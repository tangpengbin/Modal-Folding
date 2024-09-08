#include "SOBendingModeShellStructure.h"
#include <SimLib/Core/SOUtils.h>

#include <fstream>
#include <stdlib.h>
#include <random>

#include <sys/types.h>
#include <sys/stat.h>
#ifdef _WIN32
#include <direct.h>
#elif __APPLE__ || __linux__
#include <sys/stat.h>
#include <unistd.h>
#endif
#include <filesystem>

#include <iostream>
#include <iomanip>

#include "imgui.h"

//igl
#include <igl/edges.h>
#include <igl/per_vertex_normals.h>
#include <igl/readOBJ.h>
#include <igl/writeOBJ.h>
#include <igl/in_element.h>
#include <igl/adjacency_list.h>
#include <igl/boundary_loop.h>

#include <SimLib/Core/UnconstrainedLBFGS.h>
#include <SimLib/Core/colormap.h>

//ipc
#include "ipc/ipc.hpp"
#include "ipc/barrier/barrier.hpp"
#include <ipc/utils/world_bbox_diagonal_length.hpp>
#include <ipc/barrier/adaptive_stiffness.hpp>
#include <ipc/utils/local_to_global.hpp>
#include <ipc/distance/distance_type.hpp>
#include <ipc/distance/point_triangle.hpp>

#include <Spectra/SymEigsShiftSolver.h>
#include <Spectra/MatOp/SparseSymShiftSolve.h>

#include <Spectra/SymGEigsSolver.h>
#include <Spectra/MatOp/DenseSymMatProd.h>
#include <Spectra/MatOp/SparseCholesky.h>

#include "./Codegen/computeHingeEdgeAlighmentEnergy.h"
#include "./Codegen/computeHingeAngleAlighmentEnergy.h"
#include "./Codegen/computedTheta.h"
#include "./Codegen/computeTriangleDeformationGradient.h"

#include <thread>

#define USE_TBB
using namespace soutil;
const int dimension = 3;

SOBendingModeShellStructure::SOBendingModeShellStructure()
{
	//basic setting
	USE_IPC = true;
	USE_GRAVITY = false;
	USE_MAKEPD = true;
	m_useStaticSolver = false;

	m_newtonSolveResidual = 1e-7;//1e-9
	m_newtonMaxIterations = 100;
	simulationLog = false;

	m_cstThickness = 0.0005;//0.001
	time_passed = 0.0;
	m_dt = 0.01;

	m_fix_stiffness = 1e7;

	m_V.setZero();

	m_isRenderStrainColor = false;

	m_useShapeOperator = false;
	// clear all elements
	m_cstShellElements.clear();
	m_hingeElements.clear();

	m_useBothWireframeAndmesh = false;
	m_drawWireframe = false;

	//set initial green lagrange strain space variables
	m_newtonMaxLinesearchSteps = 20;
	m_drawFixingVertices = false;
	m_drawHardEdges = false;

	m_drawInteractivePairs = false;
	m_drawMaxCurvature = false;
	m_enableInteractiveSelection = false;

	m_updateModesToVisualization = true;

	m_freeLinearModesComputed = false;
	m_strainSpaceModesComputed = false;
	m_rotationStrainSpaceModesComputed = false;
	m_animateModes = false;
	m_useLinearModes = false;
	m_useStrainSpaceModes = false;
	m_useDevelopableStrainSpaceModes = false;
	m_modesAnimationID = 6;
	m_modesAmplitude = 25.0;

	m_computeModeLength = 1.0;
	m_computeModeStepLength = 1.0;

	m_visualizeSensitivityExploration = false;
	m_sensitivityMatrixUpdated = false;
	m_sensitivityExplorationProjected = false;

	m_visualizeTargetFoldingShape = false;
	m_showTargetSphere = false;
	m_targetSphereCenter = Eigen::Vector3d(0,0,-0.0);
	m_targetSphereRadius = 0.1;
}

SOBendingModeShellStructure::SOBendingModeShellStructure(const SOBendingModeShellStructure& other)
{
	(*this) = other;
	m_triMesh = new TriMesh;
	soutil::copyTriMesh(other.m_triMesh, m_triMesh);
	initial_mesh = other.initial_mesh;

	reconstructElementVector();

}

SOBendingModeShellStructure::~SOBendingModeShellStructure()
{
}

void SOBendingModeShellStructure::init(const TriMesh* pTriMesh, SOMaterial* mat, const TriMesh* targetpTriMesh,
	double thickness, 
	const std::vector<std::pair<int, int>>& fixVertexOnLinePairs, 
	const std::vector<std::pair<int, int>>& fixingVerticesDofs, 
	const std::vector<std::pair<int, Eigen::Vector3d>>& vertexIndex_force,
	bool useShapeOperator,
	string outputPath)
{
	m_useShapeOperator = useShapeOperator;
	m_outputPath = outputPath;
	m_material = mat;
	m_cstThickness = thickness;
	m_triMesh = new TriMesh;
	soutil::copyTriMesh(pTriMesh, m_triMesh);
	if (targetpTriMesh)
	{
		m_targetTriMesh = new TriMesh;
		soutil::copyTriMesh(targetpTriMesh, m_targetTriMesh);
	}
	initial_mesh = pTriMesh;
	m_fixVertexOnLinePairs = fixVertexOnLinePairs;
	for (int i = 0; i < fixingVerticesDofs.size(); i++)
		m_fixingDofs.push_back(fixingVerticesDofs[i].first * dimension + fixingVerticesDofs[i].second);
	m_vertexIndex_force = vertexIndex_force;
	m_meshDrawer->setMesh(m_triMesh);

	computeEdges();

	{
		m_nv = m_triMesh->numVertices();
		soutil::copyTriMesh(pTriMesh, m_triMesh);
		//copy positions from mesh
		m_x = m_triMesh->vertices();    //this is only valide for the initial, after simulation the vertices will change
		m_X = m_x;
		m_para_X = m_X;
		m_V = Eigen::VectorXd(3 * m_nv); //set for initial velocity
		m_V.setZero();			 //!!!!!!now this velocity is zero, TODO: parameters
	}
	//m_mesh = new DSMesh;
	//m_mesh->init(m_triMesh);

	m_faceColors = m_triMesh->faceColors();

	m_elements.clear();
	m_cstShellElements.clear();
	m_hingeElements.clear();

	initElements(m_para_X);

	reconstructMassStuff(m_para_X);

	m_v = Eigen::VectorXd(3 * m_nv); m_v.setZero();
	m_appliedForces = Eigen::VectorXd(3 * m_nv); m_appliedForces.setZero();

	m_sensitivityAnalysisProjectedState = m_para_X;

	/*for (int i = 0; i < m_x.size() / 3; i++)
	{
		if (m_x.segment<3>(3 * i).norm() < 0.02)
			m_x[3 * i + 1] += 1e-6;
	}*/
	setVisualizationConfiguration(m_x);
	
	
	initIPCConfig();
	spdlog::set_level(spdlog::level::info);
}

void SOBendingModeShellStructure::initIPCConfig()
{
	std::vector<int> faces = m_triMesh->faces();
	int numFaces = faces.size() / 3;
	Eigen::MatrixXi selfContact_triFaces(numFaces, 3);
	for (unsigned int i = 0; i < numFaces; ++i)
	{
		const unsigned int& indexA = m_triMesh->faces()[i * 3 + 0];
		const unsigned int& indexB = m_triMesh->faces()[i * 3 + 1];
		const unsigned int& indexC = m_triMesh->faces()[i * 3 + 2];

		selfContact_triFaces.row(i) = Eigen::RowVector3i(indexA, indexB, indexC);
	}

	//so we use CCD in IPC
	Eigen::MatrixXi selfContact_triEdges;
	if (selfContact_triFaces.size())
		igl::edges(selfContact_triFaces, selfContact_triEdges);

	Eigen::MatrixXd selfContact_V = (Eigen::Map<const Eigen::MatrixXd>(m_x.data(), dimension, m_x.size() / dimension)).transpose();
	m_ipc_collisionMesh = ipc::CollisionMesh(selfContact_V, selfContact_triEdges, selfContact_triFaces);
	//auto can_collide = [this](size_t vertex1, size_t vertex2) {
	//	return true;
	//};
	//m_ipc_collisionMesh.can_collide = can_collide;
	ipc::logger().set_level(spdlog::level::warn);

	if (ipc::has_intersections(m_ipc_collisionMesh, selfContact_V))
	{
		spdlog::error("The initial green lagrange strain space has intersection");
	}

	//initialize ipc adptive stiffness
	double bbox_diagonal_length = ipc::world_bbox_diagonal_length(selfContact_V);
	Eigen::VectorXd zeroEnergyGrad(m_x);
	zeroEnergyGrad.setZero();

	double averageMass = m_mass.sum() / (dimension * m_nv);// *1000.0;

	m_ipc_dhat = 1e-4;
	m_ipc_adaptiveStiffness = ipc::initial_barrier_stiffness(bbox_diagonal_length, m_ipc_dhat, averageMass, zeroEnergyGrad, zeroEnergyGrad, m_ipc_upperBoundStiffness);
	//m_ipc_adaptiveStiffness *= 1e2;
	m_ipc_upperBoundStiffness *= 1e4;
	spdlog::info("The initial ipc adaptive stiffness is {} with upper bound of {}", m_ipc_adaptiveStiffness, m_ipc_upperBoundStiffness);

}

void SOBendingModeShellStructure::reconstructMassStuff(const Eigen::VectorXd& vX)
{
	computeMasses(vX);
	m_gravity = Eigen::VectorXd(m_nv * dimension); m_gravity.setZero();
	if (USE_GRAVITY)
	{
		for (int i = 0; i < m_nv; i++)
		{
			m_gravity[3 * i + 1] = -m_g * m_mass[3 * i + 1];//gravity only set to vertices at y direction
		}
	}
}

void SOBendingModeShellStructure::initElements(const Eigen::VectorXd& vX)
{
	//assert(vX == m_X);
	m_cstShellElements.clear();

	{
		const std::vector<int>& faces = m_triMesh->faces();
		int nf = faces.size() / 3;
		for (int i = 0; i < nf; i++)
		{
			std::vector<int> inds;
			for (int j = 0; j < 3; j++)
				inds.push_back(faces[3 * i + j]);
			//cout << " triangle " << inds[0] << " " << inds[1] << " " << inds[2] << endl;
			m_cstShellElements.push_back(CSTShellElement());
			CSTShellElement& elem = m_cstShellElements.back();

			std::vector<int> dofIndices;
			for (int dof_i = 0; dof_i < inds.size(); dof_i++)
			{
				for (int j = 0; j < dimension; j++)
					dofIndices.push_back(inds[dof_i] * dimension + j);
			}

			elem.init(dofIndices, dofIndices, vX, m_material);
			elem.setThickness(m_cstThickness);
			elem.setElementProperty(InternalElement);

		}
	}

	if(m_useShapeOperator)
	{
		m_koiterShapeOperatorRestPoseBendingElements.clear();

		//for this trimesh, we find what are surrounding vertices for a triangle in the mesh
		m_triMesh->buildTriangleSurroundingTriangleStructure();
		//then, build the elements

		const auto& triangle_surroundingVertices = m_triMesh->face_ordered_FaceV_surroundingV();
		for (int face_i = 0; face_i < triangle_surroundingVertices.size(); face_i++)
		{
			const auto& tri_v_indices = triangle_surroundingVertices[face_i].first;
			const auto& tri_surroundingVerticesIndices = triangle_surroundingVertices[face_i].second;

			//the order is following, we will duplicate vertices
			//		 x4
			//     /   \
			//    x3---x2
			//   /  \ /  \
			//  x5---x1---x6
			std::vector<int> dofIndices = {
				3 * tri_v_indices[0], 3 * tri_v_indices[0] + 1, 3 * tri_v_indices[0] + 2,
				3 * tri_v_indices[1], 3 * tri_v_indices[1] + 1, 3 * tri_v_indices[1] + 2,
				3 * tri_v_indices[2], 3 * tri_v_indices[2] + 1, 3 * tri_v_indices[2] + 2 };
			for (int v_i = 0; v_i < tri_surroundingVerticesIndices.size(); v_i++)
			{
				if (tri_surroundingVerticesIndices[v_i] != -1)
					for (int d = 0; d < 3; d++)
						dofIndices.push_back(tri_surroundingVerticesIndices[v_i] * 3 + d);
				else
					break;
			}

			m_koiterShapeOperatorRestPoseBendingElements.push_back(koiterShapeOperatorRestPoseBendingElement());
			koiterShapeOperatorRestPoseBendingElement& elem = m_koiterShapeOperatorRestPoseBendingElements.back();

			elem.init(dofIndices, dofIndices, vX, m_material, m_cstThickness);
			//double energy = elem.computeEnergy(vX, vX);
			//Eigen::Matrix2d shapeOperator = elem.getCurrentShapeOperator(vX, vX);
			//spdlog::info("triangle {} shape operator bending energy {} with shape operator norm {}", face_i, energy, shapeOperator.norm());
			//std::vector<Eigen::Matrix2d> dShapeOperatordx;
			//elem.computedShapeOperatordx(vX, vX, dShapeOperatordx);

			/*if (dofIndices.size() == 18)
			{
				double theta1, theta2, theta3;
				elem.computeAngleAndShapeOperator(vX, vX, theta1, theta2, theta3);
			}*/
		}
	}
	else
	{
		m_triMesh->loadFromFile_objHardEdge();

		std::vector<Eigen::Vector2i> hard_edges = m_triMesh->hardHinge_vIndices();// { {0, 4}, { 5,12 }, { 6,14 }, { 0,15 }, { 14,15 }, { 8,14 }, { 3,15 }, { 11,13 }, { 9,13 } };
		
		m_hingeElements.clear();

		//Bending elements for bending
		//m_triMesh->buildCoHingeStructure();
		m_triMesh->buildHingeStructure();
		auto hinges = m_triMesh->hinges();
		for (int i = 0; i < hinges.size(); i++)
		{
			auto hinge = hinges[i];
			if ((hinge.tris[0] == -1) || (hinge.tris[1] == -1))
				continue; //skip boundary edges

			std::vector<int> inds(4);
			
			inds[0] = hinge.edge[0];
			inds[1] = hinge.edge[1];
			inds[2] = hinge.flaps[0];
			inds[3] = hinge.flaps[1];
			//printf("edge %d %d face %d %d\n", inds[0], inds[1], inds[2], inds[3]);

			m_hingeElements.push_back(DSHingeFixedRestPoseElement());
			DSHingeFixedRestPoseElement& elem = m_hingeElements.back();
			
			std::vector<int> dofIndices;
			for (int dof_i = 0; dof_i < inds.size(); dof_i++)
			{
				for (int j = 0; j < dimension; j++)
					dofIndices.push_back(inds[dof_i] * dimension + j);
			}

			bool isHardEdge = false;
			if (std::find(hard_edges.begin(), hard_edges.end(), Eigen::Vector2i(inds[0], inds[1])) != hard_edges.end()
				|| std::find(hard_edges.begin(), hard_edges.end(), Eigen::Vector2i(inds[1], inds[0])) != hard_edges.end())
				isHardEdge = true;
			if (isHardEdge)
				spdlog::info("Harder Edge {}-{}", inds[0], inds[1]);
			elem.init(dofIndices, dofIndices, vX, m_material, m_cstThickness, isHardEdge);// 
			//elem.setThickness(m_cstThickness);
			elem.setElementProperty(InternalElement);
			
			//spdlog::info("Theta of hinge is {}", elem.getCurrentAngle(vX));
		}

	}

	reconstructElementVector();

	//printf("init %d with shell: %d  bending: %d  motor: %d\n", m_elements.size(), m_cstShellElements.size(), m_dsBendingElements.size(), m_dihedralElements.size());
}

void SOBendingModeShellStructure::reconstructElementVector()
{
	m_elements.clear();

	for (int i = 0; i < (int)m_cstShellElements.size(); i++)
		m_elements.push_back(&m_cstShellElements[i]);
	if(m_useShapeOperator)
		for (int i = 0; i < (int)m_koiterShapeOperatorRestPoseBendingElements.size(); i++)
			m_elements.push_back(&m_koiterShapeOperatorRestPoseBendingElements[i]);
	else
		for (int i = 0; i < (int)m_hingeElements.size(); i++)
			m_elements.push_back(&m_hingeElements[i]);
	
	for (int i = 0; i < (int)m_fixDoFElements.size(); i++)
		m_elements.push_back(&m_fixDoFElements[i]);
}


void SOBendingModeShellStructure::computeEdges()
{
	//create mid edge nodes
	//1. collect all edges
	const std::vector<int>& faces = m_triMesh->faces();
	std::map<std::pair<int, int>, int> edgeMap;
	//number of 'linear' vertices
	int nlv = m_triMesh->numVertices();

	int iEdge = 0;
	int nFaces = faces.size() / 3;
	for (int i = 0; i < nFaces; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			int i1 = faces[3 * i + j];
			int i2 = faces[3 * i + ((j + 1) % 3)];
			if (i2 > i1)
				std::swap(i1, i2);

			auto eit = edgeMap.find(std::make_pair(i1, i2));
			int ie = -1;
			if (eit == edgeMap.end())
			{
				ie = iEdge++;
				edgeMap[std::make_pair(i1, i2)] = ie;
			}
			else
				ie = eit->second;
		}
	}
	m_edges.resize(2 * edgeMap.size());
	auto eit = edgeMap.begin();
	int edgeCount = 0;
	for (eit; eit != edgeMap.end(); eit++)
	{
		auto indPair = eit->first;
		m_edges[2 * edgeCount] = eit->first.first;
		m_edges[2 * edgeCount + 1] = eit->first.second;
		edgeCount++;
	}
}

void SOBendingModeShellStructure::computeMasses(const Eigen::VectorXd& vX)
{
	//const P3DIndexArray* pIB = m_triMesh->GetIndexArray();
	int nTriangles = m_triMesh->numTriangles();
	m_mass = Eigen::VectorXd(m_nv * 3); m_mass.setZero();

	for (int i = 0; i < nTriangles; i++)
	{
		V3D p[3];
		for (int j = 0; j < 3; j++)
			p[j] = get3D(m_triMesh->faces()[3 * i + j], vX);
		double area = 0.5 * crossProd3D((p[1] - p[0]), (p[2] - p[0])).norm();
		//A += area;
		//printf("%d %d %d %.15f\n", m_triMesh->faces()[3 * i + 0], m_triMesh->faces()[3 * i + 1], m_triMesh->faces()[3 * i + 2], area);
		for (int j = 0; j < 3; j++)
			for (int k = 0; k < 3; k++)
				m_mass[m_triMesh->faces()[3 * i + j] * 3 + k] += 1.0 / 3.0 * area * m_cstThickness * m_material->rho;
	}

	spdlog::info("Overall mass is {}", m_mass.sum() / 3.0);

}


bool SOBendingModeShellStructure::staticSolve_newton(Eigen::VectorXd& vx, const Eigen::VectorXd& vX)
{

	auto stepBeginningFunction = [this, vX](const Eigen::VectorXd& vx, int currentIteration)
	{
		if (USE_IPC)
		{
			ipc::Constraints constraint_set;

			Eigen::MatrixXd V = (Eigen::Map<const Eigen::MatrixXd>(vx.data(), dimension, vx.size() / dimension)).transpose();
			constraint_set.build(m_ipc_collisionMesh, V, m_ipc_dhat);

			//in the gradient computation, we will update barrier stiffness
			//since objective and hessian are computed afterward in the same step
			double min_distance = ipc::compute_minimum_distance(m_ipc_collisionMesh, V, constraint_set);
			spdlog::debug("Minimum distance is {} with number of constraints {}", min_distance, constraint_set.size());

			double bbox_diagonal_length = ipc::world_bbox_diagonal_length(V);
			m_ipc_adaptiveStiffness = ipc::update_barrier_stiffness(m_ipc_minDistance, min_distance,
				m_ipc_upperBoundStiffness, m_ipc_adaptiveStiffness, bbox_diagonal_length, 1e-5);
			spdlog::debug("Update barrier stiffness to {} with upper bound {}", m_ipc_adaptiveStiffness, m_ipc_upperBoundStiffness);
			m_ipc_minDistance = min_distance;
		}
		visualize_m_x = vx;
	};

	auto computeStaticObjective = [this, vX](const Eigen::VectorXd& vx)
	{
		//internal energy
		tbb::enumerable_thread_specific<double> storage(0);
		tbb::parallel_for(
			tbb::blocked_range<size_t>(size_t(0), m_elements.size()),
			[&](const tbb::blocked_range<size_t>& r) {
				auto& local_Eint = storage.local();
				for (size_t i = r.begin(); i < r.end(); i++)
					//for (size_t i = 0; i < m_elements.size(); i++)
				{
					local_Eint += m_elements[i]->computeEnergy(vx, vX);
					//elementEnergy[i] = m_elements[i]->computeEnergy(vx, vX);
				}
			});
		double Eint = 0.0;
		for (const auto& local_Eint : storage) {
			Eint += local_Eint;
		}

		//gravity energy
		double Egrav = 0.0;
		if (USE_GRAVITY)
		{
			Egrav += -(computeGravityForce().dot(vx));
		}

		double EappliedForces = 0.0;
		EappliedForces -= m_appliedForces.dot(vx - vX);


		//global contact barrier
		double EselfContact = 0.0;
		if(USE_IPC)
		{
			ipc::Constraints constraint_set;
			Eigen::MatrixXd V = (Eigen::Map<const Eigen::MatrixXd>(vx.data(), dimension, vx.size() / dimension)).transpose();
			constraint_set.build(m_ipc_collisionMesh, V, m_ipc_dhat);
			double barrierPotential = ipc::compute_barrier_potential(m_ipc_collisionMesh, V, constraint_set, m_ipc_dhat);

			EselfContact += m_ipc_adaptiveStiffness * barrierPotential;
		}

		double E = Eint + Egrav + EappliedForces + EselfContact;

		return E;

	};

	auto computeStaticGradient = [this, vX](const Eigen::VectorXd& vx, Eigen::VectorXd& grad)
	{
		grad = Eigen::VectorXd(vx.size());
		grad.setZero();

		//compute elements gradient
		tbb::enumerable_thread_specific<Eigen::VectorXd> storage(
			Eigen::VectorXd::Zero(vx.size()));

		tbb::parallel_for(
			tbb::blocked_range<size_t>(size_t(0), m_elements.size()),
			[&](const tbb::blocked_range<size_t>& r) {
				auto& local_grad = storage.local();
				for (size_t i = r.begin(); i < r.end(); i++)
				{
					Eigen::VectorXd element_grad;
					m_elements[i]->computeGradient(vx, vX, element_grad);

					const auto element_dofs = m_elements[i]->getDofIndices();
					for (int dof_i = 0; dof_i < element_dofs.size(); dof_i++)
					{
						local_grad[element_dofs[dof_i]] += element_grad[dof_i];
					}
				}
			});
		for (const auto& local_grad : storage) {
			grad += local_grad;
		}

		if (USE_GRAVITY)
		{
			grad += -computeGravityForce();
		}
		grad += -m_appliedForces;


		//Eigen::VectorXd strainLimitingGrad(vx.size());
		//computeAllStrainLimitingPotentialGradient(vx, vX, strainLimitingGrad);
		Eigen::VectorXd EselfContactGrad(vx.size()); EselfContactGrad.setZero();
		if(USE_IPC)
		{
			ipc::Constraints constraint_set;
			Eigen::MatrixXd V = (Eigen::Map<const Eigen::MatrixXd>(vx.data(), dimension, vx.size() / dimension)).transpose();
			constraint_set.build(m_ipc_collisionMesh, V, m_ipc_dhat);
			Eigen::VectorXd barrierPotentialGradient = ipc::compute_barrier_potential_gradient(m_ipc_collisionMesh, V, constraint_set, m_ipc_dhat);

			EselfContactGrad += m_ipc_adaptiveStiffness * barrierPotentialGradient;
		}

		grad += EselfContactGrad;

		filter(grad);
	};

	auto computeStaticHessian = [this, vX](const Eigen::VectorXd& vx, SpMat& Hessian)
	{
		Hessian = SpMat(vx.size(), vx.size());
		Hessian.setZero();

		tbb::enumerable_thread_specific<TripVec> storage;
		tbb::parallel_for(
			tbb::blocked_range<size_t>(size_t(0), m_elements.size()),
			[&](const tbb::blocked_range<size_t>& r) {
				auto& local_hess_triplets = storage.local();

				for (size_t i = r.begin(); i < r.end(); i++)
				{
					Eigen::MatrixXd element_hess;
					m_elements[i]->computeHessian(vx, vX, element_hess);

					soutil::makePD<double, Eigen::Dynamic>(element_hess, USE_MAKEPD);

					const std::vector<int>& element_dofs = m_elements[i]->getDofIndices();
					for (int dof_i = 0; dof_i < element_dofs.size(); dof_i++)
					{
						for (int dof_j = 0; dof_j < element_dofs.size(); dof_j++)
						{
							local_hess_triplets.emplace_back(element_dofs[dof_i], element_dofs[dof_j],
								element_hess(dof_i, dof_j));
						}
					}
				}
			});
		for (const auto& local_hess_triplets : storage) {
			Eigen::SparseMatrix<double> local_hess(vx.size(), vx.size());
			local_hess.setFromTriplets(
				local_hess_triplets.begin(), local_hess_triplets.end());
			Hessian += local_hess;
		}


		SpMat EselfContactHessian(vx.size(), vx.size());
		if(USE_IPC)
		{
			ipc::Constraints constraint_set;
			Eigen::MatrixXd V = (Eigen::Map<const Eigen::MatrixXd>(vx.data(), dimension, vx.size() / dimension)).transpose();
			constraint_set.build(m_ipc_collisionMesh, V, m_ipc_dhat);
			SpMat barrierPotentialHessian = ipc::compute_barrier_potential_hessian(m_ipc_collisionMesh, V, constraint_set, m_ipc_dhat, USE_MAKEPD);

			EselfContactHessian += m_ipc_adaptiveStiffness * barrierPotentialHessian;
		}

		Hessian += EselfContactHessian;

		filterMat(Hessian);
	};

	auto computeLargestStepSize = [this, vX](const Eigen::VectorXd& vx, const Eigen::VectorXd& dx)
	{
		spdlog::debug("Try to find a largest step size");


		double largestStepSize = 1.0;
		if(USE_IPC)
		{
			//we use CCD to find the largest step size
			Eigen::MatrixXd V0 = (Eigen::Map<const Eigen::MatrixXd>(vx.data(), dimension, vx.size() / dimension)).transpose();

			Eigen::MatrixXd V1 = (Eigen::Map<const Eigen::MatrixXd>((vx + largestStepSize * dx).eval().data(), dimension, vx.size() / dimension)).transpose();

			//constrain initial step size
			//while (strainSpaceAllVertices_x1.bottomRows(m_cstShellElements.size()).maxCoeff() > fixedVertices.maxCoeff())
			/*while (V1.rowwise().norm().maxCoeff() > V0.rowwise().norm().maxCoeff())
			{
				largestStepSize *= 0.5;
				V1 = (Eigen::Map<const Eigen::MatrixXd>((vx + largestStepSize * dx).eval().data(), dimension, vx.size() / dimension)).transpose();
			}*/
			while ((largestStepSize * dx).norm() > ipc::world_bbox_diagonal_length(V0))
			{
				//spdlog::debug("Search dir length {} and bbox diagonal length {}", (largestStepSize* dx).norm(), ipc::world_bbox_diagonal_length(V0));
				largestStepSize *= 0.5;
				V1 = (Eigen::Map<const Eigen::MatrixXd>((vx + largestStepSize * dx).eval().data(), dimension, vx.size() / dimension)).transpose();
			}

			spdlog::debug("Use CCD to find a largest step size with initial step size {}, in-plane", largestStepSize);

			double ipc_collisionFreeLargestStepSize = ipc::compute_collision_free_stepsize(m_ipc_collisionMesh,
				V0, V1);
			largestStepSize *= ipc_collisionFreeLargestStepSize;
			/*V1 = (Eigen::Map<const Eigen::MatrixXd>((vx + largestStepSize * dx).eval().data(), dimension, vx.size() / dimension)).transpose();

			//one more step to make sure no collision within this step
			while (!ipc::is_step_collision_free(m_ipc_collisionMesh, V0, V1))
			{
				largestStepSize *= 0.5;
				V1 = (Eigen::Map<const Eigen::MatrixXd>((vx + largestStepSize * dx).eval().data(), dimension, vx.size() / dimension)).transpose();
			}*/
			spdlog::debug("Compute largest step size done with step size {}, in-plane", largestStepSize);
		}
		return largestStepSize;
	};


	std::unique_ptr<SimOpt::UnconstrainedNewtonLinearSolver> linearSolver(new SimOpt::UnconstrainedNewtonLinearSolverLLT(computeStaticHessian));
	SimOpt::UnconstrainedNewton solver(std::move(linearSolver));
	solver.setX(vx);
	solver.setObjective(computeStaticObjective);
	solver.setObjectiveGradient(computeStaticGradient);
	solver.setLargestLineSearchStepFunctor(computeLargestStepSize);
	solver.setStepBeginningFunctor(stepBeginningFunction);
	solver.setMaxIter(5000);//m_newtonMaxIterations
	//solver.setLinesearchMaxIter(m_newtonMaxLinesearchSteps);
	solver.setGradientThreshold(1e-7);//m_newtonSolveResidual
	solver.solve();

	if (solver.state() != SimOpt::UnconstrainedNewton::SOLVED)
	{
		printf("Static solve didn't converge to %.10f\n", m_newtonSolveResidual);//, the rhs is %.10f , solver.getLastGradientNorm()
	}
	
	//after taking a step, we update states
	m_x = solver.getX();
	vx = m_x;
	return true;
}

void SOBendingModeShellStructure::computeLinearModes()
{
	if (m_freeLinearModesComputed)
		return;
	printf("Computing linear modes... ");
	//int nmodes = m_x.size();
	int nmodes = 100 < m_x.size() ? 100 : m_x.size();
	eigen_modal_analysis.resize(nmodes);
	for (int i = 0; i < nmodes; i++)
		eigen_modal_analysis[i].second = Eigen::VectorXd(m_x.size());


	computeHessian(m_x, m_para_X, m_A);
	m_A *= -1.0;
	spdlog::debug("Norm of hessian {}", m_A.norm());

	bool USE_SPECTRA = m_x.size() < 100 ? false : true;
	bool considerMass = true;

	if (considerMass)
	{
		if (!USE_SPECTRA)
			//if(true)
		{
				
			Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> generalizedEigSolver(m_A.toDense(), computeMassMatrix().toDense());

			for (int i = 0; i < nmodes; i++)
				//for (int i = nmodes-1; i >=0 ; i--)
			{
				eigen_modal_analysis[i].second = Eigen::VectorXd(generalizedEigSolver.eigenvectors().col(i));
				eigen_modal_analysis[i].first = generalizedEigSolver.eigenvalues()[i];
				//m_modes[nmodes-i-1]= dVector(eigSolver.eigenvectors().col(i));
				//m_lambdas[nmodes-i-1] = eigSolver.eigenvalues()[i];
			}
			printf("done\n");
			for (int i = 0; i < (15>nmodes ? nmodes : 15); i++)
				printf("Eigenvalue %d=%e\n", i, eigen_modal_analysis[i].first);
		}
		else
		{
			//nmodes = 15;
			Spectra::SparseSymMatProd<double> op(m_A);
			Spectra::SparseCholesky<double>  Bop(computeMassMatrix());
			Spectra::SymGEigsSolver<Spectra::SparseSymMatProd<double>, Spectra::SparseCholesky<double>, Spectra::GEigsMode::Cholesky>
				geigs(op, Bop, nmodes, 2 * nmodes);
			
			geigs.init();
			geigs.compute(Spectra::SortRule::SmallestAlge);
			if (geigs.info() == Spectra::CompInfo::Successful)
			{
				Eigen::VectorXd evalues = geigs.eigenvalues();
			}
			int neigs = geigs.eigenvalues().size();
			for (int i = 0; i < nmodes; i++)
			{
				eigen_modal_analysis[nmodes - i - 1].first = geigs.eigenvalues()[i];
				eigen_modal_analysis[nmodes - i - 1].second = Eigen::VectorXd(geigs.eigenvectors().col(i));
				spdlog::info("Eigen mode {} with eigenvalue {}", nmodes - i - 1, eigen_modal_analysis[nmodes - i - 1].first);
			}
		}
	}
	else
	{
		if (!USE_SPECTRA)
			//if(true)
		{
			Eigen::SelfAdjointEigenSolver<SpMat> eigSolver(m_A);

			for (int i = 0; i < nmodes; i++)
				//for (int i = nmodes-1; i >=0 ; i--)
			{
				eigen_modal_analysis[i].second = Eigen::VectorXd(eigSolver.eigenvectors().col(i));
				eigen_modal_analysis[i].first = eigSolver.eigenvalues()[i];
				//m_modes[nmodes-i-1]= dVector(eigSolver.eigenvectors().col(i));
				//m_lambdas[nmodes-i-1] = eigSolver.eigenvalues()[i];
			}
			printf("done\n");
			//for (int i = 0; i < (15>nmodes ? nmodes : 15); i++)
			//	printf("Eigenvalue %d=%e\n", i, eigen_modal_analysis[i].first);
		}
		else
		{
			//nmodes = 15;
			Spectra::SparseSymShiftSolve<double, Eigen::Upper> op(m_A);
			double shift = -1e-1;
			Spectra::SymEigsShiftSolver<Spectra::SparseSymShiftSolve<double, Eigen::Upper> > eigs(op, nmodes, 2 * nmodes, shift);
			eigs.init();
			eigs.compute();
			if (eigs.info() == Spectra::CompInfo::Successful)
			{
				Eigen::VectorXd evalues = eigs.eigenvalues();
				//m_lambdas = evalues;
				//std::cout << "Eigenvalues found:\n" << evalues << std::endl;
			}
			int neigs = eigs.eigenvalues().size();
			for (int i = 0; i < nmodes; i++)
			{
				eigen_modal_analysis[nmodes - i - 1].first = eigs.eigenvalues()[i];
				eigen_modal_analysis[nmodes - i - 1].second = Eigen::VectorXd(eigs.eigenvectors().col(i));
				spdlog::info("Eigen mode {} with eigenvalue {}", nmodes - i - 1, eigen_modal_analysis[nmodes - i - 1].first);
			}
		}
	}
	printf("done\n");


	for (int i = 6; i < nmodes; i++)
	{
		//im_filterRigidModes(m_x, m_modes[i]);
		//dVector vBndMode(3 * m_nv);
		//vBndMode.setZero();
		//for (int j = 0; j < bndinds.size(); j++)
		//	set3D(bndinds[j], get3D(bndinds[j], m_modes[i]), vBndMode);
		//m_modes[i] = vBndMode;
		eigen_modal_analysis[i].second.normalize();
	}
	m_freeLinearModesComputed = true;
}

void SOBendingModeShellStructure::convertToStrainSpaceMode(const Eigen::VectorXd& positionSpaceMode, Eigen::VectorXd& strainSpaceMode)
{
	Eigen::VectorXd pm = positionSpaceMode;// .normalized();
	strainSpaceMode = Eigen::VectorXd(m_hingeElements.size()); strainSpaceMode.setZero();
	//for (const auto& hinge : m_hingeElements)
	for (int i = 0; i < m_hingeElements.size(); i++)
	{
		Eigen::VectorXd dthetadx(12);
		dthetadx.setZero();
		m_hingeElements[i].computedThetadx(m_X, m_X, dthetadx);
		//std::cout << "hinge " << i << " dthetadx norm " << dthetadx.norm() << std::endl;
		const auto& inds = m_hingeElements[i].getDofIndices();
		
		Eigen::VectorXd hingeElement_position = m_hingeElements[i].gatherVector(pm);
		for (int j = 0; j < 4; j++)
		{
			V3D vj = hingeElement_position.segment<3>(j * 3);
			strainSpaceMode[i] += vj.dot(get3D(j, dthetadx));
		}
	}
}


void SOBendingModeShellStructure::convertToShapeOperatorStrainSpaceMode(const Eigen::VectorXd& positionSpaceMode, std::vector<Eigen::Matrix2d>& strainSpaceMode)
{
	strainSpaceMode.resize(m_koiterShapeOperatorRestPoseBendingElements.size());
	//for (const auto& hinge : m_hingeElements)
	for (int i = 0; i < m_koiterShapeOperatorRestPoseBendingElements.size(); i++)
	{
		std::vector<Eigen::Matrix2d> dShapeOperatordx;
		m_koiterShapeOperatorRestPoseBendingElements[i].computedShapeOperatordx(m_X, m_X, dShapeOperatordx);
		//std::cout << "hinge " << i << " dthetadx norm " << dthetadx.norm() << std::endl;
		const auto& inds = m_koiterShapeOperatorRestPoseBendingElements[i].getDofIndices();

		Eigen::VectorXd hingeElement_positionMode = m_koiterShapeOperatorRestPoseBendingElements[i].gatherVector(positionSpaceMode);
		strainSpaceMode[i].setZero();
		for (int j = 0; j < inds.size(); j++)
		{
			strainSpaceMode[i] += hingeElement_positionMode[j] * dShapeOperatordx[j];
		}
	}
}



void SOBendingModeShellStructure::animateStrainSpaceMode()
{
	static bool m_modesNeedsUpdate = true;
	if (m_modesNeedsUpdate)
	{
		//computeNatureModes(m_x, m_para_X);
		computeLinearModes();
		m_modesNeedsUpdate = false;
	}
	const Eigen::VectorXd& vm = eigen_modal_analysis[m_modesAnimationID].second;
	//dVector vm(m_modes[m_modesAnimationID]);
	//vm += m_modes[m_modesAnimationID + 1];
	//vm.normalize();
	if(m_useShapeOperator)
	{
		std::vector<Eigen::Matrix2d> ssm;
		convertToShapeOperatorStrainSpaceMode(vm, ssm);
		double alpha = sin(time_passed) * m_modesAmplitude;
		for (int m_i = 0; m_i < ssm.size(); m_i++)
			ssm[m_i] *= alpha;

		//ssm_setRestAngles(ssm);
		for (int i = 0; i < m_koiterShapeOperatorRestPoseBendingElements.size(); i++)
		{
			//std::cout << "hinge " << i << m_hingeElements[i].getRestAngle();
			Eigen::Matrix2d shapeOperator = m_koiterShapeOperatorRestPoseBendingElements[i].getCurrentShapeOperator(m_x,m_para_X);
			shapeOperator += ssm[i];
			m_koiterShapeOperatorRestPoseBendingElements[i].setRestShapeOperator(shapeOperator);
			//std::cout << " update theta0 as " << theta << " ssm " << ssm[i] << std::endl;
		}
	}
	else
	{
		Eigen::VectorXd ssm;
		convertToStrainSpaceMode(vm, ssm);
		double alpha = sin(time_passed) * m_modesAmplitude;
		ssm *= alpha;
		
		//ssm_setRestAngles(ssm);
		for (int i = 0; i < m_hingeElements.size(); i++)
		{
			//std::cout << "hinge " << i << m_hingeElements[i].getRestAngle();
			double theta = m_hingeElements[i].getCurrentAngle(m_x);
			theta += ssm[i];
			m_hingeElements[i].setRestAngle(theta);
			//std::cout << " update theta0 as " << theta << " ssm " << ssm[i] << std::endl;
		}
	}

	projectToManifold();
	setVisualizationConfiguration(m_x);

	time_passed += m_dt;
}

void SOBendingModeShellStructure::projectToManifold()
{
	bool use_gravity = USE_GRAVITY;
	USE_GRAVITY = false;
	if (!staticSolve_newton(getCurrentPositions(), getParamPosition()))
		cout << "The static solve doesn't reach equilibrim state, the clamped modes will have some problems" << endl;
	else
		cout << "The static solve has been correctly done!" << endl;

	USE_GRAVITY = use_gravity;

	/*if (m_fixDoFElements.size() == 0)
	{
		//move to center
		V3D center(0, 0, 0);
		for (int i = 0; i < m_x.size() / 3; i++)
		{
			center += m_x.segment<3>(i * 3);
		}
		center /= (double)(m_x.size() / 3);
		//V3D center = m_x.segment<3>(2 * 3);

		for (int i = 0; i < m_x.size() / 3; i++)
			m_x.segment<3>(i * 3) -= center;
	}*/
	//setVisualizationConfiguration(m_x);

	spdlog::info("Projection done");

	//{
	//	//compute in-plane energy
	//	double energy = 0.0;
	//	Eigen::VectorXd maxStrain(m_cstShellElements.size());
	//	for (int i = 0; i < m_cstShellElements.size(); i++)
	//	{
	//		energy += m_cstShellElements[i].computeEnergy(m_x, m_X);
	//		m_cstShellElements[i].computeMaxStrain(m_x, m_X, maxStrain[i]);
	//	}
	//	spdlog::info("After projection, the in-plane energy is {} and max streching over all triangle is {}", energy, maxStrain.maxCoeff());
	//}
}

void SOBendingModeShellStructure::computeBatchStrainSpaceMode()
{
	if (m_strainSpaceModesComputed)
		return;
	static bool useConstrainedMinEnergy = false;
	if(m_useShapeOperator)
	{
		computeLinearModes();
		m_shapeOperatorStrainSpaceModes.resize(eigen_modal_analysis.size());
		
		for (int i = 6; i < eigen_modal_analysis.size(); i++)
		{//compute all strain space modes
			
			m_shapeOperatorStrainSpaceModes[i].resize(m_koiterShapeOperatorRestPoseBendingElements.size());
			convertToShapeOperatorStrainSpaceMode(eigen_modal_analysis[i].second.normalized(), m_shapeOperatorStrainSpaceModes[i]);
			
			double normOfShapeOperator = 0;
			for (int j = 0; j < m_shapeOperatorStrainSpaceModes[i].size(); j++)
			{
				normOfShapeOperator += m_shapeOperatorStrainSpaceModes[i][j].norm();
				//std::cout << "elem " << j << "eigen values " << m_shapeOperatorStrainSpaceModes[i][j].eigenvalues().transpose() << std::endl;
			}

			spdlog::info("comptue strain space mode {} with norm {}", i, normOfShapeOperator);
		}
	}
	else
	{
		computeLinearModes();
		m_strainSpaceModes.resize(eigen_modal_analysis.size());
		
		for (int i = 0; i < eigen_modal_analysis.size(); i++)
		{//compute all strain space modes
			m_strainSpaceModes[i].resize(m_hingeElements.size());
			
			convertToStrainSpaceMode(eigen_modal_analysis[i].second.normalized(), m_strainSpaceModes[i]);
			//m_strainSpaceModes[i].normalize();
			
			spdlog::info("comptue strain space mode {} with norm of mdoe {}", i, m_strainSpaceModes[i].norm());
		}

		//std::cout << "angle space\n" << angleSpace << std::endl;

		//m_x = output_x;
		//setVisualizationConfiguration(m_x);
		//return;
	}
	m_strainSpaceModesComputed = true;
}

void SOBendingModeShellStructure::computeStrainSpaceModes()
{
	computeBatchStrainSpaceMode();

	//{
	//	if (m_x.size() == 3721 * 3)
	//	{
	//		spdlog::info("Fix center triangle for nonlinear compliant modes");
	//		m_fixingDofs = { 1860 * 3 + 0, 1860 * 3 + 1, 1860 * 3 + 2,
	//		1799 * 3 + 0, 1799 * 3 + 1,// m_triMesh->faces()[1] * 3 + 2,
	//		1859 * 3 + 0//, m_triMesh->faces()[2] * 3 + 1, m_triMesh->faces()[2] * 3 + 2
	//		};
	//	}
	//	m_fixDoFElements.clear();
	//	for (int i = 0; i < m_fixingDofs.size(); i++)
	//	{
	//		m_fixDoFElements.push_back(FixDoFElement());
	//		FixDoFElement& elem = m_fixDoFElements.back();

	//		elem.init(m_fixingDofs[i], m_fixingDofs[i], m_fix_stiffness, true);
	//		elem.setElementProperty(ExternalElement);

	//	}
	//	for (int i = 0; i < m_fixDoFElements.size(); i++)
	//	{
	//		m_elements.push_back(&m_fixDoFElements[i]);
	//	}
	//}

	if (m_computeModeLength <= 0.0)
	{
		spdlog::warn("The input length should larger than 0.0");
		return;
	}
	auto computeStrainSpaceMode = [this](double length, const Eigen::VectorXd& hingeElementRestAngle, const Eigen::VectorXd& initialGuess)
	{
		if (m_useStrainSpaceModes)
		{
			if (m_strainSpaceModes.size() >= m_modesAnimationID && m_strainSpaceModes[m_modesAnimationID].size() > 0)
			{
				for (int i = 0; i < m_hingeElements.size(); i++)
				{
					m_hingeElements[i].setRestAngle(hingeElementRestAngle[i]);
					//std::cout << "hinge " << i << " restAngle " << m_hingeElements[i].getRestAngle() << std::endl;
				}
				spdlog::info("Max angle is {}", hingeElementRestAngle.cwiseAbs().maxCoeff());

				m_x = initialGuess;
				projectToManifold();
				//move to center
				V3D center(0, 0, 0);
				for (int i = 0; i < m_x.size() / 3; i++)
				{
					center += m_x.segment<3>(i * 3);
				}
				center /= (double)(m_x.size() / 3);
				//V3D center = m_x.segment<3>(2 * 3);

				for (int i = 0; i < m_x.size() / 3; i++)
					m_x.segment<3>(i * 3) -= center;
				setVisualizationConfiguration(m_x);
			}
		}
		else if (m_useLinearModes)
		{
			if (eigen_modal_analysis.size() >= m_modesAnimationID)
			{
				m_x = m_para_X + length * eigen_modal_analysis[m_modesAnimationID].second;

				setVisualizationConfiguration(m_x);
			}
		}
	};



	Eigen::VectorXd initialGuess = m_para_X;
	double initialLength = 0.0;
	if (m_strainSpaceModes_index_to_lengthState.find(m_modesAnimationID) != m_strainSpaceModes_index_to_lengthState.end()
		&& m_strainSpaceModes_index_to_lengthState.at(m_modesAnimationID).size() > 0)
	{
		//initialLength = m_strainSpaceModes_index_to_lengthState[m_modesAnimationID].back().first;
		//initialGuess = m_strainSpaceModes_index_to_lengthState[m_modesAnimationID].back().second;
			
		//we first compute how many step needs be added
		double lastLength = m_strainSpaceModes_index_to_lengthState[m_modesAnimationID].back().first;
		double lastLengthStep = lastLength / double(m_strainSpaceModes_index_to_lengthState[m_modesAnimationID].size());

		int addSteps = lastLengthStep / m_computeModeStepLength - 1;
		std::vector<std::pair<double, Eigen::VectorXd>> new_length_states;

		Eigen::VectorXd angleStep(m_hingeElements.size());
		for (int i = 0; i < m_hingeElements.size(); i++)
		{
			double theta_next = m_hingeElements[i].getCurrentAngle(m_strainSpaceModes_index_to_lengthState[m_modesAnimationID][0].second);
			double theta_previous = m_hingeElements[i].getInitialRestAngle();
			angleStep[i] = (theta_next - theta_previous) / (addSteps + 1);
			//std::cout << "hinge " << i << " restAngle " << m_hingeElements[i].getRestAngle() << std::endl;
		}
		if (m_strainSpaceModes_index_to_lengthState[m_modesAnimationID][0].first > 0.0)
		{
			//start from 0 - first point
			for (int add_step_i = 0; add_step_i < addSteps; add_step_i++)
			{
				double length = 0.0 + (add_step_i + 1) * m_computeModeStepLength;// std::min(initialLength + i * m_computeModeStepLength, m_computeModeLength);
				spdlog::info("Compute mode {} with strain space length {}", m_modesAnimationID, length);

				Eigen::VectorXd targetAngle(angleStep.size());
				for (int i = 0; i < m_hingeElements.size(); i++)
				{
					double theta_previous = m_hingeElements[i].getInitialRestAngle();
					targetAngle[i] = theta_previous + (add_step_i + 1) * angleStep[i];
					//std::cout << "hinge " << i << " restAngle " << m_hingeElements[i].getRestAngle() << std::endl;
				}

				computeStrainSpaceMode(length, targetAngle, initialGuess);

				new_length_states.emplace_back(length, m_x);
			}
		}
		//from first point to the last second one
		for (int i = 0; i < m_strainSpaceModes_index_to_lengthState[m_modesAnimationID].size() - 1; i++)
		{
			initialLength = m_strainSpaceModes_index_to_lengthState[m_modesAnimationID][i].first;
			initialGuess = m_strainSpaceModes_index_to_lengthState[m_modesAnimationID][i].second;
			new_length_states.emplace_back(initialLength, initialGuess);

			for (int elem_i = 0; elem_i < m_hingeElements.size(); elem_i++)
			{
				double theta_next = m_hingeElements[elem_i].getCurrentAngle(m_strainSpaceModes_index_to_lengthState[m_modesAnimationID][i + 1].second);
				double theta_previous = m_hingeElements[elem_i].getCurrentAngle(m_strainSpaceModes_index_to_lengthState[m_modesAnimationID][i].second);
				angleStep[elem_i] = (theta_next - theta_previous) / (addSteps + 1);
				//std::cout << "hinge " << i << " restAngle " << m_hingeElements[i].getRestAngle() << std::endl;
			}

			for (int add_step_i = 0; add_step_i < addSteps; add_step_i++)
			{
				double length = initialLength + (add_step_i + 1) * m_computeModeStepLength;// std::min(initialLength + i * m_computeModeStepLength, m_computeModeLength);
				spdlog::info("Compute mode {} with strain space length {}", m_modesAnimationID, length);

				Eigen::VectorXd targetAngle(angleStep.size());
				for (int elem_i = 0; elem_i < m_hingeElements.size(); elem_i++)
				{
					double theta_previous = m_hingeElements[elem_i].getCurrentAngle(m_strainSpaceModes_index_to_lengthState[m_modesAnimationID][i].second);
					targetAngle[elem_i] = theta_previous + (add_step_i + 1) * angleStep[elem_i];
					//std::cout << "hinge " << i << " restAngle " << m_hingeElements[i].getRestAngle() << std::endl;
				}

				computeStrainSpaceMode(length, targetAngle, initialGuess);

				new_length_states.emplace_back(length, m_x);
			}
		}
		new_length_states.emplace_back(m_strainSpaceModes_index_to_lengthState[m_modesAnimationID].back().first, m_strainSpaceModes_index_to_lengthState[m_modesAnimationID].back().second);

		m_strainSpaceModes_index_to_lengthState[m_modesAnimationID] = new_length_states;
	}
	else
	{
		int numSteps = (m_computeModeLength - initialLength) / m_computeModeStepLength;
		for (int i = 1; i <= numSteps; i++)
		{
			double length = initialLength + i * m_computeModeStepLength;// std::min(initialLength + i * m_computeModeStepLength, m_computeModeLength);
			spdlog::info("Compute mode {} with strain space length {}", m_modesAnimationID, length);

			Eigen::VectorXd ssm = m_strainSpaceModes[m_modesAnimationID]  * length;// * m_computeModeStepLength;
			Eigen::VectorXd targetRestAngle(ssm.size());
			for (int i = 0; i < m_hingeElements.size(); i++)
			{
				double theta = m_hingeElements[i].getCurrentAngle(initialGuess);
				theta += ssm[i];
				targetRestAngle[i] = theta;
				//std::cout << "hinge " << i << " restAngle " << m_hingeElements[i].getRestAngle() << std::endl;
			}
			
			computeStrainSpaceMode(length, targetRestAngle, initialGuess);

			//SpMat dxdtheta_rest;
			//computeStrainSpaceSensitiveMatrix(m_x, m_para_X, dxdtheta_rest);

			m_strainSpaceModes_index_to_lengthState[m_modesAnimationID].emplace_back(length, m_x);
			initialGuess = m_x;
		}
	}


	//compute energy
	auto computeEnergy = [this](const Eigen::VectorXd& vx)
	{
		//internal energy
		tbb::enumerable_thread_specific<double> storage(0);
		tbb::parallel_for(
			tbb::blocked_range<size_t>(size_t(0), m_elements.size()),
			[&](const tbb::blocked_range<size_t>& r) {
				auto& local_Eint = storage.local();
				for (size_t i = r.begin(); i < r.end(); i++)
					//for (size_t i = 0; i < m_elements.size(); i++)
				{
					local_Eint += m_elements[i]->computeEnergy(vx, m_para_X);
					//elementEnergy[i] = m_elements[i]->computeEnergy(vx, vX);
				}
			});
		double Eint = 0.0;
		for (const auto& local_Eint : storage) {
			Eint += local_Eint;
		}
		return Eint;
	};

	for (int i = 0; i < m_hingeElements.size(); i++)
	{
		m_hingeElements[i].setRestAngle(m_hingeElements[i].getInitialRestAngle());
	}
	std::cout << "Internal Energy along with displacement and modal displacement" << std::endl;
	int numStates = m_strainSpaceModes_index_to_lengthState[m_modesAnimationID].size();
	Eigen::VectorXd internalEnergy(numStates);
	for (size_t i = 0; i < numStates; i++)
	{
		double maxStrain = 0;
		for (int cst_i = 0; cst_i < m_cstShellElements.size(); cst_i++)
		{
			double max_strain_i = 0;
			m_cstShellElements[cst_i].computeMaxStrain(m_strainSpaceModes_index_to_lengthState[m_modesAnimationID][i].second, m_para_X, max_strain_i);
			if (maxStrain < max_strain_i)
				maxStrain = max_strain_i;
		}

		internalEnergy[i] = computeEnergy(m_strainSpaceModes_index_to_lengthState[m_modesAnimationID][i].second);
		double displacement = (m_strainSpaceModes_index_to_lengthState[m_modesAnimationID][i].second - m_para_X).norm();
		std::cout << internalEnergy[i] << " " << displacement << " " << m_strainSpaceModes_index_to_lengthState[m_modesAnimationID][i].first << " " << maxStrain << std::endl;
	}
}

void SOBendingModeShellStructure::blendStrainSpaceModes(const std::vector<int>& selectedModes, const std::vector<double>& eachModeLength, double stepLength)
{
	computeBatchStrainSpaceMode();

	for (int i = 0; i < eachModeLength.size(); i++)
	{
		if (eachModeLength[i] <= 0.0)
		{
			spdlog::warn("The input length should larger than 0.0");
			return;
		}
	}
	
	Eigen::VectorXd initialGuess = m_para_X;
	
	std::vector<int> eachModeNumSteps(eachModeLength.size());
	std::vector<int> eachModeSliceNum(eachModeLength.size());
	int numSteps = 1;
	for (int i = 0; i < eachModeLength.size(); i++)
	{
		eachModeNumSteps[i] = int(eachModeLength[i] / stepLength) + 1;
		numSteps *= eachModeNumSteps[i];
		eachModeSliceNum[i] = numSteps;
	}
	m_selectedStrainSpce_to_eachModeLengthState[selectedModes] = std::make_pair(stepLength, std::map<std::vector<int>, Eigen::VectorXd>());

	auto sampleStrainSpace = [this,&selectedModes](std::vector<int>& modeLengthSampleIndices, double stepSize)
	{
		Eigen::VectorXd modeLength(modeLengthSampleIndices.size());
		for (int i = 0; i < modeLengthSampleIndices.size(); i++)
		{
			//we will accumulate length
			modeLength[i] = 0.0;
			for (int j = 1; j <= modeLengthSampleIndices[i]; j++)
			{
				modeLength[i] += j * stepSize;
			}
		}
		Eigen::VectorXd ssm(m_strainSpaceModes[6].size()); ssm.setZero();
			
		for (int i = 0; i < selectedModes.size(); i++)
		{
			if (m_strainSpaceModes.size() <= selectedModes[i] || selectedModes[i] < 0)
				return;

			ssm += m_strainSpaceModes[selectedModes[i]] * modeLength[i];
		}
			
		//find last state as the state initial guess
		Eigen::VectorXd initialGuess = m_para_X;
		for (int i = 0; i < modeLengthSampleIndices.size(); i++)
		{
			auto guessIndices = modeLengthSampleIndices;
			guessIndices[i] -= 1;
			const auto guessFinder = m_selectedStrainSpce_to_eachModeLengthState[selectedModes].second.find(guessIndices);
			if (guessFinder != m_selectedStrainSpce_to_eachModeLengthState[selectedModes].second.end())
			{
				std::cout << "Guess indices " << Eigen::Map<Eigen::VectorXi>(guessIndices.data(), guessIndices.size()).transpose() << std::endl;
					
				initialGuess = guessFinder->second;
				break;
			}
		}
			
		for (int i = 0; i < m_hingeElements.size(); i++)
		{
			double theta = m_hingeElements[i].getCurrentAngle(m_para_X);//here we always start from planar state
			theta += ssm[i];
			m_hingeElements[i].setRestAngle(theta);
		}

		m_x = m_para_X;
		projectToManifold();
		setVisualizationConfiguration(m_x);
			
		

		m_selectedStrainSpce_to_eachModeLengthState[selectedModes].second[modeLengthSampleIndices] = m_x;
	};

	for (int i = 1; i <= numSteps-1; i++)
	{
		std::vector<int> sample_i;
		int remain = i;
		for (int j = eachModeLength.size() - 2; j >=0; j--)
		{
			int index = remain / eachModeSliceNum[j];
			remain = remain % eachModeSliceNum[j];

			sample_i.push_back(index);
			if (j == 0)
				sample_i.push_back(remain);
		}
		
		std::cout << "Sampele "  << Eigen::Map<Eigen::VectorXi>(sample_i.data(), sample_i.size()).transpose() << std::endl;
		sampleStrainSpace(sample_i, stepLength);
	}


	//compute energy
	auto computeEnergy = [this](const Eigen::VectorXd& vx)
	{
		//internal energy
		tbb::enumerable_thread_specific<double> storage(0);
		tbb::parallel_for(
			tbb::blocked_range<size_t>(size_t(0), m_elements.size()),
			[&](const tbb::blocked_range<size_t>& r) {
				auto& local_Eint = storage.local();
				for (size_t i = r.begin(); i < r.end(); i++)
					//for (size_t i = 0; i < m_elements.size(); i++)
				{
					local_Eint += m_elements[i]->computeEnergy(vx, m_para_X);
					//elementEnergy[i] = m_elements[i]->computeEnergy(vx, vX);
				}
			});
		double Eint = 0.0;
		for (const auto& local_Eint : storage) {
			Eint += local_Eint;
		}
		return Eint;
	};

	for (int i = 0; i < m_hingeElements.size(); i++)
	{
		m_hingeElements[i].setRestAngle(m_hingeElements[i].getInitialRestAngle());
	}
	std::cout << "Internal Energy along with displacement and modal displacement" << std::endl;
	int numStates = m_selectedStrainSpce_to_eachModeLengthState[selectedModes].second.size();
	std::vector<double> internalEnergy; internalEnergy.reserve(numStates);
	for (const auto& state_i : m_selectedStrainSpce_to_eachModeLengthState[selectedModes].second)
	{
		internalEnergy.push_back(computeEnergy(state_i.second));
		double displacement = (state_i.second - m_para_X).norm();
		std::cout << internalEnergy.back() << " " << displacement << " sample id " << Eigen::Map<const Eigen::VectorXi>(state_i.first.data(), state_i.first.size()).transpose() << std::endl;
	}
}


void SOBendingModeShellStructure::computeStrainSpaceSensitiveMatrix(const Eigen::VectorXd& vx, const Eigen::VectorXd& vX, SpMat& dxdtheta_rest)
{//note: before computing this sensitivity analysis, we have to make sure that the hessian of the static equlibirum is positive definite!!!!!!!!!!!!!
	m_sensitivityAnalysisState = vx;

	//compute hessian
	auto computeStaticHessian = [this, vX](const Eigen::VectorXd& vx, SpMat& Hessian)
	{
		Hessian = SpMat(vx.size(), vx.size());
		Hessian.setZero();

		tbb::enumerable_thread_specific<TripVec> storage;
		tbb::parallel_for(
			tbb::blocked_range<size_t>(size_t(0), m_elements.size()),
			[&](const tbb::blocked_range<size_t>& r) {
				auto& local_hess_triplets = storage.local();

				for (size_t i = r.begin(); i < r.end(); i++)
				{
					Eigen::MatrixXd element_hess;
					m_elements[i]->computeHessian(vx, vX, element_hess);

					const std::vector<int>& element_dofs = m_elements[i]->getDofIndices();
					for (int dof_i = 0; dof_i < element_dofs.size(); dof_i++)
					{
						for (int dof_j = 0; dof_j < element_dofs.size(); dof_j++)
						{
							local_hess_triplets.emplace_back(element_dofs[dof_i], element_dofs[dof_j],
								element_hess(dof_i, dof_j));
						}
					}
				}
			});
		for (const auto& local_hess_triplets : storage) {
			Eigen::SparseMatrix<double> local_hess(vx.size(), vx.size());
			local_hess.setFromTriplets(
				local_hess_triplets.begin(), local_hess_triplets.end());
			Hessian += local_hess;
		}


		SpMat EselfContactHessian(vx.size(), vx.size());
		if (USE_IPC)
		{
			ipc::Constraints constraint_set;
			Eigen::MatrixXd V = (Eigen::Map<const Eigen::MatrixXd>(vx.data(), dimension, vx.size() / dimension)).transpose();
			constraint_set.build(m_ipc_collisionMesh, V, m_ipc_dhat);
			SpMat barrierPotentialHessian = ipc::compute_barrier_potential_hessian(m_ipc_collisionMesh, V, constraint_set, m_ipc_dhat, false);
			spdlog::info("Number contact {}", constraint_set.size());
			EselfContactHessian += m_ipc_adaptiveStiffness * barrierPotentialHessian;
		}

		Hessian += EselfContactHessian;

		filterMat(Hessian);
	};
	SpMat hessian;
	computeStaticHessian(vx, hessian);
	
	/*auto computeStaticGradient = [this, vX](const Eigen::VectorXd& vx, Eigen::VectorXd& grad)
		{
			grad = Eigen::VectorXd(vx.size());
			grad.setZero();

			//compute elements gradient
			tbb::enumerable_thread_specific<Eigen::VectorXd> storage(
				Eigen::VectorXd::Zero(vx.size()));

			tbb::parallel_for(
				tbb::blocked_range<size_t>(size_t(0), m_elements.size()),
				[&](const tbb::blocked_range<size_t>& r) {
					auto& local_grad = storage.local();
					for (size_t i = r.begin(); i < r.end(); i++)
					{
						Eigen::VectorXd element_grad;
						m_elements[i]->computeGradient(vx, vX, element_grad);

						const auto element_dofs = m_elements[i]->getDofIndices();
						for (int dof_i = 0; dof_i < element_dofs.size(); dof_i++)
						{
							local_grad[element_dofs[dof_i]] += element_grad[dof_i];
						}
					}
				});
			for (const auto& local_grad : storage) {
				grad += local_grad;
			}

			grad += -m_appliedForces;


			//Eigen::VectorXd strainLimitingGrad(vx.size());
			//computeAllStrainLimitingPotentialGradient(vx, vX, strainLimitingGrad);
			Eigen::VectorXd EselfContactGrad(vx.size()); EselfContactGrad.setZero();
			if (USE_IPC)
			{
				ipc::Constraints constraint_set;
				Eigen::MatrixXd V = (Eigen::Map<const Eigen::MatrixXd>(vx.data(), dimension, vx.size() / dimension)).transpose();
				constraint_set.build(m_ipc_collisionMesh, V, m_ipc_dhat);
				Eigen::VectorXd barrierPotentialGradient = ipc::compute_barrier_potential_gradient(m_ipc_collisionMesh, V, constraint_set, m_ipc_dhat);

				EselfContactGrad += m_ipc_adaptiveStiffness * barrierPotentialGradient;
			}

			grad += EselfContactGrad;

			filter(grad);
		};

	//Eigen::VectorXd f;
	//computeStaticGradient(vx, f);
	//std::cout << "Force norm " << f.norm() << std::endl;
	SimOpt::DerivativeTester tester;
	tester.testJacobian(vx, hessian, computeStaticGradient, 7, 1e-3);
	exit(0);*/

	//compute dfdtheta_rest
	SpMat dfdtheta_rest(vx.size(), m_hingeElements.size());
	dfdtheta_rest.setZero();
	{
		tbb::enumerable_thread_specific<TripVec> storage;
		tbb::parallel_for(
			tbb::blocked_range<size_t>(size_t(0), m_hingeElements.size()),
			[&](const tbb::blocked_range<size_t>& r) {
				auto& local_jacobian_triplets = storage.local();

				for (size_t i = r.begin(); i < r.end(); i++)
		//		for(int i =0;i<m_hingeElements.size();i++)
				{
					//d^2E/dxdtheta_rest
					Eigen::VectorXd element_gradientParameterJacobian;
					m_hingeElements[i].computeGradientParameterJacobian(vx, vX, element_gradientParameterJacobian);

					const std::vector<int>& element_dofs = m_hingeElements[i].getDofIndices();
					for (int dof_i = 0; dof_i < element_dofs.size(); dof_i++)
					{
						local_jacobian_triplets.emplace_back(element_dofs[dof_i], i,
							element_gradientParameterJacobian[dof_i]);
					}
				}
			});
		for (const auto& local_jacobian_triplets : storage) {
			Eigen::SparseMatrix<double> local_jacobian(vx.size(), m_hingeElements.size());
			local_jacobian.setFromTriplets(
				local_jacobian_triplets.begin(), local_jacobian_triplets.end());
			dfdtheta_rest += local_jacobian;
		}
	}

	/*{
		auto computeForce_theta = [this, vx, vX](const Eigen::VectorXd& vtheta, Eigen::VectorXd& f)
		{
			for (int i = 0; i < m_hingeElements.size(); i++)
			{
				m_hingeElements[i].setRestAngle(vtheta[i]);
			}

			f.setZero(vx.size());
			//compute elements gradient
			tbb::enumerable_thread_specific<Eigen::VectorXd> storage(
				Eigen::VectorXd::Zero(vx.size()));

			tbb::parallel_for(
				tbb::blocked_range<size_t>(size_t(0), m_elements.size()),
				[&](const tbb::blocked_range<size_t>& r) {
					auto& local_grad = storage.local();
					for (size_t i = r.begin(); i < r.end(); i++)
					{
						Eigen::VectorXd element_grad;
						m_elements[i]->computeGradient(vx, vX, element_grad);

						const auto element_dofs = m_elements[i]->getDofIndices();
						for (int dof_i = 0; dof_i < element_dofs.size(); dof_i++)
						{
							local_grad[element_dofs[dof_i]] += element_grad[dof_i];
						}
					}
				});
			for (const auto& local_grad : storage) {
				f += local_grad;
			}
			//f *= -1.0;
		};
		Eigen::VectorXd v_theta(m_hingeElements.size());
		for (int i = 0; i < m_hingeElements.size(); i++)
		{
			v_theta[i] = m_hingeElements[i].getRestAngle(); // m_hingeElements[i].getCurrentAngle(vx);
		}
		SimOpt::DerivativeTester tester;
		tester.testJacobian(v_theta, dfdtheta_rest, computeForce_theta, 7, 1e-4);
		exit(0);
	}*/


	// f(x(theta_rest)) = 0
	// df/dx * dx/dtheta_rest + df/dtheta_rest = 0
	// dx/dtheta_rest = (df/dx)^-1 * -df/dtheta_rest
	//then, we solve the system

	//solve for dzdw
	//Eigen::PardisoLDLT<SpMat> solver(hessian);
	Eigen::SimplicialCholesky<SpMat> solver(hessian);

	//test if this is a sime-positive definite matrix
	if (solver.info() != Eigen::Success)
	{
		/// decomposition failed
		spdlog::warn("dfdx decomposition failed when computing predictor");
	}
	dxdtheta_rest = solver.solve(-dfdtheta_rest);
	//check factorization
	if (solver.info() != Eigen::Success)
	{
		/// decomposition failed
		spdlog::warn("dfdx solve failed when computing predictor");
	}

	//check redsiduum
	SpMat r = hessian * dxdtheta_rest + dfdtheta_rest;
	if (r.norm() > 1e-7) {
		spdlog::warn(" residual too large when computing dxdtheta_rest {}", r.norm());
	}

	//we need to filter fixing dofs here
	/*for (int k = 0; k < dxdtheta_rest.outerSize(); ++k)
		for (SpMat::InnerIterator it(dxdtheta_rest, k); it; ++it)
		{
			int vertexIndex = it.row();
			if (m_constraints.find(vertexIndex) != m_constraints.end())
			{
				it.valueRef() = 0.0; //no matter how the thetas change, the displacement of fixing vertices will always be zero
			}
		}*/


}

void SOBendingModeShellStructure::step()
{
	return;
}

void SOBendingModeShellStructure::renderImGUIViewerWindow()
{
	ImGui::Begin("Mesh Window");
	{
		if (ImGui::Button("Reset simulation"))
		{
			reset();
		}
		
		ImGui::Checkbox("Draw wire frame", &m_drawWireframe);
		ImGui::Checkbox("Draw fixing vertices", &m_drawFixingVertices);
		ImGui::Checkbox("Use static solver", &m_useStaticSolver);
		ImGui::Checkbox("Draw both wire frame and mesh", &m_useBothWireframeAndmesh);
		ImGui::Checkbox("Draw max strain color", &m_isRenderStrainColor);
		ImGui::Checkbox("Draw hard edges", &m_drawHardEdges);
		
		ImGui::Separator();

		ImGui::InputInt("Animate Mode Index", &m_modesAnimationID);
		ImGui::InputDouble("Mode Amplitude", &m_modesAmplitude);


		ImGui::Separator();
		ImGui::InputDouble("Compute Mode Length", &m_computeModeLength);
		ImGui::InputDouble("Compute Mode step length", &m_computeModeStepLength);

		ImGui::Separator();
		if (ImGui::CollapsingHeader("Strain Space Modes"))
		{
			ImGui::Indent(45.0);
			/*if (ImGui::Button("Compute All Modes"))
			{
				computeBatchStrainSpaceMode();
			}*/

			ImGui::Checkbox("Animate Modes", &m_animateModes);
			ImGui::Checkbox("Use Linear Modes", &m_useLinearModes);
			ImGui::Checkbox("Use Strain Space Modes", &m_useStrainSpaceModes);
			ImGui::Checkbox("Use Developable Strain Space Modes", &m_useDevelopableStrainSpaceModes);
			/*
			if (m_animateModes)
			{
				if (m_useLinearModes && m_useStrainSpaceModes)
				{
					m_useLinearModes = true;
					m_useStrainSpaceModes = false;
				}
			}
			else
			{
				m_useLinearModes = false;
				m_useStrainSpaceModes = false;
			}*/

			if (ImGui::Button("Compute Modes State with Length"))
			{
				computeStrainSpaceModes();
			}
			ImGui::Separator();


			static int fromModeIndex = 6;
			ImGui::InputInt("From Mode Index", &fromModeIndex);
			if (ImGui::Button("Compute A Range of Modes"))
			{
				//if (fromModeIndex < 6)
				//	fromModeIndex = 6;
				int toModeIndex = m_modesAnimationID;
				for (size_t i = fromModeIndex; i <= toModeIndex; i++)
				{
					m_modesAnimationID = i;
					spdlog::info("Compute Mode {}", i);
					computeStrainSpaceModes();
					ouputStrainSpaceModesAtIndex(i);
				}
				//outputAllStrainSpaceModes();
			}
			if (ImGui::Button("Output All Modes"))
			{
				outputAllStrainSpaceModes();
			}
			static char strainSpaceModesFolder[1000] = "output/";
			ImGui::InputText("Strain Space Modes Data Folder", strainSpaceModesFolder, IM_ARRAYSIZE(strainSpaceModesFolder));

			if (ImGui::Button("Load Modes"))
			{
				loadAllStrainSpaceModes(strainSpaceModesFolder);
			}
			ImGui::Checkbox("Update modes states", &m_updateModesToVisualization);
			if (m_strainSpaceModes_index_to_lengthState.find(m_modesAnimationID) != m_strainSpaceModes_index_to_lengthState.end()
				&& m_updateModesToVisualization)
			{
				static int animateStrainSpaceModeIndex = 0;

				ImGui::InputInt("Animte Nonlinear Mode Index", &animateStrainSpaceModeIndex);//,0, 100
				//if (animateStrainSpaceModeIndex >= m_strainSpaceModes_index_to_lengthState[m_modesAnimationID].size())
				//	animateStrainSpaceModeIndex = m_strainSpaceModes_index_to_lengthState[m_modesAnimationID].size() - 1;
				//else if (animateStrainSpaceModeIndex <= -2)
				//	animateStrainSpaceModeIndex = -1;

				if(animateStrainSpaceModeIndex == -1)
					setVisualizationConfiguration(m_para_X);
				else if(animateStrainSpaceModeIndex>=0 && animateStrainSpaceModeIndex <= m_strainSpaceModes_index_to_lengthState[m_modesAnimationID].size() - 1)
				{
					ImGui::Text("The length of this mode is %f", m_strainSpaceModes_index_to_lengthState[m_modesAnimationID].at(animateStrainSpaceModeIndex).first);
					setVisualizationConfiguration(m_strainSpaceModes_index_to_lengthState[m_modesAnimationID].at(animateStrainSpaceModeIndex).second);
				}
			}
			ImGui::Unindent(45.0);
		}
		ImGui::Separator();

		//directly change this two parameters will not influence bending stiffness
		ImGui::InputDouble("Young's modulus", &m_material->k1);
		ImGui::InputDouble("Poisson's ratio", &m_material->k2);
		ImGui::InputDouble("Bending stiffness", &m_material->kB);
	}
	ImGui::Separator();

	/*if (ImGui::CollapsingHeader("Visualize strain space"))
	{
		ImGui::Indent(16.0f);

		//load tests
		static char simulationMeshFolder[1000] = "output/simulationTriangleMesh";
		ImGui::InputText("Simulation mesh data folder", simulationMeshFolder, IM_ARRAYSIZE(simulationMeshFolder));

		if (ImGui::Button("Load all mesh data"))
		{
			bool loadSuccess1 = false;
			bool loadSuccess2 = false;
			bool loadSuccess3 = false;
			loadSuccess1 = loadFolderMeshesWithKeyWords(simulationMeshFolder, "outputMesh", m_simulationMeshes);
		}

		ImGui::Separator();
		if (m_simulationMeshes.size() > 0)
		{
			//ImGui::SetNextItemWidth(250);
			ImGui::SliderInt("Choose visualization step index", &m_chooseStrainSpaceMeshIndex, 0, m_strainSpaceMeshes.size() - 1);
		}
		bool old_animateStrainSpaceMesh = m_animateStrainSpaceMesh;
		ImGui::Checkbox("Animate strain space mesh", &m_animateStrainSpaceMesh);
		if (!old_animateStrainSpaceMesh && m_animateStrainSpaceMesh)
		{
			m_chooseStrainSpaceMeshIndex = 0;
			m_animateStartClock = clock();
		}
		ImGui::Unindent(16.0f);
	}*/

	ImGui::End();

	//for video
	//ImGui::Begin("Iteractive selection");
	//{
	//	ImGui::Checkbox("Enable interactive selection", &m_enableInteractiveSelection);
	//	if (ImGui::Button("Delete last added vertex"))
	//	{
	//		if (m_selectedPairs.size() > 0)
	//		{
	//			if (m_selectedPairs.back().second == -1)
	//			{
	//				m_selectedPairs.erase(m_selectedPairs.begin() + m_selectedPairs.size() - 1);
	//				m_selectedPairsColor.erase(m_selectedPairsColor.begin() + m_selectedPairsColor.size() - 1);
	//			}
	//			else
	//			{
	//				m_selectedPairs.back().second = -1;
	//			}
	//		}
	//	}

	//	static Eigen::VectorXd staticSolverState;
	//	static bool ready_staticSolverState = false;
	//	static bool show_staticSolveState = false;
	//	if (ImGui::Button("Compute attached spring static solver"))
	//	{
	//		m_updateModesToVisualization = false;

	//		double stiffness = 1e3;

	//		//we first add spring elements
	//		for (int i = 0; i < m_selectedPairs.size(); i++)
	//		{
	//			int v0 = m_selectedPairs[i].first;
	//			int v1 = m_selectedPairs[i].second;
	//			if (v0 == -1 || v1 == -1)
	//				continue;


	//			m_pairAttachmentElements.push_back(pairAttachmentElement());
	//			pairAttachmentElement& elem = m_pairAttachmentElements.back();

	//			std::vector<int> inds = { v0,v1 };

	//			std::vector<int> dofIndices;
	//			for (int dof_i = 0; dof_i < inds.size(); dof_i++)
	//			{
	//				for (int j = 0; j < dimension; j++)
	//					dofIndices.push_back(inds[dof_i] * dimension + j);
	//			}

	//			elem.init(dofIndices, dofIndices, stiffness, true);
	//			elem.setElementProperty(ExternalElement);
	//		}

	//		for (int i = 0; i < m_pairAttachmentElements.size(); i++)
	//		{
	//			m_elements.push_back(&m_pairAttachmentElements[i]);
	//		}

	//		for (int i = 0; i < m_hingeElements.size(); i++)
	//		{
	//			m_hingeElements[i].setRestAngle(m_hingeElements[i].getInitialRestAngle());
	//		}
	//		//then solve static problem
	//		m_x = visualize_m_x;
	//		projectToManifold();
	//		setVisualizationConfiguration(m_x);

	//		//we remove all spring elements
	//		for (int i = m_pairAttachmentElements.size() - 1; i >= 0; i--)
	//		{
	//			if (m_elements[m_elements.size() - 1] != &m_pairAttachmentElements[i])
	//				spdlog::error("the element order has a problem");

	//			m_elements.erase(m_elements.begin() + m_elements.size() - 1);
	//		}


	//		//save obj
	//		{
	//			string outputMeshName = m_outputPath + "/outputMesh_interactive_staticSolve.obj";
	//			spdlog::info("Output simulation triangle mesh with name: {}", outputMeshName);

	//			m_triMesh->vertices() = m_x;
	//			m_triMesh->writeFile_tinyobj(outputMeshName);
	//		}
	//		staticSolverState = m_x;
	//		ready_staticSolverState = true;
	//		show_staticSolveState = true;
	//	}


	//	static int strainSpaceModeIndex = 0;
	//	if (m_strainSpaceModes_index_to_lengthState.find(m_modesAnimationID) != m_strainSpaceModes_index_to_lengthState.end())
	//	{


	//		ImGui::InputInt("Mode State Index", &strainSpaceModeIndex);//,0, 100
	//		if (strainSpaceModeIndex >= m_strainSpaceModes_index_to_lengthState[m_modesAnimationID].size())
	//			strainSpaceModeIndex = m_strainSpaceModes_index_to_lengthState[m_modesAnimationID].size() - 1;
	//		else if (strainSpaceModeIndex <= 0)
	//			strainSpaceModeIndex = 0;

	//		//if (animateStrainSpaceModeIndex == -1)
	//		//	setVisualizationConfiguration(m_para_X);
	//		if (strainSpaceModeIndex >= 0 && strainSpaceModeIndex <= m_strainSpaceModes_index_to_lengthState[m_modesAnimationID].size() - 1 && m_updateModesToVisualization)
	//		{
	//			//ImGui::Text("The length of this mode is %f", m_strainSpaceModes_index_to_lengthState[m_modesAnimationID].at(animateStrainSpaceModeIndex).first);
	//			setVisualizationConfiguration(m_strainSpaceModes_index_to_lengthState[m_modesAnimationID].at(strainSpaceModeIndex).second);
	//		}

	//		if (ready_staticSolverState)
	//		{
	//			ImGui::Checkbox("Show Static state", &show_staticSolveState);
	//			if(show_staticSolveState)
	//				setVisualizationConfiguration(staticSolverState);
	//			else
	//			{
	//				if (strainSpaceModeIndex >= 0 && strainSpaceModeIndex <= m_strainSpaceModes_index_to_lengthState[m_modesAnimationID].size() - 1)
	//				{
	//					//ImGui::Text("The length of this mode is %f", m_strainSpaceModes_index_to_lengthState[m_modesAnimationID].at(animateStrainSpaceModeIndex).first);
	//					setVisualizationConfiguration(m_strainSpaceModes_index_to_lengthState[m_modesAnimationID].at(strainSpaceModeIndex).second);
	//				}
	//			}
	//		}
	//	}
	//}
	//ImGui::End();

	ImGui::Begin("Iteractive Pinning-debug");
	{
		ImGui::Checkbox("Show contact pairs", &m_drawInteractivePairs);
		static double contactPairtThresold = 0.01;
		ImGui::InputDouble("Contact pair threshold", &contactPairtThresold);
		if (ImGui::Button("Compute current contact pairs"))
		{
			m_pinningPairConstraint_set.clear();
			Eigen::MatrixXd V = (Eigen::Map<const Eigen::MatrixXd>(visualize_m_x.data(), dimension, visualize_m_x.size() / dimension)).transpose();
			m_pinningPairConstraint_set.build(m_ipc_collisionMesh, V, m_ipc_dhat);

			//constraints filter


			m_pinningPairConstraintColor.resize(m_pinningPairConstraint_set.size());
			for (int i = 0; i < m_pinningPairConstraint_set.size(); i++)
			{
				m_pinningPairConstraintColor[i] = Eigen::Vector4d(((double)rand() / (RAND_MAX)), ((double)rand() / (RAND_MAX)), ((double)rand() / (RAND_MAX)), 1.0);
			}
		}
		ImGui::Separator();

		ImGui::Checkbox("Enable interactive selection", &m_enableInteractiveSelection);
		if (ImGui::Button("Delete last added vertex"))
		{
			if (m_selectedPairs.size() > 0)
			{
				if (m_selectedPairs.back().second == -1)
				{
					m_selectedPairs.erase(m_selectedPairs.begin() + m_selectedPairs.size() - 1);
					m_selectedPairsColor.erase(m_selectedPairsColor.begin() + m_selectedPairsColor.size() - 1);
				}
				else
				{
					m_selectedPairs.back().second = -1;
				}
			}
		}

		if (ImGui::Button("Compute attached spring static solver"))
		{
			m_updateModesToVisualization = false;

			double stiffness = 1e3;

			//we first add spring elements
			for (int i = 0; i < m_selectedPairs.size(); i++)
			{
				int v0 = m_selectedPairs[i].first;
				int v1 = m_selectedPairs[i].second;
				if (v0 == -1 || v1 == -1)
					continue;


				m_pairAttachmentElements.push_back(pairAttachmentElement());
				pairAttachmentElement& elem = m_pairAttachmentElements.back();

				std::vector<int> inds = { v0,v1 };

				std::vector<int> dofIndices;
				for (int dof_i = 0; dof_i < inds.size(); dof_i++)
				{
					for (int j = 0; j < dimension; j++)
						dofIndices.push_back(inds[dof_i] * dimension + j);
				}

				elem.init(dofIndices, dofIndices, stiffness, true);
				elem.setElementProperty(ExternalElement);
			}

			for (int i = 0; i < m_pairAttachmentElements.size(); i++)
			{
				m_elements.push_back(&m_pairAttachmentElements[i]);
			}

			for (int i = 0; i < m_hingeElements.size(); i++)
			{
				m_hingeElements[i].setRestAngle(m_hingeElements[i].getInitialRestAngle());
			}
			//then solve static problem
			m_x = visualize_m_x;
			projectToManifold();
			setVisualizationConfiguration(m_x);

			//we remove all spring elements
			for (int i = m_pairAttachmentElements.size() - 1; i >= 0; i--)
			{
				if (m_elements[m_elements.size() - 1] != &m_pairAttachmentElements[i])
					spdlog::error("the element order has a problem");

				m_elements.erase(m_elements.begin() + m_elements.size() - 1);
			}


			//save obj
			{
				string outputMeshName = m_outputPath + "/outputMesh_interactive_staticSolve.obj";
				spdlog::info("Output simulation triangle mesh with name: {}", outputMeshName);

				m_triMesh->vertices() = m_x;
				m_triMesh->writeFile_tinyobj(outputMeshName);
			}
		}

	}
	ImGui::End();


	ImGui::Begin("Blend strain space modes");
	{
		static std::vector<int> selectedModes;
		static std::vector<float> selectModeLength;
		static int addModeIndex = 6;
		ImGui::InputInt("Add mode index", &addModeIndex);
		if (ImGui::Button("Add mode"))
		{
			if (std::find(selectedModes.begin(), selectedModes.end(), addModeIndex) == selectedModes.end())
			{
				selectedModes.push_back(addModeIndex);
				selectModeLength.push_back(0.1);
			}
		}
		if (ImGui::Button("Delete last one"))
		{
			if(selectedModes.size()>0)
				selectedModes.pop_back();
			if (selectModeLength.size() > 0)
				selectModeLength.pop_back();
		}
		static double blendLengthStepSize = 0.1;
		//ImGui::InputDouble("Length step size", &blendLengthStepSize);
		for (int i=0; i<selectedModes.size(); i++)
		{
			string name = "mode " + to_string(selectedModes[i]) + " length";
			ImGui::SliderFloat(name.c_str(), &selectModeLength[i], 0, 1.0, "%5f");
		}

		ImGui::Separator();
		static int numMaxFastExplorationIterations = 1;
		ImGui::InputInt("Num max iterations", &numMaxFastExplorationIterations);
		static bool useOnlyTranslation = true;
		ImGui::Checkbox("Use only translation", &useOnlyTranslation);

		ImGui::Separator();

		ImGui::Checkbox("Update selected modes states", &m_updateModesToVisualization);
		if (m_selectedStrainSpce_to_eachModeLengthState.find(selectedModes) != m_selectedStrainSpce_to_eachModeLengthState.end()
			&& m_updateModesToVisualization)
		{
			static int animateStrainSpaceSelectedModeIndex = 0;

			ImGui::InputInt("Animate Selected Mode Index", &animateStrainSpaceSelectedModeIndex);//,0, 100
			
			if (animateStrainSpaceSelectedModeIndex == -1)
				setVisualizationConfiguration(m_para_X);
			else if (animateStrainSpaceSelectedModeIndex >= 0 && animateStrainSpaceSelectedModeIndex <= m_selectedStrainSpce_to_eachModeLengthState[selectedModes].second.size() - 1)
			{
				auto selectedModeAnimateIter = m_selectedStrainSpce_to_eachModeLengthState[selectedModes].second.begin();
				for (int i = 0; i < animateStrainSpaceSelectedModeIndex; i++)
					selectedModeAnimateIter++;
				Eigen::VectorXi sampleIndices = Eigen::Map<const Eigen::VectorXi>(selectedModeAnimateIter->first.data(), selectedModeAnimateIter->first.size());
				double stepSize = m_selectedStrainSpce_to_eachModeLengthState[selectedModes].first;
				Eigen::VectorXd modeLength(sampleIndices.size());
				string outputText = "The length of this mode is ";
				for (int i = 0; i < sampleIndices.size(); i++)
				{
					modeLength[i] = sampleIndices[i] * stepSize;
					outputText += (to_string(modeLength[i]) + " ");
				}

				ImGui::Text(outputText.c_str());
				setVisualizationConfiguration(selectedModeAnimateIter->second);
			}
			
		}



		ImGui::Separator();
		//sensitivity analysis
		if (ImGui::Button("Compute sensitivity matrix"))
		{
			/*Eigen::VectorXd theta(m_hingeElements.size());
			theta.setZero();
			theta[0] = 1e-3;
			for (int i = 0; i < m_hingeElements.size(); i++)
			{
				m_hingeElements[i].setRestAngle(theta[i]);
			}*/
			projectToManifold();
			m_sensitivityAnalysisProjectedState = m_x;

			computeStrainSpaceSensitiveMatrix(m_sensitivityAnalysisProjectedState, m_para_X, m_dxdtheta_rest);
			m_sensitivityMatrixUpdated = true;
			m_sensitivityExplorationProjected = false;
		}
		ImGui::Checkbox("Visualize Sensitivity Exploration", &m_visualizeSensitivityExploration);
		if (m_visualizeSensitivityExploration)
		{
			computeBatchStrainSpaceMode();
			Eigen::VectorXd sensitivityExplorationState(m_strainSpaceModes[0].size());
			if (m_sensitivityMatrixUpdated && m_dxdtheta_rest.size() > 0 && selectedModes.size() > 0 && selectModeLength.size() > 0 && m_strainSpaceModes.size()>0)
			{
				sensitivityExplorationState.setZero();
				for (int i = 0; i < selectedModes.size(); i++)
				{
					int modeIndex = selectedModes[i];
					double modeLength = selectModeLength[i];
					sensitivityExplorationState += m_strainSpaceModes[modeIndex] * modeLength;
				}
				//m_dxdtheta_rest
				setVisualizationConfiguration(m_sensitivityAnalysisState + m_dxdtheta_rest * sensitivityExplorationState);

			}

			if (ImGui::Button("Reproject"))
			{
				if (m_sensitivityMatrixUpdated && !m_sensitivityExplorationProjected)
				{
					//bool temp = USE_IPC;
					//USE_IPC = false;
					//update initial guess
					m_x = m_sensitivityAnalysisState + m_dxdtheta_rest * sensitivityExplorationState;
					Eigen::MatrixXd V = (Eigen::Map<const Eigen::MatrixXd>(m_x.data(), dimension, m_x.size() / dimension)).transpose();
					if(ipc::has_intersections(m_ipc_collisionMesh, V))
						m_x =  m_sensitivityAnalysisState;//if it has intersection, we will use the last state
					tbb::parallel_for(
						tbb::blocked_range<size_t>(size_t(0), m_hingeElements.size()),
						[&](const tbb::blocked_range<size_t>& r) {
							for (size_t i = r.begin(); i < r.end(); i++)
							{
								//set to anhle
								double theta_i = m_hingeElements[i].getCurrentAngle(m_sensitivityAnalysisProjectedState) +
									sensitivityExplorationState[i];
								m_hingeElements[i].setRestAngle(theta_i);
							}
						});
					/*for (int i = 0; i < m_hingeElements.size(); i++)
					{
						std::cout << "hinge " << i << " restAngle " << m_hingeElements[i].getRestAngle() << std::endl;
					}*/
					std::cout << endl;

					std::thread t1([this] {
						spdlog::info("Projecting...");
						m_sensitivityMatrixUpdated = false;
						projectToManifold();
						m_sensitivityAnalysisProjectedState = m_x;
						visualize_m_x = m_x;
						m_sensitivityExplorationProjected = true;
						//USE_IPC = temp;

						});
					t1.detach();
					
				}
			}
		}
		

		ImGui::Separator();

		auto computeAllStrainSpaceModes = [&]()
		{
			int startModeIndex = 0;
			int numModes = m_strainSpaceModes.size() - startModeIndex;
			Eigen::MatrixXd allStrainSpaceModes(m_hingeElements.size(), numModes);
			for (int i = startModeIndex; i < m_strainSpaceModes.size(); i++)
			{
				allStrainSpaceModes.col(i - startModeIndex) = m_strainSpaceModes[i].normalized();
			}
			return allStrainSpaceModes;
		};
		auto computeStatesAlongParameters = [&](int stepsAlongParameters, const Eigen::VectorXd& modesCoefficients)
		{
			m_statesAlongModesCoefficients.clear();
			Eigen::MatrixXd allStrainSpaceModes = computeAllStrainSpaceModes();

			Eigen::VectorXd last_state = m_para_X;
			Eigen::VectorXd mode_0(m_hingeElements.size());
			for (int i = 0; i < m_hingeElements.size(); i++)
				mode_0[i] = m_hingeElements[i].getInitialRestAngle();

			m_statesAlongModesCoefficients.push_back(m_para_X);
			for (int i = 1; i <= stepsAlongParameters; i++)
			{
				Eigen::VectorXd step_i_parameters = modesCoefficients / double(stepsAlongParameters) * double(i);

				Eigen::VectorXd mode_p = mode_0 + allStrainSpaceModes * step_i_parameters.head(allStrainSpaceModes.cols());

				for (int i = 0; i < m_hingeElements.size(); i++)
				{
					m_hingeElements[i].setRestAngle(mode_p[i]);
				}
				m_x = last_state;
				projectToManifold();

				last_state = m_x;

				if (modesCoefficients.size() == allStrainSpaceModes.cols() + 6)
				{
					//apply rigid transformation
					Eigen::Vector3d translation = step_i_parameters.tail(3);
					Eigen::Vector3d rotation_vector = step_i_parameters.tail(6).head(3);
					Eigen::Matrix3d rotationMatrix = soutil::construct_rotation_matrix(rotation_vector);

					const Eigen::MatrixXd& m_x_map = Eigen::Map<const Eigen::MatrixXd>(m_x.data(), 3, m_x.size() / 3);

					Eigen::MatrixXd transformed_m_x = rotationMatrix * m_x_map;
					transformed_m_x.colwise() += translation;
					m_statesAlongModesCoefficients.push_back(Eigen::Map<Eigen::VectorXd>(transformed_m_x.data(), transformed_m_x.size()));
				}
				else
					m_statesAlongModesCoefficients.push_back(m_x);
			}
		};

		bool before = m_visualizeTargetFoldingShape;

		ImGui::Checkbox("Show target folding shape", &m_visualizeTargetFoldingShape);
		static Eigen::VectorXd before_changed;
		if (before != m_visualizeTargetFoldingShape && m_visualizeTargetFoldingShape)
		{
			before_changed = visualize_m_x;
		}
		else if(before_changed.size()>0 && before != m_visualizeTargetFoldingShape && !m_visualizeTargetFoldingShape)
		{
			visualize_m_x = before_changed;
		}


		static bool computedFoldingParameters = false;
		ImGui::InputDouble("Target sphere radius", &m_targetSphereRadius);
		ImGui::Checkbox("Show target sphere", &m_showTargetSphere);
		if (ImGui::Button("Compute States in Sphere"))
		{
			static bool initial = false;
			if (!initial)
			{
				m_statesAlongModesCoefficients.clear();
				m_statesAlongModesCoefficients.push_back(m_para_X);//rest state is the first state
				initial = true;
			}
			//optimizePositionsToTargetSphere(m_para_X, m_targetSphereCenter, m_targetSphereRadius, m_statesAlongModesCoefficients);
			double radius = 0.1;
			while (radius > m_targetSphereRadius)
			{
				radius -= 1e-3;
				spdlog::info("Target sphere radius is {}", radius);
				optimizePositionsToTargetIPCSphere(m_para_X, m_targetSphereCenter, radius, m_statesAlongModesCoefficients);
			}
			computedFoldingParameters = true;
		}
		if (ImGui::Button("Compute Folding Parameters"))
		{
			if (m_targetTriMesh != nullptr)
			{
				computeBatchStrainSpaceMode();

				Eigen::MatrixXd allStrainSpaceModes = computeAllStrainSpaceModes();

				int numModes = allStrainSpaceModes.cols();
				
				if (m_modesCoefficients.size() != numModes + 6)
				{
					m_modesCoefficients = Eigen::VectorXd(numModes + 6);
					m_modesCoefficients.setZero();
					m_modesCoefficients[6] = 0.01;
					m_statesAlongModesCoefficients.clear();
					m_statesAlongModesCoefficients.push_back(m_para_X);//rest state is the first state
				}
				//Eigen::VectorXd targetShape_x = m_targetTriMesh->vertices();
				
				m_modesCoefficients = optimizeOneTimeFoldingParametersWithinSphereWithRigidTranformation(m_modesCoefficients, allStrainSpaceModes, m_targetSphereCenter, m_targetSphereRadius, m_statesAlongModesCoefficients);
			
				computedFoldingParameters = true;
			}
		}

		ImGui::Separator();

		if (computedFoldingParameters)
		{
			{
				static int showStateIndex = 0;
				ImGui::DragInt("Show state index", &showStateIndex, 1, 0, m_statesAlongModesCoefficients.size()-1);

				setVisualizationConfiguration(m_statesAlongModesCoefficients[showStateIndex]);
			}
		}


	}
	ImGui::End();

}

void SOBendingModeShellStructure::ouputNonlinearCompliantModesAtIndex(int index)
{
	if (m_nonlinearCompliantModes_index_to_lengthState.find(index) == m_nonlinearCompliantModes_index_to_lengthState.end())
		return;

	int modeIndex = index;

	const auto& modeStates = m_nonlinearCompliantModes_index_to_lengthState[index];
	//build a folder
	string subFolderName = "/nonlinearCompliantMode_" + to_string(modeIndex);
	struct stat st_2 = { 0 };
	if (stat((m_outputPath + subFolderName).c_str(), &st_2) == -1) {
#ifdef _WIN32
            mkdir((m_outputPath + subFolderName).c_str());
#elif __APPLE__ || __linux__ 
            mkdir((m_outputPath + subFolderName).c_str(),0777);
#endif	
	}

	//we first output rest state
	{
		double modeLength = 0.0;
		const Eigen::VectorXd modeState = m_para_X;

		string outputMeshName = m_outputPath + subFolderName + "/outputMesh_Mode_" + to_string(modeIndex) + "_" + to_string_with_precision(modeLength, 5) + ".obj";
		spdlog::info("Output simulation triangle mesh with name: {}", outputMeshName);

		m_triMesh->vertices() = modeState;
		/*Eigen::MatrixXd outputTriangleMeshVertices;
		Eigen::MatrixXi outputTriangleMeshFaces;
		loadFromMesh(m_triMesh, outputTriangleMeshVertices, outputTriangleMeshFaces);
		outputTriangleMeshVertices.transposeInPlace();
		outputTriangleMeshFaces.transposeInPlace();

		igl::writeOBJ(outputMeshName, outputTriangleMeshVertices, outputTriangleMeshFaces);*/
		m_triMesh->writeFile_tinyobj(outputMeshName);
	}

	for (int i = 0; i < modeStates.size(); i++)
	{
		double modeLength = modeStates[i].first;
		const Eigen::VectorXd modeState = modeStates[i].second;

		string outputMeshName = m_outputPath + subFolderName + "/outputMesh_Mode_" + to_string(modeIndex) + "_" + to_string_with_precision(modeLength, 5) + ".obj";
		spdlog::info("Output simulation triangle mesh with name: {}", outputMeshName);

		m_triMesh->vertices() = modeState;
		/*Eigen::MatrixXd outputTriangleMeshVertices;
		Eigen::MatrixXi outputTriangleMeshFaces;
		loadFromMesh(m_triMesh, outputTriangleMeshVertices, outputTriangleMeshFaces);
		outputTriangleMeshVertices.transposeInPlace();
		outputTriangleMeshFaces.transposeInPlace();

		igl::writeOBJ(outputMeshName, outputTriangleMeshVertices, outputTriangleMeshFaces);*/
		m_triMesh->writeFile_tinyobj(outputMeshName);

	}
}

void SOBendingModeShellStructure::ouputRotationStrainSpaceModesAtIndex(int index)
{
	if (m_rotationStrainSpaceModes_index_to_lengthState.find(index) == m_rotationStrainSpaceModes_index_to_lengthState.end())
		return;

	int modeIndex = index;

	const auto& modeStates = m_rotationStrainSpaceModes_index_to_lengthState[index];
	//build a folder
	string subFolderName = "/rotationStrainSpaceMode_" + to_string(modeIndex);
	struct stat st_2 = { 0 };
	if (stat((m_outputPath + subFolderName).c_str(), &st_2) == -1) {
#ifdef _WIN32
            mkdir((m_outputPath + subFolderName).c_str());
#elif __APPLE__ || __linux__ 
            mkdir((m_outputPath + subFolderName).c_str(),0777);
#endif	
	}

	//we first output rest state
	{
		double modeLength = 0.0;
		const Eigen::VectorXd modeState = m_para_X;

		string outputMeshName = m_outputPath + subFolderName + "/outputMesh_Mode_" + to_string(modeIndex) + "_" + to_string_with_precision(modeLength, 5) + ".obj";
		spdlog::info("Output simulation triangle mesh with name: {}", outputMeshName);

		m_triMesh->vertices() = modeState;
		/*Eigen::MatrixXd outputTriangleMeshVertices;
		Eigen::MatrixXi outputTriangleMeshFaces;
		loadFromMesh(m_triMesh, outputTriangleMeshVertices, outputTriangleMeshFaces);
		outputTriangleMeshVertices.transposeInPlace();
		outputTriangleMeshFaces.transposeInPlace();

		igl::writeOBJ(outputMeshName, outputTriangleMeshVertices, outputTriangleMeshFaces);*/
		m_triMesh->writeFile_tinyobj(outputMeshName);
	}

	for (int i = 0; i < modeStates.size(); i++)
	{
		double modeLength = modeStates[i].first;
		const Eigen::VectorXd modeState = modeStates[i].second;

		string outputMeshName = m_outputPath + subFolderName + "/outputMesh_Mode_" + to_string(modeIndex) + "_" + to_string_with_precision(modeLength, 5) + ".obj";
		spdlog::info("Output simulation triangle mesh with name: {}", outputMeshName);

		m_triMesh->vertices() = modeState;
		/*Eigen::MatrixXd outputTriangleMeshVertices;
		Eigen::MatrixXi outputTriangleMeshFaces;
		loadFromMesh(m_triMesh, outputTriangleMeshVertices, outputTriangleMeshFaces);
		outputTriangleMeshVertices.transposeInPlace();
		outputTriangleMeshFaces.transposeInPlace();

		igl::writeOBJ(outputMeshName, outputTriangleMeshVertices, outputTriangleMeshFaces);*/
		m_triMesh->writeFile_tinyobj(outputMeshName);

	}
}


void SOBendingModeShellStructure::ouputStrainSpaceModesAtIndex(int index)
{
	if (m_strainSpaceModes_index_to_lengthState.find(index) == m_strainSpaceModes_index_to_lengthState.end())
		return;

	int modeIndex = index;

	const auto& modeStates = m_strainSpaceModes_index_to_lengthState[index];
	//build a folder
	string subFolderName = "/Mode_" + to_string(modeIndex);
	struct stat st_2 = { 0 };
	if (stat((m_outputPath + subFolderName).c_str(), &st_2) == -1) {
#ifdef _WIN32
            mkdir((m_outputPath + subFolderName).c_str());
#elif __APPLE__ || __linux__ 
            mkdir((m_outputPath + subFolderName).c_str(),0777);
#endif	
	}

	//we first output rest state
	{
		double modeLength = 0.0;
		const Eigen::VectorXd modeState = m_para_X;

		string outputMeshName = m_outputPath + subFolderName + "/outputMesh_Mode_" + to_string(modeIndex) + "_" + to_string_with_precision(modeLength, 5) + ".obj";
		spdlog::info("Output simulation triangle mesh with name: {}", outputMeshName);

		m_triMesh->vertices() = modeState;
		/*Eigen::MatrixXd outputTriangleMeshVertices;
		Eigen::MatrixXi outputTriangleMeshFaces;
		loadFromMesh(m_triMesh, outputTriangleMeshVertices, outputTriangleMeshFaces);
		outputTriangleMeshVertices.transposeInPlace();
		outputTriangleMeshFaces.transposeInPlace();

		igl::writeOBJ(outputMeshName, outputTriangleMeshVertices, outputTriangleMeshFaces);*/
		m_triMesh->writeFile_tinyobj(outputMeshName);
	}

	for (int i = 0; i < modeStates.size(); i++)
	{
		double modeLength = modeStates[i].first;
		const Eigen::VectorXd modeState = modeStates[i].second;

		string outputMeshName = m_outputPath + subFolderName + "/outputMesh_Mode_" + to_string(modeIndex) + "_" + to_string_with_precision(modeLength, 5) + ".obj";
		spdlog::info("Output simulation triangle mesh with name: {}", outputMeshName);

		m_triMesh->vertices() = modeState;
		/*Eigen::MatrixXd outputTriangleMeshVertices;
		Eigen::MatrixXi outputTriangleMeshFaces;
		loadFromMesh(m_triMesh, outputTriangleMeshVertices, outputTriangleMeshFaces);
		outputTriangleMeshVertices.transposeInPlace();
		outputTriangleMeshFaces.transposeInPlace();

		igl::writeOBJ(outputMeshName, outputTriangleMeshVertices, outputTriangleMeshFaces);*/
		m_triMesh->writeFile_tinyobj(outputMeshName);

	}
}

void SOBendingModeShellStructure::outputAllStrainSpaceModes()
{
	for (const auto& item : m_strainSpaceModes_index_to_lengthState)
	{

		int modeIndex = item.first;
		ouputStrainSpaceModesAtIndex(modeIndex);
	}
}


bool SOBendingModeShellStructure::loadStrainSpaceMeshFolder(const std::string& testSetFolder, std::vector<std::pair<std::string, double>>& fileName_length)
{
	struct stat info;
	if (stat((testSetFolder).c_str(), &info) != 0)
	{
		spdlog::info("Directory not exists {:s}\n", testSetFolder.c_str());
		return false;
	}

	bool loadSuccess = true;
	for (const auto& entry : std::filesystem::directory_iterator(testSetFolder))
	{
		std::string entry_name = entry.path().string();
		//test if this is a directory or not
		struct stat info;
		if (stat(entry_name.c_str(), &info) != 0)
		{
			spdlog::info("cannot access {:s}", entry_name.c_str());
			return false;
		}
		else if (info.st_mode & S_IFDIR)
		{
			spdlog::info("{:s} is a directory", entry_name.c_str());
		}
		else
		{
			//const std::string preffix = entry_name.substr(0, entry_name.find_last_of('.'));
			//spdlog::info("{:s} is no directory, with preffix {:s}", entry_name.c_str(), preffix.c_str());//this is a file
			//fileName_length.emplace_back(entry_name, 0.0);
			const std::string preffix = entry_name.substr(0, entry_name.find_last_of('.'));
			const std::string ext = entry_name.substr(entry_name.find_last_of('.')+1, entry_name.size() - 1);
			if (ext == "obj")
			{
				spdlog::info("{:s} is no directory, with preffix {:s} and ext {:s}", entry_name.c_str(), preffix.c_str(), ext.c_str());//this is a file
				fileName_length.emplace_back(entry_name, 0.0);
			}
		}
	}

	auto myCmp = [](std::pair<string, double>& s1, std::pair<string,double>& s2)
	{
		std::size_t const s1_length_startIndex = s1.first.find_last_of("_") + 1;
		std::size_t const s1_timeStep_endIndex = s1.first.find_last_of(".obj") - 1;
		double s1_length = std::stod(s1.first.substr(s1_length_startIndex, (s1_timeStep_endIndex - s1_length_startIndex + 1)));


		std::size_t const s2_length_startIndex = s2.first.find_last_of("_") + 1;
		std::size_t const s2_timeStep_endIndex = s2.first.find_last_of(".obj") - 1;
		double s2_length = std::stod(s2.first.substr(s2_length_startIndex, (s2_timeStep_endIndex - s2_length_startIndex + 1)));

		s1.second = s1_length;
		s2.second = s2_length;

		if (s1_length < s2_length)
			return true;
		else
			return false;
	};

	std::sort(fileName_length.begin(), fileName_length.end(), myCmp);
	return loadSuccess;
}


bool SOBendingModeShellStructure::loadFolderMeshes(const string& strainSpaceMeshFolder, std::vector<std::pair<TriMesh,double>>& meshSet_length)
{

	std::vector<std::pair<std::string,double>> fileName;
	bool loadSuccess = loadStrainSpaceMeshFolder(strainSpaceMeshFolder, fileName);
	if (loadSuccess)
	{
		spdlog::info("Load test file and folder success!");

		meshSet_length.clear();
		meshSet_length.resize(fileName.size());
		//then, we can load .txt file and folder .obj parallelly
#ifdef USE_TBB
		tbb::parallel_for(0, (int)fileName.size(), 1, [&](int i)
#else
		for (int i = 0; i < fileName.size(); i++)
#endif
		{
			std::string strainSpaceMeshFileName = fileName[i].first;

			//we according each row in file to determine the name of obj in the folder
			std::string errorMessage;
			meshSet_length[i].first.loadFromFile_tinyobj(strainSpaceMeshFileName.c_str(), errorMessage);
			meshSet_length[i].second = fileName[i].second;
		}
#ifdef USE_TBB
		);
#endif
	}
	return loadSuccess;
}

bool SOBendingModeShellStructure::loadAllNonlinearCompliantModes(const string& strainSpaceModesFolder)
{
	struct stat info;
	if (stat((strainSpaceModesFolder).c_str(), &info) != 0)
	{
		spdlog::info("Directory not exists {:s}\n", strainSpaceModesFolder.c_str());
		return false;
	}

	bool loadSuccess = true;
	std::vector<std::string> folderName;
	for (const auto& entry : std::filesystem::directory_iterator(strainSpaceModesFolder))
	{
		std::string entry_name = entry.path().string();
		//test if this is a directory or not
		struct stat info;
		if (stat(entry_name.c_str(), &info) != 0)
		{
			spdlog::info("cannot access {:s}", entry_name.c_str());
			return false;
		}
		else if (info.st_mode & S_IFDIR)
		{
			spdlog::info("{:s} is a directory", entry_name.c_str());
			folderName.push_back(entry_name);
		}
		else
		{
			const std::string preffix = entry_name.substr(0, entry_name.find_last_of('.'));
			spdlog::info("{:s} is no directory, with preffix {:s}", entry_name.c_str(), preffix.c_str());//this is a file
		}
	}


	m_nonlinearCompliantModes_index_to_lengthState.clear();
	for (int i = 0; i < folderName.size(); i++)
	{
		std::size_t const mode_startIndex = folderName[i].find_last_of("_") + 1;
		std::size_t const mode_endIndex = folderName[i].size() - 1;

		int modeIndex = stoi(folderName[i].substr(mode_startIndex, mode_endIndex - mode_startIndex + 1));

		std::vector<std::pair<TriMesh, double>> meshSet_length;
		loadFolderMeshes(folderName[i], meshSet_length);

		for (int j = 0; j < meshSet_length.size(); j++)
		{
			m_nonlinearCompliantModes_index_to_lengthState[modeIndex].emplace_back(meshSet_length[j].second, meshSet_length[j].first.vertices());
			/*if (j == meshSet_length.size() - 1)
			{
				for (int e_i = 0; e_i < m_hingeElements.size(); e_i++)
				{
					double theta = m_hingeElements[e_i].getCurrentAngle(meshSet_length[j].first.vertices());
					std::cout << "hinge " << e_i << " restAngle " << theta << std::endl;
				}
			}*/
		}
	}
}

bool SOBendingModeShellStructure::loadAllStrainSpaceModes(const string& strainSpaceModesFolder)
{
	struct stat info;
	if (stat((strainSpaceModesFolder).c_str(), &info) != 0)
	{
		spdlog::info("Directory not exists {:s}\n", strainSpaceModesFolder.c_str());
		return false;
	}

	bool loadSuccess = true;
	std::vector<std::string> folderName;
	for (const auto& entry : std::filesystem::directory_iterator(strainSpaceModesFolder))
	{
		std::string entry_name = entry.path().string();
		//test if this is a directory or not
		struct stat info;
		if (stat(entry_name.c_str(), &info) != 0)
		{
			spdlog::info("cannot access {:s}", entry_name.c_str());
			return false;
		}
		else if (info.st_mode & S_IFDIR)
		{
			spdlog::info("{:s} is a directory", entry_name.c_str());
			folderName.push_back(entry_name);
		}
		else
		{
			const std::string preffix = entry_name.substr(0, entry_name.find_last_of('.'));
			spdlog::info("{:s} is no directory, with preffix {:s}", entry_name.c_str(), preffix.c_str());//this is a file
		}
	}


	m_strainSpaceModes_index_to_lengthState.clear();
	for (int i = 0; i < folderName.size(); i++)
	{
		std::size_t const mode_startIndex = folderName[i].find_last_of("_") + 1;
		std::size_t const mode_endIndex = folderName[i].size() - 1;

		int modeIndex = stoi(folderName[i].substr(mode_startIndex, mode_endIndex - mode_startIndex + 1));

		std::vector<std::pair<TriMesh,double>> meshSet_length;
		loadFolderMeshes(folderName[i], meshSet_length);

		for (int j = 0; j < meshSet_length.size(); j++)
		{
			m_strainSpaceModes_index_to_lengthState[modeIndex].emplace_back(meshSet_length[j].second, meshSet_length[j].first.vertices());
			
			/*for (int e_i = 0; e_i < m_hingeElements.size(); e_i++)
			{
				double theta = m_hingeElements[e_i].getCurrentAngle();
				std::cout << "hinge " << e_i << " restAngle " << theta << std::endl;
			}*/
			Eigen::VectorXd maxStrainVec(m_cstShellElements.size());
			for (int i = 0; i < m_cstShellElements.size(); i++)
			{
				double maxStrain = 0.0;
				m_cstShellElements[i].computeMaxStrain(meshSet_length[j].first.vertices(), m_para_X, maxStrain);
				maxStrainVec[i] = maxStrain;
			}

			spdlog::info("Mode {} length {} max Strain {}", modeIndex, meshSet_length[j].second, maxStrainVec.maxCoeff());
		}
	}
}

void SOBendingModeShellStructure::computeMaxStrainFacesColorMap(const double colorMap_UpperBound, const Eigen::VectorXd& vx, const Eigen::VectorXd& vX, Eigen::MatrixXd& facesColor)
{
	facesColor = Eigen::MatrixXd(3, m_cstShellElements.size());

	Eigen::VectorXd maxStrainVec(m_cstShellElements.size());
	for (int i = 0; i < m_cstShellElements.size(); i++)
	{
		double maxStrain = 0.0;
		m_cstShellElements[i].computeMaxStrain(vx, vX, maxStrain);
		maxStrainVec[i] = maxStrain;
		//cout << "maxStrain " << maxStrain << endl;

		V3D rgb;
		colormap::colormapFromScalar(maxStrain, 0.0, colorMap_UpperBound, rgb[0], rgb[1], rgb[2], colormap::COOL_COLORMAP);

		//m_triMesh->setFaceColor(rgb, j);
		facesColor.col(i) = rgb;
	}

	cout << "max coeff of max strain " << maxStrainVec.maxCoeff() << endl;
}


void SOBendingModeShellStructure::render()
{
	//m_primitiveDrawer->drawSphere(0.05, Eigen::Vector3d(0, 0, 0), Eigen::Vector4d(0.7, 0.7, 0.7, 0.2));
	if (m_visualizeTargetFoldingShape)
	{
		visualize_m_x = m_targetTriMesh->vertices();
	}

	{
		if (m_isRenderStrainColor)
		{ 
			double maximumStrain = 0.023;// 1.37;
			Eigen::MatrixXd facesStrainColors;
			computeMaxStrainFacesColorMap(maximumStrain, m_triMesh->vertices(), m_para_X, facesStrainColors);

			//since our drawer is vertex-color based render, we need to convert to isolated triangle mesh
			Eigen::MatrixXd isolatedVertices;
			Eigen::MatrixXi isolatedFaces;
			Eigen::MatrixXd isolatedVerticesColors;
			m_triMesh->faceColors() = Eigen::Map<Eigen::VectorXd>(facesStrainColors.data(), facesStrainColors.size());
			m_triMesh->computeIsolatedMesh(isolatedVertices, isolatedFaces, isolatedVerticesColors);

			m_meshDrawer->setMesh(isolatedVertices, isolatedFaces);
			m_meshDrawer->draw(isolatedVerticesColors);
			//m_meshDrawer->drawWireframe(Eigen::Vector4d(0.0, 0.0, 0.0, 1.0));
		}
		else
		{
			//linear mode
			//Eigen::Vector4d color(69.0 / 255.0, 179.0 / 255.0, 255.0 / 255.0, 1.0);
			//nonlinear compliant mode
			//Eigen::Vector4d color(247.0 / 255.0, 171.0 / 255.0, 104.0 / 255.0, 1.0);
			//strain space mode
			//Eigen::Vector4d color(135.0 / 255.0, 250.0 / 255.0, 135.0 / 255.0, 1.0);
			//RS mode
			//Eigen::Vector4d color(126.0 / 255.0, 47.0 / 255.0, 142.0 / 255.0, 1.0);
			Eigen::Vector4d color(0.8, 0.8, 0.8, 1.0);
			if (visualize_m_x != m_triMesh->vertices())
			{
				setVisualizationConfiguration(visualize_m_x);
			}
			m_meshDrawer->setMesh(m_triMesh);
			if (m_useBothWireframeAndmesh)
			{
				m_meshDrawer->drawWireframe(Eigen::Vector4d(0.0, 0.0, 0.0, 1.0));
				m_meshDrawer->draw(color);
			}
			else if (m_drawWireframe)
			{
				m_meshDrawer->drawWireframe(Eigen::Vector4d(0.0, 0.0, 0.0, 1.0));
			}
			else if (m_drawMaxCurvature && m_vertexMaxCurvature.size() >0)
			{
				Eigen::MatrixXd maxCurvatureColors(3, m_vertexMaxCurvature.size());
				Eigen::VectorXd abs_vertexGaussianCurvature = m_vertexMaxCurvature;// .cwiseAbs();
				double maxCoeff = abs_vertexGaussianCurvature.maxCoeff();
				double minCoeff = abs_vertexGaussianCurvature.minCoeff();
				for (int v_i = 0; v_i < m_vertexMaxCurvature.size(); v_i++)
				{
					V3D rgb;
					colormap::colormapFromScalar(m_vertexMaxCurvature[v_i], minCoeff, maxCoeff, rgb[0], rgb[1], rgb[2], colormap::JET_COLORMAP);//colormap::COOL_COLORMAP
					maxCurvatureColors.col(v_i) = rgb;
				}

				m_meshDrawer->draw(maxCurvatureColors);
			}
			else
				m_meshDrawer->draw(color);

		}
	}

	if (m_drawHardEdges)
	{
		double radius = 0.001;
		for (int i = 0; i < m_hingeElements.size(); i++)
		{
			if (m_hingeElements[i].getHardEdge())
			{
				const auto hingeDofs = m_hingeElements[i].getDofIndices();
				int e0 = hingeDofs[0] / 3;
				int e1 = hingeDofs[3] / 3;
				Eigen::Vector3d v0 = visualize_m_x.segment<3>(e0 * 3);
				Eigen::Vector3d v1 = visualize_m_x.segment<3>(e1 * 3);
				m_primitiveDrawer->drawCylinder(radius, v0, v1, Eigen::Vector4d(0.7, 0.7, 0.7, 1.0));
			}
		}
	}

	if (m_drawFixingVertices)
	{
		//draw fixed vertices
		std::set<int> fixedVertexIndices;
		for (int i = 0; i < m_fixingDofs.size(); i++)
			fixedVertexIndices.insert(int(m_fixingDofs[i] / dimension));
		for (auto iter = fixedVertexIndices.begin(); iter != fixedVertexIndices.end(); iter++)
		{
			//double radius = 0.01;
			double radius = 0.001;
			m_primitiveDrawer->drawSphere(radius, visualize_m_x.segment<3>(3 * (*iter)), Eigen::Vector4d(0.8, 0.8, 0.8, 1.0));
		}

		for (int i = 0; i < m_fixVertexOnLinePairs.size(); i++)
		{
			int v0 = m_fixVertexOnLinePairs[i].first;
			int v1 = m_fixVertexOnLinePairs[i].second;
			double radius = 0.001;
			m_primitiveDrawer->drawSphere(radius, visualize_m_x.segment<3>(3 * v0), Eigen::Vector4d(0.8, 0.8, 0.8, 1.0));
			m_primitiveDrawer->drawSphere(radius, visualize_m_x.segment<3>(3 * v1), Eigen::Vector4d(0.8, 0.8, 0.8, 1.0));
		}
	}

	//draw external force arrow
	for (int i = 0; i < m_vertexIndex_force.size(); i++)
	{
		//if (m_vertexIndex_force[i].first == 5 || m_vertexIndex_force[i].first == 115)//torusKnot
		{
			if (m_vertexIndex_force[i].second.norm() == 0.0)
				continue;
			double radius = 0.0008;
			double arrowLength = 0.008;
			Eigen::Vector4d forceColor(0.8, 0, 0, 1);


			m_primitiveDrawer->drawArrow(radius,
				visualize_m_x.segment<3>(3 * m_vertexIndex_force[i].first) - m_vertexIndex_force[i].second.normalized() * arrowLength,
				visualize_m_x.segment<3>(3 * m_vertexIndex_force[i].first),
				forceColor);
		}
	}

	if (m_drawInteractivePairs)
	{
		/*for (int i = 0; i < m_pinningPairConstraint_set.size(); i++)
		{
			const auto& c_indices = m_pinningPairConstraint_set[i].vertex_indices(m_ipc_collisionMesh.edges(), m_ipc_collisionMesh.faces());
			for (int j = 0; j < c_indices.size(); j++)
			{
				m_primitiveDrawer->drawSphere(0.001, visualize_m_x.segment<3>(c_indices[j] * 3), m_pinningPairConstraintColor[i]);
			}
		}*/
	}

	if (m_enableInteractiveSelection)
	{
		for (int i = 0; i < m_selectedPairs.size(); i++)
		{
			if (m_selectedPairs[i].first != -1)
				m_primitiveDrawer->drawSphere(0.001, visualize_m_x.segment<3>(m_selectedPairs[i].first * 3), m_pinningPairConstraintColor[i]);
			if (m_selectedPairs[i].second != -1)
				m_primitiveDrawer->drawSphere(0.001, visualize_m_x.segment<3>(m_selectedPairs[i].second * 3), m_pinningPairConstraintColor[i]);
		}
	}

	if (m_showTargetSphere)
	{
		m_primitiveDrawer->drawSphere(m_targetSphereRadius, m_targetSphereCenter, Eigen::Vector4d(0.7, 0.7, 0.7, 1.0), 2);
	}

	/*if(eigen_modal_analysis.size()>0)
		for (int i = 0; i < m_nv; i++)
		{
			Eigen::Vector3d dir = eigen_modal_analysis[m_modesAnimationID].second.segment<3>(3 * i);
			double radius = 0.0008;
			double arrowLength = 0.15;
			m_primitiveDrawer->drawArrow(radius,
				visualize_m_x.segment<3>(3 * i),
				visualize_m_x.segment<3>(3 * i) + dir * arrowLength,
				Eigen::Vector4d(0, 1, 0, 1));
		}*/
}

void SOBendingModeShellStructure::reset()
{
	time_passed = 0.0;
	Eigen::VectorXd s;
	reset(s, false, false);
}

void SOBendingModeShellStructure::reset(Eigen::VectorXd& new_status)
{
	time_passed = 0.0;
	reset(new_status, false, false);
}

void SOBendingModeShellStructure::reset(Eigen::VectorXd& new_status, bool build_graph, bool ls)
{
	cout << "--------------ShellStructure has been reset--------------" << endl;

	//basic reset
	SODeformableObject::reset();
	time_passed = 0.0;

	m_V.setZero();

	m_sensitivityAnalysisProjectedState = m_para_X;
	reInit();
}

void SOBendingModeShellStructure::reInit()
{
	m_V.setZero();

	m_faceColors = m_triMesh->faceColors();


	m_v.setZero();
	m_appliedForces.setZero();

	setVisualizationConfiguration(m_x);
}


void SOBendingModeShellStructure::mouseLeftClick(const Eigen::Vector3d& src, const Eigen::Vector3d& ray)
{
	if (m_enableInteractiveSelection)
	{
		spdlog::info("src {} {} {}, ray {} {} {}", src[0], src[1], src[2], ray[0], ray[1], ray[2]);
		if (m_selectMeshVertices.size() == 0 || m_selectMeshVertices != m_triMesh->vertices())
		{//build a new tree
			if (m_selectMeshAABBTree != nullptr)
			{
				delete m_selectMeshAABBTree;
				m_selectMeshAABBTree = nullptr;
			}

			m_selectMeshVertices = m_triMesh->vertices();
			std::vector<int> faces = m_triMesh->faces();

			int numVertices = m_selectMeshVertices.size() / 3;
			int numFaces = faces.size() / 3;


			Eigen::MatrixXi limitingSpace_triFaces(numFaces, 3);


			std::vector<SimOpt::collisionElement*> triangleList;
			triangleList.reserve(numFaces);

			for (unsigned int i = 0; i < numFaces; ++i)
			{
				const unsigned int& indexA = faces[i * 3 + 0];
				const unsigned int& indexB = faces[i * 3 + 1];
				const unsigned int& indexC = faces[i * 3 + 2];

				limitingSpace_triFaces.row(i) = Eigen::RowVector3i(indexA, indexB, indexC);

				const Eigen::Vector3d v1 = m_selectMeshVertices.segment<3>(indexA * 3);
				const Eigen::Vector3d v2 = m_selectMeshVertices.segment<3>(indexB * 3);
				const Eigen::Vector3d v3 = m_selectMeshVertices.segment<3>(indexC * 3);

				Eigen::Vector3d face_norm = ((v2 - v1).cross(v3 - v1)).normalized();

				SimOpt::CTriangle* tri = new SimOpt::CTriangle(face_norm, v1, v2, v3, limitingSpace_triFaces.row(i), i);
				//tri.setIndex(i);

				triangleList.push_back(tri);

			}

			// setup aabb parameters
			SimOpt::Heuristic heurdata(14, 1, 1, 0.01f);
			m_selectMeshAABBTree = new SimOpt::AABBTree();
			m_selectMeshAABBTree->build_tree(triangleList, 0, heurdata);
		}


		if (m_selectMeshAABBTree)
		{
			SimOpt::CRay cray(src, src + ray);
			SimOpt::HitRecord hitRecords; hitRecords.set_zero();
			std::vector<const SimOpt::collisionElement*> allIntersectionTriangleList;
			const SimOpt::collisionElement* frontCollidingTriangle = nullptr;
			m_selectMeshAABBTree->get_root().CollidingTriangles(cray, allIntersectionTriangleList, hitRecords, frontCollidingTriangle);

			if (allIntersectionTriangleList.size() > 0 && hitRecords.t != DBL_MAX && frontCollidingTriangle != nullptr)
			{

				Eigen::VectorXi triangleVericesIndice = frontCollidingTriangle->getVerticesIndices();
				int closestVertexIndex = -1;
				double closestDistance = 1e7;
				for (int i = 0; i < triangleVericesIndice.size(); i++)
				{
					double dis = (m_selectMeshVertices.segment<3>(triangleVericesIndice[i] * 3) - hitRecords.position).norm();
					if (dis < closestDistance)
					{
						closestDistance = dis;
						closestVertexIndex = triangleVericesIndice[i];
					}
				}

				//add to m_selectedPairs
				std::pair<int, int>* addToPair = nullptr;
				if (m_selectedPairs.size() > 0)
				{
					std::pair<int, int>* last_pair = &m_selectedPairs.back();
					if (last_pair->first == -1 || last_pair->second == -1)
						addToPair = last_pair;
					else
					{
						m_selectedPairs.emplace_back(-1, -1);
						addToPair = &m_selectedPairs.back();

						m_pinningPairConstraintColor.push_back(Eigen::Vector4d(((double)rand() / (RAND_MAX)), ((double)rand() / (RAND_MAX)), ((double)rand() / (RAND_MAX)), 1.0));
					}
				}
				else
				{
					m_selectedPairs.emplace_back(-1, -1);
					addToPair = &m_selectedPairs.back();

					m_pinningPairConstraintColor.push_back(Eigen::Vector4d(((double)rand() / (RAND_MAX)), ((double)rand() / (RAND_MAX)), ((double)rand() / (RAND_MAX)), 1.0));
				}


				if (addToPair->first == -1)
					addToPair->first = closestVertexIndex;
				else
					addToPair->second = closestVertexIndex;

				spdlog::info("Add vertex index {}, pair is {}, {}", closestVertexIndex, addToPair->first, addToPair->second);
			}

		}
	}
}

Eigen::VectorXd SOBendingModeShellStructure::optimizeOneTimeFoldingParametersWithinSphereWithRigidTranformation(const Eigen::VectorXd& initialParameters, const Eigen::MatrixXd& strainSpaceModes, Eigen::Vector3d& sphere_center, double sphere_radius, std::vector<Eigen::VectorXd>& intermediateStates)
{
	//intermediateStates.clear();
	assert(m_statesAlongModesCoefficients.size() > 0);

	Eigen::VectorXd last_state = m_statesAlongModesCoefficients.back();// m_para_X;
	Eigen::VectorXd mode_0(m_hingeElements.size());
	for (int i = 0; i < m_hingeElements.size(); i++)
		mode_0[i] = m_hingeElements[i].getInitialRestAngle();

	auto computeObjective = [&](const Eigen::VectorXd& p)
		{
			Eigen::VectorXd mode_p = mode_0 + strainSpaceModes * p.head(p.size() - 6);

			for (int i = 0; i < m_hingeElements.size(); i++)
			{
				m_hingeElements[i].setRestAngle(mode_p[i]);
			}

			m_x = last_state;
			projectToManifold();

			//compute rigid transformation for m_x
			Eigen::Vector3d translation = p.tail(3);
			Eigen::Vector3d rotation_vector = p.tail(6).head(3);
			Eigen::Matrix3d rotationMatrix = soutil::construct_rotation_matrix(rotation_vector);

			const Eigen::MatrixXd& m_x_map = Eigen::Map<const Eigen::MatrixXd>(m_x.data(), 3, m_x.size() / 3);

			Eigen::MatrixXd transformed_m_x = rotationMatrix * m_x_map;
			transformed_m_x.colwise() += translation;

			Eigen::MatrixXd distanceToCenter = (transformed_m_x.colwise() - sphere_center);
			double objective = 0.0;
			for (int v_i = 0; v_i < distanceToCenter.cols(); v_i++)
			{
				double dis = distanceToCenter.col(v_i).squaredNorm() - sphere_radius * sphere_radius;
				if (dis > 0)
					objective += dis;
			}

			return objective;
		};
	auto computeObjectiveGradient = [&](const Eigen::VectorXd& p, Eigen::VectorXd& grad)
		{
			Eigen::VectorXd mode_p = mode_0 + strainSpaceModes * p.head(p.size() - 6);

			for (int i = 0; i < m_hingeElements.size(); i++)
			{
				m_hingeElements[i].setRestAngle(mode_p[i]);
			}

			m_x = last_state;
			projectToManifold();
			last_state = m_x;

			//compute rigid transformation for m_x
			Eigen::Vector3d translation = p.tail(3);
			Eigen::Vector3d rotation_vector = p.tail(6).head(3);
			Eigen::Matrix3d rotationMatrix = soutil::construct_rotation_matrix(rotation_vector);

			const Eigen::MatrixXd& m_x_map = Eigen::Map<const Eigen::MatrixXd>(m_x.data(), 3, m_x.size() / 3);

			Eigen::MatrixXd transformed_m_x = rotationMatrix * m_x_map;
			transformed_m_x.colwise() += translation;
			intermediateStates.push_back(Eigen::Map<Eigen::VectorXd>(transformed_m_x.data(), transformed_m_x.size()));

			//objective gradient
			Eigen::MatrixXd distanceToCenter = (transformed_m_x.colwise() - sphere_center);
			double objective = 0.0;
			Eigen::VectorXd objective_grad(distanceToCenter.size()); objective_grad.setZero();
			for (int v_i = 0; v_i < distanceToCenter.cols(); v_i++)
			{
				double dis = distanceToCenter.col(v_i).squaredNorm() - sphere_radius * sphere_radius;
				if (dis > 0)
				{
					objective += dis;
					objective_grad.segment<3>(v_i * 3) = 2.0 * distanceToCenter.col(v_i);
				}
			}

			//sensitvity analysis to compute dxd p.head
			SpMat dxdtheta_rest;
			computeStrainSpaceSensitiveMatrix(m_x, m_para_X, dxdtheta_rest);
			auto dxdphead = dxdtheta_rest * strainSpaceModes;


			std::vector<Eigen::Matrix3d> rotationMatrix_jacobian;
			soutil::construct_rotation_matrix_jacobian(rotation_vector, rotationMatrix_jacobian);


			tbb::enumerable_thread_specific<TripVec> storage;
			tbb::parallel_for(
				tbb::blocked_range<size_t>(size_t(0), transformed_m_x.cols()),
				[&](const tbb::blocked_range<size_t>& r) {
					auto& local_dxdp_triplets = storage.local();

					for (size_t v_i = r.begin(); v_i < r.end(); v_i++)
						//for (int v_i = 0; v_i < transformed_m_x.cols(); v_i++)
					{
						// dR / dr * x
						Eigen::Vector3d transformed_vi = transformed_m_x.col(v_i).head(3);
						Eigen::Matrix3d dR_dr_mul_vi;
						soutil::tensorMulVector1<3, 3, 3, 3>(rotationMatrix_jacobian, transformed_vi, dR_dr_mul_vi);
						for (int d_i = 0; d_i < 3; d_i++)
						{
							for (int d_j = 0; d_j < 3; d_j++)
							{
								if (dR_dr_mul_vi(d_i, d_j) != 0.0)
									local_dxdp_triplets.emplace_back(v_i * 3 + d_i, p.size() - 6 + d_j, dR_dr_mul_vi(d_i, d_j));
							}
						}

						// R * dx/dp
						Eigen::MatrixXd R_dvidphead = rotationMatrix * dxdphead.middleRows(3 * v_i, 3);
						//std::cout << "nonzero rate for v_i " << v_i << "/" << transformed_m_x.cols() << ": " << double(R_dvidphead.nonZeros()) / double(R_dvidphead.size()) << std::endl;
						for (int d_i = 0; d_i < 3; d_i++)
						{
							for (int p_j = 0; p_j < p.size() - 6; p_j++)
							{
								if (R_dvidphead(d_i, p_j) != 0.0)
									local_dxdp_triplets.emplace_back(v_i * 3 + d_i, p_j, R_dvidphead(d_i, p_j));
							}
						}

						//dT / dp
						for (int d_i = 0; d_i < 3; d_i++)
						{
							local_dxdp_triplets.emplace_back(3 * v_i + d_i, p.size() - 3 + d_i, 1.0);
						}
					}
				});
			SpMat dxdp(transformed_m_x.cols() * 3, p.size());
			for (const auto& local_dxdp_triplets : storage) {
				Eigen::SparseMatrix<double> local_dxdp(transformed_m_x.cols() * 3, p.size());
				local_dxdp.setFromTriplets(
					local_dxdp_triplets.begin(), local_dxdp_triplets.end());
				dxdp += local_dxdp;
			}
			grad = dxdp.transpose() * objective_grad;

			spdlog::info("Objective is {} with its gradient norm {}", objective, grad.norm());
		};

	/*{
		SimOpt::DerivativeTester tester;
		Eigen::VectorXd grad;
		computeObjectiveGradient(initialParameters, grad);
		tester.testGradient(initialParameters, grad, computeObjective, 10, 1e-3);
		exit(0);
	}*/

	auto computeLargestStepSize = [&](const Eigen::VectorXd& p, const Eigen::VectorXd& dp)
		{
			double alpha = 1.0;
			Eigen::VectorXd dtheta = strainSpaceModes * dp.head(dp.size() - 6);
			while (dtheta.norm() * alpha > M_PI / 4.0)
			{
				alpha *= 0.5;
			}
			std::cout << "largest alpha is " << alpha;
			return alpha;
		};

	SimOpt::UnconstrainedLBFGS lbfgs_solver;
	lbfgs_solver.setObjective(computeObjective);
	lbfgs_solver.setObjectiveGradient(computeObjectiveGradient);
	lbfgs_solver.setLargestLineSearchStepFunctor(computeLargestStepSize);
	lbfgs_solver.setX(initialParameters);
	lbfgs_solver.setMaxIter(5000);//m_newtonMaxIterations
	lbfgs_solver.setGradientThreshold(1e-5);//m_newtonSolveResidual
	lbfgs_solver.solve();

	double finalObjective = computeObjective(lbfgs_solver.getX());
	spdlog::info("Final objective {}", finalObjective);
	std::cout << "with parameters " << lbfgs_solver.getX().transpose() << std::endl;

	if (lbfgs_solver.state() != SimOpt::UnconstrainedLBFGS::SOLVED)
	{
		printf("Static solve didn't converge to %.10f\n", m_newtonSolveResidual);//, the rhs is %.10f , solver.getLastGradientNorm()

		return lbfgs_solver.getX();
	}
	else
	{

		return lbfgs_solver.getX();
	}
}

Eigen::VectorXd SOBendingModeShellStructure::optimizePositionsToTargetIPCSphere(Eigen::VectorXd& vX, Eigen::Vector3d& sphere_center, double sphere_radius, std::vector<Eigen::VectorXd>& intermediateStates)
{
	for (auto& elem : m_fixDoFElements)
	{
		elem.updateApplyElement(false);
	}
	Eigen::MatrixXd vX_matrix = Eigen::Map<Eigen::MatrixXd>(vX.data(), 3, vX.size() / 3);
	
	double radius = vX_matrix.colwise().norm().maxCoeff() + 1e-3;

	ipc::CollisionMesh mesh_withSphere;
	Eigen::MatrixXd sphereVertices;
	Eigen::MatrixXi sphereTriangles;
	SimpleTriMesh::icosphere(3).toMatrices(sphereVertices, sphereTriangles);
	sphereVertices.transposeInPlace();
	sphereTriangles.transposeInPlace();
	{
		std::vector<int> faces = m_triMesh->faces();
		int numFaces = faces.size() / 3;
		Eigen::MatrixXi selfContact_triFaces(numFaces, 3);
		for (unsigned int i = 0; i < numFaces; ++i)
		{
			const unsigned int& indexA = m_triMesh->faces()[i * 3 + 0];
			const unsigned int& indexB = m_triMesh->faces()[i * 3 + 1];
			const unsigned int& indexC = m_triMesh->faces()[i * 3 + 2];

			selfContact_triFaces.row(i) = Eigen::RowVector3i(indexA, indexB, indexC);
		}
		sphereTriangles.rowwise() += Eigen::RowVector3i(m_x.size() / 3, m_x.size() / 3, m_x.size() / 3);
		Eigen::MatrixXi collisionTriangles(numFaces + sphereTriangles.rows(), 3);
		collisionTriangles.topRows(numFaces) = selfContact_triFaces;
		collisionTriangles.bottomRows(sphereTriangles.rows()) = sphereTriangles;

		//so we use CCD in IPC
		Eigen::MatrixXi collisionEdges;
		if (collisionTriangles.size())
			igl::edges(collisionTriangles, collisionEdges);

		Eigen::MatrixXd selfContact_V = (Eigen::Map<const Eigen::MatrixXd>(m_x.data(), dimension, m_x.size() / dimension)).transpose();
		Eigen::MatrixXd collisonVertices(selfContact_V.rows() + sphereVertices.rows(), 3);
		collisonVertices.topRows(selfContact_V.rows()) = selfContact_V;
		collisonVertices.bottomRows(sphereVertices.rows()) = sphereVertices * radius;
		mesh_withSphere = ipc::CollisionMesh(collisonVertices, collisionEdges, collisionTriangles);
		//auto can_collide = [this](size_t vertex1, size_t vertex2) {
		//	return true;
		//};
		//m_ipc_collisionMesh.can_collide = can_collide;
		ipc::logger().set_level(spdlog::level::warn);

		if (ipc::has_intersections(mesh_withSphere, collisonVertices))
		{
			spdlog::error("The initial green lagrange strain space has intersection");
		}
	}

	//intermediateStates.clear();
	assert(m_statesAlongModesCoefficients.size() > 0);

	Eigen::VectorXd last_state = m_statesAlongModesCoefficients.back();// m_para_X;

	Eigen::VectorXd _x(last_state.size() + 1);
	_x.head(last_state.size()) = last_state;
	_x[_x.size() - 1] = radius;

	auto computeCollisonVertices = [this, &mesh_withSphere, sphereVertices, sphere_center](const Eigen::VectorXd& vx)
		{
			Eigen::MatrixXd V_vx = (Eigen::Map<const Eigen::MatrixXd>(vx.data(), dimension, vx.size() / dimension)).transpose();
			Eigen::MatrixXd V = mesh_withSphere.vertices_at_rest();
			V.topRows(V_vx.rows()) = V_vx;
			double r = vx[vx.size() - 1];
			V.bottomRows(sphereVertices.rows()) = ((sphereVertices * r).rowwise() + sphere_center.transpose());
			return V;
		};
	auto stepBeginningFunction = [this, vX, &mesh_withSphere, computeCollisonVertices](const Eigen::VectorXd& vx, int currentIteration)
		{
			if (USE_IPC)
			{
				ipc::Constraints constraint_set;

				Eigen::MatrixXd V = computeCollisonVertices(vx);
				constraint_set.build(mesh_withSphere, V, m_ipc_dhat);

				//in the gradient computation, we will update barrier stiffness
				//since objective and hessian are computed afterward in the same step
				double min_distance = ipc::compute_minimum_distance(mesh_withSphere, V, constraint_set);
				spdlog::debug("Minimum distance is {} with number of constraints {}", min_distance, constraint_set.size());

				double bbox_diagonal_length = ipc::world_bbox_diagonal_length(V);
				m_ipc_adaptiveStiffness = ipc::update_barrier_stiffness(m_ipc_minDistance, min_distance,
					m_ipc_upperBoundStiffness, m_ipc_adaptiveStiffness, bbox_diagonal_length, 1e-5);
				spdlog::debug("Update barrier stiffness to {} with upper bound {}", m_ipc_adaptiveStiffness, m_ipc_upperBoundStiffness);
				m_ipc_minDistance = min_distance;
			}
			visualize_m_x = vx;
		};

	auto computeStaticObjective = [this, vX, &mesh_withSphere, sphere_radius, computeCollisonVertices](const Eigen::VectorXd& vx)
		{
			//internal energy
			tbb::enumerable_thread_specific<double> storage(0);
			tbb::parallel_for(
				tbb::blocked_range<size_t>(size_t(0), m_elements.size()),
				[&](const tbb::blocked_range<size_t>& r) {
					auto& local_Eint = storage.local();
					for (size_t i = r.begin(); i < r.end(); i++)
						//for (size_t i = 0; i < m_elements.size(); i++)
					{
						local_Eint += m_elements[i]->computeEnergy(vx, vX);
						//elementEnergy[i] = m_elements[i]->computeEnergy(vx, vX);
					}
				});
			double Eint = 0.0;
			for (const auto& local_Eint : storage) {
				Eint += local_Eint;
			}

			//gravity energy
			double Egrav = 0.0;
			if (USE_GRAVITY)
			{
				Egrav += -(computeGravityForce().dot(vx));
			}

			double EappliedForces = 0.0;
			EappliedForces -= m_appliedForces.dot(vx.head(vX.size()) - vX);


			//global contact barrier
			double EselfContact = 0.0;
			if (USE_IPC)
			{
				ipc::Constraints constraint_set;
				Eigen::MatrixXd V = computeCollisonVertices(vx);
				constraint_set.build(mesh_withSphere, V, m_ipc_dhat);
				double barrierPotential = ipc::compute_barrier_potential(mesh_withSphere, V, constraint_set, m_ipc_dhat);

				EselfContact += m_ipc_adaptiveStiffness * barrierPotential;
			}

			//sphere radius target
			double E_sphere = 0.5 * m_fix_stiffness * (vx[vx.size() - 1] - sphere_radius) * (vx[vx.size() - 1] - sphere_radius);

			double E = Eint + Egrav + EappliedForces + EselfContact + E_sphere;

			return E;

		};

	auto computeStaticGradient = [this, vX, &intermediateStates, &mesh_withSphere, sphereVertices, sphere_radius, computeCollisonVertices, computeStaticObjective](const Eigen::VectorXd& vx, Eigen::VectorXd& grad)
		{
			intermediateStates.push_back(vx);

			grad = Eigen::VectorXd(vx.size());
			grad.setZero();

			//compute elements gradient
			tbb::enumerable_thread_specific<Eigen::VectorXd> storage(
				Eigen::VectorXd::Zero(vx.size()));

			tbb::parallel_for(
				tbb::blocked_range<size_t>(size_t(0), m_elements.size()),
				[&](const tbb::blocked_range<size_t>& r) {
					auto& local_grad = storage.local();
					for (size_t i = r.begin(); i < r.end(); i++)
					{
						Eigen::VectorXd element_grad;
						m_elements[i]->computeGradient(vx, vX, element_grad);

						const auto element_dofs = m_elements[i]->getDofIndices();
						for (int dof_i = 0; dof_i < element_dofs.size(); dof_i++)
						{
							local_grad[element_dofs[dof_i]] += element_grad[dof_i];
						}
					}
				});
			for (const auto& local_grad : storage) {
				grad += local_grad;
			}

			if (USE_GRAVITY)
			{
				grad += -computeGravityForce();
			}
			grad.head(m_appliedForces.size()) += -m_appliedForces;


			//Eigen::VectorXd strainLimitingGrad(vx.size());
			//computeAllStrainLimitingPotentialGradient(vx, vX, strainLimitingGrad);
			Eigen::VectorXd EselfContactGrad(vx.size()); EselfContactGrad.setZero();
			if (USE_IPC)
			{
				ipc::Constraints constraint_set;
				Eigen::MatrixXd V = computeCollisonVertices(vx);

				constraint_set.build(mesh_withSphere, V, m_ipc_dhat);
				Eigen::VectorXd dEcontactdx = ipc::compute_barrier_potential_gradient(mesh_withSphere, V, constraint_set, m_ipc_dhat);

				Eigen::VectorXd barrierPotentialGradient(vx.size()); barrierPotentialGradient.setZero();
				barrierPotentialGradient.head(vx.size() - 1) = dEcontactdx.head(vx.size() - 1);
				Eigen::MatrixXd sphereVertices_transpose = sphereVertices.transpose();
				Eigen::VectorXd sphereVertices_vector = Eigen::Map<Eigen::VectorXd>(sphereVertices_transpose.data(), sphereVertices_transpose.size());
				barrierPotentialGradient[barrierPotentialGradient.size() - 1] = dEcontactdx.tail(sphereVertices.rows() * 3).transpose() * sphereVertices_vector;
				EselfContactGrad += m_ipc_adaptiveStiffness * barrierPotentialGradient;
			}
			grad += EselfContactGrad;

			//sphere radius target
			grad[vx.size() - 1] += m_fix_stiffness * (vx[vx.size() - 1] - sphere_radius);

			filter(grad);

		};

	auto computeStaticHessian = [this, vX, &mesh_withSphere, sphereVertices, sphere_radius, computeCollisonVertices, computeStaticGradient](const Eigen::VectorXd& vx, SpMat& Hessian)
		{
			Hessian = SpMat(vx.size(), vx.size());
			Hessian.setZero();

			tbb::enumerable_thread_specific<TripVec> storage;
			tbb::parallel_for(
				tbb::blocked_range<size_t>(size_t(0), m_elements.size()),
				[&](const tbb::blocked_range<size_t>& r) {
					auto& local_hess_triplets = storage.local();

					for (size_t i = r.begin(); i < r.end(); i++)
					{
						Eigen::MatrixXd element_hess;
						m_elements[i]->computeHessian(vx, vX, element_hess);

						//soutil::makePD<double, Eigen::Dynamic>(element_hess, USE_MAKEPD);

						const std::vector<int>& element_dofs = m_elements[i]->getDofIndices();
						for (int dof_i = 0; dof_i < element_dofs.size(); dof_i++)
						{
							for (int dof_j = 0; dof_j < element_dofs.size(); dof_j++)
							{
								local_hess_triplets.emplace_back(element_dofs[dof_i], element_dofs[dof_j],
									element_hess(dof_i, dof_j));
							}
						}
					}
				});
			for (const auto& local_hess_triplets : storage) {
				Eigen::SparseMatrix<double> local_hess(vx.size(), vx.size());
				local_hess.setFromTriplets(
					local_hess_triplets.begin(), local_hess_triplets.end());
				Hessian += local_hess;
			}


			SpMat EselfContactHessian(vx.size(), vx.size());
			if (USE_IPC)
			{
				ipc::Constraints constraint_set;
				Eigen::MatrixXd V = computeCollisonVertices(vx);

				constraint_set.build(mesh_withSphere, V, m_ipc_dhat);
				Eigen::VectorXd dEcontactdx = ipc::compute_barrier_potential_gradient(mesh_withSphere, V, constraint_set, m_ipc_dhat);
				SpMat d2Econtactdx2 = ipc::compute_barrier_potential_hessian(mesh_withSphere, V, constraint_set, m_ipc_dhat, false);

				//d^2 E / d vx^2
				SpMat d2Edvx2 = d2Econtactdx2.topLeftCorner(vx.size() - 1, vx.size() - 1);
				//d^2 E / d vx dvs
				SpMat d2Edvxdvs = d2Econtactdx2.topRightCorner(vx.size() - 1, sphereVertices.rows() * 3);
				//d^2 E / dvs^2
				SpMat d2Edvs2 = d2Econtactdx2.bottomRightCorner(sphereVertices.rows() * 3, sphereVertices.rows() * 3);

				//compute final hessians
				Eigen::MatrixXd sphereVertices_transpose = sphereVertices.transpose();
				Eigen::VectorXd sphereVertices_vector = Eigen::Map<Eigen::VectorXd>(sphereVertices_transpose.data(), sphereVertices_transpose.size());

				Eigen::VectorXd d2Edvxdr = d2Edvxdvs * sphereVertices_vector;
				double d2Edr2 = sphereVertices_vector.transpose() * d2Edvs2 * sphereVertices_vector;
				

				TripVec barrierPotentialHessian_triplets;
				for (int k = 0; k < d2Edvx2.outerSize(); ++k)
					for (SpMat::InnerIterator it(d2Edvx2, k); it; ++it)
					{
						barrierPotentialHessian_triplets.emplace_back(it.row(), it.col(), it.value());
					}
				for (int i = 0; i < d2Edvxdr.size(); i++)
				{
					barrierPotentialHessian_triplets.emplace_back(i, vx.size() - 1, d2Edvxdr[i]);
					barrierPotentialHessian_triplets.emplace_back(vx.size() - 1, i, d2Edvxdr[i]);
				}
				barrierPotentialHessian_triplets.emplace_back(vx.size() - 1, vx.size() - 1, d2Edr2);
				SpMat barrierPotentialHessian(vx.size(), vx.size());
				barrierPotentialHessian.setFromTriplets(barrierPotentialHessian_triplets.begin(), barrierPotentialHessian_triplets.end());

				EselfContactHessian += m_ipc_adaptiveStiffness * barrierPotentialHessian;
			}

			Hessian += EselfContactHessian;

			//sphere radius target
			Hessian.coeffRef(vx.size() - 1, vx.size() - 1) += m_fix_stiffness;

			filterMat(Hessian);

		};

	auto computeLargestStepSize = [this, vX, &mesh_withSphere, computeCollisonVertices](const Eigen::VectorXd& vx, const Eigen::VectorXd& dx)
		{
			spdlog::debug("Try to find a largest step size");


			double largestStepSize = 1.0;
			if (USE_IPC)
			{
				//we use CCD to find the largest step size
				Eigen::MatrixXd V0 = computeCollisonVertices(vx);

				Eigen::MatrixXd V1 = computeCollisonVertices(vx + largestStepSize * dx);

				//constrain initial step size
				//while (strainSpaceAllVertices_x1.bottomRows(m_cstShellElements.size()).maxCoeff() > fixedVertices.maxCoeff())
				/*while (V1.rowwise().norm().maxCoeff() > V0.rowwise().norm().maxCoeff())
				{
					largestStepSize *= 0.5;
					V1 = (Eigen::Map<const Eigen::MatrixXd>((vx + largestStepSize * dx).eval().data(), dimension, vx.size() / dimension)).transpose();
				}*/
				while ((largestStepSize * dx).norm() > ipc::world_bbox_diagonal_length(V0))
				{
					//spdlog::debug("Search dir length {} and bbox diagonal length {}", (largestStepSize* dx).norm(), ipc::world_bbox_diagonal_length(V0));
					largestStepSize *= 0.5;
				}
				V1 = computeCollisonVertices(vx + largestStepSize * dx);

				spdlog::debug("Use CCD to find a largest step size with initial step size {}, in-plane", largestStepSize);

				double ipc_collisionFreeLargestStepSize = ipc::compute_collision_free_stepsize(mesh_withSphere,
					V0, V1);
				largestStepSize *= ipc_collisionFreeLargestStepSize;
				/*V1 = (Eigen::Map<const Eigen::MatrixXd>((vx + largestStepSize * dx).eval().data(), dimension, vx.size() / dimension)).transpose();

				//one more step to make sure no collision within this step
				while (!ipc::is_step_collision_free(m_ipc_collisionMesh, V0, V1))
				{
					largestStepSize *= 0.5;
					V1 = (Eigen::Map<const Eigen::MatrixXd>((vx + largestStepSize * dx).eval().data(), dimension, vx.size() / dimension)).transpose();
				}*/
				spdlog::debug("Compute largest step size done with step size {}, in-plane", largestStepSize);
			}
			return largestStepSize;
		};

	/*{
		SimOpt::DerivativeTester tester;
		Eigen::VectorXd grad;
		computeStaticGradient(_x, grad);
		tester.testGradient(_x, grad, computeStaticObjective,10);
		exit(0);
		SpMat hessian;
		computeStaticHessian(_x, hessian);
		tester.testJacobian(_x, hessian, computeStaticGradient,10);
		exit(0);
	}*/


	std::unique_ptr<SimOpt::UnconstrainedNewtonLinearSolver> linearSolver(new SimOpt::UnconstrainedNewtonLinearSolverLLT(computeStaticHessian));
	SimOpt::UnconstrainedNewton solver(std::move(linearSolver));
	solver.setX(_x);
	solver.setObjective(computeStaticObjective);
	solver.setObjectiveGradient(computeStaticGradient);
	solver.setLargestLineSearchStepFunctor(computeLargestStepSize);
	solver.setStepBeginningFunctor(stepBeginningFunction);
	solver.setMaxIter(5000);//m_newtonMaxIterations
	//solver.setLinesearchMaxIter(m_newtonMaxLinesearchSteps);
	solver.setGradientThreshold(1e-7);//m_newtonSolveResidual
	solver.solve();

	if (solver.state() != SimOpt::UnconstrainedNewton::SOLVED)
	{
		printf("Static solve didn't converge to %.10f\n", m_newtonSolveResidual);//, the rhs is %.10f , solver.getLastGradientNorm()
	}

	/*{
		SpMat testHessian;
		bool temp = USE_MAKEPD;
		USE_MAKEPD = false;
		computeStaticHessian(solver.getX(), testHessian);
		USE_MAKEPD = temp;
		//add a very small regularization
		for (int i = 0; i < m_x.size(); i++)
			testHessian.coeffRef(i, i) += 1e-6;
		Eigen::PardisoLLT<SpMat> sp_solver;
		sp_solver.compute(testHessian);
		if (sp_solver.info() != Eigen::Success)
			spdlog::error("Final hessian is not PSD");
		else
			spdlog::info("Final hessian is PSD");
	}*/

	//after taking a step, we update states
	m_x = solver.getX();


	for (auto& elem : m_fixDoFElements)
	{
		elem.updateApplyElement(true);
	}
	return m_x;
}

Eigen::VectorXd SOBendingModeShellStructure::applyMassMatrix(const Eigen::VectorXd& vf)
{
	Eigen::VectorXd result(vf.size()); result.setZero();
	for (int i = 0; i < m_mass.size(); i++)
	{
		result[i] = m_mass[i] * vf[i];
	}
	return result;
}
Eigen::VectorXd SOBendingModeShellStructure::applyInvMassMatrix(const Eigen::VectorXd& vf)
{
	Eigen::VectorXd result(vf.size()); result.setZero();
	for (int i = 0; i < m_mass.size(); i++)
	{
		result[i] = vf[i] / m_mass[i];
	}
	return result;
}
void SOBendingModeShellStructure::addMassMatrix(Eigen::SparseMatrix<double>& mat)
{
	for (int i = 0; i < m_mass.size(); i++)
	{
		mat.coeffRef(i, i) += m_mass[i];
	}
}
SpMat SOBendingModeShellStructure::computeMassMatrix()
{
	int N = getDoFCount();
	SpMat mass_matrix(N, N);
	mass_matrix.reserve(Eigen::VectorXi::Constant(N, 1)); //reserve 1 non zero entry per column
	for (int i = 0; i < N; i++)
	{
		mass_matrix.insert(i, i) = m_mass[i];
	}
	mass_matrix.makeCompressed();

	return mass_matrix;
}
Eigen::VectorXd SOBendingModeShellStructure::computeGravityForce()
{
	return m_gravity;
}
