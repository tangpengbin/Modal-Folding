#include "SOReflectionBendingModeShellStructure.h"
#include <SimLib/Core/SOUtils.h>
#include <SimLib/Geom/Primitive.h>
#include <SimLib/Viewer/colormap.h>
#include <SimLib/Geom/SimpleTriMesh.h>

#include <fstream>
#include <stdlib.h>

#include <sys/types.h>
#include <sys/stat.h>
#ifdef _WIN32
#include <direct.h>
#elif __APPLE__ || __linux__
//#include <sys/stat.h>
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

#define USE_TBB
//ipc
#include "ipc/ipc.hpp"
#include "ipc/barrier/barrier.hpp"
#include <ipc/utils/world_bbox_diagonal_length.hpp>
#include <ipc/barrier/adaptive_stiffness.hpp>
#include <ipc/utils/local_to_global.hpp>
#include <ipc/distance/distance_type.hpp>
#include <ipc/distance/point_triangle.hpp>


#include <SimLib/Core/UnconstrainedNewton.h>


#include <Spectra/SymEigsShiftSolver.h>
#include <Spectra/MatOp/SparseSymShiftSolve.h>

#include "./Codegen/computeHingeEdgeAlighmentEnergy.h"
#include "./Codegen/computeHingeAngleAlighmentEnergy.h"
#include "./Codegen/computedTheta.h"

using namespace soutil;
const int dimension = 3;
SOReflectionBendingModeShellStructure::SOReflectionBendingModeShellStructure()
{
	m_fourCornerIndices = Eigen::Vector4i(-1, -1, -1, -1);
	m_numExtraDofs = 4;

	m_drawReflectionPatches = false;
}

SOReflectionBendingModeShellStructure::SOReflectionBendingModeShellStructure(const SOReflectionBendingModeShellStructure& other)
{
	(*this) = other;
	m_triMesh = new TriMesh;
	soutil::copyTriMesh(other.m_triMesh, m_triMesh);
	initial_mesh = other.initial_mesh;

	reconstructElementVector();

}

SOReflectionBendingModeShellStructure::~SOReflectionBendingModeShellStructure()
{
}

void SOReflectionBendingModeShellStructure::init(const TriMesh* pTriMesh, SOMaterial* mat,
	double thickness, 
	const std::vector<std::pair<int, int>>& fixVertexOnLinePairs, 
	const std::vector<std::pair<int, int>>& fixingVerticesDofs, 
	const std::vector<std::pair<int, Eigen::Vector3d>>& vertexIndex_force,
	const Eigen::Vector4d& fourReflectionPlane_rest,
	string outputPath)
{
	m_fourReflectionPlane_rest = fourReflectionPlane_rest;

	m_outputPath = outputPath;
	m_material = mat;
	m_cstThickness = thickness;
	m_triMesh = new TriMesh;
	soutil::copyTriMesh(pTriMesh, m_triMesh);
	initial_mesh = pTriMesh;
	m_fixVertexOnLinePairs = fixVertexOnLinePairs;
	for (int i = 0; i < fixingVerticesDofs.size(); i++)
		m_fixingDofs.push_back(fixingVerticesDofs[i].first * dimension + fixingVerticesDofs[i].second);
	m_vertexIndex_force = vertexIndex_force;
	m_meshDrawer->setMesh(m_triMesh);
	m_faceColors = m_triMesh->faceColors();


	m_elements.clear();
	m_cstShellElements.clear();
	m_hingeElements.clear();

	//copy positions from mesh
	m_X = m_triMesh->vertices();//m_X will keep the initial state
	m_nv = m_X.size()/3;
	processSystemDoFs();

	initElements(m_para_X);

	reconstructMassStuff(m_para_X);

	//init_newmark(m_x);

	
	initIPCConfig();

	setVisualizationConfiguration(m_x);

	spdlog::set_level(spdlog::level::info);
}

void SOReflectionBendingModeShellStructure::initIPCConfig()
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

	Eigen::VectorXd allDofs = convertSystemDofsToAllDofs(m_x);

	Eigen::MatrixXd selfContact_V = (Eigen::Map<const Eigen::MatrixXd>(allDofs.data(), dimension, allDofs.size() / dimension)).transpose();
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
	Eigen::VectorXd zeroEnergyGrad(allDofs);
	zeroEnergyGrad.setZero();

	double averageMass = m_mass.sum() / (dimension * m_nv);// *1000.0;

	m_ipc_dhat = 1e-4;
	m_ipc_adaptiveStiffness = ipc::initial_barrier_stiffness(bbox_diagonal_length, m_ipc_dhat, averageMass, zeroEnergyGrad, zeroEnergyGrad, m_ipc_upperBoundStiffness);
	m_ipc_adaptiveStiffness *= 1e2;
	m_ipc_upperBoundStiffness *= 1e6;
	spdlog::info("The initial ipc adaptive stiffness is {} with upper bound of {}", m_ipc_adaptiveStiffness, m_ipc_upperBoundStiffness);

}

void SOReflectionBendingModeShellStructure::processSystemDoFs()
{
	double numericalError = 1e-7;
	//   -----y1----
	//	 |			|
	//	 x1			x2
	//	 |			|
	//	 -----y2----
	double x1 = m_fourReflectionPlane_rest[0];//with order x1,y1,x2,y2
	double y1 = m_fourReflectionPlane_rest[1];
	double x2 = m_fourReflectionPlane_rest[2];
	double y2 = m_fourReflectionPlane_rest[3];
	int numSystemDofs = 0;
	for (int i = 0; i < m_nv; i++)
	{
		double pos_x = m_X(i * 3 + 0);
		double pos_y = m_X(i * 3 + 1);
		double diffx1 = abs(pos_x - x1);
		double diffy1 = abs(pos_y - y1);
		double diffx2 = abs(pos_x - x2);
		double diffy2 = abs(pos_y - y2);

		if (diffx1 < numericalError && diffy1 < numericalError)
		{//topleft corner
			boundaryVertices[i] = BoundaryVertexType::topleft;
			m_fourCornerIndices[0] = i;
			m_realVertex_to_systemDofAndNumDofs[i] = std::make_pair(numSystemDofs, 1);
			m_systemDof_to_realVertexAndNumDofs[numSystemDofs] = std::make_pair(i, 1);
			numSystemDofs += 1;
		}
		else if (diffx2 < numericalError && diffy1 < numericalError)
		{//topright corner
			boundaryVertices[i] = BoundaryVertexType::topright;
			m_fourCornerIndices[1] = i;
			m_realVertex_to_systemDofAndNumDofs[i] = std::make_pair(numSystemDofs, 1);
			m_systemDof_to_realVertexAndNumDofs[numSystemDofs] = std::make_pair(i, 1);
			numSystemDofs += 1;
		}
		else if (diffx1 < numericalError && diffy2 < numericalError)
		{//bottomleft corner
			boundaryVertices[i] = BoundaryVertexType::bottomleft;
			m_fourCornerIndices[2] = i;
			m_realVertex_to_systemDofAndNumDofs[i] = std::make_pair(numSystemDofs, 1);
			m_systemDof_to_realVertexAndNumDofs[numSystemDofs] = std::make_pair(i,1);
			numSystemDofs += 1;
		}
		else if (diffx2 < numericalError && diffy2 < numericalError)
		{//bottomright corner
			boundaryVertices[i] = BoundaryVertexType::bottomright;
			m_fourCornerIndices[3] = i;
			m_realVertex_to_systemDofAndNumDofs[i] = std::make_pair(numSystemDofs, 1);
			m_systemDof_to_realVertexAndNumDofs[numSystemDofs] = std::make_pair(i,1);
			numSystemDofs += 1;
		}
		else if (diffx1 < numericalError)
		{//left
			boundaryVertices[i] = BoundaryVertexType::left;
			m_realVertex_to_systemDofAndNumDofs[i] = std::make_pair(numSystemDofs, 2);
			m_systemDof_to_realVertexAndNumDofs[numSystemDofs] = std::make_pair(i,2);
			numSystemDofs += 2;
		}
		else if (diffx2 < numericalError)
		{//right
			boundaryVertices[i] = BoundaryVertexType::right;
			m_realVertex_to_systemDofAndNumDofs[i] = std::make_pair(numSystemDofs, 2);
			m_systemDof_to_realVertexAndNumDofs[numSystemDofs] = std::make_pair(i,2);
			numSystemDofs += 2;
		}
		else if (diffy1 < numericalError)
		{//top
			boundaryVertices[i] = BoundaryVertexType::top;
			m_realVertex_to_systemDofAndNumDofs[i] = std::make_pair(numSystemDofs, 2);
			m_systemDof_to_realVertexAndNumDofs[numSystemDofs] = std::make_pair(i,2);
			numSystemDofs += 2;
		}
		else if (diffy2 < numericalError)
		{//bottom
			boundaryVertices[i] = BoundaryVertexType::bottom;
			m_realVertex_to_systemDofAndNumDofs[i] = std::make_pair(numSystemDofs, 2);
			m_systemDof_to_realVertexAndNumDofs[numSystemDofs] = std::make_pair(i,2);
			numSystemDofs += 2;
		}
		else
		{//inner
			m_realVertex_to_systemDofAndNumDofs[i] = std::make_pair(numSystemDofs, 3);
			m_systemDof_to_realVertexAndNumDofs[numSystemDofs] = std::make_pair(i,3);
			numSystemDofs += 3;
		}
	}
	m_numSystemDofs = numSystemDofs + m_numExtraDofs;

	/*convert all dofs to system dofs*/
	m_x = convertAllDofsToSystemDofs(m_X);
	

	Eigen::VectorXd allDofs = convertSystemDofsToAllDofs(m_x);
	spdlog::info("Converted system dofs vs. m_X is {}", (allDofs - m_X).norm());


	m_para_X = m_x;
	m_V.setZero(m_x.size());			 //!!!!!!now this velocity is zero, TODO: parameters

	m_v = Eigen::VectorXd(m_x.size()); m_v.setZero();
	m_appliedForces = Eigen::VectorXd(m_x.size()); m_appliedForces.setZero();
}

Eigen::VectorXd SOReflectionBendingModeShellStructure::convertAllDofsToSystemDofs(const Eigen::VectorXd& allDofs)
{
	Eigen::VectorXd systemDofs;
	systemDofs.setZero(m_numSystemDofs);//plus translations

	for (int i = 0; i < m_nv; i++)
	{
		Eigen::Vector3d vi = allDofs.segment<3>(3 * i);
		int systemDofStartIndex = m_realVertex_to_systemDofAndNumDofs.at(i).first;
		int numSystemDof = m_realVertex_to_systemDofAndNumDofs.at(i).second;
		if (numSystemDof == 1)
		{//corners, only use z
			systemDofs[systemDofStartIndex] = vi[2];
		}
		else if (numSystemDof == 2)
		{//left, right, top, bottom boundaries
			auto type = boundaryVertices.at(i);
			if (type == BoundaryVertexType::left || type == BoundaryVertexType::right)
			{//we will use y and z
				systemDofs.segment<2>(systemDofStartIndex) = Eigen::Vector2d(vi[1], vi[2]);
			}
			else
			{//top and bottom, we will use x and z
				systemDofs.segment<2>(systemDofStartIndex) = Eigen::Vector2d(vi[0], vi[2]);
			}
		}
		else
		{//3
			systemDofs.segment<3>(systemDofStartIndex) = vi;
		}
	}
	//we use corners to get last four dofs x1,y1,x2,y2
	Eigen::Vector3d topleftVertex = allDofs.segment<3>(m_fourCornerIndices[0]*3);//topleft
	Eigen::Vector3d bottomRightVertex = allDofs.segment<3>(m_fourCornerIndices[3]*3);//bottomright
	systemDofs.tail(4) = Eigen::Vector4d(topleftVertex[0], topleftVertex[1], bottomRightVertex[0], bottomRightVertex[1]);

	return systemDofs;
}

Eigen::VectorXd SOReflectionBendingModeShellStructure::convertSystemDofsToAllDofs(const Eigen::VectorXd& systemDofs)
{
	Eigen::VectorXd allDofs(m_nv*3);
	Eigen::Vector4d fourReflectionPlane = systemDofs.tail(4);
	double x1 = fourReflectionPlane[0];
	double y1 = fourReflectionPlane[1];
	double x2 = fourReflectionPlane[2];
	double y2 = fourReflectionPlane[3];

	for (const auto& system_to_realVertex : m_systemDof_to_realVertexAndNumDofs)
	{
		int systemDofStartIndex = system_to_realVertex.first;
		int realVertexIndex = system_to_realVertex.second.first;
		int numSystemDof = system_to_realVertex.second.second;
		if (numSystemDof == 1)
		{//corners, only use z
			if (m_fourCornerIndices[0] == realVertexIndex)//topleft
				allDofs.segment<3>(realVertexIndex * 3) = Eigen::Vector3d(x1, y1, systemDofs[systemDofStartIndex]);
			else if(m_fourCornerIndices[1] == realVertexIndex)//topright
				allDofs.segment<3>(realVertexIndex * 3) = Eigen::Vector3d(x2, y1, systemDofs[systemDofStartIndex]);
			else if(m_fourCornerIndices[2] == realVertexIndex)//botttomleft
				allDofs.segment<3>(realVertexIndex * 3) = Eigen::Vector3d(x1, y2, systemDofs[systemDofStartIndex]);
			else//(m_fourCornerIndices[3] == realVertexIndex)
				allDofs.segment<3>(realVertexIndex * 3) = Eigen::Vector3d(x2, y2, systemDofs[systemDofStartIndex]);
		}
		else if (numSystemDof == 2)
		{//left, right, top, bottom boundaries
			auto type = boundaryVertices.at(realVertexIndex);
			if (type == BoundaryVertexType::left)
			{//we will use y and z
				Eigen::Vector2d yz = systemDofs.segment<2>(systemDofStartIndex);
				allDofs.segment<3>(realVertexIndex * 3) = Eigen::Vector3d(x1, yz[0], yz[1]);
			}
			else if (type == BoundaryVertexType::right)
			{
				Eigen::Vector2d yz = systemDofs.segment<2>(systemDofStartIndex);
				allDofs.segment<3>(realVertexIndex * 3) = Eigen::Vector3d(x2, yz[0], yz[1]);
			}
			else if(type == BoundaryVertexType::top)
			{//top and bottom, we will use x and z
				Eigen::Vector2d xz = systemDofs.segment<2>(systemDofStartIndex);
				allDofs.segment<3>(realVertexIndex * 3) = Eigen::Vector3d(xz[0], y1, xz[1]);
			}
			else
			{//bottom
				Eigen::Vector2d xz = systemDofs.segment<2>(systemDofStartIndex);
				allDofs.segment<3>(realVertexIndex * 3) = Eigen::Vector3d(xz[0], y2, xz[1]);
			}
		}
		else
		{//3
			allDofs.segment<3>(realVertexIndex * 3) = systemDofs.segment<3>(systemDofStartIndex);
		}
	}

	return allDofs;
}

SpMat SOReflectionBendingModeShellStructure::convertSystemDofsToAllDofsJacobian(const Eigen::VectorXd& systemDofs)
{
	SpMat jacobian(3 * m_nv, systemDofs.size());
	TripVec jacobian_triplets;

	//Eigen::VectorXd allDofs(m_nv * 3);
	//Eigen::Vector4d fourReflectionPlane = systemDofs.tail(4);
	//double x1 = fourReflectionPlane[0]; 
	//double y1 = fourReflectionPlane[1]; 
	//double x2 = fourReflectionPlane[2]; 
	//double y2 = fourReflectionPlane[3]; 
	int x1_dofIndex = m_numSystemDofs - 4;
	int y1_dofIndex = m_numSystemDofs - 3;
	int x2_dofIndex = m_numSystemDofs - 2;
	int y2_dofIndex = m_numSystemDofs - 1;

	for (const auto& system_to_realVertex : m_systemDof_to_realVertexAndNumDofs)
	{
		int systemDofStartIndex = system_to_realVertex.first;
		int realVertexIndex = system_to_realVertex.second.first;
		int numSystemDof = system_to_realVertex.second.second;
		if (numSystemDof == 1)
		{//corners, only use z
			if (m_fourCornerIndices[0] == realVertexIndex)//topleft
			{
				//allDofs.segment<3>(realVertexIndex * 3) = Eigen::Vector3d(x1, y1, systemDofs[systemDofStartIndex]);
				jacobian_triplets.emplace_back(realVertexIndex * 3 + 0, x1_dofIndex, 1.0);
				jacobian_triplets.emplace_back(realVertexIndex * 3 + 1, y1_dofIndex, 1.0);
			}
			else if (m_fourCornerIndices[1] == realVertexIndex)//topright
			{
				//allDofs.segment<3>(realVertexIndex * 3) = Eigen::Vector3d(x2, y1, systemDofs[systemDofStartIndex]);
				jacobian_triplets.emplace_back(realVertexIndex * 3 + 0, x2_dofIndex, 1.0);
				jacobian_triplets.emplace_back(realVertexIndex * 3 + 1, y1_dofIndex, 1.0);
			}
			else if (m_fourCornerIndices[2] == realVertexIndex)//botttomleft
			{
				//allDofs.segment<3>(realVertexIndex * 3) = Eigen::Vector3d(x1, y2, systemDofs[systemDofStartIndex]);
				jacobian_triplets.emplace_back(realVertexIndex * 3 + 0, x1_dofIndex, 1.0);
				jacobian_triplets.emplace_back(realVertexIndex * 3 + 1, y2_dofIndex, 1.0);
			}
			else//(m_fourCornerIndices[3] == realVertexIndex)
			{
				//allDofs.segment<3>(realVertexIndex * 3) = Eigen::Vector3d(x2, y2, systemDofs[systemDofStartIndex]); 
				jacobian_triplets.emplace_back(realVertexIndex * 3 + 0, x2_dofIndex, 1.0);
				jacobian_triplets.emplace_back(realVertexIndex * 3 + 1, y2_dofIndex, 1.0);
			}
			jacobian_triplets.emplace_back(realVertexIndex * 3 + 2, systemDofStartIndex, 1.0);
		}
		else if (numSystemDof == 2)
		{//left, right, top, bottom boundaries
			auto type = boundaryVertices.at(realVertexIndex);
			if (type == BoundaryVertexType::left)
			{//we will use y and z
				//Eigen::Vector2d yz = systemDofs.segment<2>(systemDofStartIndex);
				//allDofs.segment<3>(realVertexIndex * 3) = Eigen::Vector3d(x1, yz[0], yz[1]);

				jacobian_triplets.emplace_back(realVertexIndex * 3 + 0, x1_dofIndex, 1.0);
				jacobian_triplets.emplace_back(realVertexIndex * 3 + 1, systemDofStartIndex, 1.0);
				jacobian_triplets.emplace_back(realVertexIndex * 3 + 2, systemDofStartIndex+1, 1.0);
			}
			else if (type == BoundaryVertexType::right)
			{
				//Eigen::Vector2d yz = systemDofs.segment<2>(systemDofStartIndex);
				//allDofs.segment<3>(realVertexIndex * 3) = Eigen::Vector3d(x2, yz[0], yz[1]);

				jacobian_triplets.emplace_back(realVertexIndex * 3 + 0, x2_dofIndex, 1.0);
				jacobian_triplets.emplace_back(realVertexIndex * 3 + 1, systemDofStartIndex, 1.0);
				jacobian_triplets.emplace_back(realVertexIndex * 3 + 2, systemDofStartIndex + 1, 1.0);
			}
			else if (type == BoundaryVertexType::top)
			{//top and bottom, we will use x and z
				//Eigen::Vector2d xz = systemDofs.segment<2>(systemDofStartIndex);
				//allDofs.segment<3>(realVertexIndex * 3) = Eigen::Vector3d(xz[0], y1, xz[1]);

				jacobian_triplets.emplace_back(realVertexIndex * 3 + 0, systemDofStartIndex, 1.0);
				jacobian_triplets.emplace_back(realVertexIndex * 3 + 1, y1_dofIndex, 1.0);
				jacobian_triplets.emplace_back(realVertexIndex * 3 + 2, systemDofStartIndex + 1, 1.0);
			}
			else
			{//bottom
				//Eigen::Vector2d xz = systemDofs.segment<2>(systemDofStartIndex);
				//allDofs.segment<3>(realVertexIndex * 3) = Eigen::Vector3d(xz[0], y2, xz[1]);

				jacobian_triplets.emplace_back(realVertexIndex * 3 + 0, systemDofStartIndex, 1.0);
				jacobian_triplets.emplace_back(realVertexIndex * 3 + 1, y2_dofIndex, 1.0);
				jacobian_triplets.emplace_back(realVertexIndex * 3 + 2, systemDofStartIndex + 1, 1.0);
			}
		}
		else
		{//3
			//allDofs.segment<3>(realVertexIndex * 3) = systemDofs.segment<3>(systemDofStartIndex);

			jacobian_triplets.emplace_back(realVertexIndex * 3 + 0, systemDofStartIndex, 1.0);
			jacobian_triplets.emplace_back(realVertexIndex * 3 + 1, systemDofStartIndex + 1, 1.0);
			jacobian_triplets.emplace_back(realVertexIndex * 3 + 2, systemDofStartIndex + 2, 1.0);
		}
	}

	jacobian.setFromTriplets(jacobian_triplets.begin(), jacobian_triplets.end());

	return jacobian;
}


void SOReflectionBendingModeShellStructure::initElements(const Eigen::VectorXd& vX)
{
	Eigen::VectorXd allDofs = convertSystemDofsToAllDofs(vX);

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

			elem.init(dofIndices, dofIndices, allDofs, m_material);
			elem.setThickness(m_cstThickness);
			elem.setElementProperty(InternalElement);

		}
	}


	{
		m_hingeElements.clear();
		m_reflectionDSHingeFixedRestPoseElement.clear();

		//Bending elements for bending
		//m_triMesh->buildCoHingeStructure();
		m_triMesh->buildHingeStructure();
		auto hinges = m_triMesh->hinges();
		for (int i = 0; i < hinges.size(); i++)
		{
			auto hinge = hinges[i];
			if ((hinge.tris[0] == -1) || (hinge.tris[1] == -1))
			{//boundary edges
				std::vector<int> inds(3);
				inds[0] = hinge.edge[0];
				inds[1] = hinge.edge[1];
				inds[2] = hinge.flaps[0] == -1 ? hinge.flaps[1] : hinge.flaps[0];
				//std::cout << "build reflection hinge with indices " << inds[0] << ", " << inds[1] << ", " << inds[2] << std::endl;

				m_reflectionDSHingeFixedRestPoseElement.push_back(reflectionDSHingeFixedRestPoseElement());
				reflectionDSHingeFixedRestPoseElement& elem = m_reflectionDSHingeFixedRestPoseElement.back();

				std::vector<int> dofIndices;
				for (int dof_i = 0; dof_i < inds.size(); dof_i++)
				{
					for (int j = 0; j < dimension; j++)
						dofIndices.push_back(inds[dof_i] * dimension + j);
				}

				auto e0_type = boundaryVertices.at(inds[0]);
				auto e1_type = boundaryVertices.at(inds[1]);
				bool useReflectionXPlane;
				if ((e0_type == BoundaryVertexType::left || e1_type == BoundaryVertexType::left) ||
					(e0_type == BoundaryVertexType::right || e1_type == BoundaryVertexType::right))
					useReflectionXPlane = true;
				else
					useReflectionXPlane = false;

				elem.init(dofIndices, dofIndices, allDofs, m_material, m_cstThickness, useReflectionXPlane);
				//elem.setThickness(m_cstThickness);
				elem.setElementProperty(InternalElement);
			}
			else
			{
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

				elem.init(dofIndices, dofIndices, allDofs, m_material, m_cstThickness);// 
				//elem.setThickness(m_cstThickness);
				elem.setElementProperty(InternalElement);

				//spdlog::info("Theta of hinge is {}", elem.getCurrentAngle(vX));
			}
		}

	}

	reconstructElementVector();

	//printf("init %d with shell: %d  bending: %d  motor: %d\n", m_elements.size(), m_cstShellElements.size(), m_dsBendingElements.size(), m_dihedralElements.size());
}

void SOReflectionBendingModeShellStructure::reconstructElementVector()
{
	m_elements.clear();

	for (int i = 0; i < (int)m_cstShellElements.size(); i++)
		m_elements.push_back(&m_cstShellElements[i]);
	for (int i = 0; i < (int)m_hingeElements.size(); i++)
		m_elements.push_back(&m_hingeElements[i]);
	for (int i = 0; i < (int)m_reflectionDSHingeFixedRestPoseElement.size(); i++)
		m_elements.push_back(&m_reflectionDSHingeFixedRestPoseElement[i]);
	
	for (int i = 0; i < (int)m_fixDoFElements.size(); i++)
		m_elements.push_back(&m_fixDoFElements[i]);
	
}


void SOReflectionBendingModeShellStructure::computeMasses(const Eigen::VectorXd& vX)
{
	//const P3DIndexArray* pIB = m_triMesh->GetIndexArray();
	int nTriangles = m_triMesh->numTriangles();
	m_mass = Eigen::VectorXd(vX.size()); m_mass.setZero();

	Eigen::VectorXd allDofs = convertSystemDofsToAllDofs(vX);
	for (int i = 0; i < nTriangles; i++)
	{
		V3D p[3];
		for (int j = 0; j < 3; j++)
			p[j] = get3D(m_triMesh->faces()[3 * i + j], allDofs);
		double area = 0.5 * crossProd3D((p[1] - p[0]), (p[2] - p[0])).norm();
		//A += area;
		//printf("%d %d %d %.15f\n", m_triMesh->faces()[3 * i + 0], m_triMesh->faces()[3 * i + 1], m_triMesh->faces()[3 * i + 2], area);
		for (int j = 0; j < 3; j++)
		{
			int realVertexIndex = m_triMesh->faces()[3 * i + j];
			bool isBoundaryVertex = false;
			auto finder = boundaryVertices.find(realVertexIndex);
			if (finder != boundaryVertices.end())
				isBoundaryVertex = true;

			auto systemDofAndNumDofs = m_realVertex_to_systemDofAndNumDofs.at(realVertexIndex);
			int sysDofStartIndex = systemDofAndNumDofs.first;
			int numSystemDofs = systemDofAndNumDofs.second;

			double massMultiplier = 1.0;
			if (numSystemDofs == 1)
				massMultiplier = 4.0;//corners mass need multiply 4
			else if (numSystemDofs == 2)//other boundaries need multiply 2
				massMultiplier = 2.0;
			else
				massMultiplier = 1.0;//inner multiply 1
			

			for (int k = 0; k < numSystemDofs; k++)
				m_mass[sysDofStartIndex + k] += massMultiplier * 1.0 / 3.0 * area * m_cstThickness * m_material->rho;
		}
	}

	spdlog::info("Overall mass is {}", m_mass.sum() / 3.0);

}

void SOReflectionBendingModeShellStructure::reconstructMassStuff(const Eigen::VectorXd& vX)
{
	computeMasses(vX);
	m_gravity = Eigen::VectorXd(m_x.size()); m_gravity.setZero();
	if (USE_GRAVITY)
	{
		for (int i = 0; i < (m_gravity.size() - m_numExtraDofs) / 3; i++)
		{
			m_gravity[3 * i + 1] = -m_g * m_mass[3 * i + 1];//gravity only set to vertices at y direction
		}
	}
}

void SOReflectionBendingModeShellStructure::computeHessian(const Eigen::VectorXd& vx, const Eigen::VectorXd& vX, SpMat& A)
{
	A = SpMat(m_x.size(), m_x.size());
	A.setZero();


	Eigen::VectorXd allDofs_X = convertSystemDofsToAllDofs(vX);
	Eigen::VectorXd allDofs = convertSystemDofsToAllDofs(vx);
	SpMat jacobian = convertSystemDofsToAllDofsJacobian(vx);

	SpMat hessian(allDofs.size(), allDofs.size());
	//Column Major Matrix
	TripVec triplets;

	//add triplets for all edges
	for (int i = 0; i < m_elements.size(); i++)
	{
		m_elements[i]->addTriplets(allDofs, allDofs_X, triplets);
	}
	hessian.setFromTriplets(triplets.begin(), triplets.end());
	//Convert to CCS
	hessian.makeCompressed();

	A = jacobian.transpose() * hessian * jacobian;
}

bool SOReflectionBendingModeShellStructure::staticSolve_newton(Eigen::VectorXd& vx, const Eigen::VectorXd& vX)
{

	auto stepBeginningFunction = [this](const Eigen::VectorXd& vx, int currentIteration)
	{
		ipc::Constraints constraint_set;

		Eigen::VectorXd allDofs = convertSystemDofsToAllDofs(vx);
		Eigen::MatrixXd V = (Eigen::Map<const Eigen::MatrixXd>(allDofs.data(), dimension, allDofs.size() / dimension)).transpose();
		constraint_set.build(m_ipc_collisionMesh, V, m_ipc_dhat);

		//in the gradient computation, we will update barrier stiffness
		//since objective and hessian are computed afterward in the same step
		double min_distance = ipc::compute_minimum_distance(m_ipc_collisionMesh, V, constraint_set);
		spdlog::debug("Minimum distance is {} with number of constraints {}", min_distance, constraint_set.size());

		double bbox_diagonal_length = ipc::world_bbox_diagonal_length(V);
		m_ipc_adaptiveStiffness = ipc::update_barrier_stiffness(m_ipc_minDistance, min_distance,
			m_ipc_upperBoundStiffness, m_ipc_adaptiveStiffness, bbox_diagonal_length, 1e-4);
		spdlog::debug("Update barrier stiffness to {} with upper bound {}", m_ipc_adaptiveStiffness, m_ipc_upperBoundStiffness);
		m_ipc_minDistance = min_distance;

	};

	auto computeStaticObjective = [this](const Eigen::VectorXd& vx)
	{
		Eigen::VectorXd allDofs = convertSystemDofsToAllDofs(vx);
		//internal energy
		tbb::enumerable_thread_specific<double> storage(0);
		tbb::parallel_for(
			tbb::blocked_range<size_t>(size_t(0), m_elements.size()),
			[&](const tbb::blocked_range<size_t>& r) {
				auto& local_Eint = storage.local();
				for (size_t i = r.begin(); i < r.end(); i++)
					//for (size_t i = 0; i < m_elements.size(); i++)
				{
					local_Eint += m_elements[i]->computeEnergy(allDofs, m_X);
					//elementEnergy[i] = m_elements[i]->computeEnergy(vx, vX);
				}
			});
		double Eint = 0.0;
		for (const auto& local_Eint : storage) {
			Eint += local_Eint;
		}

		//global contact barrier
		double EselfContact = 0.0;
		{
			ipc::Constraints constraint_set;

			Eigen::VectorXd allDofs = convertSystemDofsToAllDofs(vx);
			Eigen::MatrixXd V = (Eigen::Map<const Eigen::MatrixXd>(allDofs.data(), dimension, allDofs.size() / dimension)).transpose();
			constraint_set.build(m_ipc_collisionMesh, V, m_ipc_dhat);
			double barrierPotential = ipc::compute_barrier_potential(m_ipc_collisionMesh, V, constraint_set, m_ipc_dhat);

			EselfContact += m_ipc_adaptiveStiffness * barrierPotential;
		}

		double E = Eint + EselfContact;

		return E;

	};

	auto computeStaticGradient = [this](const Eigen::VectorXd& vx, Eigen::VectorXd& grad)
	{
		grad = Eigen::VectorXd(vx.size());
		grad.setZero();

		Eigen::VectorXd allDofs = convertSystemDofsToAllDofs(vx);
		SpMat allDofs_jacobian = convertSystemDofsToAllDofsJacobian(vx);

		//compute elements gradient
		tbb::enumerable_thread_specific<Eigen::VectorXd> storage(
			Eigen::VectorXd::Zero(allDofs.size()));

		tbb::parallel_for(
			tbb::blocked_range<size_t>(size_t(0), m_elements.size()),
			[&](const tbb::blocked_range<size_t>& r) {
				auto& local_grad = storage.local();
				for (size_t i = r.begin(); i < r.end(); i++)
		//		for (size_t i = 0; i < m_elements.size(); i++)
				{
					Eigen::VectorXd element_grad;
					m_elements[i]->computeGradient(allDofs, m_X, element_grad);

					const auto element_dofs = m_elements[i]->getDofIndices();
					for (int dof_i = 0; dof_i < element_dofs.size(); dof_i++)
					{
						local_grad[element_dofs[dof_i]] += element_grad[dof_i];
					}
					if(isnan(element_grad.norm()))
						spdlog::debug("grad element {} norm {}", i, element_grad.norm());
				}
			});
		Eigen::VectorXd grad_int;
		grad_int.setZero(allDofs.size());
		for (const auto& local_grad : storage) {
			grad_int += local_grad;
		}
		grad = allDofs_jacobian.transpose() * grad_int;

		//Eigen::VectorXd strainLimitingGrad(vx.size());
		//computeAllStrainLimitingPotentialGradient(vx, vX, strainLimitingGrad);
		Eigen::VectorXd EselfContactGrad(vx.size()); EselfContactGrad.setZero();
		{
			ipc::Constraints constraint_set;
			Eigen::MatrixXd V = (Eigen::Map<const Eigen::MatrixXd>(allDofs.data(), dimension, allDofs.size() / dimension)).transpose();
			constraint_set.build(m_ipc_collisionMesh, V, m_ipc_dhat);
			Eigen::VectorXd barrierPotentialGradient = ipc::compute_barrier_potential_gradient(m_ipc_collisionMesh, V, constraint_set, m_ipc_dhat);

			EselfContactGrad += m_ipc_adaptiveStiffness * allDofs_jacobian.transpose() * barrierPotentialGradient;
		}


		grad += EselfContactGrad;

	};

	auto computeStaticHessian = [this](const Eigen::VectorXd& vx, SpMat& Hessian)
	{
		Hessian = SpMat(vx.size(), vx.size());
		Hessian.setZero();

		Eigen::VectorXd allDofs = convertSystemDofsToAllDofs(vx);
		SpMat allDofs_jacobian = convertSystemDofsToAllDofsJacobian(vx);

		tbb::enumerable_thread_specific<TripVec> storage;
		tbb::parallel_for(
			tbb::blocked_range<size_t>(size_t(0), m_elements.size()),
			[&](const tbb::blocked_range<size_t>& r) {
				auto& local_hess_triplets = storage.local();

				for (size_t i = r.begin(); i < r.end(); i++)
				{
					Eigen::MatrixXd element_hess;
					m_elements[i]->computeHessian(allDofs, m_X, element_hess);

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
		SpMat hessian_int(allDofs.size(), allDofs.size());
		for (const auto& local_hess_triplets : storage) {
			Eigen::SparseMatrix<double> local_hess(allDofs.size(), allDofs.size());
			local_hess.setFromTriplets(
				local_hess_triplets.begin(), local_hess_triplets.end());
			hessian_int += local_hess;
		}
		Hessian = allDofs_jacobian.transpose() * hessian_int * allDofs_jacobian;


		SpMat EselfContactHessian(vx.size(), vx.size());
		{
			ipc::Constraints constraint_set;
			Eigen::VectorXd allDofs = convertSystemDofsToAllDofs(vx);
			SpMat allDofs_jacobian = convertSystemDofsToAllDofsJacobian(vx);
			Eigen::MatrixXd V = (Eigen::Map<Eigen::MatrixXd>(allDofs.data(), dimension, allDofs.size() / dimension)).transpose();
			constraint_set.build(m_ipc_collisionMesh, V, m_ipc_dhat);
			SpMat barrierPotentialHessian = ipc::compute_barrier_potential_hessian(m_ipc_collisionMesh, V, constraint_set, m_ipc_dhat, false);

			EselfContactHessian += m_ipc_adaptiveStiffness * allDofs_jacobian.transpose() * barrierPotentialHessian * allDofs_jacobian;
		}

		Hessian += EselfContactHessian;
	};

	auto computeLargestStepSize = [this](const Eigen::VectorXd& vx, const Eigen::VectorXd& dx)
	{

		spdlog::debug("Try to find a largest step size");


		double largestStepSize = 1.0;
		{
			//we use CCD to find the largest step size
			Eigen::VectorXd allDofs0 = convertSystemDofsToAllDofs(vx);
			Eigen::MatrixXd V0 = (Eigen::Map<const Eigen::MatrixXd>(allDofs0.data(), dimension, allDofs0.size() / dimension)).transpose();


			Eigen::VectorXd allDofs1 = convertSystemDofsToAllDofs((vx + largestStepSize * dx).eval());
			Eigen::MatrixXd V1 = (Eigen::Map<const Eigen::MatrixXd>(allDofs1.data(), dimension, allDofs1.size() / dimension)).transpose();

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
				allDofs1 = convertSystemDofsToAllDofs((vx + largestStepSize * dx).eval());
				V1 = (Eigen::Map<const Eigen::MatrixXd>(allDofs1.data(), dimension, allDofs1.size() / dimension)).transpose();
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

	/*{
		SimOpt::DerivativeTester tester;
		Eigen::VectorXd grad;
		computeStaticGradient(vx, grad);
		tester.testGradient(vx, grad, computeStaticObjective,10);
		SpMat hessian;
		computeStaticHessian(vx, hessian);
		tester.testJacobian(vx, hessian, computeStaticGradient, 10);
		exit(0);
	}*/

	std::unique_ptr<SimOpt::UnconstrainedNewtonLinearSolver> linearSolver(new SimOpt::UnconstrainedNewtonLinearSolverLLT(computeStaticHessian));
	SimOpt::UnconstrainedNewton solver(std::move(linearSolver));
	solver.setX(vx);
	solver.setObjective(computeStaticObjective);
	solver.setObjectiveGradient(computeStaticGradient);
	solver.setLargestLineSearchStepFunctor(computeLargestStepSize);
	solver.setStepBeginningFunctor(stepBeginningFunction);
	solver.setMaxIter(2000);//m_newtonMaxIterations
	//solver.setLinesearchMaxIter(m_newtonMaxLinesearchSteps);
	solver.setGradientThreshold(1e-6);//m_newtonSolveResidual
	solver.solve();

	if (solver.state() != SimOpt::UnconstrainedNewton::SOLVED)
	{
		printf("Static solve didn't converge to %.10f\n", m_newtonSolveResidual);//, the rhs is %.10f , solver.getLastGradientNorm()
		return false;
	}

	//after taking a step, we update states
	m_x = solver.getX();
	vx = m_x;
	return true;
}

void SOReflectionBendingModeShellStructure::projectToManifold()
{
	bool use_gravity = USE_GRAVITY;
	USE_GRAVITY = false;
	solveStaticState();
	USE_GRAVITY = use_gravity;

	//convert from system dofs to all dofs
	//setVisualizationConfiguration(m_x);
	V3D centerAllDofs(0, 0, 0);
	for (int i = 0; i < visualize_m_x.size() / 3; i++)
	{
		centerAllDofs += visualize_m_x.segment<3>(i * 3);
	}
	centerAllDofs /= (double)(visualize_m_x.size() / 3);


	Eigen::Vector4d fourReflectionPlane = m_x.tail(4);
	double x1 = fourReflectionPlane[0];
	double y1 = fourReflectionPlane[1];
	double x2 = fourReflectionPlane[2];
	double y2 = fourReflectionPlane[3];
	//move to center
	V3D center((x1 + x2) / 2.0, (y1 + y2) / 2.0, centerAllDofs[2]);

	for (int i = 0; i < visualize_m_x.size() / 3; i++)
		visualize_m_x.segment<3>(i * 3) -= center;
	SOBendingModeShellStructure::setVisualizationConfiguration(visualize_m_x);


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

void SOReflectionBendingModeShellStructure::convertToStrainSpaceMode(const Eigen::VectorXd& positionSpaceMode, Eigen::VectorXd& strainSpaceMode)
{
	Eigen::VectorXd pm = positionSpaceMode;// .normalized();
	strainSpaceMode = Eigen::VectorXd(m_hingeElements.size() + m_reflectionDSHingeFixedRestPoseElement.size()); strainSpaceMode.setZero();


	SpMat all_dthetadallx(m_hingeElements.size() + m_reflectionDSHingeFixedRestPoseElement.size(), 3*m_nv);
	TripVec all_dthetadallx_triplets;
	for (int i = 0; i < m_hingeElements.size(); i++)
	{
		Eigen::VectorXd dthetadx(16);
		dthetadx.setZero();
		m_hingeElements[i].computedThetadx(m_X, m_X, dthetadx);
		//std::cout << "hinge " << i << " dthetadx norm " << dthetadx.norm() << std::endl;
		const auto& inds = m_hingeElements[i].getDofIndices();

		for (int dof_i = 0; dof_i < inds.size(); dof_i++)
		{
			all_dthetadallx_triplets.emplace_back(i, inds[dof_i], dthetadx[dof_i]);
		}
	}
	for (int i = 0; i < m_reflectionDSHingeFixedRestPoseElement.size(); i++)
	{
		Eigen::VectorXd dthetadx(9);
		dthetadx.setZero();
		m_reflectionDSHingeFixedRestPoseElement[i].computedThetadx(m_X, m_X, dthetadx);
		//std::cout << "hinge " << i << " dthetadx norm " << dthetadx.norm() << std::endl;
		const auto& inds = m_reflectionDSHingeFixedRestPoseElement[i].getDofIndices();

		for (int dof_i = 0; dof_i < inds.size(); dof_i++)
		{
			all_dthetadallx_triplets.emplace_back(m_hingeElements.size() + i, inds[dof_i], dthetadx[dof_i]);
		}
	}

	all_dthetadallx.setFromTriplets(all_dthetadallx_triplets.begin(), all_dthetadallx_triplets.end());

	SpMat jacobian = convertSystemDofsToAllDofsJacobian(m_para_X);

	strainSpaceMode = (all_dthetadallx * jacobian * pm).normalized();
}

void SOReflectionBendingModeShellStructure::animateStrainSpaceMode()
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
	for (int i = 0; i < m_reflectionDSHingeFixedRestPoseElement.size(); i++)
	{
		//std::cout << "hinge " << i << m_hingeElements[i].getRestAngle();
		double theta = m_reflectionDSHingeFixedRestPoseElement[i].getCurrentAngle(m_x);
		theta += ssm[m_hingeElements.size() + i];
		m_reflectionDSHingeFixedRestPoseElement[i].setRestAngle(theta);
		//std::cout << " update theta0 as " << theta << " ssm " << ssm[i] << std::endl;
	}

	projectToManifold();

	time_passed += m_dt;
}

void SOReflectionBendingModeShellStructure::computeStrainSpaceModes()
{
	computeBatchStrainSpaceMode();

	if (m_computeModeLength <= 0.0)
	{
		spdlog::warn("The input length should larger than 0.0");
		return;
	}
	Eigen::VectorXd initialGuess = m_para_X;
	double initialLength = 0.0;
	if (m_strainSpaceModes_index_to_lengthState.find(m_modesAnimationID) != m_strainSpaceModes_index_to_lengthState.end())
	{
		if (m_strainSpaceModes_index_to_lengthState[m_modesAnimationID].size() > 0)
		{
			initialLength = m_strainSpaceModes_index_to_lengthState[m_modesAnimationID].back().first;
			initialGuess = m_strainSpaceModes_index_to_lengthState[m_modesAnimationID].back().second;
		}
	}

	//double maxStrainSpaceCoeff = m_strainSpaceModes[m_modesAnimationID].maxCoeff();
	//m_computeModeStepLength = 1.0 / 20.0 / maxStrainSpaceCoeff;
	//spdlog::info("Step length is {}", m_computeModeStepLength);

	int numSteps = (m_computeModeLength - initialLength) / m_computeModeStepLength;
	for (int i = 1; i <= numSteps; i++)
	{
		double length = initialLength + i * m_computeModeStepLength;// std::min(initialLength + i * m_computeModeStepLength, m_computeModeLength);
		spdlog::info("Compute mode {} with strain space length {}", m_modesAnimationID, length);
		if (m_useStrainSpaceModes)
		{
			if (m_strainSpaceModes.size() >= m_modesAnimationID && m_strainSpaceModes[m_modesAnimationID].size() > 0)
			{
				Eigen::VectorXd ssm = m_strainSpaceModes[m_modesAnimationID] * length;
				Eigen::VectorXd initialGuess_allDofs = convertSystemDofsToAllDofs(initialGuess);
				for (int i = 0; i < m_hingeElements.size(); i++)
				{
					double theta = m_hingeElements[i].getCurrentAngle(initialGuess_allDofs);
					theta += ssm[i];
					m_hingeElements[i].setRestAngle(theta);
					//std::cout << "hinge " << i << " restAngle " << m_hingeElements[i].getRestAngle() << std::endl;
				}
				for (int i = 0; i < m_reflectionDSHingeFixedRestPoseElement.size(); i++)
				{
					double theta = m_reflectionDSHingeFixedRestPoseElement[i].getCurrentAngle(initialGuess_allDofs);
					theta += ssm[m_hingeElements.size() + i];
					m_reflectionDSHingeFixedRestPoseElement[i].setRestAngle(theta);
					//std::cout << "hinge " << i << " restAngle " << m_hingeElements[i].getRestAngle() << std::endl;
				}

				m_x = initialGuess;
				projectToManifold();

				/*double stateDiff = (visualize_m_x - m_X).norm();
				if (i == numSteps && stateDiff < 1.0)
				{
					numSteps++;
					if (1.0 / stateDiff > 50)
					{
						m_computeModeStepLength *= int(1.0 / stateDiff);
						spdlog::info("The state diff {}, so update length to {}", stateDiff, m_computeModeStepLength);
					}
				}*/
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


		m_strainSpaceModes_index_to_lengthState[m_modesAnimationID].emplace_back(length, m_x);
		initialGuess = m_x;
	}


	//compute energy
	auto computeEnergy = [this](const Eigen::VectorXd& vx)
	{
		Eigen::VectorXd allDofs = convertSystemDofsToAllDofs(vx);
		//internal energy
		tbb::enumerable_thread_specific<double> storage(0);
		tbb::parallel_for(
			tbb::blocked_range<size_t>(size_t(0), m_elements.size()),
			[&](const tbb::blocked_range<size_t>& r) {
				auto& local_Eint = storage.local();
				for (size_t i = r.begin(); i < r.end(); i++)
					//for (size_t i = 0; i < m_elements.size(); i++)
				{
					local_Eint += m_elements[i]->computeEnergy(allDofs, m_X);
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
	for (int i = 0; i < m_reflectionDSHingeFixedRestPoseElement.size(); i++)
	{
		m_reflectionDSHingeFixedRestPoseElement[i].setRestAngle(m_reflectionDSHingeFixedRestPoseElement[i].getInitialRestAngle());
	}
	std::cout << "Internal Energy along with displacement and modal displacement" << std::endl;
	int numStates = m_strainSpaceModes_index_to_lengthState[m_modesAnimationID].size();
	Eigen::VectorXd internalEnergy(numStates);
	for (size_t i = 0; i < numStates; i++)
	{
		internalEnergy[i] = computeEnergy(m_strainSpaceModes_index_to_lengthState[m_modesAnimationID][i].second);
		double displacement = (m_strainSpaceModes_index_to_lengthState[m_modesAnimationID][i].second - m_para_X).norm();
		std::cout << internalEnergy[i] << " " << displacement << " " << m_strainSpaceModes_index_to_lengthState[m_modesAnimationID][i].first << std::endl;
	}
}


void SOReflectionBendingModeShellStructure::ouputStrainSpaceModesAtIndex(int index)
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


	for (int i = 0; i < modeStates.size(); i++)
	{
		double modeLength = modeStates[i].first;
		const Eigen::VectorXd modeState = modeStates[i].second;

		string outputMeshName = m_outputPath + subFolderName + "/outputMesh_Mode_" + to_string(modeIndex) + "_" + to_string_with_precision(modeLength, 5) + ".obj";
		spdlog::info("Output simulation triangle mesh with name: {}", outputMeshName);

		Eigen::VectorXd modeState_allDofs = convertSystemDofsToAllDofs(modeState);
		m_triMesh->vertices() = modeState_allDofs;
		Eigen::MatrixXd outputTriangleMeshVertices;
		Eigen::MatrixXi outputTriangleMeshFaces;
		loadFromMesh(m_triMesh, outputTriangleMeshVertices, outputTriangleMeshFaces);
		outputTriangleMeshVertices.transposeInPlace();
		outputTriangleMeshFaces.transposeInPlace();

		igl::writeOBJ(outputMeshName, outputTriangleMeshVertices, outputTriangleMeshFaces);

	}
}


bool SOReflectionBendingModeShellStructure::loadAllStrainSpaceModes(const string& strainSpaceModesFolder)
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

		std::vector<std::pair<TriMesh, double>> meshSet_length;
		loadFolderMeshes(folderName[i], meshSet_length);

		for (int j = 0; j < meshSet_length.size(); j++)
		{
			Eigen::VectorXd allDofs = meshSet_length[j].first.vertices();
			Eigen::VectorXd systemDofs = convertAllDofsToSystemDofs(allDofs);
			m_strainSpaceModes_index_to_lengthState[modeIndex].emplace_back(meshSet_length[j].second, systemDofs);
		}
	}
}


void SOReflectionBendingModeShellStructure::computeReflectionMeshes(int startRows, int numRows, int startCols, int numCols, std::vector<Eigen::MatrixXd>& vertices_vector, std::vector<Eigen::MatrixXi>& faces_vector, std::vector<bool>& flipNorm_vector)
{

	Eigen::MatrixXi faces = Eigen::Map<Eigen::MatrixXi>(m_triMesh->faces().data(), 3, m_triMesh->faces().size() / 3);

	double x1 = visualize_m_x.segment<3>(m_fourCornerIndices[0] * 3)[0];
	double y1 = visualize_m_x.segment<3>(m_fourCornerIndices[0] * 3)[1];
	double x2 = visualize_m_x.segment<3>(m_fourCornerIndices[3] * 3)[0];
	double y2 = visualize_m_x.segment<3>(m_fourCornerIndices[3] * 3)[1];
	Eigen::Vector3d lefttoright(x2 - x1, 0, 0);
	Eigen::Vector3d toptobottom(0, y2 - y1, 0);

	//compute a 3by3 tiling
	std::map<Eigen::Vector2i, Eigen::MatrixXd, compareVector2<int>> tiling;
	for (int i = -1; i <= 1; i++)
	{
		for (int j = -1; j <= 1; j++)
		{
			if (i == 0 && j == 0)
			{
				Eigen::MatrixXd transVertices_xd = Eigen::Map<Eigen::MatrixXd>(visualize_m_x.data(), 3, visualize_m_x.size() / 3);
				tiling[Eigen::Vector2i(i, j)] = transVertices_xd;
				continue;
			}
			auto computeReflectionXPlane = [](double reflectionX, const Eigen::Vector3d& x)
			{
				double reflected_x = x[0] + 2.0 * (reflectionX - x[0]);
				return Eigen::Vector3d(reflected_x, x[1], x[2]);
			};
			auto computeReflectionYPlane = [](double reflectionY, const Eigen::Vector3d& x)
			{
				double reflected_y = x[1] + 2.0 * (reflectionY - x[1]);
				return Eigen::Vector3d(x[0], reflected_y, x[2]);
			};
			bool flipNorm = false;
			bool useReflectionXPlane = false;
			bool useReflectionYPlane = false;
			double reflectionPlaneX = 0.0;
			double reflectionPlaneY = 0.0;
			if (i == -1)
			{
				reflectionPlaneX = x1;
				useReflectionXPlane = true;
			}
			else if (i == 1)//
			{
				reflectionPlaneX = x2;
				useReflectionXPlane = true;
			}

			if (j == -1)
			{
				reflectionPlaneY = y1;
				useReflectionYPlane = true;
			}
			else if (j == 1)//
			{
				reflectionPlaneY = y2;
				useReflectionYPlane = true;
			}

			if ((abs(i) + abs(j)) % 2 == 1)
				flipNorm = true;

			Eigen::VectorXd transVertices = visualize_m_x;
			for (int v_i = 0; v_i < m_nv; v_i++)
			{
				if (useReflectionXPlane)
					transVertices.segment<3>(3 * v_i) = computeReflectionXPlane(reflectionPlaneX, transVertices.segment<3>(3 * v_i));
				if (useReflectionYPlane)
					transVertices.segment<3>(3 * v_i) = computeReflectionYPlane(reflectionPlaneY, transVertices.segment<3>(3 * v_i));
			}

			Eigen::MatrixXd transVertices_xd = Eigen::Map<Eigen::MatrixXd>(transVertices.data(), 3, transVertices.size() / 3);

			tiling[Eigen::Vector2i(i, j)] = transVertices_xd;


		}
	}

	for (int i = startRows; i <= startRows+numRows; i++)
	{
		for (int j = startCols; j <= startCols+numCols; j++)
		{
			//if (i == 0 && j == 0)
			//	continue;
			int reference_i = i;
			int reference_j = j;
			while (reference_i < -1)
				reference_i += 2;
			while (reference_i > 1)
				reference_i -= 2;
			while (reference_j < -1)
				reference_j += 2;
			while (reference_j > 1)
				reference_j -= 2;

			Eigen::MatrixXd vertices = tiling.at(Eigen::Vector2i(reference_i, reference_j));
			vertices.colwise() += ((i - reference_i) * lefttoright + (j - reference_j) * toptobottom);

			bool flipNorm = false;
			if ((abs(i) + abs(j)) % 2 == 1)
				flipNorm = true;

			vertices_vector.push_back(vertices);
			faces_vector.push_back(faces);
			flipNorm_vector.push_back(flipNorm);

		}
	}
}

void SOReflectionBendingModeShellStructure::renderImGUIViewerWindow()
{
	ImGui::Begin("Mesh Window");
	{
		if (ImGui::Button("Reset simulation"))
		{
			reset();
		}
		
		ImGui::Checkbox("Draw reflection patches", &m_drawReflectionPatches);
		if (ImGui::Button("Ouput reflection patches"))
		{
			int startRows = -2;
			int startCols = -2;
			int numRows = 4;
			int numCols = 4;
			std::vector<Eigen::MatrixXd> vertices_vector;
			std::vector<Eigen::MatrixXi> faces_vector;
			std::vector<bool> flipNorm_vector;
			computeReflectionMeshes(startRows, numRows, startCols, numCols, vertices_vector, faces_vector, flipNorm_vector);

			Eigen::MatrixXd combined_vertices;
			Eigen::MatrixXi combined_faces;
			TriMesh::combineMeshes(vertices_vector, faces_vector,
				combined_vertices, combined_faces);
			igl::writeOBJ(m_outputPath + "/reflectionMesh.obj", combined_vertices.transpose(), combined_faces.transpose());
		}

		ImGui::Checkbox("Draw wire frame", &m_drawWireframe);
		ImGui::Checkbox("Draw fixing vertices", &m_drawFixingVertices);
		ImGui::Checkbox("Use static solver", &m_useStaticSolver);
		ImGui::Checkbox("Draw both wire frame and mesh", &m_useBothWireframeAndmesh);
		ImGui::Checkbox("Draw max strain color", &m_isRenderStrainColor);
		
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
			ImGui::Checkbox("Use Strain Space Modes", &m_useStrainSpaceModes);/*
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
				if (fromModeIndex < 6)
					fromModeIndex = 6;
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
				else if(animateStrainSpaceModeIndex>=0)// && animateStrainSpaceModeIndex <= m_strainSpaceModes_index_to_lengthState[m_modesAnimationID].size() - 1)
				{
					int temp = animateStrainSpaceModeIndex >= m_strainSpaceModes_index_to_lengthState[m_modesAnimationID].size() ? m_strainSpaceModes_index_to_lengthState[m_modesAnimationID].size() - 1 : animateStrainSpaceModeIndex;
					ImGui::Text("The length of this mode is %f", m_strainSpaceModes_index_to_lengthState[m_modesAnimationID].at(temp).first);
					setVisualizationConfiguration(m_strainSpaceModes_index_to_lengthState[m_modesAnimationID].at(temp).second);
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

	ImGui::End();


	ImGui::Begin("Iteractive Pinning");
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

			//we remove all spring elements
			for (int i = m_pairAttachmentElements.size() - 1; i >= 0; i--)
			{
				if (m_elements[m_elements.size() - 1] != &m_pairAttachmentElements[i])
					spdlog::error("the element order has a problem");

				m_elements.erase(m_elements.begin() + m_elements.size() - 1);
			}
		}


		ImGui::Separator();
		static bool withAttachments = true;
		ImGui::Checkbox("With attachments", &withAttachments);
		if (ImGui::Button("Output colored curvature mesh"))
		{
			if (m_vertexMaxCurvature.size() > 0)
			{
				Eigen::MatrixXd V = (Eigen::Map<const Eigen::MatrixXd>(visualize_m_x.data(), dimension, visualize_m_x.size() / dimension));
				Eigen::MatrixXi F = m_ipc_collisionMesh.faces().transpose();

				Eigen::MatrixXd colors(3, m_vertexMaxCurvature.size());
				Eigen::VectorXd abs_vertexGaussianCurvature = m_vertexMaxCurvature;// .cwiseAbs();
				double maxCoeff = abs_vertexGaussianCurvature.maxCoeff();
				double minCoeff = abs_vertexGaussianCurvature.minCoeff();
				for (int v_i = 0; v_i < m_vertexMaxCurvature.size(); v_i++)
				{
					V3D rgb;
					colormap::colormapFromScalar(m_vertexMaxCurvature[v_i], minCoeff, maxCoeff, rgb[0], rgb[1], rgb[2], colormap::JET_COLORMAP);
					colors.col(v_i) = rgb;
				}

				//add attachments
				Eigen::MatrixXd sphereVertices;
				Eigen::MatrixXi sphereTriangles;
				SimpleTriMesh::icosphere(1).toMatrices(sphereVertices, sphereTriangles);
				if (withAttachments)
				{
					double radius = 0.001;
					for (int i = 0; i < m_selectedPairs.size(); i++)
					{

						for (int pair_i = 0; pair_i < 2; pair_i++)
						{
							if (pair_i == 0 && m_selectedPairs[i].first == -1)
							{//m_primitiveDrawer->drawSphere(0.001, visualize_m_x.segment<3>(m_selectedPairs[i].first * 3), m_pinningPairConstraintColor[i]);
								continue;
							}
							if (pair_i == 1 && m_selectedPairs[i].second == -1)
							{//m_primitiveDrawer->drawSphere(0.001, visualize_m_x.segment<3>(m_selectedPairs[i].second * 3), m_pinningPairConstraintColor[i]);
								continue;
							}

							//v
							Eigen::Vector3d translation = pair_i == 0 ?
								visualize_m_x.segment<3>(m_selectedPairs[i].first * 3)
								: visualize_m_x.segment<3>(m_selectedPairs[i].second * 3);
							Eigen::MatrixXd transformedSphereVertices = sphereVertices * radius;
							transformedSphereVertices.colwise() += translation;

							int oldVSize = V.cols();
							V.conservativeResize(3, oldVSize + transformedSphereVertices.cols());
							V.rightCols(transformedSphereVertices.cols()) = transformedSphereVertices;

							//f
							int oldFSize = F.cols();
							Eigen::MatrixXi transformedSphereFaces = sphereTriangles;
							transformedSphereFaces.colwise() += Eigen::Vector3i(oldVSize, oldVSize, oldVSize);
							F.conservativeResize(3, oldFSize + sphereTriangles.cols());
							F.rightCols(transformedSphereFaces.cols()) = transformedSphereFaces;


							//color
							int oldCSize = colors.cols();
							colors.conservativeResize(3, oldCSize + transformedSphereVertices.cols());
							colors.rightCols(transformedSphereVertices.cols()).colwise() = m_pinningPairConstraintColor[i];
						}


					}
				}
				//SimOpt::saveOBJ((m_outputPath + "/coloredMesh.obj").c_str(), V, F, colors);
			}
		}
	}
	ImGui::End();

}

void SOReflectionBendingModeShellStructure::render()
{
	{
		if (m_isRenderStrainColor)
		{ 
			double maximumStrain = 0.7;
			Eigen::MatrixXd facesStrainColors;
			computeMaxStrainFacesColorMap(maximumStrain, m_triMesh->vertices(), m_X, facesStrainColors);

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
			Eigen::Vector4d color(0.8, 0.8, 0.8, 1.0);
			//Eigen::Vector4d color(65.0 / 255.0, 113.0 / 255.0, 156.0 / 255.0, 1.0);
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


			if (m_drawReflectionPatches)
			{//draw periodic patches

				int startRows = -2;
				int startCols = -2;
				int numRows = 4;
				int numCols = 4;
				std::vector<Eigen::MatrixXd> vertices_vector;
				std::vector<Eigen::MatrixXi> faces_vector;
				std::vector<bool> flipNorm_vector;
				computeReflectionMeshes(startRows, numRows, startCols, numCols, vertices_vector, faces_vector, flipNorm_vector);

				for (int m_i = 0; m_i < vertices_vector.size(); m_i++)
				{
					m_meshDrawer->setMesh(vertices_vector[m_i], faces_vector[m_i], flipNorm_vector[m_i]);

					if (m_useBothWireframeAndmesh)
					{
						m_meshDrawer->drawWireframe(Eigen::Vector4d(0.0, 0.0, 0.0, 1.0));
						m_meshDrawer->draw(color);
					}
					else if (m_drawWireframe)
						m_meshDrawer->drawWireframe(Eigen::Vector4d(0.0, 0.0, 0.0, 1.0));
					else
						m_meshDrawer->draw(color);
				}
				
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
}
