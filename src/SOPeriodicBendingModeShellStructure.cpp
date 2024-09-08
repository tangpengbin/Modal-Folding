#include "SOPeriodicBendingModeShellStructure.h"
#include <SimLib/Core/SOUtils.h>
#include <SimLib/Viewer/colormap.h>

#include <fstream>

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

#include "tbb/tbb.h"

#include <Spectra/SymEigsShiftSolver.h>
#include <Spectra/MatOp/SparseSymShiftSolve.h>

using namespace soutil;
const int dimension = 3;
SOPeriodicBendingModeShellStructure::SOPeriodicBendingModeShellStructure()
{
	USE_MAKEPD = true;
	USE_GRAVITY = false;
	m_usePeriodicBoundaryCondition = true;
	m_drawPeriodicPathes = false;
	m_drawPeriodicTranslationPair = false;
	m_periodicPairIndex = 0;
	
	m_drawPeriodicBendingElements = false;
	m_periodicBendingElementIndex = 0;

	m_numExtraDofs = 4;
}

SOPeriodicBendingModeShellStructure::SOPeriodicBendingModeShellStructure(const SOPeriodicBendingModeShellStructure& other)
{
	(*this) = other;
	reconstructElementVector();
}

SOPeriodicBendingModeShellStructure::~SOPeriodicBendingModeShellStructure()
{
}

void SOPeriodicBendingModeShellStructure::init(const TriMesh* pTriMesh, SOMaterial* mat, 
	double thickness, 
	const std::vector<std::pair<int, int>>& fixVertexOnLinePairs, 
	const std::vector<std::pair<int, int>>& fixingVerticesDofs, 
	const std::vector<std::pair<int, Eigen::Vector3d>>& vertexIndex_force,
	const Eigen::Vector3d& periodicLeftToRight,
	const Eigen::Vector3d& periodicTopToBottom,
	string outputPath)
{
	m_periodicLeftToRight_rest = periodicLeftToRight;
	m_periodicTopToBottom_rest = periodicTopToBottom;

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

	m_periodicCSTShellElement.clear();
	m_periodicDSHingeFixedRestPoseElement.clear();



	m_nv = m_triMesh->numVertices();
	//copy positions from mesh
	m_X = m_triMesh->vertices();//m_X will keep the initial state
	computePeriodicity();
	processSystemDoFs();


	initElements(m_para_X);

	reconstructMassStuff(m_para_X);

	initIPCConfig();


	/*Eigen::VectorXd random(m_x.size()); random.setRandom(); random *= 1e-2;
	m_x += random;*/
	setVisualizationConfiguration(m_x);

	spdlog::set_level(spdlog::level::info);
}

void SOPeriodicBendingModeShellStructure::initIPCConfig()
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


void SOPeriodicBendingModeShellStructure::processSystemDoFs()
{
	/*convert all dofs to system dofs*/

	int numPeriodicPoints = m_allPeriodicVertices.size();
	int numPeriodicSystemVertices = m_nv - numPeriodicPoints + numPeriodicPoints / 2 - 1;
	m_x.setZero(numPeriodicSystemVertices * 3 + m_numExtraDofs);//plus translations
	int systemDofIndex = 0;
	m_system_to_allDofs.clear();
	m_all_to_systemDofs.clear();
	for (int i = 0; i < m_nv; i++)
	{
		if (m_allPeriodicVertices.find(i) == m_allPeriodicVertices.end())
		{//this is inner vertex
			m_x.segment<3>(systemDofIndex * 3) = m_X.segment<3>(3 * i);
			m_system_to_allDofs[systemDofIndex] = i;
			m_all_to_systemDofs[i] = systemDofIndex;
			systemDofIndex++;
		}
		else
		{
			if (m_allPeriodicVertices[i] == periodicBoundaryVertexType::topleft ||
				m_allPeriodicVertices[i] == periodicBoundaryVertexType::top ||
				m_allPeriodicVertices[i] == periodicBoundaryVertexType::left)
			{
				m_x.segment<3>(systemDofIndex * 3) = m_X.segment<3>(3 * i);
				m_system_to_allDofs[systemDofIndex] = i;
				m_all_to_systemDofs[i] = systemDofIndex;
				systemDofIndex++;
			}
		}
	}
	
	//the above process will fill all system_to_all, but vertices not belong to system will not be added to all, so we add some other pbc vertices
	for (int i = 1; i < 4; i++)
	{//set to topleft for other three corners
		int topleft_systemDofIndex = m_all_to_systemDofs.at(m_cornerVerticesIndices[0]);
		m_all_to_systemDofs[m_cornerVerticesIndices[i]] = topleft_systemDofIndex;
	}
	for (const auto& item: m_leftOrTopIndex_to_rightOrBottomWithTrans)
	{
		bool isCorners = false;
		for (int i = 0; i < 4; i++)
		{
			if (item.first == m_cornerVerticesIndices[i])
			{
				isCorners = true;
				break;
			}
		}
		if (isCorners)
			continue;

		//since others are all belong to --> top to bottom OR left to right, we directly set them
		int leftOrTop_systemDofIndex = m_all_to_systemDofs.at(item.first);
		m_all_to_systemDofs[item.second[0].first] = leftOrTop_systemDofIndex;
	}

	//periodic translations (leftToright, topTobottom)
	m_x.tail(m_numExtraDofs) = Eigen::Vector4d(m_periodicLeftToRight_rest[0], m_periodicLeftToRight_rest[1],
		m_periodicTopToBottom_rest[0], m_periodicTopToBottom_rest[1]);


	Eigen::VectorXd allDofs = convertSystemDofsToAllDofs(m_x);
	spdlog::info("Converted system dofs vs. m_X is {}",(allDofs - m_X).norm());


	m_para_X = m_x;
	m_V.setZero(m_x.size());			 //!!!!!!now this velocity is zero, TODO: parameters
	
	m_v = Eigen::VectorXd(m_x.size()); m_v.setZero();
	m_appliedForces = Eigen::VectorXd(m_x.size()); m_appliedForces.setZero();
}

Eigen::VectorXd SOPeriodicBendingModeShellStructure::convertAllDofsToSystemDofs(const Eigen::VectorXd& allDofs)
{
	Eigen::VectorXd systemDofs(m_x.size());
	int systemDofIndex = 0;
	for (int i = 0; i < m_nv; i++)
	{
		if (m_allPeriodicVertices.find(i) == m_allPeriodicVertices.end())
		{//this is inner vertex
			systemDofs.segment<3>(systemDofIndex * 3) = allDofs.segment<3>(3 * i);
			systemDofIndex++;
		}
		else
		{
			if (m_allPeriodicVertices[i] == periodicBoundaryVertexType::topleft ||
				m_allPeriodicVertices[i] == periodicBoundaryVertexType::top ||
				m_allPeriodicVertices[i] == periodicBoundaryVertexType::left)
			{
				systemDofs.segment<3>(systemDofIndex * 3) = allDofs.segment<3>(3 * i);
				systemDofIndex++;
			}
		}
	}
	Eigen::Vector2d periodicLeftToRight = allDofs.segment<2>(m_periodicReferenceLeftRight.second * 3) - allDofs.segment<2>(m_periodicReferenceLeftRight.first*3);
	Eigen::Vector2d periodicTopToBottom = allDofs.segment<2>(m_periodicReferenceTopBottom.second * 3) - allDofs.segment<2>(m_periodicReferenceTopBottom.first * 3);


	//periodic translations (leftToright, topTobottom)
	systemDofs.tail(m_numExtraDofs) = Eigen::Vector4d(periodicLeftToRight[0], periodicLeftToRight[1],
		periodicTopToBottom[0], periodicTopToBottom[1]);
	return systemDofs;
}

Eigen::VectorXd SOPeriodicBendingModeShellStructure::convertSystemDofsToAllDofs(const Eigen::VectorXd& systemDofs)
{
	//systemDofs is system dof, we convert it to the all dof
	Eigen::VectorXd allDofs(m_nv * 3);
	Eigen::Vector4d translations = systemDofs.tail(m_numExtraDofs);
	Eigen::Vector3d leftToRightTranslation; leftToRightTranslation.setZero(); leftToRightTranslation.head(2) = translations.head(2);
	Eigen::Vector3d topToBottomTranslation; topToBottomTranslation.setZero(); topToBottomTranslation.head(2) = translations.tail(2);
	int numConvertedVertices = 0;
	for (const auto& systemDof : m_system_to_allDofs)
	{
		int systemDofIndex = systemDof.first;
		int allDofIndex = systemDof.second;
		//inner, topleft, top and left
		Eigen::Vector3d thisVertexPos = systemDofs.segment<3>(systemDofIndex * 3);
		allDofs.segment<3>(allDofIndex * 3) = thisVertexPos;
		numConvertedVertices++;

		const auto& periodicVertexFinder = m_allPeriodicVertices.find(allDofIndex);
		if (periodicVertexFinder != m_allPeriodicVertices.end())
		{//this is periodic vertex (topleft, left, top)
			if (periodicVertexFinder->second == periodicBoundaryVertexType::topleft)
			{//find topright, bottomleft, bottomright
				
				int toprightvAllDofIndex = m_cornerVerticesIndices[1];
				allDofs.segment<3>(toprightvAllDofIndex * 3) = thisVertexPos + leftToRightTranslation;

				int bottomleftvAllDofIndex = m_cornerVerticesIndices[2];
				allDofs.segment<3>(bottomleftvAllDofIndex * 3) = thisVertexPos + topToBottomTranslation;

				int bottomrightvAllDofIndex = m_cornerVerticesIndices[3];
				allDofs.segment<3>(bottomrightvAllDofIndex * 3) = thisVertexPos + topToBottomTranslation + leftToRightTranslation;

				numConvertedVertices += 3;
			}
			else if (periodicVertexFinder->second == periodicBoundaryVertexType::top)
			{//find bottom
				int bottomvAllDofIndex = m_leftOrTopIndex_to_rightOrBottomWithTrans.at(allDofIndex)[0].first;
				allDofs.segment<3>(bottomvAllDofIndex * 3) = thisVertexPos + topToBottomTranslation;
				numConvertedVertices++;
			}
			else if (periodicVertexFinder->second == periodicBoundaryVertexType::left)
			{//find right
				int rightvAllDofIndex = m_leftOrTopIndex_to_rightOrBottomWithTrans.at(allDofIndex)[0].first;
				allDofs.segment<3>(rightvAllDofIndex * 3) = thisVertexPos + leftToRightTranslation;
				numConvertedVertices++;
			}
		}
	}
	assert(numConvertedVertices == m_nv);

	return allDofs;
}

SpMat SOPeriodicBendingModeShellStructure::convertSystemDofsToAllDofsJacobian(const Eigen::VectorXd& systemDofs)
{
	SpMat jacobian(3 * m_nv, systemDofs.size());
	TripVec jacobian_triplets;

	//systemDofs is system dof, we convert it to the all dof
	//Eigen::VectorXd allDofs(m_nv * 3);
	Eigen::Vector4d translations = systemDofs.tail(m_numExtraDofs);
	int leftToRightTranslation_dofStartIndex = systemDofs.size() - 4 - 1;
	int topToBottomTranslation_dofStartIndex = systemDofs.size() - 4 - 1 + 2;
	Eigen::Vector3d leftToRightTranslation; leftToRightTranslation.setZero(); leftToRightTranslation.head(2) = translations.head(2);
	Eigen::Vector3d topToBottomTranslation; topToBottomTranslation.setZero(); topToBottomTranslation.head(2) = translations.tail(2);
	for (const auto& systemDof : m_system_to_allDofs)
	{
		int systemDofIndex = systemDof.first;
		int allDofIndex = systemDof.second;
		//inner, topleft, top and left
		//Eigen::Vector3d thisVertexPos = systemDofs.segment<3>(systemDofIndex * 3);
		//allDofs.segment<3>(allDofIndex * 3) = thisVertexPos;
		for (int allDof_i = 0; allDof_i < 3; allDof_i++)
		{
			jacobian_triplets.emplace_back(allDofIndex * 3 + allDof_i, systemDofIndex * 3 + allDof_i, 1.0);
		}

		const auto& periodicVertexFinder = m_allPeriodicVertices.find(allDofIndex);
		if (periodicVertexFinder != m_allPeriodicVertices.end())
		{//this is periodic vertex (topleft, left, top)
			if (periodicVertexFinder->second == periodicBoundaryVertexType::topleft)
			{//find topright, bottomleft, bottomright

				int toprightvAllDofIndex = m_cornerVerticesIndices[1];
				//allDofs.segment<3>(toprightvAllDofIndex * 3) = thisVertexPos + leftToRightTranslation;

				int bottomleftvAllDofIndex = m_cornerVerticesIndices[2];
				//allDofs.segment<3>(bottomleftvAllDofIndex * 3) = thisVertexPos + topToBottomTranslation;

				int bottomrightvAllDofIndex = m_cornerVerticesIndices[3];
				//allDofs.segment<3>(bottomrightvAllDofIndex * 3) = thisVertexPos + topToBottomTranslation + leftToRightTranslation;

				for (int allDof_i = 0; allDof_i < 3; allDof_i++)
				{
					jacobian_triplets.emplace_back(toprightvAllDofIndex * 3 + allDof_i, systemDofIndex * 3 + allDof_i, 1.0);
					jacobian_triplets.emplace_back(bottomleftvAllDofIndex * 3 + allDof_i, systemDofIndex * 3 + allDof_i, 1.0);

					jacobian_triplets.emplace_back(bottomrightvAllDofIndex * 3 + allDof_i, systemDofIndex * 3 + allDof_i, 1.0);
				}
				for (int translationDof_i = 0; translationDof_i < 2; translationDof_i++)
				{
					jacobian_triplets.emplace_back(toprightvAllDofIndex * 3 + translationDof_i, leftToRightTranslation_dofStartIndex + translationDof_i, 1.0);
					jacobian_triplets.emplace_back(bottomleftvAllDofIndex * 3 + translationDof_i, topToBottomTranslation_dofStartIndex + translationDof_i, 1.0);

					jacobian_triplets.emplace_back(bottomrightvAllDofIndex * 3 + translationDof_i, leftToRightTranslation_dofStartIndex + translationDof_i, 1.0);
					jacobian_triplets.emplace_back(bottomrightvAllDofIndex * 3 + translationDof_i, topToBottomTranslation_dofStartIndex + translationDof_i, 1.0);
				}
			}
			else if (periodicVertexFinder->second == periodicBoundaryVertexType::top)
			{//find bottom
				int bottomvAllDofIndex = m_leftOrTopIndex_to_rightOrBottomWithTrans.at(allDofIndex)[0].first;
				//allDofs.segment<3>(bottomvAllDofIndex * 3) = thisVertexPos + topToBottomTranslation;

				for (int allDof_i = 0; allDof_i < 3; allDof_i++)
				{
					jacobian_triplets.emplace_back(bottomvAllDofIndex * 3 + allDof_i, systemDofIndex * 3 + allDof_i, 1.0);
				}
				for (int translationDof_i = 0; translationDof_i < 2; translationDof_i++)
				{
					jacobian_triplets.emplace_back(bottomvAllDofIndex * 3 + translationDof_i, topToBottomTranslation_dofStartIndex + translationDof_i, 1.0);
				}
			}
			else if (periodicVertexFinder->second == periodicBoundaryVertexType::left)
			{//find right
				int rightvAllDofIndex = m_leftOrTopIndex_to_rightOrBottomWithTrans.at(allDofIndex)[0].first;
				//allDofs.segment<3>(rightvAllDofIndex * 3) = thisVertexPos + leftToRightTranslation;

				for (int allDof_i = 0; allDof_i < 3; allDof_i++)
				{
					jacobian_triplets.emplace_back(rightvAllDofIndex * 3 + allDof_i, systemDofIndex * 3 + allDof_i, 1.0);
				}
				for (int translationDof_i = 0; translationDof_i < 2; translationDof_i++)
				{
					jacobian_triplets.emplace_back(rightvAllDofIndex * 3 + translationDof_i, leftToRightTranslation_dofStartIndex + translationDof_i, 1.0);
				}
			}
		}
	}
	jacobian.setFromTriplets(jacobian_triplets.begin(), jacobian_triplets.end());

	return jacobian;
}

void SOPeriodicBendingModeShellStructure::computePeriodicity()
{
	m_leftOrTopIndex_to_rightOrBottomWithTrans.clear();
	m_periodicReferenceTopBottom = std::make_pair(-1, -1);
	m_periodicReferenceLeftRight = std::make_pair(-1, -1);
	
	double numericalError = 1e-7;
	for (int i = 0; i < m_nv; i++)
	{
		for (int j = i + 1; j < m_nv; j++)
		{
			//top to bottom OR left to right
			Eigen::Vector3d t = m_X.segment<3>(j * 3) - m_X.segment<3>(i * 3);

			if ((t - m_periodicLeftToRight_rest).norm() < numericalError)
			{
				m_leftOrTopIndex_to_rightOrBottomWithTrans[i].emplace_back(j, &m_periodicLeftToRight_rest);
			}
			else if ((t - m_periodicTopToBottom_rest).norm() < numericalError)
			{
				m_leftOrTopIndex_to_rightOrBottomWithTrans[i].emplace_back(j, &m_periodicTopToBottom_rest);
			}
			else if ((t + m_periodicLeftToRight_rest).norm() < numericalError)
			{
				m_leftOrTopIndex_to_rightOrBottomWithTrans[j].emplace_back(i, &m_periodicLeftToRight_rest);
			}
			else if ((t + m_periodicTopToBottom_rest).norm() < numericalError)
			{
				m_leftOrTopIndex_to_rightOrBottomWithTrans[j].emplace_back(i, &m_periodicTopToBottom_rest);
			}
		}
	}


	if (m_leftOrTopIndex_to_rightOrBottomWithTrans.size() > 0)
	{
		m_usePeriodicBoundaryCondition = true;
		spdlog::info("Build periodic boundary conditions");
		m_periodicReferenceTopBottom = std::make_pair(-1, -1);
		m_periodicReferenceLeftRight = std::make_pair(-1, -1);
		for (const auto& item : m_leftOrTopIndex_to_rightOrBottomWithTrans)
		{
			if (m_periodicReferenceTopBottom.first != -1 && m_periodicReferenceTopBottom.second != -1 &&
				m_periodicReferenceLeftRight.first != -1 && m_periodicReferenceLeftRight.second != -1)
				break;
			int topOrLeftIndex = item.first;
			for (int j = 0; j < item.second.size(); j++)
			{
				int bottomOrRightIndex = item.second[j].first;

				if (item.second[j].second == &m_periodicTopToBottom_rest && m_periodicReferenceTopBottom.first == -1 && m_periodicReferenceTopBottom.second == -1)
				{
					m_periodicReferenceTopBottom.first = topOrLeftIndex;
					m_periodicReferenceTopBottom.second = bottomOrRightIndex;
				}
				else if (item.second[j].second == &m_periodicLeftToRight_rest && m_periodicReferenceLeftRight.first == -1 && m_periodicReferenceLeftRight.second == -1)
				{
					m_periodicReferenceLeftRight.first = topOrLeftIndex;
					m_periodicReferenceLeftRight.second = bottomOrRightIndex;
				}
			}
		}


		for (const auto& periodic_pairs : m_leftOrTopIndex_to_rightOrBottomWithTrans)
		{
			int topOrLeftVertexIndex = periodic_pairs.first;
			if (m_allPeriodicVertices.find(topOrLeftVertexIndex) == m_allPeriodicVertices.end())
			{
				if (periodic_pairs.second.size() == 2)
				{//this should be topleft
					//then, we directly find and set bottomleft, topright, bottomright
					int toprightVertexIndex = -1, bottomleftVertexIndex = -1;
					if (periodic_pairs.second[0].second == &m_periodicTopToBottom_rest
						&& periodic_pairs.second[1].second == &m_periodicLeftToRight_rest)
					{
						bottomleftVertexIndex = periodic_pairs.second[0].first;
						toprightVertexIndex = periodic_pairs.second[1].first;
					}
					else
					{
						bottomleftVertexIndex = periodic_pairs.second[1].first;
						toprightVertexIndex = periodic_pairs.second[0].first;
					}

					const auto& toprightPeriodicPair = m_leftOrTopIndex_to_rightOrBottomWithTrans.find(toprightVertexIndex);
					const auto& bottomleftPeriodicPair = m_leftOrTopIndex_to_rightOrBottomWithTrans.find(bottomleftVertexIndex);
					assert(toprightPeriodicPair != m_leftOrTopIndex_to_rightOrBottomWithTrans.end());
					assert(bottomleftPeriodicPair != m_leftOrTopIndex_to_rightOrBottomWithTrans.end());
					assert(toprightPeriodicPair->second.size() == 1);
					assert(bottomleftPeriodicPair->second.size() == 1);

					assert(toprightPeriodicPair->second[0].first == bottomleftPeriodicPair->second[0].first);//these two should point to the same bottomright
					int bottomrightVertexIndex = toprightPeriodicPair->second[0].first;

					//assign four corners
					m_allPeriodicVertices[topOrLeftVertexIndex] = periodicBoundaryVertexType::topleft;
					m_allPeriodicVertices[toprightVertexIndex] = periodicBoundaryVertexType::topright;
					m_allPeriodicVertices[bottomleftVertexIndex] = periodicBoundaryVertexType::bottomleft;
					m_allPeriodicVertices[bottomrightVertexIndex] = periodicBoundaryVertexType::bottomright;

					//build quick access for four corners
					m_cornerVerticesIndices = Eigen::Vector4i(topOrLeftVertexIndex, toprightVertexIndex, bottomleftVertexIndex, bottomrightVertexIndex);
				}
				else
				{//top or left
					int bottomOrRightVertexIndex = periodic_pairs.second[0].first;
					if (periodic_pairs.second[0].second == &m_periodicLeftToRight_rest)
					{//left to right
						m_allPeriodicVertices[topOrLeftVertexIndex] = periodicBoundaryVertexType::left;
						m_allPeriodicVertices[bottomOrRightVertexIndex] = periodicBoundaryVertexType::right;
					}
					else
					{//top to bottom
						m_allPeriodicVertices[topOrLeftVertexIndex] = periodicBoundaryVertexType::top;
						m_allPeriodicVertices[bottomOrRightVertexIndex] = periodicBoundaryVertexType::bottom;
					}
				}
			}
		}
	}
	else
	{
		m_usePeriodicBoundaryCondition = false;
		spdlog::info("No periodic boundary conditions");
	}
}


void SOPeriodicBendingModeShellStructure::reconstructMassStuff(const Eigen::VectorXd& vX)
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

void SOPeriodicBendingModeShellStructure::initElements(const Eigen::VectorXd& vX)
{
	auto computeTranslationMultiplier = [this](int vertexIndexInAll)
	{
		const auto& findPeriodicVertex = m_allPeriodicVertices.find(vertexIndexInAll);
		if (findPeriodicVertex != m_allPeriodicVertices.end())
		{//boundary vertices
			if (findPeriodicVertex->second == periodicBoundaryVertexType::bottom)
			{
				return Eigen::Vector2d(0.0, 1.0);
			}
			else if (findPeriodicVertex->second == periodicBoundaryVertexType::right)
			{
				return Eigen::Vector2d(1.0, 0.0);
			}
			else if (findPeriodicVertex->second == periodicBoundaryVertexType::bottomleft)
			{
				return Eigen::Vector2d(0.0, 1.0);
			}
			else if (findPeriodicVertex->second == periodicBoundaryVertexType::topright)
			{
				return Eigen::Vector2d(1.0, 0.0);
			}
			else if (findPeriodicVertex->second == periodicBoundaryVertexType::bottomright)
			{
				return Eigen::Vector2d(1.0, 1.0);
			}
			else
				return Eigen::Vector2d(0, 0);
		}
		else
			return Eigen::Vector2d(0.0, 0.0);
	};

	//assert(vX == m_X);
	m_periodicCSTShellElement.clear();
	{
		const std::vector<int>& faces = m_triMesh->faces();
		int nf = faces.size() / 3;
		for (int i = 0; i < nf; i++)
		{
			std::vector<int> inds;
			std::vector<Eigen::Vector2d> translationMultipliers;
			for (int j = 0; j < 3; j++)
			{
				int vIndex_j = faces[3 * i + j];
				Eigen::Vector2d mutiplier_j = computeTranslationMultiplier(vIndex_j);
				translationMultipliers.push_back(mutiplier_j);//multiplier for leftToright, topTobottom

				int systemVIndex_j = m_all_to_systemDofs.at(vIndex_j);
				inds.push_back(systemVIndex_j);
			}

			m_periodicCSTShellElement.push_back(periodicCSTShellElement());
			periodicCSTShellElement& elem = m_periodicCSTShellElement.back();

			std::vector<int> dofIndices;
			for (int dof_i = 0; dof_i < inds.size(); dof_i++)
			{
				for (int j = 0; j < dimension; j++)
					dofIndices.push_back(inds[dof_i] * dimension + j);
			}
			//push translation dofs to the end
			int translationDofStartIndex = m_x.size() - m_numExtraDofs;
			for (int i = 0; i < m_numExtraDofs; i++)
			{
				dofIndices.push_back(translationDofStartIndex + i);
			}

			elem.init(dofIndices, dofIndices, translationMultipliers, vX, m_material);
			elem.setThickness(m_cstThickness);
			elem.setElementProperty(InternalElement);
			double restEnergy = elem.computeEnergy(vX, vX);
			if (isnan(restEnergy))
			{
				spdlog::warn("The rest energy of cst element {} is nan", i);
			}
		}
	}

	{
		m_triMesh->buildHingeStructure();

		//elements for periodic hinge elements for bending at periodic boundaries
		std::map<int, std::vector<int>> boundaryVertex_map_boundaryHinge;
		auto hinges = m_triMesh->hinges();
		for (int i = 0; i < hinges.size(); i++)
		{
			auto hinge = hinges[i];
			if ((hinge.tris[0] == -1) || (hinge.tris[1] == -1))
			{//for boundary edges
				boundaryVertex_map_boundaryHinge[hinge.edge[0]].push_back(i);
				boundaryVertex_map_boundaryHinge[hinge.edge[1]].push_back(i);
			}
		}

		for (int i = 0; i < hinges.size(); i++)
		{
			auto hinge = hinges[i];
			if ((hinge.tris[0] == -1) || (hinge.tris[1] == -1))
			{//for boundary edges
				//we only care edges at top and left sides
				int edgeV0 = hinge.edge[0];
				int edgeV1 = hinge.edge[1];
				if (m_leftOrTopIndex_to_rightOrBottomWithTrans.find(edgeV0) == m_leftOrTopIndex_to_rightOrBottomWithTrans.end() ||
					m_leftOrTopIndex_to_rightOrBottomWithTrans.find(edgeV1) == m_leftOrTopIndex_to_rightOrBottomWithTrans.end())
					continue;//not at left or top edges


				/*we need to find flap1*/
				//we first find opposite edge vertices
				int oppositeEdge0 = -1;
				int oppositeEdge1 = -1;
				Eigen::Vector3d* oppositeTranslation = nullptr;
				for (const auto& ev0_opposite : m_leftOrTopIndex_to_rightOrBottomWithTrans[edgeV0])
				{
					for (const auto& ev1_opposite : m_leftOrTopIndex_to_rightOrBottomWithTrans[edgeV1])
					{
						if (ev0_opposite.second == ev1_opposite.second)
						{
							oppositeEdge0 = ev0_opposite.first;
							oppositeEdge1 = ev1_opposite.first;
							oppositeTranslation = ev0_opposite.second;
							break;
						}
					}
					if (oppositeEdge0 != -1 && oppositeEdge1 != -1)
						break;
				}
				int oppositeHingeIndex = -1;
				for (const auto& oppositeEdgeV0_Hinge : boundaryVertex_map_boundaryHinge.at(oppositeEdge0))
				{
					for (const auto& oppositeEdgeV1_Hinge : boundaryVertex_map_boundaryHinge.at(oppositeEdge1))
					{
						if (oppositeEdgeV0_Hinge == oppositeEdgeV1_Hinge)
						{
							oppositeHingeIndex = oppositeEdgeV0_Hinge;
							break;
						}
					}
					if (oppositeHingeIndex != -1)
						break;
				}


				int flap0 = (hinge.flaps[0] == -1 ? hinge.flaps[1] : hinge.flaps[0]); //belongs to the left and top edge
				int flap1 = (hinges[oppositeHingeIndex].flaps[0] == -1 ? hinges[oppositeHingeIndex].flaps[1] : hinges[oppositeHingeIndex].flaps[0]);
				//we only build periodic hinge elements for top and left
				std::vector<int> vinds_inAll(4);
				vinds_inAll[0] = hinge.edge[0];
				vinds_inAll[1] = hinge.edge[1];
				vinds_inAll[2] = flap0;
				vinds_inAll[3] = flap1;//opposite side flap vindex
				//inds[4] = oppositeEdge0;
				//inds[5] = oppositeEdge1;
				//printf("edge %d %d face %d %d\n", inds[0], inds[1], inds[2], inds[3]);
				std::vector<int> inds(4);
				std::vector<Eigen::Vector2d> translationMultipliers;
				for (int j = 0; j < 4; j++)
				{
					int vIndex_j = vinds_inAll[j];
					if(j < 3)
						translationMultipliers.push_back(computeTranslationMultiplier(vIndex_j));//multiplier for leftToright, topTobottom
					else
					{//for the opposite vertex, we need to translate it to the corresponding side according to the translations of edge
						if (oppositeTranslation == &m_periodicLeftToRight_rest)
							translationMultipliers.push_back(Eigen::Vector2d(-1.0, 0));
						else
							translationMultipliers.push_back(Eigen::Vector2d(0.0, -1.0));
					}
					int systemVIndex_j = m_all_to_systemDofs.at(vIndex_j);
					inds[j] = systemVIndex_j;
				}


				m_periodicDSHingeFixedRestPoseElement.push_back(periodicDSHingeFixedRestPoseElement());
				periodicDSHingeFixedRestPoseElement& elem = m_periodicDSHingeFixedRestPoseElement.back();

				std::vector<int> dofIndices;
				for (int dof_i = 0; dof_i < inds.size(); dof_i++)
				{
					for (int j = 0; j < dimension; j++)
						dofIndices.push_back(inds[dof_i] * dimension + j);
				}
				//push translation dofs to the end
				int translationDofStartIndex = m_x.size() - m_numExtraDofs;
				for (int i = 0; i < m_numExtraDofs; i++)
				{
					dofIndices.push_back(translationDofStartIndex + i);
				}

				elem.init(dofIndices, dofIndices, translationMultipliers, m_para_X, m_material, m_cstThickness);
				//elem.setThickness(m_cstThickness);
				elem.setElementProperty(InternalElement);
				//double initialEnergy = elem.computeEnergy(m_para_X, m_para_X);
				//std::cout << "initial periodic bending energy " << initialEnergy << std::endl;

				//double restEnergy = elem.computeEnergy(vX, vX);
				//if (isnan(restEnergy))
				//{
				//	spdlog::warn("The rest energy of hinge element {} is {}", i, restEnergy);
				//}
			}
			else
			{//for inner hinges
				auto hinge = hinges[i];

				std::vector<int> vinds_inAll(4);
				vinds_inAll[0] = hinge.edge[0];
				vinds_inAll[1] = hinge.edge[1];
				vinds_inAll[2] = hinge.flaps[0];
				vinds_inAll[3] = hinge.flaps[1];
				//printf("edge %d %d face %d %d\n", vinds_inAll[0], vinds_inAll[1], vinds_inAll[2], vinds_inAll[3]);

				std::vector<int> inds;
				std::vector<Eigen::Vector2d> translationMultipliers;
				for (int j = 0; j < 4; j++)
				{
					int vIndex_j = vinds_inAll[j];
					translationMultipliers.push_back(computeTranslationMultiplier(vIndex_j));//multiplier for leftToright, topTobottom

					int systemVIndex_j = m_all_to_systemDofs.at(vIndex_j);
					inds.push_back(systemVIndex_j);
				}


				m_periodicDSHingeFixedRestPoseElement.push_back(periodicDSHingeFixedRestPoseElement());
				periodicDSHingeFixedRestPoseElement& elem = m_periodicDSHingeFixedRestPoseElement.back();

				std::vector<int> dofIndices;
				for (int dof_i = 0; dof_i < inds.size(); dof_i++)
				{
					for (int j = 0; j < dimension; j++)
						dofIndices.push_back(inds[dof_i] * dimension + j);
				}
				//push translation dofs to the end
				int translationDofStartIndex = m_x.size() - m_numExtraDofs;
				for (int i = 0; i < m_numExtraDofs; i++)
				{
					dofIndices.push_back(translationDofStartIndex + i);
				}

				elem.init(dofIndices, dofIndices, translationMultipliers, vX, m_material, m_cstThickness);
				//elem.setThickness(m_cstThickness);
				elem.setElementProperty(InternalElement);

				//spdlog::info("Theta of hinge is {}", elem.getCurrentAngle(vX));

				//double restEnergy = elem.computeEnergy(vX, vX);
				//if (isnan(restEnergy))
				//{
				//	spdlog::warn("The rest energy of hinge element {} is {}", i, restEnergy);
				//}
			}
		}

	}

	/*{
		m_targetPlanarParametersElement.clear();
		std::vector<std::vector<int>> targetDofs = { 
			{(int)vX.size() - m_numExtraDofs,(int)vX.size() - m_numExtraDofs + 1},
			{(int)vX.size() - m_numExtraDofs + 2, (int)vX.size() - m_numExtraDofs + 3} };
		std::vector<Eigen::Vector2d> targetTrans = { Eigen::Vector2d(0.195, 0), Eigen::Vector2d(0, -0.195) };
		for (int i = 0; i < targetDofs.size(); i++)
		{
			m_targetPlanarParametersElement.push_back(targetPlanarParametersElement());
			targetPlanarParametersElement& elem = m_targetPlanarParametersElement.back();

			elem.init(targetDofs[i], targetDofs[i], m_fix_stiffness, targetTrans[i]);
			elem.setElementProperty(ExternalElement);

		}

	}*/

	{

		m_fixDoFElements.clear();
		m_fixingDofs = { 0,1,2 };
		for (int i = 0; i < m_fixingDofs.size(); i++)
		{
			m_fixDoFElements.push_back(FixDoFElement());
			FixDoFElement& elem = m_fixDoFElements.back();

			elem.init(m_fixingDofs[i], m_fixingDofs[i], m_fix_stiffness, true);
			elem.setElementProperty(ExternalElement);

		}
	}
	
	reconstructElementVector();
}


void SOPeriodicBendingModeShellStructure::reconstructElementVector()
{
	m_elements.clear();

	for (int i = 0; i < (int)m_periodicCSTShellElement.size(); i++)
		m_elements.push_back(&m_periodicCSTShellElement[i]);
	for (int i = 0; i < (int)m_periodicDSHingeFixedRestPoseElement.size(); i++)
		m_elements.push_back(&m_periodicDSHingeFixedRestPoseElement[i]);
	
	
	//for (int i = 0; i < (int)m_targetPlanarParametersElement.size(); i++)
	//	m_elements.push_back(&m_targetPlanarParametersElement[i]);

	//for (int i = 0; i < (int)m_sinMovementElements.size(); i++)
	//	m_elements.push_back(&m_sinMovementElements[i]);
	for (int i = 0; i < (int)m_fixDoFElements.size(); i++)
		m_elements.push_back(&m_fixDoFElements[i]);

}

void SOPeriodicBendingModeShellStructure::computeMasses(const Eigen::VectorXd& vX)
{
	//const P3DIndexArray* pIB = m_triMesh->GetIndexArray();
	int nTriangles = m_triMesh->numTriangles();
	m_mass.setZero(m_x.size());

	for (int i = 0; i < m_periodicCSTShellElement.size(); i++)
	{
		double area = m_periodicCSTShellElement[i].computeArea(m_x);
		//A += area;
		//printf("%d %d %d %.15f\n", m_triMesh->faces()[3 * i + 0], m_triMesh->faces()[3 * i + 1], m_triMesh->faces()[3 * i + 2], area);
		for (int j = 0; j < 9; j++)
		{
			m_mass[m_periodicCSTShellElement[i].getDofIndices()[j]] += 1.0 / 3.0 * area * m_cstThickness * m_material->rho;
		}
	}

	spdlog::info("Overall mass is {}", m_mass.sum() / 3.0);

}


bool SOPeriodicBendingModeShellStructure::staticSolve_newton(Eigen::VectorXd& vx, const Eigen::VectorXd& vX)
{

	auto stepBeginningFunction = [this, vX](const Eigen::VectorXd& vx, int currentIteration)
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
		{
			ipc::Constraints constraint_set;

			Eigen::VectorXd allDofs = convertSystemDofsToAllDofs(vx);
			Eigen::MatrixXd V = (Eigen::Map<const Eigen::MatrixXd>(allDofs.data(), dimension, allDofs.size() / dimension)).transpose();
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
				//for (size_t i = 0; i < m_elements.size(); i++)
				{
					Eigen::VectorXd element_grad;
					m_elements[i]->computeGradient(vx, vX, element_grad);

					const auto element_dofs = m_elements[i]->getDofIndices();
					for (int dof_i = 0; dof_i < element_dofs.size(); dof_i++)
					{
						local_grad[element_dofs[dof_i]] += element_grad[dof_i];
					}
					//spdlog::debug("grad norm {}", local_grad.norm());
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
		{
			ipc::Constraints constraint_set;
			Eigen::VectorXd allDofs = convertSystemDofsToAllDofs(vx);
			SpMat allDofs_jacobian = convertSystemDofsToAllDofsJacobian(vx);
			Eigen::MatrixXd V = (Eigen::Map<const Eigen::MatrixXd>(allDofs.data(), dimension, allDofs.size() / dimension)).transpose();
			constraint_set.build(m_ipc_collisionMesh, V, m_ipc_dhat);
			Eigen::VectorXd barrierPotentialGradient = ipc::compute_barrier_potential_gradient(m_ipc_collisionMesh, V, constraint_set, m_ipc_dhat);

			EselfContactGrad += m_ipc_adaptiveStiffness * allDofs_jacobian.transpose() * barrierPotentialGradient;
		}


		grad += EselfContactGrad;

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

	auto computeLargestStepSize = [this, vX](const Eigen::VectorXd& vx, const Eigen::VectorXd& dx)
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


void SOPeriodicBendingModeShellStructure::convertToStrainSpaceMode(const Eigen::VectorXd& positionSpaceMode, Eigen::VectorXd& strainSpaceMode)
{
	Eigen::VectorXd pm = positionSpaceMode;// .normalized();
	strainSpaceMode = Eigen::VectorXd(m_periodicDSHingeFixedRestPoseElement.size()); strainSpaceMode.setZero();
	
	SpMat all_dthetadx(m_periodicDSHingeFixedRestPoseElement.size(), m_x.size());
	TripVec all_dthetadx_triplets;
	for (int i = 0; i < m_periodicDSHingeFixedRestPoseElement.size(); i++)
	{
		Eigen::VectorXd dthetadx(16);
		dthetadx.setZero();
		m_periodicDSHingeFixedRestPoseElement[i].computedThetadx(m_para_X, m_para_X, dthetadx);
		//std::cout << "hinge " << i << " dthetadx norm " << dthetadx.norm() << std::endl;
		const auto& inds = m_periodicDSHingeFixedRestPoseElement[i].getDofIndices();
		
		for (int dof_i = 0; dof_i < inds.size(); dof_i++)
		{
			all_dthetadx_triplets.emplace_back(i, inds[dof_i], dthetadx[dof_i]);
		}
	}
	all_dthetadx.setFromTriplets(all_dthetadx_triplets.begin(), all_dthetadx_triplets.end());
	strainSpaceMode = all_dthetadx * pm;
}

void SOPeriodicBendingModeShellStructure::animateStrainSpaceMode()
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
	
	
	for (int i = 0; i < m_periodicDSHingeFixedRestPoseElement.size(); i++)
	{//std::cout << "hinge " << i << m_hingeElements[i].getRestAngle();
		double theta = m_periodicDSHingeFixedRestPoseElement[i].getCurrentAngle(m_x);
		theta += ssm[i];
		m_periodicDSHingeFixedRestPoseElement[i].setRestAngle(theta);
	}

	projectToManifold();

	time_passed += m_dt;
}

void SOPeriodicBendingModeShellStructure::projectToManifold()
{
	bool use_gravity = USE_GRAVITY;
	USE_GRAVITY = false;
	solveStaticState();
	USE_GRAVITY = use_gravity;

	//convert from system dofs to all dofs
	//setVisualizationConfiguration(m_x);

	//move to center
	V3D center(0, 0, 0);
	for (int i = 0; i < visualize_m_x.size() / 3; i++)
	{
		center += visualize_m_x.segment<3>(i * 3);
	}
	center /= (double)(visualize_m_x.size() / 3);

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


void SOPeriodicBendingModeShellStructure::computeStrainSpaceModes()
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

	int numSteps = (m_computeModeLength - initialLength) / m_computeModeStepLength;
	for (int i = 1; i <= numSteps; i++)
	{
		double length = std::min(initialLength + i * m_computeModeStepLength, m_computeModeLength);
		spdlog::info("Compute mode {} with strain space length {}", m_modesAnimationID, length);
		if (m_useStrainSpaceModes)
		{
			if (m_strainSpaceModes.size() >= m_modesAnimationID && m_strainSpaceModes[m_modesAnimationID].size() > 0)
			{
				Eigen::VectorXd ssm = m_strainSpaceModes[m_modesAnimationID] * length;

				for (int i = 0; i < m_periodicDSHingeFixedRestPoseElement.size(); i++)
				{
					double theta = m_periodicDSHingeFixedRestPoseElement[i].getCurrentAngle(initialGuess); // m_periodicDSHingeFixedRestPoseElement[i].getRestAngle();
					theta += ssm[i];
					m_periodicDSHingeFixedRestPoseElement[i].setRestAngle(theta);
				}

				m_x = initialGuess;
				projectToManifold();
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

	for (int i = 0; i < m_periodicDSHingeFixedRestPoseElement.size(); i++)
	{
		m_periodicDSHingeFixedRestPoseElement[i].setRestAngle(m_periodicDSHingeFixedRestPoseElement[i].getInitialRestAngle());
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

void SOPeriodicBendingModeShellStructure::ouputStrainSpaceModesAtIndex(int index)
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
		/*Eigen::MatrixXd outputTriangleMeshVertices;
		Eigen::MatrixXi outputTriangleMeshFaces;
		loadFromMesh(m_triMesh, outputTriangleMeshVertices, outputTriangleMeshFaces);
		outputTriangleMeshVertices.transposeInPlace();
		outputTriangleMeshFaces.transposeInPlace();

		igl::writeOBJ(outputMeshName, outputTriangleMeshVertices, outputTriangleMeshFaces);*/
		m_triMesh->writeFile_tinyobj(outputMeshName);

	}
}


bool SOPeriodicBendingModeShellStructure::loadAllStrainSpaceModes(const string& strainSpaceModesFolder)
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


void SOPeriodicBendingModeShellStructure::step()
{
	
}

void SOPeriodicBendingModeShellStructure::renderImGUIViewerWindow()
{
	ImGui::Begin("Mesh Window");
	{
		if (ImGui::Button("Reset simulation"))
		{
			reset();
		}
		
		ImGui::Checkbox("Draw periodic patches", &m_drawPeriodicPathes);
		if (ImGui::Button("Ouput periodic patches"))
		{
			std::vector<Eigen::MatrixXd> vertices_vector;
			std::vector<Eigen::MatrixXi> faces_vector;
			computePeriodicMeshes(vertices_vector, faces_vector);

			Eigen::MatrixXd combined_vertices;
			Eigen::MatrixXi combined_faces;
			TriMesh::combineMeshes(vertices_vector, faces_vector,
				combined_vertices, combined_faces);
			igl::writeOBJ(m_outputPath + "/periodicMesh.obj", combined_vertices.transpose(), combined_faces.transpose());
		}
		ImGui::Checkbox("Draw periodic pair", &m_drawPeriodicTranslationPair);
		ImGui::InputInt("Periodic pair index", &m_periodicPairIndex);

		ImGui::Checkbox("Draw periodic bending element", &m_drawPeriodicBendingElements);
		ImGui::InputInt("Periodic bending element index", &m_periodicBendingElementIndex);
		

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


	ImGui::Begin("Iteractive Pinning");

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
				m_selectedPairs.erase(m_selectedPairs.begin() + m_selectedPairs.size() -1);
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
			m_hingeElements[i].setRestAngle(0.0);
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

	ImGui::End();
}

void SOPeriodicBendingModeShellStructure::computePeriodicMeshes(std::vector<Eigen::MatrixXd>& vertices_vector, std::vector<Eigen::MatrixXi>& faces_vector)
{
	int topVertexIndex = m_periodicReferenceTopBottom.first;
	int bottomVertexIndex = m_periodicReferenceTopBottom.second;

	int leftVertexIndex = m_periodicReferenceLeftRight.first;
	int rightVertexIndex = m_periodicReferenceLeftRight.second;

	Eigen::Vector3d topToBottomTranslation = visualize_m_x.segment<3>(bottomVertexIndex * 3) - visualize_m_x.segment<3>(topVertexIndex * 3);
	Eigen::Vector3d leftToRightTranslation = visualize_m_x.segment<3>(rightVertexIndex * 3) - visualize_m_x.segment<3>(leftVertexIndex * 3);

	Eigen::MatrixXi faces = Eigen::Map<Eigen::MatrixXi>(m_triMesh->faces().data(), 3, m_triMesh->faces().size() / 3);

	for (int i = -3; i <= 3; i++)
	{
		for (int j = -3; j <= 3; j++)
		{
			//if (i == 0 && j == 0)
			//	continue;
			Eigen::VectorXd transVertices = visualize_m_x;
			Eigen::MatrixXd transVertices_xd = Eigen::Map<Eigen::MatrixXd>(transVertices.data(), 3, transVertices.size() / 3);
			transVertices_xd.colwise() += (i * topToBottomTranslation + j * leftToRightTranslation);


			vertices_vector.push_back(transVertices_xd);
			faces_vector.push_back(faces);
		}
	}
}

void SOPeriodicBendingModeShellStructure::render()
{
	{
		if (m_isRenderStrainColor)
		{ 
			double maximumStrain = 0.7;
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


			if (m_drawPeriodicPathes
				&& m_periodicReferenceTopBottom.first != -1 && m_periodicReferenceTopBottom.second != -1
				&& m_periodicReferenceLeftRight.first != -1 && m_periodicReferenceLeftRight.second != -1)
			{//draw periodic patches
				int topVertexIndex = m_periodicReferenceTopBottom.first;
				int bottomVertexIndex = m_periodicReferenceTopBottom.second;

				int leftVertexIndex = m_periodicReferenceLeftRight.first;
				int rightVertexIndex = m_periodicReferenceLeftRight.second;

				Eigen::Vector3d topToBottomTranslation = visualize_m_x.segment<3>(bottomVertexIndex * 3) - visualize_m_x.segment<3>(topVertexIndex * 3);
				Eigen::Vector3d leftToRightTranslation = visualize_m_x.segment<3>(rightVertexIndex * 3) - visualize_m_x.segment<3>(leftVertexIndex * 3);

				for (int i = -1; i <= 1; i++)
				{
					for (int j = -1; j <= 1; j++)
					{
						if (i == 0 && j == 0)
							continue;
						Eigen::VectorXd transVertices = visualize_m_x;
						Eigen::MatrixXd transVertices_xd = Eigen::Map<Eigen::MatrixXd>(transVertices.data(), 3, transVertices.size() / 3);
						transVertices_xd.colwise() += (i * topToBottomTranslation + j * leftToRightTranslation);

						m_meshDrawer->updateVertices(transVertices_xd);

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


	if (m_drawPeriodicBendingElements)
	{
		if (m_periodicBendingElementIndex >= 0 && m_periodicBendingElementIndex < m_periodicDSHingeFixedRestPoseElement.size())
		{
			const auto& dofIndices = m_periodicDSHingeFixedRestPoseElement[m_periodicBendingElementIndex].getDofIndices();

			int tri0_v0 = dofIndices[0] / 3;
			int tri0_v1 = dofIndices[3] / 3;
			int tri0_v2 = dofIndices[6] / 3;

			int tri1_v0 = dofIndices[9] / 3;
			int tri1_v1 = dofIndices[12] / 3;
			int tri1_v2 = dofIndices[15] / 3;

			Eigen::Vector4d color(1.0, 0.0, 0.0, 1.0);
			m_primitiveDrawer->drawTriangle(visualize_m_x.segment<3>(tri0_v0*3), visualize_m_x.segment<3>(tri0_v1 * 3), visualize_m_x.segment<3>(tri0_v2 * 3), color);
			m_primitiveDrawer->drawTriangle(visualize_m_x.segment<3>(tri1_v0 * 3), visualize_m_x.segment<3>(tri1_v1 * 3), visualize_m_x.segment<3>(tri1_v2 * 3), color);
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

void SOPeriodicBendingModeShellStructure::reset()
{
	time_passed = 0.0;
	Eigen::VectorXd s;
	reset(s, false, false);
}

void SOPeriodicBendingModeShellStructure::reset(Eigen::VectorXd& new_status)
{
	time_passed = 0.0;
	reset(new_status, false, false);
}

void SOPeriodicBendingModeShellStructure::reset(Eigen::VectorXd& new_status, bool build_graph, bool ls)
{
	SODeformableObject::reset();

	//we need fully reset all parts of the robot
	time_passed = 0.0;
	/*{
	printf("reInit:\n");
	printf("x[8]: %f%f%f\n", m_x[3 * 8], m_x[3 * 8 + 1], m_x[3 * 8 + 2]);
	printf("\n");
	}*/

	m_V.setZero();

	//if the initial configuration has changed
	//we need to reinitialize all the elements, but now we only change the motor speed
	reInit();
}

void SOPeriodicBendingModeShellStructure::reInit()
{
	//m_edges.clear();
	//computeEdges();

	{

		m_V.setZero();			 //!!!!!!now this velocity is zero, TODO: parameters
	}
	//m_mesh = new DSMesh;
	//m_mesh->init(m_triMesh);

	m_faceColors = m_triMesh->faceColors();


	m_v.setZero();
	//m_mass.setZero();
	m_appliedForces.setZero();

	setVisualizationConfiguration(m_x);
}
