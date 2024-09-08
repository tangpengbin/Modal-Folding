#pragma once
#ifndef HOMOGENIZED_MESH_STRUCTURE_H
#define HOMOGENIZED_MESH_STRUCTURE_H
#include <string>

//simopt
#include "./Core/SOTypes.h"
#include "./Core/SOUtils.h"
#include <nlohmann/json.hpp>
#include "SOBendingModeShellStructure.h"

class foldingShellMeshStructure
{
public:
	foldingShellMeshStructure();
	~foldingShellMeshStructure();

public:
	//this use interlock_structure_szFileName to initialize rod structure, and ipc will use this structure to initialize itself
	bool init(const std::string& structure_szFileName, PrimitiveDrawer* p_primitiveDrawer, meshDrawer* p_meshDrawer, bool m_outputRunningStates = true);
	void render();
	void renderImguiViewer();
	bool step();
	void reset();
	double getTimeStep() { return m_timeStep; }
	void mouseLeftClick(const Eigen::Vector3d& src, const Eigen::Vector3d& ray);

private://functions for simopt
	void homogenizedMeshConfigInit(const std::string& filename, PrimitiveDrawer* input_primitiveDrawer, meshDrawer* input_meshDrawer);
	void setHomogenizedMeshStructure();

	void loadHomogenizedMesh(nlohmann::json& j);
	void loadConfig(nlohmann::json& j);
	void loadHomogenizedMeshStructureConfig(const std::string& filename);


private://variables for simopt
	//SOShellStructure m_shellStructure;
	SOBendingModeShellStructure* m_shellStructure;
	bool m_useShapeOperatorShell;
	//SOProjectiveDynamicsInterlockStructure m_shellStructure;
	SpMat m_interlockStructure_massMatrix;
	//for the initialization of mesh
	string m_meshObj_str;
	string m_targetMeshObj_str;
	string m_deformationLimits_str;
	string m_bendingDeformationLimitsCurvatureEnergyDensity_str;
	string m_collisionMeshObj_str;
	string m_strainSpaceLimitingModel_str;
	string m_bendingStainSpaceLimitingModel_str;
	bool m_isPeriodic;
	Eigen::Vector3d m_periodicLeftToRight;
	Eigen::Vector3d m_periodicTopToBottom;
	bool m_isReflection;
	Eigen::Vector4d m_fourReflectionPlane; //x1,y1,x2,y2

	Eigen::VectorXd m_interlockStructure_lastState;
	int m_numMassPoints;
	double m_density;
	std::vector<double> m_crossSectionRadius;//this is the original data from jason
	double m_youngsModulus;
	double m_poissonsRatio;
	double m_thickness;
	std::vector<std::pair<int, int>> fixVertexDofs;
	std::vector<std::pair<int, Eigen::Vector3d>> vertexIndex_force;
	std::vector<std::pair<int, int>> m_fixVertexOnLinePairs;
	
private://ipc private variables
	double m_timeStep = 0.005;
	int simulationId;
	double m_self_friction = 0.1, m_dHat = 1e-3, m_epsv = 1e-3;
	bool m_turnOffGravity = false;
	string m_outputPath;


	Eigen::MatrixXd m_rendering_verticesRest;//each column is a vertex
	Eigen::MatrixXd m_rendering_vertices;//each column is a vertex
	Eigen::MatrixXi m_rendering_edges;//each column is an edge

	
	meshDrawer* m_meshDrawer;
	PrimitiveDrawer *m_primitiveDrawer;
};

#endif // HOMOGENIZED_MESH_STRUCTURE_H