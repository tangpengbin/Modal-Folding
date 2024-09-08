#ifndef SO_REFLECTION_BENDING_MODES_SHELL_STRUCTURE_H
#define SO_REFLECTION_BENDING_MODES_SHELL_STRUCTURE_H

#if _MSC_VER > 1000
#pragma once
#endif
#include "SOBendingModeShellStructure.h"
#include "./Elements/reflectionDSHingeFixedRestPoseElement.h"

class SOReflectionBendingModeShellStructure :
    public SOBendingModeShellStructure
{
public:
	SOReflectionBendingModeShellStructure();
	SOReflectionBendingModeShellStructure(const SOReflectionBendingModeShellStructure& shell);
	~SOReflectionBendingModeShellStructure();

	virtual void init(const TriMesh* pTriMesh, SOMaterial* mat, double thickness,
		const std::vector<std::pair<int, int>>& fixVertexOnLinePairs,
		const std::vector<std::pair<int, int>>& fixingVerticesDofs,
		const std::vector<std::pair<int, Eigen::Vector3d>>& vertexIndex_force,
		const Eigen::Vector4d& fourReflectionPlane_rest,
		string outputPath);

	void initIPCConfig();
	void processSystemDoFs();
	Eigen::VectorXd convertAllDofsToSystemDofs(const Eigen::VectorXd& allDofs);
	Eigen::VectorXd convertSystemDofsToAllDofs(const Eigen::VectorXd& systemDofs);
	SpMat convertSystemDofsToAllDofsJacobian(const Eigen::VectorXd& systemDofs);

	virtual void initElements(const Eigen::VectorXd& vX) override;
	virtual void reconstructElementVector() override;



	//compute mesh stuff
	void computeMasses(const Eigen::VectorXd& vX) override;
	void reconstructMassStuff(const Eigen::VectorXd& vX) override;

	void computeHessian(const Eigen::VectorXd& vx, const Eigen::VectorXd& vX, SpMat& A) override;
	bool staticSolve_newton(Eigen::VectorXd& vx, const Eigen::VectorXd& vX) override;
	void projectToManifold() override;
	virtual void convertToStrainSpaceMode(const Eigen::VectorXd& positionSpaceMode, Eigen::VectorXd& strainSpaceMode) override;
	virtual void animateStrainSpaceMode() override;
	virtual void computeStrainSpaceModes() override;

	virtual void ouputStrainSpaceModesAtIndex(int index);
	bool loadAllStrainSpaceModes(const string& strainSpaceModesFolder);

	void computeReflectionMeshes(int startRows, int numRows, int startCols, int numCols, std::vector<Eigen::MatrixXd>& vertices_vector, std::vector<Eigen::MatrixXi>& faces_vector, std::vector<bool>& flipNorm_vector);

	virtual void renderImGUIViewerWindow() override;
	virtual void render() override;


	void setVisualizationConfiguration(const Eigen::VectorXd& vx) override {
		visualize_m_x = convertSystemDofsToAllDofs(vx);
		m_triMesh->vertices() = visualize_m_x;
		m_triMesh->setDirty(true);
		m_meshDrawer->setMesh(m_triMesh);
	}

public:
	std::vector<reflectionDSHingeFixedRestPoseElement> m_reflectionDSHingeFixedRestPoseElement;

	//map from realvertex index to system dof index and number of dofs
	enum class BoundaryVertexType
	{
		topleft, top, topright, left, right, bottomleft, bottom, bottomright
	};
	std::map<int, BoundaryVertexType> boundaryVertices;
	std::map<int, std::pair<int, int>> m_realVertex_to_systemDofAndNumDofs;
	std::map<int, std::pair<int, int>> m_systemDof_to_realVertexAndNumDofs;
	int m_numSystemDofs;

	//   -----y1----
	//	 |			|
	//	 x1			x2
	//	 |			|
	//	 -----y2----
	Eigen::Vector4d m_fourReflectionPlane_rest;//with order x1,y1,x2,y2
	Eigen::Vector4i m_fourCornerIndices;//topleft, topright, bottomleft, bottomright

	int m_numExtraDofs;
	bool m_drawReflectionPatches;
};


#endif