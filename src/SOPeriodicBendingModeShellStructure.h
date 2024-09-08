#ifndef SO_PERIODIC_BENDING_MODES_SHELL_STRUCTURE_H
#define SO_PERIODIC_BENDING_MODES_SHELL_STRUCTURE_H

#if _MSC_VER > 1000
#pragma once
#endif
#include "SOBendingModeShellStructure.h"
#include "./Elements/periodicDSHingeFixedRestPoseElement.h"
#include "./Elements/periodicCSTShellElement.h"
#include "./Elements/targetPlanarParametersElement.h"

class SOPeriodicBendingModeShellStructure :
    public SOBendingModeShellStructure
{
public:
	SOPeriodicBendingModeShellStructure();
	SOPeriodicBendingModeShellStructure(const SOPeriodicBendingModeShellStructure& shell);
	~SOPeriodicBendingModeShellStructure();

	virtual void init(const TriMesh* pTriMesh, SOMaterial* mat, double thickness,
		const std::vector<std::pair<int, int>>& fixVertexOnLinePairs,
		const std::vector<std::pair<int, int>>& fixingVerticesDofs,
		const std::vector<std::pair<int, Eigen::Vector3d>>& vertexIndex_force,
		/*if the structure is periodic, we can set periodic translations. (if not these two values will ignore them)*/
		const Eigen::Vector3d& periodicLeftToRight,
		const Eigen::Vector3d& periodicTopToBottom,
		string outputPath);
	void initIPCConfig();
	void processSystemDoFs();
	Eigen::VectorXd convertAllDofsToSystemDofs(const Eigen::VectorXd& allDofs);
	Eigen::VectorXd convertSystemDofsToAllDofs(const Eigen::VectorXd& systemDofs);
	SpMat convertSystemDofsToAllDofsJacobian(const Eigen::VectorXd& systemDofs);
	void computePeriodicity();
	void reconstructMassStuff(const Eigen::VectorXd& vX) override;
	void initElements(const Eigen::VectorXd& vX) override;
	void reconstructElementVector() override;


	//compute mesh stuff
	void computeMasses(const Eigen::VectorXd& vX) override;

	bool staticSolve_newton(Eigen::VectorXd& vx, const Eigen::VectorXd& vX) override;
	void convertToStrainSpaceMode(const Eigen::VectorXd& positionSpaceMode, Eigen::VectorXd& strainSpaceMode) override;
	void animateStrainSpaceMode() override;
	void projectToManifold() override;
	void computeStrainSpaceModes() override;

	virtual void ouputStrainSpaceModesAtIndex(int index);
	bool loadAllStrainSpaceModes(const string& strainSpaceModesFolder);

	void step() override;

	virtual void renderImGUIViewerWindow() override;

	void computePeriodicMeshes(std::vector<Eigen::MatrixXd>& vertices_vector, std::vector<Eigen::MatrixXi>& faces_vector);

	virtual void render() override;
	virtual void reInit();
	virtual void reset();
	virtual void reset(Eigen::VectorXd& new_status);
	virtual void reset(Eigen::VectorXd& new_status, bool build_graph, bool ls);


	void setVisualizationConfiguration(const Eigen::VectorXd& vx) override{
		visualize_m_x = convertSystemDofsToAllDofs(vx);
		m_triMesh->vertices() = visualize_m_x;
		m_triMesh->setDirty(true);
		m_meshDrawer->setMesh(m_triMesh);
	}

public:

	//periodic boundary conditions
	bool m_usePeriodicBoundaryCondition;

	Eigen::Vector3d m_periodicLeftToRight_rest;
	Eigen::Vector3d m_periodicTopToBottom_rest;
	std::pair<int, int> m_periodicReferenceTopBottom;
	std::pair<int, int> m_periodicReferenceLeftRight;

	std::map<int, std::vector<std::pair<int, Eigen::Vector3d*>>> m_leftOrTopIndex_to_rightOrBottomWithTrans; //use rest translations (contains topelft, topright, bottomleft)
	enum class periodicBoundaryVertexType
	{
		topleft, top, topright, left, right, bottomleft, bottom, bottomright
	};
	std::map<int, periodicBoundaryVertexType> m_allPeriodicVertices;//detail periodic boundary information for alldof vertices
	Eigen::Vector4i m_cornerVerticesIndices;//this is shortcut to obtain four corners in order of (topleft, topright, bottomleft, bottomright)
	std::map<int, int> m_system_to_allDofs;
	std::map<int, int> m_all_to_systemDofs;
	int m_numExtraDofs;

	std::vector<periodicCSTShellElement> m_periodicCSTShellElement;
	std::vector<periodicDSHingeFixedRestPoseElement> m_periodicDSHingeFixedRestPoseElement;
	std::vector<targetPlanarParametersElement> m_targetPlanarParametersElement;

	//visualization
	bool m_drawPeriodicPathes;
	bool m_drawPeriodicTranslationPair;
	int m_periodicPairIndex;

	bool m_drawPeriodicBendingElements;
	int m_periodicBendingElementIndex;
};


#endif