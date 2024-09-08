#ifndef SO_BENDING_MODES_SHELL_STRUCTURE_H
#define SO_BENDING_MODES_SHELL_STRUCTURE_H

#if _MSC_VER > 1000
#pragma once
#endif

#include <SimLib/SODeformableObject.h>
#include <SimLib/Elements/FixDoFElement.h>

#include "Elements/CSTShellElement.h"
#include "Elements/DSHingeFixedRestPoseElement.h"
#include "Elements/koiterShapeOperatorRestPoseBendingElement.h"

#include "./Elements/pairAttachmentElement.h"
#include <SimLib/Collision/AABBTree.h>
#include <SimLib/Core/UnconstrainedNewton.h>
#include <SimLib/Geom/SimpleTriMesh.h>

#include <ipc/ipc.hpp>

#include "tbb/tbb.h"



class SOBendingModeShellStructure :
    public SODeformableObject
{
public:
	SOBendingModeShellStructure();
	SOBendingModeShellStructure(const SOBendingModeShellStructure& shell);
	~SOBendingModeShellStructure();

	virtual void init(const TriMesh* pTriMesh, SOMaterial* mat, const TriMesh* targetpTriMesh,
		double thickness,
		const std::vector<std::pair<int, int>>& fixVertexOnLinePairs,
		const std::vector<std::pair<int, int>>& fixingVerticesDofs,
		const std::vector<std::pair<int, Eigen::Vector3d>>& vertexIndex_force,
		bool useShapeOperator,
		string outputPath);
	void initIPCConfig();
	virtual void reconstructMassStuff(const Eigen::VectorXd& vX);
	virtual void initElements(const Eigen::VectorXd& vX);
	virtual void reconstructElementVector();

	TriMesh* getTriMesh() { return m_triMesh; }
	const TriMesh* getInitialTriMesh() { return initial_mesh; }
	const SOMaterial* getMaterial() { return m_material; }
	double get_cstThickness() { return m_cstThickness; }
	double get_rho() { return m_material->rho; }

	void setVisualizationConfiguration(const Eigen::VectorXd& vx) {
		visualize_m_x = vx;
		m_triMesh->vertices() = vx;
		m_triMesh->setDirty(true);
		m_meshDrawer->setMesh(m_triMesh);
	}
	Eigen::VectorXd getVisualizationConfiguration() {
		//if (m_visualize_m_x.size() != 0)
		//	return m_visualize_m_x;
		//else
		return m_x;
	}
	TriMesh* outputMesh() { return m_triMesh; }

	//compute mesh stuff
	void computeEdges();
	void computeMasses(const Eigen::VectorXd& vX);


	bool staticSolve_newton(Eigen::VectorXd& vx, const Eigen::VectorXd& vX) override;
	virtual void computeLinearModes();
	virtual void convertToStrainSpaceMode(const Eigen::VectorXd& positionSpaceMode, Eigen::VectorXd& strainSpaceMode);
	void convertToShapeOperatorStrainSpaceMode(const Eigen::VectorXd& positionSpaceMode, std::vector<Eigen::Matrix2d>& strainSpaceMode);

	virtual void animateStrainSpaceMode();
	virtual void projectToManifold();
	void computeBatchStrainSpaceMode();
	virtual void computeStrainSpaceModes();
	virtual void blendStrainSpaceModes(const std::vector<int>& selectedModes, const std::vector<double>& eachModeLength, double stepLength);
	
	void computeStrainSpaceSensitiveMatrix(const Eigen::VectorXd& vx, const Eigen::VectorXd& vX, SpMat &dxdtheta_rest);

	virtual void step();
	virtual void renderImGUIViewerWindow() override;

	void ouputNonlinearCompliantModesAtIndex(int index);
	void ouputRotationStrainSpaceModesAtIndex(int index);

	virtual void ouputStrainSpaceModesAtIndex(int index);
	void outputAllStrainSpaceModes();
	bool loadStrainSpaceMeshFolder(const std::string& testSetFolder, std::vector<std::pair<std::string, double>>& fileName_length);
	bool loadFolderMeshes(const string& strainSpaceMeshFolder, std::vector<std::pair<TriMesh, double>>& meshSet_length);
	bool loadAllNonlinearCompliantModes(const string& strainSpaceModesFolder);
	virtual bool loadAllStrainSpaceModes(const string& strainSpaceModesFolder);

	virtual void computeMaxStrainFacesColorMap(const double colorMap_UpperBound, const Eigen::VectorXd& vx, const Eigen::VectorXd& vX, Eigen::MatrixXd& facesColor);
	virtual void render() override;
	virtual void reInit();
	virtual void reset();
	virtual void reset(Eigen::VectorXd& new_status);
	virtual void reset(Eigen::VectorXd& new_status, bool build_graph, bool ls);
	void mouseLeftClick(const Eigen::Vector3d& src, const Eigen::Vector3d& ray);

	Eigen::VectorXd optimizeOneTimeFoldingParametersWithinSphereWithRigidTranformation(const Eigen::VectorXd& initialParameters, const Eigen::MatrixXd& strainSpaceModes, Eigen::Vector3d& sphere_center, double sphere_radius, std::vector<Eigen::VectorXd>& intermediateStates);
	Eigen::VectorXd optimizePositionsToTargetIPCSphere(Eigen::VectorXd& vX, Eigen::Vector3d& sphere_center, double sphere_radius, std::vector<Eigen::VectorXd>& intermediateStates);

	//mass
	virtual Eigen::VectorXd applyMassMatrix(const Eigen::VectorXd& v) override;
	virtual Eigen::VectorXd applyInvMassMatrix(const Eigen::VectorXd& v) override;
	virtual void addMassMatrix(Eigen::SparseMatrix<double>& mat) override;
	virtual SpMat computeMassMatrix() override;
	virtual Eigen::VectorXd computeGravityForce() override;

	virtual int getNumPoints()  const override { return m_nv; }

public:
	ipc::CollisionMesh m_ipc_collisionMesh;
	double m_ipc_adaptiveStiffness;
	double m_ipc_minDistance;
	double m_ipc_upperBoundStiffness;
	double m_ipc_dhat;


	int m_nv;
	Eigen::VectorXd m_mass;
	Eigen::VectorXd m_gravity;
	std::vector<std::pair<int, int>> m_fixVertexOnLinePairs;
	std::vector<int> m_fixingDofs;
	std::vector<std::pair<int, Eigen::Vector3d>> m_vertexIndex_force;

	double m_fix_stiffness;
	//initial mesh
	const TriMesh* initial_mesh;
	//mesh
	TriMesh* m_triMesh;
	TriMesh* m_targetTriMesh;
	bool m_isRenderStrainColor;

	std::vector<int> m_edges;
	Eigen::VectorXd m_faceColors;

	//body parameterization

	double m_cstThickness;
	bool m_useShapeOperator;
	std::vector<CSTShellElement> m_cstShellElements;
	std::vector<DSHingeFixedRestPoseElement>	m_hingeElements;
	std::vector<koiterShapeOperatorRestPoseBendingElement> m_koiterShapeOperatorRestPoseBendingElements;
	std::vector<FixDoFElement>					m_fixDoFElements;

	std::vector<pairAttachmentElement>			m_pairAttachmentElements;
	//
	bool USE_MAKEPD;
	bool USE_IPC;
	bool m_useStaticSolver;
	

	//visualization
	bool m_drawWireframe;
	bool m_useBothWireframeAndmesh;
	std::vector<TriMesh> m_simulationMeshes;
	bool m_drawFixingVertices;
	bool m_drawHardEdges;

	bool m_drawMaxCurvature;
	Eigen::VectorXd m_vertexMaxCurvature;

	bool m_updateModesToVisualization;

	//interactive 
	bool m_drawInteractivePairs;
	ipc::Constraints m_pinningPairConstraint_set;
	std::vector<Eigen::Vector4d> m_pinningPairConstraintColor;

	bool m_enableInteractiveSelection;
	std::vector<std::pair<int, int>> m_selectedPairs;
	std::vector<Eigen::Vector4d> m_selectedPairsColor;
	Eigen::VectorXd m_selectMeshVertices;


	std::vector<std::pair<double, Eigen::VectorXd>> eigen_modal_analysis;//nature modes
	//mode animation
	double time_passed;
	bool m_freeLinearModesComputed;
	bool m_strainSpaceModesComputed;
	bool m_animateModes;
	double m_modesAmplitude;
	int m_modesAnimationID;


	bool m_useLinearModes;
	bool m_useStrainSpaceModes;
	bool m_useDevelopableStrainSpaceModes;

	//compute with length
	double m_computeModeLength;
	double m_computeModeStepLength;
	//strain space modes
	std::vector<std::vector<Eigen::Matrix2d>> m_shapeOperatorStrainSpaceModes;
	std::vector<Eigen::VectorXd> m_strainSpaceModes;//each linear mode corresponds to a strain space mode
	std::map<int, std::vector<std::pair<double, Eigen::VectorXd>>> m_strainSpaceModes_index_to_lengthState;
	std::map<std::vector<int>, std::pair<double,std::map<std::vector<int>, Eigen::VectorXd>>> m_selectedStrainSpce_to_eachModeLengthState;
	std::map<std::vector<int>, std::vector<std::pair<Eigen::VectorXd, Eigen::VectorXd>>> m_selectedStrainSpce_to_modeLengthState;
	//nonlinear compliant modes
	std::map<int, std::vector<std::pair<double, Eigen::VectorXd>>> m_nonlinearCompliantModes_index_to_lengthState;
	//rotation strain space modess
	bool m_rotationStrainSpaceModesComputed;
	SpMat m_G;//discrete deformation gradient operator
	SpMat m_Vol;//volume matrix
	std::vector<Eigen::VectorXd> m_rotationStrainSpaceModes; //mode -> tet element contains (3 rotation tensor, 6 stain tensor); see "Interactive editing of deformable simulations"
	std::map<int, std::vector<std::pair<double, Eigen::VectorXd>>> m_rotationStrainSpaceModes_index_to_lengthState;

	SimOpt::AABBTree *m_selectMeshAABBTree;

	//sensitivity analysis
	bool m_visualizeSensitivityExploration;
	SpMat m_dxdtheta_rest;
	bool m_sensitivityMatrixUpdated;
	Eigen::VectorXd m_sensitivityAnalysisState;
	bool m_sensitivityExplorationProjected;
	Eigen::VectorXd m_sensitivityAnalysisProjectedState;


	//target shape optimization
	bool m_visualizeTargetFoldingShape;
	Eigen::VectorXd m_modesCoefficients;
	bool m_showTargetSphere;
	Eigen::Vector3d m_targetSphereCenter;
	double m_targetSphereRadius;
	std::vector<Eigen::VectorXd> m_statesAlongModesCoefficients;

	//output
	string m_outputPath;
};


#endif