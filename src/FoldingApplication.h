#pragma once
#ifndef FOLDING_APPLICATION
#define FOLDING_APPLICATION

#include <SimLib/Application.h>

#include <Eigen/Core>

#include "foldingShellMeshStructure.h"

class FoldingApplication :
    public Application
{
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	FoldingApplication();
	FoldingApplication(int argc, char* argv[]);
    ~FoldingApplication();

	void renderSceneObjects();	

	void renderShadowMap()
	{
		renderSceneObjects();
	}
protected:
	void render() override;
	void step();
	void runSimulation() override;
	void ImGuiSetting() override;
	bool init(int argc, char* argv[]) override;
	void beforeExit() override;
	void keyPressed(int key, int mod) override;
	virtual void mouseButtonPressed(int button, int mods, double xPos, double yPos) override;
private:
	

	void toggleOptimization(void);

private:
	bool m_running;
	bool m_saveGif;
	uint32_t GIFDelay;
	
	foldingShellMeshStructure* m_origamiShellMeshStructure_instance;
};

#endif