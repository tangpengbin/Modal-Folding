#include "FoldingApplication.h"
//basic headers
#include <stdio.h>
#include <Eigen/Geometry>
#include <iostream>

#include "imgui.h"
#include "backends/imgui_impl_glfw.h"
#include "backends/imgui_impl_opengl3.h"
//
#include <string>
#include <ctime>
#include <algorithm>

#include <spdlog/spdlog.h>
#include <CLI/CLI.hpp>

struct CLIArgs {
    int progMode;
    std::string rodStructureFileName;
    int logLevel = 0; // trace (all logs)
};


void init_cli_app(CLI::App &app, CLIArgs &args)
{
    app.add_set("progMode", args.progMode, { 0, 1, 2, 3 },
        "program mode 0 = simluation, 1 = simulation(save gif), 2 = simulation(output simulation files), 3 = simulation(output simulation files and gif),",
        true);
    app.add_option("inputRodStructureFileName", args.rodStructureFileName, "input rod file name", true);
    app.add_set("--logLevel,--log", args.logLevel, { 0, 1, 2, 3, 4, 5, 6 },
        "set log level (0=trace, 1=debug, 2=info, 3=warn, 4=error, 5=critical,"
        " 6=off)",
        true);
    
}

FoldingApplication::FoldingApplication()
{
    //m_periodicStructure_instance = nullptr;
    m_running = false;
    m_saveGif = false;
    GIFDelay = 10; //10ms
    //Application::Application(argv);
}


FoldingApplication::FoldingApplication(int argc, char* argv[]) : Application(argc, argv)
{
    //m_periodicStructure_instance = nullptr;
    m_running = false;
    m_saveGif = false;
    GIFDelay = 10; //10ms
    //Application::Application(argv);
}

FoldingApplication::~FoldingApplication()
{
    
}

void FoldingApplication::renderSceneObjects()
{
    if(m_origamiShellMeshStructure_instance !=nullptr)
        m_origamiShellMeshStructure_instance->render();
}

void FoldingApplication::render()
{
	renderSceneObjects();
}

void FoldingApplication::step()
{
   /* if (m_origamiShellMeshStructure_instance != nullptr)
    {
        if (!m_origamiShellMeshStructure_instance->step())
        {
            m_running = false;
            spdlog::error("step failed");
            return;
        }
    }
*/
}


void FoldingApplication::runSimulation()
{
    if (m_running)
    {
        step();
    }
}

void FoldingApplication::ImGuiSetting()
{
    Application::ImGuiSetting();

    if(m_origamiShellMeshStructure_instance != nullptr)
        m_origamiShellMeshStructure_instance->renderImguiViewer();
}

bool FoldingApplication::init(int argc, char* argv[])
{
#ifdef USE_TBB
#ifdef TBB_NUM_THREADS
    tbb::global_control control(
        tbb::global_control::max_allowed_parallelism, 1);
    //tbb::task_scheduler_init init(TBB_NUM_THREADS);
#endif
#endif
    CLI::App app{ "IPC" };
    CLIArgs args;
    init_cli_app(app, args);
    try {
        app.parse(argc, argv);
    }
    catch (const CLI::ParseError& e) {
        return app.exit(e);
    }

    //args.logLevel = 6;
    spdlog::set_level(static_cast<spdlog::level::level_enum>(args.logLevel));// SPDLOG_LEVEL_OFF
    if (args.logLevel == 6) {
        std::cout.setstate(std::ios_base::failbit);
    }

    //init rod structure
    if (args.rodStructureFileName == "")
        return false;
    
    bool outputRunningStates = true;
    //set simulaton mode
    switch (args.progMode) {
    case 0:
        // simulation mode
        m_running = false;
        spdlog::info("Simulation mode");
        break;

    case 1:
        // autoswitch optimization mode
        m_running = false;
        m_saveGif = true;
        spdlog::info("Simulation mode (save png)");
        break;
    case 2:
        m_running = false;
        outputRunningStates = true;
        spdlog::info("Simulation mode (output simulation files)");
        break;
    case 3:
        m_running = false;
        m_saveGif = false;
        outputRunningStates = true;
        spdlog::info("Simulation mode (output simulation files and don't save png)");
        break;
    default:
        spdlog::error("No progMode {:d}", args.progMode);
        return false;
    }
        
    //return m_periodicStructure_instance.init(args.inputFileName, args.rodStructureFileName.c_str(), &m_meshDrawer, &m_primitiveDrawer, outputRunningStates);

   /*DIMStructure* mesh_instance = new DIMStructure;
   m_homogenizedMeshStructure_instance = mesh_instance;
   if (m_homogenizedMeshStructure_instance != nullptr)
       return m_homogenizedMeshStructure_instance->init(args.rodStructureFileName.c_str(), &m_meshDrawer, &m_otherMeshDrawer, &m_primitiveDrawer, &m_rodDrawer, outputRunningStates);
   else
       return false;*/


    m_origamiShellMeshStructure_instance = new foldingShellMeshStructure;
    if (m_origamiShellMeshStructure_instance != nullptr)
    {
        if(!m_origamiShellMeshStructure_instance->init(args.rodStructureFileName.c_str(), &m_primitiveDrawer, &m_meshDrawer, outputRunningStates))
            return false;
    }
    else
        return false;

    Application::init(argc, argv);
   

}

void FoldingApplication::toggleOptimization(void)
{
    m_running = !m_running;
    if (m_running)
    {
        spdlog::info("start/resume optimization, press again to pause.");
    }
    else
    {
        spdlog::info("pause optimization, press again to resume.");
    }
}

void FoldingApplication::beforeExit()
{
}

void FoldingApplication::keyPressed(int key, int mod)
{
    Application::keyPressed(key, mod);

    switch (key)
    {
    case 'r':
    case 'R':
        //resetSimulation();
        spdlog::info("pressed R or r");
        break;
    case ' ':
        spdlog::info("Optimize for one step");
        step();
        break;
    case GLFW_KEY_ENTER://return
        spdlog::info("IPC optimize full simulation");
        toggleOptimization();

        break;
    default:
        spdlog::info("pressed key not been set yet: " + key);
        break;
    }

}

void FoldingApplication::mouseButtonPressed(int button, int mods, double xPos, double yPos)
{
    Application::mouseButtonPressed(button, mods, xPos, yPos);

    if (button == GLFW_MOUSE_BUTTON_LEFT)
    {
        bool translation = mods & GLFW_MOD_SHIFT;
        Eigen::Vector3d src, ray;
        getCamera()->computeRay(xPos,yPos, src, ray);

        //m_origamiShellMeshStructure_instance->mouseLeftClick(src, ray);
    }
}
