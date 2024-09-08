#include "foldingShellMeshStructure.h"
//system headers
#include <fstream>
#include <spdlog/spdlog.h>
#ifdef _WIN32
#include <direct.h>
#elif __APPLE__ || __linux__
//#include <sys/stat.h>
#include <unistd.h>
#endif
#include "imgui.h"

#include <sys/types.h>
#include <sys/stat.h>
#include <string>
#include <fstream>

#include <igl/readOBJ.h>

//simopt headers
#include "Core/File.h"
#include "Geom/MeshObj.h"

#include "SOPeriodicBendingModeShellStructure.h"
#include "SOReflectionBendingModeShellStructure.h"

foldingShellMeshStructure::foldingShellMeshStructure()
{
	m_numMassPoints = 0;
    m_shellStructure = nullptr;
    m_isPeriodic = false;
    
    m_useShapeOperatorShell = false;

    
    m_outputPath = "./output";
    struct stat st_2 = { 0 };
        if (stat(m_outputPath.c_str(), &st_2) == -1) {
#ifdef _WIN32
            mkdir(m_outputPath.c_str());
#elif __APPLE__ || __linux__ 
            mkdir(m_outputPath.c_str(),0777);
#endif
        }
}

foldingShellMeshStructure::~foldingShellMeshStructure()
{
    if (m_shellStructure != nullptr)
        delete m_shellStructure;
    m_shellStructure = nullptr;

}

void foldingShellMeshStructure::homogenizedMeshConfigInit(const std::string& filename, PrimitiveDrawer* input_primitiveDrawer, meshDrawer* input_meshDrawer)
{
    spdlog::info("load homogenized mesh structure config");
    loadHomogenizedMeshStructureConfig(filename);//use config to generate vertices and edges for rod structure
    std::cout << "done" << std::endl;

    auto fileType = m_meshObj_str.find_last_of(".") + 1;

    m_shellStructure = nullptr;

    if (m_isPeriodic)
        m_shellStructure = new SOPeriodicBendingModeShellStructure();
    else if (m_isReflection)
        m_shellStructure = new SOReflectionBendingModeShellStructure();
    else
    m_shellStructure = new SOBendingModeShellStructure();


    if (m_shellStructure != nullptr)
    {
        m_shellStructure->setMeshDrawer(input_meshDrawer);
        m_shellStructure->setPrimitiveDrawer(input_primitiveDrawer);
    }
    //then, we set structure after adding boundary points
    setHomogenizedMeshStructure();


    spdlog::info("rod structure initialization has been done");

}

void foldingShellMeshStructure::setHomogenizedMeshStructure()
{
    SOMaterial* pMat = new SOMaterial;
    pMat->k1 = m_youngsModulus;// young's modulus
    pMat->k2 = m_poissonsRatio;// poisson's ratio
    pMat->kB = pMat->k1 * m_thickness * m_thickness * m_thickness / (24.0 * (1.0 - pMat->k2 * pMat->k2));
    pMat->rho = m_density;
    pMat->matModel = MM_StVKMod;


    string ouputFolder = "./output";
    std::size_t const folder_startIndex = m_meshObj_str.find_last_of("/") + 1;
    std::size_t const folder_endIndex = m_meshObj_str.find_last_of(".") - 1;
    string outputSubFolder = m_meshObj_str.substr(folder_startIndex, (folder_endIndex - folder_startIndex + 1));
    ouputFolder += "/" + outputSubFolder;
    struct stat st_2 = { 0 };
    if (stat(ouputFolder.c_str(), &st_2) == -1) {
#ifdef _WIN32
            mkdir(ouputFolder.c_str());
#elif __APPLE__ || __linux__ 
            mkdir(ouputFolder.c_str(),0777);
#endif
    }


    bool use_shapeOperator = false;
    if (m_shellStructure)
    {
        //load structure or mechanisms
        TriMesh* pMesh = new TriMesh;

        /* if (!pMesh->loadFromFile(m_meshObj_str.c_str()))
         {
             printf("Couldn't open file %s for loading. Terminating..", m_meshObj_str.c_str());
             exit(1);
         }*/

        std::string error_message;
        if (!pMesh->loadFromFile_tinyobj(m_meshObj_str, error_message))
        {
            printf("Couldn't open file %s for loading. Terminating.. with error message %s", m_meshObj_str.c_str(), error_message.c_str());
            exit(1);
        }

        std::cout << "Initialize triangle mesh structure in simulation" << std::endl;

        SOPeriodicBendingModeShellStructure* m_periodicShellStructure = dynamic_cast<SOPeriodicBendingModeShellStructure*>(m_shellStructure);
        if (m_periodicShellStructure != nullptr)
            m_periodicShellStructure->init(pMesh, pMat,
                m_thickness,
                m_fixVertexOnLinePairs,
                fixVertexDofs, vertexIndex_force,
                m_periodicLeftToRight, m_periodicTopToBottom,
                ouputFolder);
        else
        {
            SOReflectionBendingModeShellStructure* m_reflectionShellStructure = dynamic_cast<SOReflectionBendingModeShellStructure*>(m_shellStructure);
            if (m_reflectionShellStructure != nullptr)
                m_reflectionShellStructure->init(pMesh, pMat,
                    m_thickness,
                    m_fixVertexOnLinePairs,
                    fixVertexDofs, vertexIndex_force,
                    m_fourReflectionPlane,
                    ouputFolder);
            else
            {
                TriMesh* target_pMesh = NULL;
                if (!m_targetMeshObj_str.empty())
                {
                    target_pMesh = new TriMesh;

                    if (!target_pMesh->loadFromFile_tinyobj(m_targetMeshObj_str, error_message))
                    {
                        printf("Couldn't open file %s for loading. Terminating.. with error message %s", m_targetMeshObj_str.c_str(), error_message.c_str());
                        exit(1);
                    }
                }
                m_shellStructure->init(pMesh, pMat, target_pMesh,
                    m_thickness,
                    m_fixVertexOnLinePairs,
                    fixVertexDofs, vertexIndex_force,
                    //m_periodicLeftToRight, m_periodicTopToBottom,
                    use_shapeOperator,
                    ouputFolder);
            }
        }
    }
    
}

bool foldingShellMeshStructure::init(const std::string& structure_szFileName, PrimitiveDrawer* p_primitiveDrawer, meshDrawer* p_meshDrawer, bool m_outputRunningStates)
{
	m_meshDrawer = p_meshDrawer;
    m_primitiveDrawer = p_primitiveDrawer;

    std::cout << "Init config from file " << structure_szFileName << std::endl;

    homogenizedMeshConfigInit(structure_szFileName.c_str(), p_primitiveDrawer, p_meshDrawer);
    
    m_interlockStructure_lastState = m_shellStructure->getParamPosition();
    
    spdlog::info("Interlocking structure has been initilized");
	return true;
}

void foldingShellMeshStructure::render()
{
    if(m_shellStructure)
        m_shellStructure->render();
}

void foldingShellMeshStructure::renderImguiViewer()
{
    if (m_shellStructure)
        m_shellStructure->renderImGUIViewerWindow();
    //materialTestImguiViewer();
}

bool foldingShellMeshStructure::step()
{
    bool optimized_flag = false;


    if (m_shellStructure)
    {
        m_interlockStructure_lastState = m_shellStructure->getCurrentPositions();
        //run simulation
        m_shellStructure->step();
        simulationId++;
        optimized_flag = true;
    }
   
    return optimized_flag;
}

void foldingShellMeshStructure::reset()
{
    //reset configuration to the initial one
    if (m_shellStructure)
        m_shellStructure->reset();
}

void foldingShellMeshStructure::mouseLeftClick(const Eigen::Vector3d& src, const Eigen::Vector3d& ray)
{
    if (m_shellStructure)
        m_shellStructure->mouseLeftClick(src, ray);
}

void foldingShellMeshStructure::loadConfig(nlohmann::json& j)
{
    auto config_json = j.find("IPCConfig");
    bool config_jsonis_array = config_json->is_array();
    if (!config_jsonis_array) throw std::invalid_argument("parse error");

    if (config_json->size() != 1) throw std::invalid_argument("parse error");
    
    double self_friction = 0.1, dHat = 1e-3, epsv = 1e-3;
    double timeStep = 0.005;
    bool turnOffGravity = false;
    for (auto& el : (*config_json))
    {
        if (!el.is_object()) throw std::invalid_argument("parse error");
        //self friction parse
        auto el_selfFric = *el.find("selfFric");
        self_friction = el_selfFric.get<double>();
        //dHat
        auto el_dHat = *el.find("dHat");
        dHat = el_dHat.get<double>();
        //epsv
        auto el_epsv = *el.find("epsv");
        epsv = el_epsv.get<double>();
        //turnOffGravity
        auto el_turnOffGravity = *el.find("turnOffGravity");
        turnOffGravity = el_turnOffGravity.get<bool>();
        //timeStep
        auto el_timeStep = *el.find("timeStep");
        timeStep = el_timeStep.get<double>();
    }

    m_self_friction = self_friction;
    m_dHat = dHat;
    m_epsv = epsv;
    m_turnOffGravity = turnOffGravity;
    m_timeStep = timeStep;
}

void foldingShellMeshStructure::loadHomogenizedMeshStructureConfig(const std::string& filename)
{
	std::string m_filename = filename;

	std::string folder = get_folder(filename);
	std::string contents = get_file_contents(filename.c_str());
    //std::cout << contents << std::endl;
	std::istringstream jsonInputStream(contents);
	nlohmann::json j;
	jsonInputStream >> j;

    loadHomogenizedMesh(j);

    //boundary conditions
    vertexIndex_force.clear();
    //std::vector<std::pair<int, Eigen::Vector3d>> vertexIndex_force;
    auto externalVertexForce_json = j.find("externalVertexForce");
    bool externalVertexForce_jsonis_array = externalVertexForce_json->is_array();
    if (!externalVertexForce_jsonis_array) throw std::invalid_argument("parse error");

    for (auto& el : (*externalVertexForce_json))
    {
        if (!el.is_object()) throw std::invalid_argument("parse error");
        //vertexIndex parse
        auto el_vertexIndex = *el.find("vertexIndex");
        int vertexIndex = el_vertexIndex.get<int>();
        //forceVector
        V3D forceVector(0.0, 0.0, 0.0);
        auto el_forceVector = *el.find("forceVector");
        if (!el_forceVector.empty())
        {
            if (!el_forceVector.is_array() || el_forceVector.size() != 3)
            {
            }
            else
            {

                forceVector = V3D(el_forceVector[0].get<double>(), el_forceVector[1].get<double>(), el_forceVector[2].get<double>());
            }
        }
        vertexIndex_force.push_back(std::make_pair(vertexIndex, forceVector));
    }


    fixVertexDofs.clear();
    auto fixVertexDofs_json = j.find("fixVertexDofs");
    if (fixVertexDofs_json != j.end())
    {
        for (auto& el : *fixVertexDofs_json)
        {
            if (!el.is_array()) throw std::invalid_argument("parse error");
            if (el.size() != 2) throw std::invalid_argument("parse error");
            fixVertexDofs.emplace_back(el[0].get<int>(), el[1].get<int>());
        }
    }

    m_fixVertexOnLinePairs.clear();
    auto fixOnline_json = j.find("fixOnLineElements");
    if (fixOnline_json != j.end())
    {
        for (auto& el : *fixOnline_json)
        {
            if (!el.is_array()) throw std::invalid_argument("parse error");
            if (el.size() != 2) throw std::invalid_argument("parse error");

            m_fixVertexOnLinePairs.emplace_back(el[0].get<int>(), el[1].get<int>());
        }
    }


	//material properties
	m_density = j.find("density")->get<double>();
	//m_crossSectionRadius = j.find("crossSectionRadius")->get<double>();
	m_youngsModulus = j.find("youngsModulus")->get<double>();
	m_poissonsRatio = j.find("poissonsRatio")->get<double>();
    m_thickness = j.find("thickness")->get<double>();
    //load ipc config
    loadConfig(j);
}

void foldingShellMeshStructure::loadHomogenizedMesh(nlohmann::json& j)
{
    auto structureParameter_json = j.find("structureParameter");
    //define variables
    for (auto& el : (*structureParameter_json))
    {
        auto meshObj_jason = *el.find("meshObj");
        m_meshObj_str = meshObj_jason.get<string>();


        auto targetMeshObj_jason = *el.find("targetMeshObj");
        m_targetMeshObj_str = targetMeshObj_jason.get<string>();

        auto m_collisionMeshObj_jason = *el.find("collisionMeshObj");
        m_collisionMeshObj_str = m_collisionMeshObj_jason.get<string>();

        
        //auto m_strainSpaceLimitingModel_jason = *el.find("strainSpaceLimitingModel");
        //m_strainSpaceLimitingModel_str = m_strainSpaceLimitingModel_jason.get<string>();

        //auto m_bendingStrainSpaceLimitingModel_jason = *el.find("bendingStrainSpaceLimitingModel");
        //m_bendingStainSpaceLimitingModel_str = m_bendingStrainSpaceLimitingModel_jason.get<string>();

        auto isPeriodic_jason = *el.find("isPeriodic");
        m_isPeriodic = isPeriodic_jason.get<bool>();

        if (m_isPeriodic)
        {
            auto periodicLeftToRight_jason = *el.find("periodicLeftToRight");
            m_periodicLeftToRight = Eigen::Vector3d(0, 0, 0);
            if (!periodicLeftToRight_jason.empty())
            {
                if (!periodicLeftToRight_jason.is_array() || periodicLeftToRight_jason.size() != 3)
                {
                }
                else
                {

                    m_periodicLeftToRight = V3D(periodicLeftToRight_jason[0].get<double>(), periodicLeftToRight_jason[1].get<double>(), periodicLeftToRight_jason[2].get<double>());
                }
            }


            auto periodicTopToBottom_jason = *el.find("periodicTopToBottom");
            m_periodicTopToBottom = Eigen::Vector3d(0, 0, 0);
            if (!periodicTopToBottom_jason.empty())
            {
                if (!periodicTopToBottom_jason.is_array() || periodicTopToBottom_jason.size() != 3)
                {
                }
                else
                {

                    m_periodicTopToBottom = V3D(periodicTopToBottom_jason[0].get<double>(), periodicTopToBottom_jason[1].get<double>(), periodicTopToBottom_jason[2].get<double>());
                }
            }
        }


        auto isReflection_jason = *el.find("isReflection");
        m_isReflection = isReflection_jason.get<bool>();
        if (m_isReflection)
        {
            auto fourReflectionPlane_jason = *el.find("fourReflectionPlane");
            m_fourReflectionPlane = Eigen::Vector4d(0, 0, 0, 0);
            if (!fourReflectionPlane_jason.empty())
            {
                if (!fourReflectionPlane_jason.is_array() || fourReflectionPlane_jason.size() != 4)
                {
                }
                else
                {
                    m_fourReflectionPlane = Eigen::Vector4d(fourReflectionPlane_jason[0].get<double>(), fourReflectionPlane_jason[1].get<double>(), 
                        fourReflectionPlane_jason[2].get<double>(), fourReflectionPlane_jason[3].get<double>());
                }
            }
        }

    }
}
