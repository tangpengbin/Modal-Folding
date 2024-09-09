cmake_minimum_required(VERSION 3.17.0)
include(FetchContent)

#add our SimLib
message(STATUS "SimLib")
FetchContent_Declare(
	SimLib
	GIT_REPOSITORY https://github.com/tangpengbin/SimLib.git
	GIT_TAG main
)
FetchContent_MakeAvailable(SimLib)


project(SimShellProject)
SET(CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake") 

set(PROJECT_SRC_DIR "${CMAKE_CURRENT_LIST_DIR}/src")

file(GLOB_RECURSE source_files
	${PROJECT_SRC_DIR}/*.h
	${PROJECT_SRC_DIR}/*.cpp
	)


#put lib into the project
add_executable(${PROJECT_NAME} ${source_files})


target_link_libraries(${PROJECT_NAME}
		PUBLIC
		SimLib
		)


source_group(TREE ${PROJECT_SRC_DIR} FILES ${source_files})

set_property(TARGET ${PROJECT_NAME} PROPERTY CXX_STANDARD 17)
