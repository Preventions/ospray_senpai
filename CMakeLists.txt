cmake_minimum_required(VERSION 3.5)
project(ospray_senpai)

if (NOT WIN32)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -pedantic")
endif()

add_definitions(-DNOMINMAX -DGLM_ENABLE_EXPERIMENTAL -DSDL_MAIN_HANDLED)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_LIST_DIR}/cmake")

find_package(Threads REQUIRED)
find_package(glm REQUIRED)
find_package(SDL2 REQUIRED)
find_package(OpenGL REQUIRED)
find_package(ospray REQUIRED)
find_package(libIS REQUIRED)

add_subdirectory(imgui)
add_subdirectory(util)

add_executable(ospray_senpai main.cpp)

set_target_properties(ospray_senpai PROPERTIES
	CXX_STANDARD 14
	CXX_STANDARD_REQUIRED ON)

target_link_libraries(ospray_senpai PUBLIC util ospray::ospray)

