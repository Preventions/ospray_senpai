cmake_minimum_required(VERSION 3.1)
project(viewer)

option(LAMMPS_IS_RENDER_WORKER "LAMMPS In Situ Rendering Worker" ON)
option(LAMMPS_IS_VIEWER "LAMMPS In Situ Remote Viewer" ON)

# Use our FindGLM module
set(CMAKE_MODULE_PATH "${viewer_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if (NOT WIN32)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -pedantic -std=c++11")
	set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_RELEASE} -O3")
endif()
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_DEBUG}")

if (WIN32)
	add_definitions(-DNOMINMAX)
endif()

find_package(ospray REQUIRED)
find_package(Threads REQUIRED)
find_package(TurboJpeg REQUIRED)
find_package(VTK REQUIRED)
include_directories(${OSPRAY_INCLUDE_DIRS} ${VTK_INCLUDE_DIRS})

add_library(lammps_is_util SHARED
	util.cpp
	image_util.cpp
	client_server.cpp
	socket_fabric.cpp
	colormap.cpp)
target_link_libraries(lammps_is_util PUBLIC
	${OSPRAY_LIBRARIES}
	${CMAKE_THREAD_LIBS_INIT}
	${VTK_LIBRARIES}
	TurboJpeg)

if (WIN32)
	target_link_libraries(lammps_is_util PUBLIC Ws2_32)
endif()


if (LAMMPS_IS_RENDER_WORKER)
	find_package(MPI REQUIRED)
	find_package(libIS REQUIRED)
	include_directories(${MPI_CXX_INCLUDE_DIRS} ${MPI_CXX_INCLUDE_PATH})

	add_executable(lammps_is_render_worker
		lammps_is_render_worker.cpp
		colormap.cpp
		query_task.cpp)
	target_link_libraries(lammps_is_render_worker PUBLIC
		lammps_is_util
		is_render
		${OSPRAY_LIBRARIES}
		${OSPRAY_MPI_LIBRARIES}
		${MPI_CXX_LIBRARIES})
endif()

if (LAMMPS_IS_VIEWER)
	add_definitions(-DGLM_FORCE_RADIANS)
	find_package(GLM REQUIRED)
	find_package(OpenGL REQUIRED)
	find_package(GLFW REQUIRED)

	add_subdirectory(gl3w)
	add_subdirectory(imgui)

	add_executable(lammps_is_viewer
		lammps_is_viewer.cpp
		arcball.cpp
		gldebug.cpp
		glshaders.cpp
		imgui_impl_glfw_gl3.cpp)
	target_include_directories(lammps_is_viewer PRIVATE
		$<BUILD_INTERFACE:${GLFW_INCLUDE_DIR}>
		$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/gl3w>
		$<BUILD_INTERFACE:${GLM_INCLUDE_DIRS}>
		$<BUILD_INTERFACE:${OPENGL_INCLUDE_DIRS}>)
	target_link_libraries(lammps_is_viewer PUBLIC
		lammps_is_util
		${GLFW_LIBRARY}
		${OPENGL_LIBRARIES}
		gl3w
		imgui)
endif()

