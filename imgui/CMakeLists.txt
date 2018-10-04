project(imgui)

add_library(imgui STATIC imgui.cpp imgui_demo.cpp imgui_draw.cpp)
set_property(TARGET imgui PROPERTY CXX_STANDARD_REQUIRED 11)
target_include_directories(imgui PUBLIC
	$<BUILD_INTERFACE:${imgui_BINARY_DIR}>
	$<BUILD_INTERFACE:${imgui_SOURCE_DIR}>
	$<INSTALL_INTERFACE:include>
)

