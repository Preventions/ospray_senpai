// ======================================================================== //
// Copyright 2018 Intel Corporation                                         //
//                                                                          //
// Licensed under the Apache License, Version 2.0 (the "License");          //
// you may not use this file except in compliance with the License.         //
// You may obtain a copy of the License at                                  //
//                                                                          //
//     http://www.apache.org/licenses/LICENSE-2.0                           //
//                                                                          //
// Unless required by applicable law or agreed to in writing, software      //
// distributed under the License is distributed on an "AS IS" BASIS,        //
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. //
// See the License for the specific language governing permissions and      //
// limitations under the License.                                           //
// ======================================================================== //

#include <iostream>
#include <iomanip>
#include <cassert>
#include <stdexcept>
#include <GL/gl3w.h>
#include "gldebug.h"

void register_debug_callback() {
	glEnable(GL_DEBUG_OUTPUT);
	glEnable(GL_DEBUG_OUTPUT_SYNCHRONOUS);
	glDebugMessageCallback(debug_callback, nullptr);
	glDebugMessageControl(GL_DONT_CARE, GL_DONT_CARE, GL_DONT_CARE, 0, NULL, GL_TRUE);
}
#ifdef _WIN32
void APIENTRY debug_callback(GLenum src, GLenum type, GLuint id, GLenum severity,
	GLsizei len, const GLchar *msg, GLvoid*)
{
	log_debug_msg(src, type, id, severity, len, msg);
}
#else
void debug_callback(GLenum src, GLenum type, GLuint id, GLenum severity,
	GLsizei len, const GLchar *msg, GLvoid*)
{
	log_debug_msg(src, type, id, severity, len, msg);
}
#endif
void log_debug_msg(GLenum src, GLenum type, GLuint, GLenum severity, GLsizei tag, const GLchar *msg){
	// Disable nvidia mapping spam
	if (tag == 300 || tag == 315 || type == GL_DEBUG_TYPE_OTHER){
		return;
	}
  std::cout << "OpenGL Debug -";
	switch (severity){
	case GL_DEBUG_SEVERITY_HIGH:
		std::cout << " High severity";
		break;
	case GL_DEBUG_SEVERITY_MEDIUM:
		std::cout << " Medium severity";
		break;
	case GL_DEBUG_SEVERITY_LOW:
		std::cout << " Low severity";
	}
	switch (src){
	case GL_DEBUG_SOURCE_API:
		std::cout << " API";
		break;
	case GL_DEBUG_SOURCE_WINDOW_SYSTEM:
		std::cout << " Window system";
		break;
	case GL_DEBUG_SOURCE_SHADER_COMPILER:
		std::cout << " Shader compiler";
		break;
	case GL_DEBUG_SOURCE_THIRD_PARTY:
		std::cout << " Third party";
		break;
	case GL_DEBUG_SOURCE_APPLICATION:
		std::cout << " Application";
		break;
	default:
		std::cout << " Other";
	}
	switch (type){
	case GL_DEBUG_TYPE_ERROR:
		std::cout << " Error";
		break;
	case GL_DEBUG_TYPE_DEPRECATED_BEHAVIOR:
		std::cout << " Deprecated behavior";
		break;
	case GL_DEBUG_TYPE_UNDEFINED_BEHAVIOR:
		std::cout << " Undefined behavior";
		break;
	case GL_DEBUG_TYPE_PORTABILITY:
		std::cout << " Portability";
		break;
	case GL_DEBUG_TYPE_PERFORMANCE:
		std::cout << " Performance";
		break;
	default:
		std::cout << " Other";
	}
	std::cout << " Tag: " << tag;
	std::cout << ":\n\t" << msg << "\n";
	// Break for a stack trace of sorts
	if (severity == GL_DEBUG_SEVERITY_HIGH && type == GL_DEBUG_TYPE_ERROR) {
		throw std::runtime_error("Fatal OpenGL Error!");
	}
}

