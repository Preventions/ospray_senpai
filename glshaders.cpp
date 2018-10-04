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
#include <stdexcept>
#include <vector>
#include <GL/gl3w.h>
#include "glshaders.h"

GLShader::GLShader(const std::string &vsrc, const std::string &fsrc)
	: shader(load_shader_program(vsrc, fsrc))
{}
GLuint GLShader::uniform_location(const std::string &unif) {
	auto fnd = uniforms.find(unif);
	if (fnd != uniforms.end()) {
		return fnd->second;
	} else {
		uniforms[unif] = glGetUniformLocation(shader, unif.c_str());
		return uniforms[unif];
	}
}
GLShader::~GLShader() {
	glDeleteProgram(shader);
}

GLuint compile_shader(const std::string &src, GLenum type) {
	GLuint shader = glCreateShader(type);
	const char *csrc = src.c_str();
	glShaderSource(shader, 1, &csrc, 0);
	glCompileShader(shader);

	GLint status;
	glGetShaderiv(shader, GL_COMPILE_STATUS, &status);
	if (status == GL_FALSE) {
		std::cout << "Shader compilation error:\n";
		GLint len;
		glGetShaderiv(shader, GL_INFO_LOG_LENGTH, &len);
		std::vector<char> log(len, 0);
		glGetShaderInfoLog(shader, log.size(), 0, log.data());
		std::cout << log.data() << "\n";
		throw std::runtime_error("Shader compilation failed");
	}
	return shader;
}
GLuint load_shader_program(const std::string &vshader_src, const std::string &fshader_src) {
	GLuint vs = compile_shader(vshader_src, GL_VERTEX_SHADER);
	GLuint fs = compile_shader(fshader_src, GL_FRAGMENT_SHADER);
	GLuint prog = glCreateProgram();
	glAttachShader(prog, vs);
	glAttachShader(prog, fs);
	glLinkProgram(prog);

	GLint status;
	glGetProgramiv(prog, GL_LINK_STATUS, &status);
	if (status == GL_FALSE) {
		std::cout << "Shader link error:\n";
		GLint len;
		glGetProgramiv(prog, GL_INFO_LOG_LENGTH, &len);
		std::vector<char> log(len, 0);
		glGetProgramInfoLog(prog, log.size(), 0, log.data());
		std::cout << log.data() << "\n";
		throw std::runtime_error("Shader link failed");
	}
	glDetachShader(prog, vs);
	glDetachShader(prog, fs);
	glDeleteShader(vs);
	glDeleteShader(fs);
	return prog;
}

