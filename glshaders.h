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

#pragma once

#include <string>
#include <unordered_map>
#include <GL/gl3w.h>

#define XSTRINGIFY(X) STRINGIFY(X)
#define STRINGIFY(X) #X

const static std::string aoParamsBlockSrc = 
"#line " XSTRINGIFY(__LINE__) R"(
layout(std140) uniform AOParams {
	int n_samples;
	int n_turns;
	float ball_radius;
	float sigma;
   	float kappa;
	int filter_scale;
	float edge_sharpness;
	int num_blur_passes;
} ao_params;
)";

const static std::string vsrc =R"(
#version 330 core
#line )" XSTRINGIFY(__LINE__) R"(

const vec4 pos[4] = vec4[4](
  vec4(-1, 1, 0.5, 1),
	vec4(-1, -1, 0.5, 1),
	vec4(1, 1, 0.5, 1),
	vec4(1, -1, 0.5, 1)
);

void main(void){
	gl_Position = pos[gl_VertexID];
}
)";


const static std::string fAOComputeSrc =R"(
#version 430 core
)" + aoParamsBlockSrc + R"(
#line )" XSTRINGIFY(__LINE__) R"(

uniform sampler2D camera_depths;

uniform vec2 viewport_dim;
uniform vec3 cam_du;
uniform vec3 cam_dv;
uniform vec3 cam_pos;
uniform vec3 cam_dir;

out vec2 ao_out;

// Compute the position of the pixel, based on its depth and the camera info
vec3 compute_pixel_position(float depth, ivec2 px) {
	vec2 p = vec2(px) / viewport_dim - vec2(0.5);
	vec3 px_dir = normalize(cam_dir + cam_du * p.x + cam_dv * p.y);
	return cam_pos + depth * px_dir;
}

//#define UINT_MAX 4294967295
#define INT_MAX 2147483647
// Pseudo-random number gen from
// http://www.reedbeta.com/blog/quick-and-easy-gpu-random-numbers-in-d3d11/
int wang_hash(int seed) {
	seed = (seed ^ 61) ^ (seed >> 16);
	seed *= 9;
	seed = seed ^ (seed >> 4);
	seed *= 0x27d4eb2d;
	seed = seed ^ (seed >> 15);
	return seed;
}

void main(void){
	const ivec2 px = ivec2(gl_FragCoord.xy);
	const float depth = texture(camera_depths, gl_FragCoord.xy / viewport_dim).x;
	if (depth > 1e10) {
		ao_out = vec2(1, depth);
		return;
	}
	vec3 pos = compute_pixel_position(depth, px);
	vec3 normal = normalize(cross(dFdx(pos), dFdy(pos)));

	// The Alchemy AO hash for random per-pixel offset
	// I'm not sure what the range we get out of this hash is, it doesn't seem to be [0, UINT_MAX], so..
	float phi = (wang_hash(int(gl_FragCoord.x + viewport_dim.x * gl_FragCoord.y)) % INT_MAX) / float(INT_MAX);
	const float TAU = 6.2831853071795864;
	float ball_radius_sqr = pow(ao_params.ball_radius, 2);
	// What's the radius of a 1m object at z = -1m to compute screen_radius properly?
	// Comments in their code mention we can compute it from the projection mat, or hardcode in like 500
	// and make the ball radius resolution dependent (as I've done currently)
	float screen_radius = abs(ao_params.ball_radius * 500.0 / depth);

	int max_mip = textureQueryLevels(camera_depths) - 1;
	float ao_value = 0;
	for (int i = 0; i < ao_params.n_samples; ++i){
		float alpha = 1.f / ao_params.n_samples * (i + 0.5);
		float h = screen_radius * alpha;
		float theta = TAU * alpha * ao_params.n_turns + phi;
		vec2 u = vec2(cos(theta), sin(theta));
		int m = clamp(findMSB(int(h)) - 4, 0, max_mip);
		ivec2 sample_pos = ivec2(h * u) + px;
		ivec2 mip_pos = clamp(sample_pos >> m, ivec2(0), textureSize(camera_depths, m) - ivec2(1));

		// TODO: I think this px pos we compute here is correct
		float q_depth = texelFetch(camera_depths, mip_pos, m).x;
		vec3 q = compute_pixel_position(q_depth, sample_pos);
		vec3 v = q - pos;
		float vv = dot(v, v);
		float vn = dot(v, normal);

		float f = max(ball_radius_sqr - vv, 0.0);
		ao_value += f * f * f * max((vn - 0.001) / (0.01 + vv), 0.0);
		//ao_value += float(vv < ball_radius_sqr) * max((vn - 0.001) / (0.01 + vv), 0.0) * ball_radius_sqr;
		//ao_value += max(0, dot(v, normal + depth * 0.0005)) / (dot(v, v) + 0.01);
	}
	// The original method in paper, from Alchemy AO
	ao_value = max(0, 1.f - 2.f * ao_params.sigma / ao_params.n_samples * ao_value);
	ao_value = pow(ao_value, ao_params.kappa);

	// Do a little bit of filtering now, respecting depth edges
	if (abs(dFdx(depth)) < 0.02) {
		ao_value -= dFdx(ao_value) * ((px.x & 1) - 0.5);
	}
	if (abs(dFdy(depth)) < 0.02) {
		ao_value -= dFdy(ao_value) * ((px.y & 1) - 0.5);
	}
	ao_value = clamp(ao_value, 0, 1);
	ao_out = vec2(ao_value, depth);
}
)";

const static std::string fBlurSrc =R"(
#version 330 core
#line )" XSTRINGIFY(__LINE__) R"(
)" + aoParamsBlockSrc + R"(

#define RADIUS 4

// Gaussian filter values from the author's blurring shader
const float gaussian[RADIUS + 1] = float[](0.153170, 0.144893, 0.122649, 0.092902, 0.062970);

uniform sampler2D ao_in;
uniform ivec2 axis;

out vec2 result;

void main(void){
	ivec2 px = ivec2(gl_FragCoord.xy);
	vec2 val = texelFetch(ao_in, px, 0).xy;
	float z_pos = val.y;

	// Compute weighting for the term at the center of the kernel
	float base = gaussian[0];
	float weight = base;
	float sum = weight * val.x;

	for (int i = -RADIUS; i <= RADIUS; ++i){
		// We handle the center pixel above so skip that case
		if (i != 0){
			// Filter scale effects how many pixels the kernel actually covers
			ivec2 p = px + axis * i * ao_params.filter_scale;
			vec2 fval = texelFetch(ao_in, p, 0).xy;
			float z = fval.y;
			float w = 0.3 + gaussian[abs(i)];
			// Decrease weight as depth difference increases. This prevents us from
			// blurring across depth discontinuities
			w *= max(0.f, 1.f - ao_params.edge_sharpness * abs(z_pos - z));
			sum += fval.x * w;
			weight += w;
		}
	}
	result = vec2(sum / (weight + 0.00001), val.y);
}
)";

const static std::string fFinalCompositeSrc = R"(
#version 330 core
#line )" XSTRINGIFY(__LINE__) R"(

uniform sampler2D camera_colors;
uniform sampler2D ao_vals;

out vec4 color;

void main(void) {
    ivec2 px = ivec2(gl_FragCoord.xy);
    color = texelFetch(camera_colors, px, 0) * texelFetch(ao_vals, px, 0).x;
	color.a = 1;
}
)";

struct GLShader {
	GLuint shader;
	std::unordered_map<std::string, GLuint> uniforms;

	GLShader(const std::string &vsrc, const std::string &fsrc);
	GLuint uniform_location(const std::string &unif);
	~GLShader();
};

GLuint compile_shader(const std::string &src, GLenum type);
GLuint load_shader_program(const std::string &vshader_src,
		const std::string &fshader_src);

