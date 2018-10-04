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

#include <cmath>
#include "arcball_camera.h"

ArcBallCamera::ArcBallCamera(const glm::mat4 &look_at, float motion_speed, float rotation_speed,
		const std::array<size_t, 2> &screen)
	: look_at(look_at), translation(glm::mat4{}), rotation(glm::quat{}), camera(look_at),
	inv_camera(glm::inverse(camera)), motion_speed(motion_speed), rotation_speed(rotation_speed),
	inv_screen({1.f / screen[0], 1.f / screen[1]})
{}
const glm::mat4& ArcBallCamera::transform() const {
	return camera;
}
const glm::mat4& ArcBallCamera::inv_transform() const {
	return inv_camera;
}
glm::vec3 ArcBallCamera::eye_pos() const {
	return glm::vec3{inv_camera * glm::vec4{0, 0, 0, 1}};
}
glm::vec3 ArcBallCamera::look_dir() const {
	return glm::vec3{inv_camera * glm::vec4{0, 0, -1, 0}};
}
glm::vec3 ArcBallCamera::up_dir() const {
	return glm::vec3{inv_camera * glm::vec4{0, 1, 0, 0}};
}
void ArcBallCamera::update_screen(const size_t screen_x, const size_t screen_y){
	inv_screen[0] = 1.0f / screen_x;
	inv_screen[1] = 1.0f / screen_y;
}
void ArcBallCamera::reset() {
	translation = glm::mat4{};
	rotation = glm::quat{};
	camera = look_at;
	inv_camera = glm::inverse(camera);
}
void ArcBallCamera::rotate(const glm::vec2 &mouse, const glm::vec2 &mouse_delta, float elapsed){
	// Compute current and previous mouse positions in clip space
	glm::vec2 mouse_cur = glm::vec2{mouse.x * 2.0 * inv_screen[0] - 1.0,
		1.0 - 2.0 * mouse.y * inv_screen[1]};
	glm::vec2 mouse_prev = glm::vec2{(mouse.x - mouse_delta.x) * 2.0 * inv_screen[0] - 1.0,
		1.0 - 2.0 * (mouse.y - mouse_delta.y) * inv_screen[1]};
	// Clamp mouse positions to stay in screen space range
	mouse_cur = glm::clamp(mouse_cur, glm::vec2{-1, -1}, glm::vec2{1, 1});
	mouse_prev = glm::clamp(mouse_prev, glm::vec2{-1, -1}, glm::vec2{1, 1});
	glm::quat mouse_cur_ball = screen_to_arcball(mouse_cur);
	glm::quat mouse_prev_ball = screen_to_arcball(mouse_prev);

	rotation = mouse_cur_ball*mouse_prev_ball*rotation;
	camera = translation * look_at * glm::mat4_cast(rotation);
	inv_camera = glm::inverse(camera);
}
void ArcBallCamera::pan(const glm::vec2 &mouse_delta, float elapsed){
	const glm::vec3 motion(mouse_delta.x * inv_screen[0], -mouse_delta.y * inv_screen[1], 0.f);
	translation = glm::translate(motion * motion_speed * elapsed * 100.f) * translation;
	camera = translation * look_at * glm::mat4_cast(rotation);
	inv_camera = glm::inverse(camera);
}
void ArcBallCamera::zoom(const float amount, const float elapsed){
	glm::vec3 motion{0.f};
	motion.z = amount;
	translation = glm::translate(motion * motion_speed * elapsed) * translation;
	camera = translation * look_at * glm::mat4_cast(rotation);
	inv_camera = glm::inverse(camera);
}
glm::quat screen_to_arcball(const glm::vec2 &p){
	float dist = glm::dot(p, p);
	// If we're on/in the sphere return the point on it
	if (dist <= 1.f){
		return glm::quat(0, p.x, p.y, std::sqrt(1.f - dist));
	}
	// otherwise we project the point onto the sphere
	else {
		const glm::vec2 unit_p = glm::normalize(p);
		return glm::quat(0, unit_p.x, unit_p.y, 0);
	}
}

