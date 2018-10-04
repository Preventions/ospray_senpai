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

#include <array>
#include <math.h>
#include <glm/glm.hpp>
#include <glm/ext.hpp>

// A simple arcball camera that moves around the camera's focal point
class ArcBallCamera {
	// We store the unmodified look at matrix along with
	// decomposed translation and rotation components
	glm::mat4 look_at, translation;
	glm::quat rotation;
	// camera is the full camera transform,
	// inv_camera is stored as well to easily compute
	// eye position and world space rotation axes
	glm::mat4 camera, inv_camera;
	// Motion and rotation speeds
	float motion_speed, rotation_speed;
	// Inverse x, y window dimensions
	std::array<float, 2> inv_screen;

public:
	/* Create an arcball camera with some look at matrix
	 * motion speed: units per second speed of panning the camera
	 * rotation speed: radians per second speed of rotation the camera
	 * screen: { WIN_X_SIZE, WIN_Y_SIZE }
	 */
	ArcBallCamera(const glm::mat4 &look_at, float motion_speed, float rotation_speed,
			const std::array<size_t, 2> &screen);
	void update_screen(const size_t screen_x, const size_t screen_y);
	// Get the camera transformation matrix
	const glm::mat4& transform() const;
	// Get the camera inverse transformation matrix
	const glm::mat4& inv_transform() const;
	// Get the eye position of the camera in world space
	glm::vec3 eye_pos() const;
	// Get the look direction of the camera
	glm::vec3 look_dir() const;
	// Get the up direction of the camera
	glm::vec3 up_dir() const;
	/* Handle keyboard events to reset the camera
	 * returns true if the camera has moved
	 */
	void reset();
	// Handle rotation events
	void rotate(const glm::vec2 &mouse, const glm::vec2 &mouse_delta, float elapsed);
	// Handle panning/zooming events
	void pan(const glm::vec2 &mouse_delta, float elapsed);
	// Zoom the camera in/out
	void zoom(const float amount, const float elapsed);
};

// Project the point in [-1, 1] screen space onto the arcball sphere
glm::quat screen_to_arcball(const glm::vec2 &p);

