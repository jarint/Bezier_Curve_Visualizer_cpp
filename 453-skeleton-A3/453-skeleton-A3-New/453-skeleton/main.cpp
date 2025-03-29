#include <glad/glad.h>

#include "Geometry.h"
#include "GLDebug.h"
#include "Log.h"
#include "ShaderProgram.h"
#include "Shader.h"
#include "Texture.h"
#include "Window.h"
#include "CurveControl.h"

const int EDITOR_SCREENS = 4;

class CurveEditor : public CallbackInterface {

public:
    bool left_click_active, right_click_active, left_click_acknowledged, 
			delete_requested, reset_requested;

    // Curve mode, scene state, wireframe
    bool bezier_mode;   // true = bezier, false = b-spline
    int scene_mode;
	bool wireframe; // true = wireframe, false = solid

	// moves
    bool move_forward, move_left, move_backward, move_right;

    glm::vec3 last_click_position;


	CurveEditor(int window_width, int window_height)
    	: left_click_active(false), right_click_active(false), left_click_acknowledged(false),
      	  delete_requested(false), reset_requested(false), bezier_mode(true), scene_mode(0), wireframe(true),
      	  move_forward(false), move_left(false), move_backward(false), move_right(false),
      	  last_click_position(0.0f), raw_mouse_delta(0), raw_mouse_position(0),
      	  screen_dimensions(window_width, window_height) 
	{}

    CurveEditor() : CurveEditor(800, 800) {}

	// keyboard input handling
	virtual void keyCallback(int key, int scancode, int action, int mods) override {
	    if (action == GLFW_PRESS) {
	        if (key == GLFW_KEY_BACKSPACE) delete_requested = true;
	        else if (key == GLFW_KEY_R) reset_requested = true;
	        else if (key == GLFW_KEY_ENTER) bezier_mode = !bezier_mode;
	        else if (key == GLFW_KEY_SPACE) scene_mode = (scene_mode + 1) % EDITOR_SCREENS;
			else if (key == GLFW_KEY_T) wireframe = !wireframe;
	        else if (key == GLFW_KEY_W) move_forward = true;
	        else if (key == GLFW_KEY_A) move_left = true;
	        else if (key == GLFW_KEY_S) move_backward = true;
	        else if (key == GLFW_KEY_D) move_right = true;
	    }
	    if (action == GLFW_RELEASE) {
	        if (key == GLFW_KEY_W) move_forward = false;
	        else if (key == GLFW_KEY_A) move_left = false;
	        else if (key == GLFW_KEY_S) move_backward = false;
	        else if (key == GLFW_KEY_D) move_right = false;
	    }
	}
    // Mouse button handling
	virtual void mouseButtonCallback(int button, int action, int mods) override {
	    if (action == GLFW_PRESS) {
	        if (button == GLFW_MOUSE_BUTTON_LEFT) { 
	            left_click_active = true; 
	            last_click_position = get_world_position_from_mouse(); 
	        } else if (button == GLFW_MOUSE_BUTTON_RIGHT) { 
	            right_click_active = true; 
	            last_click_position = get_world_position_from_mouse(); 
	        }
	    }
	    if (action == GLFW_RELEASE) {
	        if (button == GLFW_MOUSE_BUTTON_LEFT) { 
	            left_click_active = false; 
	            left_click_acknowledged = false; 
	        } else if (button == GLFW_MOUSE_BUTTON_RIGHT) { 
	            right_click_active = false; 
	        }
	    }
	}
    virtual void cursorPosCallback(double xpos, double ypos) override {
        raw_mouse_delta = glm::ivec2(xpos - raw_mouse_position.x, ypos - raw_mouse_position.y);
        raw_mouse_position = glm::vec2(xpos, ypos);
    }
    virtual void windowSizeCallback(int width, int height) override {
        CallbackInterface::windowSizeCallback(width, height);
        screen_dimensions = glm::ivec2(width, height);
    }
    glm::vec3 get_world_position_from_mouse() const {
        return glm::vec3(
            ((raw_mouse_position.x * 2.0f) / screen_dimensions.x) - 1.0f,
            -(((raw_mouse_position.y * 2.0f) / screen_dimensions.y) - 1.0f),
            0.0f
        );
    }
    glm::ivec2 consume_mouse_delta() {
        glm::ivec2 delta = raw_mouse_delta;
        raw_mouse_delta = glm::ivec2(0);
        return delta;
    }

private:
    glm::vec2 raw_mouse_delta;
    glm::vec2 raw_mouse_position;
    glm::ivec2 screen_dimensions;
};

void update_gpu_geometry(GPU_Geometry &gpu_geometry, const CPU_Geometry &cpu_geometry) {
    gpu_geometry.bind();
    gpu_geometry.setVerts(cpu_geometry.verts);
    gpu_geometry.setCols(cpu_geometry.cols);
}

void reset_points(CPU_Geometry& geometry, int& selected_point) {
    selected_point = -1;
    geometry.verts.clear();
    geometry.cols.clear();
	// starting points
    geometry.verts = {
        { -0.5f,  0.5f, 0.0f },
        { 0.5f, 0.5f, 0.0f },
        {  0.5f, -0.5f, 0.0f },
    };
}

// deCasteljau's algorithm as described in tutorial with Jeffrey
std::vector<glm::vec3> deCasteljau(std::vector<glm::vec3>& control_points) {
    std::vector<glm::vec3> curve_points;
    int point_count = control_points.size();
    for (float u = 0.0f; u <= 1.0f; u += 0.01f) {
        std::vector<glm::vec3> temp_points = control_points;
        for (int i = 1; i < point_count; i++) {
            for (int j = 0; j < point_count - i; j++) {
                temp_points[j] = (1 - u) * temp_points[j] + u * temp_points[j + 1];
            }
        }
        curve_points.push_back(temp_points[0]);
    }
    return curve_points;
}

// a single iteration of chaikin subdivision
std::vector<glm::vec3> chaikin_iteration(std::vector<glm::vec3>& control_points) {
    std::vector<glm::vec3> subdivided_points;
    int num_points = (int)control_points.size();

    subdivided_points.push_back(control_points[0]);
    subdivided_points.push_back(0.5f * control_points[0] + 0.5f * control_points[1]);

    for (int i = 1; i < num_points - 2; i++) {
        subdivided_points.push_back(0.75f * control_points[i] + 0.25f * control_points[i + 1]);
        subdivided_points.push_back(0.25f * control_points[i] + 0.75f * control_points[i + 1]);
    }

    subdivided_points.push_back(0.5f * control_points[num_points - 2] + 0.5f * control_points[num_points - 1]);
    subdivided_points.push_back(control_points[num_points - 1]);

    return subdivided_points;
}

// repeatedly applies the subdivision for the given number of levels
std::vector<glm::vec3> chaikin_subdivision(std::vector<glm::vec3>& input_points, const int& subdivision_levels) {
    std::vector<glm::vec3> subdivided_points = input_points;

    for (int i = 0; i < subdivision_levels; i++) {
        subdivided_points = chaikin_iteration(subdivided_points);
    }
    return subdivided_points;
}

// flattens a grid of control points into a single vector
std::vector<glm::vec3> flatten_surface_grid(std::vector<std::vector<glm::vec3>>& grid_points, std::vector<unsigned int>& index_buffer) {
    std::vector<glm::vec3> flattened_points;
    int num_rows = (int)grid_points.size();
    int num_cols = (int)grid_points[0].size();

    for (int col = 0; col < num_cols; col++) {
        for (int row = 0; row < num_rows; row++) {
            flattened_points.push_back(grid_points[row][col]);
        }
    }

    for (int base_index = 0; base_index < flattened_points.size() - num_cols; base_index += num_cols) {
        for (int col_idx = 0; col_idx < num_cols - 1; col_idx++) {
            unsigned int top_left = base_index + col_idx;
            unsigned int bottom_left = top_left + num_cols;
            unsigned int bottom_right = bottom_left + 1;
            unsigned int top_right = top_left + 1;

            index_buffer.push_back(top_left);
            index_buffer.push_back(bottom_left);
            index_buffer.push_back(bottom_right);

            index_buffer.push_back(top_left);
            index_buffer.push_back(bottom_right);
            index_buffer.push_back(top_right);
        }
    }

    return flattened_points;
}

// generates a surface by subdividing control points along rows and columns
std::vector<glm::vec3> generate_chaikin_surface(std::vector<glm::vec3>& control_points, std::vector<unsigned int>& index_buffer) {
    int row_count = (int)sqrt(control_points.size());
    std::vector<std::vector<glm::vec3>> subdivided_rows, subdivided_columns;
    std::vector<glm::vec3> temp_points;

    // Subdivide each column
    for (int col = 0; col < row_count; col++) {
        temp_points.clear();
        for (int row = col; row < control_points.size(); row += row_count) {
            temp_points.push_back(control_points[row]);
        }
        subdivided_rows.push_back(chaikin_subdivision(temp_points, 4));
    }

    // Subdivide each row
    for (int i = 0; i < subdivided_rows[0].size(); i++) {
        temp_points.clear();
        for (int j = 0; j < row_count; j++) {
            temp_points.push_back(subdivided_rows[j][i]);
        }
        subdivided_columns.push_back(chaikin_subdivision(temp_points, 4));
    }

    return flatten_surface_grid(subdivided_columns, index_buffer);
}


// revolve points around the Y axis
std::vector<glm::vec3> surface_revolution(std::vector<glm::vec3>& profile_points, std::vector<unsigned int>& indices) {
    std::vector<glm::vec3> revolved_vertices;
    double angle_step = M_PI / 40.0;
    int profile_size = (int)profile_points.size();

    // Generate revolved points by rotating around the Y-axis
    for (double angle = 0; angle < 2 * M_PI; angle += angle_step) {
        for (const auto& point : profile_points) {
            revolved_vertices.push_back(glm::vec3(
                point.x * cos(angle),
                point.y,
                point.x * sin(angle)
            ));
        }
    }

    int revolved_size = (int)revolved_vertices.size();

    // Creates triangle indices for connecting the revolved rings
    for (int ring_start = 0; ring_start < revolved_size; ring_start += profile_size) {
        for (int i = 0; i < profile_size - 1; i++) {
            indices.push_back((i + ring_start) % revolved_size);
            indices.push_back((i + ring_start + profile_size) % revolved_size);
            indices.push_back((i + 1 + ring_start + profile_size) % revolved_size);

            indices.push_back((i + ring_start) % revolved_size);
            indices.push_back((i + 1 + ring_start) % revolved_size);
            indices.push_back((i + 1 + ring_start + profile_size) % revolved_size);
        }
    }

    return revolved_vertices;
}

GPU_Geometry generate_curve(std::vector<glm::vec3>& control_points, bool& bezier_mode, int& draw_count, int& scene_mode, std::vector<unsigned int>& indices) {
    CPU_Geometry cpu_geometry;
    GPU_Geometry gpu_geometry;

    cpu_geometry.verts = bezier_mode ? deCasteljau(control_points) : chaikin_subdivision(control_points, 6);
    draw_count = (int)cpu_geometry.verts.size();

    switch (scene_mode) {
        case 2:
            indices.clear();
            cpu_geometry.verts = surface_revolution(cpu_geometry.verts, indices);
            draw_count = (int)cpu_geometry.verts.size();
            break;
        case 3:
            indices.clear();
            cpu_geometry.verts = generate_chaikin_surface(control_points, indices);
            break;
    }
    cpu_geometry.cols.assign(cpu_geometry.verts.size(), glm::vec3(0.0f, 0.0f, 0.0f));
    update_gpu_geometry(gpu_geometry, cpu_geometry);
    return gpu_geometry;
}

void setup_gpu_geometry(CPU_Geometry& geometry, GPU_Geometry& gpu_points, GPU_Geometry& gpu_lines, int& selected_point) {
    geometry.cols.clear();
    geometry.cols.resize(geometry.verts.size(), glm::vec3{1.0f, 0.0f, 0.0f});
    
    if (selected_point >= 0 && selected_point < geometry.verts.size()) {
        geometry.cols[selected_point] = glm::vec3{0.0f, 0.0f, 1.0f};
    }
    update_gpu_geometry(gpu_points, geometry);

	// colour reset
    geometry.cols.clear();
    geometry.cols.resize(geometry.verts.size(), glm::vec3{0.0f, 1.0f, 0.0f});
    update_gpu_geometry(gpu_lines, geometry);
}



int main() {
	Log::debug("starting main");

	// WINDOW
	glfwInit();
	Window window(800, 800, "CPSC 453: Assignment 3");

	GLDebug::enable();

	auto curve_editor = std::make_shared<CurveEditor>();
	window.setCallbacks(curve_editor);


	float camera_yaw = -90.0f;
	float camera_pitch = 0.0f;


	glm::mat4 identity_matrix(1.0f);
	glm::mat4 projection_matrix = glm::perspective(45.0f, 1.0f, 0.1f, 100.0f);
	glm::vec3 camera_position(0.0f, 0.0f, 1.5f);
	glm::vec3 camera_direction(0.0f, 0.0f, -1.0f);


	CPU_Geometry control_points;
	GPU_Geometry point_gpu_geom, line_gpu_geom, curve_gpu_geom;
	int selected_point_index;
	int curve_draw_count;
	bool current_curve_mode = false;
	int current_scene_mode = 0;


	std::vector<unsigned int> draw_indices;
	unsigned int element_buffer;
	glGenBuffers(1, &element_buffer);

	ShaderProgram shader("shaders/test.vert", "shaders/test.frag");
	GLint mvp_location = glGetUniformLocation(GLuint(shader), "projection_matrix");

	reset_points(control_points, selected_point_index);
	setup_gpu_geometry(control_points, point_gpu_geom, line_gpu_geom, selected_point_index);
	curve_gpu_geom = generate_curve(control_points.verts, curve_editor->bezier_mode, curve_draw_count, curve_editor->scene_mode, draw_indices);

	glPointSize(10.0f);

	// the given tensor surface for part IV
	std::vector<glm::vec3> tensor_surface_1{
		{-2,0,-2}, {-1,0,-2}, {0,0,-2}, {1,0,-2}, {2,0,-2},
		{-2,0,-1}, {-1,1,-1}, {0,1,-1}, {1,1,-1}, {2,0,-1},
		{-2,0,0},  {-1,1,0},  {0,-1,0}, {1,1,0},  {2,0,0},
		{-2,0,1},  {-1,1,1},  {0,1,1},  {1,1,1},  {2,0,1},
		{-2,0,2},  {-1,0,2},  {0,0,2},  {1,0,2},  {2,0,2}
	};

	// I tried to make this one look interesting, but it kinda just looks like a hill.
	// I would try to do it algorithmically but honestly, I'm kinda tired so here ya go.
	std::vector<glm::vec3> tensor_surface_2{
    	{-3,  0.4, -3},  {-2,  0.6, -3},  {-1,  0.8, -3},   {0,  0.6, -3},   {1,  0.8, -3},   {2,  0.6, -3},   {3,  0.4, -3},
    	{-3,  0.6, -2},  {-2,  0.9, -2},  {-1,  1.2, -2},   {0,  1.0, -2},   {1,  1.2, -2},   {2,  0.9, -2},   {3,  0.6, -2},
    	{-3,  0.8, -1},  {-2,  1.2, -1},  {-1,  1.5, -1},   {0,  1.4, -1},   {1,  1.5, -1},   {2,  1.2, -1},   {3,  0.8, -1},
    	{-3,  0.6,  0},  {-2,  1.0,  0},  {-1,  1.4,  0},   {0,  1.5,  0},   {1,  1.4,  0},   {2,  1.0,  0},   {3,  0.6,  0},
    	{-3,  0.4,  1},  {-2,  0.8,  1},  {-1,  1.2,  1},   {0,  1.4,  1},   {1,  1.2,  1},   {2,  0.8,  1},   {3,  0.4,  1},
    	{-3,  0.3,  2},  {-2,  0.6,  2},  {-1,  0.9,  2},   {0,  1.1,  2},   {1,  0.9,  2},   {2,  0.6,  2},   {3,  0.3,  2},
    	{-3,  0.2,  3},  {-2,  0.4,  3},  {-1,  0.6,  3},   {0,  0.8,  3},   {1,  0.6,  3},   {2,  0.4,  3},   {3,  0.2,  3}
	};
	
	while (!window.shouldClose()) {
		glfwPollEvents();

		glEnable(GL_LINE_SMOOTH);
		glEnable(GL_FRAMEBUFFER_SRGB);
		glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		shader.use();

		// reset handling
		if (curve_editor->reset_requested) {
			reset_points(control_points, selected_point_index);
			curve_editor->reset_requested = false;
			camera_position = glm::vec3(0.0f, 0.0f, 1.5f);
			camera_direction = glm::vec3(0.0f, 0.0f, -1.0f);

			if (current_scene_mode != 3) {
				setup_gpu_geometry(control_points, point_gpu_geom, line_gpu_geom, selected_point_index);
				curve_gpu_geom = generate_curve(control_points.verts, curve_editor->bezier_mode, curve_draw_count, curve_editor->scene_mode, draw_indices);
			}
		}

		// scene mode 0 (2D curve editing)
		if (curve_editor->scene_mode == 0) {
			// handle left click (add/move points)
			if (curve_editor->left_click_active) {
				if (!curve_editor->left_click_acknowledged) {
					glm::vec3 click_pos = curve_editor->last_click_position;
					selected_point_index = 0;

					for (int i = 0; i < (int)control_points.verts.size(); i++) {
						if (glm::distance(click_pos, control_points.verts[i]) < glm::distance(click_pos, control_points.verts[selected_point_index])) {
							selected_point_index = i;
						}
					}

					if (glm::distance(control_points.verts[selected_point_index], click_pos) > 0.035f) {
						selected_point_index = (int)control_points.verts.size();
						control_points.verts.push_back(click_pos);
					}
				}
				curve_editor->left_click_acknowledged = true;

				// move point
				if (selected_point_index >= 0 && selected_point_index < (int)control_points.verts.size()) {
					control_points.verts[selected_point_index] = curve_editor->get_world_position_from_mouse();
				}

				setup_gpu_geometry(control_points, point_gpu_geom, line_gpu_geom, selected_point_index);
				curve_gpu_geom = generate_curve(control_points.verts, curve_editor->bezier_mode, curve_draw_count, curve_editor->scene_mode, draw_indices);
			}

			// handle deletion with right click or delete key
			if (curve_editor->delete_requested || curve_editor->right_click_active) {
				curve_editor->delete_requested = false;

				if (curve_editor->right_click_active) {
					glm::vec3 click_pos = curve_editor->last_click_position;
					int closest_point = 0;

					for (int i = 0; i < (int)control_points.verts.size(); i++) {
						if (glm::distance(click_pos, control_points.verts[i]) < glm::distance(click_pos, control_points.verts[closest_point])) {
							closest_point = i;
						}
					}

					if (glm::distance(control_points.verts[closest_point], click_pos) <= 0.035f) {
						selected_point_index = closest_point;
					}
				}

				curve_editor->right_click_active = false;

				if (selected_point_index >= 0 && selected_point_index < (int)control_points.verts.size()) {
					control_points.verts.erase(control_points.verts.begin() + selected_point_index);
					selected_point_index = -1;

					setup_gpu_geometry(control_points, point_gpu_geom, line_gpu_geom, selected_point_index);
					curve_gpu_geom = generate_curve(control_points.verts, curve_editor->bezier_mode, curve_draw_count, curve_editor->scene_mode, draw_indices);
				}
			}

			glUniformMatrix4fv(mvp_location, 1, GL_FALSE, &identity_matrix[0][0]);
		}
		else if (curve_editor->scene_mode > 0) {
			// camera orbit control
			if (curve_editor->left_click_active) {
				curve_editor->left_click_acknowledged = true;
				glm::ivec2 mouse_delta = curve_editor->consume_mouse_delta();
				camera_yaw -= (float)mouse_delta.x / 8.0f;
				camera_pitch += (float)mouse_delta.y / 8.0f;

				camera_pitch = glm::clamp(camera_pitch, -89.0f, 89.0f);

				camera_direction = glm::normalize(glm::vec3(
					cos(glm::radians(camera_yaw)) * cos(glm::radians(camera_pitch)),
					sin(glm::radians(camera_pitch)),
					sin(glm::radians(camera_yaw)) * cos(glm::radians(camera_pitch))
				));
			}

			// WASD camera movement
			if (curve_editor->move_forward)  camera_position += camera_direction * 0.05f;
			if (curve_editor->move_left)     camera_position -= glm::normalize(glm::cross(camera_direction, glm::vec3(0, 1, 0))) * 0.05f;
			if (curve_editor->move_backward) camera_position -= camera_direction * 0.05f;
			if (curve_editor->move_right)    camera_position += glm::normalize(glm::cross(camera_direction, glm::vec3(0, 1, 0))) * 0.05f;

			glm::mat4 view_proj_matrix = projection_matrix * glm::lookAt(camera_position, camera_position + camera_direction, glm::vec3(0, 1, 0));
			glUniformMatrix4fv(mvp_location, 1, GL_FALSE, &view_proj_matrix[0][0]);
		}

		// handle scene or curve mode changes
		if (current_curve_mode != curve_editor->bezier_mode || current_scene_mode != curve_editor->scene_mode) {
			current_curve_mode = curve_editor->bezier_mode;
			current_scene_mode = curve_editor->scene_mode;

			if (current_scene_mode == 3) {
				// Choose which tensor product surface to render based on the current curve mode
				// in other words, you can control the surface rendered with 'ENTER'
				curve_gpu_geom = current_curve_mode
					? generate_curve(tensor_surface_1, current_curve_mode, curve_draw_count, current_scene_mode, draw_indices)
					: generate_curve(tensor_surface_2, current_curve_mode, curve_draw_count, current_scene_mode, draw_indices);
			} else {
				curve_gpu_geom = generate_curve(control_points.verts, curve_editor->bezier_mode, curve_draw_count, curve_editor->scene_mode, draw_indices);
			}
		}

		// render geometry
		curve_gpu_geom.bind();
		if (curve_editor->scene_mode < 2) {
			glDrawArrays(GL_LINE_STRIP, 0, curve_draw_count);
		} else {
			// check if wireframe mode toggled on or off to render
			if (curve_editor->wireframe) {
				glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
			} else {
				glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
			}
			glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, element_buffer);
			glBufferData(GL_ELEMENT_ARRAY_BUFFER, draw_indices.size() * sizeof(unsigned int), draw_indices.data(), GL_STATIC_DRAW);
			glDrawElements(GL_TRIANGLES, (GLsizei)draw_indices.size(), GL_UNSIGNED_INT, nullptr);
		}

		// draw control lines and points
		if (curve_editor->scene_mode < 3) {
			line_gpu_geom.bind();
			glDrawArrays(GL_LINE_STRIP, 0, (GLsizei)control_points.verts.size());

			point_gpu_geom.bind();
			glDrawArrays(GL_POINTS, 0, (GLsizei)control_points.verts.size());
		}

		glDisable(GL_FRAMEBUFFER_SRGB);
		window.swapBuffers();
	}

	glfwTerminate();
	return 0;
}

// main.cpp - Adapted for CPSC 453 Assignment 3 (2025)
