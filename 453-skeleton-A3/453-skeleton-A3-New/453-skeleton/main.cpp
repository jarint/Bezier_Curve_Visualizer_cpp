// #include <glad/glad.h>

// #include "GLDebug.h"
// #include "Log.h"
// #include "Window.h"

// #include "CurveControl.h"

// int main() {
//   Log::debug("Starting main");

//   // WINDOW
//   glfwInit();
//   Window window(800, 800, "CPSC 453: Assignment 3");

//   // GLDebug::enable();

//   CurveControl curveControl(window);

//   while (!window.shouldClose()) {
//     glfwPollEvents();

// 	curveControl.Update();
//     curveControl.DrawGeometry();

//     window.swapBuffers();
//   }

//   glfwTerminate();
//   return 0;
// }

// main.cpp - Adapted for CPSC 453 Assignment 3 (2025)







/*
















*/






// #include <glad/glad.h>
// #include <GLFW/glfw3.h>
// #include <iostream>
// #include <vector>
// #include <cmath>
// #include <glm/glm.hpp>
// #include <glm/gtc/matrix_transform.hpp>
// #include <glm/gtc/type_ptr.hpp>
// #include "Window.h"
// #include "ShaderProgram.h"
// #include "Geometry.h"
// #include "GLDebug.h"
// #include "Log.h"

// // Utility to update GPU geometry from CPU geometry
// void updateGPUGeometry(GPU_Geometry &gpuGeom, CPU_Geometry const &cpuGeom) {
//     gpuGeom.bind();
//     gpuGeom.setVerts(cpuGeom.verts);
//     gpuGeom.setCols(cpuGeom.cols);
// }

// // Callback class to handle input & scene states
// class Assignment3 : public CallbackInterface {
// public:
//     bool LClicked, RClicked, deletePoint, resetScene;
//     bool W, A, S, D;
//     bool bezierMode;  // true for Bezier, false for B-Spline
//     int sceneMode;    // 0: 2D Curve Editor, 1: 3D Viewing, 2: Surface of Revolution, 3: Tensor Product Surface
//     glm::vec3 clickPosition;
//     glm::ivec2 screenDimensions;
//     glm::vec2 lastMousePos, mouseDelta;

//     Assignment3(int w = 800, int h = 800)
//         : LClicked(false), RClicked(false), deletePoint(false), resetScene(false),
//           W(false), A(false), S(false), D(false), bezierMode(true), sceneMode(0),
//           clickPosition(0.f), screenDimensions(w, h), lastMousePos(0), mouseDelta(0) {}

//     virtual void keyCallback(int key, int scancode, int action, int mods) {
//         if (action == GLFW_PRESS) {
//             switch (key) {
//                 case GLFW_KEY_BACKSPACE:
//                 case GLFW_KEY_DELETE: deletePoint = true; break;
//                 case GLFW_KEY_R: resetScene = true; break;
//                 case GLFW_KEY_SPACE: sceneMode = (sceneMode + 1) % 4; break;
//                 case GLFW_KEY_TAB: bezierMode = !bezierMode; break;
//                 case GLFW_KEY_W: W = true; break;
//                 case GLFW_KEY_A: A = true; break;
//                 case GLFW_KEY_S: S = true; break;
//                 case GLFW_KEY_D: D = true; break;
//             }
//         } else if (action == GLFW_RELEASE) {
//             if (key == GLFW_KEY_W) W = false;
//             if (key == GLFW_KEY_A) A = false;
//             if (key == GLFW_KEY_S) S = false;
//             if (key == GLFW_KEY_D) D = false;
//         }
//     }

//     virtual void mouseButtonCallback(int button, int action, int mods) {
//         if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_PRESS) {
//             LClicked = true;
//             clickPosition = getWorldPos(lastMousePos);
//         } else if (button == GLFW_MOUSE_BUTTON_RIGHT && action == GLFW_PRESS) {
//             RClicked = true;
//             clickPosition = getWorldPos(lastMousePos);
//         } else if (action == GLFW_RELEASE) {
//             LClicked = false;
//             RClicked = false;
//         }
//     }

//     virtual void cursorPosCallback(double xpos, double ypos) {
//         mouseDelta = glm::vec2(xpos - lastMousePos.x, ypos - lastMousePos.y);
//         lastMousePos = glm::vec2(xpos, ypos);
//     }

//     virtual void windowSizeCallback(int width, int height) {
//         CallbackInterface::windowSizeCallback(width, height);
//         screenDimensions = glm::ivec2(width, height);
//     }

//     glm::vec3 getWorldPos(glm::vec2 mousePos) {
//         return glm::vec3(
//             (mousePos.x / (float)screenDimensions.x) * 2.0f - 1.0f,
//             -((mousePos.y / (float)screenDimensions.y) * 2.0f - 1.0f), 0.f
//         );
//     }

//     glm::vec2 consumeMouseDelta() {
//         glm::vec2 delta = mouseDelta;
//         mouseDelta = glm::vec2(0);
//         return delta;
//     }
// };

// // The next section will contain geometry generation, drawing logic, and the main render loop (continued in next part).
// // Curve generation functions
// glm::vec3 lerp(const glm::vec3 &a, const glm::vec3 &b, float t) {
//     return (1 - t) * a + t * b;
// }

// std::vector<glm::vec3> deCasteljau(const std::vector<glm::vec3> &points) {
//     std::vector<glm::vec3> result;
//     for (float u = 0; u <= 1.0f; u += 0.01f) {
//         std::vector<glm::vec3> tempPoints = points;
//         for (int i = 1; i < (int)points.size(); i++) {
//             for (int j = 0; j < (int)points.size() - i; j++) {
//                 tempPoints[j] = lerp(tempPoints[j], tempPoints[j + 1], u);
//             }
//         }
//         result.push_back(tempPoints[0]);
//     }
//     return result;
// }

// std::vector<glm::vec3> chaikinSubdivision(const std::vector<glm::vec3> &points, int depth = 4) {
//     std::vector<glm::vec3> result = points;
//     for (int d = 0; d < depth; d++) {
//         std::vector<glm::vec3> newPoints;
//         newPoints.push_back(result[0]);
//         for (int i = 0; i < (int)result.size() - 1; i++) {
//             glm::vec3 Q = 0.75f * result[i] + 0.25f * result[i + 1];
//             glm::vec3 R = 0.25f * result[i] + 0.75f * result[i + 1];
//             newPoints.push_back(Q);
//             newPoints.push_back(R);
//         }
//         newPoints.push_back(result.back());
//         result = newPoints;
//     }
//     return result;
// }

// std::vector<glm::vec3> revolve(const std::vector<glm::vec3> &curvePoints, std::vector<unsigned int> &indices) {
//     std::vector<glm::vec3> revolvedPoints;
//     const int slices = 80;
//     const float step = 2 * M_PI / slices;

//     for (int i = 0; i < slices; i++) {
//         float angle = i * step;
//         for (const auto &p : curvePoints) {
//             revolvedPoints.emplace_back(p.x * cos(angle), p.y, p.x * sin(angle));
//         }
//     }

//     int profileSize = (int)curvePoints.size();
//     for (int ring = 0; ring < slices; ring++) {
//         int nextRing = (ring + 1) % slices;
//         for (int i = 0; i < profileSize - 1; i++) {
//             indices.push_back(ring * profileSize + i);
//             indices.push_back(nextRing * profileSize + i);
//             indices.push_back(nextRing * profileSize + i + 1);

//             indices.push_back(ring * profileSize + i);
//             indices.push_back(nextRing * profileSize + i + 1);
//             indices.push_back(ring * profileSize + i + 1);
//         }
//     }

//     return revolvedPoints;
// }

// // Next part: Tensor product surface generation, main loop, camera setup, and rendering logic coming up.
// // Tensor product surface generation
// std::vector<glm::vec3> flattenSurface(const std::vector<std::vector<glm::vec3>> &grid, std::vector<unsigned int> &indices) {
//     std::vector<glm::vec3> flatPoints;
//     int rows = (int)grid.size();
//     int cols = (int)grid[0].size();

//     for (const auto &row : grid) {
//         flatPoints.insert(flatPoints.end(), row.begin(), row.end());
//     }

//     for (int i = 0; i < rows - 1; i++) {
//         for (int j = 0; j < cols - 1; j++) {
//             unsigned int topLeft = i * cols + j;
//             unsigned int topRight = topLeft + 1;
//             unsigned int bottomLeft = topLeft + cols;
//             unsigned int bottomRight = bottomLeft + 1;

//             indices.insert(indices.end(), { topLeft, bottomLeft, bottomRight });
//             indices.insert(indices.end(), { topLeft, bottomRight, topRight });
//         }
//     }

//     return flatPoints;
// }

// std::vector<glm::vec3> generateTensorSurface(const std::vector<std::vector<glm::vec3>> &controlGrid, std::vector<unsigned int> &indices) {
//     std::vector<std::vector<glm::vec3>> refinedRows;

//     for (const auto &row : controlGrid) {
//         refinedRows.push_back(chaikinSubdivision(row, 4));
//     }

//     std::vector<std::vector<glm::vec3>> transposed;
//     int refinedCols = (int)refinedRows[0].size();
//     int refinedRowsCount = (int)refinedRows.size();
//     transposed.resize(refinedCols, std::vector<glm::vec3>(refinedRowsCount));

//     for (int i = 0; i < refinedRowsCount; i++) {
//         for (int j = 0; j < refinedCols; j++) {
//             transposed[j][i] = refinedRows[i][j];
//         }
//     }

//     std::vector<std::vector<glm::vec3>> refinedBoth;
//     for (const auto &col : transposed) {
//         refinedBoth.push_back(chaikinSubdivision(col, 4));
//     }

//     std::vector<std::vector<glm::vec3>> finalSurface;
//     int finalRows = (int)refinedBoth[0].size();
//     int finalCols = (int)refinedBoth.size();

//     for (int i = 0; i < finalRows; i++) {
//         std::vector<glm::vec3> row;
//         for (int j = 0; j < finalCols; j++) {
//             row.push_back(refinedBoth[j][i]);
//         }
//         finalSurface.push_back(row);
//     }

//     return flattenSurface(finalSurface, indices);
// }

// // Next: camera setup, shader initialization, and the main render loop logic.


// // Camera setup and main render loop

// int main() {
//     Log::debug("Starting Assignment 3 (2025)");
//     glfwInit();
//     Window window(800, 800, "CPSC 453 Assignment 3");
//     GLDebug::enable();

//     auto a3 = std::make_shared<Assignment3>();
//     window.setCallbacks(a3);

//     ShaderProgram shader("shaders/test.vert", "shaders/test.frag");
//     GLint mvpLocation = glGetUniformLocation(GLuint(shader), "mvp");

//     CPU_Geometry cpuPoints;
//     GPU_Geometry gpuPoints, gpuLines, gpuCurve, gpuSurface;
//     std::vector<unsigned int> surfaceIndices;
//     unsigned int ebo;
//     glGenBuffers(1, &ebo);

//     std::vector<std::vector<glm::vec3>> tensorGrid = {
//         { {-2,0,-2}, {-1,0,-2}, {0,0,-2}, {1,0,-2}, {2,0,-2} },
//         { {-2,0,-1}, {-1,1,-1}, {0,1,-1}, {1,1,-1}, {2,0,-1} },
//         { {-2,0,0}, {-1,1,0}, {0,-1,0}, {1,1,0}, {2,0,0} },
//         { {-2,0,1}, {-1,1,1}, {0,1,1}, {1,1,1}, {2,0,1} },
//         { {-2,0,2}, {-1,0,2}, {0,0,2}, {1,0,2}, {2,0,2} }
//     };

//     int selectedPoint = -1;

//     glm::vec3 cameraPos(0, 0, 4);
//     glm::vec3 cameraFront(0, 0, -1);
//     glm::vec3 cameraUp(0, 1, 0);
//     float yaw = -90.f, pitch = 0.f;
//     float fov = 45.f;

//     while (!window.shouldClose()) {
//         glfwPollEvents();

//         if (a3->resetScene) {
//             cpuPoints.verts = { {-0.5f, -0.5f, 0}, {0.5f, -0.5f, 0}, {0.5f, 0.5f, 0}, {-0.5f, 0.5f, 0} };
//             selectedPoint = -1;
//             a3->resetScene = false;
//         }

//         if (a3->LClicked) {
//             glm::vec3 clicked = a3->clickPosition;
//             float minDist = 0.05f;
//             int closest = -1;
//             for (int i = 0; i < (int)cpuPoints.verts.size(); i++) {
//                 if (glm::distance(clicked, cpuPoints.verts[i]) < minDist) {
//                     minDist = glm::distance(clicked, cpuPoints.verts[i]);
//                     closest = i;
//                 }
//             }
//             if (closest == -1) {
//                 cpuPoints.verts.push_back(clicked);
//             } else {
//                 cpuPoints.verts[closest] = clicked;
//             }
//             a3->LClicked = false;
//         }

//         if (a3->deletePoint && !cpuPoints.verts.empty()) {
//             cpuPoints.verts.pop_back();
//             a3->deletePoint = false;
//         }

//         // Camera orbit handling
//         if (a3->sceneMode > 0) {
//             glm::vec2 delta = a3->consumeMouseDelta();
//             yaw += delta.x * 0.1f;
//             pitch -= delta.y * 0.1f;
//             pitch = glm::clamp(pitch, -89.f, 89.f);

//             cameraFront = glm::normalize(glm::vec3(
//                 cos(glm::radians(yaw)) * cos(glm::radians(pitch)),
//                 sin(glm::radians(pitch)),
//                 sin(glm::radians(yaw)) * cos(glm::radians(pitch))
//             ));

//             if (a3->W) cameraPos += 0.05f * cameraFront;
//             if (a3->S) cameraPos -= 0.05f * cameraFront;
//             if (a3->A) cameraPos -= glm::normalize(glm::cross(cameraFront, cameraUp)) * 0.05f;
//             if (a3->D) cameraPos += glm::normalize(glm::cross(cameraFront, cameraUp)) * 0.05f;
//         }

//         glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
//         glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

//         shader.use();
//         glm::mat4 view = glm::lookAt(cameraPos, cameraPos + cameraFront, cameraUp);
//         glm::mat4 projection = glm::perspective(glm::radians(fov), 1.f, 0.1f, 100.f);
//         glm::mat4 model(1.f);
//         glm::mat4 mvp = projection * view * model;

//         glUniformMatrix4fv(mvpLocation, 1, GL_FALSE, glm::value_ptr(mvp));

//         // Next section will cover drawing each scene based on sceneMode
//         // Scene drawing logic
//         if (a3->sceneMode == 0) { // 2D curve editor
//             // Update points and curves
//             cpuPoints.cols.assign(cpuPoints.verts.size(), glm::vec3(1, 0, 0));
//             updateGPUGeometry(gpuPoints, cpuPoints);

//             CPU_Geometry lineGeom;
//             lineGeom.verts = cpuPoints.verts;
//             lineGeom.cols.assign(cpuPoints.verts.size(), glm::vec3(0, 1, 0));
//             updateGPUGeometry(gpuLines, lineGeom);

//             CPU_Geometry curveGeom;
//             curveGeom.verts = a3->bezierMode ? deCasteljau(cpuPoints.verts) : chaikinSubdivision(cpuPoints.verts);
//             curveGeom.cols.assign(curveGeom.verts.size(), glm::vec3(0, 0, 0));
//             updateGPUGeometry(gpuCurve, curveGeom);

//             // Draw points, lines, and curve
//             gpuPoints.bind();
//             glPointSize(10.0f);
//             glDrawArrays(GL_POINTS, 0, (GLsizei)cpuPoints.verts.size());

//             gpuLines.bind();
//             glDrawArrays(GL_LINE_STRIP, 0, (GLsizei)cpuPoints.verts.size());

//             gpuCurve.bind();
//             glDrawArrays(GL_LINE_STRIP, 0, (GLsizei)curveGeom.verts.size());

//         } else if (a3->sceneMode == 1) { // 3D viewing of curve
//             std::vector<glm::vec3> curve3D = a3->bezierMode ? deCasteljau(cpuPoints.verts) : chaikinSubdivision(cpuPoints.verts);
//             CPU_Geometry curve3DGeom;
//             curve3DGeom.verts = curve3D;
//             curve3DGeom.cols.assign(curve3D.size(), glm::vec3(0, 0, 0));
//             updateGPUGeometry(gpuCurve, curve3DGeom);

//             gpuCurve.bind();
//             glDrawArrays(GL_LINE_STRIP, 0, (GLsizei)curve3D.size());

//         } else if (a3->sceneMode == 2) { // Surface of revolution
//             surfaceIndices.clear();
//             auto revolvedSurface = revolve(cpuPoints.verts, surfaceIndices);
//             CPU_Geometry surfaceGeom;
//             surfaceGeom.verts = revolvedSurface;
//             surfaceGeom.cols.assign(revolvedSurface.size(), glm::vec3(0.5, 0.5, 0.5));
//             updateGPUGeometry(gpuSurface, surfaceGeom);

//             glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ebo);
//             glBufferData(GL_ELEMENT_ARRAY_BUFFER, surfaceIndices.size() * sizeof(unsigned int), surfaceIndices.data(), GL_STATIC_DRAW);

//             gpuSurface.bind();
//             glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
//             glDrawElements(GL_TRIANGLES, (GLsizei)surfaceIndices.size(), GL_UNSIGNED_INT, nullptr);
//             glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

//         } else if (a3->sceneMode == 3) { // Tensor product surface
//             surfaceIndices.clear();
//             auto tensorSurface = generateTensorSurface(tensorGrid, surfaceIndices);
//             CPU_Geometry tensorGeom;
//             tensorGeom.verts = tensorSurface;
//             tensorGeom.cols.assign(tensorSurface.size(), glm::vec3(0.6, 0.6, 0.6));
//             updateGPUGeometry(gpuSurface, tensorGeom);

//             glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ebo);
//             glBufferData(GL_ELEMENT_ARRAY_BUFFER, surfaceIndices.size() * sizeof(unsigned int), surfaceIndices.data(), GL_STATIC_DRAW);

//             gpuSurface.bind();
//             glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
//             glDrawElements(GL_TRIANGLES, (GLsizei)surfaceIndices.size(), GL_UNSIGNED_INT, nullptr);
//             glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
//         }

//         window.swapBuffers();
//     }

//     glfwTerminate();
//     return 0;
// }



/*
















*/

//#include <GL/glew.h>
//#include <GLFW/glfw3.h>
#include <glad/glad.h>

#include <iostream>
#include <string>
#include <list>
#include <vector>
#include <limits>
#include <functional>

#include "Geometry.h"
#include "GLDebug.h"
#include "Log.h"
#include "ShaderProgram.h"
#include "Shader.h"
#include "Texture.h"
#include "Window.h"

#include "glm/glm.hpp"
#include "glm/gtc/type_ptr.hpp"

const int SCENE_NUM = 4;


// We gave this code in one of the tutorials, so leaving it here too
void updateGPUGeometry(GPU_Geometry &gpuGeom, CPU_Geometry const &cpuGeom) {
	gpuGeom.bind();
	gpuGeom.setVerts(cpuGeom.verts);
	gpuGeom.setCols(cpuGeom.cols);
}

// EXAMPLE CALLBACKS
class Assignment3 : public CallbackInterface {

public:
	bool LClicked;
	bool RClicked;
	bool clickAck;
	bool deleteB;
	bool resetflag;
	bool curvemode;
	int scenemode;

	bool W;
	bool A;
	bool S;
	bool D;

	glm::vec3 mouseClickPos;

	Assignment3(int x, int y) :
		LClicked(false),
		RClicked(false),
		clickAck(false),
		deleteB(false),
		resetflag(false),
		curvemode(true),
		scenemode(0),
		W(false),
		A(false),
		S(false),
		D(false),
		mousePosRawDelta(0),
		mouseClickPos(0.f),
		mousePosRaw(0),
		screenDim(x,y)
	{}

	Assignment3() : Assignment3(800,800){}

	virtual void keyCallback(int key, int scancode, int action, int mods) {
		if(action == GLFW_PRESS){
			switch(key){
				case GLFW_KEY_BACKSPACE:
				case GLFW_KEY_DELETE:
					deleteB = true;
					break;
				case GLFW_KEY_W:
					W = true;
					break;
				case GLFW_KEY_A:
					A = true;
					break;
				case GLFW_KEY_S:
					S = true;
					break;
				case GLFW_KEY_D:
					D = true;
					break;

				case GLFW_KEY_R:
					resetflag = true;
					break;
				case GLFW_KEY_TAB:
					curvemode = !curvemode;
					break;
				case GLFW_KEY_SPACE:
					scenemode = (scenemode + 1) % SCENE_NUM;
					break;

				default:
					break;
			}
		}
		if(action == GLFW_RELEASE){
			switch(key){
				case GLFW_KEY_W:
					W = false;
					break;
				case GLFW_KEY_A:
					A = false;
					break;
				case GLFW_KEY_S:
					S = false;
					break;
				case GLFW_KEY_D:
					D = false;
					break;
				
				default:
					break;
			}
		}
	}
	virtual void mouseButtonCallback(int button, int action, int mods) {
		if(action == GLFW_PRESS){
			switch (button){
				case GLFW_MOUSE_BUTTON_LEFT:
					LClicked = true;
					mouseClickPos = translateMousePos();
					break;
				case GLFW_MOUSE_BUTTON_RIGHT:
					RClicked = true;
					mouseClickPos = translateMousePos();
					break;
				
				default:
					break;
			}
		}
		if(action == GLFW_RELEASE){
			switch (button){
				case GLFW_MOUSE_BUTTON_LEFT:
					LClicked = false;
					clickAck = false;
					break;
				case GLFW_MOUSE_BUTTON_RIGHT:
					RClicked = false;
					break;
				
				default:
					break;
			}
		}
	}
	virtual void cursorPosCallback(double xpos, double ypos) {
		mousePosRawDelta = glm::ivec2(xpos - mousePosRaw.x,ypos - mousePosRaw.y);
		mousePosRaw.x = xpos;
		mousePosRaw.y = ypos;
	}
	// virtual void scrollCallback(double xoffset, double yoffset) {
	// }
	virtual void windowSizeCallback(int width, int height) {
		// The CallbackInterface::windowSizeCallback will call glViewport for us
		CallbackInterface::windowSizeCallback(width,  height);
		screenDim.x = width;
		screenDim.y = height;
	}

	// void acknowledge(){
	// 	clicked = false;
	// 	deleteB = false;
	// }

	glm::vec3 translateMousePos(){
		return glm::vec3(((float)mousePosRaw.x * 2.f / screenDim.x) - 1.f,-(((float)mousePosRaw.y* 2.f / screenDim.y) - 1.f), 0.f);
	}

	glm::ivec2 mouseDeltaRaw(){
		glm::ivec2 temp = mousePosRawDelta;
		mousePosRawDelta = glm::ivec2(0);
		return temp;
	}

private:
	glm::vec2 mousePosRawDelta;
	glm::vec2 mousePosRaw;
	glm::ivec2 screenDim;

};

void reset(CPU_Geometry& points, int& selPoint){
	selPoint = -1;
	points.cols.clear();
	points.verts.clear();
	points.verts.push_back(glm::vec3{-0.5, 0.5, 0});
	points.verts.push_back(glm::vec3{-0.5, -0.5, 0});
	points.verts.push_back(glm::vec3{0.5, -0.5, 0});
	points.verts.push_back(glm::vec3{0.5, 0.5, 0});
}

void setupGPU(CPU_Geometry& points, GPU_Geometry& pointsGPUGeom, GPU_Geometry& linesGPUGeom, int& selPoint){
	points.cols.clear();
	points.cols.resize(points.verts.size(), glm::vec3{1.0, 0.0, 0.0});
	if(selPoint >= 0 && selPoint < points.verts.size()){
		points.cols[selPoint] = glm::vec3{0.0, 0.0, 1.0};
	}
	updateGPUGeometry(pointsGPUGeom, points);

	// Reset the colors to green
	points.cols.clear();
	points.cols.resize(points.verts.size(), glm::vec3{0.0, 1.0, 0.0});
	updateGPUGeometry(linesGPUGeom, points);
}

std::vector<glm::vec3> deCasteljau(std::vector<glm::vec3>& points){
	std::vector<glm::vec3> deCastPoints;
	int d = points.size();
	for(float u = 0; u <= 1; u+=0.01){
		std::vector<glm::vec3> pointsCP = points;
		for (int i = 1; i < d; i++){
			for (int j = 0; j < d-1; j++){
				pointsCP[j] = (1-u) * pointsCP[j] + u * pointsCP[j+1];
			}
		}
		deCastPoints.push_back(pointsCP[0]);
	}
	// deCastPoints.cols.resize(deCastPoints.verts.size(), glm::vec3{0.0, 0.0, 0.0});
	return deCastPoints;
}

std::vector<glm::vec3> chaikinPoints(std::vector<glm::vec3>& points){
	std::vector<glm::vec3> chaiPoints;
	int n = points.size();
	chaiPoints.push_back(points[0]);
	chaiPoints.push_back(points[0]*0.5f+points[1]*0.5f);
	for(int i = 1; i < n-2; i++){
		chaiPoints.push_back((points[i]*0.75f)+(points[i+1]*0.25f));
		chaiPoints.push_back((points[i]*0.25f)+(points[i+1]*0.75f));
	}
	chaiPoints.push_back(((points.end()[-2])*0.5f)+(points.end()[-1]*0.5f));
	chaiPoints.push_back(points.end()[-1]);
	return chaiPoints;
}

std::vector<glm::vec3> chaikin(std::vector<glm::vec3>& points, const int& divlevel){
	std::vector<glm::vec3> cPoints;

	cPoints = points;

	for (int i = 0; i < divlevel; i++){
		cPoints = chaikinPoints(cPoints);
	}
	return cPoints;
}

std::vector<glm::vec3> flattenPoints(std::vector<std::vector<glm::vec3>>& points, std::vector<unsigned int>& indices){
	std::vector<glm::vec3> flat;
	for(int i = 0; i < points.size(); i++){
		for(int j = 0; j < points[0].size(); j++){
			flat.push_back(points[i][j]);
		}
	}
	for(int j = 0; j < flat.size()-points.size(); j+=points.size()){
		for (int i = 0; i < points.size()-1; i++){
			indices.push_back((i+j));
			indices.push_back((i+j+points.size()));
			indices.push_back((i+1+j+points.size()));
			indices.push_back((i+j));
			indices.push_back((i+1+j));
			indices.push_back((i+1+j+points.size()));
		}
	}
	return flat;
}

std::vector<glm::vec3> chaikinSurface(std::vector<glm::vec3>& points, std::vector<unsigned int>& indices){
	int n = sqrt(points.size()); 
	std::vector<std::vector<glm::vec3>> surface1;
	std::vector<std::vector<glm::vec3>> surface2;
	std::vector<glm::vec3> newPoints;
	for(int i = 0; i < n; i++){
		newPoints.clear();
		for (int j = i; j < points.size(); j+=n){
			newPoints.push_back(points[j]);
		}
		surface1.push_back(chaikin(newPoints, 4));
	}
	for(int i = 0; i < surface1[0].size(); i++){
		newPoints.clear();
		for(int j = 0; j < n; j++){
			newPoints.push_back(surface1[j][i]);
		}
		surface2.push_back(chaikin(newPoints, 4));
	}
	std::vector<glm::vec3> finalP = flattenPoints(surface2, indices);
	return finalP;
}

std::vector<glm::vec3> revolve(std::vector<glm::vec3>& points, std::vector<unsigned int>& indices){
	std::vector<glm::vec3> newPoints;
	double step = M_PI/40;
	for(double u = 0; u < M_PI * 2; u+=step){
		for(int i = 0; i < points.size(); i++){
			newPoints.push_back(glm::vec3(
				points[i].x*cos(u),
				points[i].y,
				points[i].x*sin(u)
			));
		}
	}
	int w = newPoints.size();
	for(int j = 0; j < w; j+=points.size()){
		for (int i = 0; i < points.size()-1; i++){
			indices.push_back((i+j)%w);
			indices.push_back((i+j+points.size())%w);
			indices.push_back((i+1+j+points.size())%w);
			indices.push_back((i+j)%w);
			indices.push_back((i+1+j)%w);
			indices.push_back((i+1+j+points.size())%w);
		}
	}
	
	return newPoints;
}

GPU_Geometry drawCurve(std::vector<glm::vec3>& points, bool& mode, int& drawcount, int& scenemode,std::vector<unsigned int>& indices){
	CPU_Geometry CPU;
	GPU_Geometry GPU;
	if(mode){
		drawcount = 100;
		CPU.verts = deCasteljau(points);
	}else{
		CPU.verts = chaikin(points, 6);
		drawcount = CPU.verts.size();
	}
	if(scenemode==2){
		indices.clear();
		CPU.verts = revolve(CPU.verts,indices);
		drawcount = CPU.verts.size();
	}
	if(scenemode==3){
		indices.clear();
		CPU.verts = chaikinSurface(points, indices);
	}
	CPU.cols.resize(CPU.verts.size(), glm::vec3{0.0, 0.0, 0.0});
	updateGPUGeometry(GPU, CPU);
	return GPU;
}



int main() {
	Log::debug("Starting main");

	// WINDOW
	glfwInit();
	Window window(800, 800, "CPSC 453"); // can set callbacks at construction if desired


	GLDebug::enable();

	// CALLBACKS
	auto a3 = std::make_shared<Assignment3>();
	window.setCallbacks(a3);

	glm::mat4 defaultMat(1);
	glm::mat4 projMat = glm::perspective(45.0f,1.0f,0.1f,100.0f);
	glm::vec3 camLoc(0,0,1.5);
	glm::vec3 camLook(0,0,-1);
	float yaw = -90.f;
	float pitch = 0.f;

	ShaderProgram shader("shaders/test.vert", "shaders/test.frag");
	GLint viewProjMatLoc = glGetUniformLocation(GLuint(shader), "viewProjMat");

	std::vector<unsigned int> indices;
	unsigned int EBO;
	glGenBuffers(1, &EBO);

	std::vector<glm::vec3> pt4{
		{-2, 0,-2},{-1, 0,-2},{ 0, 0,-2},{ 1, 0,-2},{ 2, 0,-2},
		{-2, 0,-1},{-1, 1,-1},{ 0, 1,-1},{ 1, 1,-1},{ 2, 0,-1},
		{-2, 0, 0},{-1, 1, 0},{ 0,-1, 0},{ 1, 1, 0},{ 2, 0, 0},
		{-2, 0, 1},{-1, 1, 1},{ 0, 1, 1},{ 1, 1, 1},{ 2, 0, 1},
		{-2, 0, 2},{-1, 0, 2},{ 0, 0, 2},{ 1, 0, 2},{ 2, 0, 2}
	};

	std::vector<glm::vec3> pt4_2{
		{-3, 2,-3},{-1, 3,-4},{ 1, 1,-3},{ 3, 3,-2},{ 4, 3,-2},
		{-5,-1,-2},{-2, 1,-1},{ 0, 0,-1},{ 2,-2,-1},{ 3, 1,-1},
		{-4, 4,-1},{-3, 0, 0},{-1,-1, 1},{ 1, 4, 0},{ 3, 2, 0},
		{-3, 0, 1},{-1, 4, 2},{ 1, 2, 2},{ 2, 1, 1},{ 3, 0, 1},
		{-2,-1, 3},{-1, 5, 4},{ 0, 3, 3},{ 1,-3, 2},{ 4, 2, 3}
	};

	// The current CPU_Geometry and GPU_Geometry classes are defined in
	// Geometry.h/Geometry.cpp They will work for this assignment, but for some of
	// the bonuses you may have to modify them.
	CPU_Geometry points;
	GPU_Geometry pointsGPUGeom;
	GPU_Geometry linesGPUGeom;
	GPU_Geometry curvesGPUGeom;
	int selPoint;
	int drawCount;
	bool currentcurve = false;
	int currentscene = 0;
	reset(points, selPoint);
	setupGPU(points, pointsGPUGeom, linesGPUGeom, selPoint);
	curvesGPUGeom = drawCurve(points.verts,a3->curvemode,drawCount,a3->scenemode,indices);


	glPointSize(10.0f);

	// Note this call only work on some systems, unfortunately.
	// In order for them to work, you have to comment out line 60
	// If you're on a mac, you can't comment out line 60, so you
	// these will have no effect. :(
	// glLineWidth(5.0f);

	// RENDER LOOP
	while (!window.shouldClose()) {
		glfwPollEvents();

		glEnable(GL_LINE_SMOOTH);
		glEnable(GL_FRAMEBUFFER_SRGB);
		glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		shader.use();

		if(a3->resetflag){
			reset(points, selPoint);
			a3->resetflag = false;
			camLoc = glm::vec3(0,0,1.5);
			camLook = glm::vec3(0,0,-1);
			if(currentscene!=3){
				setupGPU(points, pointsGPUGeom, linesGPUGeom, selPoint);
				curvesGPUGeom = drawCurve(points.verts,currentcurve,drawCount,a3->scenemode,indices);
			}
		}

		if(a3->scenemode == 0){
			if(a3->LClicked){
				if(!a3->clickAck){
					glm::vec3 clickedPos = a3->mouseClickPos;
					selPoint = 0;
					for (int i = 0; i < points.verts.size(); i++){
						if(fabsf(glm::distance(clickedPos,points.verts[i])) < fabsf(glm::distance(clickedPos,points.verts[selPoint]))){
							selPoint = i;
						}
					}

					if(fabsf(glm::distance(points.verts[selPoint],clickedPos)) > 0.035f){
						selPoint = points.verts.size();
						points.verts.push_back(clickedPos);
					}
				}
				a3->clickAck = true;

				if(selPoint >= 0 && selPoint < points.verts.size()){
					points.verts[selPoint] = a3->translateMousePos();
				}

				setupGPU(points,pointsGPUGeom,linesGPUGeom, selPoint);
				curvesGPUGeom = drawCurve(points.verts,currentcurve,drawCount,a3->scenemode,indices);
				}

			if(a3->deleteB || a3->RClicked){
				a3->deleteB = false;
				if(a3->RClicked){
					glm::vec3 clickedPos = a3->mouseClickPos;
					int localSelPoint = 0;
					for (int i = 0; i < points.verts.size(); i++){
						if(fabsf(glm::distance(clickedPos,points.verts[i])) < fabsf(glm::distance(clickedPos,points.verts[localSelPoint]))){
							localSelPoint = i;
						}
					}
					if(fabsf(glm::distance(points.verts[localSelPoint],clickedPos)) <= 0.035f){
						selPoint = localSelPoint;
					}
				}
				a3->RClicked = false;
				if(selPoint >= 0 && selPoint < points.verts.size()){
					points.verts.erase(points.verts.begin()+selPoint);
					selPoint = -1;
					setupGPU(points,pointsGPUGeom,linesGPUGeom, selPoint);
					curvesGPUGeom = drawCurve(points.verts,currentcurve,drawCount,a3->scenemode,indices);
				}
			}

			glUniformMatrix4fv(viewProjMatLoc, 1, GL_FALSE, &defaultMat[0][0]);
		}
		else if(a3->scenemode > 0){
			if(a3->LClicked){
				a3->clickAck = true;
				glm::ivec2 delta = a3->mouseDeltaRaw();
				yaw -= ((float)delta.x) / 8;
				pitch += ((float)delta.y) / 8;
				if(pitch>89.f)pitch=89.f;
				if(pitch<-89.f)pitch=-89.f;
				camLook = glm::normalize(glm::vec3(cos(glm::radians(yaw)) * cos(glm::radians(pitch)), sin(glm::radians(pitch)), sin(glm::radians(yaw)) * cos(glm::radians(pitch))));
			}
			if(a3->W){
				camLoc += camLook * 0.05f;
			}
			if(a3->S){
				camLoc -= camLook * 0.05f;
			}
			if(a3->A){
				camLoc -= glm::normalize(glm::cross(camLook, glm::vec3(0,1,0))) * 0.05f;
			}
			if(a3->D){
				camLoc += glm::normalize(glm::cross(camLook, glm::vec3(0,1,0))) * 0.05f;
			}
			glm::mat4 viewProjMat = projMat * glm::lookAt(camLoc, camLoc + camLook, glm::vec3(0,1,0));
			glUniformMatrix4fv(viewProjMatLoc, 1, GL_FALSE, &viewProjMat[0][0]);
		}

		if(currentcurve!=a3->curvemode || currentscene!=a3->scenemode){
			currentcurve = a3->curvemode;
			currentscene =a3->scenemode;
			if(currentscene==3){
				if(currentcurve)curvesGPUGeom = drawCurve(pt4,currentcurve,drawCount,a3->scenemode,indices);
				else curvesGPUGeom = drawCurve(pt4_2,currentcurve,drawCount,a3->scenemode,indices);
			}
			else
				curvesGPUGeom = drawCurve(points.verts,currentcurve,drawCount,a3->scenemode,indices);
		}



		curvesGPUGeom.bind();
		if(a3->scenemode<2)
			glDrawArrays(GL_LINE_STRIP, 0, drawCount);
		else{
			glPolygonMode( GL_FRONT_AND_BACK, GL_LINE );
			glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
			glBufferData(GL_ELEMENT_ARRAY_BUFFER, indices.size()*sizeof(unsigned int), indices.data(), GL_STATIC_DRAW);
			glDrawElements(GL_TRIANGLES, indices.size(), GL_UNSIGNED_INT, 0);
			// glDrawArrays(GL_POINTS, 0, drawCount);
		}


		if(a3->scenemode<3){
			linesGPUGeom.bind();
			glDrawArrays(GL_LINE_STRIP, 0, GLsizei(points.verts.size()));

			pointsGPUGeom.bind();
			glDrawArrays(GL_POINTS, 0, GLsizei(points.verts.size()));
		}
		


		glDisable(GL_FRAMEBUFFER_SRGB); // disable sRGB for things like imgui

		window.swapBuffers();
	}

	glfwTerminate();
	return 0;
}
