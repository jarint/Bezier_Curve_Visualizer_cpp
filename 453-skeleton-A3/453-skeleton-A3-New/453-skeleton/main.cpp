#include <glad/glad.h>

#include "GLDebug.h"
#include "Log.h"
#include "Window.h"

#include "CurveControl.h"

int main() {
  Log::debug("Starting main");

  // WINDOW
  glfwInit();
  Window window(800, 800, "CPSC 453: Assignment 3");

  // GLDebug::enable();

  CurveControl curveControl(window);

  while (!window.shouldClose()) {
    glfwPollEvents();

	curveControl.Update();
    curveControl.DrawGeometry();

    window.swapBuffers();
  }

  glfwTerminate();
  return 0;
}
