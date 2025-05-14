# Bezier Curve visualizer

Author: Jarin Thundathil

## To Build
in main folder:
    1. mkdir build
    2. cd build
    3. cmake ..
    4. make
    5. ./bezier_visualizer

on visual studio:
    1. hit play once CMakeLists.txt is linked.

## Controls
**General**
    1.'SPACE' to toggle between four scenes: 2D Curve Editor, 2D Curve Viewer (orbital camera), 3D Curve Viewer, Tensor Product Surfaces.
    2.'ENTER' to toggle between the Bezier Curve and B-Spline. In the tensor Product Surfaces, press 'ENTER' to switch between surfaces.
    3.'WASD' moves the camera in scene 2, 3, and 4
    4.'LEFT-CLICK + DRAG' to change the camera orientation in scenes 2, 3, and 4
    5.'T' toggles wireframe on and off for Scenes 3 and 4. Use this to toggle between wireframe and solid surface.
    6. Press 'R' at any point to restart the program.

**Curve Editor**
    1.'LEFT-CLICK' on blank space to add new points. 'LEFT-CLICK' on a point and drag it to move it around.
    2.'RIGHT-CLICK' on a point to delete that point.
    3.'BACKSPACE' on a slected point to delete that point.

## Specifications
Compiler: clang
OS: developed on MacOS
