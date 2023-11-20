#include <CanvasTriangle.h>
#include <DrawingWindow.h>
#include <Utils.h>
#include <fstream>
#include <vector>
#include <glm/glm.hpp>
#include <CanvasPoint.h>
#include <Colour.h>
#include <TextureMap.h>
#include <TexturePoint.h>
#include <iostream>
#include <fstream>
#include <ModelTriangle.h>
#include <string>
#include <algorithm>
#include <unordered_map>
#include <cmath>
#include <SDL.h>
#include <RayTriangleIntersection.h>

#define WIDTH 320 * 2
#define HEIGHT 240 * 2
#define HALFWIDTH 160 * 2
#define HALFHEIGHT 120 * 2
#define OBJfilename "cornell-box.obj"
#define MTLfilename "cornell-box.mtl"

#ifndef MYRAY_HPP
#define MYRAY_HPP

/*
you know where the camera is (you have its position in a vec3); 
you know where the image plane is 
(it's a fixed distance from the camera and its orientation is fixed - 
relative to the camera). 
For a particular pixel on the image plane, 
it's fairly straight-forward to calculate the vector from the camera to that pixel. 
Below is an image that shows the image plane as a bunch of big square pixels 
- it's top-down and only in 2D, but you can use the same principle to get a y coordinate.

as George points out - Raytracing is fundamentally different 
(it turns rasterising on its head in many ways). 
There will need to be some scaling and WIDTH/2 and HEIGHT/2 
shifting at some point - often the inverse of what you did with rasterising. 
The best advice is to do things step-by-step: 
(1) get the small render showing in the top-left corner 
(2) scale it to make it bigger 
(3) shift it to the centre of the image plane. 
At each stage you will need to tinker with your code to make it work appropriately 
- if you try everything at once you will probably get into a pickle.

you shoot a ray from the camera, 
through a pixel in the image plane and into the model 
then (if you hit something) you then shoot another ray from the surface 
that was hit (you can calculate its location in 3D space) towards the light. 
If that second ray hits something closer than the light, 
the point on the surface is in shadow.

could be that the calculation of the intersected surface points are incorrect. 
When you are doing you shadow calculation, 
make sure that your surface is in the world coordinate system. 
Use both formulae in week 6 task 2 (just above the triangle) 
to calculate the x,y,z position of the intersection in world space 
- make sure the two match. 
Other things to try is a bit of "visual debugging" 
- when rendering, colour each pixel the distance that 
surface is from the origin (to check that all seems correct). 
If that looks OK, try colouring the pixels the distance they are from the light.
*/

// glm::vec3 e0 = triangle.vertices[1] - triangle.vertices[0];
// glm::vec3 e1 = triangle.vertices[2] - triangle.vertices[0];
// glm::vec3 SPVector = cameraPosition - triangle.vertices[0];
// glm::mat3 DEMatrix(-rayDirection, e0, e1);
// glm::vec3 possibleSolution = glm::inverse(DEMatrix) * SPVector;

class MyRay {
public:
    /*
    {   8 pixels   }
    [][][][][][][][] [][][][][][][][]
    \           |
     \          |
      \         |
       \        | Focal length: 2
        \       |
         \      |
          \     |
           \    |
            [camera]

    [][][][]
    [][][][]
    [][][][]
    [][][][]
    */
};

#endif
