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

Steps in the process:
	
1.Loop through all the pixels on the SDL canvas one at a time. For each pixel:
	
2.Calculate the vector from the camera to the current pixel 
(we know the camera position, we know the image plane is fixed relative to the camera, 
we know the distance to the image plane, we know the x and y on the image place)
	
3.Loop through each triangle in the loaded model
	
4.Use the code given in week 6 task 2 to see if the camera-to-pixel vector 
intersects with the current triangle
	
5.Validate the intersection as indicated in week 6 task 3
	
6.Continue through all the triangles in the loaded model to find the closest valid intersection
	
7.Use the colour of the closest validly intersected triangle to colour that current pixel
    */
    RayTriangleIntersection getClosestIntersection(glm::vec3 cameraPosition, glm::vec3 rayDirection, std::vector<ModelTriangle> vecModel) {
        RayTriangleIntersection rayTriangle;
        rayTriangle.distanceFromCamera = std::numeric_limits<float>::infinity(); // positive infinity
        // -std::numeric_limits<float>::infinity(); // negative infinity
        glm::vec3 ray = cameraPosition - rayDirection; // ray: camera to vertex
        ray = ray * cameraOrientation; // camera orientation applies to ray
        for (size_t i = 0; i < vecModel.size(); i++) {
            /*
            [t]   [-dx][e0x][e1x]-1    [sx - p0x]
            [u] = [-dy][e0y][e1y]   *  [sy - p0y]
            [v]   [-dz][e0z][e1z]      [sz - p0z]
            Another form:
            [tuv] = inverse(DEMatrix) * SPVector
            */
            // glm::vec3 e0 = triangle.vertices[1] - triangle.vertices[0];
            // glm::vec3 e1 = triangle.vertices[2] - triangle.vertices[0];
            // glm::vec3 SPVector = cameraPosition - triangle.vertices[0];
            // glm::mat3 DEMatrix(-rayDirection, e0, e1);
            // glm::vec3 possibleSolution = glm::inverse(DEMatrix) * SPVector;
            glm::vec3 e0 = vecModel[i].vertices[1] - vecModel[i].vertices[0];
            glm::vec3 e1 = vecModel[i].vertices[2] - vecModel[i].vertices[0];
            glm::vec3 SPVector = cameraPosition - vecModel[i].vertices[0];
            glm::mat3 DEMatrix(-ray, e0, e1); // (-rayDirection, e0, e1)
            glm::vec3 possibleSolution = glm::inverse(DEMatrix) * SPVector; // distance along ray and the coords on the triangle
            float t = possibleSolution.x;
            float u = possibleSolution.y;
            float v = possibleSolution.z;
            /*
                    [p1]

                    r       [p2]

            [p0]
            r = p0 + u * (p1 - p0) + v * (p2 - p0)
            */
            //--------------------
            // (u >= 0.0) && (u <= 1.0)
            // (v >= 0.0) && (v <= 1.0)
            // (u + v) <= 1.0
            //--------------------
            // position = startpoint + scalar * direction
            if (((u >= 0.0) && (u <= 1.0)) && ((v >= 0.0) && (v <= 1.0)) && ((u + v) <= 1.0)) { // validation process
                if ((rayTriangle.distanceFromCamera > t) && (t > 0)) {
                    rayTriangle.distanceFromCamera = t;
                    rayTriangle.intersectedTriangle = vecModel[i];
                    rayTriangle.triangleIndex = i;
                    glm::vec3 intersection = vecModel[i].vertices[0] + u * e0 + v * e1;
                    rayTriangle.intersectionPoint = intersection;
                }
            }
        }
        return rayTriangle;
    }
    bool isShadow(RayTriangleIntersection rayTriangle, std::vector<ModelTriangle> vecModel) {
        // (0.0, 0.8, 0.0) - (x, y, z)
        glm::vec3 shadowRay = lightPosition - rayTriangle.intersectionPoint; // light - intersection point
        for (size_t i = 0; i < vecModel.size(); i++) {

            glm::vec3 e0 = vecModel[i].vertices[1] - vecModel[i].vertices[0];
            glm::vec3 e1 = vecModel[i].vertices[2] - vecModel[i].vertices[0];
            glm::vec3 SPVector = rayTriangle.intersectionPoint - vecModel[i].vertices[0];
            glm::mat3 DEMatrix(-glm::normalize(shadowRay), e0, e1); // (-ray, e0, e1)
            glm::vec3 possibleSolution = glm::inverse(DEMatrix) * SPVector; // distance along ray and the coords on the triangle
            float t = possibleSolution.x;
            float u = possibleSolution.y;
            float v = possibleSolution.z;

            if (((u >= 0.0) && (u <= 1.0)) && ((v >= 0.0) && (v <= 1.0)) && ((u + v) <= 1.0)) { // validation process
                if ((t < glm::length(shadowRay)) && (t > 0.01) && (i != rayTriangle.triangleIndex)) { // light intersects with canvas
                    return true;
                }
            }
        }
        return false;
    }
    void drawRayTrace(DrawingWindow& window, std::vector<ModelTriangle> vecModel) {
        for (size_t y = 0; y < window.height; y++) {
            for (size_t x = 0; x < window.width; x++) {
                /* Index:
                [x] [y] [focalLength]
                rayDirection = -x, y, focalLength
                */
                glm::vec3 rayDirection = glm::vec3(-(x - HALFWIDTH), y - HALFHEIGHT, focalLength);
                RayTriangleIntersection rayTriangle = getClosestIntersection(cameraPosition, rayDirection, vecModel);
                if (!isinf(rayTriangle.distanceFromCamera)) { // distance from camera is not infinite
                    if (isShadow(rayTriangle, vecModel)) { // isShadow(rayTriangle, vecModel) == true
                        // Dark shadow
                        uint32_t colour = (255 << 24) + (0 << 16) + (0 << 8) + 0; // Alpha + Red + Green + Blue
                        // size_t stdX = std::round(x);
                        // size_t stdY = std::round(y);
                        // if ((stdX < WIDTH) && (stdY < HEIGHT)) { // Fix the range
                        //     window.setPixelColour(stdX, stdY, colour);
                        // }
                        window.setPixelColour(x, y, colour);
                    } else { // isShadow(rayTriangle, vecModel) != true
                        // Normal rasterised render output...
                    }
                }
            }
        }
    }
};

#endif
