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
#include <map>
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

    /*
    In ray trace task, it simulates the physical phenomenon of light reflection & light diffraction.
    All scene we see is a reflection of light. Light carries colour, and it reflects to human eyes.
    1. Find ray direction from camera/eye to the canvas. Inverse to getCanvasIntersectionPoint.
    2. Find the intersection point between camera and canvas.
    3. According to triangles of intersection point, draw it on the canvas.
    */

    std::vector<ModelTriangle> getModelTriangle(std::vector<glm::vec3> vertices, std::vector<std::size_t> faces, std::vector<ModelTriangle> vecModel, Colour c) {
        for (size_t i = 0; i < faces.size(); i+=3) { // 96 / 3 = 32
            glm::vec3 f1 = vertices[faces[i] - 1];
            glm::vec3 f2 = vertices[faces[i + 1] - 1];
            glm::vec3 f3 = vertices[faces[i + 2] - 1];
            vecModel.push_back(ModelTriangle(f1, f2, f3, c));
        }
        return vecModel; // [0] - [31]
    }
    std::vector<ModelTriangle> readOBJ(std::vector<ModelTriangle> vecModel) {
        std::ifstream inputStream(OBJfilename);
        std::string nextLine;
        std::vector<glm::vec3> vertices;
        std::vector<std::size_t> faces;
        Colour c;
        // If OBJ file load failed
        if (!inputStream.is_open()) {
            std::cerr << "Failed to open OBJ file!" << std::endl;
        }
        // Use a while loop together with the getline() function to read the file line by line
        while (std::getline(inputStream, nextLine)) {
            auto line = split(nextLine, ' '); // std::vector<std::string>
            for (size_t i = 0; i < line.size(); i++) { // iterate all line to locate line[i]
                if (line[i] == "v") { // 0
                    glm::vec3 tmp;
                    tmp.x = std::stof(line[i+1]) * 0.35; // 1
                    tmp.y = std::stof(line[i+2]) * 0.35; // 2
                    tmp.z = std::stof(line[i+3]) * 0.35; // 3
                    vertices.push_back(tmp);
                }
                if (line[i] == "f") { // 0
                    std::size_t tmp_1;
                    std::size_t tmp_2;
                    std::size_t tmp_3;
                    tmp_1 = std::stoi(line[i+1]); // 1
                    tmp_2 = std::stoi(line[i+2]); // 2
                    tmp_3 = std::stoi(line[i+3]); // 3
                    faces.push_back(tmp_1);
                    faces.push_back(tmp_2);
                    faces.push_back(tmp_3);
                }
            }
        }
        // Update vecModel
        vecModel = getModelTriangle(vertices, faces, vecModel, c);
        // Hardcoding to insert colour name inside ModelTriangle...
        vecModel[0].colour.name = "White";
        vecModel[1].colour.name = "White";
        vecModel[2].colour.name = "Grey";
        vecModel[3].colour.name = "Grey";
        vecModel[4].colour.name = "Cyan";
        vecModel[5].colour.name = "Cyan";
        vecModel[6].colour.name = "Green";
        vecModel[7].colour.name = "Green";
        vecModel[8].colour.name = "Magenta";
        vecModel[9].colour.name = "Magenta";
        vecModel[10].colour.name = "Yellow";
        vecModel[11].colour.name = "Yellow";
        for (size_t i = 12; i < 22; i++) {
            vecModel[i].colour.name = "Red";
        }
        for (size_t i = 22; i < vecModel.size(); i++) {
            vecModel[i].colour.name = "Blue";
        }
        // Close the file
        inputStream.close();
        return vecModel;
    }
    // std::unordered_map<std::string, uint32_t> myMap
    std::unordered_map<std::string, std::vector<float>> readMTL(std::unordered_map<std::string, std::vector<float>> myMap) {
        std::ifstream inputStream(MTLfilename);
        std::string nextLine;
        std::vector<float> colourBoard;
        std::vector<std::string> colourName;
        // If MTL file load failed
        if (!inputStream.is_open()) {
            std::cerr << "Failed to open MTL file!" << std::endl;
        }
        // Use a while loop together with the getline() function to read the file line by line
        while (std::getline(inputStream, nextLine)) {
            auto line = split(nextLine, ' '); // std::vector<std::string>
            for (size_t i = 0; i < line.size(); i++) {
                if (line[i] == "newmtl") { // 0
                    colourName.push_back(line[i+1]); // 1
                }
                if (line[i] == "Kd") { // 0
                    colourBoard.push_back(std::stof(line[i+1]) * 255); // 1
                    colourBoard.push_back(std::stof(line[i+2]) * 255); // 2
                    colourBoard.push_back(std::stof(line[i+3]) * 255); // 3
                }
            }
        }
        // Insert index and corresponding colour into the hashmap
        for (size_t i = 0; i < colourName.size(); i++) {
            for (size_t j = i * 3; j < (i * 3) + 3; j+=3) {
                // uint32_t colour = (255 << 24) + (int(colourBoard[j]) << 16) + (int(colourBoard[j+1]) << 8) + int(colourBoard[j+2]);
                // myMap[colourName[i]] = colour;
                myMap[colourName[i]] = {colourBoard[j], colourBoard[j+1], colourBoard[j+2]};
            }
        }
        // Close the file
        inputStream.close();
        return myMap;
    }
    std::vector<ModelTriangle> updateModelTriangleColour(std::unordered_map<std::string, std::vector<float>> myMap, std::vector<ModelTriangle> vecModel) {
        for (size_t i = 0; i < vecModel.size(); i++) {
            for (auto pair : myMap) { // const auto& pair : myMap
                if (vecModel[i].colour.name == pair.first) { // Colour name matched
                    // Insert RGB value into colour
                    vecModel[i].colour.red = int(pair.second[0]);
                    vecModel[i].colour.green = int(pair.second[1]);
                    vecModel[i].colour.blue = int(pair.second[2]);
                }
            }
        }
        return vecModel;
    }
    glm::vec3 getRayDirection(glm::vec3 cameraPosition, size_t y, size_t x) {
        // Reverse implementation function of getCanvasIntersectionPoint(glm::vec3 cameraPosition, glm::vec3 vertexPosition, float focalLength) function
        // But... it iterates all pixels of canvas, which means it loops through all pixels of canvas...
        // In this function, we need to find out the actual position of vertex, instead of canvas point
        CanvasPoint point;
        point.y = y;
        point.x = x;
        point.depth = 1.0f; // Assign the depth to 1: Buffer [0] to [1] (normalise) (top)
        glm::vec3 vertexPosition;
        // x: positive, y: negative (inverse implementation)
        vertexPosition.x = ((point.x - HALFWIDTH) / (80 * focalLength)) * point.depth; // Reverse formula
        vertexPosition.y = - ((point.y - HALFHEIGHT) / (80 * focalLength)) * point.depth; // Reverse formula
        vertexPosition.z = 0.0f * point.depth; // Initial value of z: 0.0f, reverse formula
        vertexPosition = vertexPosition - cameraPosition;
        vertexPosition = glm::normalize(vertexPosition); // Normalise to stabilise the range between 0 and 1
        return vertexPosition;
    }
    RayTriangleIntersection getClosestValidIntersection(glm::vec3 cameraPosition, glm::vec3 rayDirection, ModelTriangle triangle) {
        RayTriangleIntersection rayTriangle;
        return rayTriangle;
    }
    RayTriangleIntersection verifyClosestValidIntersection(glm::vec3 cameraPosition, glm::vec3 rayDirection, std::vector<ModelTriangle> vecModel) {
        RayTriangleIntersection rayTriangle;
        return rayTriangle;
    }
    void drawRayTrace(DrawingWindow& window, std::vector<ModelTriangle> vecModel) {
        for (size_t y = 0; y < window.height; y++) {
            for (size_t x = 0; x < window.width; x++) {
                // Draw ray trace scene test
            }
        }
    }
    void drawRasterisedSceneRayTrace(DrawingWindow& window) {
        std::unordered_map<std::string, std::vector<float>> myMap;
        std::vector<ModelTriangle> vecModel; // 32
        myMap = readMTL(myMap); // Read MTL file
        vecModel = readOBJ(vecModel); // Read OBJ file
        vecModel = updateModelTriangleColour(myMap, vecModel); // Map colour to vecModel
        std::cout << "Test ray trace is starting" << std::endl;
        drawRayTrace(window, vecModel);
        std::cout << "Test ray trace is completed" << std::endl;
    }

};

#endif
