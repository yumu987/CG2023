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

// using namespace std;
// using namespace glm;

#define WIDTH 320 * 2
#define HEIGHT 240 * 2
#define HALFWIDTH 160 * 2
#define HALFHEIGHT 120 * 2
#define OBJfilename "cornell-box.obj"
#define MTLfilename "cornell-box.mtl"

// Global variable: depthBuffer
float depthBuffer[WIDTH][HEIGHT];
// Global variable: cameraPosition (0, 0, 4)
glm::vec3 cameraPosition (0.0f, 0.0f, 4.0f);
// Global variable: focalLength 2
float focalLength = 2.0f;
// Global variable: targetPoint (0, 0, 0)
glm::vec3 targetPoint (0.0f, 0.0f, 0.0f);
// Global variable: upVector (0, 1, 0)
glm::vec3 upVector (0.0f, 1.0f, 0.0f);

glm::vec3 upMatrix = glm::vec3(0.0f, -0.1f, 0.0f); // y-
glm::vec3 downMatrix = glm::vec3(0.0f, 0.1f, 0.0f); // y+
glm::vec3 rightMatrix = glm::vec3(-0.1f, 0.0f, 0.0f); // x-
glm::vec3 leftMatrix = glm::vec3(0.1f, 0.0f, 0.0f); // x+

double angleInDegrees = 0.1;
double angleInRadians = angleInDegrees * (M_PI / 180); // radians = degrees * (M_PI / 180)

glm::mat3 rotationMatrixX = glm::mat3( // X axis, left
	glm::vec3(1.0f, 0.0f, 0.0f),
	glm::vec3(0.0f, std::cos(angleInRadians), -std::sin(angleInRadians)),
	glm::vec3(0.0f, std::sin(angleInRadians), std::cos(angleInRadians))
);
glm::mat3 rotationMatrixXX = glm::mat3( // X axis, right
	glm::vec3(1.0f, 0.0f, 0.0f),
	glm::vec3(0.0f, std::cos(angleInRadians), std::sin(angleInRadians)),
	glm::vec3(0.0f, -std::sin(angleInRadians), std::cos(angleInRadians))
);
glm::mat3 rotationMatrixY = glm::mat3( // Y axis, up
	glm::vec3(std::cos(angleInRadians), 0.0f, std::sin(angleInRadians)),
	glm::vec3(0.0f, 1.0f, 0.0f),
	glm::vec3(-std::sin(angleInRadians), 0.0f, std::cos(angleInRadians))
);
glm::mat3 rotationMatrixYY = glm::mat3( // Y axis, down
	glm::vec3(std::cos(angleInRadians), 0.0f, -std::sin(angleInRadians)),
	glm::vec3(0.0f, 1.0f, 0.0f),
	glm::vec3(std::sin(angleInRadians), 0.0f, std::cos(angleInRadians))
);
glm::mat3 rotationMatrixZ = glm::mat3( // Z axis, but no need to implement
	glm::vec3(std::cos(angleInRadians), -std::sin(angleInRadians), 0.0f),
	glm::vec3(std::sin(angleInRadians), std::cos(angleInRadians), 0.0f),
	glm::vec3(0.0f, 0.0f, 1.0f)
);
glm::mat3 innerProductIdentityMatrix = glm::mat3(
	glm::vec3(1.0f, 1.0f, 1.0f),
	glm::vec3(1.0f, 1.0f, 1.0f),
	glm::vec3(1.0f, 1.0f, 1.0f)
);
glm::mat3 outerProductIdentityMatrix = glm::mat3( // glm::mat3(1.0f)
	glm::vec3(1.0f, 0.0f, 0.0f),
	glm::vec3(0.0f, 1.0f, 0.0f),
	glm::vec3(0.0f, 0.0f, 1.0f)
);

// Global variable: cameraOrientation
glm::mat3 cameraOrientation = outerProductIdentityMatrix;

float zoomIn = 0.1f;
float zoomOut = -0.1f;

// Global variable: lightPosition (0, 0.5, 0.5)
glm::vec3 lightPosition = glm::vec3(0.0f, 0.5f, 0.5f); // (0.0f, 0.8f, 0.0f)

// --------------------
// Global rotation matrix for ray trace render
// --------------------

glm::vec3 upMatrixRay = glm::vec3(0.0f, -0.5f, 0.0f); // y-
glm::vec3 downMatrixRay = glm::vec3(0.0f, 0.5f, 0.0f); // y+
glm::vec3 rightMatrixRay = glm::vec3(-0.5f, 0.0f, 0.0f); // x-
glm::vec3 leftMatrixRay = glm::vec3(0.5f, 0.0f, 0.0f); // x+

double angleInDegreesRay = 2.5;
double angleInRadiansRay = angleInDegreesRay * (M_PI / 180); // radians = degrees * (M_PI / 180)

glm::mat3 rotationMatrixXRay = glm::mat3( // X axis, left
	glm::vec3(1.0f, 0.0f, 0.0f),
	glm::vec3(0.0f, std::cos(angleInRadiansRay), -std::sin(angleInRadiansRay)),
	glm::vec3(0.0f, std::sin(angleInRadiansRay), std::cos(angleInRadiansRay))
);
glm::mat3 rotationMatrixXXRay = glm::mat3( // X axis, right
	glm::vec3(1.0f, 0.0f, 0.0f),
	glm::vec3(0.0f, std::cos(angleInRadiansRay), std::sin(angleInRadiansRay)),
	glm::vec3(0.0f, -std::sin(angleInRadiansRay), std::cos(angleInRadiansRay))
);
glm::mat3 rotationMatrixYRay = glm::mat3( // Y axis, up
	glm::vec3(std::cos(angleInRadiansRay), 0.0f, std::sin(angleInRadiansRay)),
	glm::vec3(0.0f, 1.0f, 0.0f),
	glm::vec3(-std::sin(angleInRadiansRay), 0.0f, std::cos(angleInRadiansRay))
);
glm::mat3 rotationMatrixYYRay = glm::mat3( // Y axis, down
	glm::vec3(std::cos(angleInRadiansRay), 0.0f, -std::sin(angleInRadiansRay)),
	glm::vec3(0.0f, 1.0f, 0.0f),
	glm::vec3(std::sin(angleInRadiansRay), 0.0f, std::cos(angleInRadiansRay))
);

#ifndef MYMODEL_HPP
#define MYMODEL_HPP

class MyModel {
public:
    void resetModel() {
        cameraPosition = glm::vec3(0.0f, 0.0f, 4.0f); // Reset position
        focalLength = 2.0f; // Reset focal length
        cameraOrientation = outerProductIdentityMatrix; // Reset orientation
        lightPosition = glm::vec3(0.0f, 0.5f, 0.5f); // Reset light
    }
    void resetDepthBuffer() {
        for (size_t x = 0; x < WIDTH; x++) {
            for (size_t y = 0; y < HEIGHT; y++) {
                depthBuffer[x][y] = 0.0f;
            }
        }
    }
    std::vector<CanvasPoint> interpolateElements(CanvasPoint from, CanvasPoint to, float numberOfSteps) {
        std::vector<CanvasPoint> coordinates;
        float xDiff = to.x - from.x;
        float yDiff = to.y - from.y;
        float zDiff = to.depth - from.depth;
        float xStepSize = xDiff / numberOfSteps;
        float yStepSize = yDiff / numberOfSteps;
        float zStepSize = zDiff / numberOfSteps;
        for (float i = 0.0; i < numberOfSteps; i++) {
            CanvasPoint point;
            point.x = from.x + (xStepSize * i);
            point.y = from.y + (yStepSize * i);
            point.depth = from.depth + (zStepSize * i);
            coordinates.push_back(point);
        }
        return coordinates;
    }
    void drawPoint(DrawingWindow &window, CanvasPoint point) {
        uint32_t colour = (255 << 24) + (255 << 16) + (255 << 8) + 255; // White colour
        window.setPixelColour(point.x, point.y, colour);
    }
    void drawLine(DrawingWindow& window, CanvasPoint from, CanvasPoint to, Colour c) {
        float xDiff = to.x - from.x;
        float yDiff = to.y - from.y;
        float numberOfSteps = std::fmax(std::abs(xDiff), std::abs(yDiff));
        float xStepSize = xDiff / numberOfSteps;
        float yStepSize = yDiff / numberOfSteps;
        uint32_t colour = (255 << 24) + (c.red << 16) + (c.green << 8) + c.blue;
        for (float i = 0.0; i < numberOfSteps; i++) {
            float x = from.x + (xStepSize * i);
            float y = from.y + (yStepSize * i);
            size_t stdX = std::round(x);
            size_t stdY = std::round(y);
            if ((stdX < WIDTH) && (stdY < HEIGHT)) { // Fix the range to prevent out of range / segmentation fault
                // if ((stdX >= 0) && (stdY >= 0)) {}
                window.setPixelColour(stdX, stdY, colour);
            } // else {continue;}
        }
    }
    CanvasTriangle bubbleSort(CanvasTriangle triangle) {
        /*
        [v0]: top
        [v1]: middle
        [v2]: bottom
        */
        // sort height only
        if (triangle.v0().y > triangle.v1().y) {
            std::swap(triangle.v0(), triangle.v1());
        }
        if (triangle.v1().y > triangle.v2().y) {
            std::swap(triangle.v1(), triangle.v2());
            // After new v1 has been updated, we need to check v0 and v1 again for double confirmation...
            if (triangle.v0().y > triangle.v1().y) {
                std::swap(triangle.v0(), triangle.v1());
            }
        }
        return triangle;
    }
    CanvasPoint updateExtraPoint(CanvasPoint from, CanvasPoint to, CanvasPoint extra) {
        // from: triangle.v0(), to: triangle.v2()
        // Line equation: y = mx + b to update extra.x
        float xDiff = to.x - from.x;
        float yDiff = to.y - from.y;
        float m = yDiff / xDiff;
        float b = to.y - (m * to.x);
        if (xDiff == 0) { // Slope is 0
            // extra.x = from.x;
            extra.x = to.x;
        } else if (yDiff == 0) { // Slope is 0
            // extra.x = from.x;
            extra.x = to.x;
        } else if (xDiff == 0 && yDiff == 0) { // Slope is both 0
            // extra.x = from.x;
            extra.x = to.x;
        } else { // xDiff != 0 && yDiff != 0
            extra.x = (extra.y - b) / m;
        }
        // Proportion method to update extra.depth
        /*
                [v0]
        [extra] [v1] (extra.y = v1.y)
                [v2]
        yDiff = v2.y - v0.y
        t = v2.y - extra.y
        Formula: extra.depth = (t/yDiff) * v0.depth + (1 - (t/yDiff)) * v2.depth;
        The reason to use 1 (1.0f) is: t/yDiff to save upper part, (1-t/yDiff) to save lower part
        */
        float t = to.y - extra.y;
        t /= yDiff; // t = t / yDiff;
        extra.depth = t * from.depth + (1.0f - t) * to.depth;
        return extra;
    }
    void drawLineDepthBuffer(DrawingWindow& window, CanvasPoint from, CanvasPoint to, Colour c) {
        float xDiff = to.x - from.x;
        float yDiff = to.y - from.y;
        float zDiff = to.depth - from.depth;
        float numberOfSteps = std::fmax(std::abs(xDiff), std::abs(yDiff));
        float xStepSize = xDiff / numberOfSteps;
        float yStepSize = yDiff / numberOfSteps;
        float zStepSize = zDiff / numberOfSteps;
        uint32_t colour = (255 << 24) + (c.red << 16) + (c.green << 8) + c.blue;
        for (float i = 0.0; i < numberOfSteps; i++) {
            float x = from.x + (xStepSize * i);
            float y = from.y + (yStepSize * i);
            float z = from.depth + (zStepSize * i);
            // Ensure that depth buffer has the same index implementation with draw pixel
            size_t stdX = std::round(x);
            size_t stdY = std::round(y);
            if ((stdX < WIDTH) && (stdY < HEIGHT)) { // Fix the range to prevent out of range (overstep the boundary) / segmentation fault
                // if ((stdX >= 0) && (stdY >= 0)) {}
                if (z >= depthBuffer[stdX][stdY]) { // Depth is larger than recorded depth
                    depthBuffer[stdX][stdY] = z; // Record depth buffer
                    window.setPixelColour(stdX, stdY, colour); // Draw
                } // else {continue;}
            }
        }
    }
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
        // Print out vecModel which is used to contain model's data (debug method)
        // for (size_t i = 0; i < vecModel.size(); i++) {
        // 	std::cout << "vecModel[i]: " << i << std::endl;
        // 	std::cout << vecModel[i] << std::endl; // [0] - [31]
        // }
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
        // Iterate through the unordered_map (debug method)
        // for (const auto& pair : myMap) {
        //     std::cout << "Colour: " << pair.first << ", Red Value: " << pair.second[0]<< ", Green Value: " << pair.second[1]<< ", Blue Value: " << pair.second[2] << std::endl;
        // }
        // Close the file
        inputStream.close();
        return myMap;
    }
    CanvasPoint getCanvasIntersectionPoint(glm::vec3 cameraPosition, glm::vec3 vertexPosition, float focalLength) {
        /*
        Formulae:
        U^i = f * (X^i / Z^i) + W / 2
        V^i = f * (Y^i / Z^i) + H / 2
        Triangle:
            []
            [][] h^v
        []h^i[][]
        d^i
        d^v
        h^v / d^v = h^i / d^i
        Multiplier: 240 (It can be adjusted)
        A single minus (-) is showed to change the direction of x
        Minus (-) will reverse the direction
        Visual system:
        [camera] -> [plane] -> [vertex]
        distance    close       far
        */
        CanvasPoint point;
        float tmpX = vertexPosition.x - cameraPosition.x;
        float tmpY = vertexPosition.y - cameraPosition.y;
        float tmpZ = vertexPosition.z - cameraPosition.z;
        point.x = 160 * focalLength * ( - tmpX / tmpZ) + HALFWIDTH; // 160: constant factor
        point.y = 160 * focalLength * (tmpY / tmpZ) + HALFHEIGHT; // 160: constant factor
        // If z(1/Z) is greater, it is close to camera
        // If z(1/Z) is smaller, it is far from camera
        point.depth = 1 / - tmpZ;
        return point;
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
    void fillTopModelTriangle(DrawingWindow& window, CanvasPoint top, CanvasPoint middle, CanvasPoint extra, Colour c) {
        std::vector<CanvasPoint> firstLine;
        firstLine = interpolateElements(top, middle, middle.y - top.y);
        std::vector<CanvasPoint> secondLine;
        secondLine = interpolateElements(top, extra, extra.y - top.y);
        for (size_t i = 0; i < firstLine.size(); i++) { // i < secondLine.size()
            drawLineDepthBuffer(window, firstLine[i], secondLine[i], c);
        }
    }
    void fillBottomModelTriangle(DrawingWindow& window, CanvasPoint bottom, CanvasPoint middle, CanvasPoint extra, Colour c) {
        std::vector<CanvasPoint> firstLine;
        firstLine = interpolateElements(middle, bottom, bottom.y - middle.y);
        std::vector<CanvasPoint> secondLine;
        secondLine = interpolateElements(extra, bottom, bottom.y - extra.y);
        for (size_t i = 0; i < firstLine.size(); i++) { // i < secondLine.size()
            drawLineDepthBuffer(window, firstLine[i], secondLine[i], c);
        }
    }
    void filledModelTriangle(DrawingWindow& window, CanvasTriangle triangle, Colour c) {
        // Initialisation of points and sorting algorithm
        CanvasPoint extra;
        triangle = bubbleSort(triangle); // Sort top, middle and bottom
        extra.y = triangle.v1().y; // Update y
        extra = updateExtraPoint(triangle.v0(), triangle.v2(), extra); // Update x and depth
        // Draw top triangle
        fillTopModelTriangle(window, triangle.v0(), triangle.v1(), extra, c); // Top point, middle point and extra point
        // Draw bottom triangle
        fillBottomModelTriangle(window, triangle.v2(), triangle.v1(), extra, c); // Bottom point, middle point and extra point
    }
    void wireFrameRender(DrawingWindow& window) {
        // std::unordered_map<std::string, uint32_t> myMap;
        std::unordered_map<std::string, std::vector<float>> myMap;
        std::vector<ModelTriangle> vecModel; // 32
        std::vector<CanvasPoint> vecPoint; // 96
        std::vector<Colour> vecColour; // 32
        myMap = readMTL(myMap);
        vecModel = readOBJ(vecModel);
        vecModel = updateModelTriangleColour(myMap, vecModel); // Map colour to vecModel
        for (size_t i = 0; i < vecModel.size(); i++) {
            for (size_t j = 0; j < vecModel[i].vertices.size(); j++) {
                glm::vec3 vertexPosition = vecModel[i].vertices[j];
                CanvasPoint point = getCanvasIntersectionPoint(cameraPosition, vertexPosition, focalLength);
                vecPoint.push_back(point); // Store canvas point
                // drawPoint(window, point);
            }
        }
        for (size_t i = 0; i < vecModel.size(); i++) {
            vecColour.push_back(vecModel[i].colour); // Store colour
        }
        for (size_t i = 0; i < vecPoint.size(); i+=3) {
            CanvasTriangle triangle;
            Colour c;
            // Extract points
            triangle.v0() = vecPoint[i];
            triangle.v1() = vecPoint[i+1];
            triangle.v2() = vecPoint[i+2];
            // Extract colour
            c = vecColour[i/3];
            // Draw edge
            drawLine(window, triangle.v0(), triangle.v1(), c);
            drawLine(window, triangle.v1(), triangle.v2(), c);
            drawLine(window, triangle.v2(), triangle.v0(), c);
            // drawLineDepthBuffer(window, triangle.v0(), triangle.v1(), c);
            // drawLineDepthBuffer(window, triangle.v1(), triangle.v2(), c);
            // drawLineDepthBuffer(window, triangle.v2(), triangle.v0(), c);
        }
    }
    void rasterisedRender(DrawingWindow& window) {
        // std::unordered_map<std::string, uint32_t> myMap;
        std::unordered_map<std::string, std::vector<float>> myMap;
        std::vector<ModelTriangle> vecModel; // 32
        std::vector<CanvasPoint> vecPoint; // 96
        std::vector<Colour> vecColour; // 32
        myMap = readMTL(myMap);
        vecModel = readOBJ(vecModel);
        vecModel = updateModelTriangleColour(myMap, vecModel); // Map colour to vecModel
        for (size_t i = 0; i < vecModel.size(); i++) {
            for (size_t j = 0; j < vecModel[i].vertices.size(); j++) {
                glm::vec3 vertexPosition = vecModel[i].vertices[j];
                CanvasPoint point = getCanvasIntersectionPoint(cameraPosition, vertexPosition, focalLength);
                vecPoint.push_back(point); // Store canvas point
                // drawPoint(window, point);
            }
        }
        for (size_t i = 0; i < vecModel.size(); i++) {
            // Store colour
            vecColour.push_back(vecModel[i].colour);
        }
        for (size_t i = 0; i < vecPoint.size(); i+=3) {
            CanvasTriangle triangle;
            Colour c;
            // Extract points
            triangle.v0() = vecPoint[i];
            triangle.v1() = vecPoint[i+1];
            triangle.v2() = vecPoint[i+2];
            // Extract colour
            c = vecColour[i/3];
            // Draw triangle
            filledModelTriangle(window, triangle, c);
            // Draw edge
            drawLineDepthBuffer(window, triangle.v0(), triangle.v1(), c);
            drawLineDepthBuffer(window, triangle.v1(), triangle.v2(), c);
            drawLineDepthBuffer(window, triangle.v2(), triangle.v0(), c);
        }
    }
    //--------------------
    // Orbit
    //--------------------
    /*
    glm::normalise: The normalize function is used to ensure that the vector has a length of 1.
    Normalization is applied to ensure a unit length.
    Normalise: Fix the range (length) -> [0] - [1]

    glm::cross: The cross product gives a vector that is perpendicular to both input vectors.
    Cross product:                          Visual version of diagram of cross product:
            ^ perpendicular vector               ^
            |                                    |    perpendicular vector
            |                                    |
        <-      vector 1                <--------|    vector 1
            |                                   /
            v vector 2                         /      vector 2
                                              v
    rotation matrix: right, up, forward
                ^ up
                |
            <- forward
                /
               v right
    */
    glm::mat3 lookAt(glm::vec3 eye, glm::vec3 target, glm::vec3 up) { // const glm::vec3& eye, const glm::vec3& target, const glm::vec3& up
        glm::vec3 forwardVector = glm::normalize(eye - target); // Forward
        glm::vec3 rightVector = glm::normalize(glm::cross(up, forwardVector)); // Right
        glm::vec3 upVector = glm::normalize(glm::cross(forwardVector, rightVector)); // Up
        // glm::mat3 rotationMatrix;
        // rotationMatrix[0] = rightV;
        // rotationMatrix[1] = upV;
        // rotationMatrix[2] = forwardV;
        glm::mat3 rotationMatrix = glm::mat3(
            glm::vec3(rightVector.x, rightVector.y, rightVector.z),
            glm::vec3(upVector.x, upVector.y, upVector.z),
            glm::vec3(forwardVector.x, forwardVector.y, forwardVector.z)
        );
        return rotationMatrix;
    }
    CanvasPoint getCanvasIntersectionPointOrbit(glm::vec3 cameraPosition, glm::vec3 vertexPosition, float focalLength) {
        /*
        cameraToVertex = vertex - camera
        adjustedVector = cameraToVertex * cameraOrientation
        Cross product (*)
        The cross product gives a vector that is perpendicular to both input vectors
        [camera] -> [plane] -> [vertex]
        */
        CanvasPoint point;
        glm::vec3 cameraToVertex (vertexPosition.x - cameraPosition.x, vertexPosition.y - cameraPosition.y, vertexPosition.z - cameraPosition.z); // The distance between vertex and camera
        glm::vec3 adjustedVector = cameraToVertex * cameraOrientation; // orientation of camera applied
        point.x = 160 * focalLength * ( - adjustedVector.x / adjustedVector.z) + HALFWIDTH; // constant factor
        point.y = 160 * focalLength * (adjustedVector.y / adjustedVector.z) + HALFHEIGHT; // constant factor
        point.depth = 1 / - adjustedVector.z;
        return point;
    }
    void wireFrameRenderOrbit(DrawingWindow& window) {
        // std::unordered_map<std::string, uint32_t> myMap;
        std::unordered_map<std::string, std::vector<float>> myMap;
        std::vector<ModelTriangle> vecModel; // 32
        std::vector<CanvasPoint> vecPoint; // 96
        std::vector<Colour> vecColour; // 32
        myMap = readMTL(myMap);
        vecModel = readOBJ(vecModel);
        vecModel = updateModelTriangleColour(myMap, vecModel); // Map colour to vecModel
        for (size_t i = 0; i < vecModel.size(); i++) {
            for (size_t j = 0; j < vecModel[i].vertices.size(); j++) {
                glm::vec3 vertexPosition = vecModel[i].vertices[j];
                CanvasPoint point = getCanvasIntersectionPointOrbit(cameraPosition, vertexPosition, focalLength);
                vecPoint.push_back(point); // Store canvas point
                // drawPoint(window, point);
            }
        }
        for (size_t i = 0; i < vecModel.size(); i++) {
            vecColour.push_back(vecModel[i].colour); // Store colour
        }
        for (size_t i = 0; i < vecPoint.size(); i+=3) {
            CanvasTriangle triangle;
            Colour c;
            // Extract points
            triangle.v0() = vecPoint[i];
            triangle.v1() = vecPoint[i+1];
            triangle.v2() = vecPoint[i+2];
            // Extract colour
            c = vecColour[i/3];
            // Draw edge
            drawLine(window, triangle.v0(), triangle.v1(), c);
            drawLine(window, triangle.v1(), triangle.v2(), c);
            drawLine(window, triangle.v2(), triangle.v0(), c);
            // drawLineDepthBuffer(window, triangle.v0(), triangle.v1(), c);
            // drawLineDepthBuffer(window, triangle.v1(), triangle.v2(), c);
            // drawLineDepthBuffer(window, triangle.v2(), triangle.v0(), c);
        }
    }
    void rasterisedRenderOrbit(DrawingWindow& window) {
        // std::unordered_map<std::string, uint32_t> myMap;
        std::unordered_map<std::string, std::vector<float>> myMap;
        std::vector<ModelTriangle> vecModel; // 32
        std::vector<CanvasPoint> vecPoint; // 96
        std::vector<Colour> vecColour; // 32
        myMap = readMTL(myMap);
        vecModel = readOBJ(vecModel);
        vecModel = updateModelTriangleColour(myMap, vecModel); // Map colour to vecModel
        for (size_t i = 0; i < vecModel.size(); i++) {
            for (size_t j = 0; j < vecModel[i].vertices.size(); j++) {
                glm::vec3 vertexPosition = vecModel[i].vertices[j];
                CanvasPoint point = getCanvasIntersectionPointOrbit(cameraPosition, vertexPosition, focalLength);
                vecPoint.push_back(point); // Store canvas point
                // drawPoint(window, point);
            }
        }
        for (size_t i = 0; i < vecModel.size(); i++) {
            // Store colour
            vecColour.push_back(vecModel[i].colour);
        }
        for (size_t i = 0; i < vecPoint.size(); i+=3) {
            CanvasTriangle triangle;
            Colour c;
            // Extract points
            triangle.v0() = vecPoint[i];
            triangle.v1() = vecPoint[i+1];
            triangle.v2() = vecPoint[i+2];
            // Extract colour
            c = vecColour[i/3];
            // Draw triangle
            filledModelTriangle(window, triangle, c);
            // Draw edge
            drawLineDepthBuffer(window, triangle.v0(), triangle.v1(), c);
            drawLineDepthBuffer(window, triangle.v1(), triangle.v2(), c);
            drawLineDepthBuffer(window, triangle.v2(), triangle.v0(), c);
        }
    }

    //--------------------
    // Ray trace
    //--------------------

    //--------------------
    // Notes of ray trace
    //--------------------
    // glm::vec3 e0 = triangle.vertices[1] - triangle.vertices[0];
    // glm::vec3 e1 = triangle.vertices[2] - triangle.vertices[0];
    // glm::vec3 SPVector = cameraPosition - triangle.vertices[0];
    // glm::mat3 DEMatrix(-rayDirection, e0, e1);
    // glm::vec3 possibleSolution = glm::inverse(DEMatrix) * SPVector;
    //--------------------
    /*
    [t]   [-dx][e0x][e1x]-1    [sx - p0x]
    [u] = [-dy][e0y][e1y]   *  [sy - p0y]
    [v]   [-dz][e0z][e1z]      [sz - p0z]
    Another form:
    [tuv] = inverse(DEMatrix) * SPVector
    */
    //--------------------
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
    //--------------------

    // RayTriangleIntersection getClosestIntersection(glm::vec3 cameraPosition, glm::vec3 rayDirection, ModelTriangle modelTriangle) {
    //     // glm::vec3 e0 = triangle.vertices[1] - triangle.vertices[0];
    //     // glm::vec3 e1 = triangle.vertices[2] - triangle.vertices[0];
    //     // glm::vec3 SPVector = cameraPosition - triangle.vertices[0];
    //     // glm::mat3 DEMatrix(-rayDirection, e0, e1);
    //     // glm::vec3 possibleSolution = glm::inverse(DEMatrix) * SPVector;
    //     RayTriangleIntersection rayTriangle;
    //     // Direction: to - from
    //     // [t]          ([-d e0 e1])   [s - p0] (x)
    //     // [u] = inverse([-d e0 e1]) * [s - p0] (y)
    //     // [v]          ([-d e0 e1])   [s - p0] (z)
    //     glm::vec3 e0 = modelTriangle.vertices[1] - modelTriangle.vertices[0]; // direction e0
    //     glm::vec3 e1 = modelTriangle.vertices[2] - modelTriangle.vertices[0]; // direction e1
    //     glm::vec3 SPVector = cameraPosition - modelTriangle.vertices[0];
    //     glm::mat3 DEMatrix(-rayDirection, e0, e1);
    //     glm::vec3 possibleSolution = glm::inverse(DEMatrix) * SPVector;
    //     // (t, u, v) = (x, y, z)
    //     float t = possibleSolution.x; // the absolute distance along the ray from the camera to the intersection point
    //     float u = possibleSolution.y; // the proportional distance along the triangle's first edge that the intersection point occurs
    //     float v = possibleSolution.z; // the proportional distance along the triangle's second edge that the intersection point occurs
    //     // r = p0 + u(p1-p0) + v(p2-p0)
    //     // position = startpoint + scalar * direction
    //     rayTriangle.intersectionPoint = modelTriangle.vertices[0] + u * e0 + v * e1;
    //     rayTriangle.distanceFromCamera = t;
    //     rayTriangle.triangleIndex = t;
    //     return rayTriangle;
    // }

    // Main
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
        vertexPosition.x = (((point.x - HALFWIDTH) / (45 * focalLength)) * point.depth) - cameraPosition.x; // Reverse formula
        vertexPosition.y = (- ((point.y - HALFHEIGHT) / (45 * focalLength)) * point.depth) - cameraPosition.y; // Reverse formula
        vertexPosition.z = (0.0f * point.depth) - cameraPosition.z; // Initial value of z: 0.0f, reverse formula
        vertexPosition = glm::normalize(vertexPosition); // Normalise to stabilise the range between 0 and 1
        return vertexPosition;
    }
    RayTriangleIntersection getClosestValidIntersection(glm::vec3 cameraPosition, glm::vec3 rayDirection, ModelTriangle triangle) {
        // glm::vec3 e0 = triangle.vertices[1] - triangle.vertices[0];
        // glm::vec3 e1 = triangle.vertices[2] - triangle.vertices[0];
        // glm::vec3 SPVector = cameraPosition - triangle.vertices[0];
        // glm::mat3 DEMatrix(-rayDirection, e0, e1);
        // glm::vec3 possibleSolution = glm::inverse(DEMatrix) * SPVector;

        RayTriangleIntersection rayTriangle;
        
        // Direction: to - from
        // [t]          ([-d e0 e1])   [s - p0] (x)
        // [u] = inverse([-d e0 e1]) * [s - p0] (y)
        // [v]          ([-d e0 e1])   [s - p0] (z)
        glm::vec3 e0 = triangle.vertices[1] - triangle.vertices[0]; // direction e0
        glm::vec3 e1 = triangle.vertices[2] - triangle.vertices[0]; // direction e1
        glm::vec3 SPVector = cameraPosition - triangle.vertices[0];
        glm::mat3 DEMatrix(-rayDirection, e0, e1);
        glm::vec3 possibleSolution = glm::inverse(DEMatrix) * SPVector;

        // (t, u, v) = (x, y, z)
        float t = possibleSolution.x; // the absolute distance along the ray from the camera to the intersection point
        float u = possibleSolution.y; // the proportional distance along the triangle's first edge that the intersection point occurs
        float v = possibleSolution.z; // the proportional distance along the triangle's second edge that the intersection point occurs

        // position = startpoint + scalar * direction
        // rayTriangle.intersectionPoint = cameraPosition + t * rayDirection;

        // r = p0 + u(p1-p0) + v(p2-p0)
        // position = startpoint + scalar * direction + scalar * direction
        // intersection point = startpoint + first edge distance * first direction + second edge distance * second direction
        rayTriangle.intersectionPoint = triangle.vertices[0] + u * e0 + v * e1;
        rayTriangle.distanceFromCamera = t;
        rayTriangle.intersectedTriangle = triangle;
        // rayTriangle.triangleIndex = 0;
        rayTriangle.triangleIndex = std::numeric_limits<size_t>::infinity(); // positive infinity

        // (u >= 0.0) && (u <= 1.0)
        // (v >= 0.0) && (v <= 1.0)
        // (u + v) <= 1.0
        if (((u >= 0.0) && (u <= 1.0)) && ((v >= 0.0) && (v <= 1.0)) && ((u + v) <= 1.0)) { // validation process
            if (t > 0) { // the absolute distance along the ray from the camera to the intersection point should be positive
                rayTriangle.triangleIndex = 1; // It indicates that ray is intersected with triangle
            }
        }
        
        return rayTriangle;
    }
    RayTriangleIntersection verifyClosestValidIntersection(glm::vec3 cameraPosition, glm::vec3 rayDirection, std::vector<ModelTriangle> vecModel) {
        RayTriangleIntersection validTriangle;
        // RayTriangleIntersection validTriangle = RayTriangleIntersection(glm::vec3(), INT_MAX, ModelTriangle(), 0);
        validTriangle.distanceFromCamera = std::numeric_limits<float>::infinity(); // positive infinity
        RayTriangleIntersection rayTriangle;
        for (size_t i = 0; i < vecModel.size(); i++) {
            rayTriangle = getClosestValidIntersection(cameraPosition, rayDirection, vecModel[i]);
            if ((validTriangle.distanceFromCamera > rayTriangle.distanceFromCamera) && (rayTriangle.triangleIndex == 1)) { // rayTriangle.triangleIndex == 1 indicates valid intersected triangle
                validTriangle = rayTriangle;
                validTriangle.triangleIndex = i; // Remove this line will remove shadow
            }
        }
        return validTriangle;
    }
    void drawRayTrace(DrawingWindow& window, std::vector<ModelTriangle> vecModel) {
        float red = 0.0f;
        float green = 0.0f;
        float blue = 0.0f;
        uint32_t black = (255 << 24) + (int(red) << 16) + (int(green) << 8) + int(blue);
        uint32_t colour;
        for (size_t y = 0; y < window.height; y++) { // HEIGHT
            for (size_t x = 0; x < window.width; x++) { // WIDTH
                glm::vec3 rayDirection = getRayDirection(cameraPosition, y, x);
                RayTriangleIntersection rayTriangle = verifyClosestValidIntersection(cameraPosition, rayDirection, vecModel);
                RayTriangleIntersection lightTriangle = verifyClosestValidIntersection(lightPosition, glm::normalize(rayTriangle.intersectionPoint - lightPosition), vecModel);
                // The match of index of triangles means that the light will reflect the colour of triangles
                if (rayTriangle.triangleIndex == lightTriangle.triangleIndex) {
                    // uint32_t colour = (255 << 24) + (int(rayTriangle.intersectedTriangle.colour.red) << 16) + (int(rayTriangle.intersectedTriangle.colour.green) << 8) + int(rayTriangle.intersectedTriangle.colour.blue);
                    // size_t stdX = std::round(x);
                    // size_t stdY = std::round(y);
                    // if ((stdX >= 0) && (stdX < WIDTH) && (stdY >= 0) && (stdY < HEIGHT)) { // Fix the range of window
                    //     window.setPixelColour(stdX, stdY, colour);
                    // }
                    colour = (255 << 24) + (int(rayTriangle.intersectedTriangle.colour.red) << 16) + (int(rayTriangle.intersectedTriangle.colour.green) << 8) + int(rayTriangle.intersectedTriangle.colour.blue);
                } else { // The mismatch of index of triangles means that the light will not reflect the colour of light (shadow/darkness)
                    // size_t stdX = std::round(x);
                    // size_t stdY = std::round(y);
                    // if ((stdX >= 0) && (stdX < WIDTH) && (stdY >= 0) && (stdY < HEIGHT)) { // Fix the range of window
                    //     window.setPixelColour(stdX, stdY, black);
                    // }
                    colour = black;
                }
                window.setPixelColour(x, y, colour);
            }
        }
    }
    void drawRasterisedScene(DrawingWindow& window) {
        std::unordered_map<std::string, std::vector<float>> myMap;
        std::vector<ModelTriangle> vecModel; // 32
        myMap = readMTL(myMap); // Read MTL file
        vecModel = readOBJ(vecModel); // Read OBJ file
        vecModel = updateModelTriangleColour(myMap, vecModel); // Map colour to vecModel
        
        // for (size_t i = 0; i < vecModel.size(); i++) {
        //     std::cout << i << std::endl;
        //     std::cout << vecModel[i].colour.name << std::endl;
        //     std::cout << "Red: " << vecModel[i].colour.red << std::endl;
        //     std::cout << "Green: " << vecModel[i].colour.green << std::endl;
        //     std::cout << "Blue: " << vecModel[i].colour.blue << std::endl;
        //     std::cout << vecModel[i] << std::endl;
        //     // for (size_t j = 0; j < vecModel[i].vertices.size(); j++) {
        //     //     std::cout << vecModel[i].vertices[j].x << std::endl;
        //     //     std::cout << vecModel[i].vertices[j].y << std::endl;
        //     //     std::cout << vecModel[i].vertices[j].z << std::endl;
        //     // }
        // }

        std::cout << "Ray trace is drawing" << std::endl;
        drawRayTrace(window, vecModel);
        std::cout << "Ray trace drawing is completed" << std::endl;
        
    }

    // --------------------
    // ray trace + orbit + lookAt
    // --------------------

    glm::vec3 getRayDirectionOrbit(glm::vec3 cameraPosition, glm::mat3 cameraOrientation, size_t y, size_t x) {
        // Reverse implementation function of getCanvasIntersectionPoint(glm::vec3 cameraPosition, glm::vec3 vertexPosition, float focalLength) function
        // But... it iterates all pixels of canvas, which means it loops through all pixels of canvas...
        // In this function, we need to find out the actual position of vertex, instead of canvas point
        CanvasPoint point;
        point.y = y;
        point.x = x;
        point.depth = 1.0f; // Assign the depth to 1: Buffer [0] to [1] (normalise) (top)
        glm::vec3 vertexPosition;

        cameraPosition = cameraPosition * cameraOrientation; // Camera moves towards left...

        // x: positive, y: negative (inverse implementation)
        vertexPosition.x = (((point.x - HALFWIDTH) / (40 * focalLength)) * point.depth) - (cameraPosition.x); // Reverse formula
        vertexPosition.y = (- ((point.y - HALFHEIGHT) / (40 * focalLength)) * point.depth) - (cameraPosition.y); // Reverse formula
        vertexPosition.z = (0.0f * point.depth) - (cameraPosition.z); // Initial value of z: 0.0f, reverse formula
        vertexPosition = glm::normalize(vertexPosition); // Normalise to stabilise the range between 0 and 1

        // vertexPosition = vertexPosition * cameraOrientation;

        return vertexPosition;
    }
    RayTriangleIntersection getClosestValidIntersectionOrbit(glm::vec3 cameraPosition, glm::vec3 rayDirection, glm::mat3 cameraOrientation, ModelTriangle triangle) {
        RayTriangleIntersection rayTriangle;
        
        // Direction: to - from
        // [t]          ([-d e0 e1])   [s - p0] (x)
        // [u] = inverse([-d e0 e1]) * [s - p0] (y)
        // [v]          ([-d e0 e1])   [s - p0] (z)
        glm::vec3 e0 = triangle.vertices[1] - triangle.vertices[0]; // direction e0
        glm::vec3 e1 = triangle.vertices[2] - triangle.vertices[0]; // direction e1
        glm::vec3 SPVector = cameraPosition - triangle.vertices[0];
        glm::mat3 DEMatrix(-rayDirection, e0, e1);
        glm::vec3 possibleSolution = glm::inverse(DEMatrix) * SPVector;

        // (t, u, v) = (x, y, z)
        float t = possibleSolution.x; // the absolute distance along the ray from the camera to the intersection point
        float u = possibleSolution.y; // the proportional distance along the triangle's first edge that the intersection point occurs
        float v = possibleSolution.z; // the proportional distance along the triangle's second edge that the intersection point occurs

        // position = startpoint + scalar * direction
        // rayTriangle.intersectionPoint = cameraPosition + t * rayDirection;

        // r = p0 + u(p1-p0) + v(p2-p0)
        // position = startpoint + scalar * direction + scalar * direction...
        // intersection point = startpoint + first edge distance * first direction + second edge distance * second direction
        rayTriangle.intersectionPoint = triangle.vertices[0] + u * e0 + v * e1;
        rayTriangle.distanceFromCamera = t;
        rayTriangle.intersectedTriangle = triangle;
        // rayTriangle.triangleIndex = 0; // rayTriangle.triangleIndex = -1; // Invalid (ray is not intersected with triangle)
        rayTriangle.triangleIndex = std::numeric_limits<size_t>::infinity(); // positive infinity
        // (u >= 0.0) && (u <= 1.0)
        // (v >= 0.0) && (v <= 1.0)
        // (u + v) <= 1.0
        if (((u >= 0.0) && (u <= 1.0)) && ((v >= 0.0) && (v <= 1.0)) && ((u + v) <= 1.0)) { // validation process
            if (t > 0) { // the absolute distance along the ray from the camera to the intersection point should be positive
                rayTriangle.triangleIndex = 1; // It indicates that ray is intersected with triangle
            }
        }
        return rayTriangle;
    }
    RayTriangleIntersection verifyClosestValidIntersectionOrbit(glm::vec3 cameraPosition, glm::vec3 rayDirection, glm::mat3 cameraOrientation, std::vector<ModelTriangle> vecModel) {
        RayTriangleIntersection validTriangle;
        validTriangle.distanceFromCamera = std::numeric_limits<float>::infinity(); // positive infinity
        RayTriangleIntersection rayTriangle;
        for (size_t i = 0; i < vecModel.size(); i++) {
            rayTriangle = getClosestValidIntersectionOrbit(cameraPosition, rayDirection, cameraOrientation, vecModel[i]);
            if ((validTriangle.distanceFromCamera > rayTriangle.distanceFromCamera) && (rayTriangle.triangleIndex == 1)) { // rayTriangle.triangleIndex == 1 indicates valid intersected triangle
                validTriangle = rayTriangle;
                validTriangle.triangleIndex = i; // Remove this line will remove shadow
            }
        }
        return validTriangle;
    }
    void drawRayTraceOrbit(DrawingWindow& window, glm::vec3 cameraPosition, glm::mat3 cameraOrientation, std::vector<ModelTriangle> vecModel) {
        float red = 0.0f;
        float green = 0.0f;
        float blue = 0.0f;
        uint32_t black = (255 << 24) + (int(red) << 16) + (int(green) << 8) + int(blue);
        for (size_t y = 0; y < window.height; y++) { // HEIGHT
            for (size_t x = 0; x < window.width; x++) { // WIDTH
                glm::vec3 rayDirection = getRayDirectionOrbit(cameraPosition, cameraOrientation, y, x);
                RayTriangleIntersection rayTriangle = verifyClosestValidIntersectionOrbit(cameraPosition, rayDirection, cameraOrientation, vecModel);
                glm::vec3 lightDirection = glm::normalize(rayTriangle.intersectionPoint - lightPosition);
                RayTriangleIntersection lightTriangle = verifyClosestValidIntersectionOrbit(lightPosition, lightDirection, cameraOrientation, vecModel);
                // The match of index of triangles means that the light will reflect the colour of triangles
                if (rayTriangle.triangleIndex == lightTriangle.triangleIndex) {
                    uint32_t colour = (255 << 24) + (int(rayTriangle.intersectedTriangle.colour.red) << 16) + (int(rayTriangle.intersectedTriangle.colour.green) << 8) + int(rayTriangle.intersectedTriangle.colour.blue);
                    size_t stdX = std::round(x);
                    size_t stdY = std::round(y);
                    if ((stdX >= 0) && (stdX < WIDTH) && (stdY >= 0) && (stdY < HEIGHT)) { // Fix the range of window
                        window.setPixelColour(stdX, stdY, colour);
                    }
                } else { // The mismatch of index of triangles means that the light will not reflect the colour of light (shadow/darkness)
                    // The scene behind the model is blocked by the model, so light cannot reach the scene and reflect the light to camera / eyes 
                    size_t stdX = std::round(x);
                    size_t stdY = std::round(y);
                    if ((stdX >= 0) && (stdX < WIDTH) && (stdY >= 0) && (stdY < HEIGHT)) { // Fix the range of window
                        window.setPixelColour(stdX, stdY, black);
                    }
                }
            }
        }
    }
    void drawRasterisedSceneOrbit(DrawingWindow& window) {
        std::unordered_map<std::string, std::vector<float>> myMap;
        std::vector<ModelTriangle> vecModel; // 32
        myMap = readMTL(myMap); // Read MTL file
        vecModel = readOBJ(vecModel); // Read OBJ file
        vecModel = updateModelTriangleColour(myMap, vecModel); // Map colour to vecModel
        std::cout << "Ray trace is rendering!" << std::endl;
        drawRayTraceOrbit(window, cameraPosition, cameraOrientation, vecModel);
        std::cout << "Ray trace rendering is completed!" << std::endl;
    }
    
};

#endif
