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

#define WIDTH 320 * 2
#define HEIGHT 240 * 2
#define HALFWIDTH 160 * 2
#define HALFHEIGHT 120 * 2
#define OBJfilename "cornell-box.obj"
#define MTLfilename "cornell-box.mtl"

float depthBuffer[WIDTH][HEIGHT];
glm::vec3 cameraPosition (0.0f, 0.0f, 4.0f);
float focalLength = 2.0f;

#ifndef MYMODEL_HPP
#define MYMODEL_HPP

class MyModel {
public:
    void resetModel() {
        cameraPosition = glm::vec3(0.0f, 0.0f, 4.0f);
        focalLength = 2.0f;
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
            window.setPixelColour(x, y, colour);
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
                if (z >= depthBuffer[stdX][stdY]) {
                    depthBuffer[stdX][stdY] = z; // Record depth buffer
                    window.setPixelColour(stdX, stdY, colour);
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
        */
        CanvasPoint point;
        point.x = 160 * focalLength * ( - (vertexPosition.x - cameraPosition.x) / (vertexPosition.z - cameraPosition.z)) + HALFWIDTH;
        point.y = 160 * focalLength * ((vertexPosition.y - cameraPosition.y) / (vertexPosition.z - cameraPosition.z)) + HALFHEIGHT;
        // If z(1/Z) is greater, it is close to camera
        // If z(1/Z) is smaller, it is far from camera
        point.depth = 1 / - (vertexPosition.z - cameraPosition.z);
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
};

#endif
