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

#ifndef MYTEXTURE_HPP
#define MYTEXTURE_HPP

class MyTexture {
public:
    void drawLine(DrawingWindow& window, CanvasPoint from, CanvasPoint to, Colour c) {
        float xDiff = to.x - from.x;
        float yDiff = to.y - from.y;
        float numberOfSteps = fmax(abs(xDiff), abs(yDiff));
        float xStepSize = xDiff / numberOfSteps;
        float yStepSize = yDiff / numberOfSteps;
        uint32_t colour = (255 << 24) + (c.red << 16) + (c.green << 8) + c.blue; // Pack colour into uint32_t package
        for (float i = 0.0; i < numberOfSteps; i++) {
            float x = from.x + (xStepSize * i);
            float y = from.y + (yStepSize * i);
            window.setPixelColour(x, y, colour);
        }
    }
    void drawWhiteEdgeTriangle(DrawingWindow& window, CanvasTriangle triangle, Colour c) {
        c.red = 255;
        c.green = 255;
        c.blue = 255;
        drawLine(window, triangle.v0(), triangle.v1(), c);
        drawLine(window, triangle.v1(), triangle.v2(), c);
        drawLine(window, triangle.v2(), triangle.v0(), c);
    }
    CanvasPoint interpolateLine(CanvasPoint from, CanvasPoint to, CanvasPoint extra) {
        // Line equation: y = mx + b
        float xDiff = to.x - from.x;
        float yDiff = to.y - from.y;
        float m = yDiff / xDiff;
        float b = to.y - (m * to.x);
        if (xDiff == 0) { // Slope is 0
            extra.x = to.x; // extra.x = from.x
        } else if (yDiff == 0) { // Slope is 0
            extra.x = to.x; // extra.x = from.x
        } else if (xDiff == 0 && yDiff == 0) { // Slope is both 0
            extra.x = to.x; // extra.x = from.x
        } else { // xDiff != 0 && yDiff != 0
            extra.x = (extra.y - b) / m;
        }
        return extra;
    }
    std::vector<CanvasPoint> interpolationCanvasPoint(CanvasPoint from, CanvasPoint to, int numberOfSteps) {
        std::vector<CanvasPoint> coordinates;

        float xDiff = to.x - from.x;
        float yDiff = to.y - from.y;
        float xStepSize = xDiff / numberOfSteps;
        float yStepSize = yDiff / numberOfSteps;
        
        // float xDiffT = to.texturePoint.x - from.texturePoint.x;
        // float yDiffT = to.texturePoint.y - from.texturePoint.y;
        // float xStepSizeT = xDiffT / numberOfSteps;
        // float yStepSizeT = yDiffT / numberOfSteps;

        for (float i = 0.0; i < numberOfSteps; i++) {
            CanvasPoint point;
            float x = from.x + (xStepSize * i);
            float y = from.y + (yStepSize * i);
            float xT = from.texturePoint.x + (xStepSize * i);
            float yT = from.texturePoint.y + (yStepSize * i);
            point.x = x;
            point.y = y;
            point.texturePoint.x = xT;
            point.texturePoint.y = yT;
            coordinates.push_back(point);
        }

        return coordinates;
    }
    void drawTextureLine(DrawingWindow& window, CanvasPoint from, CanvasPoint to, TextureMap textMap) {
        /*
        Notes for triangle's texture mapping:
        The pixel data is just stored linearly...
        Assume we got a 3x3 matrix:
        width: 3, height: 3
            0	1	2
        0	[]	[]	[]
        1	[]	[]	[]
        2	[]	[X]	[]
        We need to find X's index...
        As we mentioned before, pixel data are stored linearly...
        So, we can get:
            [0]	[1]	[2]
            [3]	[4]	[5]
            [6]	[7]	[8]
        A 3x3 matrix stores data like this.
        If we got width, height and X's coordinate: (1, 2) = (x, y), how could we get X's index?
        Use index = y * height + x. / index = y * width + x.
        Substitute into our equation:
        Index: 7 = 2 * 3 + 1 = 6 + 1 = 7
        */
        std::vector<CanvasPoint> canvasCoordinates;
        canvasCoordinates = interpolationCanvasPoint(from, to, std::round(to.x - from.x));

        // for (const CanvasPoint& canPoint : canvasCoordinates) {
        // 	window.setPixelColour(std::round(canPoint.x), std::round(canPoint.y), textMap.pixels[(canPoint.texturePoint.y * textMap.height) + canPoint.texturePoint.x]);
        // }

        for (size_t i = 0; i < canvasCoordinates.size(); i++) {
            window.setPixelColour(canvasCoordinates[i].x, canvasCoordinates[i].y, textMap.pixels[(canvasCoordinates[i].texturePoint.y * textMap.height) + canvasCoordinates[i].texturePoint.x]); // textMap.width
        }
    }
    CanvasTriangle generateVisualVerificationVertices(CanvasTriangle triangle) { // Debug function
        triangle.v0().x = 160;
        triangle.v0().y = 10;
        triangle.v1().x = 300;
        triangle.v1().y = 230;
        triangle.v2().x = 10;
        triangle.v2().y = 150;
        return triangle;
    }
    TexturePoint interpolateTextureLine(TexturePoint from, TexturePoint to, TexturePoint extra) {
        // Line equation: y = mx + b
        float xDiff = to.x - from.x;
        float yDiff = to.y - from.y;
        float m = yDiff / xDiff;
        float b = to.y - (m * to.x);
        if (xDiff == 0) { // Slope is 0
            extra.x = to.x; // extra.x = from.x
        }
        else if (yDiff == 0) { // Slope is 0
            extra.x = to.x; // extra.x = from.x
        }
        else if (xDiff == 0 && yDiff == 0) { // Slope is both 0
            extra.x = to.x; // extra.x = from.x
        }
        else { // xDiff != 0 && yDiff != 0
            extra.x = (extra.y - b) / m;
        }
        return extra;
    }
    void fillTopTextureTriangle(DrawingWindow& window, CanvasPoint top, CanvasPoint middle, CanvasPoint extra, TextureMap textMap) {
        std::vector<CanvasPoint> firstLine;
        firstLine = interpolationCanvasPoint(top, middle, std::round(middle.y - top.y));

        std::vector<CanvasPoint> secondLine;
        secondLine = interpolationCanvasPoint(top, extra, std::round(extra.y - top.y));

        for (size_t i = 0; i < firstLine.size(); i++) { // i < secondLine.size()
            drawTextureLine(window, firstLine[i], secondLine[i], textMap);
        }
    }
    void fillBottomTextureTriangle(DrawingWindow& window, CanvasPoint bottom, CanvasPoint middle, CanvasPoint extra, TextureMap textMap) {
        std::vector<CanvasPoint> firstLine;
        firstLine = interpolationCanvasPoint(middle, bottom, std::round(bottom.y - middle.y));

        std::vector<CanvasPoint> secondLine;
        secondLine = interpolationCanvasPoint(extra, bottom, std::round(bottom.y - extra.y));

        for (size_t i = 0; i < firstLine.size(); i++) { // i < secondLine.size()
            drawTextureLine(window, firstLine[i], secondLine[i], textMap);
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
    std::vector<CanvasPoint> mapPoint(CanvasPoint v0, CanvasPoint v1, CanvasPoint v2) {
        std::vector<CanvasPoint> canVector;
        // CanvasPoint initialisation
        // v0: (160, 10)
        v0.x = 160;
        v0.y = 10;
        // v1: (10, 150)
        v1.x = 10;
        v1.y = 150;
        // v2: (300, 230)
        v2.x = 300;
        v2.y = 230;
        // TexturePoint mapping
        // v0: (195, 5)
        v0.texturePoint.x = 195;
        v0.texturePoint.y = 5;
        // v1: (65, 330)
        v1.texturePoint.x = 65;
        v1.texturePoint.y = 330;
        // v2: (395, 380)
        v2.texturePoint.x = 395;
        v2.texturePoint.y = 380;
        // Push back vertices
        canVector.push_back(v0);
        canVector.push_back(v1);
        canVector.push_back(v2);
        return canVector;
    }
    void fillTexture(DrawingWindow& window, CanvasTriangle textureTriangle, TextureMap textMap) {
        // Initialisation of points and sorting algorithm
        CanvasPoint v0;
        CanvasPoint v1;
        CanvasPoint v2;
        CanvasPoint extra;
        std::vector<CanvasPoint> canVector;
        std::vector<CanvasPoint> lineExtra;
        canVector = mapPoint(v0, v1, v2);
        textureTriangle.v0() = canVector[0];
        textureTriangle.v1() = canVector[1];
        textureTriangle.v2() = canVector[2];
        textureTriangle = bubbleSort(textureTriangle); // Sort top, middle and bottom
        // Update extra point in canvas
        extra.y = textureTriangle.v1().y;
        extra = interpolateLine(textureTriangle.v0(), textureTriangle.v2(), extra);
        // Update extra point in texture
        lineExtra = interpolationCanvasPoint(textureTriangle.v0(), textureTriangle.v2(), textureTriangle.v2().y - textureTriangle.v0().y);
        for (size_t i = 0; i < lineExtra.size(); i++) {
            if (extra.x == lineExtra[i].x && extra.y == lineExtra[i].y) {
                extra.texturePoint.x = lineExtra[i].texturePoint.x;
                extra.texturePoint.y = lineExtra[i].texturePoint.y;
            }
        }
        /*----------*/
        std::cout << "textureTriangle.v0(): " << textureTriangle.v0() << std::endl;
        std::cout << "textureTriangle.v1(): " << textureTriangle.v1() << std::endl;
        std::cout << "textureTriangle.v2(): " << textureTriangle.v2() << std::endl;
        std::cout << "extra: " << extra << std::endl;
        std::cout << "textureTriangle.v0().texture: " << textureTriangle.v0().texturePoint << std::endl;
        std::cout << "textureTriangle.v1().texture: " << textureTriangle.v1().texturePoint << std::endl;
        std::cout << "textureTriangle.v2().texture: " << textureTriangle.v2().texturePoint << std::endl;
        std::cout << "extra.texturePoint: " << extra.texturePoint << std::endl;
        /*----------*/
        // Divide triangle
        drawTextureLine(window, textureTriangle.v1(), extra, textMap); // Middle point with extra
        // Draw top
        fillTopTextureTriangle(window, textureTriangle.v0(), textureTriangle.v1(), extra, textMap); // Top point, middle point and extra
        // Draw bottom
        fillBottomTextureTriangle(window, textureTriangle.v2(), textureTriangle.v1(), extra, textMap); // Bottom point, middle point and extra
    }
    void drawTexture(DrawingWindow& window) {
        CanvasTriangle triangle; // Debug triangle for white edge
        Colour c; // Debug colour
        CanvasTriangle textureTriangle;
        TextureMap textMap = TextureMap("texture.ppm"); // Load texture.ppm
        fillTexture(window, textureTriangle, textMap);
        // Debug white edge below to ensure texture is correct
        triangle = generateVisualVerificationVertices(triangle);
        triangle = bubbleSort(triangle);
        drawWhiteEdgeTriangle(window, triangle, c);
    }
};

#endif
