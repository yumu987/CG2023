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

#ifndef MYSDL_HPP
#define MYSDL_HPP

class MySDL {
public:
    CanvasTriangle generateThreeRandomVertices(CanvasTriangle triangle) {
	    triangle.v0().x = rand() % WIDTH;
	    triangle.v0().y = rand() % HEIGHT;
        triangle.v1().x = rand() % WIDTH;
        triangle.v1().y = rand() % HEIGHT;
        triangle.v2().x = rand() % WIDTH;
        triangle.v2().y = rand() % HEIGHT;
        return triangle;
    }
    Colour generateRandomColour(Colour c) {
        c.red = rand() % 256;
        c.green = rand() % 256;
        c.blue = rand() % 256;
        return c;
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
    CanvasPoint interpolateLine(CanvasPoint from, CanvasPoint to, CanvasPoint extra) {
        // Line equation: y = mx + b
        float xDiff = to.x - from.x;
        float yDiff = to.y - from.y;
        float m = yDiff / xDiff;
        float b = to.y - (m * to.x);
        if (xDiff == 0) { // Slope is 0
            extra.x = to.x; // extra.x = from.x
        }
        else if (yDiff == 0) { // Slope is 0
            extra.x = to.x;
        }
        else if (xDiff == 0 && yDiff == 0) { // Slope is both 0
            extra.x = to.x;
        }
        else { // xDiff != 0 && yDiff != 0
            extra.x = (extra.y - b) / m;
        }
        return extra;
    }
    void drawWhiteEdgeTriangle(DrawingWindow& window, CanvasTriangle triangle, Colour c) {
        c.red = 255;
        c.green = 255;
        c.blue = 255;
        drawLine(window, triangle.v0(), triangle.v1(), c);
        drawLine(window, triangle.v1(), triangle.v2(), c);
        drawLine(window, triangle.v2(), triangle.v0(), c);
    }
    void strokedTriangle(DrawingWindow& window, CanvasTriangle triangle, Colour c) {
        triangle = generateThreeRandomVertices(triangle);
        c = generateRandomColour(c);
        drawLine(window, triangle.v0(), triangle.v1(), c);
        drawLine(window, triangle.v1(), triangle.v2(), c);
        drawLine(window, triangle.v2(), triangle.v0(), c);
    }
    void fillTopTriangle(DrawingWindow& window, CanvasPoint top, CanvasPoint middle, CanvasPoint extra, Colour c) {
        float xDiffF = middle.x - top.x;
        float yDiffF = middle.y - top.y;
        float mF = yDiffF / xDiffF;
        float bF = middle.y - (mF * middle.x);
        float xDiffS = extra.x - top.x;
        float yDiffS = extra.y - top.y;
        float mS = yDiffS / xDiffS;
        float bS = extra.y - (mS * extra.x);
        // Fill triangle
        if (xDiffF == 0 && yDiffF == 0) { // middle = top
        }
        else if (xDiffS == 0 && yDiffS == 0) { // extra = top
        }
        else {
            if (xDiffF == 0) { // middle.x = top.x
                CanvasPoint firstPoint;
                CanvasPoint secondPoint;
                firstPoint.x = top.x;
                for (float i = top.y; i < middle.y; i++) {
                    firstPoint.y = i;
                    secondPoint.y = i;
                    secondPoint.x = (secondPoint.y - bS) / mS;
                    drawLine(window, firstPoint, secondPoint, c);
                }
            }
            else if (xDiffS == 0) { // extra.x = top.x
                CanvasPoint firstPoint;
                CanvasPoint secondPoint;
                secondPoint.x = top.x;
                for (float i = top.y; i < middle.y; i++) {
                    firstPoint.y = i;
                    secondPoint.y = i;
                    firstPoint.x = (firstPoint.y - bF) / mF;
                    drawLine(window, firstPoint, secondPoint, c);
                }
            }
            else if (yDiffF == 0) { // middle.y = top.y
                CanvasPoint firstPoint;
                CanvasPoint secondPoint;
                firstPoint.y = top.y;
                for (float i = top.x; i < middle.x; i++) {
                    firstPoint.x = i;
                    secondPoint.x = i;
                    secondPoint.y = (mS * secondPoint.x) + bS;
                    drawLine(window, firstPoint, secondPoint, c);
                }
            }
            else if (yDiffS == 0) { // extra.y = top.y
                CanvasPoint firstPoint;
                CanvasPoint secondPoint;
                secondPoint.y = top.y;
                for (float i = top.x; i < middle.x; i++) {
                    firstPoint.x = i;
                    secondPoint.x = i;
                    firstPoint.y = (mF * firstPoint.x) + bF;
                    drawLine(window, firstPoint, secondPoint, c);
                }
            }
            else {
                CanvasPoint firstPoint;
                CanvasPoint secondPoint;
                for (float i = top.y; i < middle.y; i++) {
                    firstPoint.y = i;
                    secondPoint.y = i;
                    firstPoint.x = (firstPoint.y - bF) / mF;
                    secondPoint.x = (secondPoint.y - bS) / mS;
                    drawLine(window, firstPoint, secondPoint, c);
                }
            }
        }
    }
    void fillBottomTriangle(DrawingWindow& window, CanvasPoint bottom, CanvasPoint middle, CanvasPoint extra, Colour c) {
        float xDiffF = bottom.x - middle.x;
        float yDiffF = bottom.y - middle.y;
        float mF = yDiffF / xDiffF;
        float bF = middle.y - (mF * middle.x);
        float xDiffS = bottom.x - extra.x;
        float yDiffS = bottom.y - extra.y;
        float mS = yDiffS / xDiffS;
        float bS = extra.y - (mS * extra.x);
        // Fill triangle
        if (xDiffF == 0 && yDiffF == 0) { // bottom = middle
        }
        else if (xDiffS == 0 && yDiffS == 0) { // bottom = extra
        }
        else {
            if (xDiffF == 0) { // bottom.x = middle.x
                CanvasPoint firstPoint;
                CanvasPoint secondPoint;
                firstPoint.x = bottom.x;
                for (float i = middle.y; i < bottom.y; i++) {
                    firstPoint.y = i;
                    secondPoint.y = i;
                    secondPoint.x = (secondPoint.y - bS) / mS;
                    drawLine(window, firstPoint, secondPoint, c);
                }
            }
            else if (xDiffS == 0) { // bottom.x = extra.x
                CanvasPoint firstPoint;
                CanvasPoint secondPoint;
                secondPoint.x = bottom.x;
                for (float i = middle.y; i < bottom.y; i++) {
                    firstPoint.y = i;
                    secondPoint.y = i;
                    firstPoint.x = (firstPoint.y - bF) / mF;
                    drawLine(window, firstPoint, secondPoint, c);
                }
            }
            else if (yDiffF == 0) { // bottom.y = middle.y
                CanvasPoint firstPoint;
                CanvasPoint secondPoint;
                firstPoint.y = bottom.y;
                for (float i = middle.x; i < bottom.x; i++) {
                    firstPoint.x = i;
                    secondPoint.x = i;
                    secondPoint.y = (mS * secondPoint.x) + bS;
                    drawLine(window, firstPoint, secondPoint, c);
                }
            }
            else if (yDiffS == 0) { // bottom.y = extra.y
                CanvasPoint firstPoint;
                CanvasPoint secondPoint;
                secondPoint.y = bottom.y;
                for (float i = middle.x; i < bottom.x; i++) {
                    firstPoint.x = i;
                    secondPoint.x = i;
                    firstPoint.y = (mF * firstPoint.x) + bF;
                    drawLine(window, firstPoint, secondPoint, c);
                }
            }
            else {
                CanvasPoint firstPoint;
                CanvasPoint secondPoint;
                for (float i = middle.y; i < bottom.y; i++) {
                    firstPoint.y = i;
                    secondPoint.y = i;
                    firstPoint.x = (firstPoint.y - bF) / mF;
                    secondPoint.x = (secondPoint.y - bS) / mS;
                    drawLine(window, firstPoint, secondPoint, c);
                }
            }
        }
    }
    void filledTriangle(DrawingWindow& window, CanvasTriangle triangle, Colour c) {
        triangle = generateThreeRandomVertices(triangle);
        c = generateRandomColour(c);
        CanvasPoint topPoint; // Close to the top edge
        CanvasPoint middlePoint; // Close to the middle edge
        CanvasPoint bottomPoint; // Close to the bottom edge
        CanvasPoint extraPoint; // It has the same y coordinate with middlePoint
        // Sort vertices by vertical position (from top to bottom)
        // If y != y && x != x...
        if (triangle.v0().y > triangle.v1().y && triangle.v0().y > triangle.v2().y) {
            bottomPoint = triangle.v0();
            if (triangle.v1().y > triangle.v2().y) {
                middlePoint = triangle.v1();
                topPoint = triangle.v2();
            }
            else { // triangle.v2().y > triangle.v1().y
                middlePoint = triangle.v2();
                topPoint = triangle.v1();
            }
        }
        else if (triangle.v1().y > triangle.v0().y && triangle.v1().y > triangle.v2().y) {
            bottomPoint = triangle.v1();
            if (triangle.v0().y > triangle.v2().y) {
                middlePoint = triangle.v0();
                topPoint = triangle.v2();
            }
            else { // triangle.v2().y > triangle.v0().y
                middlePoint = triangle.v2();
                topPoint = triangle.v0();
            }
        }
        else if (triangle.v2().y > triangle.v0().y && triangle.v2().y > triangle.v1().y) {
            bottomPoint = triangle.v2();
            if (triangle.v0().y > triangle.v1().y) {
                middlePoint = triangle.v0();
                topPoint = triangle.v1();
            }
            else { // triangle.v1().y > triangle.v0().y
                middlePoint = triangle.v1();
                topPoint = triangle.v0();
            }
        }
        // Divide triangle into 2 "flat-bottomed" triangles
        extraPoint.y = middlePoint.y;
        extraPoint = interpolateLine(topPoint, bottomPoint, extraPoint);
        // If y == y && x == x...
        if (triangle.v0().y == triangle.v1().y && triangle.v0().x == triangle.v2().x) {
            topPoint = triangle.v2();
            middlePoint = triangle.v1();
            bottomPoint = triangle.v0();
            extraPoint = bottomPoint;
        }
        else if (triangle.v1().y == triangle.v2().y && triangle.v1().x == triangle.v0().x) {
            topPoint = triangle.v0();
            middlePoint = triangle.v2();
            bottomPoint = triangle.v1();
            extraPoint = bottomPoint;
        }
        else if (triangle.v2().y == triangle.v1().y && triangle.v2().x == triangle.v0().x) {
            topPoint = triangle.v0();
            middlePoint = triangle.v1();
            bottomPoint = triangle.v2();
            extraPoint = bottomPoint;
        }
        else if (triangle.v0().y == triangle.v2().y && triangle.v0().x == triangle.v1().x) {
            topPoint = triangle.v1();
            middlePoint = triangle.v2();
            bottomPoint = triangle.v0();
            extraPoint = bottomPoint;
        }
        else if (triangle.v1().y == triangle.v0().y && triangle.v1().x == triangle.v2().x) {
            topPoint = triangle.v2();
            middlePoint = triangle.v0();
            bottomPoint = triangle.v1();
            extraPoint = bottomPoint;
        }
        else if (triangle.v2().y == triangle.v0().y && triangle.v2().x == triangle.v1().x) {
            topPoint = triangle.v1();
            middlePoint = triangle.v0();
            bottomPoint = triangle.v2();
            extraPoint = bottomPoint;
        }
        else {
            // If y == y && x != x...
            if (triangle.v0().y == triangle.v1().y) {
                topPoint = triangle.v2();
                middlePoint = triangle.v0();
                bottomPoint = triangle.v1();
                extraPoint = bottomPoint;
            }
            else if (triangle.v1().y == triangle.v2().y) {
                topPoint = triangle.v0();
                middlePoint = triangle.v1();
                bottomPoint = triangle.v2();
                extraPoint = bottomPoint;
            }
            else if (triangle.v2().y == triangle.v0().y) {
                topPoint = triangle.v1();
                middlePoint = triangle.v2();
                bottomPoint = triangle.v0();
                extraPoint = bottomPoint;
            }
            // If y != y && x == x...
            // They are handled by previous sort vertices methods
        }
        drawLine(window, middlePoint, extraPoint, c);
        fillTopTriangle(window, topPoint, middlePoint, extraPoint, c);
        fillBottomTriangle(window, bottomPoint, middlePoint, extraPoint, c);
        drawWhiteEdgeTriangle(window, triangle, c);
    }
};

#endif
