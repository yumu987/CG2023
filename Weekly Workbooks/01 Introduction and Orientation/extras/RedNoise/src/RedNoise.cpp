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

#define WIDTH 320
#define HEIGHT 240

/*
Diagram: WIDTH * HEIGHT (x * y)
[(0, 0)   (320, 0)  ]
[                   ]
[(0, 240) (320, 240)]
*/

std::vector<float> interpolateSingleFloats(float from, float to, int numberOfValues) {
	std::vector<float> v;
	float tmp = to - from; // 'tmp' could be positive or negative
	int num = numberOfValues - 1;
	int numtmp = numberOfValues - 2; // 'from' and 'to' have been pushed back to v
	float arithmetic = tmp / num;
	v.push_back(from);
	for (int i = 0; i < numtmp; i++) {
		from = from + arithmetic; // 'from' + positive / negative number
		v.push_back(from);
	}
	v.push_back(to);
	return v;
}

std::vector<glm::vec3> interpolateThreeElementValues(glm::vec3 from, glm::vec3 to, int numberOfValues) {
	std::vector<glm::vec3> v;
	std::vector<float> resultX;
	std::vector<float> resultY;
	std::vector<float> resultZ;
	resultX = interpolateSingleFloats(from.x, to.x, numberOfValues);
	resultY = interpolateSingleFloats(from.y, to.y, numberOfValues);
	resultZ = interpolateSingleFloats(from.z, to.z, numberOfValues);
	for (int i = 0; i < numberOfValues; i++) {
		glm::vec3 tmpA;
		tmpA.x = resultX[i];
		tmpA.y = resultY[i];
		tmpA.z = resultZ[i];
		v.push_back(tmpA);
	}
	return v;
}

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

/*
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

std::vector<CanvasPoint> interpolationCanvasPoint(CanvasPoint from, CanvasPoint to, int numberOfSteps) {
	std::vector<CanvasPoint> coordinates;

	float xDiff = to.x - from.x;
	float yDiff = to.y - from.y;
	float xStepSize = xDiff / numberOfSteps;
	float yStepSize = yDiff / numberOfSteps;
	
	float xDiffT = to.texturePoint.x - from.texturePoint.x;
	float yDiffT = to.texturePoint.y - from.texturePoint.y;
	float xStepSizeT = xDiffT / numberOfSteps;
	float yStepSizeT = yDiffT / numberOfSteps;

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
	std::vector<CanvasPoint> canvasCoordinates;
	canvasCoordinates = interpolationCanvasPoint(from, to, std::round(to.x - from.x));

	// for (const CanvasPoint& canPoint : canvasCoordinates) {
	// 	window.setPixelColour(std::round(canPoint.x), std::round(canPoint.y), textMap.pixels[(canPoint.texturePoint.y * textMap.height) + canPoint.texturePoint.x]);
	// }

	for (size_t i = 0; i < canvasCoordinates.size(); i++) {
		window.setPixelColour(canvasCoordinates[i].x, canvasCoordinates[i].y, textMap.pixels[(canvasCoordinates[i].texturePoint.y * textMap.height) + canvasCoordinates[i].texturePoint.x]); // textMap.width
	}
}

CanvasPoint interpolateLine(CanvasPoint from, CanvasPoint to, CanvasPoint extra) {
	// Line equation: y = mx + b
	float xDiff = to.x - from.x;
	float yDiff = to.y - from.y;
	float m = yDiff / xDiff;
	float b = to.y - (m * to.x);
	if (xDiff == 0) { // Slope is 0
		std::cout << "xDiff == 0 interpolateLine" << std::endl;
		extra.x = to.x; // extra.x = from.x
	}
	else if (yDiff == 0) { // Slope is 0
		std::cout << "yDiff == 0 interpolateLine" << std::endl;
		extra.x = to.x;
	}
	else if (xDiff == 0 && yDiff == 0) { // Slope is both 0
		std::cout << "xDiff == 0 && yDiff == 0 interpolateLine" << std::endl;
		extra.x = to.x;
	}
	else { // xDiff != 0 && yDiff != 0
		extra.x = (extra.y - b) / m;
	}
	return extra;
}

void strokedTriangle(DrawingWindow& window, CanvasTriangle triangle, Colour c) {
	c = generateRandomColour(c);
	drawLine(window, triangle.v0(), triangle.v1(), c);
	drawLine(window, triangle.v1(), triangle.v2(), c);
	drawLine(window, triangle.v2(), triangle.v0(), c);
}

void drawWhiteEdgeTriangle(DrawingWindow& window, CanvasTriangle triangle, Colour c) {
	c.red = 255;
	c.green = 255;
	c.blue = 255;
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
		std::cout << "xDiffF == 0 && yDiffF == 0 fillTopTriangle" << std::endl;
		std::cout << "pass" << std::endl;
	}
	else if (xDiffS == 0 && yDiffS == 0) { // extra = top
		std::cout << "xDiffS == 0 && yDiffS == 0 fillTopTriangle" << std::endl;
		std::cout << "pass" << std::endl;
	}
	else {
		if (xDiffF == 0) { // middle.x = top.x
			std::cout << "xDiffF == 0 fillTopTriangle" << std::endl;
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
			std::cout << "xDiffS == 0 fillTopTriangle" << std::endl;
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
			std::cout << "yDiffF == 0 fillTopTriangle" << std::endl;
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
			std::cout << "yDiffS == 0 fillTopTriangle" << std::endl;
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
		std::cout << "xDiffF == 0 && yDiffF == 0 fillBottomTriangle" << std::endl;
		std::cout << "pass" << std::endl;
	}
	else if (xDiffS == 0 && yDiffS == 0) { // bottom = extra
		std::cout << "xDiffS == 0 && yDiffS == 0 fillBottomTriangle" << std::endl;
		std::cout << "pass" << std::endl;
	}
	else {
		if (xDiffF == 0) { // bottom.x = middle.x
			std::cout << "xDiffF == 0 fillBottomTriangle" << std::endl;
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
			std::cout << "xDiffS == 0 fillBottomTriangle" << std::endl;
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
			std::cout << "yDiffF == 0 fillBottomTriangle" << std::endl;
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
			std::cout << "yDiffS == 0 fillBottomTriangle" << std::endl;
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
	// Random colour
	c = generateRandomColour(c);
	// Points
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
	// Draw a line to divide triangle
	drawLine(window, middlePoint, extraPoint, c);
	// Fill triangle
	fillTopTriangle(window, topPoint, middlePoint, extraPoint, c);
	fillBottomTriangle(window, bottomPoint, middlePoint, extraPoint, c);
	// Draw white edge
	drawWhiteEdgeTriangle(window, triangle, c);
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
		std::cout << "xDiff == 0 interpolateTextureLine" << std::endl;
		extra.x = to.x; // extra.x = from.x
	}
	else if (yDiff == 0) { // Slope is 0
		std::cout << "yDiff == 0 interpolateTextureLine" << std::endl;
		extra.x = to.x;
	}
	else if (xDiff == 0 && yDiff == 0) { // Slope is both 0
		std::cout << "xDiff == 0 && yDiff == 0 interpolateTextureLine" << std::endl;
		extra.x = to.x;
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
	firstLine = interpolationCanvasPoint(middle, bottom, bottom.y - middle.y);

	std::vector<CanvasPoint> secondLine;
	secondLine = interpolationCanvasPoint(extra, bottom, bottom.y - extra.y);

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

void draw(DrawingWindow& window) {
	window.clearPixels();
	for (size_t y = 0; y < window.height; y++) {
		for (size_t x = 0; x < window.width; x++) {
			// Test Black Background
			float red = 0.0;
			float green = 0.0;
			float blue = 0.0;
			uint32_t colour = (255 << 24) + (int(red) << 16) + (int(green) << 8) + int(blue);
			window.setPixelColour(x, y, colour);
		}
	}
}

void handleEvent(SDL_Event event, DrawingWindow& window) {
	if (event.type == SDL_KEYDOWN) {
		if (event.key.keysym.sym == SDLK_LEFT) std::cout << "LEFT" << std::endl;
		else if (event.key.keysym.sym == SDLK_RIGHT) std::cout << "RIGHT" << std::endl;
		else if (event.key.keysym.sym == SDLK_UP) std::cout << "UP" << std::endl;
		else if (event.key.keysym.sym == SDLK_DOWN) std::cout << "DOWN" << std::endl;
		else if (event.key.keysym.sym == SDLK_u) { // u: Draw random stroked (unfilled) triangle
			std::cout << "u" << "/" << "U" << ":" << "Stroked triangle" << std::endl;
			CanvasTriangle triangle;
			Colour c;
			triangle = generateThreeRandomVertices(triangle);
			strokedTriangle(window, triangle, c);
		}
		else if (event.key.keysym.sym == SDLK_f) { // f: Draw random filled triangle
			std::cout << "f" << "/" << "F" << ":" << "Filled triangle" << std::endl;
			CanvasTriangle triangle;
			Colour c;
			triangle = generateThreeRandomVertices(triangle);
			filledTriangle(window, triangle, c);
		}
	}
	else if (event.type == SDL_MOUSEBUTTONDOWN) {
		window.savePPM("output.ppm");
		window.saveBMP("output.bmp");
	}
}

int main(int argc, char* argv[]) {
	DrawingWindow window = DrawingWindow(WIDTH, HEIGHT, false);
	SDL_Event event;
	/*----------*/
	CanvasTriangle triangle; // Debug triangle for white edge
	Colour c; // Debug colour
	CanvasTriangle textureTriangle;
	TextureMap textMap = TextureMap("texture.ppm"); // Load texture.ppm
	/*----------*/
	draw(window); // Fixed black background
	/*----------*/
	fillTexture(window, textureTriangle, textMap);
	/*----------*/
	// Debug white edge below to ensure texture is correct
	triangle = generateVisualVerificationVertices(triangle);
	// triangle = bubbleSort(triangle);
	drawWhiteEdgeTriangle(window, triangle, c);
	/*----------*/
	while (true) {
		// We MUST poll for events - otherwise the window will freeze !
		if (window.pollForInputEvents(event)) handleEvent(event, window);
		// Need to render the frame at the end, or nothing actually gets shown on the screen !
		window.renderFrame();
	}
}
