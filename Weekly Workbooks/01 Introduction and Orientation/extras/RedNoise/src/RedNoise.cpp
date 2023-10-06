#include <CanvasTriangle.h>
#include <DrawingWindow.h>
#include <Utils.h>
#include <fstream>
#include <vector>
#include <glm/glm.hpp>
#include <CanvasPoint.h>
#include <Colour.h>

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
	/*--------------------*/
	// std::vector<glm::vec3> v;
	// int num = numberOfValues - 1;
	// int numtmp = numberOfValues - 2;
	// v.push_back(from);
	// for(int i = 0; i < numtmp; i++) {
	// 	if(from.x > to.x) {
	// 		float tmpFirst = from.x - to.x;
	// 		float decreFirst = tmpFirst / num;
	// 		from.x = from.x - decreFirst;
	// 	} else if(from.x < to.x) {
	// 		float tmpSecond = to.x - from.x;
	// 		float increSecond = tmpSecond / num;
	// 		from.x = from.x + increSecond;
	// 	} else {
	// 		from.x = from.x;
	// 	}
	// 	if(from.y > to.y) {
	// 		float tmpThird = from.y - to.y;
	// 		float decreThird = tmpThird / num;
	// 		from.y = from.y - decreThird;
	// 	} else if(from.y < to.y) {
	// 		float tmpFourth = to.y - from.y;
	// 		float increFourth = tmpFourth / num;
	// 		from.y = from.y + increFourth;
	// 	} else {
	// 		from.y = from.y;
	// 	}
	// 	if(from.z > to.z) {
	// 		float tmpFifth = from.z - to.z;
	// 		float decreFifth = tmpFifth / num;
	// 		from.z = from.z - decreFifth;
	// 	} else if(from.z < to.z) {
	// 		float tmpSixth = to.z - from.z;
	// 		float increSixth = tmpSixth / num;
	// 		from.z = from.z + increSixth;
	// 	} else {
	// 		from.z = from.z;
	// 	}
	// 	v.push_back(from);
	// }
	// v.push_back(to);
	// return v;
	/*--------------------*/
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
	// Random three vertices (points)
	triangle.v0().x = rand() % WIDTH;
	triangle.v0().y = rand() % HEIGHT;
	triangle.v1().x = rand() % WIDTH;
	triangle.v1().y = rand() % HEIGHT;
	triangle.v2().x = rand() % WIDTH;
	triangle.v2().y = rand() % HEIGHT;
	return triangle;
}

Colour generateRandomColour(Colour c) {
	// Random colour
	c.red = rand() % 256;
	c.green = rand() % 256;
	c.blue = rand() % 256;
	return c;
}

void drawLine(DrawingWindow &window, CanvasPoint from, CanvasPoint to, Colour c) {
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

CanvasPoint flatLine(CanvasPoint from, CanvasPoint to, CanvasPoint extra) {
	// Line equation: y = mx + b
	// m = yDiff / xDiff // m = (y2 - y1) / (x2 - x1)
	// b = y - mx // Substitute one exact point (x1, y1) OR (x2, y2) to get constant b
	// x = (y - b) / m // Use exact y value of point to get x
	float xDiff = to.x - from.x;
	float yDiff = to.y - from.y;
	float m = yDiff / xDiff;
	float b = to.y - (m * to.x); // Substitute one exact point to get constant b
	extra.x = (extra.y - b) / m; // Substitute extraPoint's y
	return extra;
}

void strokedTriangle(DrawingWindow &window, CanvasTriangle triangle, Colour c) {
	// Random colour
	c = generateRandomColour(c);
	// Draw edge
	drawLine(window, triangle.v0(), triangle.v1(), c);
	drawLine(window, triangle.v1(), triangle.v2(), c);
	drawLine(window, triangle.v2(), triangle.v0(), c);
}

void fillTopTriangle() {
	// Fill top triangle by drawing lines
}

void fillBottomTriangle() {
	// Fill bottom triangle by drawing lines
}

void filledTriangle(DrawingWindow &window, CanvasTriangle triangle, Colour c) {
	// Random colour
	c = generateRandomColour(c);
	// Draw edge
	drawLine(window, triangle.v0(), triangle.v1(), c);
	drawLine(window, triangle.v1(), triangle.v2(), c);
	drawLine(window, triangle.v2(), triangle.v0(), c);
	CanvasPoint topPoint; // Close to the top edge
	CanvasPoint middlePoint; // Close to the middle edge
	CanvasPoint bottomPoint; // Close to the bottom edge
	CanvasPoint extraPoint; // Extra point has the same y-coordinate with middle point
	// Sort vertices by vertical position (top to bottom)
	if (triangle.v0().y > triangle.v1().y && triangle.v0().y > triangle.v2().y) {
		bottomPoint = triangle.v0();
		if (triangle.v1().y > triangle.v2().y) {
			middlePoint = triangle.v1();
			topPoint = triangle.v2();
		} else { // triangle.v2().y > triangle.v1().y
			middlePoint = triangle.v2();
			topPoint = triangle.v1();
		}
	} else if (triangle.v1().y > triangle.v0().y && triangle.v1().y > triangle.v2().y) {
		bottomPoint = triangle.v1();
		if (triangle.v0().y > triangle.v2().y) {
			middlePoint = triangle.v0();
			topPoint = triangle.v2();
		} else { // triangle.v2().y > triangle.v0().y
			middlePoint = triangle.v2();
			topPoint = triangle.v0();
		}
	} else if (triangle.v2().y > triangle.v0().y && triangle.v2().y > triangle.v1().y) {
		bottomPoint = triangle.v2();
		if (triangle.v0().y > triangle.v1().y) {
			middlePoint = triangle.v0();
			topPoint = triangle.v1();
		} else { // triangle.v1().y > triangle.v0().y
			middlePoint = triangle.v1();
			topPoint = triangle.v0();
		}
	}
	// Divide triangle into 2 "flat-bottomed" triangles
	extraPoint.y = middlePoint.y;
	extraPoint = flatLine(topPoint, bottomPoint, extraPoint);
	// Draw a line to divide triangle
	drawLine(window, middlePoint, extraPoint, c);
	// Fill top triangle
	// ......
	// Fill bottom triangle
	// ......
}

void draw(DrawingWindow &window) {
	window.clearPixels();
	/*--------------------*/
	// Black & White
	// float startColour = 255; // White colour
	// float endColour = 0; // Black colour
	/*--------------------*/
	// Colour
	// glm::vec3 topLeft(255, 0, 0);        // red 
	// glm::vec3 topRight(0, 0, 255);       // blue 
	// glm::vec3 bottomRight(0, 255, 0);    // green 
	// glm::vec3 bottomLeft(255, 255, 0);   // yellow
	// std::vector<glm::vec3> rightColumn;
	// std::vector<glm::vec3> leftColumn;
	// rightColumn = interpolateThreeElementValues(topLeft, bottomLeft, HEIGHT);
	// leftColumn = interpolateThreeElementValues(topRight, bottomRight, HEIGHT);
	/*--------------------*/
	for (size_t y = 0; y < window.height; y++) {
		/*--------------------*/
		// Colour
		// std::vector<glm::vec3> rowVector;
		// rowVector = interpolateThreeElementValues(rightColumn[y], leftColumn[y], WIDTH);
		/*--------------------*/
		for (size_t x = 0; x < window.width; x++) {
			/*--------------------*/
			// RedNoise
			// float red = rand() % 256;
			// float green = 0.0;
			// float blue = 0.0;
			// uint32_t colour = (255 << 24) + (int(red) << 16) + (int(green) << 8) + int(blue);
			// window.setPixelColour(x, y, colour);
			/*--------------------*/
			// Greyscale
			// float red = (startColour + (endColour - startColour) * x / WIDTH);
			// float green = (startColour + (endColour - startColour) * x / WIDTH);
			// float blue = (startColour + (endColour - startColour) * x / WIDTH);
			// uint32_t colour = (255 << 24) + (int(red) << 16) + (int(green) << 8) + int(blue);
			// window.setPixelColour(x, y, colour);
			/*--------------------*/
			// Greyscale
			// std::vector<float> resultFor;
			// resultFor = interpolateSingleFloats(255, 0, WIDTH); // split from white to black for a whole window
			// // iterate from strat index to end index
			// float red = resultFor[x];
			// float green = resultFor[x];
			// float blue = resultFor[x];
			// uint32_t colour = (255 << 24) + (int(red) << 16) + (int(green) << 8) + int(blue);
			// window.setPixelColour(x, y, colour);
			/*--------------------*/
			// Colour
			// float red = rowVector[x].x;
			// float green = rowVector[x].y;
			// float blue = rowVector[x].z;
			// uint32_t colour = (255 << 24) + (int(red) << 16) + (int(green) << 8) + int(blue);
			// window.setPixelColour(x, y, colour);
			/*--------------------*/
			// Test Black Background
			float red = 0.0;
			float green = 0.0;
			float blue = 0.0;
			uint32_t colour = (255 << 24) + (int(red) << 16) + (int(green) << 8) + int(blue);
			window.setPixelColour(x, y, colour);
		}
	}
}

void handleEvent(SDL_Event event, DrawingWindow &window) {
	if (event.type == SDL_KEYDOWN) {
		if (event.key.keysym.sym == SDLK_LEFT) std::cout << "LEFT" << std::endl;
		else if (event.key.keysym.sym == SDLK_RIGHT) std::cout << "RIGHT" << std::endl;
		else if (event.key.keysym.sym == SDLK_UP) std::cout << "UP" << std::endl;
		else if (event.key.keysym.sym == SDLK_DOWN) std::cout << "DOWN" << std::endl;
		else if (event.key.keysym.sym == SDLK_u) { // u: Draw random stroked (unfilled) triangle
			std::cout << "u" << "/" << "U" << std::endl;
			CanvasTriangle triangle;
			Colour c;
			triangle = generateThreeRandomVertices(triangle);
			strokedTriangle(window, triangle, c);
		} else if (event.key.keysym.sym == SDLK_f) { // f: Draw random filled triangle
			std::cout << "f" << "/" << "F" << std::endl;
			CanvasTriangle triangle;
			Colour c;
			triangle = generateThreeRandomVertices(triangle);
			filledTriangle(window, triangle, c);
		}
	} else if (event.type == SDL_MOUSEBUTTONDOWN) {
		window.savePPM("output.ppm");
		window.saveBMP("output.bmp");
	}
}

int main(int argc, char *argv[]) {
	DrawingWindow window = DrawingWindow(WIDTH, HEIGHT, false);
	SDL_Event event;
	draw(window); // Fixed black background
	/*--------------------*/
	// std::cout << std::endl; // Debug copy
	/*--------------------*/
	while (true) {
		// We MUST poll for events - otherwise the window will freeze !
		if (window.pollForInputEvents(event)) handleEvent(event, window);
		// Need to render the frame at the end, or nothing actually gets shown on the screen !
		window.renderFrame();
	}
}
