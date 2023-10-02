#include <CanvasTriangle.h>
#include <DrawingWindow.h>
#include <Utils.h>
#include <fstream>
#include <vector>
#include <glm/glm.hpp>

#define WIDTH 320
#define HEIGHT 240

std::vector<float> interpolateSingleFloats(float from, float to, int numberOfValues) {
	std::vector<float> v;
	float tmp = to - from; // tmp could be positive or negative
	int num = numberOfValues - 1;
	int numtmp = numberOfValues - 2; // from and to have been pushed back to v
	float arithmetic = tmp / num;
	v.push_back(from);
	for(int i = 0; i < numtmp; i++) {
		from = from + arithmetic; // from + positive/negative number
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
	for(int i = 0; i < numberOfValues; i++) {
		glm::vec3 tmpA;
		tmpA.x = resultX[i];
		tmpA.y = resultY[i];
		tmpA.z = resultZ[i];
		v.push_back(tmpA);
	} 
	return v;
}

void draw(DrawingWindow &window) {
	window.clearPixels();
	// Black & White
	// float startColour = 255; // white colour
	// float endColour = 0; // black colour
	// Colour
	glm::vec3 topLeft(255, 0, 0);        // red 
	glm::vec3 topRight(0, 0, 255);       // blue 
	glm::vec3 bottomRight(0, 255, 0);    // green 
	glm::vec3 bottomLeft(255, 255, 0);   // yellow
	std::vector<glm::vec3> rightColumn;
	std::vector<glm::vec3> leftColumn;
	rightColumn = interpolateThreeElementValues(topLeft, bottomLeft, HEIGHT);
	leftColumn = interpolateThreeElementValues(topRight, bottomRight, HEIGHT);
	for (size_t y = 0; y < window.height; y++) {
		std::vector<glm::vec3> rowVector;
		rowVector = interpolateThreeElementValues(rightColumn[y], leftColumn[y], WIDTH);
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
			float red = rowVector[x].x;
			float green = rowVector[x].y;
			float blue = rowVector[x].z;
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
	} else if (event.type == SDL_MOUSEBUTTONDOWN) {
		window.savePPM("output.ppm");
		window.saveBMP("output.bmp");
	}
}

int main(int argc, char *argv[]) {
	DrawingWindow window = DrawingWindow(WIDTH, HEIGHT, false);
	SDL_Event event;
	/*--------------------*/
	std::vector<float> result;
	result = interpolateSingleFloats(2.2, 8.5, 7);
	for(size_t i=0; i<result.size(); i++) std::cout << result[i] << " ";
	std::cout << std::endl;
	/*--------------------*/
	glm::vec3 from(1.0, 4.0, 9.2);
	glm::vec3 to(4.0, 1.0, 9.8);
	std::vector<glm::vec3> newResult;
	newResult = interpolateThreeElementValues(from, to, 4);
	for(size_t i=0; i<newResult.size(); i++) std::cout << "(" << newResult[i].x << ", " << newResult[i].y << ", " << newResult[i].z << ")" << std::endl;
	/*--------------------*/
	while (true) {
		// We MUST poll for events - otherwise the window will freeze !
		if (window.pollForInputEvents(event)) handleEvent(event, window);
		draw(window);
		// Need to render the frame at the end, or nothing actually gets shown on the screen !
		window.renderFrame();
	}
}
