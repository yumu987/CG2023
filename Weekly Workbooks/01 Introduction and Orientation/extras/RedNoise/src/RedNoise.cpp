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
#include "mysdl.hpp"
#include "mymodel.hpp"

#define WIDTH 320 * 2
#define HEIGHT 240 * 2
#define HALFWIDTH 160 * 2
#define HALFHEIGHT 120 * 2
#define OBJfilename "cornell-box.obj"
#define MTLfilename "cornell-box.mtl"

// float depthBuffer[WIDTH][HEIGHT];
// glm::vec3 cameraPosition (0.0f, 0.0f, 4.0f);
// float focalLength = 2.0f;
// glm::vec3 cameraToVertex (3.0f, 4.0f, 5.0f);
// glm::vec3 targetPoint(-0.04830505f, 0.0039152836f, -0.0014743528f); // (0.0f, 0.0f, 0.0f)
glm::vec3 targetPoint(0.0f, 0.0f, 0.0f);
glm::vec3 upVector(0.0f, 1.0f, 0.0f);

/*
Notes for drawing 2D diagram:
Coordinates System:
Diagram: WIDTH * HEIGHT / (x * y)
[(0, 0)   (320, 0)  ]
[                   ]
[(0, 240) (320, 240)]
------------------------------
[(0, 0)   		(WIDTH, 0)  ]
[                       	]
[(0, HEIGHT) (WIDTH, HEIGHT)]
*/

/*
Notes for camera position:
camera position goes reversed direction with vertex positions.
vertex positions are fixed, but camera position can be adjusted.
camera position also follows the diagram's rule.
y:
[-]
[...]
[+]
x:
[-][...][+]
*/

/*
Notes for center of mass in 3D's triangle:
Three vertices of triangle:
A(x1, y1, z1)
B(x2, y2, z2)
C(x3, y3, z3)
------------------------------
Gx = (x1 + x2 + x3) / 3
Gy = (y1 + y2 + y3) / 3
Gz = (z1 + z2 + z3) / 3
The coordinate of center of mass of triangle: (Gx, Gy, Gz).
*/

bool key = true; // Initial switch

glm::vec3 upMatrix = glm::vec3(0.0f, -0.1f, 0.0f); // y-
glm::vec3 downMatrix = glm::vec3(0.0f, 0.1f, 0.0f); // y+
glm::vec3 rightMatrix = glm::vec3(-0.1f, 0.0f, 0.0f); // x-
glm::vec3 leftMatrix = glm::vec3(0.1f, 0.0f, 0.0f); // x+

int angleInDegrees = 1;
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
/*
[1][2][3]
[1][2][3] * (inner product) innerProductIdentityMatrix = inputMatrix (mapping multiplication)
[1][2][3]
*/
glm::mat3 innerProductIdentityMatrix = glm::mat3(
	glm::vec3(1.0f, 1.0f, 1.0f),
	glm::vec3(1.0f, 1.0f, 1.0f),
	glm::vec3(1.0f, 1.0f, 1.0f)
);
/*
[1][2][3]
[1][2][3] (*) (outer product) outerProductIdentityMatrix = inputMatrix
[1][2][3]
                  -
([1][2][3])      |1| 0 0   [1] [] []
[1][2][3]   (*)  |0| 1 0 = [] [] []
[1][2][3]        |0| 0 1   [] [] []
                  -
*/
glm::mat3 outerProductIdentityMatrix = glm::mat3( // glm::mat3(1.0f)
	glm::vec3(1.0f, 0.0f, 0.0f),
	glm::vec3(0.0f, 1.0f, 0.0f),
	glm::vec3(0.0f, 0.0f, 1.0f)
);

float zoomIn = 0.1;
float zoomOut = -0.1;

std::vector<float> interpolateSingleFloats(float from, float to, float numberOfValues) {
	std::vector<float> v;
	float tmp = to - from; // 'tmp' could be positive or negative
	int num = numberOfValues - 1;
	int numtmp = numberOfValues - 2; // 'from' and 'to' have been pushed back to v
	float arithmetic = tmp / num;
	v.push_back(from);
	for (float i = 0.0; i < numtmp; i++) {
		from = from + arithmetic; // 'from' + positive / negative number
		v.push_back(from);
	}
	v.push_back(to);
	return v;
}

std::vector<glm::vec3> interpolateThreeElementValues(glm::vec3 from, glm::vec3 to, float numberOfValues) {
	std::vector<glm::vec3> v;
	std::vector<float> resultX;
	std::vector<float> resultY;
	std::vector<float> resultZ;
	resultX = interpolateSingleFloats(from.x, to.x, numberOfValues);
	resultY = interpolateSingleFloats(from.y, to.y, numberOfValues);
	resultZ = interpolateSingleFloats(from.z, to.z, numberOfValues);
	for (float i = 0.0; i < numberOfValues; i++) {
		glm::vec3 tmpA;
		tmpA.x = resultX[i];
		tmpA.y = resultY[i];
		tmpA.z = resultZ[i];
		v.push_back(tmpA);
	}
	return v;
}

glm::mat3 lookAt(glm::vec3 eye, glm::vec3 target, glm::vec3 up) { // const glm::vec3& eye, const glm::vec3& target, const glm::vec3& up
    glm::vec3 zaxis = glm::normalize(eye - target); // Forward
    glm::vec3 xaxis = glm::normalize(glm::cross(up, zaxis)); // Right
    glm::vec3 yaxis = glm::cross(zaxis, xaxis); // Up
    glm::mat3 rotationMatrix = glm::mat3(
        xaxis.x, yaxis.x, zaxis.x,
        xaxis.y, yaxis.y, zaxis.y,
        xaxis.z, yaxis.z, zaxis.z
    );
    return rotationMatrix;
}

void draw(DrawingWindow& window) {
	window.clearPixels();
	for (size_t y = 0; y < window.height; y++) {
		for (size_t x = 0; x < window.width; x++) {
			// Test black background
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
		if (event.key.keysym.sym == SDLK_LEFT) { // Camera left
			std::cout << "LEFT" << std::endl;
			window.clearPixels();
			MyModel myModel;
			myModel.resetDepthBuffer();
			cameraPosition = cameraPosition + leftMatrix;
			myModel.rasterisedRender(window);
		} else if (event.key.keysym.sym == SDLK_RIGHT) {// Camera right
			std::cout << "RIGHT" << std::endl;
			window.clearPixels();
			MyModel myModel;
			myModel.resetDepthBuffer();
			cameraPosition = cameraPosition + rightMatrix;
			myModel.rasterisedRender(window);
		} else if (event.key.keysym.sym == SDLK_UP) { // Camera up
			std::cout << "UP" << std::endl;
			window.clearPixels();
			MyModel myModel;
			myModel.resetDepthBuffer();
			cameraPosition = cameraPosition + upMatrix;
			myModel.rasterisedRender(window);
		} else if (event.key.keysym.sym == SDLK_DOWN) { // Camera down
			std::cout << "DOWN" << std::endl;
			window.clearPixels();
			MyModel myModel;
			myModel.resetDepthBuffer();
			cameraPosition = cameraPosition + downMatrix;
			myModel.rasterisedRender(window);
		} else if (event.key.keysym.sym == SDLK_u) { // u: Draw random stroked (unfilled) triangle
			std::cout << "u" << "/" << "U" << ":" << "Stroked triangle" << std::endl;
			MySDL mySDL;
			CanvasTriangle triangle;
			Colour c;
			mySDL.strokedTriangle(window, triangle, c);
		} else if (event.key.keysym.sym == SDLK_f) { // f: Draw random filled triangle
			std::cout << "f" << "/" << "F" << ":" << "Filled triangle" << std::endl;
			MySDL mySDL;
			CanvasTriangle triangle;
			Colour c;
			mySDL.filledTriangle(window, triangle, c);
		} else if (event.key.keysym.sym == SDLK_o) { // o: Orbit model
			std::cout << "o" << "/" << "O" << ":" << "Orbit" << std::endl;
			// A bool can be used to judge orbit or not...
			// key = -key;
			/*
							[camera]->
				^
				|	
			[camera]		[box]    	[camera]
											|
											v
						<-[camera]

			*/
			window.clearPixels();
			MyModel myModel;
			myModel.resetDepthBuffer();
			cameraPosition = cameraPosition * rotationMatrixYY; // Model goes right, camera goes left
			glm::mat3 rotationMatrix = lookAt(cameraPosition, targetPoint, upVector);
			cameraPosition = cameraPosition * rotationMatrix;
			myModel.rasterisedRender(window);
		} else if (event.key.keysym.sym == SDLK_g) { // g
			std::cout << "g" << "/" << "G" << ":" << "Test: Z axis" << std::endl;
			window.clearPixels();
			MyModel myModel;
			myModel.resetDepthBuffer();
			cameraPosition = cameraPosition * rotationMatrixZ;
			myModel.rasterisedRender(window);
		} else if (event.key.keysym.sym == SDLK_r) { // r: Reset model
			std::cout << "r" << "/" << "R" << ":" << "Reset" << std::endl;
			window.clearPixels();
			MyModel myModel;
			myModel.resetDepthBuffer();
			myModel.resetModel();
			myModel.rasterisedRender(window);
		} else if (event.key.keysym.sym == SDLK_z) { // z: Zoom in
			std::cout << "z" << "/" << "Z" << ":" << "Zoom in" << std::endl;
			window.clearPixels();
			MyModel myModel;
			myModel.resetDepthBuffer();
			focalLength = focalLength + zoomIn;
			myModel.rasterisedRender(window);
		} else if (event.key.keysym.sym == SDLK_x) { // x: Zoom out
			std::cout << "x" << "/" << "X" << ":" << "Zoom out" << std::endl;
			window.clearPixels();
			MyModel myModel;
			myModel.resetDepthBuffer();
			focalLength = focalLength + zoomOut;
			if (focalLength >= 0.0f) { // Normal model
				myModel.rasterisedRender(window);
			} else { // Disappeared model / Reversed model disabled
				std::cout << "Disappeared model disabled | Reversed model disabled" << std::endl;
			}
		} else if (event.key.keysym.sym == SDLK_a) { // Rotating the camera in the Y axis (panning)
			std::cout << "a" << "/" << "A" << ":" << "Panning: Y axis LEFT" << std::endl;
			window.clearPixels();
			MyModel myModel;
			myModel.resetDepthBuffer();
			cameraPosition = cameraPosition * rotationMatrixY; // Left
			myModel.rasterisedRender(window);
		} else if (event.key.keysym.sym == SDLK_d) { // Rotating the camera in the Y axis (panning)
			std::cout << "d" << "/" << "D" << ":" << "Panning: Y axis RIGHT" << std::endl;
			window.clearPixels();
			MyModel myModel;
			myModel.resetDepthBuffer();
			cameraPosition = cameraPosition * rotationMatrixYY; // Right
			myModel.rasterisedRender(window);
		} else if (event.key.keysym.sym == SDLK_w) { // Rotating the camera in the X axis (tilting)
			std::cout << "w" << "/" << "W" << ":" << "Tilting: X axis UP" << std::endl;
			window.clearPixels();
			MyModel myModel;
			myModel.resetDepthBuffer();
			cameraPosition = cameraPosition * rotationMatrixX; // Up
			myModel.rasterisedRender(window);
		} else if (event.key.keysym.sym == SDLK_s) { // Rotating the camera in the X axis (tilting)
			std::cout << "s" << "/" << "S" << ":" << "Tilting: X axis DOWN" << std::endl;
			window.clearPixels();
			MyModel myModel;
			myModel.resetDepthBuffer();
			cameraPosition = cameraPosition * rotationMatrixXX; // Down
			myModel.rasterisedRender(window);
		}
	} else if (event.type == SDL_MOUSEBUTTONDOWN) {
		window.savePPM("output.ppm");
		window.saveBMP("output.bmp");
	}
}

int main(int argc, char* argv[]) {
	DrawingWindow window = DrawingWindow(WIDTH, HEIGHT, false);
	SDL_Event event;
	draw(window);
	MyModel myModel;
	// myModel.wireFrameRender(window); // Wireframe render
	myModel.rasterisedRender(window); // Rasterised render
	while (true) {
		// if (key == true) {
		// 	// orbit
		// } else { // key == false
		// 	// stop orbit
		// }
		// We MUST poll for events - otherwise the window will freeze !
		if (window.pollForInputEvents(event)) handleEvent(event, window);
		// Need to render the frame at the end, or nothing actually gets shown on the screen !
		window.renderFrame();
	}
}
