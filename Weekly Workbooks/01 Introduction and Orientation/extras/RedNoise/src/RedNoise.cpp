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
#include "mysdl.hpp"
#include "mymodel.hpp"
#include "myray.hpp"
#include "mytexture.hpp"

#define WIDTH 320 * 2
#define HEIGHT 240 * 2
#define HALFWIDTH 160 * 2
#define HALFHEIGHT 120 * 2
#define OBJfilename "cornell-box.obj"
#define MTLfilename "cornell-box.mtl"

// Finite state machine
enum class State {
    INIT, // initialRender
    WIREFRAME, // wireFrameRender
    RASTERISED, // rasterisedRender
	RAYTRACING // ray tracing & shadows
};

// Initial state
State mode = State::INIT;

// float depthBuffer[WIDTH][HEIGHT];
// glm::vec3 cameraPosition (0.0f, 0.0f, 4.0f);
// float focalLength = 2.0f;
// glm::vec3 cameraToVertex (3.0f, 4.0f, 5.0f);
// glm::vec3 targetPoint (0.0f, 0.0f, 0.0f);
// glm::vec3 upVector (0.0f, 1.0f, 0.0f);

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

// glm::vec3 upMatrix = glm::vec3(0.0f, -0.1f, 0.0f); // y-
// glm::vec3 downMatrix = glm::vec3(0.0f, 0.1f, 0.0f); // y+
// glm::vec3 rightMatrix = glm::vec3(-0.1f, 0.0f, 0.0f); // x-
// glm::vec3 leftMatrix = glm::vec3(0.1f, 0.0f, 0.0f); // x+

// int angleInDegrees = 1;
// double angleInRadians = angleInDegrees * (M_PI / 180); // radians = degrees * (M_PI / 180)

// glm::mat3 rotationMatrixX = glm::mat3( // X axis, left
// 	glm::vec3(1.0f, 0.0f, 0.0f),
// 	glm::vec3(0.0f, std::cos(angleInRadians), -std::sin(angleInRadians)),
// 	glm::vec3(0.0f, std::sin(angleInRadians), std::cos(angleInRadians))
// );
// glm::mat3 rotationMatrixXX = glm::mat3( // X axis, right
// 	glm::vec3(1.0f, 0.0f, 0.0f),
// 	glm::vec3(0.0f, std::cos(angleInRadians), std::sin(angleInRadians)),
// 	glm::vec3(0.0f, -std::sin(angleInRadians), std::cos(angleInRadians))
// );
// glm::mat3 rotationMatrixY = glm::mat3( // Y axis, up
// 	glm::vec3(std::cos(angleInRadians), 0.0f, std::sin(angleInRadians)),
// 	glm::vec3(0.0f, 1.0f, 0.0f),
// 	glm::vec3(-std::sin(angleInRadians), 0.0f, std::cos(angleInRadians))
// );
// glm::mat3 rotationMatrixYY = glm::mat3( // Y axis, down
// 	glm::vec3(std::cos(angleInRadians), 0.0f, -std::sin(angleInRadians)),
// 	glm::vec3(0.0f, 1.0f, 0.0f),
// 	glm::vec3(std::sin(angleInRadians), 0.0f, std::cos(angleInRadians))
// );

// glm::mat3 rotationMatrixZ = glm::mat3( // Z axis, but no need to implement
// 	glm::vec3(std::cos(angleInRadians), -std::sin(angleInRadians), 0.0f),
// 	glm::vec3(std::sin(angleInRadians), std::cos(angleInRadians), 0.0f),
// 	glm::vec3(0.0f, 0.0f, 1.0f)
// );
/*
[1][2][3]
[1][2][3] * (inner product) innerProductIdentityMatrix = inputMatrix (mapping multiplication)
[1][2][3]
*/
// glm::mat3 innerProductIdentityMatrix = glm::mat3(
// 	glm::vec3(1.0f, 1.0f, 1.0f),
// 	glm::vec3(1.0f, 1.0f, 1.0f),
// 	glm::vec3(1.0f, 1.0f, 1.0f)
// );
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
// glm::mat3 outerProductIdentityMatrix = glm::mat3( // glm::mat3(1.0f)
// 	glm::vec3(1.0f, 0.0f, 0.0f),
// 	glm::vec3(0.0f, 1.0f, 0.0f),
// 	glm::vec3(0.0f, 0.0f, 1.0f)
// );
// float zoomIn = 0.1;
// float zoomOut = -0.1;

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

void draw(DrawingWindow& window) {
	window.clearPixels();
	for (size_t y = 0; y < window.height; y++) {
		for (size_t x = 0; x < window.width; x++) {
			// Test black background
			float red = 0.0;
			float green = 0.0;
			float blue = 0.0;
			// Pack colour into uint32_t (channel type: alpha + RGB) type: int
			uint32_t colour = (255 << 24) + (int(red) << 16) + (int(green) << 8) + int(blue); // Alpha + Red + Green + Blue
			window.setPixelColour(x, y, colour);
		}
	}
}

void handleEvent(SDL_Event event, DrawingWindow& window) {
	if (event.type == SDL_KEYDOWN) {
		if (event.key.keysym.sym == SDLK_LEFT) { // Camera left
			std::cout << "LEFT" << std::endl;
			if (mode == State::RASTERISED) {
				window.clearPixels();
				MyModel myModel;
				myModel.resetDepthBuffer();
				cameraPosition = cameraPosition + leftMatrix;
				myModel.rasterisedRenderOrbit(window);
			} else if (mode == State::WIREFRAME) {
				window.clearPixels();
				MyModel myModel;
				myModel.resetDepthBuffer();
				cameraPosition = cameraPosition + leftMatrix;
				myModel.wireFrameRenderOrbit(window);
			}
		} else if (event.key.keysym.sym == SDLK_RIGHT) {// Camera right
			std::cout << "RIGHT" << std::endl;
			if (mode == State::RASTERISED) {
				window.clearPixels();
				MyModel myModel;
				myModel.resetDepthBuffer();
				cameraPosition = cameraPosition + rightMatrix;
				myModel.rasterisedRenderOrbit(window);
			} else if (mode == State::WIREFRAME) {
				window.clearPixels();
				MyModel myModel;
				myModel.resetDepthBuffer();
				cameraPosition = cameraPosition + rightMatrix;
				myModel.wireFrameRenderOrbit(window);
			}
		} else if (event.key.keysym.sym == SDLK_UP) { // Camera up
			std::cout << "UP" << std::endl;
			if (mode == State::RASTERISED) {
				window.clearPixels();
				MyModel myModel;
				myModel.resetDepthBuffer();
				cameraPosition = cameraPosition + upMatrix;
				myModel.rasterisedRenderOrbit(window);
			} else if (mode == State::WIREFRAME) {
				window.clearPixels();
				MyModel myModel;
				myModel.resetDepthBuffer();
				cameraPosition = cameraPosition + upMatrix;
				myModel.wireFrameRenderOrbit(window);
			}
		} else if (event.key.keysym.sym == SDLK_DOWN) { // Camera down
			std::cout << "DOWN" << std::endl;
			if (mode == State::RASTERISED) {
				window.clearPixels();
				MyModel myModel;
				myModel.resetDepthBuffer();
				cameraPosition = cameraPosition + downMatrix;
				myModel.rasterisedRenderOrbit(window);
			} else if (mode == State::WIREFRAME) {
				window.clearPixels();
				MyModel myModel;
				myModel.resetDepthBuffer();
				cameraPosition = cameraPosition + downMatrix;
				myModel.wireFrameRenderOrbit(window);
			}
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
			/*
							[camera]->
				^
				|	
			[camera]		[box]    	[camera]
											|
											v
						<-[camera]

			*/
			if (mode == State::RASTERISED) {
				window.clearPixels();
				MyModel myModel;
				myModel.resetDepthBuffer();
				cameraPosition = cameraPosition * rotationMatrixYY; // Move camera towards left of the box (box goes right)
				cameraOrientation = myModel.lookAt(cameraPosition, targetPoint, upVector); // Let camera focus on the centre
				myModel.rasterisedRenderOrbit(window);
			} else if (mode == State::WIREFRAME) {
				window.clearPixels();
				MyModel myModel;
				myModel.resetDepthBuffer();
				cameraPosition = cameraPosition * rotationMatrixYY; // Move camera towards left of the box
				cameraOrientation = myModel.lookAt(cameraPosition, targetPoint, upVector); // Let camera focus on the centre
				myModel.wireFrameRenderOrbit(window);
			}
		} else if (event.key.keysym.sym == SDLK_g) { // g: reset ALL!
			std::cout << "g" << "/" << "G" << ":" << "Reset all to black background!" << std::endl;
			draw(window);
			mode = State::INIT;
		} else if (event.key.keysym.sym == SDLK_r) { // r: Reset model to rasterised
			std::cout << "r" << "/" << "R" << ":" << "Reset to rasterised" << std::endl;
			window.clearPixels();
			MyModel myModel;
			myModel.resetDepthBuffer();
			myModel.resetModel();
			myModel.rasterisedRenderOrbit(window);
			mode = State::RASTERISED;
		} else if (event.key.keysym.sym == SDLK_t) { // t: Reset model to wireframe
			std::cout << "t" << "/" << "T" << ":" << "Reset to wireframe" << std::endl;
			window.clearPixels();
			MyModel myModel;
			myModel.resetDepthBuffer();
			myModel.resetModel();
			myModel.wireFrameRenderOrbit(window);
			mode = State::WIREFRAME;
		} else if (event.key.keysym.sym == SDLK_z) { // z: Zoom in
			std::cout << "z" << "/" << "Z" << ":" << "Zoom in" << std::endl;
			if (mode == State::RASTERISED) {
				window.clearPixels();
				MyModel myModel;
				myModel.resetDepthBuffer();
				focalLength = focalLength + zoomIn;
				myModel.rasterisedRenderOrbit(window);
			} else if (mode == State::WIREFRAME) {
				window.clearPixels();
				MyModel myModel;
				myModel.resetDepthBuffer();
				focalLength = focalLength + zoomIn;
				myModel.wireFrameRenderOrbit(window);
			}
		} else if (event.key.keysym.sym == SDLK_x) { // x: Zoom out
			std::cout << "x" << "/" << "X" << ":" << "Zoom out" << std::endl;
			if (mode == State::RASTERISED) {
				window.clearPixels();
				MyModel myModel;
				myModel.resetDepthBuffer();
				focalLength = focalLength + zoomOut;
				if (focalLength >= 0.0f) { // Normal model
					myModel.rasterisedRenderOrbit(window);
				} else { // Disappeared model / Reversed model disabled
					std::cout << "Disappeared model disabled | Reversed model disabled" << std::endl;
				}
			} else if (mode == State::WIREFRAME) {
				window.clearPixels();
				MyModel myModel;
				myModel.resetDepthBuffer();
				focalLength = focalLength + zoomOut;
				if (focalLength >= 0.0f) { // Normal model
					myModel.wireFrameRenderOrbit(window);
				} else { // Disappeared model / Reversed model disabled
					std::cout << "Disappeared model disabled | Reversed model disabled" << std::endl;
				}
			}
		} else if (event.key.keysym.sym == SDLK_a) { // Rotating the camera in the Y axis (panning)
			std::cout << "a" << "/" << "A" << ":" << "Panning: Y axis LEFT" << std::endl;
			if (mode == State::RASTERISED) {
				window.clearPixels();
				MyModel myModel;
				myModel.resetDepthBuffer();
				cameraPosition = cameraPosition * rotationMatrixY; // Left
				myModel.rasterisedRenderOrbit(window);
			} else if (mode == State::WIREFRAME) {
				window.clearPixels();
				MyModel myModel;
				myModel.resetDepthBuffer();
				cameraPosition = cameraPosition * rotationMatrixY; // Left
				myModel.wireFrameRenderOrbit(window);
			}
		} else if (event.key.keysym.sym == SDLK_d) { // Rotating the camera in the Y axis (panning)
			std::cout << "d" << "/" << "D" << ":" << "Panning: Y axis RIGHT" << std::endl;
			if (mode == State::RASTERISED) {
				window.clearPixels();
				MyModel myModel;
				myModel.resetDepthBuffer();
				cameraPosition = cameraPosition * rotationMatrixYY; // Right
				myModel.rasterisedRenderOrbit(window);
			} else if (mode == State::WIREFRAME) {
				window.clearPixels();
				MyModel myModel;
				myModel.resetDepthBuffer();
				cameraPosition = cameraPosition * rotationMatrixYY; // Right
				myModel.wireFrameRenderOrbit(window);
			}
		} else if (event.key.keysym.sym == SDLK_w) { // Rotating the camera in the X axis (tilting)
			std::cout << "w" << "/" << "W" << ":" << "Tilting: X axis UP" << std::endl;
			if (mode == State::RASTERISED) {
				window.clearPixels();
				MyModel myModel;
				myModel.resetDepthBuffer();
				cameraPosition = cameraPosition * rotationMatrixX; // Up
				myModel.rasterisedRenderOrbit(window);
			} else if (mode == State::WIREFRAME) {
				window.clearPixels();
				MyModel myModel;
				myModel.resetDepthBuffer();
				cameraPosition = cameraPosition * rotationMatrixX; // Up
				myModel.wireFrameRenderOrbit(window);
			}
		} else if (event.key.keysym.sym == SDLK_s) { // Rotating the camera in the X axis (tilting)
			std::cout << "s" << "/" << "S" << ":" << "Tilting: X axis DOWN" << std::endl;
			if (mode == State::RASTERISED) {
				window.clearPixels();
				MyModel myModel;
				myModel.resetDepthBuffer();
				cameraPosition = cameraPosition * rotationMatrixXX; // Down
				myModel.rasterisedRenderOrbit(window);
			} else if (mode == State::WIREFRAME) {
				window.clearPixels();
				MyModel myModel;
				myModel.resetDepthBuffer();
				cameraPosition = cameraPosition * rotationMatrixXX; // Down
				myModel.wireFrameRenderOrbit(window);
			}
		} else if (event.key.keysym.sym == SDLK_1) {
			std::cout << "1" << ":" << "WireFrameRender" << std::endl;
			window.clearPixels();
			MyModel myModel;
			myModel.resetDepthBuffer();
			myModel.wireFrameRenderOrbit(window);
			mode = State::WIREFRAME;
		} else if (event.key.keysym.sym == SDLK_2) {
			std::cout << "2" << ":" << "RasterisedRender" << std::endl;
			window.clearPixels();
			MyModel myModel;
			myModel.resetDepthBuffer();
			myModel.rasterisedRenderOrbit(window);
			mode = State::RASTERISED;
		}
	} else if (event.type == SDL_MOUSEBUTTONDOWN) {
		std::cout << "Saved!" << std::endl;
		window.savePPM("output.ppm");
		window.saveBMP("output.bmp");
	}
}

int main(int argc, char* argv[]) {
	DrawingWindow window = DrawingWindow(WIDTH, HEIGHT, false);
	SDL_Event event;
	draw(window); // Initial render
	MyModel myModel;
	// Reset depth buffer and model
	myModel.resetDepthBuffer();
	myModel.resetModel();
	// Wireframe render
	// myModel.wireFrameRenderOrbit(window); // Wireframe render with orbit
	// mode = State::WIREFRAME;
	// Rasterised render
	myModel.rasterisedRenderOrbit(window); // Rasterised render with orbit
	mode = State::RASTERISED;
	while (true) {
		// We MUST poll for events - otherwise the window will freeze !
		if (window.pollForInputEvents(event)) handleEvent(event, window);
		// Need to render the frame at the end, or nothing actually gets shown on the screen !
		window.renderFrame();
	}
	// return 0; // Not reachable
}
