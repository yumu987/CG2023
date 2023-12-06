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

// Finite state machine of rendering modes
enum class State {
    INIT, // initialRender
    WIREFRAME, // wireFrameRender
    RASTERISED, // rasterisedRender
	RAYTRACING // ray tracing & shadows
};

// Finite state machine of orbit
enum class Orb {
	OFF, // off
	ON // on
};

// Initial state
State mode = State::INIT;

// Initial state
Orb orbMode = Orb::OFF;

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
			} else if (mode == State::RAYTRACING) {
				// Ray tracing
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
			} else if (mode == State::RAYTRACING) {
				// Ray tracing
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
			} else if (mode == State::RAYTRACING) {
				// Ray tracing
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
			} else if (mode == State::RAYTRACING) {
				// Ray tracing
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
		} else if (event.key.keysym.sym == SDLK_o) { // o: Orbit model (switch to control ON / OFF)
			/*
							[camera]->
				^
				|	
			[camera]		[box]    	[camera]
											|
											v
						<-[camera]

			*/
			if (orbMode == Orb::ON) {
				std::cout << "o" << "/" << "O" << ":" << "Orbit: turn off!" << std::endl;
				orbMode = Orb::OFF;
			} else if (orbMode == Orb::OFF) {
				std::cout << "o" << "/" << "O" << ":" << "Orbit: turn on!" << std::endl;
				orbMode = Orb::ON;
			}
		} else if (event.key.keysym.sym == SDLK_g) { // g: reset ALL!
			std::cout << "g" << "/" << "G" << ":" << "Reset all to black background!" << std::endl;
			draw(window);
			mode = State::INIT;
			orbMode = Orb::OFF;
		} else if (event.key.keysym.sym == SDLK_r) { // r: Reset model to rasterised
			std::cout << "r" << "/" << "R" << ":" << "Reset to rasterised" << std::endl;
			window.clearPixels();
			MyModel myModel;
			myModel.resetDepthBuffer();
			myModel.resetModel();
			myModel.rasterisedRenderOrbit(window);
			mode = State::RASTERISED;
			orbMode = Orb::OFF;
		} else if (event.key.keysym.sym == SDLK_t) { // t: Reset model to wireframe
			std::cout << "t" << "/" << "T" << ":" << "Reset to wireframe" << std::endl;
			window.clearPixels();
			MyModel myModel;
			myModel.resetDepthBuffer();
			myModel.resetModel();
			myModel.wireFrameRenderOrbit(window);
			mode = State::WIREFRAME;
			orbMode = Orb::OFF;
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
			} else if (mode == State::RAYTRACING) {
				// Ray tracing
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
			} else if (mode == State::RAYTRACING) {
				// Ray tracing
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
			} else if (mode == State::RAYTRACING) {
				// Ray tracing
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
			} else if (mode == State::RAYTRACING) {
				// Ray tracing
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
			} else if (mode == State::RAYTRACING) {
				// Ray tracing
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
			} else if (mode == State::RAYTRACING) {
				// Ray tracing
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
		} else if (event.key.keysym.sym == SDLK_3) {
			std::cout << "3" << ":" << "RayTraceRender" << std::endl;
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
	// myModel.rasterisedRenderOrbit(window); // Rasterised render with orbit
	// mode = State::RASTERISED;
	// Raytrace render
	myModel.drawRasterisedScene(window); // Ray trace render with orbit
	mode = State::RAYTRACING;
	while (true) {
		// Orbit
		if (orbMode == Orb::ON) {
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
			} else if (mode == State::RAYTRACING) {
				// Ray tracing
			}
		}
		// We MUST poll for events - otherwise the window will freeze !
		if (window.pollForInputEvents(event)) handleEvent(event, window);
		// Need to render the frame at the end, or nothing actually gets shown on the screen !
		window.renderFrame();
	}
	// return 0; // Not reachable
}

// #include <CanvasTriangle.h>
// #include <DrawingWindow.h>
// #include <Utils.h>
// #include <fstream>
// #include <vector>
// #include <glm/glm.hpp>
// #include <CanvasPoint.h>
// #include <Colour.h>
// #include <math.h>
// #include <CanvasTriangle.h>
// #include <TextureMap.h>
// #include <ModelTriangle.h>
// #include <iostream>
// #include <Utils.h>
// #include <map>
// #include <RayTriangleIntersection.h>

// #define WIDTH 320
// #define HEIGHT 240
// #define objFileName "cornell-box.obj"
// #define mtlFileName "cornell-box.mtl"
 

// using namespace std;
// using namespace glm;
// vector<float> depthInit = vector<float>(WIDTH* HEIGHT, 0.0);
// float move_X = 0.0;
// float move_Y = 0.0;
// float roate_Y = 0.0;
// vec3 camPos = vec3(0, 0, 4);
// mat3 canROt = glm::mat3(-1, 0, 0, 0, 1, 0, 0, 0, 1);
// vec3 lightSource = vec3(0.0, 0.5, 0.5);




// vector<float> interpolateSingleFloats(float from, float to, int numberOfValues) {
// 	vector<float> result(numberOfValues);
// 	for (int i = 0; i < numberOfValues; i++)
// 	{
// 		result[i] = from + ((to - from) / (numberOfValues - 1)) * i;
// 	}
// 	return result;
// }

// vector<CanvasPoint> interpolateSingleFloatsTexture(CanvasPoint from, CanvasPoint to,float numberOfValues) {
// 	vector<CanvasPoint> result;
//     float xDiff = to.x - from.x;
// 	float yDiff = to.y - from.y;
// 	float xStepSize = xDiff / (numberOfValues-1);
// 	float yStepSize = yDiff / (numberOfValues-1);

// 	float xDiffTexture = to.texturePoint.x - from.texturePoint.x;
// 	float yDiffTexture = to.texturePoint.y - from.texturePoint.y;
//     float xStepSizeTexture = xDiffTexture / (numberOfValues-1);
// 	float yStepSizeTexture = yDiffTexture / (numberOfValues-1);

// 	for (float i = 0; i < numberOfValues; i++)
// 	{
// 		float x = from.x + (xStepSize * i);
// 		float y = from.y + (yStepSize * i);
// 		float xTexture = from.texturePoint.x + (xStepSizeTexture * i);
// 		float yTexture = from.texturePoint.y + (yStepSizeTexture * i);
// 		CanvasPoint canvansPoint(round(x),round(y));
// 		canvansPoint.texturePoint = TexturePoint(round(xTexture),round(yTexture));
// 		result.push_back(canvansPoint);
// 	}
// 	return result;
// }

// vector<vec3> interpolateThreeElementValues(vec3 from, vec3 to, int numberOfValues) {
// 	vector<vec3> result(numberOfValues);
// 	vector<float> resultX = interpolateSingleFloats(from[0], to[0], numberOfValues);
// 	vector<float> resultY = interpolateSingleFloats(from[1], to[1], numberOfValues);
// 	vector<float> resultZ = interpolateSingleFloats(from[2], to[2], numberOfValues);
// 	for (int i = 0; i < numberOfValues; i++)
// 	{   
// 		result[i].x = resultX[i];
// 		result[i].y = resultY[i];
// 		result[i].z = resultZ[i];
// 	}
// 	return result;
// }

// void drawLine(CanvasPoint from, CanvasPoint to, DrawingWindow &window, Colour input_colour) {
// 	float xDiff = to.x - from.x;
// 	float yDiff = to.y - from.y;
// 	float numberOfSteps = std::max(abs(xDiff),abs(yDiff));
// 	float xStepSize = xDiff / numberOfSteps;
// 	float yStepSize = yDiff / numberOfSteps;
// 	float zStepSize = (to.depth - from.depth) / numberOfSteps;
// 	uint32_t colour = (255 << 24) + (input_colour.red << 16) + (input_colour.green << 8) + input_colour.blue;
// 	for (float i = 0.0; i < numberOfSteps; i++) {
// 		float x = from.x + (xStepSize * i);
// 		float y = from.y + (yStepSize * i);
// 		float z = from.depth + (zStepSize * i);
// 		if ((round(x) >= 0 && round(x) < WIDTH)&&(round(y)>=0 && round(y) < HEIGHT))
// 		{
// 			if (z < depthInit[round(y) * WIDTH + round(x)])
// 			{
// 				continue;
// 			}
// 			else {
// 				float a = round(y) * WIDTH + round(x);
// 				depthInit[round(y) * WIDTH + round(x)] = z;
// 				window.setPixelColour(round(x), round(y), colour);
// 			}
// 		}
// 		else {
// 			continue;
// 		}
// 	}
// }

// void drawTriangle(CanvasTriangle triangleVertice, DrawingWindow &window, Colour colour) {
// 	drawLine(triangleVertice.v0(), triangleVertice.v1(), window, colour);
// 	drawLine(triangleVertice.v0(), triangleVertice.v2(), window, colour);
// 	drawLine(triangleVertice.v1(), triangleVertice.v2(), window, colour
// );
// }

// CanvasTriangle sortTriangle(CanvasTriangle triangleVertice){
//     if (triangleVertice.v0().y > triangleVertice.v1().y) {
//         swap(triangleVertice.v0(), triangleVertice.v1());
//     }
//     if (triangleVertice.v1().y > triangleVertice.v2().y) {
//         swap(triangleVertice.v1(), triangleVertice.v2());
//     }
//     if (triangleVertice.v0().y > triangleVertice.v1().y) {
//         swap(triangleVertice.v0(), triangleVertice.v1());
//     }
// 	return triangleVertice;
// }

// CanvasPoint findXDistance(CanvasPoint v0, CanvasPoint v1,float y){
// 	float xDiff = v1.x - v0.x;
// 	float yDiff = v1.y - v0.y;
// 	float numberOfSteps = std::max(abs(xDiff),abs(yDiff));
// 	float xStepSize = xDiff / numberOfSteps;
// 	float yStepSize = yDiff / numberOfSteps;
// 	float zStepSize = (v1.depth - v0.depth) / numberOfSteps;
// 	float i = (y - v0.y) / yStepSize;
// 	float z = v0.depth + (zStepSize * i);
//     CanvasPoint result(round(v0.x + (xStepSize * i)),round(y),z);
// 	return result;
// }

// CanvasPoint findForthPoint(CanvasTriangle triangleVertice){
// 	CanvasPoint point = findXDistance(triangleVertice.v0(),triangleVertice.v2(),triangleVertice.v1().y);
// 	return point;
// }

// void fillHalfTriangle(CanvasPoint from,CanvasPoint to,CanvasPoint findXFrom,CanvasPoint FindXTo, DrawingWindow &window,Colour colour){
// 	float xDiff = to.x - from.x;
// 	float yDiff = to.y - from.y;
// 	float numberOfSteps = std::max(abs(xDiff),abs(yDiff));
// 	float xStepSize = xDiff/ numberOfSteps;
// 	float yStepSize = yDiff / numberOfSteps;
// 	float zStepSzie = (to.depth - from.depth) / numberOfSteps;
// 	for (float i = 0; i < numberOfSteps; i++){
// 		float x = from.x + (xStepSize * i);
// 		float y = from.y + (yStepSize * i);
// 		float z = from.depth + (zStepSzie * i);
// 		CanvasPoint x1 = findXDistance(findXFrom,FindXTo,round(y));
// 		drawLine(CanvasPoint(round(x),round(y),z),x1,window,colour);
// 	}
// }

// void fillTriangle(CanvasTriangle triangleVertice,DrawingWindow &window, Colour colour){
// 	CanvasTriangle sortedTriangle = sortTriangle(triangleVertice);
// 	CanvasPoint ForthPoint = findForthPoint(sortedTriangle);
// 	//fill half 
// 	fillHalfTriangle(sortedTriangle.v0(),ForthPoint,sortedTriangle.v0(),sortedTriangle.v1(),window,colour);
// 	//fill other half
//     fillHalfTriangle(ForthPoint,sortedTriangle.v2(),sortedTriangle.v1(),sortedTriangle.v2(),window,colour);
// 	//draw white line
	
// 	drawTriangle(triangleVertice,window, colour);
// }

// void textureLineMapper(TextureMap texture,CanvasPoint from,CanvasPoint to,DrawingWindow &window){
// 	vector<CanvasPoint> points = interpolateSingleFloatsTexture(from,to,abs(to.x-from.x)+1);
// 	if (points.size() == 1) {
// 		uint32_t colour(texture.pixels[from.texturePoint.y * texture.width + from.texturePoint.x]);
// 		window.setPixelColour(from.x, from.y, colour);
// 		return;
// 	}	

// 	for (auto i : points)
// 	{
// 		uint32_t colour(texture.pixels[i.texturePoint.y*texture.width+i.texturePoint.x]);
// 		window.setPixelColour(round(i.x), round(i.y), colour);
// 	}
// }

// void fillHalfTriangleTexture(CanvasPoint from,CanvasPoint to,CanvasPoint findXFrom,CanvasPoint findXTo, DrawingWindow &window,TextureMap texture){
// 	float numberOfValue1 = abs(to.y - from.y);
// 	float numberOfValue2 = abs(findXTo.y - findXFrom.y);
// 	vector<CanvasPoint> points1 = interpolateSingleFloatsTexture(from, to,numberOfValue1+1);
// 	vector<CanvasPoint> points2 = interpolateSingleFloatsTexture(findXFrom, findXTo,numberOfValue2+1);
// 	for (float i = 0; i < numberOfValue1+1; i++)
// 	{  
// 		textureLineMapper(texture, points1[i], points2[i], window);
// 	}
// }

// CanvasPoint findForthPointTexture(CanvasPoint Forthpoint, CanvasTriangle storedVertircs) {
// 	vector<CanvasPoint> result; 
// 		result = interpolateSingleFloatsTexture(storedVertircs.v0(), storedVertircs.v2(), storedVertircs.v2().y - storedVertircs.v0().y + 1);
// 		for (auto i : result) {
// 			if (i.x == Forthpoint.x && i.y == Forthpoint.y) {
// 				Forthpoint.texturePoint = i.texturePoint;
// 			}
// 		}	
// 		return Forthpoint;
// }

// void fillTriangleTexture(CanvasTriangle triangleVertice,DrawingWindow &window, TextureMap texture){
// 	CanvasTriangle sortedTriangle = sortTriangle(triangleVertice);
// 	CanvasPoint ForthPoint = findForthPoint(sortedTriangle);
// 	ForthPoint = findForthPointTexture(ForthPoint, sortedTriangle);
// 	//fill half left to right
// 	fillHalfTriangleTexture(sortedTriangle.v0(),ForthPoint,sortedTriangle.v0(),sortedTriangle.v1(),window,texture);
// 	//fill other half
//     fillHalfTriangleTexture(ForthPoint,sortedTriangle.v2(),sortedTriangle.v1(),sortedTriangle.v2(),window,texture);
// 	//draw white line
// 	 /*Colour white(255,255,255);
// 	 drawTriangle(triangleVertice,window,white);*/
// }

// void drawTexture(DrawingWindow& window) {
// 	CanvasPoint point_1(160, 10);
// 	CanvasPoint point_2(300, 230);
// 	CanvasPoint point_3(10, 150);
// 	point_1.texturePoint = TexturePoint(195, 5);
// 	point_2.texturePoint = TexturePoint(395, 380);
// 	point_3.texturePoint = TexturePoint(65, 330);
// 	CanvasTriangle textureTriangle(point_1, point_2, point_3);
// 	TextureMap textureMap("texture.ppm");
// 	fillTriangleTexture(textureTriangle, window, textureMap);
// }

// //week4 task 3
// ModelTriangle getModelTriangle(string nextLine, vector<vec3> vertices, Colour colour) {
// 	auto splitList = split(nextLine, ' ');
// 	vec3 f1 = vertices[stoi(splitList[1])-1];
// 	vec3 f2 = vertices[stoi(splitList[2])-1];
// 	vec3 f3 = vertices[stoi(splitList[3])-1];
// 	return ModelTriangle(f1, f2, f3, colour);
// }
// //week task 2
// map<string, Colour> readMaterialFiles() {
// 	ifstream inputStream(mtlFileName);
// 	string nextLine;
// 	map<string, Colour> colourMap;
// 	string currentColour;
// 	if (!inputStream.is_open())
// 	{
// 		std::cout << "the file open fail" << endl;
// 	}

// 	while (!inputStream.eof())
// 	{
// 		getline(inputStream, nextLine);
// 		if (nextLine.find("newmtl") != string::npos)
// 		{
// 			auto splitList = split(nextLine, ' ');
// 			currentColour=splitList[1];
// 		}
// 		else if (nextLine.size() == 0)
// 		{
// 			continue;
// 		}
// 		else if (nextLine.find("Kd") != string::npos) {
// 			auto splitList = split(nextLine, ' ');
// 			int red = stof(splitList[1])*255;
// 			int green = stof(splitList[2]) * 255;
// 			int blue = stof(splitList[3]) * 255;
// 			colourMap[currentColour] = Colour(red, green, blue);
// 		}
// 	}
// 	inputStream.close();
// 	return colourMap;
// }

// vector<ModelTriangle> readGeometryFiles() {
// 	ifstream inputStream(objFileName);
// 	string nextLine;
// 	vector<ModelTriangle> modelTriangles;
// 	vector<string> objNames;
// 	string colours;
// 	vector<vec3> vertices;
// 	string currentColour;
// 	map<string, Colour> colourMap = readMaterialFiles();

	
// 	if (!inputStream.is_open())
// 	{
// 		std::cout << "the file open fail" << endl;
// 	}

// 	while (!inputStream.eof())
// 	{
// 		getline(inputStream, nextLine);
// 		if (nextLine.find("mtllib") != string :: npos)
// 		{
// 			continue;
// 		}
// 		else if (nextLine.size() == 0)
// 		{
// 			continue;
// 		}
// 		else if (nextLine.at(0) == 'o'){
// 			auto splitList = split(nextLine, ' ');
// 			objNames.push_back(splitList[1]);
// 		}
// 		else if (nextLine.find("usemtl") != string::npos)
// 		{
// 			auto splitList = split(nextLine, ' ');
// 			currentColour = splitList[1];
// 		}
// 		else if (nextLine.at(0) == 'v')
// 		{
// 			auto splitList = split(nextLine, ' ');
// 			float v1 = stof(splitList[1]) * 0.35;
// 			float v2 = stof(splitList[2]) * 0.35;
// 			float v3 = stof(splitList[3]) * 0.35;
// 			vec3 tmp(v1, v2, v3);
// 			vertices.push_back(tmp);
// 		}
// 		else if (nextLine.at(0) == 'f')
// 		{  
// 			auto splitList = split(nextLine, ' ');
// 			Colour usedColour;
// 			if (colourMap.find(currentColour) != colourMap.end()) {
// 				usedColour = colourMap[currentColour];
// 			}
// 			modelTriangles.push_back(getModelTriangle(nextLine,vertices,usedColour));
// 		}
// 	}
// 	inputStream.close();
// 	return modelTriangles;
// }
// //week4 task5
// CanvasPoint getCanvasIntersectionPoint(vec3 cameraPosition, vec3 vertexPosition, float focalLength) {
// 	vec3 relativePostion = vertexPosition - cameraPosition;
// 	relativePostion = canROt * relativePostion;
// 	float u = (focalLength * relativePostion.x / relativePostion.z) * 175 + WIDTH /2;
// 	float v = (focalLength * relativePostion.y / relativePostion.z)  * 175 + HEIGHT /2;
// 	float depth =  1/ -relativePostion.z;
// 		return CanvasPoint(round(u),round(v), depth);
// }
// //week4 task6
// void drawPointcloudRender(DrawingWindow& windows) {
// 	vector<ModelTriangle> modelTriangles = readGeometryFiles();
// 	for (float i = 0; i < modelTriangles.size(); i++)
// 	{
// 		for (float j = 0; j < modelTriangles[i].vertices.size(); j++) {
// 			CanvasPoint tmp = getCanvasIntersectionPoint(vec3(0, 0, 4), modelTriangles[i].vertices[j], 2);
// 			windows.setPixelColour(round(tmp.x), round(tmp.y), UINT32_MAX);
// 		}
// 	}
// }
// //week4 task7
// void drawWireframeRender(DrawingWindow& windows) {
// 	vector<ModelTriangle> modelTriangles = readGeometryFiles();
// 	for (float i = 0; i < modelTriangles.size(); i++)
// 	{
// 		CanvasPoint tmp1 = getCanvasIntersectionPoint(vec3(0, 0, 4), modelTriangles[i].vertices[0], 2);
// 		CanvasPoint tmp2 = getCanvasIntersectionPoint(vec3(0, 0, 4), modelTriangles[i].vertices[1], 2);
// 		CanvasPoint tmp3 = getCanvasIntersectionPoint(vec3(0, 0, 4), modelTriangles[i].vertices[2], 2);
// 		drawTriangle(CanvasTriangle(tmp1, tmp2, tmp3), windows, modelTriangles[i].colour);
// 	}
// }
// //week4 task8
// void drawRasterisedRender(DrawingWindow &windows) {
// 	vector<ModelTriangle> modelTriangles = readGeometryFiles();
// 	for (float i = 0; i < modelTriangles.size(); i++)
// 	{
// 		CanvasPoint tmp1 = getCanvasIntersectionPoint(vec3(0, 0, 4), modelTriangles[i].vertices[0], 2);
// 		CanvasPoint tmp2 = getCanvasIntersectionPoint(vec3(0, 0, 4), modelTriangles[i].vertices[1], 2);
// 		CanvasPoint tmp3 = getCanvasIntersectionPoint(vec3(0, 0, 4), modelTriangles[i].vertices[2], 2);
// 		fillTriangle(CanvasTriangle(tmp1, tmp2, tmp3), windows, modelTriangles[i].colour);
// 	}
// }


// void translateCameraPostion(DrawingWindow& windows) {
// 	windows.clearPixels();
// 	vector<ModelTriangle> modelTriangles = readGeometryFiles();
// 	for (float i = 0; i < modelTriangles.size(); i++)
// 	{
// 		CanvasPoint tmp1 = getCanvasIntersectionPoint(camPos, modelTriangles[i].vertices[0], 2);
// 		CanvasPoint tmp2 = getCanvasIntersectionPoint(camPos, modelTriangles[i].vertices[1], 2);
// 		CanvasPoint tmp3 = getCanvasIntersectionPoint(camPos, modelTriangles[i].vertices[2], 2);
// 		fillTriangle(CanvasTriangle(tmp1, tmp2, tmp3), windows, modelTriangles[i].colour);
// 	}
// }

// void rotateX(DrawingWindow& windows,float x, float y) {
// 	float degrees_X = 0.01;
// 	//float tmpDegree = degrees_X * M_PI / 180.0;
// 	mat3 x_mat3 = mat3(1, 0, 0, 0, cos(degrees_X), -sin(degrees_X), 0, sin(degrees_X), cos(degrees_X));
// 	//cameraOrientation = cameraOrientation * x_mat3;
// 	camPos = camPos * x_mat3;
// 	windows.clearPixels();
// 	vector<ModelTriangle> modelTriangles = readGeometryFiles();
// 	for (float i = 0; i < modelTriangles.size(); i++)
// 	{
// 		CanvasPoint tmp1 = getCanvasIntersectionPoint(camPos, modelTriangles[i].vertices[0], 2);
// 		CanvasPoint tmp2 = getCanvasIntersectionPoint(camPos, modelTriangles[i].vertices[1], 2);
// 		CanvasPoint tmp3 = getCanvasIntersectionPoint(camPos, modelTriangles[i].vertices[2], 2);
// 		fillTriangle(CanvasTriangle(tmp1, tmp2, tmp3), windows, modelTriangles[i].colour);
// 	}

// }
// void rotateY(DrawingWindow& windows,float x, float y) {
// 	float degrees_Y = 0.01;
// 	//float tmpDegree = degrees_Y * M_PI / 180.0;
// 	mat3 y_mat3 = mat3(cos(degrees_Y), 0, sin(degrees_Y), 0, 1, 0, -sin(degrees_Y), 0, cos(degrees_Y));
// 	//cameraOrientation = cameraOrientation * y_mat3;
// 	camPos = camPos * y_mat3;
// 	windows.clearPixels();
// 	vector<ModelTriangle> modelTriangles = readGeometryFiles();
// 	for (float i = 0; i < modelTriangles.size(); i++)
// 	{
// 		CanvasPoint tmp1 = getCanvasIntersectionPoint(camPos, modelTriangles[i].vertices[0], 2);
// 		CanvasPoint tmp2 = getCanvasIntersectionPoint(camPos, modelTriangles[i].vertices[1], 2);
// 		CanvasPoint tmp3 = getCanvasIntersectionPoint(camPos, modelTriangles[i].vertices[2], 2);
// 		fillTriangle(CanvasTriangle(tmp1, tmp2, tmp3), windows, modelTriangles[i].colour);
// 	}
// }
// mat3 lookAtMid(glm::vec3 cameraPos) {
// 	vec3 forward = glm::normalize(cameraPos - vec3(0, 0, 0));
// 	vec3 right = glm::normalize(glm::cross(forward,vec3(0, 1, 0)));
// 	vec3 up = glm::normalize(glm::cross(right,forward));

// 	return mat3(right, up, forward);
// }

// void Orbit(DrawingWindow &window) {
// 	    window.clearPixels();
// 	    float degree = M_PI/180;
// 		vec3 relativePosition = camPos - vec3(0,0,0);

// 		float cosA = cos(degree);
// 		float sinA = sin(degree);

// 		vec3 newPosition;
// 		newPosition.x = relativePosition.x * cosA - relativePosition.z * sinA;
// 		newPosition.y = relativePosition.y;
// 		newPosition.z = relativePosition.x * sinA + relativePosition.z * cosA;
// 		camPos = newPosition;

// 		canROt = lookAtMid(camPos);

// 		vector<ModelTriangle> modelTriangles = readGeometryFiles();
// 		for (float i = 0; i < modelTriangles.size(); i++)
// 		{
// 			CanvasPoint tmp1 = getCanvasIntersectionPoint(camPos, modelTriangles[i].vertices[0], 2);
// 			CanvasPoint tmp2 = getCanvasIntersectionPoint(camPos, modelTriangles[i].vertices[1], 2);
// 			CanvasPoint tmp3 = getCanvasIntersectionPoint(camPos, modelTriangles[i].vertices[2], 2);
// 			fillTriangle(CanvasTriangle(tmp1, tmp2, tmp3), window, modelTriangles[i].colour);
// 		}

// }

// //week06 task02 
// RayTriangleIntersection getClosestVaildIntersection(vec3 cameraPostion, vec3 rayDirection, ModelTriangle triangle) {
// 	    RayTriangleIntersection result;
// 		vec3 e0 = triangle.vertices[1] - triangle.vertices[0];
// 		vec3 e1 = triangle.vertices[2] - triangle.vertices[0];
// 		vec3 SPVector = cameraPostion - triangle.vertices[0];
// 		mat3 DEMatrix(-rayDirection, e0, e1);
// 		vec3 possibleSolution = inverse(DEMatrix) * SPVector;
// 		float t = possibleSolution.x;
// 		float u = possibleSolution.y;
// 		float v = possibleSolution.z;
// 		vec3 intersectionPoint = cameraPostion + rayDirection * possibleSolution.x;
// 		float distanceFromCamera = t;
// 		float index = 0;
// 		if (t>0 && (u >= 0.0) && (u <= 1.0) && (v >= 0.0) && (v <= 1.0) && (u + v) <= 1.0)
// 		{
// 			index = 1;
// 		}
// 		result = RayTriangleIntersection(intersectionPoint, distanceFromCamera, triangle, index);
// 		//cout << result << endl;
// 		return result;
// }

// RayTriangleIntersection getClosestVaildIntersection(vec3 cameraPostion, vec3 rayDirection, vector<ModelTriangle> triangles) {
// 	RayTriangleIntersection closetVaildIntersection = RayTriangleIntersection(vec3(), INT_MAX, ModelTriangle(), 0);
// 	RayTriangleIntersection tmp;
// 	for (int i = 0; i < triangles.size(); i++)
// 	{
// 		tmp = getClosestVaildIntersection(cameraPostion, rayDirection, triangles[i]);
// 		if (closetVaildIntersection.distanceFromCamera > tmp.distanceFromCamera && tmp.triangleIndex == 1)
// 		{
// 			closetVaildIntersection = tmp;
// 			closetVaildIntersection.triangleIndex = i;
// 			closetVaildIntersection.intersectedTriangle.normal = triangles[i].normal;
// 		}
// 	}
// 	//cout << closetVaildIntersection << endl;
// 	return closetVaildIntersection;
// }
// vec3 convert(vec3 camera, vec2 point) {
// 	vec3 ray;
// 	vec3 vec3Point;
// 	int z = 1;
// 	vec3Point.x = z * (point.x - WIDTH / 2) / (45 * 2);
// 	vec3Point.y = -z * (point.y - HEIGHT / 2) / (45 * 2);
// 	vec3Point.z = 0;
// 	ray = vec3Point - camera;
// 	ray = normalize(ray);
// 	return ray;
// }


// //week06 task 04
// void drawRay(DrawingWindow& window) {
// 	vector<ModelTriangle> modelTriangles = readGeometryFiles();
	
// 	for (size_t i = 0; i < modelTriangles.size(); i++) {
// 		std::cout << i << std::endl;
// 		std::cout << modelTriangles[i].colour.name << std::endl;
// 		std::cout << "Red: " << modelTriangles[i].colour.red << std::endl;
// 		std::cout << "Green: " << modelTriangles[i].colour.green << std::endl;
// 		std::cout << "Blue: " << modelTriangles[i].colour.blue << std::endl;
// 		std::cout << modelTriangles[i] << std::endl;
// 		// for (size_t j = 0; j < vecModel[i].vertices.size(); j++) {
// 		//     std::cout << vecModel[i].vertices[j].x << std::endl;
// 		//     std::cout << vecModel[i].vertices[j].y << std::endl;
// 		//     std::cout << vecModel[i].vertices[j].z << std::endl;
// 		// }
// 	}

// 	uint32_t black = (255 << 24) + (int(0) << 16) + (int(0) << 8) + int(0);
// 	for (int j = 0; j < HEIGHT; j++)
// 	{
// 		for (int i = 0; i < WIDTH; i++) {
// 			uint32_t colour = 0;
// 			vec3 converter = convert(camPos, vec2(i, j));
// 			RayTriangleIntersection intersectionPoint = getClosestVaildIntersection(camPos,converter,modelTriangles);
// 			RayTriangleIntersection tmp = getClosestVaildIntersection(lightSource, normalize(intersectionPoint.intersectionPoint - lightSource), modelTriangles);
// 			if (intersectionPoint.triangleIndex == tmp.triangleIndex) {
// 				colour = (255 << 24) + (int(intersectionPoint.intersectedTriangle.colour.red) << 16) + (int(intersectionPoint.intersectedTriangle.colour.green) << 8) + int(intersectionPoint.intersectedTriangle.colour.blue);
// 			}
// 			else
// 			{
// 				colour = black;
// 			}
// 			window.setPixelColour(i, j, colour);
// 		}
// 	}
// 	cout << "end" << endl;
// }

// //week07 task 02
// float calculateIntensity(vec3 lightPosition, vec3 pointPosition) {
// 	vec3 distanceVec = lightPosition - pointPosition;
// 	float distance = length(distanceVec); 
// 	float intensity = 1.0f / (4.0f * M_PI * distance * distance); 
// 	return intensity;
// }

// Colour adjustPixelBrightness(Colour pixelColor, float intensity) {
// 	Colour result;
// 	result.red = pixelColor.red * intensity;
// 	result.green = pixelColor.green * intensity;
// 	result.blue = pixelColor.blue * intensity;
// 	return result;
// }

// vec3 calculateSurfaceNormal(const vec3& v0, const vec3& v1, const vec3& v2) {
// 	vec3 edge1 = v1 - v0;
// 	vec3 edge2 = v2 - v0;
// 	vec3 normalVector = cross(edge1, edge2);
// 	return normalVector;
// }

// float calculateAngleOfIncidence(vec3 normalVec, vec3 lightVec) {
// 	float reuslt = dot(normalVec, lightVec);
// 	return reuslt;
// }

// vector<ModelTriangle> MoelTriangleNormal() {
// 	vector<ModelTriangle> modelTriangles = readGeometryFiles();
// 	for (float i = 0; i < modelTriangles.size(); i++)
// 	{
// 		modelTriangles[i].normal = calculateSurfaceNormal(modelTriangles[i].vertices[0], modelTriangles[i].vertices[1], modelTriangles[i].vertices[2]);
// 	}
// 	return modelTriangles;
// }

// void drawLighting(DrawingWindow& window){
// 	vector<ModelTriangle> modelTriangles = MoelTriangleNormal();
// 	/*for (float i = 0; i < modelTriangles.size(); i++) {
// 		CanvasPoint tmp1 = getCanvasIntersectionPoint(camPos, modelTriangles[i].vertices[0], 2);
// 		CanvasPoint tmp2 = getCanvasIntersectionPoint(camPos, modelTriangles[i].vertices[1], 2);
// 		CanvasPoint tmp3 = getCanvasIntersectionPoint(camPos, modelTriangles[i].vertices[2], 2);
// 		float incidence = calculateIntensity(lightSource, );
// 		Colour drawColour = adjustPixelBrightness(modelTriangles[i].colour, incidence);
// 		if (incidence > 0.2){ 
// 			modelTriangles[i].colour = drawColour; 
// 		}
// 	}*/

// 	for (int j = 0; j < HEIGHT; j++)
// 	{
// 		for (int i = 0; i < WIDTH; i++) {
// 			uint32_t colour = 0;
// 			vec3 converter = convert(camPos, vec2(i, j));
// 			RayTriangleIntersection intersectionPoint = getClosestVaildIntersection(camPos, converter, modelTriangles);
// 			float intensity = calculateIntensity(lightSource, intersectionPoint.intersectionPoint);
// 			float incidence = calculateAngleOfIncidence(intersectionPoint.intersectedTriangle.normal, lightSource);
// 			float diffues = incidence * intensity;
// 			cout << diffues << endl;
// 			int red = 0, green = 0 , blue = 0;
// 			RayTriangleIntersection tmp = getClosestVaildIntersection(lightSource, normalize(intersectionPoint.intersectionPoint - lightSource), modelTriangles);
// 			if (intersectionPoint.triangleIndex == tmp.triangleIndex) {
// 				red = clamp(int(intersectionPoint.intersectedTriangle.colour.red * diffues),0,255);
// 				green =clamp(int(intersectionPoint.intersectedTriangle.colour.green * diffues), 0, 255);
// 				blue = clamp(int(intersectionPoint.intersectedTriangle.colour.blue * diffues), 0, 255);
// 			}
// 			else
// 			{
// 				red = int(intersectionPoint.intersectedTriangle.colour.red);
// 				green = int(intersectionPoint.intersectedTriangle.colour.green);
// 				blue = int(intersectionPoint.intersectedTriangle.colour.blue);
// 			}
// 			colour = (255 << 24) + (red << 16) + (green << 8) + blue;
// 			window.setPixelColour(i, j, colour);
// 		}
// 	}
// 	cout << "end" << endl;
// }




// void handleEvent(SDL_Event event, DrawingWindow &window) {
// 	if (event.type == SDL_KEYDOWN) {
// 		if (event.key.keysym.sym == SDLK_LEFT) {
// 			std::cout << "LEFT" << std::endl;
// 			fill(depthInit.begin(), depthInit.end(), 0.0);
// 			camPos.x = camPos.x + 1;
// 			translateCameraPostion(window);
// 		}
// 		else if (event.key.keysym.sym == SDLK_RIGHT) {
// 			std::cout << "RIGHT" << std::endl;
// 			fill(depthInit.begin(), depthInit.end(), 0.0);
// 			camPos.x = camPos.x - 1;
// 			translateCameraPostion(window);
// 		}
// 		else if (event.key.keysym.sym == SDLK_UP) {
// 			std::cout << "UP" << std::endl;
// 			fill(depthInit.begin(), depthInit.end(), 0.0);
// 			camPos.y = camPos.y - 1;
// 			translateCameraPostion(window);
// 		}
// 		else if (event.key.keysym.sym == SDLK_DOWN) {
// 			std::cout << "DOWN" << std::endl;
// 			fill(depthInit.begin(), depthInit.end(), 0.0);
// 			camPos.y = camPos.y + 1;
// 			translateCameraPostion(window);
// 		}
// 		else if (event.key.keysym.sym == SDLK_x)
// 		{
// 			std::cout << "rotate_X" << std::endl;
// 			fill(depthInit.begin(), depthInit.end(), 0.0);
// 			rotateX(window, move_X, move_Y);
// 		}
// 		else if (event.key.keysym.sym == SDLK_y)
// 		{
// 			std::cout << "rotate_Y" << std::endl;
// 			fill(depthInit.begin(), depthInit.end(), 0.0);
// 			rotateY(window, move_X, move_Y);
// 		}
// 		else if (event.key.keysym.sym == SDLK_o) {
// 			std::cout << "0" << std::endl;
// 			fill(depthInit.begin(), depthInit.end(), 0.0);
// 			Orbit(window);
// 		}
// 		else if (event.key.keysym.sym == SDLK_l) {
// 			std::cout << "l" << std::endl;
// 			fill(depthInit.begin(), depthInit.end(), 0.0);
// 			drawLighting(window);
// 		}
// 		else if (event.key.keysym.sym == SDLK_u) {
// 			CanvasPoint v0(rand() % window.width, rand() % window.height, 1), v1(rand() % window.width, rand() % window.height, 1), v2(rand() % window.width, rand() % window.height, 1);
// 			CanvasTriangle triangleVertice(v0, v1, v2);
// 			Colour colour(rand() % 256, rand() % 256, rand() % 256);
// 			drawTriangle(triangleVertice, window, colour);
// 			std::cout << "Draw triangle" << std::endl;
// 		}
// 		else if (event.key.keysym.sym == SDLK_f) {
// 			CanvasPoint v0(rand() % window.width, rand() % window.height, 1), v1(rand() % window.width, rand() % window.height, 1), v2(rand() % window.width, rand() % window.height, 1);
// 			CanvasTriangle triangleVertice(v0, v1, v2);
// 			Colour colour(rand() % 256, rand() % 256, rand() % 256);
// 			fillTriangle(triangleVertice, window, colour);
// 			std::cout << "Fill triangle" << std::endl;
// 		}
// 		else if (event.key.keysym.sym == SDLK_q)
// 		{
// 			window.clearPixels();
// 			camPos = vec3(0, 0, 4);
// 			fill(depthInit.begin(), depthInit.end(), 0.0);
// 			canROt = glm::mat3(-1, 0, 0, 0, 1, 0, 0, 0, 1);
// 			drawWireframeRender(window);
// 		}
// 		else if (event.key.keysym.sym == SDLK_w)
// 		{
// 			window.clearPixels();
// 			camPos = vec3(0, 0, 4);
// 			fill(depthInit.begin(), depthInit.end(), 0.0);
// 			canROt = glm::mat3(-1, 0, 0, 0, 1, 0, 0, 0, 1);
// 			drawRasterisedRender(window);
// 		}
// 		else if (event.key.keysym.sym == SDLK_r)
// 		{
// 			window.clearPixels();
// 			//camPos = vec3(0, 0, 4);
// 			fill(depthInit.begin(), depthInit.end(), 0.0);
// 			//canROt = glm::mat3(-1, 0, 0, 0, 1, 0, 0, 0, 1);
// 			drawRay(window);
// 		}
// 	}
// 	else if (event.type == SDL_MOUSEBUTTONDOWN) {
// 		window.savePPM("output.ppm");
// 		window.saveBMP("output.bmp");
// 	}
// }


// void drawRasterisedScene(DrawingWindow& window) {
// 	window.clearPixels();

// 	//Single Dimension Greyscale Interpolation
// 	//vector<float> result_single(window.width);
// 	//Orbit(window);
// 	// Two Dimensional Colour Interpolation
// 	vec3 topLeft(255, 0, 0);        // red 
// 	vec3 topRight(0, 0, 255);       // blue 
// 	vec3 bottomRight(0, 255, 0);    // green 
// 	vec3 bottomLeft(255, 255, 0);   // yellow

// 	vector<vec3> result_leftTop2Bottom(window.height);
// 	vector<vec3> result_rightTop2Bottom(window.height);
// 	result_leftTop2Bottom = interpolateThreeElementValues(topLeft, bottomLeft, window.height);
// 	result_rightTop2Bottom = interpolateThreeElementValues(topRight, bottomRight, window.height);

// 	vector<vec3> result_threeElement(window.width);

// 	for (size_t y = 0; y < window.height; y++) {
// 		//Single Dimension Greyscale Interpolation
// 		//result_single = interpolateSingleFloats(255, 0, window.width);

// 		//Two Dimensional Colour Interpolation
// 		result_threeElement = interpolateThreeElementValues(result_leftTop2Bottom[y], result_rightTop2Bottom[y], window.width);

// 		for (size_t x = 0; x < window.width; x++) {
// 			//Two Dimensional Colour Interpolation
// 			/*float red = result_threeElement[x].x;
// 			float green = result_threeElement[x].y;
// 			float blue = result_threeElement[x].z;*/
// 			//float red = 0;
// 			//float green = 0;
// 			//float blue = 0; 

// 			//uint32_t colour = (255 << 24) + (int(red) << 16) + (int(green) << 8) + int(blue);
// 			//window.setPixelColour(x, y, colour);
// 		}
// 	}
// }

// int main(int argc, char *argv[]) {
// 	DrawingWindow window = DrawingWindow(WIDTH, HEIGHT, false);
// 	SDL_Event event;
// 	int window_width = window.width - 1;
// 	int window_height = window.height - 1;

// 	//test drawLine ponnits
// 	CanvasPoint from_top_left_corner(0, 0);
// 	CanvasPoint to_center(int(window_width/2), int(window_height/2));
// 	CanvasPoint from_top_right(window_width, 0);
// 	CanvasPoint from_top_center(int(window_width/2), 0);
// 	CanvasPoint to_bottom_center(int(window_width /2), window_height);
// 	CanvasPoint from_test(int(window_width/3),int(window_height/2));
// 	CanvasPoint to_test(int(window_width-window.width/3),int(window_height/2));

    
// 	//drawTexture(window);
// 	//drawPointcloudRender(window);
// 	//drawRasterisedScene(window);
	
	
	
// 	drawRay(window);
	
// 	while (true) {
// 		// We MUST poll for events - otherwise the window will freeze !
// 		if (window.pollForInputEvents(event)) handleEvent(event, window);
// 		// Need to render the frame at the end, or nothing actually gets shown on the screen !
// 		window.renderFrame();
// 	}
// }
