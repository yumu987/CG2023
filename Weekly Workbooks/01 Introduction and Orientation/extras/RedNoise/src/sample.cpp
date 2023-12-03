#include <iostream>
#include <vector>
#include <cmath>
#include <limits>

// Define a simple 3D vector class
struct Vec3 {
    float x, y, z;

    Vec3() : x(0), y(0), z(0) {}
    Vec3(float x, float y, float z) : x(x), y(y), z(z) {}

    Vec3 operator+(const Vec3& other) const {
        return Vec3(x + other.x, y + other.y, z + other.z);
    }

    Vec3 operator-(const Vec3& other) const {
        return Vec3(x - other.x, y - other.y, z - other.z);
    }

    Vec3 operator*(float scalar) const {
        return Vec3(x * scalar, y * scalar, z * scalar);
    }

    float dot(const Vec3& other) const {
        return x * other.x + y * other.y + z * other.z;
    }

    Vec3 cross(const Vec3& other) const {
        return Vec3(y * other.z - z * other.y, z * other.x - x * other.z, x * other.y - y * other.x);
    }

    float length() const {
        return std::sqrt(x * x + y * y + z * z);
    }

    Vec3 normalized() const {
        float len = length();
        return len > 0 ? *this * (1.0f / len) : *this;
    }
};

// Define a simple sphere class
struct Sphere {
    Vec3 center;
    float radius;

    Sphere(const Vec3& center, float radius) : center(center), radius(radius) {}

    // Ray-sphere intersection test
    bool intersect(const Vec3& rayOrigin, const Vec3& rayDirection, float& t) const {
        Vec3 oc = rayOrigin - center;
        float a = rayDirection.dot(rayDirection);
        float b = 2.0f * oc.dot(rayDirection);
        float c = oc.dot(oc) - radius * radius;
        float discriminant = b * b - 4 * a * c;

        if (discriminant < 0) {
            return false; // No intersection
        }

        // Find the closest intersection point
        float sqrtDiscriminant = std::sqrt(discriminant);
        float t1 = (-b - sqrtDiscriminant) / (2.0f * a);
        float t2 = (-b + sqrtDiscriminant) / (2.0f * a);

        t = (t1 < t2) ? t1 : t2;
        return true;
    }
};

// Simple ray tracing function
Vec3 traceRay(const Vec3& rayOrigin, const Vec3& rayDirection, const std::vector<Sphere>& spheres) {
    Vec3 color(0, 0, 0);
    float minDistance = std::numeric_limits<float>::max();

    // Check for intersection with each sphere
    for (const Sphere& sphere : spheres) {
        float t;
        if (sphere.intersect(rayOrigin, rayDirection, t) && t < minDistance) {
            minDistance = t;
            color = Vec3(1, 0, 0); // Set color to red for simplicity
        }
    }

    return color;
}

int main() {
    // Define the camera parameters
    Vec3 cameraPosition(0, 0, -5);
    Vec3 cameraDirection(0, 0, 1);
    Vec3 cameraUp(0, 1, 0);
    Vec3 cameraRight = cameraDirection.cross(cameraUp).normalized();
    Vec3 cameraUpCorrected = cameraRight.cross(cameraDirection).normalized();

    // Define the image resolution
    int imageWidth = 800;
    int imageHeight = 600;

    // Create a list of spheres in the scene
    std::vector<Sphere> spheres = {
        Sphere(Vec3(0, 0, 3), 1),
        Sphere(Vec3(2, 0, 4), 1),
        Sphere(Vec3(-2, 0, 4), 1),
    };

    // Simple ray tracing loop
    for (int y = 0; y < imageHeight; ++y) {
        for (int x = 0; x < imageWidth; ++x) {
            // Compute the ray direction for each pixel
            float u = static_cast<float>(x) / imageWidth;
            float v = static_cast<float>(y) / imageHeight;
            Vec3 rayDirection = (cameraDirection + (cameraRight * u) + (cameraUpCorrected * v)).normalized();

            // Trace the ray and get the color
            Vec3 pixelColor = traceRay(cameraPosition, rayDirection, spheres);

            // Output the pixel color (in this case, just print the color values)
            std::cout << pixelColor.x << " " << pixelColor.y << " " << pixelColor.z << " ";
        }
        std::cout << std::endl; // Move to the next line after each row
    }

    return 0;
}
