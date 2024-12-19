// // // #include <SFML/Graphics.hpp>
// // // #include <complex>
// // // #include <vector>
// // // #include <thread>
// // // #include <mutex>

// // // class MandelbrotRenderer {
// // // private:
// // //     static const int MAX_THREADS = 8;
// // //     std::mutex imageMutex;
// // //     sf::Image image;
// // //     sf::Texture texture;
// // //     sf::Sprite sprite;
    
// // //     // View parameters
// // //     double centerX = -0.5;
// // //     double centerY = 0.0;
// // //     double zoom = 1.0;
// // //     int colorShift = 0;
    
// // //     // Constants
// // //     const int width;
// // //     const int height;
// // //     const int maxIterations;

// // //     int mandelbrot(std::complex<double> c) const {
// // //         std::complex<double> z = 0;
// // //         int iterations = 0;
        
// // //         // Optimization: Skip computation for main cardioid and period-2 bulb
// // //         double q = std::pow(c.real() - 0.25, 2) + std::pow(c.imag(), 2);
// // //         if (q * (q + (c.real() - 0.25)) <= 0.25 * std::pow(c.imag(), 2) ||
// // //             std::pow(c.real() + 1.0, 2) + std::pow(c.imag(), 2) <= 0.0625) {
// // //             return maxIterations;
// // //         }

// // //         while (std::abs(z) <= 2.0 && iterations < maxIterations) {
// // //             z = z * z + c;
// // //             ++iterations;
// // //         }
// // //         return iterations;
// // //     }

// // //     sf::Color getColor(int iterations) const {
// // //         if (iterations == maxIterations) return sf::Color::Black;

// // //         // Smooth coloring algorithm
// // //         double t = (double)iterations / maxIterations;
// // //         int hue = (int)(360 * t + colorShift) % 360;
        
// // //         // Convert HSV to RGB
// // //         double h = hue / 60.0;
// // //         double x = 1.0 - std::abs(std::fmod(h, 2.0) - 1.0);
        
// // //         double r = 0, g = 0, b = 0;
// // //         if (h < 1) { r = 1; g = x; }
// // //         else if (h < 2) { r = x; g = 1; }
// // //         else if (h < 3) { g = 1; b = x; }
// // //         else if (h < 4) { g = x; b = 1; }
// // //         else if (h < 5) { r = x; b = 1; }
// // //         else { r = 1; b = x; }

// // //         return sf::Color(
// // //             static_cast<sf::Uint8>(r * 255),
// // //             static_cast<sf::Uint8>(g * 255),
// // //             static_cast<sf::Uint8>(b * 255)
// // //         );
// // //     }

// // //     void renderSegment(int startY, int endY) {
// // //         double aspectRatio = (double)width / height;
// // //         double scale = 2.0 / (zoom * height);

// // //         for (int y = startY; y < endY; ++y) {
// // //             for (int x = 0; x < width; ++x) {
// // //                 double real = (x - width/2) * scale * aspectRatio + centerX;
// // //                 double imag = (y - height/2) * scale + centerY;
                
// // //                 int iterations = mandelbrot(std::complex<double>(real, imag));
// // //                 sf::Color color = getColor(iterations);

// // //                 std::lock_guard<std::mutex> lock(imageMutex);
// // //                 image.setPixel(x, y, color);
// // //             }
// // //         }
// // //     }

// // // public:
// // //     MandelbrotRenderer(int w, int h, int maxIter) 
// // //         : width(w), height(h), maxIterations(maxIter) {
// // //         image.create(width, height);
// // //         texture.create(width, height);
// // //     }

// // //     void render() {
// // //         std::vector<std::thread> threads;
// // //         int segmentHeight = height / MAX_THREADS;

// // //         for (int i = 0; i < MAX_THREADS; ++i) {
// // //             int startY = i * segmentHeight;
// // //             int endY = (i == MAX_THREADS-1) ? height : (i+1) * segmentHeight;
// // //             threads.emplace_back(&MandelbrotRenderer::renderSegment, this, startY, endY);
// // //         }

// // //         for (auto& thread : threads) {
// // //             thread.join();
// // //         }

// // //         texture.loadFromImage(image);
// // //         sprite.setTexture(texture, true);
// // //     }

// // //     void handleZoom(sf::Vector2i mousePos, bool zoomIn) {
// // //         double zoomFactor = zoomIn ? 1.5 : 0.67;
        
// // //         // Convert mouse position to complex plane coordinates before zoom
// // //         double aspectRatio = (double)width / height;
// // //         double scale = 2.0 / (zoom * height);
// // //         double real = (mousePos.x - width/2) * scale * aspectRatio + centerX;
// // //         double imag = (mousePos.y - height/2) * scale + centerY;
        
// // //         zoom *= zoomFactor;
        
// // //         // Adjust center to keep mouse position fixed
// // //         centerX = real - (mousePos.x - width/2) * scale * aspectRatio / zoomFactor;
// // //         centerY = imag - (mousePos.y - height/2) * scale / zoomFactor;
// // //     }

// // //     void shiftColors() {
// // //         colorShift = (colorShift + 30) % 360;
// // //     }

// // //     const sf::Sprite& getSprite() const { return sprite; }
// // // };

// // // int main() {
// // //     const int width = 800;
// // //     const int height = 800;
// // //     const int maxIterations = 1000;

// // //     sf::RenderWindow window(sf::VideoMode(width, height), "Enhanced Mandelbrot Set");
// // //     window.setFramerateLimit(60);

// // //     MandelbrotRenderer renderer(width, height, maxIterations);
// // //     renderer.render();

// // //     while (window.isOpen()) {
// // //         sf::Event event;
// // //         bool needsUpdate = false;

// // //         while (window.pollEvent(event)) {
// // //             if (event.type == sf::Event::Closed) {
// // //                 window.close();
// // //             }
// // //             else if (event.type == sf::Event::MouseButtonPressed) {
// // //                 if (event.mouseButton.button == sf::Mouse::Left) {
// // //                     renderer.handleZoom(sf::Mouse::getPosition(window), true);
// // //                     needsUpdate = true;
// // //                 }
// // //                 else if (event.mouseButton.button == sf::Mouse::Right) {
// // //                     renderer.handleZoom(sf::Mouse::getPosition(window), false);
// // //                     needsUpdate = true;
// // //                 }
// // //             }
// // //             else if (event.type == sf::Event::KeyPressed) {
// // //                 if (event.key.code == sf::Keyboard::Space) {
// // //                     renderer.shiftColors();
// // //                     needsUpdate = true;
// // //                 }
// // //             }
// // //         }

// // //         if (needsUpdate) {
// // //             renderer.render();
// // //         }

// // //         window.clear();
// // //         window.draw(renderer.getSprite());
// // //         window.display();
// // //     }

// // //     return 0;
// // // }


// // // #include <iostream>
// // // #include <complex>
// // // #include <string>
// // // #include <vector>
// // // #include <chrono>
// // // #include <thread>

// // // class ASCIIMandelbrot {
// // // private:
// // //     const int width;
// // //     const int height;
// // //     const int maxIterations;
// // //     const std::string densityMap = " .'`^\",:;Il!i><~+_-?][}{1)(|\\/tfjrxnuvczXYUJCLQ0OZmwqpdbkhao*#MW&8%B@$";

// // //     int mandelbrot(std::complex<double> c) const {
// // //         std::complex<double> z = 0;
// // //         int iterations = 0;
        
// // //         // Optimization: Skip computation for main cardioid and period-2 bulb
// // //         double q = std::pow(c.real() - 0.25, 2) + std::pow(c.imag(), 2);
// // //         if (q * (q + (c.real() - 0.25)) <= 0.25 * std::pow(c.imag(), 2) ||
// // //             std::pow(c.real() + 1.0, 2) + std::pow(c.imag(), 2) <= 0.0625) {
// // //             return maxIterations;
// // //         }

// // //         while (std::abs(z) <= 2.0 && iterations < maxIterations) {
// // //             z = z * z + c;
// // //             ++iterations;
// // //         }
// // //         return iterations;
// // //     }

// // //     char getASCIIChar(int iterations) const {
// // //         if (iterations == maxIterations) return densityMap.back();
        
// // //         // Map iterations to density character
// // //         double t = static_cast<double>(iterations) / maxIterations;
// // //         int index = static_cast<int>(t * (densityMap.length() - 1));
// // //         return densityMap[index];
// // //     }

// // // public:
// // //     ASCIIMandelbrot(int w, int h, int maxIter) 
// // //         : width(w), height(h), maxIterations(maxIter) {}

// // //     void render(double centerX = -0.5, double centerY = 0.0, double zoom = 1.0) const {
// // //         // Adjust for console character aspect ratio (characters are taller than wide)
// // //         double aspectRatio = 2.5; // Typical console font aspect ratio
// // //         double scale = 2.0 / (zoom * height);

// // //         // Clear screen (ANSI escape code)
// // //         std::cout << "\033[2J\033[H";
        
// // //         for (int y = 0; y < height; ++y) {
// // //             for (int x = 0; x < width; ++x) {
// // //                 double real = (x - width/2) * scale * aspectRatio + centerX;
// // //                 double imag = (y - height/2) * scale + centerY;
                
// // //                 int iterations = mandelbrot(std::complex<double>(real, imag));
// // //                 char c = getASCIIChar(iterations);
                
// // //                 // Print character twice to make the aspect ratio more square-like
// // //                 std::cout << c << c;
// // //             }
// // //             std::cout << '\n';
// // //         }
// // //     }

// // //     void animate() const {
// // //         double zoom = 1.0;
// // //         double centerX = -0.5;
// // //         double centerY = 0.0;
        
// // //         while (true) {
// // //             render(centerX, centerY, zoom);
// // //             zoom *= 1.1;
            
// // //             // Sleep for a short duration
// // //             std::this_thread::sleep_for(std::chrono::milliseconds(100));
            
// // //             // Break animation after certain zoom level
// // //             if (zoom > 50.0) break;
// // //         }
// // //     }
// // // };

// // // // ANSI escape codes for colors
// // // struct ANSIColor {
// // //     static std::string reset() { return "\033[0m"; }
// // //     static std::string fg(int code) { return "\033[38;5;" + std::to_string(code) + "m"; }
// // //     static std::string bg(int code) { return "\033[48;5;" + std::to_string(code) + "m"; }
// // // };

// // // int main() {
// // //     // Get terminal size (you might want to adjust these values based on your terminal)
// // //     const int width = 40;
// // //     const int height = 20;
// // //     const int maxIterations = 100;

// // //     ASCIIMandelbrot mandelbrot(width, height, maxIterations);
    
// // //     // Option 1: Static render
// // //     //mandelbrot.render();
    
// // //     // Option 2: Animated zoom
// // //      mandelbrot.animate();
    
// // //     return 0;
// // // }

// // #include<iostream>
// // #include<string>
// // #include<thread>
// // #include<map>
// // #include <cmath>
// // #include <cstring>
// // using namespace std;
// // int k;
// // int main()
// // {
// // 	float A = 0;
// // 	float B = 0;
// // 	float i;
// // 	float j;
// // 	float z[1760];
// // 	char b[1760];
// // 	cout<<("\x1b[2J\033[34m");
// // 	for (int x=0; x<=1000; x++)
// // 	{
// // 		memset(b, 32, 1760);
// // 		memset(z, 0, 7040);
// // 		for (j = 0; 6.28 > j; j += 0.07)
// // 			for (i = 0; 6.28 > i; i += 0.02)
// // 			{	float c = sin(i);
// // 				float d = cos(j);
// // 				float e = sin(A);
// // 				float f = sin(j);
// // 				float g = cos(A);
// // 				float h = d +  2;
// // 				float D = 1 / (c * h * e + f * g + 5);
// // 				float l = cos(i);
// // 				float m = cos(B);
// // 				float n = sin(B);
// // 				float t = c * h * g - f * e;

// // 				int x = 40 + 30 * D * (l * h * m - t * n);
// // 				int y = 12 + 15 * D * (l * h * n + t * m);
// // 				int o = x + 80 * y;
// // 				int N = 8 * ((f * e - c * d * g) * m - c * d * e - f * g - l * d * n);

// // 				if (22 > y && y > 0 && x > 0 && 80 > x && D > z[o])
// // 				{	z[o] = D;
// // 					;
// // 					;
// // 					b[o] =".,-~:;-+/=%?!^*#$@"[N > 0 ? N : 0]; } } 

// // 		cout<<("\x1b[H");
// // 		for (k = 0; 1761 > k; k++)
// // 			putchar(k % 80 ? b[k] : 10);cout.flush();
// // 		A += 0.05504;
// // 		B +=  0.0552;

		  
// // 	}  }

// #include <iostream>
// #include <string>
// #include <thread>
// #include <chrono>
// #include <cmath>
// #include <vector>
// #include <array>

// class DonutRenderer {
// private:
//     static constexpr int SCREEN_WIDTH = 80;
//     static constexpr int SCREEN_HEIGHT = 24;
//     static constexpr double THETA_SPACING = 0.07;
//     static constexpr double PHI_SPACING = 0.02;
//     static constexpr int R1 = 1;  // Radius of the tube
//     static constexpr int R2 = 2;  // Distance from center to tube
//     static constexpr double K2 = 5;  // Distance from viewer to donut
//     static constexpr double K1 = SCREEN_HEIGHT * K2 * 3 / (8 * (R1 + R2));  // Scale factor

//     // Prettier ASCII density map
//     const std::string densityMap = " .,-~:;=!*#$@";
    
//     struct Point3D {
//         double x, y, z;
//     };

//     double A = 0;  // Rotation about the X-axis
//     double B = 0;  // Rotation about the Z-axis

//     std::vector<char> output;
//     std::vector<double> zBuffer;

// public:
//     DonutRenderer() : 
//         output(SCREEN_WIDTH * SCREEN_HEIGHT, ' '),
//         zBuffer(SCREEN_WIDTH * SCREEN_HEIGHT, 0) {}

//     void render() {
//         // Clear buffers
//         std::fill(output.begin(), output.end(), ' ');
//         std::fill(zBuffer.begin(), zBuffer.end(), 0);

//         // Precompute trigonometric values
//         double cosA = cos(A), sinA = sin(A);
//         double cosB = cos(B), sinB = sin(B);

//         // Render the donut
//         for (double theta = 0; theta < 2 * M_PI; theta += THETA_SPACING) {
//             double cosTheta = cos(theta), sinTheta = sin(theta);

//             for (double phi = 0; phi < 2 * M_PI; phi += PHI_SPACING) {
//                 double cosPhi = cos(phi), sinPhi = sin(phi);

//                 // 3D coordinates of the point on the torus
//                 Point3D p = {
//                     (R2 + R1 * cosTheta) * (cosPhi * cosB + sinA * sinPhi * sinB) - R1 * sinTheta * cosA * sinB,
//                     (R2 + R1 * cosTheta) * (cosPhi * sinB - sinA * sinPhi * cosB) + R1 * sinTheta * cosA * cosB,
//                     K2 + cosA * (R2 + R1 * cosTheta) * sinPhi + R1 * sinTheta * sinA
//                 };

//                 // Project 3D point to 2D screen coordinates
//                 int x = static_cast<int>(SCREEN_WIDTH/2 + K1 * p.x / p.z);
//                 int y = static_cast<int>(SCREEN_HEIGHT/2 - K1 * p.y / p.z);

//                 // Calculate luminance
//                 double L = cosPhi * cosTheta * sinB - cosA * cosTheta * sinPhi - 
//                           sinA * sinTheta + sinB * (cosA * sinTheta - cosTheta * sinA * sinPhi);

//                 // Only render if point is visible and in front of previous points
//                 if (y >= 0 && y < SCREEN_HEIGHT && x >= 0 && x < SCREEN_WIDTH && p.z > zBuffer[y * SCREEN_WIDTH + x]) {
//                     zBuffer[y * SCREEN_WIDTH + x] = p.z;
//                     int luminanceIndex = static_cast<int>((L + 1) * (densityMap.length() - 1) / 2);
//                     luminanceIndex = std::max(0, std::min(static_cast<int>(densityMap.length() - 1), luminanceIndex));
//                     output[y * SCREEN_WIDTH + x] = densityMap[luminanceIndex];
//                 }
//             }
//         }

//         // Move cursor to top-left corner and clear screen
//         std::cout << "\x1b[H\x1b[2J";
        
//         // Draw the frame
//         for (int y = 0; y < SCREEN_HEIGHT; y++) {
//             for (int x = 0; x < SCREEN_WIDTH; x++) {
//                 std::cout << output[y * SCREEN_WIDTH + x];
//             }
//             std::cout << '\n';
//         }

//         // Update rotation angles
//         A += 0.07;
//         B += 0.03;
//     }

//     void animate(int frames = -1) {
//         int frameCount = 0;
//         while (frames == -1 || frameCount < frames) {
//             render();
//             std::this_thread::sleep_for(std::chrono::milliseconds(50));
//             frameCount++;
//         }
//     }

//     // Add color support
//     void enableColor(bool rainbow = false) {
//         densityMap = rainbow ? 
//             "\x1b[31m.\x1b[32m,\x1b[33m-\x1b[34m~\x1b[35m:\x1b[36m;\x1b[37m=\x1b[91m!\x1b[92m*\x1b[93m#\x1b[94m$\x1b[95m@\x1b[0m" :
//             "\x1b[36m .-~=+*#%@\x1b[0m";
//     }
// };

// int main() {
//     DonutRenderer donut;
    
//     // Optional: Enable color output
//     // donut.enableColor(true);  // Set to true for rainbow mode
    
//     // Start the animation
//     donut.animate();
    
//     return 0;
// }


#include <iostream>
#include <complex>
#include <string>
#include <vector>

class FractalGenerator {
private:
    static const int MAX_ITERATIONS = 100;
    int width;
    int height;
    std::vector<std::vector<char>> canvas;

    // Function to map iteration count to ASCII characters
    char mapToAscii(int iterations) {
        const std::string charset = " .:-=+*#%@";
        if (iterations == MAX_ITERATIONS) return ' ';
        return charset[iterations % charset.length()];
    }

    // Mandelbrot set calculation
    int mandelbrotIterations(std::complex<double> c) {
        std::complex<double> z(0, 0);
        int iterations = 0;
        
        while (abs(z) <= 2.0 && iterations < MAX_ITERATIONS) {
            z = z * z + c;
            iterations++;
        }
        
        return iterations;
    }

    // Sierpinski triangle recursive generation
    void drawSierpinskiTriangle(int x, int y, int size, int level) {
        if (level == 0) {
            if (x >= 0 && x < width && y >= 0 && y < height) {
                canvas[y][x] = '*';
            }
            return;
        }

        int newSize = size / 2;
        drawSierpinskiTriangle(x, y - newSize, newSize, level - 1);
        drawSierpinskiTriangle(x - newSize, y + newSize, newSize, level - 1);
        drawSierpinskiTriangle(x + newSize, y + newSize, newSize, level - 1);
    }

public:
    FractalGenerator(int w, int h) : width(w), height(h) {
        canvas.resize(height, std::vector<char>(width, ' '));
    }

    // Generate Mandelbrot set
    void generateMandelbrot() {
        for (int y = 0; y < height; y++) {
            for (int x = 0; x < width; x++) {
                double real = (x - width/2.0) * 4.0/width;
                double imag = (y - height/2.0) * 4.0/height;
                std::complex<double> c(real, imag);
                
                int iterations = mandelbrotIterations(c);
                canvas[y][x] = mapToAscii(iterations);
            }
        }
    }

    // Generate Sierpinski Triangle
    void generateSierpinski(int levels) {
        // Clear canvas
        for (auto& row : canvas) {
            std::fill(row.begin(), row.end(), ' ');
        }
        
        int size = std::min(width, height) / 2;
        int centerX = width / 2;
        int centerY = height / 2;
        
        drawSierpinskiTriangle(centerX, centerY, size, levels);
    }

    // Custom recursive tree fractal
    void drawTree(int x, int y, int length, double angle, int depth) {
        if (depth == 0) return;

        int endX = x + int(length * cos(angle));
        int endY = y + int(length * sin(angle));

        // Draw line
        drawLine(x, y, endX, endY);

        // Recursive branches
        drawTree(endX, endY, length * 0.7, angle - 0.5, depth - 1);
        drawTree(endX, endY, length * 0.7, angle + 0.5, depth - 1);
    }

    // Helper function to draw a line using Bresenham's algorithm
    void drawLine(int x1, int y1, int x2, int y2) {
        int dx = abs(x2 - x1);
        int dy = abs(y2 - y1);
        int sx = (x1 < x2) ? 1 : -1;
        int sy = (y1 < y2) ? 1 : -1;
        int err = dx - dy;

        while (true) {
            if (x1 >= 0 && x1 < width && y1 >= 0 && y1 < height) {
                canvas[y1][x1] = '*';
            }

            if (x1 == x2 && y1 == y2) break;
            int e2 = 2 * err;
            if (e2 > -dy) { err -= dy; x1 += sx; }
            if (e2 < dx) { err += dx; y1 += sy; }
        }
    }

    // Generate fractal tree
    void generateTree() {
        // Clear canvas
        for (auto& row : canvas) {
            std::fill(row.begin(), row.end(), ' ');
        }
        
        drawTree(width/2, height-1, height/3, -M_PI/2, 9);
    }

    // Display the fractal
    void display() {
        for (const auto& row : canvas) {
            for (char c : row) {
                std::cout << c << c; // Double each character for better aspect ratio
            }
            std::cout << '\n';
        }
    }
};

int main() {
    FractalGenerator mandelbrot(60, 30);
    std::cout << "Mandelbrot Set:\n\n";
    mandelbrot.generateMandelbrot();
    mandelbrot.display();

    std::cout << "\nSierpinski Triangle:\n\n";
    FractalGenerator sierpinski(60, 30);
    sierpinski.generateSierpinski(6);
    sierpinski.display();

    std::cout << "\nFractal Tree:\n\n";
    FractalGenerator tree(60, 30);
    tree.generateTree();
    tree.display();

    return 0;
}