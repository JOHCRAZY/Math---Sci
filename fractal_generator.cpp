#include <iostream>
#include <complex>
#include <string>
#include <vector>
#include <algorithm>

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
    bool loop = true;
    int width = 90;
    int height = 90;
        
    while (loop) {
        std::cout << "\n=== Fractal Generator  ===\n"
                  << "1. Mandelbrot Set\n"
                  << "2. Sierpinski Triangle\n"
                  << "3. Fractal Tree\n"
                  << "4. Set width and Height\n"
                  << "5. Exit\n"
                  << "Enter your choice (1-5): ";

        int choice;
        std::cin >> choice;

        if (std::cin.fail()) {
            std::cin.clear();
            std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            std::cout << "Invalid input. Please enter a number.\n";
            continue;
        }

        
        switch (choice) {
            case 1: {
                FractalGenerator mandelbrot(width, height);
                mandelbrot.generateMandelbrot();
                
                std::cout << "\nMandelbrot Set:\n\n";
                mandelbrot.display();
                break;
            }
            
            case 2: {

                FractalGenerator sierpinski(width, height);
                sierpinski.generateSierpinski(6);
                
                std::cout << "\nSierpinski Triangle:\n\n";
                sierpinski.display();
                break;
            }
            
            case 3: {
                FractalGenerator tree(width, height);
                tree.generateTree();
                
                std::cout << "\nFractal Tree:\n\n";
                tree.display();
                break;
            }
            
            case 4: {
                std::cout << "\nChange Dimensions:\n"
                          << "Enter width (20 - 90): ";
                std::cin >> width;
                width = std::clamp(width, 20, 90);
                
                std::cout << "Enter height (20 - 90): ";
                std::cin >> height;
                height = std::clamp(height, 20, 90);
                
                std::cout << "\nNew dimensions set to " << width << "x" << height << "\n";
                break;
            }
            
            case 5:
                std::cout << "\nThank you for using the Fractal Generator!\n";
                loop = false;
                break;
                
            default:
                std::cout << "Invalid choice. Please select 1-5.\n";
        }
    }
    return 0;
}