// #include <iostream>
// #include <iomanip> // For controlling decimal precision
// #include <vector>

// using namespace std;

// // Define the first-order ODE: dy/dx = f(x, y)
// double f(double x, double y) {
//     return x + y; // Example: dy/dx = x + y
// }

// // Runge-Kutta 4th Order Method (RK4)
// vector<pair<double, double>> rungeKutta4(double x0, double y0, double xEnd, double h) {
//     vector<pair<double, double>> results; // Store (x, y) pairs
//     results.emplace_back(x0, y0);

//     double x = x0;
//     double y = y0;

//     while (x < xEnd) {
//         // Ensure the final step does not overshoot xEnd
//         if (x + h > xEnd) h = xEnd - x;

//         // Calculate the Runge-Kutta coefficients
//         double k1 = h * f(x, y);
//         double k2 = h * f(x + h / 2.0, y + k1 / 2.0);
//         double k3 = h * f(x + h / 2.0, y + k2 / 2.0);
//         double k4 = h * f(x + h, y + k3);

//         // Update the value of y using RK4 formula
//         y += (k1 + 2 * k2 + 2 * k3 + k4) / 6.0;
        
//         // Update x to the next step
//         x += h;

//         // Store the current (x, y) pair
//         results.emplace_back(x, y);
//     }

//     return results;
// }

// int main() {
//     // Initial conditions
//     double x0 = 0.0; // Initial value of x
//     double y0 = 1.0; // Initial value of y

//     // Range of x
//     double xEnd = 2.0; // Solve from x0 to xEnd

//     // Step size
//     double h = 0.1;

//     // Solve the differential equation using RK4
//     vector<pair<double, double>> results = rungeKutta4(x0, y0, xEnd, h);

//     // Print the results
//     cout << fixed << setprecision(6); // Set precision for output
//     cout << "x\t\ty" << endl;
//     for (const auto& [x, y] : results) {
//         cout << x << "\t" << y << endl;
//     }

//     return 0;
// }


#include <iostream>
#include <vector>
#include <functional>
#include <cmath>

class DifferentialSolver {
private:
    // Function type for differential equations: dy/dx = f(x,y)
    using DiffFunction = std::function<double(double, double)>;
    
    // RK4 method implementation
    static double rk4Step(DiffFunction f, double x, double y, double h) {
        double k1 = f(x, y);
        double k2 = f(x + h/2, y + h*k1/2);
        double k3 = f(x + h/2, y + h*k2/2);
        double k4 = f(x + h, y + h*k3);
        
        return y + (h/6) * (k1 + 2*k2 + 2*k3 + k4);
    }

public:
    // Solve the differential equation using RK4 method
    static std::vector<std::pair<double, double>> solve(
        DiffFunction f,     // The differential equation dy/dx = f(x,y)
        double x0,          // Initial x value
        double y0,          // Initial y value
        double xEnd,        // Final x value
        double stepSize     // Step size (h)
    ) {
        std::vector<std::pair<double, double>> solution;
        
        double x = x0;
        double y = y0;
        
        // Store initial point
        solution.push_back({x, y});
        
        // Solve until we reach xEnd
        while (x < xEnd) {
            // Adjust step size if we would overshoot xEnd
            double h = std::min(stepSize, xEnd - x);
            
            // Calculate next y using RK4
            y = rk4Step(f, x, y, h);
            x += h;
            
            // Store point
            solution.push_back({x, y});
        }
        
        return solution;
    }
};

// Example usage
int main() {
    // Example 1: dy/dx = x + y
    auto diff_eq1 = [](double x, double y) { return x + y; };
    
    // Solve from x=0 to x=1 with initial condition y(0)=1
    auto solution1 = DifferentialSolver::solve(diff_eq1, 0.0, 1.0, 1.0, 0.1);
    
    std::cout << "Solution for dy/dx = x + y, y(0) = 1:\n";
    std::cout << "x\t\ty\n";
    for (const auto& point : solution1) {
        std::cout << point.first << "\t\t" << point.second << "\n";
    }
    
    // Example 2: dy/dx = -y (exponential decay)
    auto diff_eq2 = [](double x, double y) { return -y; };
    
    // Solve from x=0 to x=2 with initial condition y(0)=1
    auto solution2 = DifferentialSolver::solve(diff_eq2, 0.0, 1.0, 2.0, 0.1);
    
    std::cout << "\nSolution for dy/dx = -y, y(0) = 1:\n";
    std::cout << "x\t\ty\n";
    for (const auto& point : solution2) {
        std::cout << point.first << "\t\t" << point.second << "\n";
    }
    
    return 0;
}