#include <iostream>
#include <functional>
#include <vector>
#include <cmath>
#include <stdexcept>
#include <iomanip>

class NumericalIntegrator {
private:
    // Function to validate input parameters
    static void validateInput(double a, double b, int n) {
        if (n <= 0) {
            throw std::invalid_argument("Number of intervals must be positive");
        }
        if (a >= b) {
            throw std::invalid_argument("Upper bound must be greater than lower bound");
        }
    }

public:
    // Rectangular method (Midpoint rule)
    static double rectangular(std::function<double(double)> f, double a, double b, int n) {
        validateInput(a, b, n);
        double h = (b - a) / n;
        double result = 0.0;
        
        for (int i = 0; i < n; i++) {
            double x_mid = a + (i + 0.5) * h;
            result += f(x_mid);
        }
        
        return h * result;
    }

    // Trapezoidal rule
    static double trapezoidal(std::function<double(double)> f, double a, double b, int n) {
        validateInput(a, b, n);
        double h = (b - a) / n;
        double result = (f(a) + f(b)) / 2.0;
        
        for (int i = 1; i < n; i++) {
            result += f(a + i * h);
        }
        
        return h * result;
    }

    // Simpson's 1/3 rule
    static double simpsons(std::function<double(double)> f, double a, double b, int n) {
        validateInput(a, b, n);
        if (n % 2 != 0) {
            throw std::invalid_argument("Number of intervals must be even for Simpson's rule");
        }

        double h = (b - a) / n;
        double result = f(a) + f(b);
        
        for (int i = 1; i < n; i++) {
            double x = a + i * h;
            result += (i % 2 == 0 ? 2 : 4) * f(x);
        }
        
        return h * result / 3.0;
    }

    // Simpson's 3/8 rule
    static double simpsons38(std::function<double(double)> f, double a, double b, int n) {
        validateInput(a, b, n);
        if (n % 3 != 0) {
            throw std::invalid_argument("Number of intervals must be divisible by 3 for Simpson's 3/8 rule");
        }

        double h = (b - a) / n;
        double result = f(a) + f(b);
        
        for (int i = 1; i < n; i++) {
            double x = a + i * h;
            result += ((i % 3 == 0) ? 2 : 3) * f(x);
        }
        
        return 3.0 * h * result / 8.0;
    }

    // Boole's rule
    static double booles(std::function<double(double)> f, double a, double b, int n) {
        validateInput(a, b, n);
        if (n % 4 != 0) {
            throw std::invalid_argument("Number of intervals must be divisible by 4 for Boole's rule");
        }

        double h = (b - a) / n;
        double result = 7 * (f(a) + f(b));
        
        for (int i = 1; i < n; i++) {
            double x = a + i * h;
            int coef;
            switch (i % 4) {
                case 0: coef = 14; break;
                case 1: case 3: coef = 32; break;
                case 2: coef = 12; break;
                default: coef = 0;
            }
            result += coef * f(x);
        }
        
        return 2.0 * h * result / 45.0;
    }

    // Romberg Integration
    static double romberg(std::function<double(double)> f, double a, double b, int maxOrder) {
        validateInput(a, b, maxOrder);
        std::vector<std::vector<double>> R(maxOrder + 1, std::vector<double>(maxOrder + 1));
        
        // Calculate R(0,0)
        R[0][0] = (b - a) * (f(a) + f(b)) / 2.0;
        
        // Calculate subsequent rows
        for (int i = 1; i <= maxOrder; i++) {
            // Calculate R(i,0)
            double h = (b - a) / std::pow(2, i);
            double sum = 0.0;
            for (int k = 1; k <= std::pow(2, i - 1); k++) {
                sum += f(a + (2 * k - 1) * h);
            }
            R[i][0] = R[i-1][0] / 2.0 + h * sum;
            
            // Calculate R(i,j)
            for (int j = 1; j <= i; j++) {
                double coef = std::pow(4, j);
                R[i][j] = (coef * R[i][j-1] - R[i-1][j-1]) / (coef - 1);
            }
        }
        
        return R[maxOrder][maxOrder];
    }

    // Adaptive quadrature using Simpson's rule
    static double adaptiveSimpson(
        std::function<double(double)> f,
        double a,
        double b,
        double tolerance = 1e-8,
        int maxDepth = 50
    ) {
        auto simpsonRule = [&](double a, double b) {
            double c = (a + b) / 2;
            return (b - a) / 6 * (f(a) + 4 * f(c) + f(b));
        };

        std::function<double(double, double, double, double, int)> adaptive;
        adaptive = [&](double a, double b, double fa, double fb, int depth) {
            double c = (a + b) / 2;
            double fc = f(c);
            double s = simpsonRule(a, b);
            double s1 = simpsonRule(a, c);
            double s2 = simpsonRule(c, b);
            double diff = std::abs(s1 + s2 - s);

            if (depth > maxDepth) {
                return s1 + s2;
            }

            if (diff < 15 * tolerance) {
                return s1 + s2 + diff / 15;
            }

            return adaptive(a, c, fa, fc, depth + 1) +
                   adaptive(c, b, fc, fb, depth + 1);
        };

        return adaptive(a, b, f(a), f(b), 0);
    }

    // Error estimation
    static double estimateError(
        std::function<double(double)> f,
        double a,
        double b,
        double actualValue,
        const std::string& method,
        int n
    ) {
        double numericalResult;
        
        if (method == "trapezoidal") {
            numericalResult = trapezoidal(f, a, b, n);
        } else if (method == "simpsons") {
            numericalResult = simpsons(f, a, b, n);
        } else if (method == "rectangular") {
            numericalResult = rectangular(f, a, b, n);
        } else {
            throw std::invalid_argument("Unknown method for error estimation");
        }
        
        return std::abs(actualValue - numericalResult);
    }
};

// Example usage
int main() {
    // Test function: f(x) = x^2
    auto f = [](double x) { return x * x; };
    double a = 0.0;
    double b = 1.0;
    int n = 36;

    // Known analytical result for integral of x^2 from 0 to 1 is 1/3
    double actual = 1.0 / 3.0;

    try {
        std::cout << std::fixed << std::setprecision(10);
        std::cout << "Integrating f(x) = x^2 from " << a << " to " << b << "\n\n";

        // Compare different methods
        std::cout << "Rectangular method: " 
                  << NumericalIntegrator::rectangular(f, a, b, n) << "\n";
        
        std::cout << "Trapezoidal rule:  " 
                  << NumericalIntegrator::trapezoidal(f, a, b, n) << "\n";
        
        std::cout << "Simpson's rule:    " 
                  << NumericalIntegrator::simpsons(f, a, b, n) << "\n";
        
        std::cout << "Simpson's 3/8:     " 
                  << NumericalIntegrator::simpsons38(f, a, b, n) << "\n";
        
        std::cout << "Boole's rule:      " 
                  << NumericalIntegrator::booles(f, a, b, n) << "\n";
        
        std::cout << "Romberg (order 4): " 
                  << NumericalIntegrator::romberg(f, a, b, 4) << "\n";
        
        std::cout << "Adaptive Simpson:  " 
                  << NumericalIntegrator::adaptiveSimpson(f, a, b) << "\n";
        
        std::cout << "\nActual value:      " << actual << "\n";

        // Error estimation
        std::cout << "\nError estimation:\n";
        std::cout << "Trapezoidal error: " 
                  << NumericalIntegrator::estimateError(f, a, b, actual, "trapezoidal", n) << "\n";
        std::cout << "Simpson's error:   " 
                  << NumericalIntegrator::estimateError(f, a, b, actual, "simpsons", n) << "\n";

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }

    return 0;
}