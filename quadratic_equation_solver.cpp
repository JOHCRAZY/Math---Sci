#include <iostream>
#include <cmath>
#include <complex>
#include <iomanip>
#include <limits>
#include <vector>

class QuadraticSolver {
public:
    // Enum to represent different types of roots
    enum class RootType {
        REAL_DISTINCT,
        REAL_EQUAL,
        COMPLEX
    };

    // Structure to hold solving results
    struct SolutionResult {
        RootType type;
        std::complex<double> root1;
        std::complex<double> root2;
    };

    // Solve quadratic equation ax^2 + bx + c = 0
    static SolutionResult solve(double a, double b, double c) {
        // Check for invalid input (a cannot be zero)
        if (std::abs(a) < std::numeric_limits<double>::epsilon()) {
            throw std::invalid_argument("Coefficient 'a' cannot be zero for a quadratic equation.");
        }

        // Calculate discriminant
        double discriminant = b*b - 4*a*c;

        SolutionResult result;

        // Determine root type based on discriminant
        if (discriminant > 0) {
            // Two distinct real roots
            result.type = RootType::REAL_DISTINCT;
            double sqrt_disc = std::sqrt(discriminant);
            result.root1 = std::complex<double>((-b + sqrt_disc) / (2*a), 0);
            result.root2 = std::complex<double>((-b - sqrt_disc) / (2*a), 0);
        }
        else if (std::abs(discriminant) < std::numeric_limits<double>::epsilon()) {
            // One real root (repeated)
            result.type = RootType::REAL_EQUAL;
            result.root1 = result.root2 = std::complex<double>(-b / (2*a), 0);
        }
        else {
            // Complex roots
            result.type = RootType::COMPLEX;
            double real_part = -b / (2*a);
            double imag_part = std::sqrt(-discriminant) / (2*a);
            result.root1 = std::complex<double>(real_part, imag_part);
            result.root2 = std::complex<double>(real_part, -imag_part);
        }

        return result;
    }

    // Pretty print the solution
    static void printSolution(const SolutionResult& solution) {
        std::cout << std::fixed << std::setprecision(4);

        switch (solution.type) {
            case RootType::REAL_DISTINCT:
                std::cout << "Two distinct real roots:" << std::endl;
                std::cout << "x1 = " << solution.root1.real() << std::endl;
                std::cout << "x2 = " << solution.root2.real() << std::endl;
                break;

            case RootType::REAL_EQUAL:
                std::cout << "One real root (repeated):" << std::endl;
                std::cout << "x = " << solution.root1.real() << std::endl;
                break;

            case RootType::COMPLEX:
                std::cout << "Two complex roots:" << std::endl;
                std::cout << "x1 = " << solution.root1.real() << " + " 
                          << solution.root1.imag() << "i" << std::endl;
                std::cout << "x2 = " << solution.root2.real() << " - " 
                          << std::abs(solution.root2.imag()) << "i" << std::endl;
                break;
        }
    }
};

int main() {
    try {
        // Test cases covering different scenarios
        std::vector<std::pair<double, std::pair<double, double>>> test_cases = {
            // a, b, c
            {1, {-3, 2}},    // Two distinct real roots
            {1, {-2, 1}},    // One real (repeated) root
            {1, {0, 1}}      // Complex roots
        };

        for (const auto& test : test_cases) {
            double a = test.first;
            double b = test.second.first;
            double c = test.second.second;

            std::cout << "\nSolving: " 
                      << a << "x² + " 
                      << b << "x + " 
                      << c << " = 0" << std::endl;

            // Solve and print
            auto solution = QuadraticSolver::solve(a, b, c);
            QuadraticSolver::printSolution(solution);
        }

        // Interactive mode
        double a, b, c;
        std::cout << "\nEnter coefficients (a b c): ";
        std::cin >> a >> b >> c;

        std::cout << "\nSolving: " 
                      << a << "x² + " 
                      << b << "x + " 
                      << c << " = 0" << std::endl;

        auto solution = QuadraticSolver::solve(a, b, c);
        QuadraticSolver::printSolution(solution);

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}