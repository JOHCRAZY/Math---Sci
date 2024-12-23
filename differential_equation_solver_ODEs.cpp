#include <iostream>
#include <vector>
#include <functional>
#include <cmath>
#include <string>
#include <map>
#include <limits>
#include <iomanip>
#include <muParser.h>
#include <memory>


/*
use: 
     g++ differential_equation_solver_ODEs.cpp -o differential_equation_solver_ODEs -lmuparser
*/
class EquationParser {
private:
    std::unique_ptr<mu::Parser> parser;
    double x_val, y_val;

public:
    EquationParser() : x_val(0), y_val(0) {
        parser = std::make_unique<mu::Parser>();
        
        // Define variables
        parser->DefineVar("x", &x_val);
        parser->DefineVar("y", &y_val);
        
        // Define constants
        parser->DefineConst("pi", M_PI);
        parser->DefineConst("e", M_E);
        
        // Define common functions
        parser->DefineFun("sin", sin);
        parser->DefineFun("cos", cos);
        parser->DefineFun("tan", tan);
        parser->DefineFun("exp", exp);
        parser->DefineFun("log", log);
        parser->DefineFun("sqrt", sqrt);
    }

    std::function<double(double, double)> parseEquation(const std::string& eqStr) {
        try {
            parser->SetExpr(eqStr);
            
            return [this](double x, double y) {
                this->x_val = x;
                this->y_val = y;
                return parser->Eval();
            };
        }
        catch (mu::Parser::exception_type& e) {
            std::cerr << "Parser error: " << e.GetMsg() << std::endl;
            throw;
        }
    }

    static std::string getAvailableFunctions() {
        return "Available functions:\n"
               "sin(x), cos(x), tan(x), exp(x), log(x), sqrt(x)\n"
               "Constants: pi, e\n"
               "Operators: +, -, *, /, ^, (, )\n";
    }
};

class DifferentialSolver {
private:
    using DiffFunction = std::function<double(double, double)>;
    
    static double rk4Step(DiffFunction f, double x, double y, double h) {
        double k1 = f(x, y);
        double k2 = f(x + h/2, y + h*k1/2);
        double k3 = f(x + h/2, y + h*k2/2);
        double k4 = f(x + h, y + h*k3);
        
        return y + (h/6) * (k1 + 2*k2 + 2*k3 + k4);
    }

public:
    static std::vector<std::pair<double, double>> solve(
        DiffFunction f,
        double x0,
        double y0,
        double xEnd,
        double stepSize
    ) {
        std::vector<std::pair<double, double>> solution;
        double x = x0;
        double y = y0;
        solution.push_back({x, y});
        
        while (x < xEnd) {
            double h = std::min(stepSize, xEnd - x);
            y = rk4Step(f, x, y, h);
            x += h;
            solution.push_back({x, y});
        }
        
        return solution;
    }
};

class UserInterface {
private:
    EquationParser parser;
    std::map<std::string, std::string> predefinedEquations = {
        {"Linear", "x + y"},
        {"Decay", "-y"},
        {"Harmonic", "-x"},
        {"Growth", "y"},
        {"Nonlinear", "x^2 + y^2"},
        {"Trigonometric", "sin(x) * y"}
    };

    void clearInput() {
        std::cin.clear();
        std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    }

    double getNumberInput(const std::string& prompt, double minVal, double maxVal) {
        double value;
        while (true) {
            std::cout << prompt;
            if (std::cin >> value && value >= minVal && value <= maxVal) {
                clearInput();
                return value;
            }
            std::cout << "Invalid input. Please enter a number between " 
                      << minVal << " and " << maxVal << ".\n";
            clearInput();
        }
    }

    void displayResults(const std::vector<std::pair<double, double>>& solution, 
                       bool detailed = true) {
        std::cout << std::fixed << std::setprecision(6);
        std::cout << "\nSolution:\n";
        std::cout << "x\t\ty\n";
        
        size_t step = detailed ? 1 : std::max(size_t(1), solution.size() / 10);
        
        for (size_t i = 0; i < solution.size(); i += step) {
            std::cout << solution[i].first << "\t\t" << solution[i].second << "\n";
        }
    }

public:
    void run() {
        while (true) {
            try {
                std::cout << "\n=== Differential Equation Solver ===\n";
                std::cout << "1. Use predefined equation\n";
                std::cout << "2. Enter custom equation\n";
                std::cout << "3. Show available functions\n";
                std::cout << "4. Exit\n";
                std::cout << "Choose option (1-4): ";

                int choice;
                std::cin >> choice;
                clearInput();

                if (choice == 4) break;
                if (choice == 3) {
                    std::cout << EquationParser::getAvailableFunctions() << std::endl;
                    continue;
                }

                std::function<double(double, double)> equation;

                if (choice == 1) {
                    std::cout << "\nAvailable equations (dy/dx =):\n";
                    int i = 1;
                    for (const auto& eq : predefinedEquations) {
                        std::cout << i++ << ". " << eq.first << ": " << eq.second << "\n";
                    }

                    int eqChoice = getNumberInput("Select equation (1-" + 
                        std::to_string(predefinedEquations.size()) + "): ", 
                        1, predefinedEquations.size());

                    auto it = predefinedEquations.begin();
                    std::advance(it, eqChoice - 1);
                    equation = parser.parseEquation(it->second);
                }
                else if (choice == 2) {
                    std::cout << "\nEnter equation (e.g., 'sin(x) + y^2'): ";
                    std::string eqStr;
                    std::getline(std::cin, eqStr);
                    equation = parser.parseEquation(eqStr);
                }
                else {
                    std::cout << "Invalid choice.\n";
                    continue;
                }

                // Get solving parameters
                double x0 = getNumberInput("Enter initial x value: ", -1000.0, 1000.0);
                double y0 = getNumberInput("Enter initial y value: ", -1000.0, 1000.0);
                double xEnd = getNumberInput("Enter final x value: ", x0, 1000.0);
                double stepSize = getNumberInput("Enter step size (0.001-1.0): ", 0.001, 1.0);

                // Solve the equation
                auto solution = DifferentialSolver::solve(equation, x0, y0, xEnd, stepSize);

                // Display options
                std::cout << "\nDisplay options:\n";
                std::cout << "1. Show all points\n";
                std::cout << "2. Show summary (10 points)\n";
                int displayChoice = getNumberInput("Choose display option (1-2): ", 1, 2);

                displayResults(solution, displayChoice == 1);

                // Calculate and display error estimates
                auto refinedSolution = DifferentialSolver::solve(equation, x0, y0, xEnd, stepSize/2);
                
                double maxError = 0.0;
                for (size_t i = 0; i < solution.size(); i++) {
                    double error = std::abs(solution[i].second - 
                        refinedSolution[i*2].second);
                    maxError = std::max(maxError, error);
                }

                std::cout << "\nEstimated maximum error: " << maxError << "\n";
            }
            catch (const mu::Parser::exception_type& e) {
                std::cout << "Error in equation: " << e.GetMsg() << "\n";
            }
            catch (const std::exception& e) {
                std::cout << "Error: " << e.what() << "\n";
            }
        }
    }
};

int main() {
    UserInterface ui;
    ui.run();
    return 0;
}



// #include <iostream>
// #include <vector>
// #include <functional>
// #include <cmath>

// class DifferentialSolver {
// private:
//     // Function type for differential equations: dy/dx = f(x,y)
//     using DiffFunction = std::function<double(double, double)>;
    
//     // RK4 method implementation
//     static double rk4Step(DiffFunction f, double x, double y, double h) {
//         double k1 = f(x, y);
//         double k2 = f(x + h/2, y + h*k1/2);
//         double k3 = f(x + h/2, y + h*k2/2);
//         double k4 = f(x + h, y + h*k3);
        
//         return y + (h/6) * (k1 + 2*k2 + 2*k3 + k4);
//     }

// public:
//     // Solve the differential equation using RK4 method
//     static std::vector<std::pair<double, double>> solve(
//         DiffFunction f,     // The differential equation dy/dx = f(x,y)
//         double x0,          // Initial x value
//         double y0,          // Initial y value
//         double xEnd,        // Final x value
//         double stepSize     // Step size (h)
//     ) {
//         std::vector<std::pair<double, double>> solution;
        
//         double x = x0;
//         double y = y0;
        
//         // Store initial point
//         solution.push_back({x, y});
        
//         // Solve until we reach xEnd
//         while (x < xEnd) {
//             // Adjust step size if we would overshoot xEnd
//             double h = std::min(stepSize, xEnd - x);
            
//             // Calculate next y using RK4
//             y = rk4Step(f, x, y, h);
//             x += h;
            
//             // Store point
//             solution.push_back({x, y});
//         }
        
//         return solution;
//     }
// };

// // Example usage
// int main() {
//     // Example 1: dy/dx = x + y
//     auto diff_eq1 = [](double x, double y) { return x + y; };
    
//     // Solve from x=0 to x=1 with initial condition y(0)=1
//     auto solution1 = DifferentialSolver::solve(diff_eq1, 0.0, 1.0, 1.0, 0.1);
    
//     std::cout << "Solution for dy/dx = x + y, y(0) = 1:\n";
//     std::cout << "x\t\ty\n";
//     for (const auto& point : solution1) {
//         std::cout << point.first << "\t\t" << point.second << "\n";
//     }
    
//     // Example 2: dy/dx = -y (exponential decay)
//     auto diff_eq2 = [](double x, double y) { return -y; };
    
//     // Solve from x=0 to x=2 with initial condition y(0)=1
//     auto solution2 = DifferentialSolver::solve(diff_eq2, 0.0, 1.0, 2.0, 0.1);
    
//     std::cout << "\nSolution for dy/dx = -y, y(0) = 1:\n";
//     std::cout << "x\t\ty\n";
//     for (const auto& point : solution2) {
//         std::cout << point.first << "\t\t" << point.second << "\n";
//     }
    
//     return 0;
// }