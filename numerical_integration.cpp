#include <iostream>
#include <functional>
#include <vector>
#include <cmath>
#include <stdexcept>
#include <iomanip>
#include <map>
#include <string>
#include <stack>
#include <cctype>
#include <sstream>
#include <muParser.h> //TODO: install muPerser (muparser-2.3.5.zip)
#include <chrono>

/*
use: 
     g++ numerical_integration.cpp -o numerical_integration -lmuparser
*/
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
            result += 2.0 * f(a + i * h);
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

};


// Function to parse mathematical expressions using muParser
std::function<double(double)> parseFunction(const std::string& expression) {
    return [expression](double x) {
        mu::Parser parser;
        parser.SetExpr(expression);
        parser.DefineVar("x", &x); // Define 'x' as the variable in the expression
        
        try {
            return parser.Eval();
        } catch (mu::ParserError& e) {
            throw std::invalid_argument(e.GetMsg());
        }
    };
}

class UI {
private:
    static void clearScreen() {
        #ifdef _WIN32
            system("cls");
        #else
            system("clear");
        #endif
    }

    static void showProgressBar(int progress, int total) {
        const int barWidth = 50;
        float ratio = static_cast<float>(progress) / total;
        int pos = static_cast<int>(barWidth * ratio);

        std::cout << "[";
        for (int i = 0; i < barWidth; ++i) {
            if (i < pos) std::cout << "=";
            else if (i == pos) std::cout << ">";
            else std::cout << " ";
        }
        std::cout << "] " << int(ratio * 100.0) << "%\r";
        std::cout.flush();
    }

public:
    static void displayHeader() {
        clearScreen();
        std::cout << "╔════════════════════════════════════════════╗\n";
        std::cout << "║       Numerical Integration                ║\n";
        std::cout << "╚════════════════════════════════════════════╝\n\n";
    }

    static void displayMethodMenu() {
        std::cout << "\nAvailable Integration Methods:\n";
        std::cout << "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n";
        std::cout << "1. Trapezoidal Rule\n";
        std::cout << "2. Rectangular Rule\n";
        std::cout << "3. Simpson's Rule 1/8\n";
        std::cout << "4. Simpson's Rule 3/8\n";
        std::cout << "5. Boole's Rule\n";
        std::cout << "6. Romberg Integration\n";
        std::cout << "0. Exit Program\n";
        std::cout << "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n";
    }

    static void displayResult(std::string func,double result, const std::string& method, 
                            double a, double b, int n, double executionTime) {
        std::cout << "\n╔════════════════════════════════════════════╗\n";
        std::cout << "║              Integration Result            ║\n";
        std::cout << "╚════════════════════════════════════════════╝\n\n";
        std::cout << "Method: " << method << "\n";
        std::cout << "Function f(x): " << func << "\n";
        std::cout << "Interval: [" << a << ", " << b << "]\n";
        std::cout << "Subintervals: " << n << "\n";
        std::cout << "Result: " << std::fixed << std::setprecision(8) << result << "\n";
        std::cout << "Execution time: " << std::fixed << std::setprecision(3) 
                 << executionTime << " ms\n\n";
    }

    static void showCalculationProgress(int current, int total) {
        showProgressBar(current, total);
    }
};

class InputValidator {
public:
    static double getValidDouble(const std::string& prompt) {
        double value;
        while (true) {
            std::cout << prompt;
            if (std::cin >> value) {
                std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                return value;
            }
            std::cin.clear();
            std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            std::cout << "Invalid input. Please enter a valid number.\n";
        }
    }

    static int getValidInt(const std::string& prompt, int min = 1, int max = std::numeric_limits<int>::max()) {
        int value;
        while (true) {
            std::cout << prompt;
            if (std::cin >> value && value >= min && value <= max) {
                std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                return value;
            }
            std::cin.clear();
            std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            std::cout << "Invalid input. Please enter a number between " << min << " and " << max << ".\n";
        }
    }

    static std::string getValidFunction() {
        std::string funcExpr;
        std::cout << "Enter a mathematical function using x as variable\n";
        std::cout << "Examples: sin(x), x^2 + 2*x, exp(x)\n";
        std::cout << "Function: ";
        std::getline(std::cin, funcExpr);
        std::cout<<std::endl;
        return funcExpr;
    }
};

int main() {
    while (true) {
        try {
            UI::displayHeader();
            
            // Get function expression
            std::string funcExpr = InputValidator::getValidFunction();
            
            // Get integration bounds
            double a = InputValidator::getValidDouble("Enter lower bound (a): ");
            double b = InputValidator::getValidDouble("Enter upper bound (b): ");
            
            if (a >= b) {
                throw std::runtime_error("Lower bound must be less than upper bound!");
            }
            
            // Get number of subintervals
            int n = InputValidator::getValidInt("Enter number of subintervals (1-1000000): ", 1, 1000000);
            
            // Display method menu and get choice
            UI::displayMethodMenu();
            int methodChoice = InputValidator::getValidInt("Select method (0-6): ", 0, 6);
            
            if (methodChoice == 0) {
                std::cout << "Exit!\n";
                break;
            }
            
            // Parse the function
            auto f = parseFunction(funcExpr);
            
            // Create integrator and get method
            NumericalIntegrator integrator;
            std::string methodName;
            std::function<double(std::function<double(double)>, double, double, int)> method;
            
            switch (methodChoice) {
                case 1:
                    method = integrator.trapezoidal;
                    methodName = "Trapezoidal Rule";
                    break;
                case 2:
                    method = integrator.rectangular;
                    methodName = "Rectangular Rule";
                    break;
                case 3:
                    method = integrator.simpsons;
                    methodName = "Simpson's Rule 1/8";
                    break;
                case 4:
                    method = integrator.simpsons38;
                    methodName = "Simpson's Rule 3/8";
                    break;
                case 5:
                    method = integrator.booles;
                    methodName = "Boole's Rule";
                    break;
                case 6:
                    method = integrator.romberg;
                    methodName = "Romberg Integration";
                    break;
                default:
                    throw std::runtime_error("Invalid method choice");
            }
            
            // Calculate with progress indication
            std::cout << "\nCalculating...\n";
            auto start = std::chrono::high_resolution_clock::now();
            
            double result = method(f, a, b, n);
            
            auto end = std::chrono::high_resolution_clock::now();
            double executionTime = std::chrono::duration<double, std::milli>(end - start).count();
            
            UI::displayResult(funcExpr,result, methodName, a, b, n, executionTime);
            
            // Ask to continue
            char continueChoice;
            std::cout << "Would you like to continue? (y/n): ";
            std::cin >> continueChoice;
            std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            
            if (std::tolower(continueChoice) != 'y') {
                std::cout << "Quit!\n";
                break;
            }

        } catch (const std::exception& e) {
            std::cerr << "\nError: " << e.what() << "\n\n";
            std::cout << "Press Enter to continue...";
            std::cin.get();
        }
    }

    return 0;
}
