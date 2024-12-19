#include <iostream>
#include <cmath>
#include <iomanip>
#include <string>
#include <stdexcept>
#include <limits>
#include <vector>

// Expression tokenizer and evaluator
class ExpressionEvaluator {
private:
    std::string expr;
    size_t pos;

    // Skip whitespace
    void skipWhitespace() {
        while (pos < expr.length() && std::isspace(expr[pos])) pos++;
    }

    // Get next token
    bool isOperator(char c) const {
        return c == '+' || c == '-' || c == '*' || c == '/' || c == '^';
    }

    // Parse number
    double parseNumber() {
        skipWhitespace();
        size_t start = pos;
        
        // Handle negative numbers
        if (expr[pos] == '-') {
            pos++;
        }
        
        // Parse digits before decimal
        while (pos < expr.length() && (std::isdigit(expr[pos]) || expr[pos] == '.')) {
            pos++;
        }
        
        std::string numStr = expr.substr(start, pos - start);
        return std::stod(numStr);
    }

    // Parse function
    std::string parseFunction() {
        skipWhitespace();
        size_t start = pos;
        while (pos < expr.length() && std::isalpha(expr[pos])) {
            pos++;
        }
        return expr.substr(start, pos - start);
    }

    // Get operator precedence
    int getPrecedence(char op) const {
        switch (op) {
            case '+': case '-': return 1;
            case '*': case '/': return 2;
            case '^': return 3;
            default: return 0;
        }
    }

    // Apply operator
    double applyOperator(double a, double b, char op) const {
        switch (op) {
            case '+': return a + b;
            case '-': return a - b;
            case '*': return a * b;
            case '/': 
                if (b == 0) throw std::invalid_argument("Division by zero");
                return a / b;
            case '^': return std::pow(a, b);
            default: throw std::invalid_argument("Invalid operator");
        }
    }

public:
    explicit ExpressionEvaluator(const std::string& expression) : expr(expression), pos(0) {}

    double evaluate() {
        std::vector<double> values;
        std::vector<char> operators;

        while (pos < expr.length()) {
            skipWhitespace();
            
            if (pos >= expr.length()) break;

            char currentChar = expr[pos];

            if (std::isdigit(currentChar) || currentChar == '-') {
                values.push_back(parseNumber());
            }
            else if (std::isalpha(currentChar)) {
                std::string func = parseFunction();
                skipWhitespace();
                
                if (pos >= expr.length() || expr[pos] != '(') {
                    throw std::invalid_argument("Expected '(' after function");
                }
                pos++; // skip '('
                
                // Recursively evaluate function argument
                size_t start = pos;
                int parentheses = 1;
                while (pos < expr.length() && parentheses > 0) {
                    if (expr[pos] == '(') parentheses++;
                    else if (expr[pos] == ')') parentheses--;
                    pos++;
                }
                
                if (parentheses > 0) {
                    throw std::invalid_argument("Mismatched parentheses");
                }
                
                std::string argExpr = expr.substr(start, pos - start - 1);
                ExpressionEvaluator argEval(argExpr);
                double arg = argEval.evaluate();
                
                // Apply function
                if (func == "sin") values.push_back(std::sin(arg));
                else if (func == "cos") values.push_back(std::cos(arg));
                else if (func == "tan") values.push_back(std::tan(arg));
                else if (func == "log") values.push_back(std::log10(arg));
                else if (func == "ln") values.push_back(std::log(arg));
                else if (func == "sqrt") values.push_back(std::sqrt(arg));
                else throw std::invalid_argument("Unknown function: " + func);
            }
            else if (currentChar == '(') {
                operators.push_back(currentChar);
                pos++;
            }
            else if (currentChar == ')') {
                while (!operators.empty() && operators.back() != '(') {
                    double b = values.back(); values.pop_back();
                    double a = values.back(); values.pop_back();
                    char op = operators.back(); operators.pop_back();
                    values.push_back(applyOperator(a, b, op));
                }
                if (operators.empty()) {
                    throw std::invalid_argument("Mismatched parentheses");
                }
                operators.pop_back(); // Remove '('
                pos++;
            }
            else if (isOperator(currentChar)) {
                while (!operators.empty() && operators.back() != '(' && 
                       getPrecedence(operators.back()) >= getPrecedence(currentChar)) {
                    double b = values.back(); values.pop_back();
                    double a = values.back(); values.pop_back();
                    char op = operators.back(); operators.pop_back();
                    values.push_back(applyOperator(a, b, op));
                }
                operators.push_back(currentChar);
                pos++;
            }
            else {
                throw std::invalid_argument("Invalid character in expression");
            }
        }

        while (!operators.empty()) {
            if (operators.back() == '(' || operators.back() == ')') {
                throw std::invalid_argument("Mismatched parentheses");
            }
            double b = values.back(); values.pop_back();
            double a = values.back(); values.pop_back();
            char op = operators.back(); operators.pop_back();
            values.push_back(applyOperator(a, b, op));
        }

        if (values.size() != 1) {
            throw std::invalid_argument("Invalid expression");
        }

        return values[0];
    }
};

class ScientificCalculator {
private:
    // Constants
    static constexpr double PI = 3.14159265358979323846;
    static constexpr double E = 2.71828182845904523536;
    
    // Current angle mode (degrees or radians)
    bool isDegreeMode = true;

    // Convert degrees to radians
    double toRadians(double angle) const {
        return isDegreeMode ? angle * PI / 180.0 : angle;
    }

    // Convert radians to degrees
    double toDegrees(double angle) const {
        return isDegreeMode ? angle * 180.0 / PI : angle;
    }

    // Calculate factorial
    double factorial(double n) const {
        if (n < 0) {
            throw std::invalid_argument("Factorial is not defined for negative numbers");
        }
        if (n > 170) {
            throw std::overflow_error("Factorial result too large");
        }
        if (floor(n) != n) {
            throw std::invalid_argument("Factorial is only defined for integers");
        }
        
        double result = 1.0;
        for (int i = 2; i <= n; ++i) {
            result *= i;
        }
        return result;
    }

public:
    // Constructor
    ScientificCalculator() = default;

    // Toggle between degrees and radians
    void toggleAngleMode() {
        isDegreeMode = !isDegreeMode;
        std::cout << "Angle mode changed to: " << (isDegreeMode ? "Degrees" : "Radians") << std::endl;
    }

    // Basic arithmetic operations
    double add(double a, double b) const { return a + b; }
    double subtract(double a, double b) const { return a - b; }
    double multiply(double a, double b) const { return a * b; }
    
    double divide(double a, double b) const {
        if (std::abs(b) < std::numeric_limits<double>::epsilon()) {
            throw std::invalid_argument("Division by zero");
        }
        return a / b;
    }

    // Power and root functions
    double power(double base, double exponent) const {
        return std::pow(base, exponent);
    }

    double squareRoot(double number) const {
        if (number < 0) {
            throw std::invalid_argument("Square root of negative number");
        }
        return std::sqrt(number);
    }

    // Exponential and logarithmic functions
    double exp(double x) const { return std::exp(x); }
    
    double ln(double x) const {
        if (x <= 0) {
            throw std::invalid_argument("Logarithm of non-positive number");
        }
        return std::log(x);
    }

    double lne() const {
        return std::log(E);
    }
    
    double log10(double x) const {
        if (x <= 0) {
            throw std::invalid_argument("Logarithm of non-positive number");
        }
        return std::log10(x);
    }

    double logBase(double base, double x) const {
        if (x <= 0 || base <= 0 || base == 1) {
            throw std::invalid_argument("Invalid arguments for logarithm");
        }
        return std::log(x) / std::log(base);
    }

    // Trigonometric functions
    double sin(double angle) const { return std::sin(toRadians(angle)); }
    double cos(double angle) const { return std::cos(toRadians(angle)); }
    double tan(double angle) const {
        double rad = toRadians(angle);
        if (std::abs(std::cos(rad)) < std::numeric_limits<double>::epsilon()) {
            throw std::invalid_argument("Tangent undefined at this angle");
        }
        return std::tan(rad);
    }

    // Inverse trigonometric functions
    double arcsin(double x) const {
        if (x < -1 || x > 1) {
            throw std::invalid_argument("Arcsin argument must be between -1 and 1");
        }
        return toDegrees(std::asin(x));
    }

    double arccos(double x) const {
        if (x < -1 || x > 1) {
            throw std::invalid_argument("Arccos argument must be between -1 and 1");
        }
        return toDegrees(std::acos(x));
    }

    double arctan(double x) const {
        return toDegrees(std::atan(x));
    }

    // Hyperbolic functions
    double sinh(double x) const { return std::sinh(x); }
    double cosh(double x) const { return std::cosh(x); }
    double tanh(double x) const { return std::tanh(x); }

    // Factorial and combination functions
    double fact(double n) const { return factorial(n); }
    
    double nCr(double n, double r) const {
        if (r > n) {
            throw std::invalid_argument("r cannot be greater than n in nCr");
        }
        return factorial(n) / (factorial(r) * factorial(n - r));
    }

    double nPr(double n, double r) const {
        if (r > n) {
            throw std::invalid_argument("r cannot be greater than n in nPr");
        }
        return factorial(n) / factorial(n - r);
    }
};

// Helper function to print a dividing line
void printLine() {
    std::cout << "\n----------------------------------------\n";
}

// Interactive calculator function
void runInteractiveCalculator() {
    ScientificCalculator calc;
    std::string input;
    
    std::cout << "\nScientific Calculator" << std::endl;
    std::cout << "Enter expressions to evaluate (type 'exit' to quit)" << std::endl;
    std::cout << "Supported operations: +, -, *, /, ^" << std::endl;
    std::cout << "Supported functions: sin, cos, tan, log, ln, sqrt" << std::endl;
    std::cout << "Example: sin(45) + 2 * sqrt(16)" << std::endl;
    
    while (true) {
        std::cout << "\n> ";
        std::getline(std::cin, input);
        
        if (input == "exit" || input == "quit") {
            break;
        }
        
        try {
            ExpressionEvaluator evaluator(input);
            double result = evaluator.evaluate();
            std::cout << "= " << std::fixed << std::setprecision(6) << result << std::endl;
        }
        catch (const std::exception& e) {
            std::cerr << "Error: " << e.what() << std::endl;
        }
    }
}

int main() {
    ScientificCalculator calc;
    std::cout << std::fixed << std::setprecision(6);

    try {
        // Demonstrate basic arithmetic
        std::cout << "Basic Arithmetic Examples:" << std::endl;
        std::cout << "12 + 3 = " << calc.add(12, 3) << std::endl;
        std::cout << "12 - 3 = " << calc.subtract(12, 3) << std::endl;
        std::cout << "12 * 3 = " << calc.multiply(12, 3) << std::endl;
        std::cout << "12 / 3 = " << calc.divide(12, 3) << std::endl;

        printLine();

        // Demonstrate trigonometric functions (in degrees by default)
        std::cout << "Trigonometric Functions (in degrees):" << std::endl;
        std::cout << "sin(30) = " << calc.sin(30) << std::endl;
        std::cout << "cos(60) = " << calc.cos(60) << std::endl;
        std::cout << "tan(45) = " << calc.tan(45) << std::endl;

        printLine();

        // Switch to radians and demonstrate
        calc.toggleAngleMode();
        std::cout << "Trigonometric Functions (in radians):" << std::endl;
        std::cout << "sin(π/6) = " << calc.sin(M_PI/6) << std::endl;
        std::cout << "cos(π/3) = " << calc.cos(M_PI/3) << std::endl;
        std::cout << "tan(π/4) = " << calc.tan(M_PI/4) << std::endl;

        printLine();

        // Demonstrate logarithmic functions
        std::cout << "Logarithmic Functions:" << std::endl;
        std::cout << "ln(e) = " << calc.lne() << std::endl;
        std::cout << "log10(100) = " << calc.log10(100) << std::endl;
        std::cout << "log2(8) = " << calc.logBase(2, 8) << std::endl;

        printLine();

        // Demonstrate exponential and power functions
        std::cout << "Exponential and Power Functions:" << std::endl;
        std::cout << "e^2 = " << calc.exp(2) << std::endl;
        std::cout << "2^3 = " << calc.power(2, 3) << std::endl;
        std::cout << "√16 = " << calc.squareRoot(16) << std::endl;

        printLine();

        // Demonstrate factorial and combinations
        std::cout << "Factorial and Combinations:" << std::endl;
        std::cout << "5! = " << calc.fact(5) << std::endl;
        std::cout << "5C2 = " << calc.nCr(5, 2) << std::endl;
        std::cout << "5P2 = " << calc.nPr(5, 2) << std::endl;

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }

    // Run interactive calculator
    std::cout << "\nStarting interactive calculator mode..." << std::endl;
    runInteractiveCalculator();
    
    return 0;
}