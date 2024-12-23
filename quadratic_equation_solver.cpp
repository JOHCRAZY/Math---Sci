#include <iostream>
#include <cmath>
#include <complex>
#include <iomanip>
#include <limits>
#include <vector>
#include <chrono>
#include <cstring>
#include <string>
class QuadraticSolver {
public:
    enum class RootType {
        REAL_DISTINCT,
        REAL_EQUAL,
        COMPLEX
    };

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
};


class InteractiveQuadraticSolver : public QuadraticSolver{
    private:
        struct EquationHistory {
            double a, b, c;
            SolutionResult solution;
            std::string timestamp;
        };
        
        std::vector<EquationHistory> history;
        
        static std::string getCurrentTimestamp() {
            auto now = std::chrono::system_clock::now();
            auto time = std::chrono::system_clock::to_time_t(now);
            std::string timestamp = std::ctime(&time);
            timestamp.pop_back();
            return timestamp;
        }
    
       static  void printBox(const std::string& message, const char* decoratorChar = "═") {
            std::vector<std::string> lines;
            std::istringstream stream(message);
            std::string line;
            
            while (std::getline(stream, line)) {
                lines.push_back("║ " + line);
            }
            
            size_t maxLength = 0;
            for (const auto& line : lines) {
                maxLength = std::max(maxLength, line.length());
            }
            
            maxLength += 2; 
            
            std::cout << "╔" << std::string(maxLength - 2, '=') << "╗\n";
                        for (const auto& line : lines) {
                std::cout << line << std::string(maxLength - line.length(),' ') << "  ║\n";
            }
            std::cout << "╚" << std::string(maxLength - 2, '=') << "╝\n";
        }
    
        static void printFormattedMenu() {
            std::vector<std::string> menuItems = {
                "solve     - Solve a new quadratic equation",
                "history   - View solution history",
                "clear     - Clear solution history",
                "help      - Show this menu",
                "exit      - Exit the program"
            };
    
            const size_t boxWidth = 50;
            std::string title = "Quadratic Equation Solver";
            size_t titlePadding = (boxWidth - title.length()) / 2;
    
            std::cout << "╔" << std::string(boxWidth, '=') << "╗\n";
            std::cout << "║" << std::string(titlePadding, ' ') << title 
                      << std::string(boxWidth - titlePadding - title.length(), ' ') << "║\n";
            std::cout << "╠" << std::string(boxWidth, '=') << "╣\n";
    
            for (const auto& item : menuItems) {
                std::cout << "║  " << item 
                          << std::string(boxWidth - item.length() - 3, ' ') << " ║\n";
            }
    
            std::cout << "╚" << std::string(boxWidth, '=') << "╝\n";
        }
    
        static void printEquation(double a, double b, double c) {
            std::ostringstream eqn;
            eqn << std::fixed << std::setprecision(2);
            
            eqn << a << "x² ";
            if (b >= 0) eqn << "+ " << b << "x ";
            else eqn << "- " << -b << "x ";
            
            if (c >= 0) eqn << "+ " << c;
            else eqn << "- " << -c;
            
            eqn << " = 0";
            
            printBox(eqn.str(), "-");
        }
    
        static void printSolution(const SolutionResult& solution) {
            std::cout << std::fixed << std::setprecision(4);
            std::vector<std::string> lines;
                switch (solution.type) {
                case RootType::REAL_DISTINCT: {
                    lines.push_back("Two distinct real roots:");
                    lines.push_back("────────────────────────");
                    std::ostringstream x1, x2;
                    x1 << "x₁ = " << solution.root1.real();
                    x2 << "  x₂ = " << solution.root2.real();
                    lines.push_back(x1.str());
                    lines.push_back(x2.str());
                    break;
                }
                case RootType::REAL_EQUAL: {
                    lines.push_back("One repeated real root:");
                    lines.push_back("────────────────────────");
                    std::ostringstream x;
                    x << "  x = " << solution.root1.real();
                    lines.push_back(x.str());
                    break;
                }
                case RootType::COMPLEX: {
                    lines.push_back("Two complex roots:");
                    lines.push_back("───────────────────");
                    std::ostringstream x1, x2;
                    x1 << "  x₁ = " << solution.root1.real() << " + " 
                       << solution.root1.imag() << "i";
                    x2 << "  x₂ = " << solution.root2.real() << " - " 
                       << std::abs(solution.root2.imag()) << "i";
                    lines.push_back(x1.str());
                    lines.push_back(x2.str());
                    break;
                }
            }
            
            size_t maxLength = 0;
            for (const auto& line : lines) {
                maxLength = std::max(maxLength, line.length());
            }
            
            size_t boxWidth = maxLength;  
            
            std::cout << "╔" << std::string(boxWidth - 20, '=') << "╗\n";
            for (const auto& line : lines) {
                std::cout << "  " << line 
                         << std::string(boxWidth - line.length(), ' ') << "\n";
            }
            std::cout << "╚" << std::string(boxWidth - 20, '=') << "╝\n";
        }
    
        void solveEquation() {
            double a, b, c;
            std::cout << "\nEnter coefficients (a b c): ";
            std::string input;
            std::getline(std::cin, input);
            std::istringstream iss(input);
        
            if (!(iss >> a >> b >> c)) {
                printBox(" Invalid input! Please enter three numbers.", "!");
                return;
            }
    
            if (std::abs(a) < 1e-10) {
                printBox(" 'a' coefficient cannot be zero!", "!");
                return;
            }
    
            std::cout << "\n Solving.....\n";
            printEquation(a, b, c);
            
            auto solution = QuadraticSolver::solve(a, b, c);
            printSolution(solution);
            
            // Save to history
            history.push_back({a, b, c, solution, getCurrentTimestamp()});
        }
    
        void showHistory() const {
            if (history.empty()) {
                printBox(" No solutions in history", "=");
                return;
            }
    
            for (size_t i = 0; i < history.size(); ++i) {
                const auto& entry = history[i];
                std::cout << "\n[ " << entry.timestamp << " ]\n";
                printEquation(entry.a, entry.b, entry.c);
                printSolution(entry.solution);
                std::cout << "\n";
            }
        }
    
    public:
        void run() {
            printFormattedMenu();
            
            std::string command;
            while (true) {
                std::cout << "\n> ";
                std::getline(std::cin, command);
                
                if (command == "exit") {
                    printBox("Goodbye!", "=");
                    break;
                }
                else if (command == "solve") {
                    solveEquation();
                }
                else if (command == "history") {
                    showHistory();
                }
                else if (command == "clear") {
                    history.clear();
                    printBox("History cleared", "-");
                }
                else if (command == "help") {
                    printFormattedMenu();
                }
                else if (!command.empty()) {
                    printBox("Unknown command. Type 'help' for menu.", "!");
                }
            }
        }
    };
    
    int main() {
        try {
            InteractiveQuadraticSolver solver;
            solver.run();
        }
        catch (const std::exception& e) {
            std::cerr << "\n╔════════════════════════════════════╗\n";
            std::cerr << "║ Error: " << e.what() << std::string(30 - strlen(e.what()), ' ') << " ║\n";
            std::cerr << "╚════════════════════════════════════╝\n";
            return 1;
        }
    
        return 0;
    }
