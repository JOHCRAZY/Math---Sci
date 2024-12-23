#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>
#include <map>
#include <sstream>
#include <regex>
#include <cstring>

class Polynomial {
private:
  // Vector to store coefficients, index represents the power
  std::vector<double> coeffs;

  // Trim leading zero coefficients
  void normalize() {
    while (!coeffs.empty() && std::abs(coeffs.back()) < 1e-10) {
      coeffs.pop_back();
    }
    if (coeffs.empty()) {
      coeffs.push_back(0);
    }
  }

public:
  // Constructors
  Polynomial() : coeffs({0}) {}

  // Constructor from a vector of coefficients
  Polynomial(const std::vector<double> &coefficients) : coeffs(coefficients) {
    normalize();
  }

  // Constructor from a single coefficient
  Polynomial(double coeff) : coeffs({coeff}) {}

  // Degree of the polynomial
  int degree() const { return coeffs.size() - 1; }

  // Get coefficient at a specific power
  double getCoeff(int power) const {
    return (power >= 0 && power < coeffs.size()) ? coeffs[power] : 0;
  }

  double evaluate(double x) const {
    double result = 0.0;
    for (size_t i = 0; i < coeffs.size(); ++i) {
        result += coeffs[i] * std::pow(x, i);
    }
    return result;
}

  // Addition operator
  Polynomial operator+(const Polynomial &other) const {
    std::vector<double> result_coeffs(
        std::max(coeffs.size(), other.coeffs.size()), 0.0);

    for (size_t i = 0; i < coeffs.size(); ++i) {
      result_coeffs[i] += coeffs[i];
    }

    for (size_t i = 0; i < other.coeffs.size(); ++i) {
      result_coeffs[i] += other.coeffs[i];
    }

    return Polynomial(result_coeffs);
  }

  // Subtraction operator
  Polynomial operator-(const Polynomial &other) const {
    std::vector<double> result_coeffs(
        std::max(coeffs.size(), other.coeffs.size()), 0.0);

    for (size_t i = 0; i < coeffs.size(); ++i) {
      result_coeffs[i] += coeffs[i];
    }

    for (size_t i = 0; i < other.coeffs.size(); ++i) {
      result_coeffs[i] -= other.coeffs[i];
    }

    return Polynomial(result_coeffs);
  }

  // Multiplication operator
  Polynomial operator*(const Polynomial &other) const {
    if (degree() < 0 || other.degree() < 0) {
      return Polynomial(0);
    }

    std::vector<double> result_coeffs(coeffs.size() + other.coeffs.size() - 1,
                                      0.0);

    for (size_t i = 0; i < coeffs.size(); ++i) {
      for (size_t j = 0; j < other.coeffs.size(); ++j) {
        result_coeffs[i + j] += coeffs[i] * other.coeffs[j];
      }
    }

    return Polynomial(result_coeffs);
  }

  // Print polynomial
  void print() const {
    if (coeffs.empty() || (coeffs.size() == 1 && std::abs(coeffs[0]) < 1e-10)) {
      std::cout << "0";
      return;
    }

    bool first_term = true;
    for (int i = coeffs.size() - 1; i >= 0; --i) {
      // Skip terms with effectively zero coefficients
      if (std::abs(coeffs[i]) < 1e-10)
        continue;

      // Sign and spacing
      if (!first_term) {
        std::cout << (coeffs[i] > 0 ? " + " : " - ");
      } else if (coeffs[i] < 0) {
        std::cout << "-";
      }
      first_term = false;

      // Coefficient
      double abs_coeff = std::abs(coeffs[i]);
      if (abs_coeff != 1 || i == 0) {
        std::cout << abs_coeff;
      }

      // Variable and power
      if (i > 0) {
        std::cout << "x";
        if (i > 1) {
          std::cout << "^" << i;
        }
      }
    }
  }
};

void menu() {
  // Define box drawing characters
  const char* TOP_LEFT = "╔";
  const char* TOP_RIGHT = "╗";
  const char* BOTTOM_LEFT = "╚";
  const char* BOTTOM_RIGHT = "╝";
  const char* HORIZONTAL = "═";
  const char* VERTICAL = "║";
  const char* TITLE_LEFT = "╠";
  const char* TITLE_RIGHT = "╣";

  // Store menu items in a vector for easy processing
  std::vector<std::string> menuItems = {
      "create <name> <coefficients> - Create a new polynomial (e.g., create p1 1 2 3)",
      "eval <name> <x>            - Evaluate polynomial at x",
      "add <p1> <p2>             - Add two polynomials",
      "subtract <p1> <p2>        - Subtract two polynomials",
      "multiply <p1> <p2>        - Multiply two polynomials",
      "list                      - List all stored polynomials",
      "history                   - Show operation history",
      "clear                     - Clear all stored polynomials",
      "help                      - Show this help message",
      "exit                      - Exit the program"
  };

  // Find the longest line to determine box width
  size_t maxLength = 0;
  for (const auto& item : menuItems) {
      maxLength = std::max(maxLength, item.length());
  }
  
  // Add padding to width
  size_t boxWidth = maxLength + 4;  // 2 spaces on each side
  std::string title = "Polynomial Solver";
  
  // Center the title
  size_t titlePadding = (boxWidth - title.length()) / 2;

  // Print top border
  std::cout << TOP_LEFT;
  for (size_t i = 0; i < boxWidth; i++) std::cout << HORIZONTAL;
  std::cout << TOP_RIGHT << "\n";

  // Print title
  std::cout << VERTICAL;
  std::cout << std::string(titlePadding, ' ') << title;
  std::cout << std::string(boxWidth - titlePadding - title.length(), ' ');
  std::cout << VERTICAL << "\n";

  // Print title separator
  std::cout << TITLE_LEFT;
  for (size_t i = 0; i < boxWidth; i++) std::cout << HORIZONTAL;
  std::cout << TITLE_RIGHT << "\n";

  // Print menu items
  for (const auto& item : menuItems) {
      std::cout << VERTICAL << "  " << item;
      std::cout << std::string(boxWidth - item.length() - 2, ' ') << VERTICAL << "\n";
  }

  // Print bottom border
  std::cout << BOTTOM_LEFT;
  for (size_t i = 0; i < boxWidth; i++) std::cout << HORIZONTAL;
  std::cout << BOTTOM_RIGHT << "\n";
}


int main() {
  std::vector<Polynomial> history;
  std::map<std::string, Polynomial> stored_polynomials;


  menu();

  std::string command;
  while (true) {
      std::cout << "\n -$ ";
      std::string line;
      std::getline(std::cin, line);
      std::istringstream iss(line);
      iss >> command;

      if (command == "exit") {
          break;
      }
      else if (command == "create") {
          std::string name;
          iss >> name;
          std::vector<double> coeffs;
          int coeff;
          while (iss >> coeff) {
              coeffs.push_back(coeff);
          }

          if (!coeffs.empty()) {
              stored_polynomials[name] = Polynomial(coeffs);
              std::cout << "Created polynomial " << name << ": ";
              stored_polynomials[name].print();
              std::cout << std::endl;
              history.push_back(stored_polynomials[name]);
          }
      }
      else if (command == "eval") {
          std::string name;
          double x;
          iss >> name >> x;
          if (stored_polynomials.count(name)) {
              double result = stored_polynomials[name].evaluate(x);
              std::cout << name << " at x = " << x << " equals: " << result << std::endl;
          }
          else {
              std::cout << "Polynomial " << name << " not found!" << std::endl;
          }
      }
      else if (command == "add") {
          std::string p1_name, p2_name;
          iss >> p1_name >> p2_name;
          if (stored_polynomials.count(p1_name) && stored_polynomials.count(p2_name)) {
              Polynomial result = stored_polynomials[p1_name] + stored_polynomials[p2_name];
              std::cout << "Result of " << p1_name << " + " << p2_name << " = ";
              result.print();
              stored_polynomials[(p1_name + "+" + p2_name)] = result;
              std::cout << std::endl;
              history.push_back(result);
          }
          else {
              std::cout << "One or both polynomials not found!" << std::endl;
          }
      }
      else if (command == "subtract") {
          std::string p1_name, p2_name;
          iss >> p1_name >> p2_name;
          if (stored_polynomials.count(p1_name) && stored_polynomials.count(p2_name)) {
              Polynomial result = stored_polynomials[p1_name] - stored_polynomials[p2_name];
              std::cout << "Result of " << p1_name << " - " << p2_name << " = ";
              result.print();
              stored_polynomials[(p1_name + "-" + p2_name)] = result;
              std::cout << std::endl;
              history.push_back(result);
          }
          else {
              std::cout << "One or both polynomials not found!" << std::endl;
          }
      }
      else if (command == "multiply") {
          std::string p1_name, p2_name;
          iss >> p1_name >> p2_name;
          if (stored_polynomials.count(p1_name) && stored_polynomials.count(p2_name)) {
              Polynomial result = stored_polynomials[p1_name] * stored_polynomials[p2_name];
              std::cout << "Result of " << p1_name << " * " << p2_name << " = ";
              result.print();
              stored_polynomials[(p1_name + "*" + p2_name)] = result;
              std::cout << std::endl;
              history.push_back(result);
          }
          else {
              std::cout << "One or both polynomials not found!" << std::endl;
          }
      }
      else if (command == "list") {
          if (stored_polynomials.empty()) {
              std::cout << "No polynomials stored." << std::endl;
          }
          else {
              std::cout << "Stored polynomials:" << std::endl;
              for (const auto& pair : stored_polynomials) {
                  std::cout << pair.first << ": ";
                  pair.second.print();
                  std::cout << std::endl;
              }
          }
      }
      else if (command == "history") {
          if (history.empty()) {
              std::cout << "No operations in history." << std::endl;
          }
          else {
              std::cout << "Operation history:" << std::endl;
              for (size_t i = 0; i < history.size(); ++i) {
                  std::cout << i + 1 << ": ";
                  history[i].print();
                  std::cout << std::endl;
              }
          }
      }
      else if (command == "clear") {
          stored_polynomials.clear();
          history.clear();
          std::cout << "All polynomials and history cleared." << std::endl;
      }
      else if (command == "help") {
        menu();
      }
      else {
          std::cout << "Unknown command. Type 'help' for available commands." << std::endl;
      }
  }

  return 0;
}