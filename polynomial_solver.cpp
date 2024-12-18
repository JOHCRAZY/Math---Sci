#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

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

// Main function to demonstrate polynomial operations
int main() {
  // Create some polynomials
  Polynomial p1({1, 2, 3}); // 3x^2 + 2x + 1
  Polynomial p2({4, 5});    // 5x + 4

  std::cout << "Polynomial 1: ";
  p1.print();
  std::cout << std::endl;

  std::cout << "Polynomial 2: ";
  p2.print();
  std::cout << std::endl;

  // Addition
  std::cout << "Addition (p1 + p2): ";
  (p1 + p2).print();
  std::cout << std::endl;

  // Subtraction
  std::cout << "Subtraction (p1 - p2): ";
  (p1 - p2).print();
  std::cout << std::endl;

  // Multiplication
  std::cout << "Multiplication (p1 * p2): ";
  (p1 * p2).print();
  std::cout << std::endl;

  return 0;
}