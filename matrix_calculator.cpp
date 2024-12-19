#include <iostream>
#include <vector>
#include <stdexcept>
#include <iomanip>
#include <cmath>

class Matrix {
private:
    std::vector<std::vector<double>> data;
    int rows, cols;

    // Helper function for Gaussian elimination with partial pivoting
    void gaussianElimination(std::vector<std::vector<double>>& augmented) {
        int n = augmented.size();
        
        for (int i = 0; i < n; ++i) {
            // Partial pivoting
            int max_row = i;
            for (int k = i + 1; k < n; ++k) {
                if (std::abs(augmented[k][i]) > std::abs(augmented[max_row][i])) {
                    max_row = k;
                }
            }
            
            // Swap maximum row with current row
            std::swap(augmented[max_row], augmented[i]);
            
            // Make all rows below this one 0 in current column
            for (int k = i + 1; k < n; ++k) {
                double factor = augmented[k][i] / augmented[i][i];
                
                for (int j = i; j < 2*n; ++j) {
                    augmented[k][j] -= factor * augmented[i][j];
                }
            }
        }
        
        // Back substitution
        for (int i = n - 1; i >= 0; --i) {
            for (int k = i - 1; k >= 0; --k) {
                double factor = augmented[k][i] / augmented[i][i];
                
                for (int j = 2*n - 1; j >= i; --j) {
                    augmented[k][j] -= factor * augmented[i][j];
                }
            }
        }
    }

public:
    // Constructor
    Matrix(int r, int c) : rows(r), cols(c), data(r, std::vector<double>(c, 0.0)) {}
    
    // Constructor from 2D vector
    Matrix(const std::vector<std::vector<double>>& input) : 
        data(input), rows(input.size()), cols(input[0].size()) {}

    // Get number of rows
    int getRows() const { return rows; }
    
    // Get number of columns
    int getCols() const { return cols; }

    // Access element
    double& operator()(int r, int c) {
        if (r < 0 || r >= rows || c < 0 || c >= cols) {
            throw std::out_of_range("Matrix index out of range");
        }
        return data[r][c];
    }

    // Const version of element access
    const double& operator()(int r, int c) const {
        if (r < 0 || r >= rows || c < 0 || c >= cols) {
            throw std::out_of_range("Matrix index out of range");
        }
        return data[r][c];
    }

    // Matrix addition
    Matrix operator+(const Matrix& other) const {
        if (rows != other.rows || cols != other.cols) {
            throw std::invalid_argument("Matrix dimensions must match for addition");
        }
        
        Matrix result(rows, cols);
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                result(i, j) = data[i][j] + other.data[i][j];
            }
        }
        return result;
    }

    // Matrix subtraction
    Matrix operator-(const Matrix& other) const {
        if (rows != other.rows || cols != other.cols) {
            throw std::invalid_argument("Matrix dimensions must match for subtraction");
        }
        
        Matrix result(rows, cols);
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                result(i, j) = data[i][j] - other.data[i][j];
            }
        }
        return result;
    }

    // Matrix multiplication
    Matrix operator*(const Matrix& other) const {
        if (cols != other.rows) {
            throw std::invalid_argument("Incompatible matrix dimensions for multiplication");
        }
        
        Matrix result(rows, other.cols);
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < other.cols; ++j) {
                for (int k = 0; k < cols; ++k) {
                    result(i, j) += data[i][k] * other.data[k][j];
                }
            }
        }
        return result;
    }

    // Matrix inversion (for square matrices)
    Matrix inverse()  {
        if (rows != cols) {
            throw std::invalid_argument("Only square matrices can be inverted");
        }

        // Create augmented matrix [A|I]
        std::vector<std::vector<double>> augmented(rows, std::vector<double>(2*rows, 0.0));
        
        // Fill left side with original matrix
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < rows; ++j) {
                augmented[i][j] = data[i][j];
            }
            
            // Fill right side with identity matrix
            augmented[i][i + rows] = 1.0;
        }

        // Perform Gaussian elimination
        gaussianElimination(augmented);

        // Normalize rows and extract inverse
        Matrix inverse_matrix(rows, rows);
        for (int i = 0; i < rows; ++i) {
            // Check if matrix is singular (determinant is zero)
            if (std::abs(augmented[i][i]) < 1e-10) {
                throw std::runtime_error("Matrix is not invertible");
            }

            for (int j = 0; j < rows; ++j) {
                inverse_matrix(i, j) = augmented[i][j + rows] / augmented[i][i];
            }
        }

        return inverse_matrix;
    }

    // Print matrix
    void print() const {
        for (int i = 0; i < rows; ++i) {
            std::cout<<"| ";
            for (int j = 0; j < cols; ++j) {
                if(j != 0)std::cout << std::fixed << std::setprecision(2)<< std::setw(10);
                std::cout  << data[i][j] << " ";
            }
            std::cout << "|" << std::endl;
        }
    }
};

int main() {
    try {
        // Example usage of matrix operations
        
        // Create matrices
        Matrix m1({{1, 2}, {3, 4},{4,5}});
        Matrix m2({{5, 6}, {7, 8}});

        std::cout << "Matrix 1:" << std::endl;
        m1.print();
        std::cout << "\nMatrix 2:" << std::endl;
        m2.print();

        // Addition
        std::cout << "\nAddition (m1 + m2):" << std::endl;
        (m1 + m2).print();

        // Subtraction
        std::cout << "\nSubtraction (m1 - m2):" << std::endl;
        (m1 - m2).print();

        // Multiplication
        std::cout << "\nMultiplication (m1 * m2):" << std::endl;
        (m1 * m2).print();

        // Inversion
        std::cout << "\nInverse of Matrix 1:" << std::endl;
        m1.inverse().print();

        // Different sized matrix for multiplication
        Matrix m3(2, 3);
        for (int i = 0; i < 2; ++i) {
            for (int j = 0; j < 3; ++j) {
                m3(i, j) = i + j;
            }
        }

        Matrix m4(3, 2);
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 2; ++j) {
                m4(i, j) = i * j;
            }
        }

        std::cout << "\nMatrix 3:" << std::endl;
        m3.print();
        std::cout << "\nMatrix 4:" << std::endl;
        m4.print();

        std::cout << "\nMultiplication (m3 * m4):" << std::endl;
        (m3 * m4).print();

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}