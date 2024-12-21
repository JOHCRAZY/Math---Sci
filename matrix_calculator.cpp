#include <algorithm>
#include <iostream>
#include <vector>
#include <stdexcept>
#include <iomanip>
#include <cmath>
#include <map>
#include <sstream>
#include <string>

std::string extractMatrixString(const std::string& command);
class Matrix {
private:
    std::vector<std::vector<double>> data;
    int rows, cols;


    static std::vector<std::vector<double>> parseMatrixString(const std::string& input) {
        std::vector<std::vector<double>> result;
        

        auto isNotSpace = [](unsigned char ch) { return !std::isspace(ch); };

    // Remove leading and trailing spaces
    auto trimmedStart = std::find_if(input.begin(), input.end(), isNotSpace);
    auto trimmedEnd = std::find_if(input.rbegin(), input.rend(), isNotSpace).base();

    if (trimmedStart >= trimmedEnd || 
        *trimmedStart != '[' || *(trimmedEnd - 1) != ']') {
        throw std::invalid_argument("Invalid matrix format: must be enclosed in []");
    }
        
        std::string cleaned = input.substr(1, input.length()-2);
        std::stringstream rows_ss(cleaned);
        std::string row_str;
        size_t expected_cols = 0;
        
        while (getline(rows_ss, row_str, ';')) {
            std::stringstream row_ss(row_str);
            std::vector<double> row_values;
            std::string value_str;
            
            while (getline(row_ss, value_str, ',')) {
                // Trim whitespace
                value_str.erase(0, value_str.find_first_not_of(" \t"));
                value_str.erase(value_str.find_last_not_of(" \t") + 1);
                
                if (value_str.empty()) {
                    throw std::invalid_argument("Empty value in matrix");
                }
                
                try {
                    double value = std::stod(value_str);
                    row_values.push_back(value);
                } catch (const std::exception&) {
                    throw std::invalid_argument("Invalid number format: " + value_str);
                }
            }
            
            if (row_values.empty()) {
                throw std::invalid_argument("Empty row in matrix");
            }
            
            if (expected_cols == 0) {
                expected_cols = row_values.size();
            } else if (row_values.size() != expected_cols) {
                throw std::invalid_argument("Inconsistent number of columns");
            }
            
            result.push_back(row_values);
        }
        
        if (result.empty()) {
            throw std::invalid_argument("Empty matrix");
        }
        
        return result;
    }

    //     // Helper function for Gaussian elimination with partial pivoting
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
    // Constructors
    Matrix() : rows(0), cols(0) {}  // Creates an empty 0x0 matrix
    Matrix(int r, int c) : rows(r), cols(c), data(r, std::vector<double>(c, 0.0)) {}
    Matrix(const std::vector<std::vector<double>>& input) : data(input), rows(input.size()), cols(input[0].size()) {}
    
    // Static method to create matrix from string
    static Matrix fromString(const std::string& input) {
        auto parsed = parseMatrixString(input);
        return Matrix(parsed);
    }

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

    // Add this to the Matrix class
double determinant() const {
    if (rows != cols) {
        throw std::invalid_argument("Determinant can only be calculated for square matrices");
    }

    // For 1x1 matrix
    if (rows == 1) {
        return data[0][0];
    }

    // For 2x2 matrix - direct calculation
    if (rows == 2) {
        return data[0][0] * data[1][1] - data[0][1] * data[1][0];
    }

    // For larger matrices, use LU decomposition
    // Create a copy of the matrix that we can modify
    std::vector<std::vector<double>> lu = data;
    double det = 1.0;
    
    // Perform LU decomposition and calculate determinant
    for (int i = 0; i < rows; i++) {
        // Check for zero pivot
        if (std::abs(lu[i][i]) < 1e-10) {
            return 0.0;  // Matrix is singular
        }

        det *= lu[i][i];  // Multiply diagonal elements
        
        // Update remaining elements
        for (int j = i + 1; j < rows; j++) {
            double factor = lu[j][i] / lu[i][i];
            for (int k = i + 1; k < rows; k++) {
                lu[j][k] -= factor * lu[i][k];
            }
        }
    }

    return det;
}

    // Enhanced print with formatting options
    void print(const std::string& name = "") const {
        if (!name.empty()) {
            std::cout << "\033[1;34m" << name << " = \033[0m\n";
        }
        for (int i = 0; i < rows; ++i) {
            std::cout << "\033[1;33m|\033[0m ";
            for (int j = 0; j < cols; ++j) {
                if(j != 0) std::cout << std::fixed << std::setprecision(2) << std::setw(10);
                std::cout << data[i][j] << " ";
            }
            std::cout << "\033[1;33m|\033[0m\n";
        }
    }
};

class MatrixCalculator {
private:
    std::map<std::string, Matrix> stored_matrices;
    
    void displayHelp() {
        std::cout << "\033[1;35m╔═══════════════════  Matrix Calculator Commands ══════════════════╗\n";
        std::cout << "║ 1. Create matrix:  create A [2,3,4; 5,6,7]                       ║\n";
        std::cout << "║ 2. Show matrix:    show A                                        ║\n";
        std::cout << "║ 3. Add:            add A B C  (C = A + B)                        ║\n";
        std::cout << "║ 4. Subtract:       sub A B C  (C = A - B)                        ║\n";
        std::cout << "║ 5. Multiply:       mul A B C  (C = A * B)                        ║\n";
        std::cout << "║ 6. Inverse:        inv A B    (B = A^-1)                         ║\n";
        std::cout << "║ 7. Determinant:     det A                                        ║\n";
        std::cout << "║ 8. Delete:          del A                                        ║\n";
        std::cout << "║ 9. List matrices:  ls                                            ║\n";
        std::cout << "║ 10. Clear screen:  clear                                         ║\n";
        std::cout << "║ 11. Help:          help                                          ║\n";
        std::cout << "║ 12. Exit:          exit                                          ║\n";
        std::cout << "╚══════════════════════════════════════════════════════════════════╝\n";
        std::cout << "\033[1;36mMatrix Format: [row1_values; row2_values; ...]\033[0m\n";
        std::cout << "Example: [1,2,3; 4,5,6; 7,8,9]\n\033[0m";
    }

    void displayPrompt() {
        std::cout << "\033[1;32mmatrix>\033[0m ";
    }

    void listMatrices() {
        if (stored_matrices.empty()) {
            std::cout << "\033[1;31mNo matrices stored.\033[0m\n";
            return;
        }
        std::cout << "\033[1;36mStored Matrices:\033[0m\n";
        for (const auto& pair : stored_matrices) {
            std::cout << pair.first << " (" << pair.second.getRows() << "x" 
                     << pair.second.getCols() << ")\n";
        }
    }

public:
    void run() {
        std::string command;
        displayHelp();

        while (true) {
            displayPrompt();
            std::getline(std::cin, command);

            try {
                std::stringstream ss(command);
                std::string cmd;
                ss >> cmd;

                if (cmd == "exit") {
                    std::cout << "\033[1;33mGoodbye!\033[0m\n";
                    break;
                }
                else if (cmd == "help") {
                    displayHelp();
                }
                else if (cmd == "clear") {
                    std::cout << "\033[2J\033[H";
                }
                else if (cmd == "ls") {
                    listMatrices();
                }
                else if (cmd == "create") {
                    std::string name;
                    std::string matrix_str;
                    ss >> name;

                    std::string matrixString = extractMatrixString(command);
                    std::getline(ss, matrix_str);
                    //std::cout <<"Matrix: "<<name<< ", Extracted Matrix String:" << matrixString <<",="<<matrix_str<< std::endl;
                    stored_matrices[name] = Matrix::fromString(matrixString);
                    std::cout << "\033[1;32mCreated matrix " << name << ":\033[0m\n";
                    stored_matrices[name].print();
                }
                else if (cmd == "show") {
                    std::string name;
                    ss >> name;
                    if (stored_matrices.count(name)) {
                        stored_matrices[name].print(name);
                    } else {
                        std::cout << "\033[1;31mMatrix " << name << " not found.\033[0m\n";
                    }
                }
                else if (cmd == "del") {
                    std::string name;
                    ss >> name;
                    if (stored_matrices.erase(name)) {
                        std::cout << "\033[1;32mMatrix " << name << " deleted.\033[0m\n";
                    } else {
                        std::cout << "\033[1;31mMatrix " << name << " not found.\033[0m\n";
                    }
                }
                else if (cmd == "det") {
                    std::string A;
                    ss >> A;
                    if (!stored_matrices.count(A)) {
                        throw std::runtime_error("Matrix not found");
                    }
                    double det = stored_matrices[A].determinant();
                    std::cout << "\033[1;32mDeterminant of " << A << " = " << det << "\033[0m\n";
                }
                else if (cmd == "add" || cmd == "sub" || cmd == "mul") {
                    std::string A, B, C;
                    ss >> A >> B >> C;
                    if (!stored_matrices.count(A) || !stored_matrices.count(B)) {
                        throw std::runtime_error("One or both input matrices not found");
                    }
                    if (cmd == "add") {
                        stored_matrices[C] = stored_matrices[A] + stored_matrices[B];
                    }
                    else if (cmd == "sub") {
                        stored_matrices[C] = stored_matrices[A] - stored_matrices[B];
                    }
                    else { // mul
                        stored_matrices[C] = stored_matrices[A] * stored_matrices[B];
                    }
                    std::cout << "\033[1;32mResult stored in " << C << ":\033[0m\n";
                    stored_matrices[C].print();
                }
                else if (cmd == "inv") {
                    std::string A, B;
                    ss >> A >> B;
                    if (!stored_matrices.count(A)) {
                        throw std::runtime_error("Matrix not found");
                    }
                    stored_matrices[B] = stored_matrices[A].inverse();
                    std::cout << "\033[1;32mInverse stored in " << B << ":\033[0m\n";
                    stored_matrices[B].print();
                }
                else {
                    std::cout << "\033[1;31mUnknown command. Type 'help' for available commands.\033[0m\n";
                }
            }
            catch (const std::exception& e) {
                std::cout << "\033[1;31mError: " << e.what() << "\033[0m\n";
            }
        }
    }
};

std::string extractMatrixString(const std::string& command) {
    // Find the first opening '[' and last closing ']'
    std::size_t start = command.find('[');
    std::size_t end = command.find_last_of(']');

    // Check if both brackets exist and are correctly positioned
    if (start == std::string::npos || end == std::string::npos || start >= end) {
        throw std::invalid_argument("Invalid command: Matrix part not found or incorrectly formatted");
    }

    // Extract substring enclosed in brackets
    return command.substr(start, end - start + 1);
}

int main() {
    MatrixCalculator calc;
    calc.run();
    return 0;
}