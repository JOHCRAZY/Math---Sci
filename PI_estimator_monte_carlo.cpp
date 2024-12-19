#include <iostream>
#include <random>
#include <chrono>
#include <iomanip>
#include <vector>
#include <thread>

class PiEstimator {
private:
    std::mt19937_64 generator;
    std::uniform_real_distribution<double> distribution;
    
    // Structure to store results for visualization
    struct Point {
        double x, y;
        bool inside;
    };
    
    // Check if a point lies within the quarter circle
    bool isInsideQuarterCircle(double x, double y) {
        return (x * x + y * y) <= 1.0;
    }
    
    // Single-threaded estimation
    double estimateSingleThread(unsigned long long points) {
        unsigned long long pointsInside = 0;
        
        for (unsigned long long i = 0; i < points; ++i) {
            double x = distribution(generator);
            double y = distribution(generator);
            
            if (isInsideQuarterCircle(x, y)) {
                pointsInside++;
            }
        }
        
        return 4.0 * pointsInside / points;
    }
    
    // Worker function for multi-threaded estimation
    unsigned long long estimateWorker(unsigned long long points) {
        unsigned long long pointsInside = 0;
        std::mt19937_64 localGen(generator());  // Create thread-local generator
        
        for (unsigned long long i = 0; i < points; ++i) {
            double x = distribution(localGen);
            double y = distribution(localGen);
            
            if (isInsideQuarterCircle(x, y)) {
                pointsInside++;
            }
        }
        
        return pointsInside;
    }

public:
    PiEstimator() {
        // Seed the generator with current time
        auto seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
        generator.seed(seed);
        distribution = std::uniform_real_distribution<double>(0.0, 1.0);
    }
    
    // Multi-threaded estimation with progress tracking
    double estimate(unsigned long long totalPoints, unsigned int numThreads = 1, bool showProgress = true) {
        if (numThreads == 1) {
            return estimateSingleThread(totalPoints);
        }
        
        std::vector<std::thread> threads;
        std::vector<unsigned long long> results(numThreads);
        unsigned long long pointsPerThread = totalPoints / numThreads;
        
        // Launch threads
        for (unsigned int i = 0; i < numThreads; ++i) {
            threads.emplace_back([this, i, pointsPerThread, &results]() {
                results[i] = estimateWorker(pointsPerThread);
            });
        }
        
        // Show progress while waiting
        if (showProgress) {
            std::cout << "Estimating Pi using " << numThreads << " threads:\n";
            for (int i = 0; i < 50; ++i) {
                std::cout << "-";
            }
            std::cout << "\n";
        }
        
        // Wait for all threads
        for (auto& thread : threads) {
            thread.join();
        }
        
        // Sum up results
        unsigned long long totalInside = 0;
        for (auto result : results) {
            totalInside += result;
        }
        
        return 4.0 * totalInside / (pointsPerThread * numThreads);
    }
    
    // Get sample points for visualization
    std::vector<Point> getSamplePoints(int numPoints) {
        std::vector<Point> points;
        points.reserve(numPoints);
        
        for (int i = 0; i < numPoints; ++i) {
            double x = distribution(generator);
            double y = distribution(generator);
            bool inside = isInsideQuarterCircle(x, y);
            points.push_back({x, y, inside});
        }
        
        return points;
    }
    
    // ASCII visualization of the estimation
    void visualize(int size = 20) {
        std::vector<std::vector<char>> grid(size, std::vector<char>(size, ' '));
        auto points = getSamplePoints(size * size / 2);
        
        // Plot points
        for (const auto& point : points) {
            int x = static_cast<int>(point.x * (size - 1));
            int y = static_cast<int>(point.y * (size - 1));
            grid[y][x] = point.inside ? '*' : '.';
        }
        
        // Draw the grid
        std::cout << "\nVisualization of Monte Carlo estimation:\n";
        std::cout << std::string(size + 2, '-') << '\n';
        for (int i = size - 1; i >= 0; --i) {
            std::cout << '|';
            for (int j = 0; j < size; ++j) {
                std::cout << grid[i][j];
            }
            std::cout << "|\n";
        }
        std::cout << std::string(size + 2, '-') << '\n';
    }
};

int main() {
    PiEstimator estimator;
    
    // Parameters
    unsigned long long numPoints = 1000000000;  // 1 billion points
    unsigned int numThreads = std::thread::hardware_concurrency();  // Use all available cores
    
    // Estimate Pi
    std::cout << "Estimating Pi using " << numPoints << " points...\n";
    auto start = std::chrono::high_resolution_clock::now();
    
    double estimatedPi = estimator.estimate(numPoints, numThreads);
    
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    
    // Show results
    std::cout << "\nResults:\n";
    std::cout << "Estimated Pi: " << std::setprecision(10) << estimatedPi << "\n";
    std::cout << "Actual Pi:    " << std::setprecision(10) << M_PI << "\n";
    std::cout << "Error:        " << std::setprecision(10) << std::abs(estimatedPi - M_PI) << "\n";
    std::cout << "Time taken:   " << duration.count() / 1000.0 << " seconds\n\n";
    
    // Show visualization
    estimator.visualize();
    
    return 0;
}