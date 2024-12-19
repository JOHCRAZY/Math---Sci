#include <iostream>
#include <random>
#include <vector>
#include <cmath>
#include <algorithm>
#include <chrono>
#include <functional>

class RandomNumberSimulator {
private:
    std::mt19937 generator;
    
public:
    // Initialize with random seed based on current time
    RandomNumberSimulator() {
        auto seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
        generator.seed(seed);
    }
    
    // Initialize with specific seed
    RandomNumberSimulator(unsigned int seed) {
        generator.seed(seed);
    }
    
    // Uniform distribution in range [min, max]
    std::vector<double> uniformDistribution(double min, double max, int count) {
        std::uniform_real_distribution<double> distribution(min, max);
        std::vector<double> numbers;
        numbers.reserve(count);
        
        for (int i = 0; i < count; ++i) {
            numbers.push_back(distribution(generator));
        }
        
        return numbers;
    }
    
    // Normal (Gaussian) distribution with mean and standard deviation
    std::vector<double> normalDistribution(double mean, double stddev, int count) {
        std::normal_distribution<double> distribution(mean, stddev);
        std::vector<double> numbers;
        numbers.reserve(count);
        
        for (int i = 0; i < count; ++i) {
            numbers.push_back(distribution(generator));
        }
        
        return numbers;
    }
    
    // Poisson distribution with given mean
    std::vector<int> poissonDistribution(double mean, int count) {
        std::poisson_distribution<int> distribution(mean);
        std::vector<int> numbers;
        numbers.reserve(count);
        
        for (int i = 0; i < count; ++i) {
            numbers.push_back(distribution(generator));
        }
        
        return numbers;
    }
    
    // Exponential distribution with given lambda (rate parameter)
    std::vector<double> exponentialDistribution(double lambda, int count) {
        std::exponential_distribution<double> distribution(lambda);
        std::vector<double> numbers;
        numbers.reserve(count);
        
        for (int i = 0; i < count; ++i) {
            numbers.push_back(distribution(generator));
        }
        
        return numbers;
    }
    
    // Binomial distribution with n trials and probability p
    std::vector<int> binomialDistribution(int n, double p, int count) {
        std::binomial_distribution<int> distribution(n, p);
        std::vector<int> numbers;
        numbers.reserve(count);
        
        for (int i = 0; i < count; ++i) {
            numbers.push_back(distribution(generator));
        }
        
        return numbers;
    }
    
    // Custom distribution using rejection sampling
    std::vector<double> customDistribution(
        std::function<double(double)> pdf,  // Probability density function
        double xMin,                        // Minimum x value
        double xMax,                        // Maximum x value
        double yMax,                        // Maximum y value of PDF
        int count                           // Number of samples
    ) {
        std::vector<double> numbers;
        numbers.reserve(count);
        std::uniform_real_distribution<double> xDist(xMin, xMax);
        std::uniform_real_distribution<double> yDist(0, yMax);
        
        while (numbers.size() < count) {
            double x = xDist(generator);
            double y = yDist(generator);
            
            if (y <= pdf(x)) {
                numbers.push_back(x);
            }
        }
        
        return numbers;
    }
    
    // Calculate basic statistics
    struct Statistics {
        double mean;
        double median;
        double stddev;
        double min;
        double max;
    };
    
    template<typename T>
    Statistics calculateStatistics(const std::vector<T>& numbers) {
        if (numbers.empty()) {
            return {0, 0, 0, 0, 0};
        }
        
        // Calculate mean
        double sum = 0;
        for (const auto& num : numbers) {
            sum += num;
        }
        double mean = sum / numbers.size();
        
        // Calculate standard deviation
        double sqSum = 0;
        for (const auto& num : numbers) {
            sqSum += (num - mean) * (num - mean);
        }
        double stddev = std::sqrt(sqSum / numbers.size());
        
        // Get min and max
        auto [minIt, maxIt] = std::minmax_element(numbers.begin(), numbers.end());
        
        // Calculate median
        std::vector<T> sorted = numbers;
        std::sort(sorted.begin(), sorted.end());
        double median;
        if (sorted.size() % 2 == 0) {
            median = (sorted[sorted.size()/2 - 1] + sorted[sorted.size()/2]) / 2.0;
        } else {
            median = sorted[sorted.size()/2];
        }
        
        return {mean, median, stddev, static_cast<double>(*minIt), static_cast<double>(*maxIt)};
    }
};

// Example usage
int main() {
    RandomNumberSimulator rng;
    
    // Generate numbers from different distributions
    auto uniform = rng.uniformDistribution(0, 1, 1000);
    auto normal = rng.normalDistribution(0, 1, 1000);
    auto poisson = rng.poissonDistribution(5, 1000);
    auto exponential = rng.exponentialDistribution(1.0, 1000);
    auto binomial = rng.binomialDistribution(10, 0.5, 1000);
    
    // Custom distribution example (triangle distribution)
    auto trianglePDF = [](double x) {
        if (x < 0 || x > 1) return 0.0;
        return 2.0 * (x < 0.5 ? x : 1.0 - x);
    };
    auto custom = rng.customDistribution(trianglePDF, 0, 1, 1, 1000);
    
    // Calculate and print statistics for normal distribution
    auto stats = rng.calculateStatistics(normal);
    std::cout << "Normal Distribution Statistics:\n";
    std::cout << "Mean: " << stats.mean << "\n";
    std::cout << "Median: " << stats.median << "\n";
    std::cout << "StdDev: " << stats.stddev << "\n";
    std::cout << "Min: " << stats.min << "\n";
    std::cout << "Max: " << stats.max << "\n";
    
    return 0;
}