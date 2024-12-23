#include <iostream>
#include <random>
#include <vector>
#include <cmath>
#include <algorithm>
#include <chrono>
#include <functional>
#include <fstream>
#include <string>
#include <climits>
#include <iomanip>
#include <cfloat>
#include <variant>


using SampleVector = std::variant<std::vector<int>, std::vector<double>>;


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
    
            // Cauchy distribution with location and scale parameters
            std::vector<double> cauchyDistribution(double location, double scale, int count) {
                std::cauchy_distribution<double> distribution(location, scale);
                std::vector<double> numbers;
                numbers.reserve(count);
                
                for (int i = 0; i < count; ++i) {
                    numbers.push_back(distribution(generator));
                }
                return numbers;
            }
            
            // Chi-squared distribution with degrees of freedom
            std::vector<double> chiSquaredDistribution(double degreesOfFreedom, int count) {
                std::chi_squared_distribution<double> distribution(degreesOfFreedom);
                std::vector<double> numbers;
                numbers.reserve(count);
                
                for (int i = 0; i < count; ++i) {
                    numbers.push_back(distribution(generator));
                }
                return numbers;
            }
            
            // Student's t-distribution with degrees of freedom
            std::vector<double> studentTDistribution(double degreesOfFreedom, int count) {
                std::student_t_distribution<double> distribution(degreesOfFreedom);
                std::vector<double> numbers;
                numbers.reserve(count);
                
                for (int i = 0; i < count; ++i) {
                    numbers.push_back(distribution(generator));
                }
                return numbers;
            }

            // Generalized Student's t-distribution with location and scale
std::vector<double> studentTDistribution(double degreesOfFreedom, double location, double scale, int count) {
    std::student_t_distribution<double> distribution(degreesOfFreedom);
    std::vector<double> numbers;
    numbers.reserve(count);
    
    for (int i = 0; i < count; ++i) {
        // Transform standard t-distribution to get location and scale
        numbers.push_back(location + scale * distribution(generator));
    }
    return numbers;
}
            
            // Fisher F-distribution with d1 and d2 degrees of freedom
            std::vector<double> fisherFDistribution(double d1, double d2, int count) {
                std::fisher_f_distribution<double> distribution(d1, d2);
                std::vector<double> numbers;
                numbers.reserve(count);
                
                for (int i = 0; i < count; ++i) {
                    numbers.push_back(distribution(generator));
                }
                return numbers;
            }
            
            // Lognormal distribution
            std::vector<double> lognormalDistribution(double m, double s, int count) {
                std::lognormal_distribution<double> distribution(m, s);
                std::vector<double> numbers;
                numbers.reserve(count);
                
                for (int i = 0; i < count; ++i) {
                    numbers.push_back(distribution(generator));
                }
                return numbers;
            }
            
            // Gamma distribution with shape and scale parameters
            std::vector<double> gammaDistribution(double alpha, double beta, int count) {
                std::gamma_distribution<double> distribution(alpha, beta);
                std::vector<double> numbers;
                numbers.reserve(count);
                
                for (int i = 0; i < count; ++i) {
                    numbers.push_back(distribution(generator));
                }
                return numbers;
            }
            
            // Weibull distribution
            std::vector<double> weibullDistribution(double shape, double scale, int count) {
                std::weibull_distribution<double> distribution(shape, scale);
                std::vector<double> numbers;
                numbers.reserve(count);
                
                for (int i = 0; i < count; ++i) {
                    numbers.push_back(distribution(generator));
                }
                return numbers;
            }
            
            
            // Geometric distribution with probability p
            std::vector<int> geometricDistribution(double p, int count) {
                std::geometric_distribution<int> distribution(p);
                std::vector<int> numbers;
                numbers.reserve(count);
                
                for (int i = 0; i < count; ++i) {
                    numbers.push_back(distribution(generator));
                }
                return numbers;
            }
            
             // Bernoulli distribution with probability p
             std::vector<int> bernoulliDistribution(double p, int count) {
                std::bernoulli_distribution distribution(p);
                std::vector<int> numbers;
                numbers.reserve(count);
                
                for (int i = 0; i < count; ++i) {
                    numbers.push_back(distribution(generator));
                }
                return numbers;
            }

            // Negative binomial distribution
            std::vector<int> negativeBinomialDistribution(int r, double p, int count) {
                std::negative_binomial_distribution<int> distribution(r, p);
                std::vector<int> numbers;
                numbers.reserve(count);
                
                for (int i = 0; i < count; ++i) {
                    numbers.push_back(distribution(generator));
                }
                return numbers;
            }
            
            // Discrete uniform distribution in range [min, max]
            std::vector<int> discreteUniformDistribution(int min, int max, int count) {
                std::uniform_int_distribution<int> distribution(min, max);
                std::vector<int> numbers;
                numbers.reserve(count);
                
                for (int i = 0; i < count; ++i) {
                    numbers.push_back(distribution(generator));
                }
                return numbers;
            }
        
            // Custom distributions
            
            // Triangular distribution with min, peak, and max values
            std::vector<double> triangularDistribution(double min, double peak, double max, int count) {
                std::vector<double> numbers;
                numbers.reserve(count);
                
                for (int i = 0; i < count; ++i) {
                    double u = std::uniform_real_distribution<double>(0, 1)(generator);
                    double f = (peak - min) / (max - min);
                    
                    if (u < f) {
                        numbers.push_back(min + std::sqrt(u * (max - min) * (peak - min)));
                    } else {
                        numbers.push_back(max - std::sqrt((1 - u) * (max - min) * (max - peak)));
                    }
                }
                return numbers;
            }
            
            // Mixture of two normal distributions
            std::vector<double> mixtureNormalDistribution(
                double mean1, double stddev1, double weight1,
                double mean2, double stddev2, int count) {
                
                std::normal_distribution<double> dist1(mean1, stddev1);
                std::normal_distribution<double> dist2(mean2, stddev2);
                std::uniform_real_distribution<double> uniform(0, 1);
                
                std::vector<double> numbers;
                numbers.reserve(count);
                
                for (int i = 0; i < count; ++i) {
                    if (uniform(generator) < weight1) {
                        numbers.push_back(dist1(generator));
                    } else {
                        numbers.push_back(dist2(generator));
                    }
                }
                return numbers;
            }
        
            // Utility function to reseed the generator
            void reseed(unsigned int seed) {
                generator.seed(seed);
            }
    
    // Calculate basic statistics
    struct Statistics {
        double mean;
        double median;
        double stddev;
        double min;
        double max;
        double skewness;
        double kurtosis;
    };
    
    template<typename T>
    Statistics calculateStatistics(const std::vector<T>& numbers) {
        if (numbers.empty()) {
            return {0, 0, 0, 0, 0, 0, 0};
        }
        
        // Calculate mean
        double sum = 0;
        for (const auto& num : numbers) {
            sum += num;
        }
        double mean = sum / numbers.size();
        
        // Calculate standard deviation and prepare for skewness/kurtosis
        double sqSum = 0;
        double cubSum = 0;  // For skewness
        double quadSum = 0; // For kurtosis
        
        for (const auto& num : numbers) {
            double diff = num - mean;
            double sqDiff = diff * diff;
            sqSum += sqDiff;
            cubSum += sqDiff * diff;
            quadSum += sqDiff * sqDiff;
        }
        
        double variance = sqSum / numbers.size();
        double stddev = std::sqrt(variance);
        
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
        
        // Calculate skewness
        // Skewness = E[(X-μ)³] / σ³
        double skewness = 0;
        if (stddev > 0) {
            skewness = (cubSum / numbers.size()) / (stddev * stddev * stddev);
        }
        
        // Calculate kurtosis
        // Kurtosis = E[(X-μ)⁴] / σ⁴ - 3
        // Note: We subtract 3 to make the kurtosis of a normal distribution = 0
        double kurtosis = 0;
        if (stddev > 0) {
            kurtosis = (quadSum / numbers.size()) / (variance * variance) - 3.0;
        }
        
        return {
            mean,
            median,
            stddev,
            static_cast<double>(*minIt),
            static_cast<double>(*maxIt),
            skewness,
            kurtosis
        };
    }
};


enum class DistributionType {
    Uniform,
    DiscreteUniform,
    Normal,
    Poisson,
    Exponential,
    Binomial,
    NegativeBinomial,
    Bernoulli,
    Cauchy,
    ChiSquared,
    Gamma,
    Geometric,
    Triangular,
    MixtureNormal,
    StudentT,
    Weibull,
    Lognormal
};

struct DistributionParams {
    std::vector<double> params;
    int sampleSize;
};

class RandomNumberSimulatorUI : public RandomNumberSimulator {
private:
    RandomNumberSimulator rng;

    void displayMenu() {
      std::cout << "\n╔════════════ Random Number Generator ════════════╗\n";
        std::cout << "║ Available Distributions:                        ║\n";
        std::cout << "║  1. Uniform                                     ║\n";
        std::cout << "║  2. Discrete Uniform                            ║\n";
        std::cout << "║  3. Normal (Gaussian)                           ║\n";
        std::cout << "║  4. Poisson                                     ║\n";
        std::cout << "║  5. Exponential                                 ║\n";
        std::cout << "║  6. Binomial                                    ║\n";
        std::cout << "║  7. Negative Binomial                           ║\n";
        std::cout << "║  8. Bernoulli                                   ║\n";
        std::cout << "║  9. Cauchy                                      ║\n";
        std::cout << "║ 10. Chi-Squared                                 ║\n";
        std::cout << "║ 11. Gamma                                       ║\n";
        std::cout << "║ 12. Geometric                                   ║\n";
        std::cout << "║ 13. Triangular                                  ║\n";
        std::cout << "║ 14. Mixture Normal                              ║\n";
        std::cout << "║ 15. Student's t                                 ║\n";
        std::cout << "║ 16. Weibull                                     ║\n";
        std::cout << "║ 17. Lognormal                                   ║\n";
        std::cout << "║  0. Exit                                        ║\n";
        std::cout << "╚═════════════════════════════════════════════════╝\n";
    }

    void printVector(const std::vector<double>& vec, int precision = 1, int width = 5) {
        // Input validation
        if (precision < 0 || width < 1) {
            throw std::invalid_argument("Invalid precision or width parameters");
        }
    
        // Save original stream state
        auto oldFlags = std::cout.flags();
        auto oldPrecision = std::cout.precision();
        
        try {
            // Set stream formatting
            std::cout << std::fixed << std::setprecision(precision);
            
            std::cout << "[";
            
            // Handle empty vector case
            if (vec.empty()) {
                std::cout << "]" << std::endl;
                return;
            }
    
            // Print first element without leading comma
            std::cout << std::setw(width) << vec.front();
            
            // Print remaining elements with leading commas
            for (auto it = std::next(vec.begin()); it != vec.end(); ++it) {
                std::cout << ",   " << std::setw(width) <<std::fixed<< *it;
                
                // Add line break for readability if vector is long
                if (std::distance(vec.begin(), it) % 18 == 0) {
                    std::cout << "\n";
                }
            }
            
            std::cout << "]" << std::endl;
        }
        catch (...) {
            // Restore stream state before re-throwing
            std::cout.flags(oldFlags);
            std::cout.precision(oldPrecision);
            throw;
        }
        
        // Restore original stream state
        std::cout.flags(oldFlags);
        std::cout.precision(oldPrecision);
    }

    double getValidatedDouble(const std::string& prompt, double min = -DBL_MAX, double max = DBL_MAX) {
        while (true) {
            std::cout << prompt;
            std::string input;
            std::getline(std::cin, input);
            
            try {
                double value = std::stod(input);
                if (value >= min && value <= max) {
                    return value;
                }
                std::cout << "Error: Value must be between " << min << " and " << max << "\n";
            }
            catch (...) {
                std::cout << "Error: Please enter a valid number\n";
            }
        }
    }

    int getValidatedInt(const std::string& prompt, int min = INT_MIN, int max = INT_MAX) {
        while (true) {
            std::cout << prompt;
            std::string input;
            std::getline(std::cin, input);
            
            try {
                int value = std::stoi(input);
                if (value >= min && value <= max) {
                    return value;
                }
                std::cout << "Error: Value must be between " << min << " and " << max << "\n";
            }
            catch (...) {
                std::cout << "Error: Please enter a valid integer\n";
            }
        }
    }

    DistributionParams getDistributionParameters(DistributionType type) {
        DistributionParams params;
        
        // Get sample size for all distributions
        params.sampleSize = getValidatedInt("Enter sample size (1-1000000): ", 1, 1000000);

        const double EPSILON = 0.000001; // Minimum value for positive parameters
    
        switch (type) {
            case DistributionType::Uniform:
                params.params.push_back(getValidatedDouble("Enter minimum value: "));
                params.params.push_back(getValidatedDouble("Enter maximum value: ", params.params[0]));
                if (params.params[0] >= params.params[1]) {
                    throw std::runtime_error("Uniform distribution: minimum must be less than maximum");
                }
                break;

            case DistributionType::DiscreteUniform:
                params.params.push_back(getValidatedInt("Enter minimum value: "));
                params.params.push_back(getValidatedInt("Enter maximum value: ", params.params[0]));
                if (params.params[0] >= params.params[1]) {
                    throw std::runtime_error("Discrete Uniform distribution: minimum must be less than maximum");
                }
                break;
    
            case DistributionType::Normal:
                params.params.push_back(getValidatedDouble("Enter mean: "));
                params.params.push_back(getValidatedDouble("Enter standard deviation (>0): ", EPSILON));
                break;
    
            case DistributionType::Poisson:
                params.params.push_back(getValidatedDouble("Enter mean (lambda > 0): ", EPSILON));
                break;
    
            case DistributionType::Exponential:
                params.params.push_back(getValidatedDouble("Enter rate parameter (lambda > 0): ", EPSILON));
                break;
    
            case DistributionType::Binomial:
                params.params.push_back(getValidatedInt("Enter number of trials (n > 0): ", 1));
                params.params.push_back(getValidatedDouble("Enter probability (0-1): ", 0, 1));
                break;
            
            case DistributionType::NegativeBinomial:
                params.params.push_back(getValidatedInt("Enter number of success (r > 0): ", 1));
                params.params.push_back(getValidatedDouble("Enter probability (0-1): ", 0, 1));
                break;
    
            case DistributionType::Cauchy:
                params.params.push_back(getValidatedDouble("Enter location parameter: "));
                params.params.push_back(getValidatedDouble("Enter scale parameter (>0): ", EPSILON));
                break;
    
            case DistributionType::ChiSquared:
                params.params.push_back(getValidatedDouble("Enter degrees of freedom (>0): ", EPSILON));
                break;
    
            case DistributionType::Gamma:
                params.params.push_back(getValidatedDouble("Enter shape parameter (k > 0): ", EPSILON));
                params.params.push_back(getValidatedDouble("Enter scale parameter (theta > 0): ", EPSILON));
                break;
    
            case DistributionType::Geometric:
                params.params.push_back(getValidatedDouble("Enter success probability (0-1): ", 0, 1));
                break;

            case DistributionType::Bernoulli:
                params.params.push_back(getValidatedDouble("Enter success probability (0-1): ", 0, 1));
                break;
    
            case DistributionType::Triangular:
                params.params.push_back(getValidatedDouble("Enter minimum value (a): "));
                params.params.push_back(getValidatedDouble("Enter most likely value (c): ", params.params[0]));
                params.params.push_back(getValidatedDouble("Enter maximum value (b): ", params.params[1]));
                if (params.params[0] >= params.params[1] || params.params[1] >= params.params[2]) {
                    throw std::runtime_error("Triangular distribution: must satisfy a < c < b");
                }
                break;
    
            case DistributionType::MixtureNormal:
                // For a mixture of two normal distributions
                params.params.push_back(getValidatedDouble("Enter mean of first component: "));
                params.params.push_back(getValidatedDouble("Enter std dev of first component (>0): ", EPSILON));
                params.params.push_back(getValidatedDouble("Enter weight of first component (0-1): ", 0, 1));

                params.params.push_back(getValidatedDouble("Enter mean of second component: "));
                params.params.push_back(getValidatedDouble("Enter std dev of second component (>0): ", EPSILON));
                break;
    
            case DistributionType::StudentT:
                params.params.push_back(getValidatedDouble("Enter degrees of freedom (>0): ", EPSILON));
                params.params.push_back(getValidatedDouble("Enter location parameter: "));
                params.params.push_back(getValidatedDouble("Enter scale parameter (>0): ", EPSILON));
                break;
    
            case DistributionType::Weibull:
                params.params.push_back(getValidatedDouble("Enter shape parameter (k > 0): ", EPSILON));
                params.params.push_back(getValidatedDouble("Enter scale parameter (lambda > 0): ", EPSILON));
                break;
    
            case DistributionType::Lognormal:
                params.params.push_back(getValidatedDouble("Enter log-scale parameter (mu): "));
                params.params.push_back(getValidatedDouble("Enter shape parameter (sigma > 0): ", EPSILON));
                break;
    
            default:
                throw std::runtime_error("Unsupported distribution type");
        }
        
        return params;

    }

    std::vector<double> generateSamples(DistributionType type, const DistributionParams& params) {
        switch (type) {
            case DistributionType::Uniform:
                std::cout<<"Uniform Distribution: "<<std::endl;
                return rng.uniformDistribution(params.params[0], params.params[1], params.sampleSize);

            case DistributionType::Normal:
               std::cout<<"Normal Distribution: "<<std::endl;
                 return rng.normalDistribution(params.params[0], params.params[1], params.sampleSize);

            case DistributionType::DiscreteUniform:{
                std::cout<<"Discrete Uniform Distribution: "<<std::endl;
                auto samples = rng.discreteUniformDistribution(params.params[0], params.params[1], params.sampleSize);
                return std::vector<double>(samples.begin(), samples.end());
        }

        case DistributionType::Poisson:{
            std::cout<<"Poison Distribution: "<<std::endl;
            auto samples = rng.poissonDistribution(params.params[0], params.sampleSize);
            return std::vector<double>(samples.begin(), samples.end());
        }

        case DistributionType::Exponential:
        std::cout<<"Exponential Distribution: "<<std::endl;
            return rng.exponentialDistribution(params.params[0], params.sampleSize);

        

        case DistributionType::Triangular:
        std::cout<<"Triangular Distribution: "<<std::endl;
                return rng.triangularDistribution(params.params[0], params.params[1], params.params[2], params.sampleSize);
        
        case DistributionType::Bernoulli:{
            auto samples =  rng.bernoulliDistribution(params.params[0], params.sampleSize);
            return std::vector<double>(samples.begin(), samples.end());
        }

        case DistributionType::Binomial:{
            std::cout<<"Binomial Distribution: "<<std::endl;
            auto samples =  rng.binomialDistribution(params.params[0], params.params[1], params.sampleSize);
            return std::vector<double>(samples.begin(), samples.end());
        }
        case DistributionType::NegativeBinomial:{
            std::cout<<"Negative Binomial Distribution: "<<std::endl;
            auto samples =  rng.negativeBinomialDistribution(params.params[0], params.params[1], params.sampleSize);
                return std::vector<double>(samples.begin(), samples.end());
        }

        case DistributionType::Cauchy:{
            std::cout<<"Cauchy Distribution: "<<std::endl;
            return rng.cauchyDistribution(params.params[0], params.params[1], params.sampleSize);
        }

        case DistributionType::ChiSquared:{
            std::cout<<"Chi Squared Distribution: "<<std::endl;
            return rng.chiSquaredDistribution(params.params[0], params.sampleSize);
        }

        case DistributionType::Gamma:{
            std::cout<<"Gamma Distribution: "<<std::endl;
            return rng.gammaDistribution(params.params[0], params.params[1], params.sampleSize);
        }

        case DistributionType::Geometric:{
            std::cout<<"Geometric Distribution: "<<std::endl;
            auto samples =  rng.geometricDistribution(params.params[0], params.sampleSize);
            return std::vector<double>(samples.begin(), samples.end());
        }

        case DistributionType::StudentT:{
            std::cout<<"Student T Distribution: "<<std::endl;
            return rng.studentTDistribution(params.params[0], params.params[1], params.params[2],params.sampleSize);
        }

        case DistributionType::MixtureNormal:{
            std::cout<<"Mixture Normal Distribution: "<<std::endl;
            return rng.mixtureNormalDistribution(params.params[0], params.params[1], params.params[2], params.params[3],params.params[4], params.sampleSize);
        }

        case DistributionType::Weibull:{
            std::cout<<"Weibull Distribution: "<<std::endl;
            return rng.weibullDistribution(params.params[0], params.params[1], params.sampleSize);
        }

        case DistributionType::Lognormal:{
            std::cout<<"Lognormal Distribution: "<<std::endl;
            return rng.lognormalDistribution(params.params[0], params.params[1], params.sampleSize);
        }

            default:
                throw std::runtime_error("Unsupported distribution type");
        }
    }

    void displayStatistics(const std::vector<double>& samples) {
        auto stats = calculateStatistics(samples);
        
      std::cout << "\n╔═════════ Statistics ══════════════╗\n";
        std::cout << "║ Mean:      " << std::setw(20) << stats.mean << "   ║\n";
        std::cout << "║ Median:    " << std::setw(20) << stats.median << "   ║\n";
        std::cout << "║ Std Dev:   " << std::setw(20) << stats.stddev << "   ║\n";
        std::cout << "║ Min:       " << std::setw(20) << stats.min << "   ║\n";
        std::cout << "║ Max:       " << std::setw(20) << stats.max << "   ║\n";
        std::cout << "║ Skewness:  " << std::setw(20) << stats.skewness << "   ║\n";
        std::cout << "║ Kurtosis:  " << std::setw(20) << stats.kurtosis << "   ║\n";
        std::cout << "╚═══════════════════════════════════╝\n";
    }

    void saveToFile(const std::vector<double>& samples, const std::string& filename) {
        std::ofstream file(filename);
        if (!file) {
            throw std::runtime_error("Could not open file for writing");
        }
        
        for (const auto& sample : samples) {
            file << sample << "\n";
        }
    }

public:
    void run() {
        while (true) {
            displayMenu();
            
            int choice = getValidatedInt("Enter your choice (0 - 17) $:  ", 0, 17);
            if (choice == 0) break;

            try {
                DistributionType type = static_cast<DistributionType>(choice - 1);
                DistributionParams params = getDistributionParameters(type);
                
                std::cout << "\nGenerating Numbers.....\n";
                auto samples = generateSamples(type, params);
                printVector(samples);



                displayStatistics(samples);

                 // Ask if user wants to save samples
                 std::cout << "\nWould you like to save the generated numbers to a file? (y/n): ";
                 std::string response;
                 std::getline(std::cin, response);
                 
                 if (response == "y" || response == "Y") {
                     std::cout << "Enter filename: ";
                     std::string filename;
                     std::getline(std::cin, filename);
                     saveToFile(samples, filename);
                     std::cout << "Numbers saved to " << filename << "\n";
                 }
                
               
            }
            catch (const std::exception& e) {
                std::cerr << "Error: " << e.what() << "\n";
            }
        }
    }
};

int main() {
    RandomNumberSimulatorUI simulator;
    simulator.run();
    return 0;
}