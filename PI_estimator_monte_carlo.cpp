#include <iostream>
#include <random>
#include <chrono>
#include <iomanip>
#include <vector>
#include <thread>
#include <mutex>
#include <condition_variable>

using namespace std;

class PiEstimator {
private:
    mt19937_64 generator;
    uniform_real_distribution<double> distribution;
    
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
        auto seed = chrono::high_resolution_clock::now().time_since_epoch().count();
        generator.seed(seed);
        distribution = uniform_real_distribution<double>(0.0, 1.0);
    }
    
    // Multi-threaded estimation with progress tracking
    double estimate(unsigned long long totalPoints, unsigned int numThreads = 1, bool showProgress = true) {
        if (numThreads == 1) {
            return estimateSingleThread(totalPoints);
        }
        
        vector<std::thread> threads;
        vector<unsigned long long> results(numThreads);
        unsigned long long pointsPerThread = totalPoints / numThreads;
        
        // Launch threads
        for (unsigned int i = 0; i < numThreads; ++i) {
            threads.emplace_back([this, i, pointsPerThread, &results]() {
                results[i] = estimateWorker(pointsPerThread);
            });
        }
        
        // Show progress while waiting
        if (showProgress) {
            cout << "Estimating Pi using " << numThreads << " threads:\n";
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
    
};

class Loading {
    private:
        bool running;
        thread animationThread;
        mutex mtx;
        condition_variable cv;
        const vector<std::string> frames = {
            "⠋", "⠙", "⠹", "⠸", "⠼", "⠴", "⠦", "⠧", "⠇", "⠏"
        };
        
    public:
        Loading() : running(false) {}
        
        void start(const string& message) {
            running = true;
            animationThread = thread([this, message]() {
                int frame = 0;
                auto lastUpdate = chrono::steady_clock::now();
                
                while (running) {
                    auto now = chrono::steady_clock::now();
                    if (chrono::duration_cast<chrono::milliseconds>(now - lastUpdate).count() >= 80) {
                        // Clear the current line
                        cout << "\r" << string(80, ' ') << "\r";
                        cout << frames[frame] << " " << message << flush;
                        
                        frame = (frame + 1) % frames.size();
                        lastUpdate = now;
                    }
                    
                    // Sleep to prevent high CPU usage
                    std::this_thread::sleep_for(std::chrono::milliseconds(10));
                }
            });
        }
        
        void stop() {
            running = false;
            if (animationThread.joinable()) {
                animationThread.join();
            }
            // Clear the animation line
            cout << "\r" << string(80, ' ') << "\r" << flush;
        }
        
        ~Loading() {
            stop();
        }
    };

// Function to format large numbers with commas
string formatNumber(unsigned long long n) {
    string number = to_string(n);
    int insertPosition = number.length() - 3;
    while (insertPosition > 0) {
        number.insert(insertPosition, ",");
        insertPosition -= 3;
    }
    return number;
}


    void calculatePi(unsigned long long numPoints = 1000000000) {
        PiEstimator estimator;
        unsigned int numThreads = thread::hardware_concurrency();
        
        
        cout << "\n╔══════════════════════════════════════════╗" << std::endl;
        cout << "║         Monte Carlo Pi Estimation        ║" << std::endl;
        cout << "╚══════════════════════════════════════════╝\n" << std::endl;

        cout<<"Enter number of points: ";
        cin>>numPoints;

        // Format the number of points with commas
        string formattedPoints = formatNumber(numPoints);
        
        cout << "Configuration:" << std::endl;
        cout << "• Points:  " << formattedPoints << std::endl;
        cout << "• Threads: " << numThreads << std::endl;
        
        // Create and start the loading animation
        Loading loader;
        std::string loadingMessage = "Estimating Pi using " + formattedPoints + " points...";
        loader.start(loadingMessage);
        
        auto start = std::chrono::high_resolution_clock::now();
        
        // Perform the estimation
        double estimatedPi = estimator.estimate(numPoints, numThreads);
        
        // Stop timing
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
        
        loader.stop();
        
        // Calculate error
        double error = abs(estimatedPi - M_PI);
        double errorPercentage = (error / M_PI) * 100;
        
        // Display results with a nice format
        cout << "\nResults:" << endl;
        cout << "╭──────────────────── ──────────────────────╮" << endl;
        cout << "│ Estimated π: " << setw(28) << setprecision(10) << fixed << estimatedPi << " │" << endl;
        cout << "│ Actual π:    " << setw(28) << M_PI << " │" << endl;
        cout << "│ Error:       " << setw(28) << error << " │" << endl;
        cout << "│ Error (%):   " << setw(27) << setprecision(4) << fixed  << errorPercentage << "% │" << std::endl;
        cout << "│ Time taken:        " << setw(14) << setprecision(2) << fixed << duration.count() / 1000.0 << " seconds │" << endl;
        cout << "╰──────────────────── ──────────────────────╯" << endl;
    }


int main() {
    calculatePi();
    return 0;
}