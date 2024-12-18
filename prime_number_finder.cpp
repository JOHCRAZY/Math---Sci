#include <iostream>
#include <vector>

void sieveOfEratosthenes(int limit) {
    // Step 1: Initialize a boolean vector and a result container
    std::vector<bool> isPrime(limit + 1, true); // Mark all as true initially
    std::vector<int> primes; // Vector to store prime numbers

    isPrime[0] = isPrime[1] = false; // 0 and 1 are not prime numbers

    // Step 2: Sieve logic
    for (int p = 2; p * p <= limit; ++p) {
        if (isPrime[p]) {
            // Mark all multiples of p as not prime
            for (int i = p * p; i <= limit; i += p) {
                isPrime[i] = false;
            }
        }
    }

    // Step 3: Collect all prime numbers into the vector
    for (int i = 2; i <= limit; ++i) {
        if (isPrime[i]) {
            primes.push_back(i);
        }
    }

    // Step 4: Output the prime numbers
    std::cout << "Prime numbers up to " << limit << " are:" << std::endl;
    for (int prime : primes) {
        std::cout << prime << " ";
    }
    std::cout << std::endl;
}

int main() {
    int limit;

    // User input for the limit
    std::cout << "Enter the limit up to which you want to find prime numbers: ";
    std::cin >> limit;

    // Call Sieve of Eratosthenes
    if (limit >= 2) {
        sieveOfEratosthenes(limit);
    } else {
        std::cout << "There are no prime numbers less than 2." << std::endl;
    }

    return 0;
}
