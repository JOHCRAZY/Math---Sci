#include <iostream>
#include <string>
#include <algorithm>
#include <set>
#include <sstream>
#include <vector>
#include <iomanip>
#include <map>

using namespace std;

class SET {
private:
    set<string> elements;
    set<int> numericElements;
    
    // Helper method to split string by delimiter
    vector<string> split(const string& str, char delim) {
        vector<string> tokens;
        stringstream ss(str);
        string token;
        while (getline(ss, token, delim)) {
            if (!token.empty()) {
                tokens.push_back(token);
            }
        }
        return tokens;
    }

    

public:
// Helper method to trim whitespace
string trim(const string& str) {
    size_t first = str.find_first_not_of(" \t\n\r");
    size_t last = str.find_last_not_of(" \t\n\r");
    if (first == string::npos || last == string::npos) return "";
    return str.substr(first, last - first + 1);
}
    set<char> setNames;  // Track available set names

    // Enhanced element insertion with error handling
    bool insertElements(const string& setStr) {
        try {
            string elementList = setStr.substr(setStr.find('=') + 1);
            vector<string> tokens = split(elementList, ',');
            
            elements.clear();
            numericElements.clear();
            
            for (const auto& token : tokens) {
                string trimmed = trim(token);
                if (!trimmed.empty()) {
                    elements.insert(trimmed);
                    // Try to convert to number if possible
                    try {
                        numericElements.insert(stoi(trimmed));
                    } catch (...) {
                        // Not a number, keep as string only
                    }
                }
            }
            return true;
        } catch (...) {
            cout << "\033[1;31mError: Invalid input format\033[0m\n";
            return false;
        }
    }

    // Set operations
    SET setUnion(const SET& other) const {
        SET result;
        set_union(
            numericElements.begin(), numericElements.end(),
            other.numericElements.begin(), other.numericElements.end(),
            inserter(result.numericElements, result.numericElements.begin())
        );
        return result;
    }

    SET setIntersection(const SET& other) const {
        SET result;
        set_intersection(
            numericElements.begin(), numericElements.end(),
            other.numericElements.begin(), other.numericElements.end(),
            inserter(result.numericElements, result.numericElements.begin())
        );
        return result;
    }

    SET setDifference(const SET& other) const {
        SET result;
        set_difference(
            numericElements.begin(), numericElements.end(),
            other.numericElements.begin(), other.numericElements.end(),
            inserter(result.numericElements, result.numericElements.begin())
        );
        return result;
    }

    // Enhanced display methods
    void displaySet(const string& name) const {
        cout << "\033[1;33m" << name << " = {";
        bool first = true;
        for (const auto& elem : elements) {
            if (!first) cout << ", ";
            cout << elem;
            first = false;
        }
        cout << "}\033[0m\n";
    }

    void displayOperation(const string& op, const SET& result) const {
        cout << "\033[1;36m" << op << " = {";
        bool first = true;
        for (const auto& elem : result.numericElements) {
            if (!first) cout << ", ";
            cout << elem;
            first = false;
        }
        cout << "}\033[0m\n";
    }

    // Enhanced prompt with better formatting
    static void displayPrompt(const set<char>& setNames) {
        cout << "\033[1;34m";
        for (char name : setNames) {
            cout << name;
            if (name != *setNames.rbegin()) cout << ":";
        }
        cout << "\033[0m\033[1;32m$ \033[0m";
    }

    // Enhanced help menu
    static void displayHelp() {
        cout << "\033[1;35m╔════════════════════════ SET OPERATIONS HELP ═══════════════════════╗\n";
        cout << "║ Commands:                                                          ║\n";
        cout << "║  1. Define set:    A=1,2,3,4                                       ║\n";
        cout << "║  2. Display set:   A?                                              ║\n";
        cout << "║  3. Union:         A+B                                             ║\n";
        cout << "║  4. Intersection:  A*B                                             ║\n";
        cout << "║  5. Difference:    A-B                                             ║\n";
        cout << "║  6. Help:          help                                            ║\n";
        cout << "║  7. Clear screen:  clear                                           ║\n";
        cout << "║  8. Exit:          exit                                            ║\n";
        cout << "╚════════════════════════════════════════════════════════════════════╝\033[0m\n";
    }
};

// Global set instances with better organization
map<char, SET> sets;

int main() {
    string input;
    SET::displayHelp();

    while (true) {
        SET::displayPrompt(sets['U'].setNames);
        
        if (!getline(cin, input)) break;
        input = SET().trim(input);
        if (input.empty()) continue;

        if (input == "exit") {
            cout << "\033[1;33mGoodbye!\033[0m\n";
            break;
        }
        
        if (input == "help") {
            SET::displayHelp();
            continue;
        }

        if (input == "clear") {
            cout << "\033[2J\033[H";
            continue;
        }

        // Process set operations
        if (input.length() >= 3 && isupper(input[0])) {
            char setName = input[0];
            char operation = input[1];

            switch (operation) {
                case '=':
                    if (sets[setName].insertElements(input)) {
                        sets['U'].setNames.insert(setName);
                    }
                    break;

                case '?':
                    if (sets.find(setName) != sets.end()) {
                        sets[setName].displaySet(string(1, setName));
                    } else {
                        cout << "\033[1;31mError: Set " << setName << " not defined\033[0m\n";
                    }
                    break;

                case '+':
                case '*':
                case '-':
                    if (input.length() == 3 && isupper(input[2])) {
                        char secondSet = input[2];
                        SET result;
                        
                        if (sets.find(setName) == sets.end() || sets.find(secondSet) == sets.end()) {
                            cout << "\033[1;31mError: One or both sets not defined\033[0m\n";
                            break;
                        }

                        switch (operation) {
                            case '+':
                                result = sets[setName].setUnion(sets[secondSet]);
                                result.displayOperation(input, result);
                                break;
                            case '*':
                                result = sets[setName].setIntersection(sets[secondSet]);
                                result.displayOperation(input, result);
                                break;
                            case '-':
                                result = sets[setName].setDifference(sets[secondSet]);
                                result.displayOperation(input, result);
                                break;
                        }
                    }
                    break;
            }
        }
    }
    return 0;
}