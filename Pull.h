#ifndef PULL_H
#define PULL_H

#include <string>
#include <cmath>

// Simple structure to calculate and store pull values
struct Pull {
    std::string name;
    double predicted;
    double observed;
    double uncertainty;
    double pullValue;
    
    Pull(const std::string& n, double pred, double obs, double unc) 
        : name(n), predicted(pred), observed(obs), uncertainty(unc) {
        pullValue = (predicted - observed) / uncertainty;
    }
    
    Pull() : predicted(0), observed(0), uncertainty(1), pullValue(0) {}
};

#endif // PULL_H
