#ifndef PDGAVERAGE_H
#define PDGAVERAGE_H

#include "dato.h"

class PDGAverage {
public:
    PDGAverage(); // default constructor
    PDGAverage(std::string name, const std::vector<dato> & data, const bool & isAngle = false); // constructor
    ~PDGAverage(); // destructor

    // methods
    void CalculateAverage();

    // setter for fData
    void setData(const std::vector<dato>& newData) {
        fData = newData;
    }

    // getter for fAverage
    double getAverage() const {
        return fAverage;
    }

    // getter for fUncertainty
    double getUncertainty() const {
        return fUncertainty;
    }

    // getter for fScaleFactor
    double getScaleFactor() const {
        return fScaleFactor;
    }

    // getter for fName
    std::string getName() const {
        return fName;
    }

    // setter for fName
    void setName(const std::string& newName) {
        fName = newName;
    }

    // setter for fIsAngle
    void setIsAngle(bool isAngle) {
        fIsAngle = isAngle;
    }

private:
    // member variables
    std::string fName;
    std::vector<dato> fData;
    double fAverage;
    double fUncertainty;
    double fScaleFactor;
    bool fIsAngle;

    // Helper: signed angular difference wrapped into (-pi, pi]
    double AngleDiff(double pred, double meas) const
    {
        return std::remainder(pred - meas, 2.0 * M_PI);
    }

};

#endif // PDGAVERAGE_H
