#include "PDGAverage.h"
#include <algorithm>

//default constructor
PDGAverage::PDGAverage() {};

// constructor
PDGAverage::PDGAverage(std::string name, const std::vector<dato> & data) : fData(data) {
    fName = name;   
    isAngle = false;
    for (std::vector<dato>::iterator it = fData.begin(); it != fData.end(); ++it) {
        if (it->getIsAngle() && !isAngle) {
            isAngle = true;
        }
        else if (!it->getIsAngle() && isAngle) {
            std::cerr << "Error: Inconsistent angle flag in data for " << fName << std::endl;
            exit(1);
        }
    }
    CalculateAverage();
}

// destructor
PDGAverage::~PDGAverage() {}

// methods
void PDGAverage::CalculateAverage() {
    double sum = 0.;
    double sum2 = 0.;
    for (std::vector<dato>::iterator it = fData.begin(); it != fData.end(); ++it) {
        sum += (it->getIsAngle() ? remainder(it->getMean(), 2.*M_PI) : it->getMean()) / it->getSigma() / it->getSigma();
        sum2 += 1. / it->getSigma() / it->getSigma();
    }
    fAverage = sum / sum2;
    fUncertainty = 1. / sqrt(sum2);
    // compute the PDG scale factor
    double chi2 = 0.;
    for (std::vector<dato>::iterator it = fData.begin(); it != fData.end(); ++it) {
        double num = (it->getIsAngle() ? remainder(it->getMean() - fAverage, 2.*M_PI) : it->getMean() - fAverage);
        chi2 += num * num / it->getSigma() / it->getSigma();
    }
    fScaleFactor = std::max(1.,sqrt(chi2 / (fData.size() - 1.)));
    // rescale the uncertainty
    fUncertainty *= fScaleFactor;
    std::cout << "PDG average for " << fName << ": " << fAverage << " +- " << fUncertainty << " with scale factor " << fScaleFactor << std::endl;
}
