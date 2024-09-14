#include "PDGAverage.h"
#include <algorithm>

//default constructor
PDGAverage::PDGAverage() {};

// constructor
PDGAverage::PDGAverage(std::string name, const std::vector<dato> & data) : fData(data) {
    fName = name;   
        CalculateAverage();
    }

    // destructor
    ~PDGAverage() {}

    // methods
    void CalculateAverage() {
        double sum = 0.;
        double sum2 = 0.;
        for (std::vector<dato>::iterator it = fData.begin(); it != fData.end(); ++it) {
            sum += it->getMean() / it->getSigma() / it->getSigma();
            sum2 += 1. / it->getSigma() / it->getSigma();
        }
        fAverage = sum / sum2;
        fUncertainty = 1. / sqrt(sum2);
        // compute the PDG scale factor
        double chi2 = 0.;
        for (std::vector<dato>::iterator it = fData.begin(); it != fData.end(); ++it) {
            chi2 += (it->getMean() - fAverage) * (it->getMean() - fAverage) / it->getSigma() / it->getSigma();
        }
        fScaleFactor = std::max(1.,sqrt(chi2 / (fData.size() - 1.)));
        // rescale the uncertainty
        fUncertainty *= fScaleFactor;
        std::cout << "PDG average for" << fName << ": " << fAverage << " +- " << fUncertainty << " with scale factor " << fScaleFactor << std::endl;
    }
