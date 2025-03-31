/* 
 * Copyright (C) 2013 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include <stdexcept>
#include <vector>
#include <fstream>
#include <sstream>
#include <math.h>
#include <boost/lexical_cast.hpp>
#include "CorrelatedGaussianParameters.h"

CorrelatedGaussianParameters::CorrelatedGaussianParameters(std::string name_i) {
    name = name_i;
    v = NULL;
    e = NULL;
    IsEOF = false;
}

CorrelatedGaussianParameters::CorrelatedGaussianParameters() {
    v = NULL;
    e = NULL;
    IsEOF = false;
}

CorrelatedGaussianParameters::CorrelatedGaussianParameters(const CorrelatedGaussianParameters& orig) : Cov(orig.Cov)
{
    Pars = orig.Pars;
    name = orig.name;
    v = new TMatrixD(*orig.v);  // Replace gslpp::matrix with TMatrixD
    e = new TVectorD(*orig.e);  // Replace gslpp::vector with TVectorD
    DiagPars = orig.DiagPars;
    IsEOF = orig.IsEOF;
}

CorrelatedGaussianParameters::~CorrelatedGaussianParameters() {
    if (v != NULL)
        delete v;  // Use TMatrixD instead of gslpp::matrix
    if (e != NULL)
        delete e;  // Use TVectorD instead of gslpp::vector
}


void CorrelatedGaussianParameters::AddPar(ModelParameter& Par_i) {
    Pars.push_back(Par_i);
}

void CorrelatedGaussianParameters::DiagonalizePars(TMatrixDSym Corr) {
    unsigned int size = Pars.size();
    if (Corr.GetNrows() != size || Corr.GetNcols() != size)
        throw std::runtime_error("The size of the correlated parameters in " + name + " does not match the size of the correlation matrix!");
    
    Cov.ResizeTo(size, size);
    for (unsigned int i = 0; i < size; i++)
        for (unsigned int j = 0; j < size; j++)
            Cov(i, j) = Pars.at(i).geterrg() * Corr(i, j) * Pars.at(j).geterrg();

    TMatrixDSymEigen SE(Cov);

    TMatrixD vv(SE.GetEigenVectors());
    TVectorD ee(SE.GetEigenValues());

    v = new TMatrixD(size, size);  // Allocate new TMatrixD for v
    e = new TVectorD(size);        // Allocate new TVectorD for e
    unsigned int EVbad = 0;

    for (unsigned int i = 0; i < size; i++) {
        (*e)(i) = ee(i);  // Copy eigenvalues
        for (unsigned int j = 0; j < size; j++) {
            (*v)(i, j) = vv(i, j);  // Copy eigenvectors
        }
        if (ee(i) <= 0.) {
            EVbad++;
        }
    }

    if (EVbad > 0) {
        std::cout << "WARNING: Covariance matrix of the correlated parameters in "<< name <<" is not a positive definite matrix!" << std::endl;
        std::cout << "("<< EVbad <<" non positive eigenvalue(s).)" << std::endl;
        sleep(2);
    }

    TVectorD ave_in(size);  // Create a TVectorD for average values

    int ind = 0;
    for (std::vector<ModelParameter>::iterator it = Pars.begin(); it != Pars.end(); it++) {
        ave_in(ind) = it->getave();  // Copy average values to ave_in
        ind++;
    }

    TVectorD ave = v->T() * ave_in;  // Perform matrix-vector multiplication

    for (unsigned int i = 0; i < size; i++) {
        std::stringstream ss;
        ss << (i + 1);
        std::string namei = name + ss.str();
        ModelParameter current(namei, ave(i), sqrt((*e)(i)), 0.);  // Use ave and eigenvalue
        current.setCgp_name(name);
        DiagPars.push_back(current);  // Push the parameter to DiagPars
    }
}


std::vector<double> CorrelatedGaussianParameters::getOrigParsValue(const std::vector<double>& DiagPars_i) const {
    if (DiagPars_i.size() != DiagPars.size()) {
        std::stringstream out;
        out << DiagPars_i.size();
        throw std::runtime_error("CorrelatedGaussianParameters::getOrigParsValue(DiagPars_i): DiagPars_i.size() = " + out.str() + " does not match the size of DiagPars");
    }

    TVectorD pars_in(DiagPars_i.size());  // Use TVectorD instead of gslpp::vector
    pars_in.Zero();  // Initialize to zero

    for (unsigned int i = 0; i < DiagPars_i.size(); i++) {
        pars_in(i) = DiagPars_i[i];  // Copy the values
    }

    TVectorD val = (*v) * pars_in;  // Matrix-vector multiplication with TMatrixD

    std::vector<double> res;
    for (unsigned int i = 0; i < DiagPars_i.size(); i++) {
        res.push_back(val(i));  // Store the result
    }

    return res;
}



int CorrelatedGaussianParameters::ParseCGPFromData(
    std::vector<ModelParameter>& ModPars,
    const std::vector<double>& meanValues,          // Means for each parameter
    const std::vector<double>& uncertainties,       // Uncertainties for each parameter
    const TMatrixD& correlations,                   // Correlation matrix
    int rank) {
    
    name = "CorrelatedGaussianParameters";  // Use a default name or pass it as an argument
    int size = meanValues.size();            // Get the size from the mean values vector
    int nlines = 0;

    // Check if the sizes of inputs match
    if (uncertainties.size() != size || correlations.GetNrows() != size || correlations.GetNcols() != size) {
        throw std::runtime_error("Size mismatch between mean values, uncertainties, or correlation matrix.");
    }

    std::vector<bool> lines(size, true);  // Assume all parameters are not fixed

    // Add parameters to the Pars vector based on the input meanValues and uncertainties
    for (int i = 0; i < size; ++i) {
        ModelParameter tmpMP;
        tmpMP.setValue(meanValues[i]);        // Set the mean value
        tmpMP.setError(uncertainties[i]);     // Set the uncertainty
        tmpMP.setName("Param" + std::to_string(i));  // Generate parameter name
        tmpMP.setCgp_name(name);              // Associate the parameter with the set name
        AddPar(tmpMP);                        // Add the parameter to Pars
        ModPars.push_back(tmpMP);             // Also add it to ModPars
        nlines++;
    }

    // Proceed to diagonalize the correlation matrix if there is more than one parameter
    if (nlines > 1) {
        TMatrixDSym mySCorr(nlines);
        for (int i = 0; i < nlines; i++) {
            for (int j = 0; j < nlines; j++) {
                mySCorr(i, j) = correlations(i, j);
            }
        }

        DiagonalizePars(mySCorr);
        ModPars.insert(ModPars.end(), getDiagPars().begin(), getDiagPars().end());
    } else {
        if (rank == 0) {
            std::cout << "\nWARNING: Correlated Gaussian Parameters " << name.c_str() << " defined with less than two correlated parameters. The set is being marked as normal Parameters." << std::endl;
        }
    }

    return nlines;
}
