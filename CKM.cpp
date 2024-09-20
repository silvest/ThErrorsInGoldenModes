#include "CKM.h"
#include <iostream>
#include <cmath>
#include <TMatrixDEigen.h>
#include <stdexcept>
#include "CorrelatedGaussianParameters.h"



// Constructor
CKMParameters::CKMParameters() 
    : meanValues(4), covariance(4), sampledValues(4), uncertainties(4), randGen(0), V(3, std::vector<std::complex<double>>(3)) {

    // Set the central values for the Wolfenstein parameters
    meanValues[0] = 0.2251;  // lambda 
    meanValues[1] = 0.828;    // A
    meanValues[2] = 0.160;    // rho_bar
    meanValues[3] = 0.346;    // eta_bar

    // Uncertainties on the Wolfenstein parameters
    uncertainties[0] = 0.0005;
    uncertainties[1] = 0.015;
    uncertainties[2] = 0.023;
    uncertainties[3] = 0.015;

    // Correlation coefficients between the Wolfenstein parameters
    correlations.ResizeTo(4, 4);  // Ensure the matrix has the correct size
    correlations(0, 0) = 1.00; correlations(0, 1) = -0.581; correlations(0, 2) = 0.223; correlations(0, 3) = -0.127;
    correlations(1, 0) = -0.581; correlations(1, 1) = 1.000; correlations(1, 2) = -0.105; correlations(1, 3) = -0.280;
    correlations(2, 0) = 0.223; correlations(2, 1) = -0.105; correlations(2, 2) = 1.000; correlations(2, 3) = 0.098;
    correlations(3, 0) = -0.127; correlations(3, 1) = -0.280; correlations(3, 2) = 0.098; correlations(3, 3) = 1.000;
}

void CKMParameters::sampleCKMParameters() {
    // Create an instance of the CorrelatedGaussianParameters class
    CorrelatedGaussianParameters CGP;
    std::vector<ModelParameter> ModPars;

    // Convert TVectorD (meanValues, uncertainties) to std::vector<double>
    std::vector<double> meanVec = TVectorDToStdVector(meanValues);
    std::vector<double> uncertaintyVec = TVectorDToStdVector(uncertainties);

    // Sample the CKM parameters (lambda, A, rho, eta) using the correlation matrix
    CGP.ParseCGPFromData(ModPars, meanVec, uncertaintyVec, correlations, 0);  // Rank = 0 for simplicity

    // Extract the sampled CKM parameters (from ModPars) and store them in sampledValues
    // Limit the loop to the size of sampledValues, which is 4
    for (size_t i = 0; i < sampledValues.GetNrows() && i < ModPars.size(); ++i) {
        sampledValues[i] = ModPars[i].getave(); 
    }


    // Compute the CKM matrix using the sampled Wolfenstein parameters
    computeCKMwithWolfenstein();
}


// Compute the CKM matrix with the Wolfenstein parameters
void CKMParameters::computeCKMwithWolfenstein() {
    double lambda = sampledValues[0];
    double A = sampledValues[1];
    double rho = sampledValues[2];
    double eta = sampledValues[3];

    // Compute CKM matrix based on sampled Wolfenstein parameters
    computeCKMfromAngles(lambda, A, rho, eta);
}

// Helper function to compute the CKM matrix from angles (Wolfenstein parameters)
void CKMParameters::computeCKMfromAngles(double lambda, double A, double rho, double eta) {
    std::complex<double> num(rho, eta);  // num = rho + i*eta
    num = num * sqrt(1.0 - pow(A, 2.0) * pow(lambda, 4.0));  // Adjust numerator
    
    std::complex<double> den = sqrt(1.0 - pow(lambda, 2.0)) * 
                               std::complex<double>(1.0 - pow(A, 2.0) * pow(lambda, 4.0) * rho, 
                                                    -pow(A, 2.0) * pow(lambda, 4.0) * eta);  // Denominator
    
    std::complex<double> ratio = num / den;  // ratio = num / den

    double rho_nb = ratio.real();  // Real part
    double eta_nb = ratio.imag();  // Imaginary part

    // Sine of the mixing angles
    double s12 = lambda;
    double s23 = A * pow(lambda, 2.0);
    double s13 = std::abs(std::complex<double>(A * pow(lambda, 3.0) * rho_nb, -A * pow(lambda, 3.0) * eta_nb));
    
    // CP-violating phase (delta)
    double delta = -std::arg(std::complex<double>(A * pow(lambda, 3.0) * rho_nb, -A * pow(lambda, 3.0) * eta_nb));

    // Cosines of the mixing angles
    double c12 = sqrt(1.0 - s12 * s12);
    double c13 = sqrt(1.0 - s13 * s13);
    double c23 = sqrt(1.0 - s23 * s23);

    // Compute the CKM matrix based on these angles
    V[0][0] = std::complex<double>(c12 * c13, 0.0);
    V[0][1] = std::complex<double>(s12 * c13, 0.0);
    V[0][2] = std::complex<double>(s13 * cos(delta), -s13 * sin(delta));

    V[1][0] = std::complex<double>(-s12 * c23 - c12 * s23 * s13 * cos(delta), -c12 * s23 * s13 * sin(delta));
    V[1][1] = std::complex<double>(c12 * c23 - s12 * s23 * s13 * cos(delta), -s12 * s23 * s13 * sin(delta));
    V[1][2] = std::complex<double>(s23 * c13, 0.0);

    V[2][0] = std::complex<double>(s12 * s23 - c12 * c23 * s13 * cos(delta), -c12 * c23 * s13 * sin(delta));
    V[2][1] = std::complex<double>(-c12 * s23 - s12 * c23 * s13 * cos(delta), -s12 * c23 * s13 * sin(delta));
    V[2][2] = std::complex<double>(c23 * c13, 0.0);
}

// Implement q/p for B0 (B_d)
std::complex<double> CKMParameters::get_q_p_Bd() const {
    std::complex<double> Vtd = getVtd();
    std::complex<double> Vtb = getVtb();
    return -std::conj(Vtd) * Vtb / (Vtd * std::conj(Vtb));
}

// Implement q/p for B_s
std::complex<double> CKMParameters::get_q_p_Bs() const {
    std::complex<double> Vts = getVts();
    std::complex<double> Vtb = getVtb();
    return -std::conj(Vts) * Vtb / (Vts * std::conj(Vtb));
}

// Implement q/p for K0 (K_S mixing)
std::complex<double> CKMParameters::get_q_p_KS() const {
    // Assuming you have CKM elements for K0 (Vcd and Vcs)
    std::complex<double> Vcd = V[1][0];  // Vcd (K0)
    std::complex<double> Vcs = V[1][1];  // Vcs (K0)
    return -std::conj(Vcd) * Vcs / (Vcd * std::conj(Vcs));
}

// Print the sampled CKM parameters
void CKMParameters::printParameters() const {
    std::cout << "Sampled CKM Parameters:" << std::endl;
    std::cout << "Lambda: " << sampledValues[0] << std::endl;
    std::cout << "A: " << sampledValues[1] << std::endl;
    std::cout << "Rho: " << sampledValues[2] << std::endl;
    std::cout << "Eta: " << sampledValues[3] << std::endl;
}
