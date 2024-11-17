#include "CKM.h"
#include <iostream>
#include <cmath>
#include <TMatrixDEigen.h>
#include <stdexcept>
#include "CorrelatedGaussianParameters.h"



// Constructor
CKMParameters::CKMParameters(bool initialize) 
    : meanValues(4), covariance(4), sampledValues(4), uncertainties(4), randGen(0), V(3, std::vector<std::complex<double>>(3)) {

    // Set the central values for the Wolfenstein parameters
    meanValues[0] = 0.22504;  // lambda 
    meanValues[1] = 0.816;    // A
    meanValues[2] = 0.163;    // rho_bar
    meanValues[3] = 0.356;    // eta_bar

    // Uncertainties on the Wolfenstein parameters
    uncertainties[0] = 0.00077;
    uncertainties[1] = 0.016;
    uncertainties[2] = 0.024;
    uncertainties[3] = 0.028;

    // Correlation coefficients between the Wolfenstein parameters
    correlations.ResizeTo(4, 4);  // Ensure the matrix has the correct size
    correlations(0, 0) = 1.00; correlations(0, 1) = -0.581; correlations(0, 2) = 0.223; correlations(0, 3) = -0.127;
    correlations(1, 0) = -0.581; correlations(1, 1) = 1.000; correlations(1, 2) = -0.105; correlations(1, 3) = -0.280;
    correlations(2, 0) = 0.223; correlations(2, 1) = -0.105; correlations(2, 2) = 1.000; correlations(2, 3) = 0.098;
    correlations(3, 0) = -0.127; correlations(3, 1) = -0.280; correlations(3, 2) = 0.098; correlations(3, 3) = 1.000;

     // Create an actual std::vector<ModelParameter> 
    std::vector<ModelParameter> modPars;


    // Initialize the CorrelatedGaussianParameters with the data
    CGP.ParseCGPFromData(modPars, TVectorDToStdVector(meanValues), TVectorDToStdVector(uncertainties), correlations, 0);
    }

//constructor
CKMParameters::CKMParameters() : V(3, std::vector<std::complex<double>>(3, std::complex<double>(0.0, 0.0))) {}


// Sample CKM parameters using the diagonalized parameters provided by BAT
void CKMParameters::sampleCKMParameters(const std::vector<double>& diagonalizedPars) {
    // Convert diagonalized parameters to original CKM parameters
    std::vector<double> originalPars = CGP.getOrigParsValue(diagonalizedPars);
    
    // Store the sampled CKM parameters
    for (int i = 0; i < 4; ++i) {
        sampledValues[i] = originalPars[i];
    }

    // Compute the CKM matrix from the sampled Wolfenstein parameters
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

void CKMParameters::computeCKM(double Vus_v, double Vcb_v, double Vub_v, double gamma_v, bool useVud)
{
   // Calculate sines and cosines for the angles
    s13 = Vub_v;
    c13 = std::sqrt(1.0 - s13 * s13);

    if (useVud) {
        c12 = Vus_v / c13;
        s12 = std::sqrt(1.0 - c12 * c12);
    } else {
        s12 = Vus_v / c13;
        c12 = std::sqrt(1.0 - s12 * s12);
    }

    s23 = Vcb_v / c13;
    c23 = std::sqrt(1.0 - s23 * s23);

    // Calculate the phase delta using gamma_v
    double a = c12 * s13 * s23 / (s12 * c23);
    if (std::fabs(gamma_v) < 1.e-10) {
        delta = 0.0;
    } else {
        delta = 2.0 * atan((1.0 + sqrt(1.0 - (a * a - 1.0) * pow(tan(gamma_v), 2.0)) * (cos(gamma_v) < 0.0 ? 1.0 : -1.0)) / ((a - 1.0) * tan(gamma_v)));
    }

    // Compute the CKM elements based on these angles and delta
    computeCKMfromAngles();

    // Wolfenstein to all orders (updating parameters Lambda, A, Rho, Eta)
    Lambda = s12;
    A = s23 / (Lambda * Lambda);

    // Calculate Rb using std::complex<double>
    std::complex<double> V_ud = std::conj(V[0][0]);
    std::complex<double> V_ub = std::conj(V[0][2]);
    std::complex<double> V_cd = std::conj(V[1][0]);
    std::complex<double> V_cb = std::conj(V[1][2]);

    std::complex<double> Rb = (V_ud * V_ub) / (V_cd * V_cb);
    Rho = -Rb.real();
    Eta = -Rb.imag();
}


// Compute CKM elements from sines, cosines, and phase angles
void CKMParameters::computeCKMfromAngles() {
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
    return std::conj(Vtb) * Vtd / (Vtb * std::conj(Vtd));
}

// Implement q/p for B_s
std::complex<double> CKMParameters::get_q_p_Bs() const {
    std::complex<double> Vts = getVts();
    std::complex<double> Vtb = getVtb();
    return std::conj(Vtb) * Vts / (Vtb * std::conj(Vts));
}

// Print the sampled CKM parameters
void CKMParameters::printParameters() const {
    std::cout << "Sampled CKM Parameters:" << std::endl;
    std::cout << "Lambda: " << sampledValues[0] << std::endl;
    std::cout << "A: " << sampledValues[1] << std::endl;
    std::cout << "Rho: " << sampledValues[2] << std::endl;
    std::cout << "Eta: " << sampledValues[3] << std::endl;
}

