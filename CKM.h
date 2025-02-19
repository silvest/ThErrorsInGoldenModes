#ifndef CKM_H
#define CKM_H
#include "CorrelatedGaussianParameters.h"
#include "ModelParameter.h"
#include <vector>
#include <complex>
#include <TMatrixD.h>
#include <TMatrixDSym.h>
#include <TRandom3.h>
#include <TVectorD.h>
#include <TDecompChol.h>

class CKMParameters {
public:
    // Constructor
  CKMParameters(); //constructor to be used if I want to work with CKM elements
  CKMParameters(bool initialize); //constructor to be used if i want to construct the matrix from Wolfenstein parameters with correlation matrix

    // Sample CKM parameters from the multivariate Gaussian
    void sampleCKMParameters(const std::vector<double>& diagonalizedPars);
     // Compute the CKM matrix based on sampled Wolfenstein parameters
    void computeCKMwithWolfenstein();
    void computeCKM(double Vus_v, double Vcb_v, double Vub_v, double gamma_v, bool useVud);
    void computeCKMfromAngles();

    // Helper function to compute the CKM matrix from angles (Wolfenstein parameters)
    void computeCKMfromAngles(double lambda, double A, double rho, double eta);

    // Print sampled CKM parameters
    void printParameters() const;
    // Access CKM matrix elements (direct values)
    std::complex<double> getVud() const { return V[0][0]; }
    std::complex<double> getVus() const { return V[0][1]; }
    std::complex<double> getVub() const { return V[0][2]; }

    std::complex<double> getVcd() const { return V[1][0]; }
    std::complex<double> getVcs() const { return V[1][1]; }
    std::complex<double> getVcb() const { return V[1][2]; }

    std::complex<double> getVtd() const { return V[2][0]; }
    std::complex<double> getVts() const { return V[2][1]; }
    std::complex<double> getVtb() const { return V[2][2]; }

    // Access conjugates of CKM matrix elements
    std::complex<double> getVdu() const { return std::conj(V[0][0]); }
    std::complex<double> getVsu() const { return std::conj(V[0][1]); }
    std::complex<double> getVbu() const { return std::conj(V[0][2]); }

    std::complex<double> getVdc() const { return std::conj(V[1][0]); }
    std::complex<double> getVsc() const { return std::conj(V[1][1]); }
    std::complex<double> getVbc() const { return std::conj(V[1][2]); }

    std::complex<double> getVdt() const { return std::conj(V[2][0]); }
    std::complex<double> getVst() const { return std::conj(V[2][1]); }
    std::complex<double> getVbt() const { return std::conj(V[2][2]); }

       // Add correlations and uncertainties as member variables
    TVectorD uncertainties;  // Store uncertainties for the CKM parameters
    TMatrixD correlations;   // Store correlation matrix for CKM parameters



   // Implement q/p for B0 (B_d) mixing
    std::complex<double> get_q_p_Bd() const;

    // Implement q/p for B_s mixing
    std::complex<double> get_q_p_Bs() const;
    std::complex<double> get_q_p_KS() const;

    std::vector<double> TVectorDToStdVector(const TVectorD& tvec) {
    return std::vector<double>(tvec.GetMatrixArray(), tvec.GetMatrixArray() + tvec.GetNrows());
}


private:
  // Create an instance of the CorrelatedGaussianParameters class
    CorrelatedGaussianParameters CGP;
    TVectorD meanValues;   // Mean values of the Wolfenstein parameters
    TVectorD sampledValues;  // Sampled Wolfenstein parameters
    TMatrixDSym covariance;  // Covariance matrix
    TMatrixD L;  // Cholesky decomposition result (lower triangular matrix)

    double Rho, Eta, Lambda, A; ///< The Wolfenstein parameters.
    double s12, s13, s23; ///< The sine of the three mixing angles
    double c12, c23, c13; ///< The cosine of the three mixing angles
    double delta; ///< The CP violating phase in the CKM matrix.
    std::vector<std::vector<std::complex<double>>> V;  // CKM matrix



};

#endif  // CKM_H
