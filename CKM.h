#ifndef CKM_H
#define CKM_H

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
    CKMParameters();

    // Sample CKM parameters from the multivariate Gaussian
    void sampleCKMParameters();

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

    // Implement q/p for K0 (K_S mixing)
    std::complex<double> get_q_p_KS() const;

    std::vector<double> TVectorDToStdVector(const TVectorD& tvec) {
    return std::vector<double>(tvec.GetMatrixArray(), tvec.GetMatrixArray() + tvec.GetNrows());
}


private:
    TVectorD meanValues;   // Mean values of the Wolfenstein parameters
    TVectorD sampledValues;  // Sampled Wolfenstein parameters
    TMatrixDSym covariance;  // Covariance matrix
    TMatrixD L;  // Cholesky decomposition result (lower triangular matrix)
    TRandom3 randGen;  // Random number generator
    std::vector<std::vector<std::complex<double>>> V;  // CKM matrix

    // Compute the CKM matrix based on sampled Wolfenstein parameters
    void computeCKMwithWolfenstein();

    // Helper function to compute the CKM matrix from angles (Wolfenstein parameters)
    void computeCKMfromAngles(double lambda, double A, double rho, double eta);
};

#endif  // CKM_H




   
