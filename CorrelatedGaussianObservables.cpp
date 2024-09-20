#include "CorrelatedGaussianObservables.h"

CorrelatedGaussianObservables::CorrelatedGaussianObservables(vector<dato> v_i,
        const TMatrixDSym& corr_i)
  {

    Int_t n = v_i.size(), i=0;
    Obs.ResizeTo(n);
    TVectorD Sig(n);

    for(vector<dato>::iterator it = v_i.begin(); it != v_i.end(); ++it) {
      Sig(i) = it->getSigma();
      Obs(i++) = it->getMean(); //equivalente a Obs(i) = it -> GetMean(); i++;
    }

    Cov.ResizeTo(n, n);
    for(int i = 0; i < n; i++)
      for(int j = i; j < n; j++)
	Cov(i,j) = Sig(i)*corr_i(i,j)*Sig(j);

    //norm = 1./sqrt(pow(2.*M_PI,n)*Cov.Determinant()); //Fattore di Normalizzazione inutile per il fattore di Bayes

    Cov.InvertFast();

  }

CorrelatedGaussianObservables::CorrelatedGaussianObservables(vector<dato> v_i,
        const TMatrixDSym& corr_1, const TMatrixDSym& corr_2)
  {
    Int_t n = v_i.size(), i=0;
    Obs.ResizeTo(n);
    TVectorD Sig1(n), Sig2(n);

    for(vector<dato>::iterator it = v_i.begin(); it != v_i.end(); ++it) {
      Sig1(i) = it->getSigma1();
      Sig2(i) = it->getSigma2();
      Obs(i++) = it->getMean();
    }

    Cov.ResizeTo(n, n);
    for(int i = 0; i < n; i++)
      for(int j = i; j < n; j++)
	Cov(i,j) = Sig1(i)*corr_1(i,j)*Sig1(j) + Sig2(i)*corr_2(i,j)*Sig2(j);

    //    norm = 1./sqrt(pow(2.*M_PI,n)*Cov.Determinant());

    Cov.InvertFast(); //Calcola l'inversa

  }

CorrelatedGaussianObservables::CorrelatedGaussianObservables(vector<dato> v_i,
        const TMatrixDSym& corr_1, const TMatrixDSym& corr_2, const TMatrixDSym& corr_3)
  {
    Int_t n = v_i.size(), i=0;
    Obs.ResizeTo(n);
    TVectorD Sig1(n), Sig2(n), Sig3(n);

    for(vector<dato>::iterator it = v_i.begin(); it != v_i.end(); ++it) {
      Sig1(i) = it->getSigma1();
      Sig2(i) = it->getSigma2();
      Sig3(i) = it->getSigma3();
      Obs(i++) = it->getMean();
    }

    Cov.ResizeTo(n, n);
    for(int i = 0; i < n; i++)
      for(int j = i; j < n; j++)
	      Cov(i,j) = Sig1(i)*corr_1(i,j)*Sig1(j) + Sig2(i)*corr_2(i,j)*Sig2(j) + Sig3(i)*corr_3(i,j)*Sig3(j);

    //    norm = 1./sqrt(pow(2.*M_PI,n)*Cov.Determinant());

    Cov.InvertFast(); //Calcola l'inversa

  }


double CorrelatedGaussianObservables::logweight(const TVectorD& v) {
    int n = Obs.GetNrows();
    double chisq = 0.;

    if (n != v.GetNrows()) {
      cout << "Errore in getWeight()" << endl;
      exit(EXIT_FAILURE);
    }
    for(int i = 0; i < n; i++)
      for(int j = 0; j < n; j++)
	chisq += (v(i)-Obs(i))*Cov(i,j)*(v(j)-Obs(j));

    return(-0.5*chisq);
  }
