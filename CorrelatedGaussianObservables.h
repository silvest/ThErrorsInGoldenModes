/*
 * File:   CorrelatedGaussianObservables.h
 * Author: silvest
 *
 * Created on January 21, 2014, 12:14 PM
 */

#ifndef CORRELATEDGAUSSIANOBSERVABLES_H
#define	CORRELATEDGAUSSIANOBSERVABLES_H

#include <TVectorD.h>
#include <TMatrixD.h>
#include <vector>
#include <cstdlib>
#include "dato.h"

using namespace std;

class CorrelatedGaussianObservables {
public:
  virtual ~CorrelatedGaussianObservables() {};
  CorrelatedGaussianObservables(vector<dato> v_i, const TMatrixDSym& corr_1, const TMatrixDSym& corr_2, const TMatrixDSym& corr_3) ;
  CorrelatedGaussianObservables(vector<dato> v_i, const TMatrixDSym& corr_1, const TMatrixDSym& corr_2) ;
  CorrelatedGaussianObservables(vector<dato> v_i, const TMatrixDSym& corr_i) ;
  CorrelatedGaussianObservables(const CorrelatedGaussianObservables& orig)
  :Obs(orig.getObs()),Cov(orig.getCov()) {};

  const TVectorD& getObs() const { return Obs; }

  const TMatrixDSym& getCov() const { return Cov; }

  double logweight(const TVectorD& v) ;

private:
  TVectorD Obs;
  TMatrixDSym Cov;
};

#endif	/* CORRELATEDGAUSSIANOBSERVABLES_H */
