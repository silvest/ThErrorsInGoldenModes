/*
 * File:   dato.h
 * Author: silvest
 *
 * Created on October 31, 2012, 3:40 PM
 */

#ifndef DATO_H
#define	DATO_H

#include <TRandom3.h>
#include <math.h>
#include <iostream>

class dato {

public:
    dato(double ave, double sig1, double sig2 = 0., double sig3 = 0.) {
        sigma1 = sig1;
        sigma2 = sig2;
        sigma3 = sig3;
        mean = ave;
        sigma = sqrt(sigma1 * sigma1 + sigma2 * sigma2 + sigma3 * sigma3);
    };

    virtual ~dato() {
    };

    double getMean() {
        return mean;
    };

    double getSigma() {
        return sigma;
    };

    double weight(double x) {
        return exp(-0.5 * (x - mean)*(x - mean) / sigma / sigma);
    };


    //Peso Gaussiano del dato
    double logweight(double x) {
        return (-0.5 * (x - mean)*(x - mean) / sigma / sigma);
    };

    double getRandom() {
        return gRandom->Gaus(mean, sigma);
    };

    double getFlat() {
        return gRandom->Uniform(mean - sigma, mean + sigma);
    };

    double getSigma1() const
    {
        return sigma1;
    }

    double getSigma2() const
    {
        return sigma2;
    }

    double getSigma3() const
    {
        return sigma3;
    }

private:
    double mean, sigma, sigma1, sigma2, sigma3;

};

#endif	/* DATO_H */
