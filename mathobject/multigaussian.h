// Class definition of MultiGaussian
// Version 1.0
#ifndef MULTIGAUSSIAN_H
#define MULTIGAUSSIAN_H

#include "miscmath.h"

#include <vector>
#include <iostream>

#include <vector>
#include <iostream>

/*! Class for generating multivariate gaussian random variables. This Class
generates a set of N gaussian random variables given a covariance matrix S.
The procedure works as follows: Find eigenvectors and eigenvalues of S, then
create the transformation matrix T for diagonalization. Generate N gaussian
random variables y with sigma^2 given by the eigenvalues of S. Then obtain
x=T*y as the required set of random variables.
A ususal call sequence would be
\code
Multigaussian mgauss(6);
mgauss.setBounds(lowBoundVector,highBoundVector);
mgauss.generateEigenvectors(covarianceMatrix);
mgauss.throwDice();
var=mgauss.getRandomValue(0);
\endcode
You have to call generateEigenvectors(vector<vector<double> >) only if you have a new covariance
matrix. throwDice(double *meanarray) generates a new set of random values centered at meanarray. Note that the
function getRandomValue(int i) is only for extracting the variables generated
by throwDice(double*), this function will not generate new variables.
*/

class MultiGaussian
{
  private:
    double Scale;
    bool Lock;

  protected:
    virtual void generateAndTransform();

    /*!The row number of the covariance matrix and number of
      random variables sought. */
    unsigned int _SIZE;

    /*!Upper and lower bounds are encoded here*/
    double *_lbounds;
    double *_hbounds;

    /*! The mean of the gaussian to sample from*/
    double *_center;

    /*! Encodes the eigenvectors/transformation matrix*/
    std::vector<std::vector<double> > _MasterMatrix;

    /*!The eigenvalues*/
    double *_eigenvalues;

    /*! The random values, generated anew if throwDice() is called*/
    std::vector<double> _generatedValues;

    Miscmath math;

  public:
    MultiGaussian(unsigned int Size);
    virtual ~MultiGaussian();
    unsigned int size() const { return _SIZE; }
    void setBounds(std::vector<double> lowBound, std::vector<double> highBound);
    void generateEigenvectors(std::vector<std::vector<double> >  covarianceMatrix, double scale);
    bool throwDice(const std::vector<double>& center);
    bool throwDice(const double* mean);
    double getRandomValue(int i) const;
    std::vector<double> getRandomValues() const;
    void printInfo(std::ostream& o);
    void printExtInfo(std::ostream& o);
    void lock() { Lock=true;}
    std::vector<std::vector<double> > copyMasterMatrix() { return _MasterMatrix; } //!< return copy of master matrix
};

#endif
