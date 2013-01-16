#include "multigaussian.h"

#include "global.h"

#include "gsl/gsl_math.h"
#include "gsl/gsl_eigen.h"

#include <iomanip>
#include <vector>
#include <iostream>
#include <iterator>

/*!Constructor. Size<=0 will give an error.
\param Size this is the rank of the (symmetric) covariance matrix and, correspondingly the
    number of random variables sought.
*/
MultiGaussian::MultiGaussian(unsigned int Size) : Scale(1.0), Lock(false) {
    if (Size==0) {
      std::cout << "Error in MultiGaussian - Size is ==0";
    }                       // initialize the various arrays and variables
    _SIZE=Size;
    _MasterMatrix.resize(_SIZE);
    for(unsigned int i=0; i< _SIZE; i++){
      _MasterMatrix[i].resize(_SIZE);
    }
    _eigenvalues=new double[_SIZE];
    _generatedValues.resize(_SIZE);
    _lbounds=new double[_SIZE];
    _hbounds=new double[_SIZE];
    _center=new double[_SIZE];
    std::fill(_eigenvalues, _eigenvalues + _SIZE, 0.);
    std::fill(_generatedValues.begin(), _generatedValues.end(), 0.);
    std::fill(_lbounds, _lbounds + _SIZE, 0.);
    std::fill(_hbounds, _hbounds + _SIZE, 0.);
    std::fill(_center, _center + _SIZE, 0.);
}

MultiGaussian::~MultiGaussian() {
  delete [] _lbounds;
  delete [] _hbounds;
  delete [] _eigenvalues;
  delete [] _center;
}

/*! Set bounds for the random variables. This needs to be called before any other
function, or unpredictive behaviour will result. In most applications, you will
probably set this only once.
\param lowBound a vector with lower Bounds for the random variables
\param highBound a vector with upper Bounds for the random variables
*/
void MultiGaussian::setBounds(std::vector<double> lowBound, std::vector<double> highBound) {
  if(lowBound.size()!=_SIZE || highBound.size()!=_SIZE) {
    std::cout << "Error in MultiGaussian::setBounds:  highBound, lowBound are not of correct size";
  }
  for(unsigned int i=0; i <_SIZE; i++ ) {
    _lbounds[i]=lowBound[i];
    _hbounds[i]=highBound[i];
  }
}

/*! Generate a set of eigenvectors and eigenvalues from the covariance matrix.
Thus, if you whish to sample from a new distribution with covariance matrix S,
you need to call this function again. This function uses the eigenvalue and
eigenvector routines of the GNU Scientific Library, so if you need this for
large matrices it is probably best to change this code and use LAPACK or
something similar.
\param covarianceMatrix the covariance matrix of the distribution you want to draw from
*/
void MultiGaussian::generateEigenvectors(std::vector<std::vector<double> > covarianceMatrix,double scale) {
  if (Lock) throw Bad_Error("MultiGaussian::generateEigenvectors() MultiGaussian is locked!");
  Scale = scale;
  double *tempArray = new double[_SIZE*_SIZE];               // reformat [][] into [] for creating matrix
    int k=0;
    for(unsigned int j=0; j< _SIZE; j++)    {            // copy into tempArray
        for(unsigned int i=0; i< _SIZE; i++) {
             tempArray[k]=covarianceMatrix[i][j];
             k++;
         }
    }
    //!!!!!!!!!!!! part from GNU Scientific Library starts
    gsl_eigen_symmv_workspace * w = gsl_eigen_symmv_alloc (_SIZE);
    gsl_matrix_view m = gsl_matrix_view_array (tempArray, _SIZE, _SIZE);
    // ordering of rows and columns is irrelevant, as this is a symmetric matrix.
    gsl_vector *eval = gsl_vector_alloc (_SIZE);
    gsl_matrix *evec = gsl_matrix_alloc (_SIZE, _SIZE);
    gsl_eigen_symmv (&m.matrix, eval, evec, w); //obtain eigenvalues and eigenvectorx
    gsl_eigen_symmv_free (w);
    gsl_eigen_symmv_sort (eval, evec, GSL_EIGEN_SORT_ABS_ASC);  //sort in ascending order
                                                            // (not strictly neccessary)
    for(unsigned int i=0; i<_SIZE; i++) {                // copy into the class member variables.
         _eigenvalues[i] = gsl_vector_get (eval, i);
	 //  gsl_vector_view evec_i = gsl_matrix_column (evec, i);
         for (unsigned int j=0;j <_SIZE; j++){
         _MasterMatrix[i][j]= gsl_matrix_get(evec,i,j);
         }
    }
    
    for(unsigned int i=0; i<_SIZE; i++) {           
      _eigenvalues[i] *= scale;
    }

    gsl_vector_free(eval);   // free memory
    gsl_matrix_free(evec);
    delete[] tempArray;
    //!!!!!!!!!!!!!Nothing more from the GNU Scientific Library
}

/*! Generate a new set of variables, will be written to the private variable _generatedValues.
Call this function if you want to draw new variables from the distribution.
 Returns true if the step is within all boundaries, false, if it fell outside the boundaries
*/
bool MultiGaussian::throwDice(const std::vector<double>& mean) {
  for(unsigned int i=0;i < _SIZE; i++) _center[i]=mean[i];

  generateAndTransform();

  //static int falseCounter=0;
  for(unsigned int i=0; i<_SIZE;i++) {
    bool exceed_low = ((_generatedValues[i]+_center[i]) < _lbounds[i]);
    bool exceed_high = ((_generatedValues[i]+_center[i]) > _hbounds[i]);

    if (exceed_low || exceed_high) {
      //std::cout << (_generatedValues[i]+_center[i]) <<  " is out of bounds for param no : "
      //        << i << " lbound: " << _lbounds[i] << " and hbound: " << _hbounds[i] << endl;
      //std:: cout << "consecutive times out of bounds: " << falseCounter++ << endl;
      return false;
    }
  }
  //falseCounter=0;
  return true;
}

bool MultiGaussian::throwDice(const double* mean) {
  std::vector<double> vec(_SIZE);
  for(unsigned int i=0;i < _SIZE; i++)
    vec[i]=mean[i];
  return throwDice(vec);
}

/*! Extract random variable i; all i together have a multivariate gaussian
 *  distribution with the given covariance matrix.
 */
double MultiGaussian::getRandomValue(int i) const {
  return _generatedValues[i]+_center[i];
}

//! returns a vector of all generated values. \see getRandomValue()
std::vector<double> MultiGaussian::getRandomValues() const {
  std::vector<double> randomVals = _generatedValues;
  for (int i=0; i < _SIZE; ++i) {
    randomVals[i] += _center[i];
  }
  return randomVals;
}


/*!Print some information (usually for debugging purposes)*/
void MultiGaussian::printInfo(std::ostream& o) {
  o.setf(ios::scientific);
  o.precision(5);
    o << "************************************" << std::endl;
    o << "MultiGaussian::printInfo:" << std::endl;
    o << "Eigenvalues (we used scale = " << Scale << "):" << std::endl;
    for(unsigned int i=0; i< _SIZE; i++){
      o << _eigenvalues[i] << "  --> sigma: " << sqrt(_eigenvalues[i]) <<  "  :: eigenvector w/o scale:  " << _eigenvalues[i]/Scale << "  --> sigma: " << sqrt(_eigenvalues[i]/Scale)  << std::endl;
    }
    o << "Eigenvectors:" << std::endl;
    for(unsigned int i=0; i < _SIZE; i++){
        for(unsigned int j=0; j < _SIZE ; j++){
          o << _MasterMatrix[i][j] << "  " ;
        }
        o << std::endl;
    }
    o << "***********************************" << std::endl;
}
/*! A bit more information than from printInfo()*/
void MultiGaussian::printExtInfo(std::ostream& o) {
    printInfo(o);
    o << "EXTENDED information: " << std::endl;
    o << "Size: " << _SIZE << std::endl;
    o << "Bounds: " << std::endl;
    for(unsigned int i=0; i < _SIZE; i++ ) {
        o << _lbounds[i] << "  " << _hbounds[i] << std::endl;
    }
    o << "Generated values currently stored: " << std::endl;
    for(unsigned int i=0; i < _SIZE; i++ ) {
      bool exceed_low = _generatedValues[i]+_center[i] < _lbounds[i];
      bool exceed_high = _generatedValues[i]+_center[i] > _hbounds[i];
      o << "generated: " << _generatedValues[i] << " return: " << (_generatedValues[i]+_center[i]);
      if (exceed_low) o << " (smaller than lower bound)";
      if (exceed_high) o << " (larger than upper bound)";
      o << std::endl;
    }
    o << "************EXTENDED information ends******" << std::endl;
}

/*!This utility function is used only by the member function throwDice()*/
void MultiGaussian::generateAndTransform(){
    double *y = new double[_SIZE];
    for(unsigned int i=0;i<_SIZE;i++) {
       y[i]= math.gaussRnd(0.0,sqrt(_eigenvalues[i])); // generate y random numbers with variances given above
    }
    for(unsigned int i=0; i <_SIZE;i++) {
      _generatedValues[i]=0;
        for(unsigned int j=0; j< _SIZE; j++) {
            _generatedValues[i]+=_MasterMatrix[i][j]*y[j];
        }
    }
    delete[] y;
}

