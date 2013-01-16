#ifndef MCUTILS_H
#define MCUTILS_H

#include "mctaskinfo.h"

#include <algorithm>
#include <numeric>
#include <cmath>
//X #include <vector>
//X #include <functional>
//X #include <iostream>
//X #include <limits>

namespace McUtils
{
  typedef std::vector<std::vector<double> > Matrix;

  std::ostream& operator<<(std::ostream& os, const Matrix& m);
  std::ostream& operator<<(std::ostream& os, const std::vector<double>& v);
  Matrix invertMatrix(const Matrix& matrix);
  Matrix matrixMatrixProd(const Matrix& m1, const Matrix& m2);
  Matrix transposed(const Matrix& m);
  Matrix eigenVectors(const Matrix& m, std::vector<double>* eVals=0);
  std::vector<double> matrixVectorProd(const Matrix& m, const std::vector<double>& v);
  McTaskInfo findBestLogLike(const McTaskInfo& startingPoint);
  Matrix computeHesseMatrix(const McTaskInfo& center, const double bestNegLoglike);
  Matrix estimateCovarianceAroundPoint(McTaskInfo& t, const double bestNegLoglike);

  template<class ForwardIterator>
  void updateCovarianceMatrix(const ForwardIterator& first, const ForwardIterator& last,
                              Matrix& covMat, std::vector<double>* mean=0)
  {
    typename ForwardIterator::value_type valueTypeInstance;
    McTaskInfo check(valueTypeInstance); // if this fails, we have been called on a container
                                         // that doesn't hold things of type McTaskInfo, but
                                         // something else
    int dim=covMat.size();
    for (int i=0; i<dim; ++i) {
      fill(covMat[i].begin(), covMat[i].end(), 0.);
    }

    ForwardIterator it=first;
    unsigned int totalSize=0;
    std::vector<double> average(dim);


    while (it!=last) {
      std::vector<double> paramValues = it->parameterVector();
      for (unsigned int k=0; k<dim; ++k) {
        average[k] += paramValues[k]*it->Multiplicity;
      }
      totalSize += it->Multiplicity; // total weight of all points
      ++it; //next point
    }

    for(unsigned int k=0; k<dim; ++k){ //normalize
      average[k]=average[k]/totalSize;
    }
    if (mean) {
      *mean=average;
    }

    it=first; // reset iterator

    while (it!=last) {
      std::vector<double> paramValues = it->parameterVector();
      for (unsigned int k=0; k<dim; ++k) {
        for (unsigned int j=0; j<dim; ++j) {
          covMat[k][j] += it->Multiplicity
                          *(paramValues[k]-average[k])
                          *(paramValues[j]-average[j]);
        }
      }
      ++it;
    }

    for(unsigned int k=0; k<dim; ++k) {    //normalize
      for(unsigned int j=0; j<dim; ++j) {
        covMat[k][j]/=(totalSize-1.);
      }
    }
  }

  std::vector<double> randomParameterPoint();
  bool insideBounds(const std::vector<double>& params);

  class BoundingEllipse
  {
    public:
      template<class ForwardIterator>
        void setPoints(const ForwardIterator& first, const ForwardIterator& last) {
          using namespace McUtils;
          using namespace std;
          typename ForwardIterator::value_type valueTypeInstance;
          McTaskInfo check(valueTypeInstance); // if this fails, we have been called on a container
          // that doesn't hold things of type McTaskInfo, but
          unsigned int dim = McTaskInfo::paramCount();
          Matrix covMat;
          covMat.resize(dim, vector<double>(dim));
          updateCovarianceMatrix(first, last, covMat, &mMean);
          eVecs = eigenVectors(covMat);
          Matrix eVecsT = transposed(eVecs);

          Matrix diag = matrixMatrixProd(covMat, eVecs);
          diag = matrixMatrixProd(eVecsT, diag);
          /*
             cout << "covMatrix:\n" << covMatrix << endl;
             cout << "evecs:\n" << eVecs << endl;
             cout << "evecsT:\n" << eVecsT << endl;
             cout << "diag:\n" << diag << endl;
             */
          mSqrtDiag.resize(dim, vector<double>(dim));
          for (int i=0; i<diag.size(); ++i) {
            mSqrtDiag[i][i] = sqrt(diag[i][i]);
          }
          std::vector<double> zeroVec(covMat.size());
          for (int i=0; i<covMat.size(); ++i) {
            if (covMat[i]==zeroVec)
              covMat[i][i]=1.;
          }
          mTrafo = matrixMatrixProd(eVecs, mSqrtDiag);
          //cout << "CovMat: " << covMat << endl;
          Matrix invCov = invertMatrix(covMat);
          //std::cout << std::endl << "inv Cov: " << std::endl << invCov << std::endl;
          //std::cout << std::endl << "should be unity: " << std::endl << matrixMatrixProd(invCov, covMat) << std::endl;
          //std::cout << std::endl << "should also be unity: " << std::endl << matrixMatrixProd(covMat, invCov) << std::endl;
          k_max = -numeric_limits<double>::infinity();
          ForwardIterator pointsIterator=first;
          for ( ; pointsIterator!=last; ++pointsIterator) {
            vector<double> pV = pointsIterator->parameterVector();
            transform(pV.begin(), pV.end(), mMean.begin(), pV.begin(), minus<double>());
            vector<double> v = matrixVectorProd(invCov, pV);
            k_max = max(k_max, inner_product(pV.begin(), pV.end(), v.begin(), 0.));
          }
//X           for ( ; pointsIterator!=last; ++pointsIterator) {
//X             vector<double> pV = pointsIterator->parameterVector();
//X             transform(pV.begin(), pV.end(), mMean.begin(), pV.begin(), minus<double>());
//X             //vector<double> v = matrixVectorProd(invCov, pV);
//X             //k_max = max(k_max, inner_product(pV.begin(), pV.end(), v.begin(), 0.));
//X             vector<double> v = matrixVectorProd(mTrafo, pV);
//X             k_max = max(k_max, *max_element(v.begin(), v.end()));
//X           }
          if (std::isinf(k_max))
            k_max=1;
        }

      std::vector<double> fromUnitCube(std::vector<double>& point) {
        using namespace std;
        //cout << "transforming: " << mTrafo <<endl;
        vector<double> res = McUtils::matrixVectorProd(mTrafo, point);
        transform(res.begin(), res.end(), res.begin(), bind1st(multiplies<double>(), sqrt(k_max)));
        //cout << "from point: " << point << " to " << res << endl;
        transform(res.begin(), res.end(), mMean.begin(), res.begin(), plus<double>());
        return res;
      }

      std::vector<std::vector<double> > boundingPoints(int xIndex, int yIndex) {
        using namespace McUtils;
        using namespace std;
        vector<vector<double> > retVec;


        struct FillSphere {
          vector<vector<double> >& mRetVec;
          BoundingEllipse& parent;
          FillSphere(vector<vector<double> >& boundingPoints, BoundingEllipse& p): mRetVec(boundingPoints), parent(p) {
            for (int i=0; i<McTaskInfo::paramCount(); ++i) {
              McTaskInfo info;
              fillAngle(i, info);
            }
          }

          void fillAngle(int no, McTaskInfo& coord) {
              for (double theta = 0; theta <= M_PI; theta += 2.*M_PI/30.) {
                if (no==McTaskInfo::paramCount()-1) {
                  for (double phi = 0; phi <= 2.*M_PI; phi += 2.*M_PI/30.) {
                    vector<double> v = coord.parameterVector();
                    v[no] = phi;
                    vector<double> boundingPoint;
                    for (int i=0; i<McTaskInfo::paramCount(); ++i) {
                      boundingPoint[i] = cos(v[i]);
                      if (i+1==McTaskInfo::paramCount()) {
                        boundingPoint[i] = sin(v[i]);
                      }
                      for (int j=0; j<i; ++j) {
                        boundingPoint[j] *= sin(v[j]);
                      }
                    }
                    mRetVec.push_back(parent.fromUnitCube(v));
                  }
                  return;
                } else {
                  vector<double> v = coord.parameterVector();
                  v[no] = theta;
                  coord.setParameterValues(v);
                  fillAngle(no+1, coord);
                  return;
                }
              }
          }
        };

        FillSphere filler(retVec, *this);
        return retVec;



//X         for (double phi = 0; phi <= 2.*M_PI; phi += 2.*M_PI/30.) {
//X           for (double theta = 0; theta <= M_PI; theta += 2.*M_PI/30.) {
//X             vector<double> coords(McTaskInfo::paramCount());
//X             coords[2] = cos(phi)*sin(theta);
//X             coords[1] = sin(phi)*sin(theta);
//X             coords[0] = cos(theta);
//X             retVec.push_back(fromUnitCube(coords));
//X           }
//X         }
//X         for (double phi = 0; phi <= 2.*M_PI; phi += 2.*M_PI/30.) {
//X           vector<double> coords(McTaskInfo::paramCount());
//X           coords[xIndex]=sin(phi);
//X           coords[yIndex]=cos(phi);
//X           retVec.push_back(fromUnitCube(coords));
//X         }
        //cout << "returning x: " << retVec.size() << endl;
        return retVec;
      }

        private:
      double k_max;
      std::vector<double> mMean;
      McUtils::Matrix mTrafo, mSqrtDiag, eVecs;
      };
};
#endif // MCUTILS_H
