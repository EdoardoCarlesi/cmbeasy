#include "mcutils.h"

#include "mcerror.h"
#include "mcmodel.h"
#include "mclikelihoodcalculator.h"
#include "mcsettings.h"

#include "gsl/gsl_multimin.h"
#include "gsl/gsl_matrix.h"
#include "gsl/gsl_linalg.h"
#include "gsl/gsl_math.h"
#include "gsl/gsl_eigen.h"

#if 0
//#warning Eigen
#include "../eigen2/Eigen/LU"
#include "../eigen2/Eigen/Array"
#endif

#include <iomanip>
#include <numeric>

using namespace std;

namespace McUtils
{
  typedef vector<vector<double> > Matrix;

  ostream& operator<<(ostream& os, const Matrix& m) {
    os << "xxx";
    for (int j=0; j<m.size(); ++j) {
      os << setw(12) << setfill(' ') << j;
    }
    for (int i=0; i<m.size(); ++i) {
      os << endl << i << ": ";
      for (int j=0; j<m.size(); ++j) {
        os << "("<< setw(12) << m[i][j] << ")";
      }
    }
    return os;
  }

  ostream& operator<<(ostream& os, const vector<double>& v) {
    if (v.empty()) {
      os << "(-empty-)";
      return os;
    }
    os << setw(12) << setfill(' ') << '(';
    for (int i=0; i<v.size()-1; ++i) {
      os << v[i] << ", ";
    }
    os << v[v.size()-1] << ')';
    return os;
  }

  Matrix invertMatrix(const Matrix& matrix) {
    int dim = matrix.size();
    if (dim<1) {
      throw McError("McUtils::invertMatrix() - cannot invert matrix with dim<1");
    }
    if (matrix[0].size()!=dim) {
      throw McError("McUtils::invertMatrix() - only quadr. matrices.");
    }

    /*
#warning check this - seems to crash sometimes, and I don't understand the crash
    Eigen::MatrixXd mat(dim, dim);
    for ( int i = 0; i < dim; ++i ) {
      for ( int j = 0; j < dim; ++j ) {
        mat(i,j) = matrix[i][j];
      }
    }
    //
    //if (!mat.lu().isInvertible()) {
    //  cout << "McUtils::invertMatrix() - matrix is not invertible." << endl;
    //}
    Eigen::MatrixXd invMat(dim, dim);
    mat.computeInverse(&invMat);
    Matrix inverse1(dim);
    for ( int i = 0; i < dim; ++i) {
      inverse1[i].resize(dim);
    }
    for ( int i = 0; i < dim; ++i ) {
      for ( int j = 0; j < dim; ++j ) {
        inverse1[i][j] = invMat(i,j);
      }
    }
    return inverse1;
    // */

    gsl_permutation* perm = gsl_permutation_alloc( dim );
    gsl_matrix* originalMatrix = gsl_matrix_alloc( dim, dim );
    gsl_matrix* inverseMatrix = gsl_matrix_alloc( dim, dim );

    for ( int i = 0; i < dim; ++i ) {
      for ( int j = 0; j < dim; ++j ) {
        gsl_matrix_set( originalMatrix, i, j, matrix[i][j] );
      }
    }


    int status, signum;
    status = gsl_linalg_LU_decomp(originalMatrix, perm, &signum);
    if (!status)
      gsl_linalg_LU_invert(originalMatrix, perm, inverseMatrix);
    else
      std::cerr << "LU decomp. failed" << endl;
    if (status) {
      std::cerr<< "failed, gsl_errno=" << status << ": " << gsl_strerror (status) << std::endl;
    }

    Matrix inverse(dim);
    for ( int i = 0; i < dim; ++i) {
      inverse[i].resize(dim);
    }


    for ( int i = 0; i < dim; ++i ) {
      for ( int j = 0; j < dim; ++j ) {
        inverse[i][j] = gsl_matrix_get( inverseMatrix, i, j );
      }
    }
    return inverse;
  }

  std::vector<double> matrixVectorProd(const Matrix& m, const std::vector<double>& v) {
    using namespace std;
    vector<double> prod(m.size());
    for (int i=0; i<m.size(); ++i) {
      prod[i] = inner_product(m[i].begin(), m[i].end(), v.begin(), 0.);
    }
    return prod;
  }

  Matrix matrixMatrixProd(const Matrix& m1, const Matrix& m2) {
    using namespace std;
    struct ExtractColumns {
      const Matrix& matrix;
      ExtractColumns(const Matrix& m): matrix(m) {}
       vector<double> operator()(int i) {
         vector<double> v(matrix.size());
         for (int c=0; c<matrix.size(); ++c) {
           v[c] = matrix[c][i];
         }
         return v;
       }
    };
    Matrix m(m1.size(), vector<double>(m1.size()));
    ExtractColumns m2Col(m2);
    for (int c=0; c<m1.size(); ++c) {
      vector<double> v = m2Col(c);
      for (int d=0; d<m1.size(); ++d) {
        m[d][c] = inner_product(m1[d].begin(), m1[d].end(), v.begin(), 0.);
      }
    }
    return m;
  }

  Matrix eigenVectors(const Matrix& m, vector<double>* eVals) {
    unsigned int dim = m.size();
    Matrix eVecs(dim, vector<double>(dim));
    if (eVals) {
      eVals->resize(dim);
    }
    double *tempArray = new double[dim*dim];
    int k=0;
    for(unsigned int j=0; j< dim; j++)    {            // copy into tempArray
      for(unsigned int i=0; i< dim; i++) {
        tempArray[k]=m[i][j];
        ++k;
      }
    }

    // from GNU Scientific Library
    gsl_eigen_symmv_workspace * w = gsl_eigen_symmv_alloc (dim);
    gsl_matrix_view mw = gsl_matrix_view_array (tempArray, dim, dim);
    gsl_vector *eval = gsl_vector_alloc (dim);
    gsl_matrix *evec = gsl_matrix_alloc (dim, dim);
    gsl_eigen_symmv (&mw.matrix, eval, evec, w);
    //gsl_eigen_symmv_sort (eval, evec, GSL_EIGEN_SORT_ABS_ASC);  //sort in ascending order
    gsl_eigen_symmv_free (w);
    for(unsigned int i=0; i<dim; ++i) {
      if (eVals) {
        eVals->at(i) = gsl_vector_get (eval, i);
      }
      for (unsigned int j=0;j<dim; ++j) {
        eVecs[i][j]=gsl_matrix_get(evec,i,j);
      }
    }

    gsl_vector_free(eval);   // free memory
    gsl_matrix_free(evec);
    delete[] tempArray;
    return eVecs;
  }

  Matrix transposed(const Matrix& m) {
    const unsigned int rows=m.size();
    const unsigned int cols=m[0].size();

    Matrix t(cols, vector<double>(rows));
    for (unsigned int i=0; i<rows; ++i) {
      for (unsigned int j=0; j<cols; ++j) {
        t[j][i]=m[i][j];
      }
    }
    return t;
  }

  double negLogLikeFct(McTaskInfo& t) {
    McModel& model = McSettings::self()->model();
    static McLikelihoodCalculator calcLike;

    model.setParameters(t);
    model.compute();
    calcLike.computeLogLike(model, t);
    return -t.totalLogLike();
  }

  double negLogLikeFct(const gsl_vector *v, void *params) {
    static McTaskInfo t;
    for (int i = 0; i<t.paramCount(); ++i) {
      t.setParameterValue(i, gsl_vector_get(v, i));
    }
    return negLogLikeFct(t);
  }

  double safeNegLogLikeFct(McTaskInfo& t) {
    double negLogLike;
    try {
      negLogLike = negLogLikeFct(t);
    } catch (Bad_Error e) {
      negLogLike = -1e100;
    }
    if (isnan(negLogLike) || isinf(negLogLike)) {
      negLogLike = -1e100;
    }
    return negLogLike;
  }

  double safeNegLogLikeFct(const gsl_vector *v, void *params) {
    static McTaskInfo t;
    for (int i = 0; i<t.paramCount(); ++i) {
      t.setParameterValue(i, gsl_vector_get(v, i));
    }
    return safeNegLogLikeFct(t);
  }

  McTaskInfo findBestLogLike(const McTaskInfo& startingPoint) {
    vector<double> startVector = startingPoint.parameterVector();
    size_t dim = startingPoint.paramCount();

    const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex;
    gsl_multimin_fminimizer *s = gsl_multimin_fminimizer_alloc (T, dim);

    double* par = 0;

    gsl_multimin_function function;

    function.n = dim;
    function.f = &safeNegLogLikeFct;
    function.params = par;

    gsl_vector *x = gsl_vector_alloc(dim);
    gsl_vector *initialStepSize = gsl_vector_alloc(dim);
    vector<double> initialSigma = McTaskInfo::initialParameterSigmas();
    for (int i = 0; i < startingPoint.paramCount(); ++i) {
      gsl_vector_set(x, i, startVector[i]);
      gsl_vector_set(initialStepSize, i, initialSigma[i]);
    }


    gsl_multimin_fminimizer_set(s, &function, x, initialStepSize);

    size_t iter = 0;
    int status;
    McTaskInfo t;
    do
    {
      iter++;
      status = gsl_multimin_fminimizer_iterate(s);

      if (status)
        break;

      double size = gsl_multimin_fminimizer_size (s);
      status = gsl_multimin_test_size (size, 1e-2);

      cout << iter << ": " << prettyPrint << *McSettings::self()->model().currentParameters() << endl;;

      if (status == GSL_SUCCESS) {
        cout << ("Minimum found at:\n");

        // the minimum isn't necessarily the last point we looked at, so we can't use
        // model().currentParameters()
        for (int i = 0; i < startingPoint.paramCount(); ++i) {
          t.setParameterValue(i, gsl_vector_get(s->x, i));
        }
        McModel& model = McSettings::self()->model();
        McLikelihoodCalculator calcLike;
        model.setParameters(t);
        model.compute();
        calcLike.computeLogLike(model, t);
      }
    }
    while (status == GSL_CONTINUE && iter < 500);

    gsl_multimin_fminimizer_free (s);
    gsl_vector_free (x);
    gsl_vector_free (initialStepSize);

    return t;
  }

  Matrix computeHesseMatrix(const McTaskInfo& center, const double bestNegLoglike) {

    vector<double> initialSigma, centerPoint, upperBound, lowerBound,
                   upperOneSigBound, lowerOneSigBound;
    upperBound = McTaskInfo::upperParameterBounds();
    lowerBound = McTaskInfo::lowerParameterBounds();
    initialSigma = McTaskInfo::initialParameterSigmas();
    centerPoint = center.parameterVector();
    McTaskInfo trialPoint;
    for (int i=0; i<center.paramCount(); ++i) {
      trialPoint = center;
      double upper = min(upperBound[i], centerPoint[i]+initialSigma[i]);
      double lower = max(lowerBound[i], centerPoint[i]-initialSigma[i]);

      int maxSteps = 20;
      int stopCountUp = 0;
      trialPoint.setParameterValue(i, upper);
      while (fabs(negLogLikeFct(trialPoint)-bestNegLoglike)<0.6) {
        upper+=initialSigma[i];
        trialPoint.setParameterValue(i, upper);
        if (++stopCountUp > maxSteps) {
          cout << "broke out upper for param: " << i << endl;
          break;
        }
      }
      int stopCountDown=0;
      trialPoint.setParameterValue(i, lower);
      while (fabs(negLogLikeFct(trialPoint)-bestNegLoglike)<0.6) {
        lower-=initialSigma[i];
        trialPoint.setParameterValue(i, lower);
        if (++stopCountDown > maxSteps) {
          cout << "broke out lower for param: " << i << endl;
          break;
        }
      }
      if (stopCountUp<=maxSteps && stopCountDown <=maxSteps) {
        struct NegLogLikeFunctor: public Mathobject {
          NegLogLikeFunctor(McTaskInfo c, const int i): mCenter(c), mIdx(i) {}
          double negLogLike(const double param) const {
            McTaskInfo& c = const_cast<McTaskInfo&>(mCenter);
            c.setParameterValue(mIdx, param);
            double retVal = negLogLikeFct(c);
            cout << "checking " << c.parameterVector()[mIdx] << " got " << retVal << endl;
            return retVal;
          }
          McTaskInfo mCenter;
          const int mIdx;
        };
        NegLogLikeFunctor f(center, i);
        try {
          double ret = Miscmath::zbrent((moSingle)&NegLogLikeFunctor::negLogLike,
              f, centerPoint[i], upper, 1e-1, bestNegLoglike+0.5);
          upperOneSigBound.push_back(ret);
          cout << "found upper 1 sigma bound for " << i << " " << upperOneSigBound[i] << endl;
          ret = Miscmath::zbrent((moSingle)&NegLogLikeFunctor::negLogLike,
              f, centerPoint[i], lower, 1e-1, bestNegLoglike+0.5);
          lowerOneSigBound.push_back(ret);
        } catch (Bad_Error) { // most likely something like "too many steps in zbrent"
          upperOneSigBound.push_back(centerPoint[i]+0.75*(upperBound[i]-centerPoint[i]));
          lowerOneSigBound.push_back(centerPoint[i]+0.75*(lowerBound[i]-centerPoint[i]));
        }
        cout << "found lower 1 sigma bound for " << i << " " << lowerOneSigBound[i] << endl;
      } else {
        upperOneSigBound.push_back(centerPoint[i]+0.75*(upperBound[i]-centerPoint[i]));
        lowerOneSigBound.push_back(centerPoint[i]+0.75*(lowerBound[i]-centerPoint[i]));
      }

    }
    cout << endl << "done finding one sigma bounds" << endl;
      Matrix deltaLnLike;
      for (int i=0; i<center.paramCount(); ++i) {
        vector<double> row(center.paramCount(), 0.);
        deltaLnLike.push_back(row);
        double iStepSize = (upperOneSigBound[i]-lowerOneSigBound[i])/2.;
        deltaLnLike[i][i] = 1./(iStepSize*iStepSize)*1.;
        for (int j=0; j<i; ++j) {
          trialPoint = center;
          trialPoint.setParameterValue(i, upperOneSigBound[i]);
          trialPoint.setParameterValue(j, upperOneSigBound[j]);
          double upUp = negLogLikeFct(trialPoint);
          trialPoint.setParameterValue(j, lowerOneSigBound[j]);
          double upDown = negLogLikeFct(trialPoint);
          trialPoint.setParameterValue(i, lowerOneSigBound[i]);
          double downDown = negLogLikeFct(trialPoint);
          trialPoint.setParameterValue(j, upperOneSigBound[j]);
          double downUp = negLogLikeFct(trialPoint);

          double jStepSize = (upperOneSigBound[j]-lowerOneSigBound[j])/2.;
          deltaLnLike[i][j] = 1./(iStepSize*jStepSize*4.)*(upUp+downDown-upDown-downUp);
          deltaLnLike[j][i] = deltaLnLike[i][j];
      }
      cout << "done with param " << i << endl;
    }
    return deltaLnLike;
  }

  Matrix estimateCovarianceAroundPoint(McTaskInfo& t, const double bestNegLoglike) {
    Matrix hesse = McUtils::computeHesseMatrix(t, bestNegLoglike);
    cout << "=========== found hess: ============\n";
    cout << hesse << endl;;
    Matrix covariance = invertMatrix(hesse);
    cout << "=========== found convariance ===========\n";
    cout << covariance << endl;
    return covariance;
  }

  vector<double> randomParameterPoint()
  {
    int paramCount=McTaskInfo::paramCount();
    vector<double> lowerBounds, upperBounds;
    lowerBounds=McTaskInfo::lowerParameterBounds();
    upperBounds=McTaskInfo::upperParameterBounds();
    vector<double> randomVals(paramCount);
    for (int i=0; i<paramCount; ++i) {
      randomVals[i]=lowerBounds[i] + Miscmath::posRnd(upperBounds[i]-lowerBounds[i]);
    }
    return randomVals;
  }

  bool insideBounds(const vector<double>& params)
  {
    unsigned int size=params.size();
    if (size!=McTaskInfo::paramCount()) {
      throw McError("McUtils::insideBounds() - size mismatch.");
    }

    const vector<double>& lower=McTaskInfo::lowerParameterBounds();
    const vector<double>& upper=McTaskInfo::upperParameterBounds();
    for (int i=0; i<size; ++i) {
      if (params[i]<lower[i] || params[i]>upper[i]) {
        return false;
      }
    }
    return true;
}

};


