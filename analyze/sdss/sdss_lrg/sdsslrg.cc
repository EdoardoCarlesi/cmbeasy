#include "sdsslrg.h"

#include "distcosmos.h"

#include <iterator>
#include <numeric>
#include <algorithm>

SdssLrg::SdssLrg(const baseCosmos& c, std::string dataDir)
                     : mCosmos(c), mDataDir(dataDir), mMaxBand(12)
{
  std::string filename = dataDir+"sdss_lrg_kbands.txt";
  ifstream kbandFile(filename.c_str());
  if (!kbandFile) throw Bad_Error("SdssLrg::SdssLrg() Could not open sdss_lrg data file: " + filename);

  std::copy(std::istream_iterator<double>(kbandFile), std::istream_iterator<double>(),
                                                            std::back_inserter(mkvals));

  filename = dataDir+"sdss_lrg_measurements.txt";
  ifstream measurements(filename.c_str());
  if (!measurements) throw Bad_Error("Could not open sdss_lrg data file: " + filename);

  int line = 0;
  while ( line++ < 2) measurements.ignore(4096, '\n'); //skip comments
  Band b;
  while ( measurements >> b.keff >> b.klow >> b.khigh >> b.obs >> b.sdev )
      mBands.push_back(b);

  mkbands = mkvals.size();

  filename = dataDir+"sdss_lrg_windows.txt";
  ifstream windows(filename.c_str());
  if (!windows) throw Bad_Error("Could not open sdss_lrg data file: " + filename);

  for (unsigned int rowNo = 0; rowNo < mBands.size(); ++rowNo)
  {
    std::vector<double> row;
    for (unsigned int col = 0; col < mkbands; ++col)
    {
      double val;
      windows >> val;
      row.push_back(val);
    }
    mMatrix.push_back(row);
  }
}

double SdssLrg::aFactor()
{
  // this should be cached, but it's not too much of a bottleneck
  distCosmos fiducialCosmos;
  fiducialCosmos.setOmega_m(0.25);
  fiducialCosmos.setOmega_v(0.75);
  fiducialCosmos.history();

  double z = 0.35;
  double d_a = mCosmos.angulardiameterDistance(z)*mCosmos.h();
  double d_afiducial = fiducialCosmos.angulardiameterDistance(z)*fiducialCosmos.h();

  return pow( pow(d_a/d_afiducial, 2) * fiducialCosmos.Z2H(z) / fiducialCosmos.h()
                        / const_cast<baseCosmos&>(mCosmos).Z2H(z) * mCosmos.h(), 1./3. );
}

std::vector<double> SdssLrg::shiftedPower(Spline& pk, double bias)
{
  std::vector<double> theoPower;
  std::vector<double>::const_iterator it,end;
  it = mkvals.begin();
  end = mkvals.end();
  double a = aFactor();
  //cout << "a scaling is: " << a << endl;
  for ( ; it != end; ++it) {
      double p = pk(*it*a);
      p *= bias*bias;
      p /=a*a*a;
      theoPower.push_back(p);
  }
  return theoPower;
}


double SdssLrg::chi2(Spline& pk, double bias)
{
  std::vector<double> theoPower = shiftedPower(pk, bias);

  Matrix m = mMatrix;
  for (unsigned int i = 0; i < mBands.size(); ++i)
      for (unsigned int j = 0; j < mkbands; ++j)
          m[i][j] /= mBands[i].sdev;

  std::vector<double> matrixp;
  for (unsigned int i = 0; i <= mMaxBand; ++i){
      matrixp.push_back( std::inner_product(m[i].begin(), m[i].end(), theoPower.begin(), 0.));
  }

  double chi2=0;
  for (unsigned int i=0; i <= mMaxBand; ++i) {
      chi2 += pow(mBands[i].obs/mBands[i].sdev-matrixp[i], 2);
  }
  return chi2;
}

double SdssLrg::chi2bestBias(Spline& pk, double firstBias, double lastBias, unsigned int steps)
{
    double db = (lastBias-firstBias)/(double)steps;
    double b = firstBias;
    std::vector<double> chi2vec;
    while (b <= lastBias){
        chi2vec.push_back(chi2(pk, b));
        b += db;
    }
    return *std::min_element(chi2vec.begin(), chi2vec.end());
}

double SdssLrg::chi2margedOverb(Spline& pk)
{
  std::vector<double> theoPower = shiftedPower(pk, 1.);

  // marginalize analytically over b^2 over bias, flat prior on b^2,
  // pretty much exactly the way CAMB does it
  std::vector<double> w;
  BandIterator bandsEnd = mBands.end();
  for ( BandIterator it = mBands.begin(); it != bandsEnd; ++it)
      w.push_back(1./pow(it->sdev, 2));

  std::vector<double> matrixp;
  for (unsigned int i = 0; i <= mMaxBand; ++i)
      matrixp.push_back( std::inner_product(mMatrix[i].begin(), mMatrix[i].end(), theoPower.begin(), 0.));

  //std::copy(matrixp.begin(), matrixp.end(), std::ostream_iterator<double>(cout, " "));

  double norm;
  for (unsigned int i = 0; i <= mMaxBand; ++i)
      norm += matrixp[i]*matrixp[i]*w[i];

  double tmp;
  for (unsigned int i = 0; i <= mMaxBand; ++i)
      tmp += matrixp[i]*mBands[i].obs*w[i]/norm;


  double chi2=0;
  for (unsigned int i=0; i <= mMaxBand; ++i)
      chi2 += mBands[i].obs*(mBands[i].obs-matrixp[i]*tmp)*w[i];
  chi2 += log(norm);

  cout << "tmp: " << tmp << endl;
  cout << "norm: " << norm << endl;
  return chi2;
}

typedef std::vector<std::vector<double> > Matrix;
static void invert2x2Matrix( Matrix& m)
{
    Matrix mInv; for (int i = 0; i <=1; ++i) mInv.push_back(vector<double>());

    double det = m[0][0]*m[1][1]-m[0][1]*m[1][0];
    mInv[0][0] = m[1][1]/det;
    mInv[1][1] = m[0][0]/det;
    mInv[0][1] = -m[1][0]/det;
    mInv[1][0] = -m[0][1]/det;
    m[0][0] = mInv[0][0];
    m[1][1] = mInv[1][1];
    m[0][1] = mInv[0][1];
    m[1][0] = mInv[1][0];
}


double SdssLrg::chi2margedOverb2Q(Spline& pk, double Ag)
{
    using namespace std;

    double a = aFactor();
    vector<double> scaledkvals = mkvals;
    vector<double> onepAgkscaled = mkvals;
    transform(onepAgkscaled.begin(), onepAgkscaled.end(), onepAgkscaled.begin(), bind1st(multiplies<double>(), Ag));
    transform(onepAgkscaled.begin(), onepAgkscaled.end(), onepAgkscaled.begin(), bind1st(plus<double>(), 1.));
    transform(scaledkvals.begin(), scaledkvals.end(), scaledkvals.begin(), bind1st(multiplies<double>(), a));

    vector<double> theoPower = shiftedPower(pk, 1.);
    vector<double> mpk_Pth = theoPower;
    transform(mpk_Pth.begin(), mpk_Pth.end(), onepAgkscaled.begin(), mpk_Pth.begin(), divides<double>());

    vector<double> mpk_k2 = scaledkvals;
    transform(mpk_k2.begin(), mpk_k2.end(), mpk_k2.begin(), mpk_k2.begin(), multiplies<double>());
    transform(mpk_k2.begin(), mpk_k2.end(), mpk_Pth.begin(), mpk_k2.begin(), multiplies<double>());

    vector<double> mpk_WPth;
    for (unsigned int i = 0; i <= mMaxBand; ++i)
        mpk_WPth.push_back( inner_product(mMatrix[i].begin(), mMatrix[i].end(), mpk_Pth.begin(), 0.));

    vector<double> mpk_WPth_k2;
    for (unsigned int i = 0; i <= mMaxBand; ++i)
        mpk_WPth_k2.push_back( inner_product(mMatrix[i].begin(), mMatrix[i].end(), mpk_k2.begin(), 0.));

    vector<double> w;
    BandIterator bandsEnd = mBands.end();
    for ( BandIterator it = mBands.begin(); it != bandsEnd; ++it)
        w.push_back(1./pow(it->sdev, 2));

    vector<double> covdat,covth,covth_k2;
    for ( unsigned int i = 0; i <= mMaxBand; ++i){
        covdat.push_back(mBands[i].obs*w[i]);
        covth.push_back(mpk_WPth[i]*w[i]);
        covth_k2.push_back(mpk_WPth_k2[i]*w[i]);
    }

    Matrix mat; for (int i = 0; i <=1; ++i) mat.push_back(vector<double>());
    mat[0][0] = inner_product(covth.begin(), covth.end(), mpk_WPth.begin(), 0.);
    mat[1][1] = inner_product(covth_k2.begin(), covth_k2.end(), mpk_WPth_k2.begin(), 0.);
    mat[0][1] = mat[1][0] = inner_product(covth.begin(), covth.end(), mpk_WPth_k2.begin(), 0.);
    double lnLike = log( mat[0][0]*mat[1][1]-mat[0][1]*mat[1][0]);

    invert2x2Matrix(mat);
    vector<double> vec(2);
    vec[0] = inner_product(covdat.begin(), covdat.end(), mpk_WPth.begin(), 0.);
    vec[1] = inner_product(covdat.begin(), covdat.end(), mpk_WPth_k2.begin(), 0.);

    double sum1 = 0;
    for ( unsigned int i = 0; i <= mMaxBand; ++i) sum1 += mBands[i].obs*covdat[i];
    double sum2 = (mat[0][0]*vec[0] + mat[0][1]*vec[1])*vec[0];
    sum2 += (mat[1][0]*vec[0] + mat[1][1]*vec[1])*vec[1];

    return sum1-sum2+lnLike;
}

