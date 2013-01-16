// vim:ts=2:sw=2:et
#include "union2.h"

#include <gsl/gsl_permutation.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>

//! reimplemented from base class (i.e. slightly modified) / adapted for union data file format
void Union2SNeData::readData(const string& FileName, bool verbose)
{
  ifstream file(FileName.c_str());
  if (!file) {
    stringstream errstrng("");
    errstrng << "SNeData::readData() Data file " << FileName << " not found";
    throw Bad_Error(errstrng.str());
  }
  string line;
  vector<string> words;
  stringstream str("");
  SNeDataEntry sndata;
  //File reading loop
  while(file) {
    getline(file,line);
    if (line[0]=='#' || line[0]==';')
      continue; //Comment line  stringwords(line,words);
    stringwords(line,words);
    if (words.size() < 4)
      continue; //Not enough entries on line
    sndata.name=words[0];
   //Yes, I'm sure there is a more elegant way to do this.
    str.str(words[1]); str.clear(); str >> sndata.zcmb;
    str.str(words[2]); str.clear(); str >> sndata.mag;
    str.str(words[3]); str.clear(); str >> sndata.dmag;
    points.push_back(sndata);
    line.clear();
  }
  if (verbose)
    cout << "about to read  " << " from " <<  FileName << endl;
  if (points.size() == 0)
    throw Bad_Error("Union2SNeData::readData() No data read");
  if (verbose)
    cout << "Read: " << points.size() << " lines from " <<  FileName << endl;
  file.close();
}


void Union2SNeData::readCovMat(const string& fileNameCovMatSys, const string& fileNameCovMatNoSys,
                                bool includeSysError)
{
  string fileName;
  if(includeSysError) {
    fileName = fileNameCovMatSys;
  } else {
    fileName = fileNameCovMatNoSys;
  }
  ifstream file(fileName.c_str());

  gsl_matrix *m = gsl_matrix_alloc (557, 557);
  //Matrix matrix;
  //matrix.resize(557, vector<double>(557, 0.));

  if (!file)  {
   stringstream errstrng("");
   errstrng << "SNeData::readMat() Data file " << fileName << " not found";
   throw Bad_Error(errstrng.str());
  }

  string line;
  unsigned int i=0;
  while (getline(file, line)) {
    if (line.size()==0 || line[0]=='#') {
      continue;
    }
    istringstream lineStream(line);
    unsigned int j=0;
    double x;
    while (lineStream>>x) {
      gsl_matrix_set (m, i, j, x);
      //matrix[i][j] = x;
      ++j;
    }
    if (j!=557 || i>556) {
        cout << line << i << " x " << j << endl;
        throw Bad_Error("Union2SNeData::readCovMat() - expected 557x557 matrix");
    }
    ++i;
  }

  gsl_permutation *p=gsl_permutation_alloc (557);
  gsl_matrix *inverse=gsl_matrix_alloc (557, 557);
  //gsl_matrix *LU=gsl_matrix_alloc (557, 557);
  int s;
  gsl_linalg_LU_decomp (m, p, &s);
  gsl_linalg_LU_invert (m, p, inverse);
  gsl_permutation_free (p);
  gsl_matrix_free (m);

  sum=0;

  mInverseMatrix.resize(557, vector<double>(557, 0.));
  for (unsigned int i=0; i<557; i++) {
    for(unsigned int j=0; j<557; j++) {
      sum+=gsl_matrix_get(inverse,i,j);
      mInverseMatrix[i][j]=gsl_matrix_get(inverse,i,j);
    }
  }
  //mInverseMatrix = McUtils::invertMatrix(matrix);
}

double Union2SNeData::chi2(const baseCosmos& cosmos)
{
  if (!mInitialized) {
      init();
  }
  Vec diffs(557);
  int count=0;
  int i,j;

  vector<SNeDataEntry>::const_iterator sn = points.begin();
  for ( ; sn != points.end(); ++sn) {
    diffs[count++]=5.*log10(cosmos.luminosityDistance(sn->zcmb))+25.-sn->mag;
  }

  Vec prod(557);
  for(i=0; i<557; i++) {
    prod[i] = 0.;
    for(j = 0; j<557; j++) {
      prod[i] += mInverseMatrix[i][j]*diffs[j];
    }
  }

  double dot1=0;

  for(count=0; count<557; count++) {
    dot1 += prod[count]*diffs[count];
  }

  double sum1=0;
  for(count=0; count<557; count++) {
    sum1+=prod[count];
  }
  double chisq=dot1-(pow(sum1,2)/sum);

  cout << "union2 chisq = " << chisq << endl;
  return chisq;
}

