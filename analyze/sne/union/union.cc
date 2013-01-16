// vim:ts=2:sw=2:et
#include "union.h"

#include <gsl/gsl_permutation.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>

//! reimplemented from base class (i.e. slightly modified) / adapted for union data file format
void UnionSNeData::readData(const string& FileName, bool verbose)
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
    throw Bad_Error("UnionSNeData::readData() No data read");
  if (verbose)
    cout << "Read: " << points.size() << " lines from " <<  FileName << endl;
  file.close();
}


void UnionSNeData::readCovMat(const string& fileNameCovMatSys, const string& fileNameCovMatNoSys,
                                bool includeSysError)
{
  string fileName;
  if(includeSysError) {
    fileName = fileNameCovMatSys;
  } else {
    fileName = fileNameCovMatNoSys;
  }
  ifstream file(fileName.c_str());

  gsl_matrix *m = gsl_matrix_alloc (307, 307);
  //Matrix matrix;
  //matrix.resize(307, vector<double>(307, 0.));

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
    if (j!=307 || i>306) {
        cout << line << i << " x " << j << endl;
        throw Bad_Error("UnionSNeData::readCovMat() - expected 306x306 matrix");
    }
    ++i;
  }

  gsl_permutation *p=gsl_permutation_alloc (307);
  gsl_matrix *inverse=gsl_matrix_alloc (307, 307);
  //gsl_matrix *LU=gsl_matrix_alloc (307, 307);
  int s;
  gsl_linalg_LU_decomp (m, p, &s);
  gsl_linalg_LU_invert (m, p, inverse);
  gsl_permutation_free (p);
  gsl_matrix_free (m);

  sum=0;

  mInverseMatrix.resize(307, vector<double>(307, 0.));
  for (unsigned int i=0; i<307; i++) {
    for(unsigned int j=0; j<307; j++) {
      sum+=gsl_matrix_get(inverse,i,j);
      mInverseMatrix[i][j]=gsl_matrix_get(inverse,i,j);
    }
  }
  //mInverseMatrix = McUtils::invertMatrix(matrix);
}

double UnionSNeData::chi2(const baseCosmos& cosmos)
{
  if (!mInitialized) {
      init();
  }
  Vec diffs(307);
  int count=0;
  int i,j;

  vector<SNeDataEntry>::const_iterator sn = points.begin();
  for ( ; sn != points.end(); ++sn) {
    diffs[count++]=5.*log10(cosmos.luminosityDistance(sn->zcmb))+25.-sn->mag;
  }

  Vec prod(307);
  for(i=0; i<307; i++) {
    prod[i] = 0.;
    for(j = 0; j<307; j++) {
      prod[i] += mInverseMatrix[i][j]*diffs[j];
    }
  }

  double dot1=0;

  for(count=0; count<307; count++) {
    dot1 += prod[count]*diffs[count];
  }

  double sum1=0;
  for(count=0; count<307; count++) {
    sum1+=prod[count];
  }
  double chisq=dot1-(pow(sum1,2)/sum);

  //cout << "union chisq = " << chisq << endl;
  return chisq;
}

