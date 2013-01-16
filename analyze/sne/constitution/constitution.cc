#include "constitution.h"
#include "analyzethis.h"
#include "data.h"

#include <gsl/gsl_permutation.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>

#include <iostream>

//! reimplemented from base class (i.e. slightly modified) / adapted for constitution data file format
void ConstitutionSNeData::readData(const string& FileName, bool verbose)
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
  Constitution_09 = new Data("");

  //File reading loop
  while(file) {
    getline(file,line);
    if (line[0]=='#' || line[0]==';')
      continue; //Comment line  stringwords(line,words);
    stringwords(line,words);
    if (words.size() < 4)
      continue; //Not enough entries on line

    //Yes, I'm sure there is a more elegant way to do this.
    str.str(words[1]); str.clear(); str >> sndata.zcmb;
    str.str(words[2]); str.clear(); str >> sndata.dz;
    str.str(words[9]); str.clear(); str >> sndata.mag;
    str.str(words[10]); str.clear(); str >> sndata.dmag;
    
    DataEntry data = DataEntry(sndata.zcmb,sndata.mag,sndata.dz,sndata.dmag);
   // Maybe it doesn't make too much sense to push_back these data
   // since we already store them into Constitution_09...
   // But I get some errors otherwise!
    points.push_back(sndata);
    Constitution_09->points.push_back(data);

    line.clear();
  }
  if (verbose)
    cout << "about to read  " << " from " <<  FileName << endl;
  if (points.size() == 0)
    throw Bad_Error("ConstitutionSNeData::readData() No data read");
  if (verbose)
    cout << "Read: " << points.size() << " lines from " <<  FileName << endl;
  file.close();
}

double ConstitutionSNeData::chi2(const baseCosmos& cosmos)
 {

  Spline lum(100, "Model::lum");
  for (double z=0; z<=2; z+=0.05) {
  lum.set(z, cosmos.luminosityDistance(z));
  }
  lum.arm();
  double chisq;
  AnalyzeThis *aT;
  Data data("ConstitutionSNeData09");

  data = *Constitution_09;
  // Calling the old likelihood-computation 
  // routine from the AnalyzeThis class
  chisq=aT->Sn1aCore(lum, data, 5.0);
  
  cout << "constitution chisq = " << chisq << endl;
  
  return chisq;
}

