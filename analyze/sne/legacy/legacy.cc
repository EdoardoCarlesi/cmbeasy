#include "legacy.h"
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

//! reimplemented from base class (i.e. slightly modified) / adapted for legacy data file format
void LegacySNeData::readData(const string& FileName, bool verbose)
{
  ifstream file(FileName.c_str());
  if (!file) {
    stringstream errstrng("");
    errstrng << "SN::readData() Data file " << FileName << " not found";
    throw Bad_Error(errstrng.str());
  }
  string line;
  vector<string> words;
  stringstream str("");
  SNeDataEntry sndata;
  Legacy = new Data("");

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
    str.str(words[3]); str.clear(); str >> sndata.dz;
    str.str(words[4]); str.clear(); str >> sndata.mag;
    str.str(words[5]); str.clear(); str >> sndata.dmag;
    
    DataEntry data = DataEntry(sndata.zcmb,sndata.mag,sndata.dz,sndata.dmag);
   // Maybe it doesn't make too much sense to push_back these data
   // since we already store them into Legacy...
   // But I get some errors otherwise!
    points.push_back(sndata);
    Legacy->points.push_back(data);

    line.clear();
  }
  if (verbose)
    cout << "about to read  " << " from " <<  FileName << endl;
  if (points.size() == 0)
    throw Bad_Error("LegacySN::readData() No data read");
  if (verbose)
    cout << "Read: " << points.size() << " lines from " <<  FileName << endl;
  file.close();
}

double LegacySNeData::chi2(const baseCosmos& cosmos)
 {

  Spline lum(100, "Model::lum");
  for (double z=0; z<=2; z+=0.05) {
  lum.set(z, cosmos.luminosityDistance(z));
  }
  lum.arm();
  double chisq;
  AnalyzeThis *aT;
  Data data("LegacySN");

  data = *Legacy;
  // Calling the old likelihood-computation 
  // routine from the AnalyzeThis class
  chisq=aT->Sn1aCore(lum, data, 5.0);
  
  cout << "legacy chisq = " << chisq << endl;
  
  return chisq;
}

