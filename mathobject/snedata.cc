#include <string>
#include <fstream>
#include <cmath>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <map>
#include <stdexcept>

#include "global.h"
#include "snedata.h"

SNeData::SNeData() : Name("") {}

SNeData::SNeData(const string& name) : Name(name) {}

SNeData::SNeData(const string& name, 
		 const string &fileName) : Name(name) {
  readData(fileName,false);
}


//Cheerfully stolen from Robert Knop
void SNeData::stringwords(const string &ins,vector<string> &words) const {
  string s,tmp;
  unsigned int i,p;
  int first,last;

  s = ins;

  // Trim spaces from beginning and end

  first=s.find_first_not_of(" ");
  if (first==-1) {
    s="";
  } else {
    last=s.find_last_not_of(" ");
    s=s.substr(first,last-first+1);
  }
  words.clear();

  p=s.find_first_not_of(" \t\r\n");
  if (p>=s.size()) return;
  while ((i=s.find_first_of(" \t\r\n",p))) {
    tmp=s.substr(p,i-p);
    words.push_back(tmp);
    p=s.find_first_not_of(" \t\r\n",i);
    if (p>=s.size()) return;
  }
  tmp=s.substr(p);
  words.push_back(tmp);
}


void SNeData::zcmbsort() {
  if (points.size() == 0) return;
  sort ( points.begin(), points.end(), SNeSortByZcmb() );
}

/*!
The data file should have the format
snname zcmb zhel dz widthpar dwidthpar colourpar dcolourpar cov_mag_width cov_mag_colourpar cov_widthpar_colourpar
Lines that begin with # are ignored

*/

void SNeData::readData(const string& FileName, bool verbose) { 
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
    if (line[0]=='#' || line[0]==';') continue; //Comment line
    stringwords(line,words);
    if (words.size() < 13) continue; //Not enough entries on line
    
    sndata.name = words[0];

    //Yes, I'm sure there is a more elegant way to do this.
    str.str(words[1]); str.clear(); str >> sndata.zcmb;
    str.str(words[2]); str.clear(); str >> sndata.zhel;
    str.str(words[3]); str.clear(); str >> sndata.dz;
    str.str(words[4]); str.clear(); str >> sndata.mag;
    str.str(words[5]); str.clear(); str >> sndata.dmag;
    str.str(words[6]); str.clear(); str >> sndata.widthpar;
    str.str(words[7]); str.clear(); str >> sndata.dwidthpar;
    str.str(words[8]); str.clear(); str >> sndata.colourpar;
    str.str(words[9]); str.clear(); str >> sndata.dcolourpar;
    str.str(words[10]); str.clear(); str >> sndata.cov_mag_widthpar;
    str.str(words[11]); str.clear(); str >> sndata.cov_mag_colourpar;
    str.str(words[12]); str.clear(); str >> sndata.cov_widthpar_colourpar;

    points.push_back(sndata);
  }
  
  if (points.size() == 0) throw Bad_Error("SNeData::readData() No data read");

  if (verbose) cout << "Read: " << points.size() << " lines from " << 
		 FileName << endl;

  file.close();
}

/*!
  Uses the same format as SNeData::readData.
  \param outfile File to write to
*/
void SNeData::writeData(const string& outfile) const {  
  //Since we are doing formatted output, and the c++ style formatted
  // output sucks out loud, use C style

  FILE *fp;
  fp = fopen(outfile.c_str(),"w");
  if (!fp) {
    stringstream errstrng("");
    errstrng <<  "SNeData::writeFile() Couldn't open file " << outfile << 
      " for writing";
    throw Bad_Error(errstrng.str());
  }

  string hdfmt="%-10s %-8s %-8s %-8s %-8s %-8s %-8s %-8s %-8s %-8s ";
  hdfmt += "%-11s %-11s %-11s";
  string fmt = "%-10s %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f ";
  fmt += "%11.8f %11.8f %11.8f";

  //Write the header
  fprintf(fp,hdfmt.c_str(),"name","zcmb","zhel","dz","mag","dmag",
	  "width","dwidth","colour","dcolour","cov_m_w","cov_m_c",
	  "cov_w_c");

  //And the data
  for (vector<SNeDataEntry>::const_iterator sn = points.begin(); 
       sn != points.end(); ++sn)
    fprintf(fp,fmt.c_str(),sn->name.c_str(),sn->zcmb,sn->zhel,sn->dz,
	    sn->mag,sn->dmag,sn->widthpar,sn->dwidthpar,
	    sn->colourpar,sn->dcolourpar,sn->cov_mag_widthpar,
	    sn->cov_mag_colourpar,sn->cov_widthpar_colourpar);

  fclose(fp);
}


/*!
  \param a List of SNe to append to this one.
  Duplicates are not prevented by this operator.
*/
SNeData& SNeData::operator+=(const SNeData& a) {

  if (&a == this) {
    //Self addition.  It's an interesting question if
    // we should respect the const reference to a or not --
    // here I do the addition anyways
    int nsn = points.size();
    points.reserve( 2 * nsn );
    //This is a case where using iterators could be a problem
    for (int i = 0; i < nsn; ++i) points[i+nsn] = points[i];
  } else {
    points.reserve( points.size() + a.points.size() );
    for (SNeData::const_iterator i = a.begin(); i != a.end(); ++i)
      points.push_back( *i );
  }

  return *this;
}

/*!
  \param indxarr Array of indicies into old array.  
  \returns A new SNeData with only the specified objects present

  No bounds checking is performed, so it is up to the caller
  to ensure that the indicies in indxarr are valid.
*/
SNeData SNeData::operator[](const vector<int>& indxarr) const {
  int nelements = indxarr.size();

  if (nelements == 0) return SNeData(Name);

  SNeData retval(Name);

  retval.points.resize( nelements );
  
  for (int i = 0, j = 0; j < nelements; ++j)
    retval.points[i++] = points[indxarr[j]];

  return retval;
}

/*!
  \param namearr Vector of strings corresponding to names
    of SN
  \returns A new SNeData with only the specified objects present

  Error checking is done, and an out_of_range exception is thrown
  if any of the passed in names are not found.
*/
SNeData SNeData::operator[](const vector<string>& namearr) const {

  int nelements = namearr.size();

  if (nelements == 0) return SNeData(Name);  

  SNeData retval(Name);
  
  retval.points.resize( nelements );

  //The indexing is a lot more complex in this case, especially
  // because we want to retain the order requested

  map<string,int> name_index_map;
  int ncurrent = points.size();
  for (int i = 0; i < ncurrent; ++i)
    name_index_map[ points[i].name ] = i;

  map<string,int>::iterator endpoint = name_index_map.end();
  map<string,int>::iterator location;
  
  for (int i = 0, j = 0; i < nelements; ++i) {
    location = name_index_map.find( namearr[i] );
    if (location == endpoint) {
      //Not found
      stringstream s("");
      s << "Name " << namearr[i] << " not found in " << Name;
      throw out_of_range(s.str());
    }
    retval.points[j++] = points[ location->second ];

  }
  return retval;
}

/*!
  \param indxarr List of indicies to remove 

  Duplicate indicies are allowed, but no range checking is performed,
  so it is up to the caller to ensure that only valid indicies are passed.
*/
void SNeData::remove( const vector<int>& indxarr ) {

  int nelements = points.size();
  int nremove = indxarr.size();

  if (nremove == 0) return;

  vector<SNeDataEntry> newpoints(nelements-nremove);

  //Make a vector to hold which elements to include
  //A bitset can't be used because we don't know the number of
  // elements ahead of time
  vector<bool> included(nelements, true);

  for (int i =0; i < nremove; ++i)
    included[ indxarr[i] ] = false;
  
  for (int i = 0, j = 0; j < nelements; ++j)
    if ( included[j] ) newpoints[i++] = points[j];
  
  points = newpoints;
}

/*!
  \param namearr List of names to remove
  \param strict Complain if passed in name isn't found in object

  Duplicate indicies are allowed.  Invalid entries in namearr
  result in a out_of_range error being thrown if strict is true
*/
void SNeData::remove( const vector<string>& namearr, 
		      bool strict) {

  int nelements = points.size();
  int nremove = namearr.size();

  if (nremove == 0) return; //Nothing to remove

  //The indexing is a lot more complex in this case, especially
  // because we want to retain the order requested, and because
  // we allow duplicates

  map<string,int> name_index_map;
  for (int i = 0; i < nelements; ++i)
    name_index_map[ points[i].name ] = i;

  map<string,int>::iterator endpoint = name_index_map.end();
  map<string,int>::iterator location;

  //Make a vector to hold which elements to include
  //A bitset can't be used because we don't know the number of
  // elements ahead of time
  vector<bool> included(nelements, true);

  for (vector<string>::const_iterator i = namearr.begin();
       i != namearr.end(); ++i) {
    location = name_index_map.find( *i );
    if (location == endpoint) {
      //Not found
      if (strict) {
	stringstream s("");
	s << "Name " << *i << " not found in " << Name;
	throw out_of_range(s.str());
      }
      //Else ignore this entry
    } else included[ location->second ] = false;
  }

  int nremaining = count( included.begin(), included.end(), true );

  vector<SNeDataEntry> newpoints(nremaining);

  for (int i = 0, j = 0; j < nelements; ++j)
    if ( included[j] ) newpoints[i++] = points[j];

  points = newpoints;
}


/*!
  \param indxarr List of indicies to remove from returned copy

  Duplicate indicies are allowed, but no range checking is performed,
  so it is up to the caller to ensure that only valid indicies are passed.
*/
SNeData SNeData::copy_remove( const vector<int>& indxarr ) const {

  SNeData retval(Name);

  int nelements = points.size();
  int nremove = indxarr.size();

  if (nremove == 0) return retval = *this; //Nothing removed

  retval.points.resize(nelements-nremove);

  //Make a vector to hold which elements to include
  //A bitset can't be used because we don't know the number of
  // elements ahead of time
  vector<bool> included(nelements, true);

  for (int i =0; i < nremove; ++i)
    included[ indxarr[i] ] = false;
  
  for (int i = 0, j = 0; j < nelements; ++j)
    if ( included[j] ) retval.points[i++] = points[j];

  return retval;
}

/*!
  \param namearr List of names to remove from returned copy
  \param strict Complain if passed in name is not found in object
  \returns Copy of SNeData with specified entries removed

  Duplicate indicies are allowed.  Invalid entries in namearr
  result in a out_of_range error being thrown unless strict is false.
*/
SNeData SNeData::copy_remove( const vector<string>& namearr,
			      bool strict ) const {

  SNeData retval(Name);

  int nelements = points.size();
  int nremove = namearr.size();

  if (nremove == 0) return retval = *this; //Nothing removed

  //The indexing is a lot more complex in this case, especially
  // because we want to retain the order requested, and because
  // we allow duplicates

  map<string,int> name_index_map;
  for (int i = 0; i < nelements; ++i)
    name_index_map[ points[i].name ] = i;

  map<string,int>::iterator endpoint = name_index_map.end();
  map<string,int>::iterator location;

  //Make a vector to hold which elements to include
  //A bitset can't be used because we don't know the number of
  // elements ahead of time
  vector<bool> included(nelements, true);

  for (vector<string>::const_iterator i = namearr.begin();
       i != namearr.end(); ++i) {
    location = name_index_map.find( *i );
    if (location == endpoint) {
      //Not found
      if (strict) {
	stringstream s("");
	s << "Name " << *i << " not found in " << Name;
	throw out_of_range(s.str());
      }
      //Otherwise ignore this entry
    } else included[ location->second ] = false;
  }

  int nremaining = count( included.begin(), included.end(), true );

  retval.points.resize(nremaining);

  for (int i = 0, j = 0; j < nelements; ++j)
    if ( included[j] ) retval.points[i++] = points[j];

  return retval;
}
