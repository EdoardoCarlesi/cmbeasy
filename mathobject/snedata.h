#ifndef snedata_h
#define snedata_h

#include <vector>
#include <string>
using namespace std;

/*!
  \brief Data structure specialized for holding supernova data.
*/
struct SNeDataEntry {
  string name; //!< Name of supernova
  double zcmb; //!< Redshift of supernova in CMB frame
  double zhel; //!< Heliocentric redshift of supernova
  double dz; //!< Error in redshift
  double widthpar; //!< Width-luminosity parameter (stretch, delta m_15, etc.)
  double dwidthpar; //!< Error in widthpar
  double mag; //!< Magnitude
  double dmag; //!< Error in magnitude
  double colourpar; //!< Colour related parameter ( E(B-V), A_V, etc. )
  double dcolourpar; //!< Error in colourpar
  double cov_mag_widthpar; //!< Covariance between magnitude and widthpar
  double cov_mag_colourpar; //!< Covariance between magnitude and colourpar
  double cov_widthpar_colourpar; //!< Covariance between widthpar and colourpar

  /*! \brief Constructor */
  SNeDataEntry(string SNNAME="",double ZCMB=0,double ZHEL=0,double DZ=0,
	       double WIDTHPAR=0, double DWIDTHPAR=0,double MAG=0,double DMAG=0, 
	       double COLOURPAR=0, double DCOLOURPAR=0, double COV_MAG_WIDTHPAR=0,
	       double COV_MAG_COLOURPAR=0, double COV_WIDTHPAR_COLOURPAR=0) :
    name(SNNAME),zcmb(ZCMB),zhel(ZHEL),dz(DZ),widthpar(WIDTHPAR),
       dwidthpar(DWIDTHPAR),mag(MAG),dmag(DMAG),colourpar(COLOURPAR),
       dcolourpar(DCOLOURPAR),cov_mag_widthpar(COV_MAG_WIDTHPAR),
       cov_mag_colourpar(COV_MAG_COLOURPAR),
       cov_widthpar_colourpar(COV_WIDTHPAR_COLOURPAR) {}; 

  /*! \brief Copy Constructor */
  SNeDataEntry(const SNeDataEntry& inval) {
    name = inval.name; zcmb = inval.zcmb; zhel = inval.zhel; dz = inval.dz;
    widthpar = inval.widthpar; dwidthpar = inval.dwidthpar; mag = inval.mag;
    dmag = inval.dmag; colourpar = inval.colourpar; 
    dcolourpar = inval.dcolourpar; cov_mag_widthpar = inval.cov_mag_widthpar;
    cov_mag_colourpar = inval.cov_mag_colourpar;
    cov_widthpar_colourpar = inval.cov_widthpar_colourpar;
  }

  /*! \brief Copy operator */
  SNeDataEntry& operator=(const SNeDataEntry& inval) {
    if (this == &inval) return *this; //Self copy protection
    name = inval.name; zcmb = inval.zcmb; zhel = inval.zhel; dz = inval.dz;
    widthpar = inval.widthpar; dwidthpar = inval.dwidthpar; mag = inval.mag;
    dmag = inval.dmag; colourpar = inval.colourpar; 
    dcolourpar = inval.dcolourpar; cov_mag_widthpar = inval.cov_mag_widthpar;
    cov_mag_colourpar = inval.cov_mag_colourpar;
    cov_widthpar_colourpar = inval.cov_widthpar_colourpar;
    return *this;
  }

  /*! \brief Magnitude corrected for stretch and colour */
  double getCorrMag(double alpha, double beta, double widthmean=1) const {
    return mag + alpha * ( widthpar - widthmean ) - 
      beta * colourpar;
  } 

  /*! \brief Retuns variance of corrected magnitude.
    Does not include peculiar velocity errors, since the exact
    effect of these depends (mildly) on the cosmological parameters.
    Also does not include any intrinsic 'dispersion'. */
  double getCorrVar(double alpha, double beta) const {
    return dmag*dmag + alpha*alpha * dwidthpar*dwidthpar + 
      beta*beta * dcolourpar*dcolourpar + 
      2 * alpha * cov_mag_widthpar - 2 * beta * cov_mag_colourpar
      - 2 * alpha * beta * cov_widthpar_colourpar;
  }

};


/*!
  \brief Sorting class for SNeDataEntry based on zcmb.
*/
class SNeSortByZcmb {
 public:
  /*! \brief Comparison operator on zcmb */
  int operator()(const SNeDataEntry& sn1, const SNeDataEntry& sn2) {
    return sn1.zcmb < sn2.zcmb;
  }
};


/*!
  \brief Collection of SNeDataEntry objects

  All data points are stored in the vector<SNeDataEntry> points, which you can
  access through the [] operator.  This class also provides serialization
  capability (i.e., it can write the data to a file and read it back).
*/
class SNeData {
 protected:
  string Name; //!< Just a handy name, potentially useful for debugging
  vector<SNeDataEntry> points; //!< The data as SNeDataEntry's

  void stringwords(const string&, vector<string>&) const; //!< Breaks input string into a vector of words

 public:

  typedef vector<SNeDataEntry>::iterator iterator;
  typedef vector<SNeDataEntry>::const_iterator const_iterator;
  typedef vector<SNeDataEntry>::reverse_iterator reverse_iterator;
  typedef vector<SNeDataEntry>::const_reverse_iterator const_reverse_iterator;

  SNeData(); //!< Default constructor
  explicit SNeData(const string& name); //!< No data read
  SNeData(const string& name, const string& fileName); //!< Create SNeData object and read in data from file 
  void readData(const string&, bool verbose=false); //!< Read in the data from file 
  void writeData(const string&) const; //!< Write the data to a file

  void setName(const string& name) { Name = name; } //!< Set name of object
  string getName() const { return Name; } //!< Return the name of the object

  //Iterators
  iterator begin() { return points.begin(); } //!< Returns read/write iterator pointing to first element
  const_iterator begin() const { return points.begin(); } //!< Returns read only iterator pointing to first element
  iterator end() { return points.end(); } //!< Returns read/write iterator pointing to one past the last element
  const_iterator end() const { return points.end(); } //!< Returns read only iterator pointing to one past the last element

  //Reverse iterators
  reverse_iterator rbegin() { return points.rbegin(); } //!< Returns read/write reverse iterator pointing to last element
  const_reverse_iterator rbegin() const { return points.rbegin(); } //!< Returns read only reverse iterator pointing to last element
  reverse_iterator rend() { return points.rend(); } //!< Returns read/write reverse iterator pointing to one before first element
  const_reverse_iterator rend() const { return points.rend(); } //!< Returns read only reverse iterator pointing to one before first element

  //Size/capacity stuff
  int capacity() const { return points.capacity(); } //!< Space allocated (but not necessarily filled) in points
  int size() const { return points.size(); } //!< Number of elements in points
  bool empty() const { return points.empty(); } //!< Is points empty?
  void resize(int newsize) { points.resize(newsize); } //!< Resize points
  void reserve(int newsize) { points.reserve(newsize); } //!< Allocate but don't initialize more space in points
  
  //sorting
  void zcmbsort(); //!< Sort supernovae by zcmb

  //Combining two lists
  SNeData& operator+=(const SNeData&); //!< Concatenate two lists.

  //indexing
  SNeData operator[](const vector<int>&) const; //!< Return new list indexed from old
  SNeData operator[](const vector<string>&) const; //!< Return new list indexed by supernova name
  SNeDataEntry& operator[](const int i) { return points[i]; } //!< Unchecked subscript
  const SNeDataEntry& operator[](const int i) const { return points[i]; } //!< Unchecked subscript
  SNeDataEntry& at(const int i) { return points.at(i); } //!<Checked subscript
  const SNeDataEntry& at(const int i) const { return points.at(i); } //!< Checked subscript

  //Removal
  void remove(const vector<int>&); //!< Removes specified elements in place
  void remove(const vector<string>&, bool strict=true); //!< Removes specified elements in place using names
  SNeData copy_remove(const vector<int>&) const; //!< Returns a new list with indexed entries removed
  SNeData copy_remove(const vector<string>&, bool strict=true) const; //!< Returns a new list with entries with names in argument removed
  
};

#endif
