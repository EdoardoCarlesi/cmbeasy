#ifndef SPLINE
#define SPLINE


#include "miscmath.h"
#include "cleanvector.h"
#include "anchor.h"

#include <map>
#include <list>
#include <string>
#include <vector>
#include <fstream>
#include <cmath>


#ifdef SPLINEDBG
#define SPLINE_DBG_MSG(m) cerr << "[SPLINE DEBUG] Spline::" << m << endl;
#else
#define SPLINE_DBG_MSG(m) 
#endif

#ifdef _OPENMP
#include "omp.h"
#endif

/*! 
  Class for cubic spline interpolation of data. 
  Supports shared x-data among splines that have
  this data in common, fast access to the
  right interpolation interval and caching of
  precalculated data.

  It is fast, efficient and stable. 
  \par Constructors
  There are several constructors:
  \code
  Spline(const int = 1000, string = "unnamed",Anchor* =0);
  Spline(Spline *, string = "unnamed",Anchor* = 0);
  Spline(const Spline&);
  Spline(ifstream &i, string = "unnamed",Anchor* =0); 
  Spline(moSingle, const Mathobject&, Spline*, string name, Anchor* =0);
  \endcode
  Common to all is that you can specify a name (which you
  should for debugging is easier with this) and an Anchor.
  Splines are AnchorEnabled and hence the destruction
  of the Anchor will destroy the Spline, if you give an
  Anchor. This makes it very easy to prevent memory
  leaks.
  The first constructor creates a spline with default 1000 data
  points. If you fill more than 1000 points, this is no problem,
  as Splines will double their size each time you set more
  data points than allocated. On calling arm(), the overhead
  memory is freed.  
  The second constructor is for constructing a child
  Spline from a mother given as the first argument.
  Children share the same x-data with their mother
  spline, in fact, the x-data is only allocated once for
  the mother and never again. 
  Using child splines provides very effecient
  caching possibilities: Say you have two splines,
  mother m and child c and you usually access both
  in the same routine at the same x-value. If you 
  then use m(x) and subsequently c(x), the c(x)
  call will notice that already several variables
  have been calculated for the mother and these
  will be used for calculating c(x). 
  The fourth constructor reads data from file. Use
  save() to write this file beforehand.
  The last constructor is for constructing child splines
  that given a function moSingle func which you
  specify returns a spline with the same x-data as
  mother (naturally) and the y-data set to func(x-data)

  \par Basic usage
  Using Splines is deceivingly simple. After you have
  created a spline, you set the data points using 
  the set(x,y) or if the spline is child spline the
  set(y) function. If you just don't bother, use
  setForce(x,y), this will decide for you.
  Setting is always in some order, either decending
  or ascending. If you set in descending order, you
  will have to call flip() to flip the data in ascending
  order, before you arm() the spline, cause 
  spline x-data comes in increasing values.
  It is not allowed to arm() twice (see below, why
  you should want to do this). Call disarm() first.
  Then you may call arm() again to get the proper
  Spline-interpolation table.
  So as an example, we create a sine - spline
  \code
  Spline example(1000,"example");
  for (double x = -3; x < 5 ; x += 0.1) 
  example.set(x, sin(x));
  \endcode
  In order to access the interpolation capabilities,
  you now have to call arm(). This builds up
  the table of second derivatives. If you try
  to access and the spline is not arm()ed, you
  will be thrown an Bad_Error().
  \code
  example.arm();
  \endcode
  This is all. If you now would like to know the
  value at x=0.12345, you call
  \code
  double y = example(0.12345);
  \endcode
  As a synonym, you can call fastY() this is
  exactly the same and only provided for 
  cases where you have a pointer to a spline.
  In this case you have two possibilities:
  either use the operator () which sometimes looks
  akward or use fastY. We also use an anchor here:
  \code 
  Anchor anchor;
  Spline *pointerSpline = new Spline(1000,"pspline",&anchor);
  for (double x = -3; x < 5 ; x += 0.1) 
  pointerSpline->set(x, sin(x));
  pointerSpline->arm();
  double y  = pointerSpline->fastY(0.12345); // 1st possib.
  double y2= (*pointerSpline)(0.12345);        // 2nd possib
  \endcode
  


  \par Arithmetics
  You can multiply or add a usual number to the 
  whole spline y-data, say we already have a spline example
  with values set, then we could
  \code 
  example *= 0.3;    // multiply all y-data by 0.3
  example += 5;      // and add a constant 5 to it
  \endcode
  If you would like to act differently on different 
  data points, you can call 
  \code
  example.mul(7,0.3);
  example.mul(8,0.5);
  \endcode
  This will multiply data point #7 by 0.3 and #8 by 0.5.
  If for instance you would like to multiply the 
  spline by a function sin(x*x), you could do:
  \code
  for (int i = 0; i < example.size(); i++) 
  example.mul(i, sin(example.x(i)*example.x(i));
  \endcode
  The size() just returns the number of data points, 
  x(i) returns x-data points number #i and hence
  you multiply by sin(x*x).
  After doing arithmetics, you will have to call
  arm() again. Keep in mind that if you have
  already called arm(), you will have to 
  disarm() first.

  \par Differentiation and Integration
  You may differentiate a given spline, say the
  spline holding sine data
  \code
  Spline Sine, Cosine;  // we dont care for names here
  .... fill Sine with sine values
  Sine.arm();
  Sine.derive(&Cosine);
  Cosine.dump("cosine");  // write y-data to file "cosine.dat"
  Sine.dump("sine");        
  \endcode
  So both splines data be written (human readably) to the files
  sine.dat and cosine.dat and you may
  take a look at them in e.g. gnuplot.

  If you like to know the integral of the Sine spline
  above ,call 
  \code
  double y  = Sine.integrate();          // whole x -range
  double y2  = Sine.integrate(-1,1);  // x = -1... 1 
  double y3 = Sine.average();     // integrate and divide by x range
  \endcode
  Notice the average() call above, returning the integral divided by the integration interval 
(in this case, no interval is given, which defaults to the  
  entire Spline).
  
  



*/


class  Spline : public Mathobject, public AnchorEnabled
{
public:
  enum childrenToo { no, all, thoseReady,yes};  //!< as a matter of fact, this just gives some "words" you can use
  string name;   //!< the name

  double *xdat, *ydat, *sdat;
  bool armed; //! if the spline data is calculated

private:
  int n;
  int maxDim;
  Spline *mother;   //!< if there is a mother, its this
  bool own;  //!< if the xdat belongs to the spline, i.e. mother ==0. Own does not say, that xdat is a valid pointer, it may be 0, if for instance, mother is killed, own will be true, but the pointer will be set to zero
  bool valid;  //!< if the xdat pointer is valid, i.e. own or mother alive And valid
  bool childExists; //!< if children is not empty
  int nvc,vc;
  int acc;  //!< Counter for operator() this is for perfomance monitoring
   int splintcC;
  static double splg[20010];        //!< splini and splder need it...
  bool splgValid;  //!< derive() needs the splg[] array, which is initialized if it is not valid

  int lastarm;  //!< keep track of last arm-position for rearm(), if you want to use this feature...

  static unsigned int nvar;  //!< The size of the convolution vector 
  static double zeroLevel;

  map <Spline*, Spline*> children;  //!< a map of children, if there are any


  struct CachedValues {
    double lastX, h, h26, a, b;
    double a3a,b3b;
    unsigned int klo, khi;
    unsigned int guess;
    bool valid;
  }; //!< caching these for children's use

  struct Cache {
      Cache() {
#ifdef _OPENMP
          threadCount = omp_get_max_threads();
#else
          threadCount = 1;
#endif
          mData = new CachedValues[threadCount];
          for(int i=0; i<threadCount; ++i) {
             mData[i].valid=false;
             mData[i].guess=0;
          }
      }
      ~Cache() { delete[] mData; }

      void initGuess() {
          for(int i=0; i<threadCount; ++i) { mData[i].guess=0; }
      }
      void invalidate() {
          for(int i=0; i<threadCount; ++i) { mData[i].valid=false; }
      }
    private:
      Cache(const Cache&);
      CachedValues *mData;
      friend class Spline;
      int threadCount;
  } mCache;

  inline CachedValues& cachedVals() {
#ifdef _OPENMP
      return mCache.mData[omp_get_thread_num()];
#else
      return mCache.mData[0];
#endif
  }

  int cC;

  int generated;

  double gen; // see function generatedAt() below
  double convResult; //!< the result of a family convolution, if requested


  static Spline *partner;   //!< convolutions partner spline
  static double shift;   //!< a shift in integration
  static double faktor;  //!< and a faktor for x in convolution
  static double *convXdat; //!< for vectorConvolution this is the pointer to the actual x-data. Either taken fom the first spline in convolutionVektor or its mother
  static Spline* convTwin; //!< for the specialised twin convolution, the second twin
  


 

  double checkSm, checkSm2;
  //int cC;
  void motherKilled();  //!< if mother is killed
  void addChild(Spline *c) {
#pragma omp critical (SplineChildren)
    {
      children[c] = c; childExists = true;
    }
  }  //!< adds another Child to the Child-map
public:
  void removeChild(Spline *c) {
#pragma omp critical (SplineChildren)
    {
      children.erase(children.find(c));
      childExists = !children.empty();
    }
  }  //!< removes a child from the childmap



  void checkSpace(int=-2);  //<! Checks whether there is space for one more data-point. If not, allocates more

  double ck(double x) { 
    if (x == 0) return 1.0; 
    return log(fabs(x));
  }

  double conv(const double);  //!< romberg-convolution integrand
  void convD(const double,const double*, double*); //!< usual convolution integrand
  // void convVektor(const double,const double*, double*); //!< vektor convoltuion integrand
  void convVektor2(const double,const double*, double*); //!< new version 
  // void convVektorTwin(const double,const double*, double*); //!< specialised to two splines 
  // void convVektor_backup(const double,const double*, double*); //!< vektor convoltuion integrand
  //void convFamily(const double,const  double*, double*); 
  double convOneDim(const double);

  double inte(const double);  //!< romberg-integration integrand
  void inteD(const double,const double*, double*);  //!< odeint integration integrand
  double inteOneDim(const double);

  childrenToo killChildrenWhenDying;
  
  
 public:
  static int WatchCount; //!< If spline is compiled with ENABLE_SPLINE_WATCH, counts number of splines alive
  static map<string,int>  WatchMap; //!< If compiled with ENABLE_SPLINE_WATCH, tracks number of splines with that name
  static void printWatchMap(); //!< output how many splines with what name are currently alive
  static vector<Spline*>  convolutionVektor;
  void dump(string,bool);
  void dump(string,childrenToo InterpolatedAlso=no);
  void dump();
  void dump(ofstream &, bool = true);
  void save(ofstream&); //!< save spline x and y data to file 
  void save(const string&); //!< save spline x and y data to a file
  void save(); //!< save spline x and y data to  file; name of the file is "spline name" + ".spline"
  Spline(const int = 1000, string = "unnamed",Anchor* =0);

  explicit Spline(Spline *, string = "unnamed",Anchor* = 0);
  explicit Spline(const Spline&, string, Anchor* =0); //!< Construct a Spline that has its own data, is un-armed, has no children and otherwise gets its x and y data from the Spline in the argument.
  explicit Spline(const Spline&); //!< Create a Spline that is an exact copy of the Spline in the argument.
  Spline(ifstream &i, string = "unnamed",Anchor* =0); //!< read data from file (i.e re-create a spline from the save() data )
  Spline(moSingle, const Mathobject&, Spline*, string name, Anchor* =0); //!< use function moSingle to create a spline at x-data of another spline
  ~Spline();

  //! create a spline from text file (e.g. re-create a spline from a file written by dump())
  static Spline* fromTextFile(const string& fileName, string="unnamed", Anchor* =0);
  static Spline* fromTextTable(ifstream &i, string="unnamed", Anchor* =0); //!< create a spline from a text x-y table (e.g. re-create a spline from a file written by dump())

  void error(string,string,int = 1234) const;  //!< shows this spline, mother and a lot of information
  bool isArmed() const { return armed; }

  
  void createChecksum();
  void checksum();
  Spline* getChild(string); //!< return child with name s. If there are more than a few, this is costy !!!
  void setName(string s) { name = s; } //!< set's the splines name
  int first() const { return 0;}  //!< index of first data-point
  int last() const { return n; } ; //!<  index of last data-point

  int indexToRight(const double x) const; //!< return the smallest i with x(i) >= x
  int indexToLeft(const double x) const; //!< return the largest i with x(i) <= x

  int size() const { return n+1;}  //!< the size i.e. n+1

  void arm(childrenToo = no, bool keepSize = false, int=0);

  //void arm(bool childrenToo = false, bool keepSize = false, int =1);  //!< creates the spline-interpolation table
  void rearm(bool keepSize = true);

  //! verifies that x-data is correctly ordered and does not contain
  //! duplicates; throws Bad_Error if that is not the case
  void proper();
  void makeProper(); //!< reshuffles the spline data so that x-data is monotonic, removes duplicates

  void disarm(childrenToo = no, bool keepSize=false); //!< deletes the spline-interpolation table 
  void set(const double, const double); //!< new data point (x,y) for a spline that owns its xdata
  void set(const double); //!< new data point (y) for a spline that shares the xdata with mother and sisters
  /*! Wrapper for set() that depending on wether
    the spline owns the xdata or not, calls the
    two possible set() functions */
  void setForce(double x, double y ) {  
    if (own) set(x,y); else set(y); 
  }

  //  template<class T> void set(const map<T,T>& m);  
  void set(const map<float,float>& m);
  void set(const map<double,double>& m);

  void setKillChildrenWhenDying(childrenToo c) { killChildrenWhenDying = c; }

  void getData(const Spline&,childrenToo=no); //!< Copy x (if necessary) and y data from another spline, set n. If spline does not have own x-data, the copy uses mothers-xdata 


  void resize(int);  //!< allocates arrays for x and y and copies existing data to these
  void add(int,double);  //!< add x to ydat[k] 
  void mul(int k, double x); //!< multiply ydat[k] by x
  void div(int k, double x) { mul(k,1.0/x);} //!< call mul(1/x)
  void setX(int,double);  //!< set xdat[index] to x
  void setX(double x) { setX(n+1,x); } //!< sets next X
  void setY(int,double); //!< sets ydat[index] to y 
  double *ptrY(int);   //!< returns pointer to element i of ydat. Also adjusts n. 
  double *ptrY() { return ptrY( n+1); } 

  void flip(); //!< flip the x axis, i.e. make first last and so on

  void setN(int k=-2); //!< Set n of this spline.
  void setChildrensN(int k =-2); //!< Call setN() for all children. All precautions of setN are true here, too

  double y(int) const;   //!< return ydat[k]
  //! return xdat[k]
  double x(int k) const {   
    if (k > n || k < 0 ) { 
      cout << "asking for " << k << endl;
      error("x","Element not available");
    }
    if (own) return xdat[k]; else return mother->xdat[k];
  };

  double start() const { return x(0);};   //! Return first x
  double stop() const { return x(n);}; //! Return last x 
  double start(int k) const { return x(k); };
  double stop(int k) const { return x(n-k); };

  double front() const { return y(0);};
  double back() const { return y(n);}; 
  double front(int k) const { return y(k); };
  double back(int k) const { return y(n-k); };

  double range(int  step=1) { return (x(n) - x(0))/step; } //!< return (last-first)/step (useful for for-loops)

  double operator() (const double);  //!< Fast and caching access to the interpolation value. Should usesually be the one to call (not splint() or fastSplint())

  double operator[] (const int n) { return y(n); }

  double fastY(const double x) { return (*this)(x); } //!< stub to operator()
  double safeY(const double x);//!< used by derive to keep within spline boundaries

  bool isWithinBounds(const double x); //!< true, if x value is within the boundaries of spline data

  double randomY(double);  //!< return interpolation at point x. use this instead of operator (), if you know that the access is random 

  double inBetweenY(double);

  //! Scalar multiply all ydat[] by const double x
  const Spline& operator *=(const double x) {
    for (int i=0; i<= n; i++) ydat[i] *=x;
    return *this;
  }

  //! add x to all ydat[]
  const Spline& operator +=(const double x) {
    for (int i=0; i<= n; i++) ydat[i] +=x;
    return *this;
  }


  // DERIVATIVES, INTEGRALS AND CONVOLUTIONS:
  void splder(Spline &s);
  void derive(Spline &s,childrenToo=no); // derive and if specified anything else but none, arm this spline

  double convolution(Spline*,double,double,double = 1e-4,double=1.0,double=1.0);
  // double rombConvolution(Spline*,double,double,double =1e-4,double=1.0,double=1.0);
  // void  familyConvolution(Spline*,double,double,double = 1e-4,double=1.0,double=1.0);
  // void  vectorConvolution(Spline*,double,double,double = 1e-4,double=1.0,double=1.0);
  static void  vectorConvolution2(Spline*,double,double,double = 1e-4,double=1.0,double=1.0,double hnext=0);


  void checkConvolutionRange(Spline *p, double *, double*, const double, const double);

  void dumpConvolution(Spline*,double,double,double,double,const char*);

  double integrate(double a, double b,double tol = 1e-6);
  double integrate(double tol = 1e-6) { return integrate( min(x(0), x(n)) , max(x(0),x(n)), tol); }

  double average(double a, double b, double to = 1e-6);
  double average(double b) { return average(start(), b); }
  double average() { return average(stop()); }


  double stepIntegrate(double a, double b);
  double stepIntegrate() { return stepIntegrate( min(x(0), x(n)) , max(x(0),x(n))); }
  double stepInt(int a, int b);
  double findZero(double a, double b, double =1e-5);
  //! return x value for which f(x)=y between x=a and x=b (proxy for Miscmath::zbrent(), same
  //! restrictions apply
  //! This is essentially the same as findZeroBetween(), but less clumsy (e.g. no class-static
  //! zeroLevel)
  double findXForY(const double y, const double a, const double b, const double tol=1e-5);
  double findZeroBetween(const int a, const int b, double =1e-5, const double zl=0.0);
  list<double>* getZeroList(const double zl=0.0, const double tol = 1e-2, int a=0,int b=0);
  vector<double> getZeroVector(const double zl=0.0, const double tol = 1e-2, int a=0,int b=0); //!< convenience wrapper for getZeroList()
  int  getZeroArray(int a, const  int b, double *, const int,const double tol = 1e-5, double zl=0.0); 
  void getMaxima(list<double>*, int a=0,int b=0,const double tol = 1e-5);
  void getMinima(list<double>*, int a=0,int b=0,const double tol = 1e-5);
  void getExtrema(list<double>*, list<double>*, int a=0,int b=0,const double tol = 1e-5);
  double maximum(); //!< return x  of the largest maximum in this spline *including* the boundaries, if they are higher
  double minimum(); //!< return x  of the smallest minimum in this spline *including* the boundaries, if they are below
  
  void merge(const double*, const int, const bool = true,const double* =0,const double keepSize = true);
  void merge(Spline &s, const bool = true);

  void smoothen(int k=1, int m=2); //!< smoothen the spline within k neigbouring points to the left and k to the right and 2*m*k + 1 points in this range


  void  printStatus() const;
  void  explizit();
  void setMother(Spline&);
 
  static void generateFamily(vector<Spline*> &, const unsigned int , const int , const string);
  static void generateSisters(vector<Spline*> &, const unsigned int , const int , const string);
  static void generateChildren(vector<Spline*> &, Spline*, const unsigned int , const string);
  static void generateOne2OneChild(vector<Spline*> &, vector<Spline*> & , const string);

  static void armVector(vector<Spline*>&);
  static void disarmVector(vector<Spline*>&);

  static void kill(Spline *);
  static  vector<Spline*>* expandVectorBetween(vector <Spline*> &, vector<Spline*> &,int between); //!< Generate "between" splines between each existing pair of spines in v
  static  vector<Spline*>* expandVector(vector <Spline*> &, vector<Spline*> &,int total); //!< Equally -spaced expansion of the spline vector
  static void  generateSpline(Spline* e, vector<Spline *> &v, vector<Spline *> &p, const double k,bool=false);
  
  void  generateSplineXXL(Spline* e, vector<Spline *> &v, vector<Spline *> &p, const double k);
  static void generateAndDump(vector<Spline *> &v, vector<Spline *> &p, string,const double k);

  double generatedAt() const ;
  double convolutionResult() { return convResult; }

  void splini(); //! initialize whatever...
  void Splder(double*, double*, const int); //! Fits cubic splines to y and interpolates first derivs at grid points dy

  double fitNeeds(const pair<double,double>);

  void fitToSingle(const pair<double,double>);

  void origSpline(double x[], double y[], int, int n, double yp1, double ypn, double y2[]);
  double splint(double xa[], double x, unsigned int *guess); //! Spline interpolation Nrecipes + guess for next

   //! Spline interpolation, modified version that uses sequential search starting from a guess, instead of interval halfing
  double fastSplint(double xa[],  double x, unsigned int *guess);


#ifdef _OPENMP
 private:
  omp_lock_t mElementAccessLock;
  omp_lock_t* lock() { omp_set_lock(&mElementAccessLock); return &mElementAccessLock; }
#endif
};


inline double Spline::y(int k) const {
    if (k > n || k < 0 ) error("y","Element not available",k); 
    return ydat[k];
}

inline double* Spline::ptrY(int k) {
  if (k >= maxDim) error("ptrY","Element not in range",k); 
  if (k > n) n = k;
  return &ydat[k];
}





/*! 
  The same as the original splint, except that the original version does a half
  step strategy in order to find the right interval and this one 
  expects a guess, and searches sequentially from that guess on for a 
  matching interval. This is of course much faster for a 
  nearby solution to the problem. 
  
  However, if after a few sequential searches, it didn't succeed, it calles
  splint() for its half-step strategy

  In additon, it caches some quantities for faster access by
  related splines.

  The guess is set to the new guess.
*/

inline  double Spline::fastSplint(double xa[], double x, unsigned int *guess) {
  CachedValues& cache = cachedVals();
  if (x >= xa[n]) {
#ifndef FORCE_SPLINE_SPEED
    if (x > 3*xa[n] - 2*xa[n-1]) { cout << x << "  bound  " << xa[n] << endl; error("fastSplint","x larger than upper bound"); }
#endif
    *guess = n-1; 
    cache.valid = false;
    return ydat[n]; 
  }
  if (x <= xa[0]) { 
#ifndef FORCE_SPLINE_SPEED
    if (x < 3*xa[0] - 2*xa[1])    { cout << x << "  bound  " << xa[0] << endl;error("fastSplint","x smaller than lower bound"); }
#endif
    *guess = 0; 
    cache.valid = false;
    return ydat[0]; 
  }


  unsigned int& klo = cache.klo;
  klo = *guess;
  unsigned int& khi = cache.khi;
  khi = klo+1;
 
  //cout << "++ " << klo << "  " << khi << endl;
#ifndef  FORCE_SPLINE_SPEED
  if (klo <0 || khi > n) error("fastsplint","guess is out of range");
#endif

  if (xa[klo] <= x) { 
    while (klo < n) {
      if (xa[khi] > x) break;
      klo++; khi++;
      if (klo - *guess > 5) return splint(xa, x, guess); //! if too many steps, do half step
    } 
  } else {
    while (klo >0) {
      if (xa[klo] < x) break;
      klo--; khi--;
      if (*guess - klo > 5) return splint(xa,x,guess); //! if too many steps, do half step
    }
  }

  //if (name == "besselfunction") cout << "abstand: " << (klo-*guess) << endl;
  *guess = klo;
  // if there are children, then we fill the cache, if not, we use faster local variables
  if (own && childExists) {
    double &h = cache.h;
    double &a = cache.a;
    double &b = cache.b;
    double &a3a = cache.a3a;
    double &b3b = cache.b3b;
    double &h26 = cache.h26;
    double &lastX = cache.lastX;
    h = xa[khi]-xa[klo];
    if (h == 0.0) throw Bad_Error("Spline::fastSplint() Bad xa input to routine splint");

    double invh=1.0/h;
    a=(xa[khi]-x)*invh;
    b=(x-xa[klo])*invh;
    a3a = a*a*a -a;
    b3b = b*b*b -b;
    h26=h*h*0.1666666666666;
    lastX = x;
    cache.valid = true;
    return a*ydat[klo]+b*ydat[khi]+(a3a*sdat[klo]+b3b*sdat[khi])*h26;
  } else {
    double h = xa[khi]-xa[klo];
    if (h == 0.0) throw Bad_Error("Spline::fastSplint() Bad xa input to routine splint");
    double invh=1.0/h;
    double a=(xa[khi]-x)*invh;
    double b=(x-xa[klo])*invh;
    return a*ydat[klo]+b*ydat[khi]+((a*a*a-a)*sdat[klo]+(b*b*b-b)*sdat[khi])*(h*h)*0.1666666666666;
  }
}

struct noConvolutionOverlap {}; //!< thrown, if the convolution would have no overlap, i.e. result is 0


/*!
  Very convenient class for handling the common
  problem of having a 2-D grid of data points which
  you would like to slice through and print out to
  some file. 

  In cosmos and cmbcalc wildly used to store and
  print transfer functions.

  Creating a SplineWeb is particularily easy, you
  just have to provide the size of the grid in x and y
  directions. 
  \code
  SplineWeb myWeb(10,20);
  \endcode
  Now you can set data points by
  \code
  myWeb.set(3,5,1.234);
  \endcode
  This will create 10 Splines pointing in Y-Direction
  and 20 Splines pointing in X-Direction.

  Of course, you then are allowed to set 10 data points
  in Y-direction and 20 Datapoints in X-direction. As
  the splines are the limiting factor.

  For convenience, you just tell the set() routine the
  real world number, not the index of the spline.
  
  Make sure that these numbers are always the same
  for each row and collumn:

  In other words, for our example, make sure that 
  you do only give 10 different X-Values and only 
  20 different Y-values.

  But as you usually will be filling the grid by a 
  nested loop, no problem here in practice.
 
*/

class SplineWeb : public AnchorEnabled {
  CleanVector<Spline*> X;
  CleanVector<Spline*> Y;

 private:
  pair<double,double> maximum_internal(double a, double b); //! for internal maximum use

 public:
  string name;
  map<double,unsigned int> xgrid, ygrid;
  int igrid, jgrid;  // count the number of points in X and Y direction

  SplineWeb(string n="unknown",Anchor *a=0,const int x=3, const int y=3);
  ~SplineWeb();

  void printStatus(std::ostream& os=std::cout);

  //! Scalar multiply all ydat[] by const double x
  // After this, you will have to arm the Web again !!!!
  const SplineWeb& operator *=(const double x) {
    for (unsigned int i=0; i< X.size(); i++) (*X[i]) *= x;
    for (unsigned int i=0; i< Y.size(); i++) (*Y[i]) *= x;
    
    return *this;
  }

  void shrinkToFit(); //!< resize to fit igrid and jgrid, in general, make it smaller to save memory
  void makeProper(); //!< calls makeProper() for each internal Spline.

  void arm(); //!< arm all vectors in the web. If already armed, this will disarm() first..
  void disarm(); //!< disarm the web, if it armed
  void set(const double x, const double y, const double z);
  void dumpAlongX(const double y, string);
  void dumpAlongY(const double x, string);
  Spline* createAlongX(const double y, string,Anchor *a=0, Spline *mother=0);
  Spline* createAlongY(const double x, string,Anchor *a=0, Spline *mother=0);
  void fitToSingleAlongX(const double y, const pair<double,double> c);
  void fitToSingleAlongY(const double x, const pair<double,double> c);
  //void resize(const int x, const int y);

  int xSize() const { return igrid; }  //!< return number of points along x-direction
  int ySize() const { return jgrid;}   //!< return number of points along y-direction

  //! return the i-th Spline of points in Y direction (i.e. at x-index i)
  Spline* xSplineAt(const int i) const { if (i > ySize()) return 0; return X[i]; }
  //! return the i-th Spline of points in X direction (i.e. at y-index i)
  Spline* ySplineAt(const int i) const { if (i > xSize()) return 0; return Y[i]; }


  pair<double,double> maximum(); //!< find the coordinates of the maximum of the web


  /*
#define WEB_RASTER_MAGIC 1.234009754e-32
  vector< vector<double> > * rasterize(int nx, int ny, double& totalInt, double xmin=WEB_RASTER_MAGIC, double xmax=WEB_RASTER_MAGIC, double ymin=WEB_RASTER_MAGIC, double ymax=WEB_RASTER_MAGIC); //!< put a raster on top of the spline web, create splines along the raster, evaluate and return a 2D vector of the function values
  
  list < pair<int,int> > * fractionOfRaster(vector< vector<double> > *, double frac); //!< return a list of indices (i,j) such that the summed up value of these raster patches gives frac of the totalInt of the whole  raster 
  */


  double x(int i); //!< return i-th x value, *ordered*, not in order of setting
  double y(int i); //!< return i-th y value, *ordered*, not in order of setting
};

#if 0
// short class to integrate an entire family of splines
class SplineFamilyIntegrator: public Mathobject
{
  public:
    SplineFamilyIntegrator(Spline* mother, double tol=1e-6);

    double integrate();
    double result(const Spline *const s) const;
    void integrationHelper(const double x, const double* y, double* dy); //!< internal

  private:
    const double mTolerance;
    Spline *mMother;
    std::vector<Spline*> mSplines;
    std::vector<double> mResult;
};
#endif

#endif
