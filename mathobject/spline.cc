#include "spline.h"
#include "global.h"


#include <iomanip>
#include <sstream>
#include <limits>
#include <algorithm>

// if you think that some memory - holes due to allocated splines may
// exist, uncomment the following line. With this, the WatchCount and
// WatchMap will keep track of each Spline created (and destroyed).
// By calling printWatchMap(), you can output this information to cout
//#define ENABLE_SPLINE_WATCH

/*! Constructor for an independent spline (owning its x-data). You may specify the initial size
  of the x and y data arrays as well as a name, which is useful if something goes wrong, cause
  error can tell you which spline actually crashed. 

  The initial size is not too important, as the spline will dynamically allocate more memory if needed,
  also, when arm()ed, the spline will free access space. 
*/

Spline::Spline(const int nvar, string nom,Anchor* a)
            : AnchorEnabled(a), n(-1), name(nom), maxDim(nvar),
              mother(0), sdat(0), own(true), valid(true), armed(false), childExists(false),
              nvc(0), vc(0), acc(0), splintcC(0), splgValid(false), lastarm(0),
              gen(std::numeric_limits<double>::infinity()),
              checkSm(std::numeric_limits<double>::infinity()),
              checkSm2(std::numeric_limits<double>::infinity()),
              killChildrenWhenDying(no)
{
#ifdef _OPENMP
  omp_init_lock(&mElementAccessLock);
#endif
  xdat = new double[nvar];
  ydat = new double[nvar];

  for (int i = 0; i < nvar; i++)   ydat[i] = 0;
  
  if (xdat == 0 || ydat == 0 ) throw Bad_Error("no memory");

#ifdef ENABLE_SPLINE_WATCH
  WatchCount++;
  WatchMap[name]++;
#endif
 
}

Spline::Spline(Spline *s, string nom,Anchor* a)
           : AnchorEnabled(a),  n(-1),  name(nom), maxDim(s->maxDim),
             mother(s), xdat(0), sdat(0), own(false), valid(s->valid), armed(false),
             childExists(false),  splgValid(false), lastarm(0),
             gen(std::numeric_limits<double>::infinity()),
             checkSm(std::numeric_limits<double>::infinity()),
             checkSm2(std::numeric_limits<double>::infinity()),
             killChildrenWhenDying(no)
{
#ifdef _OPENMP
  omp_init_lock(&mElementAccessLock);
#endif
    ydat = new double[s->maxDim];
    if (ydat == 0 ) error("Spline(Spline*)","no memory");    
    for (int i = 0; i < s->maxDim; i++)  ydat[i] = 0;

    if (s->own) s->addChild(this); 
    else {
      mother = s->mother;
      mother -> addChild(this);
    }
#ifdef ENABLE_SPLINE_WATCH
  WatchCount++;
  WatchMap[name]++;
#endif
}

Spline::Spline(const Spline& s)
           : AnchorEnabled(s.MyAnchor), n(s.n), name(s.name),
             maxDim(s.maxDim), mother(s.mother), xdat(s.xdat), sdat(s.sdat), own(s.own),
             valid(s.valid), armed(s.armed), childExists(s.childExists), splgValid(s.splgValid),
             lastarm(s.lastarm),
             gen(std::numeric_limits<double>::infinity()),
             checkSm(std::numeric_limits<double>::infinity()),
             checkSm2(std::numeric_limits<double>::infinity()),
             killChildrenWhenDying(s.killChildrenWhenDying)
{
#ifdef _OPENMP
  omp_init_lock(&mElementAccessLock);
#endif
  if (own) xdat = new double[maxDim];

  if (armed) sdat = new double[maxDim];
  ydat = new double[maxDim];
  if ( ydat == 0) error("Spline(const Spline&)","no memory");  
  getData(s,all); // get its data (even sdat if its there)
#ifdef ENABLE_SPLINE_WATCH
  WatchCount++;
  WatchMap[name]++;
#endif
}

/*!
  This constructor is in some ways superior to the Spline(const Spline&), because it only is a copy w.r.t to the data.
*/
Spline::Spline(const Spline& s, string nom, Anchor* a)
                 : AnchorEnabled(a), n(s.n), name(nom), maxDim(s.n+1),
                   mother(0), sdat(0), own(true), valid(true), armed(false), childExists(false), nvc(0),
                   vc(0), acc(0), splintcC(0), splgValid(false), lastarm(0),
                   gen(std::numeric_limits<double>::infinity()),
                   checkSm(std::numeric_limits<double>::infinity()),
                   checkSm2(std::numeric_limits<double>::infinity()),
                   killChildrenWhenDying(no)
{
#ifdef _OPENMP
  omp_init_lock(&mElementAccessLock);
#endif
  xdat = new double[maxDim];
  ydat = new double[maxDim];
  if (xdat && ydat) getData(s); 
  else throw Bad_Error("Spline::Spline(const Spline& s, string nom, Anchor* a)  out of memory");
#ifdef ENABLE_SPLINE_WATCH
  WatchCount++;
  WatchMap[name]++;
#endif
}

Spline::Spline(ifstream &in, string nom, Anchor *a)
                : AnchorEnabled(a), n(-1), name(nom), maxDim(0),
                  mother(0), sdat(0), own(true), valid(true), armed(false),
                  childExists(false), nvc(0), vc(0), acc(0), splintcC(0), splgValid(false),
                  lastarm(0),
                  gen(std::numeric_limits<double>::infinity()),
                  checkSm(std::numeric_limits<double>::infinity()),
                  checkSm2(std::numeric_limits<double>::infinity()),
                  killChildrenWhenDying(no)
{
#ifdef _OPENMP
  omp_init_lock(&mElementAccessLock);
#endif
  int FormerN;
  double x,y;
  in.read((char*)&FormerN,sizeof(int));
  //cout << "SIZE: " << FormerN << endl;
  maxDim = FormerN+1;
  xdat = new double[maxDim];
  ydat = new double[maxDim];

  if (xdat == 0 || ydat == 0 ) throw Bad_Error("no memory");
  for (int k = 0; k <= FormerN; k ++) {
    in.read((char*)&x,sizeof(double));
    in.read((char*)&y,sizeof(double));
    set(x,y);
  }
#ifdef ENABLE_SPLINE_WATCH
  WatchCount++;
  WatchMap[name]++;
#endif
}

/*!
  Create a spline as child of mother s using (naturally) mother x-data
  and y-data coming from calls of moSingle func(x-data)
  This is an ultra - convenient way of generating mappings of splines
  using arbitrary functions :-)
*/
Spline::Spline(moSingle func, const Mathobject& mobj, Spline *s, string nom,Anchor* a)
                 : AnchorEnabled(a), n(s->n), name(nom), maxDim(s->maxDim), mother(s),
                   xdat(0), sdat(0), own(false), valid(s->valid), armed(false), childExists(false), splgValid(false),
                   lastarm(0),
                   gen(std::numeric_limits<double>::infinity()),
                   checkSm(std::numeric_limits<double>::infinity()),
                   checkSm2(std::numeric_limits<double>::infinity()),
                   killChildrenWhenDying(no)
{
#ifdef _OPENMP
  omp_init_lock(&mElementAccessLock);
#endif
    ydat = new double[s->maxDim];
    if (ydat == 0 ) error("Spline(Spline*)","no memory");    
    for (int i = 0; i < s->maxDim; i++)  ydat[i] = 0;
    if (s->own) s->addChild(this); 
    else {
      mother = s->mother;
      mother -> addChild(this);
    }

    for (int i = 0; i < s->n; i++) 
      ydat[i] = (mobj.*func)(x(i));
#ifdef ENABLE_SPLINE_WATCH
  WatchCount++;
  WatchMap[name]++;
#endif  
}

Spline* Spline::fromTextFile(const string& fileName, string splineName, Anchor* a)
{
  ifstream in(fileName.c_str());
  if (in) {
    return Spline::fromTextTable(in, splineName, a);
  }
  throw Bad_Error("Spline::fromTextFile() - Could not read fileName");
}

Spline* Spline::fromTextTable(ifstream& i, string splineName, Anchor* a)
{
  typedef vector<pair<double, double> > Points;
  Points points;
  double x, y;
  while (i >> x >> y) {
    i.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    points.push_back(make_pair<double, double>(x, y));
  }
  Spline* s = new Spline(points.size(), splineName, a);
  Points::const_iterator it, end = points.end();
  for (it = points.begin(); it != end; ++it) {
    s->set(it->first, it->second);
  }
  return s;
}


Spline::~Spline() {
#ifdef _OPENMP
  omp_destroy_lock(&mElementAccessLock);
#endif
  // if (armed) checksum();

  // 
  // when we exit, we check if we have to kill children too,
  // or just information is required (in rare case you wish to construct another mother)
  // 

  if (killChildrenWhenDying == all) {
    //cout << "deleting all.." << endl;
    while (!children.empty()) {
      // the next two lines are precautions. We really
      // want the daugther spline to aks us to kill it from our
      // child - list (consistency) [see removeChild call below]
      Spline *s = children.begin()->second;
      s->mother = this;
      s->own = false;
      delete children.begin()->second;  // this will remove the spline from children map also
    }
    //cout << "finished deleting all" << endl;
  } else {
    for(map<Spline*,Spline*>::iterator i = children.begin(); i != children.end(); i++)
      i->second->motherKilled();
  }

  if (mother !=0 && !own) {
    //cout << "~Spline:: I inform mother " << name << "   " << mother << endl;
    mother -> removeChild(this); 
    //cout << "I did inform " << endl;
  }

  if (own) delete[] xdat; 
  delete[] ydat; 
  delete[] sdat; 
  xdat = ydat = sdat = 0;
  //cout << " ... accomplished" << endl;

#ifdef ENABLE_SPLINE_WATCH
  WatchCount--;
  WatchMap[name]--;
#endif
  //this->AnchorEnabled::~AnchorEnabled();
  //cout << "exiting ~Spline for " << name << "   " << this << endl;
}


void Spline::motherKilled() {
  if (xdat != 0) error("motherKilled", "strange");
  valid = false;
  armed = false;
  mother = 0;
}

Spline* Spline::getChild(string s) {
  for (map<Spline*,Spline*>::iterator i = children.begin(); i != children.end(); i++) 
    if (i->second->name == s) return i->second;
  error("getChild", "not child named " + s);
  return 0; // to silence compiler
}

void Spline::rearm(bool keepSize) {
  arm(no,keepSize,lastarm);
  lastarm = n;
}

void Spline::arm(childrenToo childtoo, bool keepSize, int anf) {
  if (armed && anf==0)
    error("arm", "Spline already armed");
  if (n <= 1)
    error("arm", "Trying to arm a spline with less than 2 points");
  proper();

  if (keepSize) {  // if the arrays should be kept
    if (sdat == 0) sdat = new double[maxDim];  // we take the size of xdat  and ydat...
  } else {
      resize(n);      // if we are not forced to keep the size of the arrays, we will adopt minimum size
      delete[] sdat;
      sdat = new double[n+1];
  }
  if (sdat == 0) throw Bad_Error("no memory");   

  if (own) {
    if (valid) {
      origSpline(xdat,ydat,anf, n,1e30,1e30,sdat);
      armed = true;
      mCache.initGuess();
    }  else throw  Bad_Error("arm() not valid");
  } else {
    if (!mother->valid)  throw  Bad_Error("mother is not not valid");     
    if (mother->n != n) error("arm", "mother and child have different number of datapoints");  
    origSpline(mother->xdat,ydat, anf, n,1e30,1e30,sdat);
    armed = true;
  }
  //createChecksum();

  // if there are children, we may have to arm those (or those with same n) too...
  if (childExists) {
    switch (childtoo) {
    case all:
      for (map< Spline*, Spline*>::iterator i = children.begin(); i != children.end(); i++) i->second->arm();
      break;
    case thoseReady:
      for (map< Spline*, Spline*>::iterator i = children.begin(); i != children.end(); i++) {
	if (i->second->n == n)  i->second->arm();
      }
      break; 
    default: ;
    }
  }
}

void Spline::disarm(childrenToo childtoo, bool keepSize) {
  if (!armed) return; 
  if (!keepSize) {  delete[] sdat; sdat = 0; }

  if (childExists && childtoo == all) {
    for (map< Spline*, Spline*>::iterator i = children.begin(); i != children.end(); i++)
      i->second->disarm();
  }
  if (childExists && childtoo == thoseReady) {
    for (map< Spline*, Spline*>::iterator i = children.begin(); i != children.end(); i++)
     if (i->second->n==n)
       i->second->disarm();
  }
  armed=false;
  if (own)
    mCache.invalidate();
}

void Spline::proper() {
  if (n==0) {
    error("proper", "only one point in spline.");
  }
  for (int i =0; i < n; i++) {
    if (x(i+1) < x(i)) {
      cout << "for i=" << i << " spline error for spline named: " << name << endl;
      cout << "hint: x(i) = " << x(i) << "  and x(i+1) : " << x(i+1) << endl;
      dump("splineNoOrder",true);
      error("proper","the x-data is not in proper order [ x(i+1) is < x(i)].\nPlease make sure x(i+1) > x(i)"); 
    }
    if (x(i+1) == x(i)) {
      cout << "for i=" << i << " spline error for spline named: " << name << endl;
      cout << "hint: x(i) = " << x(i) << endl;
      dump("splineNoOrder",true);
      error("proper","the x-data is not in proper order [ x(i+1) is = x(i)].\nPlease make sure x(i+1) > x(i)");
    }
  }
}

struct XIndexSorter
{
  XIndexSorter(double *xdat): mXdat(xdat) {}
  bool operator()(unsigned int i, unsigned int j) {
    return mXdat[i] < mXdat[j];
  }
private:
  double* mXdat;
};

struct XdatEqual
{
  XdatEqual(double *xdat): mXdat(xdat) {}
  bool operator()(unsigned int i, unsigned int j) {
    return mXdat[i] == mXdat[j];
  }
private:
  double* mXdat;
};

void Spline::makeProper()
{
  if (!own) {
    error("makeProper()", "spline must have its own x-data (i.e. no mother) for calling makeProper()", -1);
  }
  if (!valid) {
    error("makeProper()", "spline not valid", -1);
  }
  bool wasArmed = armed;
  if (armed) {
    disarm();
  }

  std::vector<unsigned long int> indices;
  for (int i =0; i <= n; ++i) {
    indices.push_back(i);
  }

  std::vector<unsigned long int>::iterator newIndexEnd;
  std::vector<unsigned long int>::iterator it = indices.begin();
  std::sort(it, indices.end(), XIndexSorter(xdat));
  newIndexEnd = std::unique(it, indices.end(), XdatEqual(xdat));
  indices.erase(newIndexEnd, indices.end());

  double *tmpX = new double[n+1];
  double *tmpY = new double[n+1];
  map<Spline*, double*> childrensTmpYdat;
  map< Spline*, Spline*>::iterator childIt, childrenEnd = children.end();
  if (childExists) {
    for (childIt = children.begin(); childIt != childrenEnd; ++childIt) {
      if (childIt->second->n==n)
        childrensTmpYdat[childIt->first] = new double[n+1];
      else
        childrensTmpYdat[childIt->first] = 0;
    }
  }
  unsigned long int i=0;
  for ( it = indices.begin(); it != newIndexEnd; ++it) {
    tmpX[i] = xdat[*it];
    tmpY[i] = ydat[*it];
    if (childExists) {
      for (childIt = children.begin(); childIt != childrenEnd; ++childIt) {
        double *childTmpY = childrensTmpYdat[childIt->first];
        if (childTmpY)
          childTmpY[i]=childIt->second->ydat[*it];
      }
    }
    ++i;
  }
  delete[] xdat;
  delete[] ydat;
  xdat = tmpX;
  ydat = tmpY;
  maxDim = indices.size();
  n = maxDim-1;
  if (childExists) {
    for (childIt = children.begin(); childIt != childrenEnd; ++childIt) {
      double *childTmpY = childrensTmpYdat[childIt->first];
      if (childTmpY) {
        delete[] childIt->second->ydat;
        childIt->second->ydat = childTmpY;
        childIt->second->maxDim = maxDim;
        childIt->second->n = n;
      }
    }
  }
  if (wasArmed) {
    arm(thoseReady);
  }
}


/*!
  Wrapper for dump(string,bool). Reason: I can never
  remember, if bool has to be true or false to get in addition
  to name.dat a file name.plt with interpolated spline values
  in addition.

  So I now use the plain-text attribute "yes" and "no"
*/
void Spline::dump(string nam,childrenToo InterpolatedAlso) {
  dump(nam,InterpolatedAlso==no);
}

/*! 
  Write spline data to file called "nam.dat".
  If bool only is false, it will also create a file "nam.plt"
  containing 10 times as many points interpolated
  from this spline 
*/
void Spline::dump(string nam, bool only) { 
  string daten = nam + ".dat";
  string plot = nam + ".plt";
  const char *na = daten.c_str();
  const char *na2 =  plot.c_str();
  ofstream d(na);
  ofstream d2;

  if (!only) {
    d2.open(na2);
    if (!d2) {
      throw Bad_Error("Spline::dump() - failed to write to " + plot);
    }
    dump(d2,false);
  }
  if (!d) {
    throw Bad_Error("Spline::dump() - failed to write to " + daten);
  }
  dump(d,true);
}

void Spline::dump() {
  if (name.empty()) {
      error("dump()", "Spline doesn't have a name, so you need to specify one if you use dump() without a filename.");
  }
  dump(name);
}

/*! 
  Write spline data to ofstream o,
  If bool only is false, 10 interpolated points per
  spline - data point are written (to give a smoother
  picture
*/
void Spline::dump(ofstream &o, bool only) {
  o.setf(ios::scientific);
  o << std::setprecision(12);
  for (int k = 0; k <= n; k++) {
    o << x(k) << "  " << ydat[k] << endl;
    if (k < n && !only) {
      for (double xx = x(k)+(x(k+1)-x(k))*0.1; xx < x(k+1); xx+= (x(k+1)-x(k))*0.1) {
	o << xx << "  " << (*this)(xx) << endl;
      }
    }
  }
}

void Spline::save(const string &name) {
  ofstream out(name.c_str());
  save(out);
}

void Spline::save() {
  if (name.empty())
    throw Bad_Error("Spline::save() - spline needs to have a name when using save() without a filename");
  save(name+".spline");
}

void Spline::save(ofstream &o) {
  static double X,Y;
  o.write((char*)&n, sizeof(int));
  for (int k =0; k<=n; k++) {
    X = x(k);
    Y = y(k);
    o.write((char*)&X,sizeof(double));
    o.write((char*)&Y,sizeof(double));
  }
}



/*! Kept for debugging purposes. If You think
  something very strange is going on, you can 
  ask the spline to build checksums of its y and s data.
  Hence you will be able to notice memory corruptions
*/
void Spline::createChecksum() {
    checkSm = 0;
    for (int k = 0; k <= n; k++) {
      if (own) checkSm += ck(xdat[k]); 
      checkSm += ck(ydat[k]);
      checkSm += ck(sdat[k]);
    }
    checkSm2 = 0;
    for (int k = 0; k <= n; k++) {
      if (own) checkSm2 += ck(xdat[k]); 
      checkSm2 += ck(ydat[k]);
    }
}

/*!
  After you have creatChecksum()ed, you can
  call checksum(). It will reevaluate and if
  the checksums do not coincide, throw a Bad_Error
*/
void Spline::checksum() {
  double chk =0,chk2 =0;
   for (int k = 0; k <= n; k++) {   
     if (own) { chk += ck(xdat[k]); chk2 += ck(xdat[k]);}
     chk += ck(ydat[k]);
     chk2 += ck(ydat[k]);
     chk +=ck(sdat[k]);
   }
   
   if (chk != checkSm || chk2 != checkSm2) {
     cout << "checksum" << checkSm << "  chk: " << chk << endl; 
     cout << "checksum2 " << checkSm2 << "  chk2: " << chk2 << endl; 
     cout << "chk - checkSm: " << chk - checkSm << endl;
     cout << "chk2 - checkSm2: " << chk2 - checkSm2 << endl;
     error("Checksum","not identical"); 
   }
}
  



void Spline::set(const double x, const double y) {
  checkSpace();
  if (!own) error("set(x,y)","add to something not owned");
  xdat[++n] = x;
  ydat[n] = y;
}

void Spline::set(const double y) {
  checkSpace();
  if (own) { printStatus(); throw Bad_Error("set:: no mother, so x data also needed"); }
  ydat[++n] = y;
}




void Spline::set(const map<float,float>& m) {
  for ( map<float,float>::const_iterator i= m.begin(); i != m.end(); i++) setForce(i->first,i->second);
} 

void Spline::set(const map<double,double>& m) {
  for ( map<double,double>::const_iterator i= m.begin(); i != m.end(); i++) setForce(i->first,i->second);
} 

void Spline::setX(int k, double x) {
  if (!own)  error("setX", "This spline doesn't own xdata, it belongs to the mother: change xdata there");
  checkSpace(k-1);  // k-1, cause checkspace assummes that k+1 is written to, see set()
  xdat[k] = x;
  if (k> n) n=k;
}
void Spline::setY(int k, double x) {
  checkSpace(k-1);// k-1, cause checkspace assummes that k+1 is written to, see set()
  ydat[k] = x;
  if (k> n) n=k;
}

/*! 
  This routine lets you set the counter of the data points of this spline EXPLICITLY.
  No need to mention that this is dangerous, not nice and should be avoided. 
  However, there are two closely related exceptions:

  (1) If the data points starting from some x value on are all zero, this can be useful,
  as by default the new y-data is always zero'd when allocated.

  (2) If the spline is not really initialized, but the daugther of some other spline
  and for the sake of simple programming it has to be accessed nevertheless (to
  avoid a lot of if thens), the access delivering zero, of course,
  then a call to SetN() WITHOUT an argument will set the
  counter to the n value of the mother which is assumed to contain the amount
  of data you wish to have controle of here 

  setN() will return error, if n is out of maxDim range. This is a safety precaution.
  use resize() explicitly, if you really need to do such strange things. In general:

  Do not use this function easily, try set, setX, setY, ptrY etc, which all have
  checks, increase n if necessary and so on and so on and so on

*/

void Spline::setN(int k) {
  if (k == -2 && !own && mother != 0) k = mother->n;
  if (k<-1) error("setN","k to small",k);
  if (k>=maxDim) error("setN","k larger than maxDim",k);
  n = k;
}

void Spline::setChildrensN(int k) {
  for (map< Spline*, Spline*>::iterator i = children.begin(); i != children.end(); i++) i->second->setN(k);
}
  

void Spline::getData(const Spline &s, childrenToo sdataToo) {
  int sizze = s.size();
  if (sizze > maxDim) error("getData","spline to copy data from is too large");
  
  for (int i = 0; i < sizze; i++) {
    if (own) xdat[i] = s.x(i);
    ydat[i] = s.y(i);
  }
  if (sdataToo == all && s.armed && s.sdat != 0) {
    for (int i = 0; i < sizze; i++) sdat[i] = s.sdat[i];
  }
  n = s.n;
}


/*!
  resize xdat and ydat to fit index k, if k > maxDim. If this spline is a mother,
  it calls resize() for all children
 */
void Spline::resize(int k) {
  //cout << "resize: " << name;
  if (k < n) error("resize","desired size too small to hold data",k);
  if (k+1 == maxDim) return; // nothing to do
  double *tmp;
  //cout << "I am resizing myself: " << name << "   " << maxDim << "  to: " << k << endl;
  //cout << "...." << endl;
  if (own && valid) {
    //cout << "own and valid :" << n << endl;
    tmp = new double[k+1];
    if (tmp == 0) throw Bad_Error("no memory i");
    for (int i =0; i<= n; i++) tmp[i]=xdat[i]; 
    //cout << "deleting xdat: " << xdat << endl;
    delete[] xdat; 
    //cout << "done" << endl;
    xdat = tmp;

    if (childExists) {
      // cout << "... children exist" << endl;
      for (map< Spline*, Spline*>::iterator i = children.begin(); i != children.end(); i++) i->second->resize(k);
    }
  }
  
  tmp = new double[k+1];
  if (tmp == 0) throw Bad_Error("no memory ii");
  for (int i =n+1; i <= k; i++) tmp[i]=0; // initialize
  for (int i =0; i<=n;i++) tmp[i]=ydat[i]; // copy old data
  delete[] ydat;
  ydat = tmp;
  maxDim = k+1;

  if (sdat != 0) {  delete[] sdat; sdat = 0;}
 
}

void Spline::checkSpace(int k) {
  if (!valid)  throw Bad_Error("Spline::checkSpace() - spline not valid");

  if (k == -2) k = n;
  if (k < -1) error("checkSpace", "k smaller than -1 is not accepted",k);
  if (k + 1 < maxDim) return;
  //cout << name <<"  doubling " << n << " maxDim: " << maxDim << endl;
  int need = maxDim*2;  // double the size ...
  if (need < k) need = k;    // if this is still to small, adopt to k
  resize(need);  // resize to fit at indices up to index need
}


double Spline::randomY(double x) {
 if (!valid)  Bad_Error("not valid");     
 if (!armed)  error("randomY","not armed");     

 unsigned int fireandforget;

  if (own) {
    return splint(xdat, x, &fireandforget);
  } else {
    // if we are valid, so must be mummy
    return splint(mother->xdat, x, &fireandforget);
  }
}

/*!
  moSingle for Miscmath::zbrent() for finding the "zeros" of 
  a spline. Zero is relative, cause the static double zeroLevel
  will be subtracted from the y(x). 
  So for instance: if y(x) = 3 and zeroLevel = 2 then  one is 
  returnde. Usually for finding really zero zero's, zeroLevel =0
*/
double Spline::inBetweenY(double x) { 
  double *xa;
  if (own) xa = xdat; else xa = mother->xdat;
  const CachedValues& cache = cachedVals();
  double h = xa[cache.khi]-xa[cache.klo];
  double a=(xa[cache.khi]-x)/h;
  double b=(x-xa[cache.klo])/h;
  return a*ydat[cache.klo]+b*ydat[cache.khi]+((a*a*a-a)*sdat[cache.klo]
                                             +(b*b*b-b)*sdat[cache.khi])
         *(h*h)/6.0 - zeroLevel;
}


void Spline::printStatus() const {

  string s[] = {"false","true"};
  cout << endl;
  cout << "============================================" << endl;
  cout <<  "Name: " << name << "  at: " << &(*this) << endl;
  cout <<  "generated at: " << generatedAt() << endl;
  cout << "own: " << s[own] << endl;
  cout << "valid: " << s[valid] << endl;
  cout << "armed: " << s[armed] << endl;
  cout << "--------------------------------------------" << endl;
  cout << "n: " << n << endl;
  cout << "maxDim:: " << maxDim << endl;
  cout << "--------------------------------------------" << endl;
  cout << "children: " << children.size() << endl;
  if (own) { 
    cout << "no mother" << endl;
    cout << "calls to () " << acc << endl;
    cout << "cache hit:   " << vc << endl;
    cout << "cache miss: " << nvc << endl;
    cout << "ratio (h/m):    " << (double)(vc)/(double)(nvc) << endl;
    cout << "calls to splint(): " << splintcC << endl;
  }
  else cout << "mother: " << mother->name << "   at: " << mother << endl;
  cout << "--------------------------------------------" << endl;
  if (own) cout << "xdat: " << xdat << endl;
  cout << "ydat: " << ydat << endl;
  cout << "sdat: " << sdat << endl;
  cout << "--------------------------------------------" << endl;
  cout << "checksum : " << checkSm << endl;
  cout << "checksum2 : " << checkSm2 << endl;
  cout << "--------------------------------------------" << endl;
  // cout << "x-range: " << start() << "  ...  " << stop() << endl;

  if (own && xdat != 0)  cout << "x[1] : " << xdat[0] <<  "   x[n]:   " << xdat[n] << endl;
  if (ydat != 0) cout << "y[1] : " << ydat[0] <<  "   y[n]:   " << ydat[n] << endl;
  if ((n>1) && (sdat != 0))  cout << "s[2] : " << sdat[1] <<  "   s[n-1]: " << sdat[n-1] << endl;

  cout << "============================================" << endl;
  cout << endl;

}

void Spline::setMother(Spline& m) {  
  if (!own) throw Bad_Error("Spline already has a mother");
  if (!m.own) throw Bad_Error("mother is itself a child");
  mother = &m;
  if (xdat != 0) { delete[] xdat; xdat =0; }
  own = false;
  armed = false;
  valid = m.valid;
  m.addChild(this);
}

void Spline::splder(Spline &s) {
  if (s.ydat != 0)  delete[] s.ydat; 
    
  s.ydat  = new double[n+1];   

  if (s.own && s.xdat !=0){ 
    throw Bad_Error("Not implemented");
  } 
  s.n = n;
  s.maxDim = n + 1;
  s.valid = valid;
  Splder(ydat,s.ydat, n);
}

/*!
  Fill spline s with data of derivatives of this spline.
  In contrast to splder(), no equally spaced points
  are assumed and the spline function is evaluated,
  nothing else....

  The x-data does not have to coincide (spline s 
  should however lie in the definition range of this 
  spline)

  If spline s has no xdata, the xdata of this spline
  is taken.
*/

void Spline::derive(Spline &s, childrenToo  ct) {
  if (s.own) {       // if it maybe has own
    delete[] s.xdat;  // copy our xdat
    s.xdat = new double[n+1];
    for (int i = 0; i <= n; i++) s.xdat[i] = x(i);
    s.valid = true;
    s.n = n;
  } else {
    s.n = n;   // if no own data, set s.n
  }

  if (s.ydat != 0) delete[] s.ydat; 
  s.ydat = new double [s.n+1];

  s.maxDim = s.n+1;

  // so until now, the spline has (if necessary) new x-data and 
  // always new y -data  fitting s.n
  // we now fill y-data with the derivative of this spline

  double dummy;

  for (int k = 1; k <  s.n; k++) {  // first the inner part, the boundaries later
    s.ydat[k] = Miscmath::dfridr( (moSingle) &Spline::safeY, *this, s.x(k) , (s.x(k+1)- s.x(k)), &dummy); 
    //cout << "f'("<<s.x(k) << ")  =  " <<      s.ydat[k] << "   with error: " <<  dummy << endl;
  }

  s.ydat[0] = Miscmath::dfridr( (moSingle) &Spline::safeY, *this, s.start() + 1e-2*(s.start(1)-s.start()), 1e-3*(s.start(1)-s.start()), &dummy); 
  s.ydat[s.n] = Miscmath::dfridr( (moSingle) &Spline::safeY, *this, s.stop() - 1e-2*(s.stop()-s.stop(1)), 1e-3*(s.stop()-s.stop(1)), &dummy); 
  
  if (ct != no) s.arm();
}


double Spline::safeY(const double x) {
  if (x > start() ) {
    if (x < stop()) return fastY(x);
    return back();
  } else return front();
}

bool Spline::isWithinBounds(const double x) {
  if (x >= start() && x <= stop()) return true;
  return false;
}


/*! 
  Convolute this Spline with a partner Spline p. The Convolution starts at start
  and ends at end, and has the tolerance tol.

  In principle this looks like:   \int_start^end  this(x)  times partner(fak * x + shft) dx

  So fak(tor) and sh(i)ft give additional freedom to convolute shifted and streched
  partners.
*/

double Spline::convolution(Spline* p,double start, double end,double tol,double fak, double shft) {
  throw Bad_Error("convoltution no support");
  checkConvolutionRange(p, &start, &end, fak, shft);
 
  double y[2];   y[1] = 0;  
  Miscmath::rungeInt(y, 1, start,end, 1e-4, (end-start)*0.01,0,  (moDerivs) &Spline::convD, *this);
  return y[1];
}

//! convolution function for single convolution type: moDerivs
void Spline::convD(const double x,const  double *y, double* dy)  {
  cC++;
  double px = x*faktor + shift;
  dy[1] = (*partner)(px) * (*this)(x);  
  throw Bad_Error("convd no support");
}
/*
double Spline::rombConvolution(Spline* p,double start, double end,double tol,double fak, double shft) {
  checkConvolutionRange(p, &start, &end, fak, shft);
  throw Bad_Error("rombconvoltution no support");
  return  Miscmath::rombint( (moSingle) &Spline::conv, this, start,end, tol);
}
 
//! convolution function for rombConvolution
double Spline::conv(const double x)  {
  double px = x*faktor + shift;
   throw Bad_Error("conv romb convoltution no support");
  return (*partner)(px) * (*this)(x);
}
*/

/*! 
  Performs convolution of a whole family with one partner spline p.
  Starting with mother, it climbs down the map of children... (in convFamily() )
*/
/*
void Spline::familyConvolution(Spline *p, double start, double end,double tol,double fak, double shft) {
  cC = 0; throw Bad_Error("family convoltution no support");
  checkConvolutionRange(p, &start, &end, fak, shft);
  double *y = new double[children.size() + 2];
  for (unsigned int i =1; i <= children.size() + 1 ; i++) y[i] =0;
  Miscmath::rungeInt(y, 1+children.size(), start,end,tol, (end-start)*0.01,0,  (moDerivs) &Spline::convFamily, *this);
  // cout << "convolution took: " << cC <<endl;

  int k =2;
  convResult = y[1];
  for (map<Spline*,Spline*>::iterator i = children.begin(); i != children.end(); k++,i++) 
    i->second->convResult = y[k];
  
  delete[] y;
}
*/

//! convolution function for family convolution type: moDerivs
/*
void Spline::convFamily(const double x,const  double *y, double* dy)  {
  cC++; throw Bad_Error("conv family convoltution no support");
  double px = x*faktor + shift;
  double py =  (*partner)(px);
  dy[1] = py * fastSplint(xdat,x, &guess);
  int k = 2;
  for (map<Spline*,Spline*>::iterator i = children.begin(); i != children.end(); k++, i++) { 
    dy[k] =py * (* (i->first))(x);
  }
}
*/



/*!
  Optimized convolution of several splines with S one partner spline P.  

  Effecient in the case that several splines have the same x-data (this need not be owned by themselves). The splines S are store in the static vector convolutionVektor[] which you have to resize and fill with the Spline S. 
  
  vektorConvolution() will then temporariliy make convolutionVektor[0] the mother of all other splines S. Hence, access to the splines is greatly sped up by caching of interpolation variables etc.
*/
void Spline::vectorConvolution2(Spline* p , double start, double end,double tol,double fak, double shft,double hnext) {  
  Spline * null = convolutionVektor[0]; 
  null->checkConvolutionRange(p, &start, &end, fak, shft);
  null->cachedVals().guess=0;

  if (hnext == 0) hnext =  (end-start)*1e-3;

  //bool was_own = null->own;
  nvar  =  convolutionVektor.size();

  // remember the old settings 
  bool was_childExists = null->childExists;
  double * was_xdat = null->xdat;
  vector<Spline*> was_mother(nvar+1);
  vector<bool> was_own(nvar+1);

  double *y = new double[nvar+1];
  for (unsigned int i =0; i < nvar ; i++) {
    if (!convolutionVektor[i]->isArmed()) convolutionVektor[i]->error("vectorConvolution2","splines need to be armed");
    y[i+1] =0; 
    was_mother[i] = convolutionVektor[i]->mother;
    was_own[i] = convolutionVektor[i]->own;
    if (i) {  // so we are not null
      convolutionVektor[i] -> mother = null;    // null is now mother
      convolutionVektor[i] -> own = false;      // and of course we do not own
    }
  }

  if (! null->own) null->xdat = null->mother->xdat;

  null->own = true;
  null->childExists = true;

  Miscmath::rungeInt(y, nvar, start,end,tol,hnext,0,  (moDerivs) &Spline::convVektor2, *null);
  for (unsigned int i =0; i < convolutionVektor.size(); i++) convolutionVektor[i]->convResult = y[i+1];  
  delete[] y;

  null->childExists = was_childExists;
  null->xdat = was_xdat; 

  for (unsigned int i =0; i < nvar ; i++) {
    convolutionVektor[i]->mother = was_mother[i];
    convolutionVektor[i]->own = was_own[i];
  }
}


void Spline::convVektor2(const double x,const  double *y, double* dy)  { 
  double px = x*faktor + shift, py;

  if (partner->own) py = partner->fastSplint(partner->xdat, px, &partner->cachedVals().guess); else py =  (*partner)(px);

  CachedValues& cache = cachedVals();
  dy[1] = py * fastSplint(xdat,x, &cache.guess); // this automatically fills the cache variables (if x is in bound)
  //printStatus();
  if (cache.valid) {
    for (unsigned int i = 2; i <= nvar; i++) { 
      double& h = cache.h;
      double& a = cache.a;
      double& b = cache.b;
      double& a3a = cache.a3a;
      double& b3b = cache.b3b;
      double& h26 = cache.h26;
      double& lastX = cache.lastX;
      unsigned int& khi = cache.khi;
      unsigned int& klo = cache.klo;
      // dy[i] = py*(*convolutionVektor[i-1])(x);
      Spline *s = convolutionVektor[i-1];
      dy[i] = (a*s->ydat[klo]+b*s->ydat[khi]+(a3a*s->sdat[klo]+ b3b*s->sdat[khi])*h26)*py;
    }
  } else {
    for (unsigned int i = 2; i <= nvar; i++) 
      dy[i]=py*convolutionVektor[i-1]->fastSplint(xdat,x, &cache.guess);
  }
}

/*!
  Special version of convVektor2() that doesn't need vektor access cause it is 
  optimized for two splines
  OBSOLETE
*/
/*
void Spline::convVektorTwin(const double x,const  double *y, double* dy)  { 
  double px = x*faktor + shift, py; 
  //cout << "deb1" << endl;
  throw Bad_Error("conv Vektor Twin  no support");
  if (partner->own) py = partner->fastSplint(partner->xdat, px, &partner->guess); else py =  (*partner)(px);
  //cout << "deb2" << endl;
  dy[1] = py * fastSplint(xdat,x, &guess); // this automatically fills the cache variables (if x is in bound)
  // cout << "deb3" << endl;
  
//   if (convTwin -> mother) {
//     cout << "there is mother" << endl;
//     convTwin->mother->printStatus();
//     if (convTwin->mother != this) throw Bad_Error("mother is not me");
//   }

  dy[2] = py * (*convTwin)(x);
  // cout << "deb4" << endl;
  }

*/

/*! 
  Fill vector <Spline*> with this family
*/
/*
void Spline::vectorConvolution(Spline *p, double start, double end,double tol,double fak, double shft) {
  cC = 0;  
  throw Bad_Error("vector convoltution no support");
  checkConvolutionRange(p, &start, &end, fak, shft);
  nvar  =  1 + children.size();
  guess = 1;
  double *y = new double[nvar+1];
  
  for (unsigned int i =1; i <= nvar ; i++) y[i] =0;
  if (convolutionVektor.size() != nvar) convolutionVektor.resize(nvar);
  
  int k = 0;
  for (map<Spline*,Spline*>::iterator i  = children.begin(); i != children.end(); i++,k++) {
    convolutionVektor[k] = i->second;
  }
  
  
  Miscmath::rungeInt(y, nvar, start,end,tol, (end-start)*1e-3,0,  (moDerivs) &Spline::convVektor, *this);
  
  
  // cout << "convolution took: " << cC <<endl;
  
  convResult = y[1];
  for (unsigned int i = 2; i <= nvar; i ++) convolutionVektor[i-2]->convResult = y[i];
  delete[] y;
}
*/




/*! 
  Print the convolution integrand of this and partner to file datei
*/
void Spline::dumpConvolution(Spline *p, double start, double end,double fak, double shft, const char* datei) {
  throw Bad_Error("dump convoltution no support");
  cout << "request: " << start << " ..." << end;
  checkConvolutionRange(p, &start, &end, fak, shft);
 
  double step = (end-start)*1e-4;
  ofstream ausgabe(datei);
  ofstream par("partner");


  for (double x = start; x < end; x+= step) {
    double px = x*fak + shft;
    ausgabe << x << "  " <<   (*partner)(px) * (*this)(x) << endl;
    par << x << " " <<   (*partner)(px) << endl;
  }
  ausgabe.flush();
  cout << "weiter? "; string weiter;
  cin >> weiter;
  
}


/*
void Spline::convVektor_backup(const double x,const  double *y, double* dy)  { 
  cC++;
   throw Bad_Error("convvecotr backup convoltution no support");
  double px = x*faktor + shift, py;
  if (partner->own) py = partner->fastSplint(partner->xdat, px, &partner->guess); else py =  (*partner)(px);

  dy[1] = py * fastSplint(xdat,x, &guess); // this automatically fills the cache variables (if x is in bound)
  for (unsigned int i = 2; i <= nvar; i++)  {
    if (validCache) {        // usually this is true, however, if out of bound we have to call operator()
      Spline *s = convolutionVektor[i-2];
      dy[i] = py *  (a * s->ydat[klo] + b* s->ydat[khi]+(a3a*s->sdat[klo]+ b3b*s->sdat[khi])*h26 );
    } else dy[i] =py * (*(convolutionVektor[i-2]))(x);	
  }
}
*/

/*
void Spline::convVektor(const double x,const  double *y, double* dy)  {
  // cC++;
  double px = x*faktor + shift,py = (*partner)(px);
  double *yptr,*sptr;
  throw Bad_Error("conv vector 22 convoltution no support");
  // the range checking has been done by checkConvolutionRange
  // so we only have to do a small version of the fastSplint() story
  dy++;
  *dy = py * fastSplint(xdat,x, &guess); // this automatically fills the cache variables (if x is in bound)
 
  for (unsigned int i = 2; i <= nvar; i++)  {
    dy++;  
    if (validCache) {        // usually this is true, however, if out of bound we have to call operator()
      Spline *s = convolutionVektor[i-2];
      yptr = &(s->ydat[klo]);
      sptr = &(s->sdat[klo]);

      *dy = py *  (a * (*yptr) + b* (*(++yptr)) +(a3a*(*sptr)+ b3b*(*(++sptr)))*h26 );

    } else *dy =py * (*(convolutionVektor[i-2]))(x);	
  }
}
*/


double Spline::convOneDim(const double x) {
  cC++;
  double px = x*faktor + shift;
  throw Bad_Error("convone dim convoltution no support");
  if (px < partner->x(1)) error("convD","px<");
  if (px > partner->x(partner->n))  error("convD","px<");
  
  return  (*partner)(px) * (*this)(x);
}



void Spline::checkConvolutionRange(Spline *p, double *start, double *end, const double fak, const double shft) {
  double partStart,partStop;
  partner = p;
  shift = shft;
  faktor = fak;
  
  cC = 0;

  if (*start == *end) {
    WARN("Spline::CheckconvolutionRange() start == end");
    return;
  }
 
  if (fak > 0) {
    partStart = (p->start() - shft)/fak;
    partStop = (p->stop() - shft)/fak; 
  } else {
    partStart = (p->stop() - shft)/fak;
    partStop = (p->start() - shft)/fak;
  } 

  *start = max(*start, x(0) );  *start = max(*start,partStart);   //finding the maximal start
  *end = min(*end, x(n) );      *end = min(*end,partStop);       // and the minimal end (i.e. the total overlap :-)
  if (*start >= *end) throw noConvolutionOverlap();
}



/*! 
  Integrate from double a to double b by performing a
  runge - kutta integration, i.e. converting the problem
  to the solution of a ordinary diff-eqn. 
*/

double Spline::integrate(double a, double b, double tol) {
  double y[2];
  y[1] = 0;
  double initialGuess = 5e-3;
  Miscmath::rungeInt(y, 1, a, b, tol, (a-b)*initialGuess, 0, (moDerivs) &Spline::inteD, *this);
//X   Miscmath::odeint(y, 1, a, b, tol, (a-b)*initialGuess, 0, (moDerivs) &Spline::inteD, *this);
  return y[1];
  return Miscmath::oneDimOdeint(a, b, 1e-5, (moSingle)&Spline::inteOneDim, *this);
}


double Spline::inteOneDim(const double x) {   throw Bad_Error("indeonedim no support"); return (*this)(x); }

void Spline::inteD(const double x ,const double* y, double* dy) {
  dy[1] = (*this)(x);
}


double Spline::stepIntegrate(double a, double b) {
  unsigned int anf, ende;
  throw Bad_Error("step integrate, not supported");
  double *xd;
  if (own) xd = xdat ; else  xd = mother->xdat;

  splint(xd,a,&anf);
  splint(xd,b,&ende);

  cout << "stepintegrate: " << anf << "  " << ende << endl;

  return stepInt(anf,ende);
}


double Spline::stepInt(int a, int b) {
    throw Bad_Error("step int, not supported ");
  double result = 0;
  double xe = x(a);
   for (int i = a; i < b; i++) {
    double xa  = xe;
    xe = x(i+1);
    double delta = xe-xa;
    result += delta*ydat[i];
   }
   return result;
}


double Spline::average(double a, double b, double tol) { return integrate(a,b,tol) / fabs(b-a); }

void Spline::flip() {
  for (int k = 0; k <= n; k++) {
    int pos = n-k;
    if (pos <= k) break;
  
    if (own) { 
      double swap_x = x(pos);
      xdat[pos] = xdat[k];
      xdat[k] = swap_x;
    }
    double swap_y = y(pos);
    ydat[pos] = ydat[k];
    ydat[k] = swap_y;
  }
}


double Spline::findZero(double a, double b, double tol) {
  SPLINE_DBG_MSG("findZero()" << a << " and " << b << "  ");
  return Miscmath::zbrent( (moSingle) &Spline::fastY, *this, a,b, fabs(a-b)*tol);
}

double Spline::findXForY(const double y, const double a, const double b, const double tol) {
  SPLINE_DBG_MSG("findXForY() " << fastY(a) << " and " << fastY(b) << " ");
  return Miscmath::zbrent( (moSingle) &Spline::fastY, *this, a,b, fabs(a-b)*tol, y);
}

double Spline::findZeroBetween(const int a, const int  b, double tol, const double zl) {
  SPLINE_DBG_MSG("findZeroBetween() " << a << " and " << b << "  ");
  cachedVals().klo = a;
  cachedVals().khi  =b;
  zeroLevel = zl;
  return Miscmath::zbrent( (moSingle) &Spline::inBetweenY, *this, x(a),x(b), fabs(x(a)-x(b))*tol);
}

/*! 
  get a list of zeros or actually crossings of zl  (zerolevel) 
  of this spline, starting from data point #a to data point #b.
  If the y-data at #a coincides with zl,  then this point is skipped. This
  happens with all subsequent points. Thus, if for instance you would
  like to know the zeros of your Spline and say from start() = -3.5 to x-data = -1.7 
  the Spline is zero at the x-data points of this inverval, then the 
  search for zeros will only start from the next x-data points higher than -1.7 on.

  Naturally, in between two data  points, there can only be one zero and hence this
  procedure will never cost you any zeros *except* trailing identical zeros.
*/
  
list<double> *Spline::getZeroList(const double zl,const double tol, int a, int b) {
  //cout << "get Zero List::::: "<<endl;
  SPLINE_DBG_MSG("getZeroList()");
  list<double> *nullen =  new list<double>();
  if (!b) b = n;   // if b has not been given (i.e. is 0), set this to n
  while (y(a) == zl && a < b) a++;
  if (a+1 >= b) return nullen;

  bool oben = (y(a) > zl);
  for (int i = a+1; i <= b; i++) {
    if (oben) { 
      if ( y(i) < zl) {
	oben = false;
	nullen->push_back( findZeroBetween( i-1, i, tol,zl) );
      }
    } else {
      if ( y(i) > zl) {
	oben = true;
	nullen->push_back( findZeroBetween( i-1, i, tol,zl) );
      }
    }
  } 
  //cout << "ende gelaende" << endl;
  return nullen;
}

vector<double> Spline::getZeroVector(const double zl,const double tol, int a, int b) {
  list<double> *lst = getZeroList(zl,tol,a,b);
  vector<double> vec(lst->size());
  int count = 0;
  for (list<double>::iterator i = lst->begin(); i != lst->end(); i++) vec[count++] = *i; 
  delete lst;
  return vec;
}

int Spline::getZeroArray(int a, const int b, double *nullen,const int max,const double tol, const double zl) {
  int k =0;
  //cout <<"GZA: from: " << x(a) << "  " << x(b) << endl;
  while (y(a) == zl && a < b) a++;
  if (a+1 >= b) return 0;
  // cout <<"GZA: from: " << x(a) << "  " << x(b) << endl;

  // cout << "first a: " << y(a);

  bool oben = (y(a) > zl);
  //cout << "  oben: " << oben << endl;
  for (int i = a+1; i <= b; i++) {
    if (oben) { 
      if ( y(i) < zl) {
	oben = false;
	k++; 
	if (k > max) error("getZeroArray","array too small to store zero's");
	nullen[k]= findZeroBetween( i-1, i,  tol,zl);
      }
    } else {
      if ( y(i) > zl) {
	oben = true;
	k++; 
	if (k > max) error("getZeroArray","array too small to store zero's");
	nullen[k]= findZeroBetween( i-1, i,  tol, zl);
      }
    }
  }
  return k;
}


void Spline::getExtrema(list<double>* pos, list<double>* snd , int a,int b,const double tol) {
  Spline&  firstDerivClts = *new Spline(this,"fdclts");
  Spline& sndDerivClts = *new  Spline(this,"sndclts");
 
  bool hasnotbeenarmed = false;

  if (!armed) { arm(); hasnotbeenarmed = true; } 
 
  derive(firstDerivClts);
  firstDerivClts.arm(); 
  firstDerivClts.derive(sndDerivClts);
  sndDerivClts.arm();
 
  list<double> *zeros = firstDerivClts.getZeroList(0.0,tol,a,b);
  for (list<double>::iterator i = zeros->begin(); i != zeros->end(); i++) {
    pos->push_back(*i);
    snd->push_back(sndDerivClts(*i));
  }
 
  delete &firstDerivClts;
  delete &sndDerivClts;
  delete zeros;

  if (hasnotbeenarmed) disarm();
}

void Spline::getMaxima(list<double>* pos, int a,int b,const double tol) {
  list<double> ext, snd;
  getExtrema(&ext,&snd,a,b,tol);
  list<double>::iterator j = snd.begin();
  for (list<double>::iterator i = ext.begin(); i != ext.end(); i++, j++) {
    if ( *j < 0) pos->push_back(*i);
  }
}

void Spline::getMinima(list<double>* pos, int a,int b,const double tol) {
  list<double> ext, snd;
  getExtrema(&ext,&snd,a,b,tol);
  list<double>::iterator j = snd.begin();
  for (list<double>::iterator i = ext.begin(); i != ext.end(); i++, j++) {
    if ( *j > 0) pos->push_back(*i);
  }
}

double Spline::maximum() {
  list<double> max;
  getMaxima(&max);
  double sup = front(); 
  double pos = start();
  for (list<double>::iterator i=max.begin(); i != max.end() ;i++) {
    if (fastY(*i) > sup) { sup = fastY(*i); pos = *i;} 
  }
  if (back() > sup)  pos = stop();
  return pos;
}

double Spline::minimum() {
  list<double> mini;
  getMinima(&mini);
  double inf = front(); 
  double pos = start();
  for (list<double>::iterator i=mini.begin(); i != mini.end() ;i++) {
    if (fastY(*i) < inf) { inf = fastY(*i); pos = *i;} 
  }
  if (back() < inf)  pos = stop();
  return pos;
}


void Spline::merge(Spline &s, const bool ignoredouble) {
  if (s.own) merge(s.xdat, s.n,ignoredouble,s.ydat); 
  else {
    if (!s.mother->valid) error("merge","mother does not have valid x-data");
    if (s.mother->n < n) error("merge","mothers x-data not sufficient for this splines y-data");
    merge(s.mother->xdat, s.n, ignoredouble,s.ydat);
  }
}

/*! 
  if keepSize, then the new xdata und ydata arrays do either have the old size
  maxDim, or if the combined size would be too large, n + max +1.

  if keepSize is false, then the new arrays will always have the smallest possible
  size, however, this may cause doubling soon after. If however, you do not want
  to add anything after merging then keepSize = false is the better choice 
*/
void Spline::merge(const double *X, const int maX, const bool ignoredouble, const double *Y, const double keepSize) {

  throw Bad_Error("temporarily unsopported, due to move to 0 as base index");
  int j = 1;
  int pos = 0;
  if (maX <1) return; // nothing to do
  if (! own) error("merge","not owned");
 
  int need = n + maX + 1;
  if (keepSize) need = max(need , maxDim);

  double *xd = new double[need];
  double *yd = new double[need];

  if (sdat != 0 && need >  maxDim +1 ) {delete[] sdat; sdat =0;}   // if there is some sdat that is too small

  maxDim = need-1;
  for (int i = 0; i <= n; i ++) {    
    if ( j <= maX) {                    // if end of s is not reached
      while (X[j] < x(i)) {           // and as long as s.x is smaller than x,
	pos++;
	xd[pos] = X[j]; 
	if (Y != 0) yd[pos] = Y[j]; else yd[pos] = 0.0;          // merge it in
	j++;
	if (j > maX) break;
      }
    }
    if ( j <= maX) {
      if (X[j] == xdat[i]) {    // if the x-data coincides
	if (ignoredouble) j++;  // if ignore, just increase s count
	else { // if not, we may have to tell there is an error -- if the ydat does not coincide
	  if (Y != 0) {
	    if (Y[j] == ydat[i]) j ++;       // if the y-data also, just skip this point for the other one 
	    else error("merge","data coincides on x but not on y (1)");
	  } else {
	    if (ydat[i] == 0.0) j++; // same as above, if Y is not given then 0 assumed
	    else {
	      cout << "i, j: " << i << "  " << j << endl;
	      cout << "x[i] : " << xdat[i] << "   X[j] = "<<X[j]<<endl;
	      cout << "y[i] : " << ydat[i] << "   Y[j] = 0"<<endl;
	      error("merge","data coincides on x but not on y (2)");
	    }
	  }
	}
      }
    }

    pos++;
    xd[pos] = xdat[i]; yd[pos] = ydat[i];         // if this spline's x ist the smaller one , take this
  }

  delete[] xdat;
  delete[] ydat;

  xdat = xd;
  ydat = yd;
  n = pos;
  
  if (armed) disarm();

  // for (int i = 1; i <= n; i ++) cout << "zus: " << xdat[i] << "  " << ydat[i] << endl;
}



/*!
  Smoothen uses a gaussian window along the spline to smoothen out wiggles.
  It starts k points to the left and k to the right and in total uses m-fold more  points to
  interpolate.

  In other words:  k sets the range over which we average and m sets the number of
  points that are considered in this range.

  It disarms but does not arm, so it is up to you to arm() again
*/

void Spline::smoothen(int k, int m) {

  m = m*k*2+1;
  
  Spline s(*this);
  disarm();  // then, we disarm ourselves
  
  for (int p = 0; p <= n; p ++) {   // for all points of this spline

    int max_right = min(n,p+k);
    int max_left = max(1,p-k); // find the lowest x-value in the range p-k that is still valid (i.e. >= 1)
  
    double max_dist = max(pow(x(max_left) - x(p),2), pow(x(max_right) - x(p),2) ); // get the maximum x-distance ^2
    double sigma2  = max_dist;  // sigma^(-2) such that at max-distance we are down by 1/e
  
    double sum =0,totalweight=0;
    double step = ( x(max_right)-x(max_left) ) / (m-1);  // use m points, i.e. left most, m-2 in the middle and then right most 
     for (double q = x(max_left); q <= x(max_right); q += step ) {
       double dist = q-x(p);
       double weight = exp(- dist*dist*sigma2);
       sum += s(q)*weight;
       totalweight += weight;
       //cout << "p: " << p << "   " << "  q: " << q << endl;
     }
     if (totalweight == 0) error("smoothen","strange");
     ydat[p] = sum / totalweight;
  }  
}
      
   

double Spline::inte(const double x)  {
  throw Bad_Error("inte(x) not supported in spline");
    return (*this)(x);
}



void Spline::explizit() {
  double *x;
  if (own) x = xdat ; else x = mother->xdat;
  for (int i = 0; i <= n;i++) {
    cout << x[i] << " " << ydat[i] << endl;
  }
}
 


/*! 
  This ist the true (well, a bit modified) NR spline interpolation, sorry for the misnomer 
  of CMBFASTS splint which is an integration routine

  If "x" is out of the range xa[1] xa[n], then the boundary values will be given.
  Meaning:: You will not get any errormessages, if you are out of the boundary,
  just the closest available value 

*/
double Spline::splint(double xa[], double x, unsigned int *guess)  {     
  int k;
  CachedValues& cache = cachedVals();
  cache.valid = false;

  double &h = cache.h;
  double &a = cache.a;
  double &b = cache.b;
  double &a3a = cache.a3a;
  double &b3b = cache.b3b;
  double &h26 = cache.h26;
  double &lastX = cache.lastX;
  unsigned int& klo = cache.klo;
  unsigned int& khi = cache.khi;

#ifndef _OPENMP
  if (own) splintcC++; else mother->splintcC++;
#endif

  if (x >= xa[n]) { 
    if (guess) *guess = n-1; 
    return ydat[n]; 
  }
  if (x <= xa[0]) {
    if (guess)  *guess = 0; 
    return ydat[0]; 
  }

  klo=0;
  khi=n;
  while (khi-klo > 1) {
    k=(khi+klo) >> 1;
    if (xa[k] > x) khi=k;
    else klo=k;
  }

  h= xa[khi] -xa[klo];
  if (h == 0.0) {
   //cout << omp_get_thread_num() << ": " << xa << " " << klo << " " << khi << endl;
   error("splint","Spline::splint() Bad xa input to routine splint");
  }

  a=(xa[khi]-x)/h;
  b=(x-xa[klo])/h;
  a3a = a*a*a -a;
  b3b = b*b*b -b;
  h26 = h*h/6.0;
  lastX = x;
  cache.valid = true;

  if (guess != 0)
    *guess = klo;
  return a*ydat[klo]+b*ydat[khi]+(a3a*sdat[klo]+b3b*sdat[khi])*h26;
}

int Spline::indexToLeft(const double x) const
{
  if (x > xdat[n] || x <= xdat[0]) {
    error("indexToLeft", "Spline::indexToLeft () - x outside x-range of spline");
  }
  int lowerIndex=0;
  int upperIndex=n;
  while (upperIndex-lowerIndex > 1) {
    int i =(upperIndex+lowerIndex) >> 1;
    if (xdat[i] == x) return i;
    if (xdat[i] > x) upperIndex=i;
    else lowerIndex=i;
  }
  return lowerIndex;
}

int Spline::indexToRight(const double x) const
{
  if (x > xdat[n] || x <= xdat[0]) {
    error("indexToRight", "Spline::indexToRight () - x outside x-range of spline");
  }
  int lowerIndex=0;
  int upperIndex=n;
  while (upperIndex-lowerIndex > 1) {
    int i =(upperIndex+lowerIndex) >> 1;
    if (xdat[i] == x) return i;
    if (xdat[i] > x) upperIndex=i;
    else lowerIndex=i;
  }
  return upperIndex;
}



/* cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
void Spline::Splder(double *y, double *dy, const int n) {
    /* System generated locals */
    int  i__1;

    /* Local variables */
    static double f[20001];
    static int i, n1;

    if (!splgValid) splini();


/*  Splder fits a cubic spline to y and returns the first derivatives at */
/*  the grid points in dy.  Dy is equivalent to a 4th-order Pade */
/*  difference formula for dy/di. */

    /* Function Body */
    n1 = n - 1;
    if (n1 > 20000) {
      cout << "Spline array overflow!!! n1= " << n1 << "  >20000" <<  endl;
      throw Bad_Error("subroutines::splder() stop");
    }
/*  Quartic fit to dy/di at boundaries, assuming d3y/di3=0. */
    f[0] = (y[0] * -10. + y[1] * 15. - y[2] * 6. + y[3]) / 6.;
    f[n - 1] = (y[n1] * 10. - y[n1-1] * 15. + y[n1 - 2] * 6. - y[n1 - 3]) / 6.;
/*  Solve the tridiagonal system */
/*  dy(i-1)+4*dy(i)+dy(i+1)=3*(y(i+1)-y(i-1)), i=2,3,...,n1, */
/*  with dy(1)=f(1), dy(n)=f(n). */
    i__1 = n1;
    for (i = 2; i <= i__1; ++i) {
	f[i - 1] = splg[i - 1] * ((y[i] - y[i-2]) * 3. - f[i - 2]);
/* L10: */
    }
    dy[n] = f[n - 1];
    for (i = n1; i >= 1; --i) {
	dy[i-1] = f[i - 1] - splg[i-1] * dy[i];
/* L20: */
    }
} /* splder_ */

void Spline::splini() {
    static int i;

    //  Splini must be called before splder to initialize array g in common. 
    splg[0] = 0.;
    for (i = 2; i <= 20001; ++i) splg[i - 1] = 1. / (4. - splg[i - 2]);
    splgValid = true;
} 

void Spline::origSpline(double x[], double y[], int anf, int n, double yp1, double ypn, double y2[])
{
  int i,k;
  double p,qn,sig,un;


  double * u = new double[n+2];
  if (yp1 > 0.99e30)
    y2[anf]=u[anf]=0.0;
  else {
    y2[anf] = -0.5;
    u[anf]=(3.0/(x[anf+1]-x[anf]))*((y[anf+1]-y[anf])/(x[anf+1]-x[anf])-yp1);
  }
  for (i=anf+1;i<n;i++) {
    sig=(x[i]-x[i-1])/(x[i+1]-x[i-1]);
    p=sig*y2[i-1]+2.0;
    y2[i]=(sig-1.0)/p;
    u[i]=(y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
    u[i]=(6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
  }
  if (ypn > 0.99e30)
    qn=un=0.0;
  else {
    qn=0.5;
    un=(3.0/(x[n]-x[n-1]))*(ypn-(y[n]-y[n-1])/(x[n]-x[n-1]));
  }
  y2[n]=(un-qn*u[n-1])/(qn*y2[n-1]+1.0);
  for (k=n-1;k>=anf;k--)
    y2[k]=y2[k]*y2[k+1]+u[k];
  delete[] u;
}


void Spline::error(string r, string s, int x) const {
  cout << endl << "****************** SPLINE ERROR ***************" << endl;
  cout << "Spline::" <<r<<"() " << endl << endl << s;
  if (x != 1234) cout << "  hint: " << x; 
  cout << endl;
  cout << "***********************************************" << endl;
  cout << "Status of spline: " << endl;
  printStatus();
  if (!own) {
    cout << "Status of mother: " << endl;
    mother->printStatus();
  }
  throw Bad_Error("Spline::" + r + "()\n" + s);
}



inline double Spline::operator() (const double x)  {
#if 0
//#ifdef _OPENMP
  omp_set_lock(&mElementAccessLock);
  struct OnScopeExit
  {
    omp_lock_t* mLock;
    OnScopeExit(omp_lock_t* l) : mLock(l) {}
    ~OnScopeExit() { omp_unset_lock(mLock); }
  };
  OnScopeExit releaseLock(&mElementAccessLock);
#endif
#ifndef FORCE_SPLINE_SPEED
 if (!valid)  throw Bad_Error("y() not valid");     
 if (!armed) error("operator()","not armed");     
#endif

 if (own) {
#ifndef FORCE_SPLINE_SPEED
#ifndef _OPENMP
   acc++;
#endif
#endif
   return fastSplint(xdat, x, &cachedVals().guess);
 } else {
#if 0
//#ifdef _OPENMP
  omp_lock_t* myMothersLock = mother->lock();
  OnScopeExit releaseLock(myMothersLock);
#endif
   // if we are valid, so must be mummy

   /*
     Children usually are accessed at around roughly the same x as their mother.
     It is assumed, that mother has already been evaluated and hence useful quantities
     have already been calculated. If mother->cacheValid is true, then
     this function will try to use mothers last klo value as its guess. Usually, this is right.
     It then calculates the interpolation. 
   */
#ifndef FORCE_SPLINE_SPEED
#ifndef _OPENMP
   mother->acc++;
#endif
#endif
   CachedValues& mothersCache = mother->cachedVals();
   if (mothersCache.valid) {
     if (mothersCache.lastX == x) {
#ifndef FORCE_SPLINE_SPEED
#ifndef _OPENMP
       mother->vc++;
#endif
#endif
       unsigned int klo = mothersCache.klo;	
       unsigned int khi = klo+1;
       return mothersCache.a*ydat[klo]+mothersCache.b*ydat[khi]
              +(mothersCache.a3a*sdat[klo]+ mothersCache.b3b*sdat[khi])*mothersCache.h26;
     }  // if we come to this, mother's cache must be wrong...
     mothersCache.valid=false;
   }
#ifndef FORCE_SPLINE_SPEED
#ifndef _OPENMP
   mother->nvc++;
#endif
#endif
   return fastSplint(mother->xdat,x, &mothersCache.guess);
 }
}





void Spline::add(int k, double y) {
  if (k > maxDim)  { printStatus();throw Bad_Error("spline::add not handled yet");}
  ydat[k] += y;
  if (k> n) n=k;
}

void Spline::mul(int k, double x) {
    if (k > n || k < 0) throw Bad_Error("Spline::mul() Element not available");
    ydat[k] *=x;
}

void Spline::kill(Spline *s) {
  for (map< Spline*, Spline*>::iterator i = s->children.begin(); i != s->children.end(); i++) {
    i->second->motherKilled();
    delete i->second;
  }
  delete s;
}


void Spline::generateFamily(vector<Spline*> &v, const unsigned int members, const int k, const string s) {
  if (v.size() != members) v.resize(members);
  v[0] = new Spline(k,s+"_mother");
  //cout << "GENERATE SPLINEFAMILY: " << s << "  at : " << v[1] << endl;
  for (unsigned int j = 1; j < members; j++) {
    stringstream num;
    num << "_" << j;
    v[j] = new Spline(v[0],s+num.str());
    v[j]->generated = j;
  }
  //cout << "I have generated: " << s << " with: " <<members << " members and actual: " << v.size() << endl;

}

void Spline::generateSisters(vector<Spline*> &v, const unsigned int members, const int k, const string s) {
  if (v.size() != members) v.resize(members);
  for (unsigned int j= 0; j < members; j++) v[j] = new Spline(k,s);
}

void Spline::generateChildren(vector<Spline*> &v, Spline*  p, const unsigned int members, const string s) {
  if (v.size() != members) v.resize(members);
  for (unsigned int j= 0; j < members; j++) { v[j] = new Spline(p,s); v[j]->generated = j;}
}

void Spline::generateOne2OneChild(vector<Spline*> &v, vector<Spline*> &p,  const string s) {
  if (v.size() != p.size()) v.resize(p.size());
  for (unsigned int j= 0; j < p.size(); j++) { v[j] = new Spline(p[j],s); v[j]->generated = j;}
}

void Spline::armVector(vector<Spline*> &v) { for (unsigned int i = 0; i < v.size(); i++) v[i]->arm();  }
void Spline::disarmVector(vector<Spline*> &v) { for (unsigned int i = 0; i < v.size(); i++) v[i]->disarm();  }


/*!
  Now, it gets tricky.

  In some applications (as in cmbfast) it is useful to calculate some function S(k,tau) on a grid
  of k and tau values.
  
  In the later calculation, we may need this function on a much denser grid. 

  Now, assume, as we don't have 2-dim Splines, but only 1-dim, that there are 2 vectors<Spline*>
  each holding splines T and K respectively that do a T_k(tau) and K_tau(k) interpolation of S(k,tau).

  The subscript _tau and _k resembles the fact that we deal with vectors T[k](tau).

  
  On calling expandVector( vector of T splines, vector of K splines, int between),
  a new vector is calculated that in between of any k value, the T splines exist (i.e. T_k),
  has "between" new Splines that are also T splines but with k values in between two
  neighbouring T splines.
  
  Now how do we know which k values to take and who will give us the right Datapoints ???

  To answer this, we assume that the Splines K_tau(k) have tau values EXACTLY at
  the tau data-points of T_k(tau), i.e. for all values that exist in T[k]->x(i), there is 
  a spline K[tau]. 

  Put in other words: K[1] is a spline along the k-direction that has fixed value tau = T[anything]->x(1)
  K[2] = T[anything]->x(2) etc. etc.

  If this sounds pretty strange then please observe that in most situation exactly these
  properties arise:

  Example:

  for (k=... )
  for (tau = ...)

  S( k, tau )   = some function

  T[k]->set(tau, S(k,tau))
  K[tau]->set(k, S(k,tau))
  
  That's it: You now have exactly these splines and on calling e.g expandVector(T, K, 3) 
  you get about two to three times the number of T[k] splines back that let you interpolate
  much better for differen k values


  Technical note:

  (1) expandVector() does new a vector<Spline*>, please delelte it after use
  (2) as you will want to know the k value these new splines belong to, call generatedAt()
  (3) vector[1] is mother, all others are children. This is true for the result as well as for the argument
  (4) all xdata is assumed to INCREASE with index i
  (5) the new k - data is EQUALLY SPACED, i.e. now logarithmic spacing in between each interval.
       if however the intervalls have already logarithmic spacing and are not to far apart, then this
       should not matter too much

*/

vector<Spline*>* Spline::expandVectorBetween(vector<Spline *> &v, vector<Spline *> &p,int between) {
 
  throw Bad_Error("expandvector between not supported");

  Spline &first = (*v[1]);
  vector<Spline *> &e = *new vector<Spline *>( (p[1]->size()-1)*(between+1)  + 2);

  e[1] = new Spline(first.size(),"generated."+v[1]->name);
  
  // first, generate the daugther splines
  int c = 2;
  if (v.size() < 3) throw Bad_Error("Spline::expandVectorBetweeen() needs at least vector of three splines");
  for (unsigned int i = 1; i < v.size()-1; i++) {            // we only have intervals until size()-1 and size(), so size()+1 !!!
    for (int j = 1; j <= between; j++) {
      e[c] = new Spline(e[1],"generated."+v[2]->name);   // between the splines, new splines
      c++;
    }
    e[c] = new Spline(e[1],"generated."+v[2]->name);   // this is the boundary spline
    c++; 
  }

  e[1]->gen = p[1]->x(1);  // generated at the beginning :-)
  e[1]->getData(*v[1]);
  // now fill those splines...
  // i assume that the first Spline in vector "p" coresponds to this tau value
  // and i also assume that I will generate q new splines in between the values p[1]->x(i) and p[1]->x(i+1)
  c = 2;
  for (unsigned int k = 1; k < v.size()-1; k++) { 
    double kmin = p[1]->x(k);
    double kmax = p[1]->x(k+1);
    
    double step = (kmax-kmin)/(between + 1);     
    double kk = kmin;
    for (int j = 1; j <= between; j++) {
      kk += step;
      generateSpline(e[c], v,p,kk);
      c++;
    }
    e[c]->gen = kmax;
    e[c]->getData(*v[k+1]);  // just copy the boundary splines that are already there
    c++;   
  }
  e[1]->arm(all);
  return &e;
}


/*!
  Like expandVectorBetween, however it does not generate a fixed number of new
  splines in between each pair of existing splines, but a equally spaced total number
  total between to left and rightmost existing spline 
*/
vector<Spline*>* Spline::expandVector(vector<Spline *> &v, vector<Spline *> &p, int total) { 
  Spline &first = (*v[1]);
  vector<Spline *> &e = *new vector<Spline *>(total+1);
  throw Bad_Error("expandvector  not supported");
  generateFamily(e, total, first.size(), "generated." +v[1]->name);

  e[1]->gen = p[1]->start();  // generated at the beginning :-)
  e[total]->gen = p[1]->stop();
  e[1]->getData(*v[1]);
  e[total]->getData(*v[v.size()-1]);

  double step = (p[1]->stop()- p[1]->start() )/(total - 1);     
  double kk =  p[1]->start();
  for (int k = 2; k < total; k++) {
    kk += step;
    generateSpline(e[k], v, p, kk); 
  }
  e[1]->arm(all);
  return &e;
}

/*!
  Give a vector<Spline*> v with splines T interpolating along some tau-axis
  and vector<Spline*> p of splines K interpolating along some k-axis,
    
  fill an empty spline e with data points at tau values that are at the x-data of
  T containing values interpolated using the K splines.

  It is assumed that the first k-spline coressponds to the first tau x-data value
  of the T-Splines. 

  See comments at expandVectorBetween()
*/

void Spline::generateSpline(Spline* e, vector<Spline *> &v, vector<Spline *> &p, const double k, bool linear) {  
  e->gen = k;  

    for (int i = 0; i < v[0]->size(); i++) {          // go through all data points of vectors first spline
      // now i assume that the first Spline in vector "p" coresponds to this tau value         
      e->setForce(v[0]->x(i), (*p[i])(k));
    }
}


/*!
  Same as generateSpline() however, it assumes that p[0] is
  this spline and that this spline is the  mother of all other p[j] Splines. 

  Hence, it can efficiently cache the access and speed things up
*/

void Spline::generateSplineXXL(Spline* e, vector<Spline *> &v, vector<Spline *> &p, const double k) {
  e->gen = k;
  CachedValues& cache = cachedVals();
  e->setForce(v[0]->x(0),  fastSplint(xdat,k, &cache.guess) );  // this fills the cache variables
  // cout << "xxl2 " << v[1]->size() << endl;
  if (cache.valid) {
    double& h = cache.h;
    double& a = cache.a;
    double& b = cache.b;
    double& a3a = cache.a3a;
    double& b3b = cache.b3b;
    double& h26 = cache.h26;
    double& lastX = cache.lastX;
    unsigned int& khi = cache.khi;
    unsigned int& klo = cache.klo;
    double y;
    Spline *s;
    for (int i = 1; i < v[0]->size(); i++) {          // go through all data points of vectors first spline
      // now i assume that the first Spline in vector "p" coresponds to this tau value
      // also, we take this splines cache value ...
      s = p[i];
      y = (a*s->ydat[klo]+b*s->ydat[khi]+(a3a*s->sdat[klo]+ b3b*s->sdat[khi])*h26);
      e->setForce(v[0]->x(i), y);
    }
  } else {
    for (int i = 1; i < v[0]->size(); i++)   e->setForce(v[0]->x(i), (*p[i])(k));
  }
}

/*!
  Easy to use convenience function. For example to generate 
  Transferfunctions at some specific time and output it to a file

  Input: vector<Spline*> &v  must be a vector that stores function values of the
  same k-position as the yet to be generated spline.

  vector<Spline*>&p has the same function values but knows the corresponding
  tau values. 

  The variable k in the function argument is then the TIME (not k :-) at which
  you would like to have a slice through this 2-D plot, so to speak. The direction
  is of course the k-direction of the plot, i.e. the x-values of the v splines.
  
  If you need examples, look at cosmos.cc
*/
void Spline::generateAndDump(vector<Spline *> &v, vector<Spline *> &p,string file,const double k) {
  Spline * tmp = new Spline(v.size(),"Spline::generateAndDump");   
  generateSpline(tmp, v, p, k);     

  tmp->arm();
  tmp->dump(file); 
  delete tmp;
}
/*! 
  If this Spline has been automatically generated via expandVector() 
  then you may very well be interested which value of "kk" of expandVector()
  belongs to this specific spline. In other words and in the language of 
  expandVector():

  We generated splines which are interpolating along the tau direction by
  taking "kk" values along the k direction in between a given spline-vector interpolating
  along the k direction at the tau values of the splines T. 

  We store thes kk values so you know to which value of k this spline belongs. 
*/

double Spline::generatedAt()  const {
  return gen;
}


// bool Spline::splgValid = false;
double Spline::splg[20010];

vector<Spline*> Spline::convolutionVektor;
unsigned int Spline::nvar;
double Spline::zeroLevel(0.0);
Spline* Spline::partner;
double Spline::shift;
double Spline::faktor;
double* Spline::convXdat;
Spline* Spline::convTwin;

double Spline::fitNeeds(const pair<double,double> c) {
  double y = fastY(c.first);
  if (y == 0.0) error("fitNeeds","Function value is zero, thus infinity-scaling would be needed to fit");
  return c.second/y;
}

void Spline::fitToSingle(const pair<double,double> c) {  *this *= fitNeeds(c); }

int Spline::WatchCount=0;
map<string,int> Spline::WatchMap;
void Spline::printWatchMap() {
  for (map<string,int>::iterator i = WatchMap.begin(); i != WatchMap.end();i++) {
    cout << i->first << "     " << i->second << endl;
  }
}






//
// Here come the splineweb class functions
//

SplineWeb::SplineWeb(string n,Anchor *a,const int x, const int y) : AnchorEnabled(a),X(x), Y(y) , name(n),  igrid(-1), jgrid(-1) {
  Spline::generateFamily(X, x, y, n+"_x");
  Spline::generateFamily(Y, y, x, n+"_y");
}

SplineWeb::~SplineWeb() {
  // no need to do anything, CleanVectors will look after themselves
}

void SplineWeb::printStatus(std::ostream& os) {
  os << "-------------- Status of SplineWeb " << name << "------------------\n";
  os << "points in x-direction: " << xSize() << endl;
  os << "points in y-direction: " << ySize() << endl;
  os << "-------------------------------------------------------------------\n";
  os << "Splines in x-direction: \n";
  for (int i = 0; i <= xSize(); ++i) {
    Spline* s = X[i];
    os << "Spline no " << i << ": size " << s->size() << ", x[0]=" << s->start() << ", x[n]=" << s->stop()
       << " y[0]=" << s->front() << ", y[n]=" << s->back() << "\n";
    if (xSize()>15 && i > 5 && i < xSize()-6) {
      cout << "                 ...\n";
      i = xSize()-6;
    }
  }
  os << "-------------------------------------------------------------------\n";
  os << "Splines in y-direction: \n";
  for (int i = 0; i <= ySize(); ++i) {
    Spline* s = Y[i];
    os << "Spline no " << i << ": size " << s->size() << ", x[0]=" << s->start() << ", x[n]=" << s->stop()
       << " y[0]=" << s->front() << ", y[n]=" << s->back() << "\n";
    if (ySize()>15 && i > 5 && i < ySize()-6) {
      cout << "                 ...\n";
      i = ySize()-6;
    }
  }
  os << "-------------------------------------------------------------------" << endl;
}


void SplineWeb::set(const double x, const double y, const double z) {
  map<double,unsigned int>::iterator i,j;
  unsigned int xpos, ypos;

  i = xgrid.find(x);
  j = ygrid.find(y);

  if (i == xgrid.end()) {
    xgrid[x] = ++igrid;
    xpos = igrid;
  } else {
    xpos = xgrid[x];
  }
  if (j == ygrid.end()) {
    ygrid[y] = ++jgrid; ypos = jgrid;
  } else {
    ypos = ygrid[y];
  }

  if (xpos == X.size()) {
    // cout << "SPLINEWEB: " << name << "  resizing x : " << X.size()*2 << endl;
    // cout << "I resize X-Vector: entry one: " << X[1] << endl;
    X.resize(X.size()*2);
    string name_x = name + "_x";
    for (unsigned int i = xpos; i < X.size();i++) {
      X[i] = new Spline(X[0],name_x);
    }
  }
  if (ypos == Y.size()) {
    //cout << "SPLINEWEB: " << name << "  resizing y : " << Y.size()*2 << endl;
    //cout << "I resize Y-Vector: entry one: " << Y[1] << endl;
    Y.resize(Y.size()*2);
    string name_y = name + "_y";
    for (unsigned int i = ypos; i < Y.size();i++) {
      Y[i] = new Spline(Y[0],name_y);
    }
    //cout << "resized y up to: " << count << endl;
  }
  //  cout << xpos<< " : " << ypos << endl;

  X[xpos]->setForce(y,z);
  Y[ypos]->setForce(x,z);
  //  cout <<"done"<<endl;
}


void SplineWeb::makeProper()
{
  if (xSize() <= 0 || ySize() <= 0)
    return;
  if (X[0]->isArmed() || Y[0]->isArmed()) {
    throw Bad_Error("SplineWeb::makeProper() - SplineWeb is already armed(); no need to call makeProper().");
  }

  Spline* oldXMother = X[0];
  Spline* oldYMother = Y[0];

  oldXMother->makeProper();
  oldYMother->makeProper();

  CleanVector<Spline*> tmpX(X.size()), tmpY(Y.size());
  map<double,unsigned int>::iterator it, end;
  it=xgrid.begin();
  end=xgrid.end();
  unsigned int xpos=0;
  for( ; it!=end; ++it) {
    tmpX[xpos] = X[it->second];
    it->second = xpos++;
  }
  it=ygrid.begin();
  end=ygrid.end();
  unsigned int ypos=0;
  for( ; it!=end; ++it) {
    tmpY[ypos] = Y[it->second];
    it->second = ypos++;
  }
  X = tmpX;
  Y = tmpY;
  fill(tmpX.begin(), tmpX.end(), (Spline*)0);
  fill(tmpY.begin(), tmpY.end(), (Spline*)0);

  vector<Spline*>::iterator newXposMother = find(X.begin(), X.end(), oldXMother);
  vector<Spline*>::iterator newYposMother = find(Y.begin(), Y.end(), oldYMother);

  swap((*X.begin())->ydat, oldXMother->ydat);
  swap(X[0], *newXposMother);
  swap((*Y.begin())->ydat, oldYMother->ydat);
  swap(Y[0], *newYposMother);

  for (xpos=igrid+1; xpos<X.size(); ++xpos) {
    string name_x = name + "_x";
    //stringstream names; names << xpos; name_x+names.str();
    X[xpos] = new Spline(X[0], name_x);
  }
  for (ypos=jgrid+1; ypos<Y.size(); ++ypos) {
    string name_y = name + "_y";
    //stringstream names; names << ypos; namey+=names.str();
    Y[ypos] = new Spline(Y[0], name_y);
  }
}

void SplineWeb::shrinkToFit() {
  X.resize(igrid+1);
  Y.resize(jgrid+1);
}


void SplineWeb::arm() {
  disarm();
  Spline::armVector(X);
  Spline::armVector(Y); 
}

void SplineWeb::disarm() {
  if (X[0]->isArmed()) Spline::disarmVector(X); 
  if (Y[0]->isArmed()) Spline::disarmVector(Y);
}

/*! 
  Dump a spline that hast const y-position y and the direction x
*/
void SplineWeb::dumpAlongX(const double y, string file) {
   if (! X[0]->isArmed()) { 
     cout << "going to arm X: " << X.size();
    X.resize(igrid+1); 
    cout << " after resize: " << X.size() << endl;
    Spline::armVector(X); }  // this also kills all splines that are not used
  if (! Y[0]->isArmed()) { 
    cout << "going to arm Y: " << Y.size();
    Y.resize(jgrid+1); 
    cout << " after resize: " << Y.size() << endl;
    Spline::armVector(Y); } 
  Spline::generateAndDump(Y, X, file, y);
}

/*! 
  Dump a spline that hast const x-position x and the direction y
*/
void SplineWeb::dumpAlongY(const double x, string file) {
  if (! X[0]->isArmed()) { 
    //cout << "going to arm: " << X.size();
    X.resize(igrid+1); 
    //cout << " after resize: " << X.size() << endl;
    Spline::armVector(X); }  // this also kills all splines that are not used
  if (! Y[0]->isArmed()) { 
    //cout << "going to arm: " << Y.size();
    Y.resize(jgrid+1); 
    //cout << " after resize: " << Y.size() << endl;
    Spline::armVector(Y); }
  Spline::generateAndDump(X, Y, file, x);
}

/*! 
  Create a spline that hast const  y-position y and the direction x
  return pointer to this new spline. It is your responsibility to delete 
  the spline, if you do not use it anymore... (However, specifyig an anchor will make this easier for you)
    If mother is null the new spline will own its data, if not, mother  
    (or if itself has a mother then this of course) 
    will  be the mother to this spline (i.e. hold the x-data) 
*/

Spline* SplineWeb::createAlongX(const double y, string nombre,Anchor* a, Spline* mother) {
  try {
    if (! X[0]->isArmed()) {
      // this also kills all splines that are not used
      X.resize(igrid+1); Spline::armVector(X);
    }
    if (! Y[0]->isArmed()) {
      Y.resize(jgrid+1); Spline::armVector(Y);
    }
    Spline *neu;
    if (mother) {
      neu = new Spline(mother,nombre,a); 
    } else  {
      neu = new Spline(Y[0]->size(), nombre,a);
    }
    X[0]->generateSplineXXL(neu,Y, X, y);
    return neu;
  } catch (...) {
    cout << "Error occured in SplineWeb::createAlongX() for SplineWeb '" << name << "'\n";
    printStatus();
    throw;
  }
  return 0;
}

/*! 
  Create a spline that hast const x-position x and the direction y
  return pointer to this new spline. It is your responsibility to delete 
  the spline, if you do not use it anymore...
  If mother is null the new spline will own its data, if not, 
  mother  (or if itself has a mother then this of course)  will
  be the mother to this spline (i.e. hold the x-data) 
*/
Spline* SplineWeb::createAlongY(const double x, string nombre,Anchor *a, Spline *mother) {
  //X   cout << X[0]->isArmed()
  //X     << Y[0]->isArmed() << " " << "Y[0]->sdat " << Y[0]->sdat<<  " is status" << endl;

  try {
    if (! X[0]->isArmed()) { 
      // this also kills all splines that are not used
      X.resize(igrid+1);
      Spline::armVector(X);
    }
    if (! Y[0]->isArmed()) {
      Y.resize(jgrid+1);
      Spline::armVector(Y);
    }
    Spline *neu;
    if (mother) {
      neu = new Spline(mother,nombre,a); 
    } else {
      neu = new Spline(X[0]->size(), nombre,a);
    }
    Y[0]->generateSplineXXL(neu,X, Y, x);
    return neu;
  } catch (...) {
    cout << "Error occured in SplineWeb::createAlongY() for SplineWeb '" << name << "'\n";
    printStatus();
    throw;
  }
  return 0;
}

void SplineWeb::fitToSingleAlongX(const double y, const pair<double,double> c) {
    if (! X[0]->isArmed()) { X.resize(igrid+1); Spline::armVector(X); }  // this also kills all splines that are not used
    if (! Y[0]->isArmed()) { Y.resize(jgrid+1); Spline::armVector(Y); }
    Spline* neu = new Spline(Y.size(), "fitToSingleX");
    Spline::generateSpline(neu,Y, X, y);
    neu->arm();
    *this *= neu->fitNeeds(c);
    //neu->dump("neuvorher");
    neu->fitToSingle(c);
    //neu->dump("neunachher");
    Spline::disarmVector(X);
    Spline::disarmVector(Y);
    delete neu;
}


void SplineWeb::fitToSingleAlongY(const double x, const pair<double,double> c) {
    if (! X[0]->isArmed()) { X.resize(igrid+1); Spline::armVector(X); }  // this also kills all splines that are not used
    if (! Y[0]->isArmed()) { Y.resize(jgrid+1); Spline::armVector(Y); }
    Spline* neu = new Spline(X.size(), "fitToSingleY");
    Spline::generateSpline(neu,X, Y, x);
    neu->arm();
    *this *= neu->fitNeeds(c);
    Spline::disarmVector(X);
    Spline::disarmVector(Y);
    delete neu;
}

/*! 
  Expensive walk-through the xgrid map. It just counts from
  1...i and iterates through xgrid until it finds it.
*/
double SplineWeb::x(int i) {
  map<double,unsigned int>::iterator j=xgrid.begin();
  for (int k=0; k < i; k++) { 
    j++;
    if (j == xgrid.end()) throw Bad_Error("SplineWeb::x() index out of range");
  }
  if (j == xgrid.end()) throw Bad_Error("SplineWeb::x() index out of range");
  return j->first;
}

double SplineWeb::y(int i) {
  map<double,unsigned int>::iterator j=ygrid.begin();
  for (int k=0; k < i; k++) { 
    j++;
    if (j == ygrid.end()) throw Bad_Error("SplineWeb::y() index out of range");
  }
  if (j == ygrid.end()) throw Bad_Error("SplineWeb::y() index out of range");
  return j->first;
}


/*!
  This gets the maximum of the spline web. 
*/
pair<double,double> SplineWeb::maximum_internal(double a, double b) {
  
  bool virgin=true;
  double sup=0, pos_x=0, pos_y=0;
 
  double step = (b-a)*2e-1 / xgrid.size(); // if equally spaced, five steps per interval
  for (double x = a; x < b; x +=step) {
      Spline *s = createAlongY(x,"maximum");
      s->arm();
      double p =  s->maximum();
      double y = s->fastY(p);
      if (y > sup || virgin) { sup = y; pos_y = p; pos_x=x;virgin = false;} 
      delete s;
  };
  
  if (step < 1e-6*(Y[0]->stop()-Y[0]->start())) return make_pair(pos_x,pos_y);
  
  a = max(pos_x - step,Y[0]->start());  // lowest x
  b = min(pos_x + step,Y[0]->stop()); // highest x
  // recursivly call this function
  return maximum_internal(a,b);
}

pair<double,double> SplineWeb::maximum() {
  return maximum_internal(Y[0]->start(),Y[0]->stop());
}


/*!
  Lays a raster of dimension nx times ny over the Spline Web. If you specify xmin etc, it uses the patch
  enclosed by the coordinates you give. Otherwise, the whole Spline is covered.
  Returns a vector< vector<double> > * that is yours (also yours to delete !). 
  It runs from 0 ... nx-1 and 0.. ny-1. 

  totalInt will be filled with the Integral over the patch. 


vector< vector<double> > * SplineWeb::rasterize(int nx, int ny, double &totalInt, double xmin, double xmax, double ymin, double ymax) {
  if (xmin == WEB_RASTER_MAGIC) xmin = x(0);
  if (ymin == WEB_RASTER_MAGIC) ymin = y(0);
  if (xmax == WEB_RASTER_MAGIC) xmax = x(igrid);
  if (ymax == WEB_RASTER_MAGIC) ymax = y(jgrid);
  
  totalInt = 0;
  
  double dx =  (xmax - xmin)/(double) (nx -1);
  double dy = (ymax - ymin)/(double)(ny-1);
  // get 2d-array and resize it appropriatly
  vector< vector<double> > * v = new vector< vector<double> >(nx);
  for (int i = 0; i <  nx; i++) (*v)[i].resize(ny);
  
  for (int i = 0; i < nx ; i++) {
    double xpos = xmin + i*dx;
    // at each xpos, we create a spline pointing in y - direction
    Spline *s = createAlongY(xpos , "rasterize");
    s->arm();
    for (int j =0; j< ny; j++) {
      // and at each ypos, we evaluate this spline 
      double ypos =  ymin + j*dy;
      double value = s->fastY(ypos);
      (*v)[i][j] = value;
      totalInt += value;
    }
    delete s;
  }
  totalInt *= dx*dy;
  return v;
} 


  Given a vector[double][double] v, and a double frac, this function
  returns a list of pair<int,int> which should be interpreted as points (x,y).
  Each two pairs in the list define start and end point of a line.
  If you draw these lines, then the total area covered will have the integrated
  fractional weight frac of the total integral over this raster.

  In other words: If you have a SplineWeb, you may call rasterize() and
  using the vector[][] returned, fractionOfRaster() will:
  1. Integrate the total vector[][]
  2. Set a threshold value such that the sum of all vector[][] entries that
  has a value > threshold makes up fraction frac of the total integral
  3. In principle, it could give back all the indices [][] of these points,
  however, it is of advantage to give it back in junks, i.e. lines, cause
  usually, several points in a row will be selected. This is why you get pairs
  of points(x,y) that are the beginning and the end of a line


  
list < pair<int,int> >* SplineWeb::fractionOfRaster(vector< vector<double> > *v, double frac) {
  int nx = (*v).size();
  int ny = (*v)[0].size();
  
  double totalSum=0;
  double min = (*v)[0][0];
  double max = (*v)[0][0];

  for (int i = 0; i < nx; i++) {
    for (int j = 0; j < ny; j++) {
      double value = (*v)[i][j];
      totalSum += value;
      if (value < min) min = value;
      if (value > max) max = value;
    }
  }

  if (totalSum == 0) throw fractionOfRasterFailed("Integrated propability was 0");

  double threshold = min + 0.5*(max - min);
  double step = 0.25*(max - min);

  for (;;) {
    double sum=0;
    for (int i = 0; i < nx; i++) {
      for (int j = 0; j < ny; j++) {
	double value = (*v)[i][j];	
	if (value > threshold) sum+= value;      
      }
    }
    // we break now, if the relative accuracy is better than 1 %
    if ( fabs ( (sum - totalSum*frac)/(totalSum*frac))  < 1e-3 || step < (max -min) *1e-20) break; 
    if (sum > totalSum*frac) threshold += step;
    if (sum < totalSum*frac) threshold -= step;
    step *= 0.5;
    //cout << "to: " << totalSum << " sum: " << sum << "  min: " << min << "   max: " << max << " thres: " << threshold << "   step size: " << step << endl;
  }


  list < pair<int,int> > * lst = new list <pair<int,int> >;

  for (int i = 0; i < nx; i++) {
    bool start = true, need_end=false;
    int links=0;  // initialization to silence compiler
    for (int j = 0; j < ny; j++) {
      double value = (*v)[i][j];	
      if (value > threshold) {
	if (start) {	  
	  lst->push_back(make_pair(i,j)); 
	  start = false; 
	  links = j;
	  need_end=true;
	} else {
	  if (j == links + 1) {  // naechster Nachbar 
	    links = j;
	  }  else {  //also es gibt start und wir haben einen sprung
	    lst->push_back(make_pair(i,links));  // end point of the last line
	    lst->push_back(make_pair(i,j));   // beginning of the new line
	    need_end=true;
	    links = j;
	  }
	}
      }
    }
    if (need_end) {  // then we necessarily need an end point
      lst->push_back(make_pair(i,links));
    }  
  }
  return lst;
}
*/

#if 0
//----------SplineFamilyIntegrator--------------
SplineFamilyIntegrator::SplineFamilyIntegrator(Spline* mother, double tol)
      :mTolerance(tol), mMother(mother)
{
      mSplines.push_back(mother);
      map<Spline*,Spline*>::iterator i = mother->children.begin();
      map<Spline*,Spline*>::iterator end = mother->children.end();
      for(; i != end; ++i) {
        mSplines.push_back(i->second);
      }
}


double SplineFamilyIntegrator::integrate()
{
  const int size=mSplines.size()+1;
  double* y=new double[size];
  std::fill(y, y+size, 0.);
  double initialGuess = 5e-3;

  Miscmath::rungeInt(y, size-1, mMother->start(), mMother->stop(), mTolerance,
                     (mMother->start()-mMother->stop())*initialGuess, 0,
                     (moDerivs)&SplineFamilyIntegrator::integrationHelper, *this);
  mResult.clear();
  std::copy(y+1, y+size, std::back_inserter(mResult));
//X   mResult = std::vector<double>();
//X  cout << mSplines.size()<< "here        ------------ " << flush;
//X   for (int j = 1; j <= mSplines.size(); ++j, cout << " " << j << endl) {
//X     cout << "going to insert: " << j << " - " << endl;
//X     cout << "that is: " << y[j] << " - " << endl;
//X     mResult.push_back(y[j]);
//X   }
//X   cout << "--------done----------------" << endl;
}

void SplineFamilyIntegrator::integrationHelper(const double x, const double* y, double* dy)
{
  const int size = mSplines.size();
  for (int i = 0; i < size; ++i)
   dy[i+1] = mSplines[i]->fastY(x);
}

double SplineFamilyIntegrator::result(const Spline *const s) const
{
  const double size = mSplines.size();
  for (int index=0; index < size; ++index)
    if (mSplines[index] == s)
      return mResult[index];

  throw Bad_Error("SplineFamilyIntegrator::result - no such Spline in the family.");
}
#endif


