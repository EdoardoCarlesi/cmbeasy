#include "allskylensing.h"

#include "cosmos.h"
#include "controlpanel.h"
#include "cl.h"

#define MAX_ANGLE ( M_PI / 32 )
#define SUM_STEP_SIZE 1
#define INTEGRATION_POINTS ( 200 )

/*
static inline double ctnorm( double l ) { return l * (l - 1) * (l + 1) * (l + 2); }     // conventional l factors
*/

// The only class-static member of AllSkyLensing
const AllSkyLensing::FactorStruct AllSkyLensing::factorial[12] =
{
  { 0,  1.0 },
  { 1,  1.0 },
  { 2,  2.0 },
  { 3,  6.0 },
  { 4,  24.0 },
  { 5,  120.0 },
  { 6,  720.0 },
  { 7,  5040.0 },
  { 8,  40320.0 },
  { 9,  362880.0 },
  { 10, 3628800.0 },
  { 11, 39916800.0 }
};


AllSkyLensing::AllSkyLensing( Cosmos& cosmos, const ControlPanel& controlPanel, const CL& cl,
           const CmbCalc& cmbcalc )
         : mCosmos( cosmos ), mCmbCalc( cmbcalc ), mControlPanel( controlPanel ),
           mCl( cl )
{
  if ( controlPanel.lensing == ControlPanel::nonelinear)
    throw Bad_Error( "Nonlinear corrections are not yet implemented for AllSkyLensing" );
}

CL * AllSkyLensing::lensedCls()
{
  DlSplines::initializeFactorTables( (int)mCl.ts[ 0 ]->stop() );

  Anchor anchor; Anchor *a = &anchor;
  CL * lensed = new CL;  // a new set of Cl's namely the lensed ones
  // Create a whole family of those splines, ts[0] being the mother
  Spline::generateFamily(lensed->ts,mCosmos.InitialPower.size(), mCl.ts[0]->size(), "clts");
  Spline::generateChildren(lensed->es,lensed->ts[0], mCosmos.InitialPower.size(), "lensedes");
  Spline::generateChildren(lensed->cs,lensed->ts[0], mCosmos.InitialPower.size(), "lensedcs");
  Spline::generateChildren(lensed->bs,lensed->ts[0], mCosmos.InitialPower.size(), "lensedbs");
  Spline::generateChildren(lensed->kk,lensed->ts[0], mCosmos.InitialPower.size(), "lensedkk");
  Spline::generateChildren(lensed->tk,lensed->ts[0], mCosmos.InitialPower.size(), "lensedtk");
  Spline::generateChildren(lensed->tt,lensed->ts[0], mCosmos.InitialPower.size(), "lensedtt");
  Spline::generateChildren(lensed->et,lensed->ts[0], mCosmos.InitialPower.size(), "lensedet");
  Spline::generateChildren(lensed->bt,lensed->ts[0], mCosmos.InitialPower.size(), "lensedbt");
  Spline::generateChildren(lensed->ct,lensed->ts[0], mCosmos.InitialPower.size(), "lensedct");
  

  for ( unsigned int n = 0; n < mCosmos.InitialPower.size(); ++n )
  {
    Spline *c =  mCosmos.createLensingPowerLimber( n, "mCpsi_l", 2, (int)mCl.ts[ 0 ]->stop(), a);
    c->arm();
    mCpsi_l.push_back( c );

    LensingDifference * diffs = new LensingDifference( mCl.ts[ n ]->size() );
    for ( double theta = 0; theta <= MAX_ANGLE; theta += MAX_ANGLE/INTEGRATION_POINTS )
    {
      oneAngle( n, theta, diffs );
    }

    /*
    const int size = mCl.ts[ n ]->size();
    Spline *tsLensed = new Spline( size, "tsLensedSpline" );
    Spline *csLensed = new Spline( size, "csLensedSpline" );
    Spline *esLensed = new Spline( size, "esLensedSpline" );
    Spline *bsLensed = new Spline( size, "bsLensedSpline" );
    */
    for ( int i = mCl.ts[ n ]->first(); i <= mCl.ts[ n ]->last(); ++i )
    {
      const double l = mCl.ts[ n ]->x( i );
      lensed->ts[n]->set( l, mCl.ts[ n ]->y( i ) );
      lensed->ts[n]->add( i, 2 * M_PI * norml2( l ) * diffs->ts->y( i ) );
      lensed->cs[n]->set(mCl.cs[ n ]->y( i ) );
      lensed->cs[n]->add( i, 2 * M_PI * ( l*( l+1 ) ) * diffs->cs->y( i ) ); //sqrt( ctnorm( l ) ) * norml2( l ) ?
      lensed->es[n]->set( mCl.es[ n ]->y( i ) );  // ctnorm( l )*norml2( l ) ?
      lensed->es[n]->add( i,  2 * M_PI * l * ( l+1 ) * diffs->e->y( i ) );
      lensed->bs[n]->set( 2 * M_PI * l * ( l+1 ) * diffs->b->y( i ) );

      // do nothing for the rest
      lensed->tt[n]->set(mCl.tt[n]->y(i));
      lensed->et[n]->set(mCl.et[n]->y(i));
      lensed->bt[n]->set(mCl.bt[n]->y(i));
      lensed->ct[n]->set(mCl.ct[n]->y(i));
      // even copy legacy stuff
      lensed->kk[n]->set(mCl.kk[n]->y(i));
      lensed->tk[n]->set(mCl.tk[n]->y(i));
    }
    delete diffs;
    
    /* Feature came close to bug :-)
    tsLensed->arm();
    csLensed->arm();
    bsLensed->arm();
    esLensed->arm();
  
    lensed->ts.push_back( tsLensed );
    lensed->cs.push_back( csLensed );
    lensed->es.push_back( esLensed );
    lensed->bs.push_back( bsLensed );
    */
  }
  return lensed;
}

void AllSkyLensing::dumpLensedCls( const string fileName )
{
  CL * lensed = lensedCls();
  lensed->ts[0]->arm(Spline::all); 
  mCmbCalc.dumpCl( mCosmos.InitialPower, *lensed, mControlPanel, fileName, std::string(), true );
  delete lensed;
}

void AllSkyLensing::oneAngle( int n, double theta, LensingDifference * diffs )
{
  Spline * const cltsn = mCl.ts[ n ];
  const int lmax = ( int ) cltsn->stop();

  static double lastTheta = 0.;
  const double step = sin( theta ) *(- lastTheta + theta ); // integrals are over dcostheta
  lastTheta = theta;

  Anchor a;

  DlSplines dls( lmax, theta, &a );
  const double s2 = sigma2( n, theta, dls );
  const double cgl2 = C_gl2( n, theta, dls );

  const Ximn ximns( theta, s2, cgl2, ( int )cltsn->start(), lmax );

  const double tsContribution = step * xiLensed( n, theta, cgl2, dls, ximns );
  const double csContribution = step * xiLensed_X( n, theta, cgl2, dls, ximns );
  const double xip = xiLensed_plus( n, theta, cgl2, dls, ximns );
  const double xim = xiLensed_minus( n, theta, cgl2, dls, ximns );


  const int last = cltsn->last();
  for ( int i = cltsn->first(); i <= last; ++i )
  {
    const double l = cltsn->x( i );
    const int intl = ( int ) l;

    diffs->ts->add( i, tsContribution * dls.dl00->y(intl) );
    diffs->cs->add( i, csContribution * dls.dl20->y(intl) );
    diffs->e->add( i, 0.5 * ( xip * dls.dl22->y( intl ) + xim * dls.dl2m2->y(intl) ) * step );
    diffs->b->add( i, 0.5 * ( xip * dls.dl22->y( intl ) - xim * dls.dl2m2->y(intl) ) * step );
  }
}

double AllSkyLensing::sigma2( int n, double beta, const DlSplines& s )
{
  Spline * const cpsi_l = mCpsi_l[n];
  const double stop = cpsi_l->stop();
  double value = 0;
  for ( double l = cpsi_l->start(); l <= stop; ++l )
  {
    value += ( ( 2.*l + 1. )*l*( l+1. )/pow(l, 4.) )*cpsi_l->fastY( l )*( 1. - s.dl11->y((int)l) );
  }
  return value/(4.*M_PI );
}

double AllSkyLensing::C_gl2( int n, double beta, const DlSplines& s )
{
  Spline * const cpsi_l = mCpsi_l[n];
  double value = 0;

  const double stop = cpsi_l->stop();
  for ( double l = cpsi_l->start(); l <= stop; ++l )
  {
    value += ( ( 2.*l + 1. )*l*( l+1. )/pow( l, 4. ) )*cpsi_l->fastY( l )*s.dl1m1->y((int) l );
  }
  return value/( 4 * M_PI );
}

double AllSkyLensing::xiLensed( int n, double beta, double cgl2, const DlSplines& s, const Ximn& x )
{
  const int lmax = ( int )mCl.ts[ n ]->stop();
  double retval = 0;
  const double * X_000 = &(x.X_000[0]);
  const double * Xprime_000 = &(x.Xprime_000[0]);
  const double * X_220 = &(x.X_220[0]);

  for (double l = mCl.ts[ n ]->start(); l <= lmax; l += SUM_STEP_SIZE )
  {
    double sum = ( pow( *X_000,2 ) - 1. )*s.dl00->y(( int )l);
    sum += 8./(l*(l+1.))*cgl2* (*Xprime_000) * (*Xprime_000 ) * s.dl1m1->y(( int )l);
    sum += pow(cgl2,2)*((*Xprime_000*( *Xprime_000 ))*s.dl00->y(( int )l)+pow( *X_220,2)*s.dl2m2->y(( int )l));
    retval += SUM_STEP_SIZE * (2.*l+1.) * mCl.ts[ n ]->fastY(l) * sum / norml2(l);
    ++X_000; ++Xprime_000; ++X_220;
  }
  return retval/( 4*M_PI );
}


double AllSkyLensing::xiLensed_plus( int n, double beta, double cgl2, const DlSplines& s, const Ximn& x )
{
  const double * X_022 = &(x.X_022[0]);
  const double * X_220 = &(x.X_220[0]);
  const double * X_132 = &(x.X_132[0]);
  const double * X_121 = &(x.X_121[0]);
  const double * X_242 = &(x.X_242[0]);
  const double * Xprime_022 = &(x.Xprime_022[0]);

  const int lmax = ( int )mCl.ts[ n ]->stop();
  double retval = 0;
  for ( int l = ( int ) mCl.ts[ n ]->start(); l <= lmax; l += SUM_STEP_SIZE )
  {
    double sum = 0;
    sum = ( pow( *X_022, 2 )-1. )*s.dl22->y(l);
    sum += 2*cgl2* (*X_132) * (*X_121) * s.dl31->y(l);
    sum += pow( cgl2, 2 )*( (*Xprime_022 ) * (* Xprime_022) * s.dl22->y(l)+ (*X_242) * (*X_220)*s.dl40->y(l));
    retval +=  SUM_STEP_SIZE * sum * ( 2*l+1 )/( l*( l+1 ) )*mCl.es[ n ]->fastY( l );
    //retval +=  SUM_STEP_SIZE * sum * ( 2*l+1 )/( ctnorm( l )*norml2( l ))*mCl.es[ n ]->fastY( l );
    ++X_022; ++X_220; ++X_132; ++X_121; ++X_242; ++Xprime_022;
  }
  return retval / ( 4 * M_PI );
}

double AllSkyLensing::xiLensed_minus( int n, double beta, double cgl2, const DlSplines& s, const Ximn& x )
{
  const double * X_022 = &(x.X_022[0]);
  const double * X_121 = &(x.X_121[0]);
  const double * X_132 = &(x.X_132[0]);
  const double * Xprime_022 = &(x.Xprime_022[0]);
  const double * X_220 = &(x.X_220[0]);
  const double * X_242 = &(x.X_242[0]);

  const int lmax = ( int )mCl.ts[ n ]->stop();
  double retval = 0;
  for ( int l = ( int )mCl.ts[ n ]->start(); l <= lmax; l += SUM_STEP_SIZE )
  {
    double sum = 0;
    sum = ( pow( *X_022, 2 )-1. )*s.dl2m2->y(l);
    sum += cgl2*( pow( *X_121, 2 ) * s.dl1m1->y(l) + pow( *X_132, 2 ) * s.dl3m3->y(l) );
    sum += 0.5 * pow( cgl2, 2 )*( 2 * pow( *Xprime_022, 2 ) * s.dl2m2->y(l)
                   + pow( *X_220, 2 ) * s.dl00->y(l)+ pow( *X_242, 2 ) * s.dl4m4->y(l) );
    retval +=  SUM_STEP_SIZE * sum * ( 2*l+1 )/( l*( l+1 ) )*mCl.es[ n ]->fastY( l );
    //retval +=  SUM_STEP_SIZE * sum * ( 2*l+1 )/(ctnorm( l )*norml2( l ) )*mCl.es[ n ]->fastY( l );
    ++X_022; ++X_121; ++X_132; ++Xprime_022; ++X_220; ++X_242;
  }
  return retval / ( 4 * M_PI );
}

double AllSkyLensing::xiLensed_X( int n, double beta, double cgl2, const DlSplines& s, const Ximn& x )
{
  const Spline * const dlm24 = s.dl4m2;

  const double * X_000 = &(x.X_000[0]);
  const double * X_022 = &(x.X_022[0]);
  const double * X_121 = &(x.X_121[0]);
  const double * X_132 = &(x.X_132[0]);
  const double * X_220 = &(x.X_220[0]);
  const double * X_242 = &(x.X_242[0]);
  const double * Xprime_022 = &(x.Xprime_022[0]);
  const double * Xprime_000 = &(x.Xprime_000[0]);

  const int lmax = (int)mCl.ts[ n ]->stop();
  double retval = 0;
  for ( int l = (int)mCl.ts[ n ]->start(); l <= lmax; l += SUM_STEP_SIZE )
  {
    double sum = 0;
    sum = ( ( (*X_022) * (*X_000) ) - 1. )*s.dl20->y(l);
    sum += cgl2* 2. *  (*Xprime_000) / sqrt( ( double )l*(( double ) l+1 ) ) *( (*X_121) * s.dl11->y(l)  // not -X_121
           + (*X_132)*s.dl3m1->y(l));
    sum += 0.5*pow( cgl2, 2 )*( ( 2* (*Xprime_022) * (*Xprime_000) + pow( *X_220, 2 ) )*s.dl20->y(l) + (*X_220)* (*X_242)*dlm24->y(l) );
    sum *= SUM_STEP_SIZE;
    sum *= (2.*l+1.)/ ( l*( l+1 ) ); //( sqrt(ctnorm( l ) )*norml2( l ) );
    retval += mCl.cs[ n ]->fastY(l)*sum;
    ++X_000; ++X_022; ++X_121; ++X_132; ++X_220; ++X_242;
    ++Xprime_022; ++Xprime_000;
  }
  return retval / ( 4 * M_PI );
}

// The class-static members of DlSplines:
typedef pair<int, int> poi;
map< poi, vector<double> >  AllSkyLensing::DlSplines::DlSplines::mFactor1;
map< poi, vector<double> >  AllSkyLensing::DlSplines::DlSplines::mFactor2;
map< poi, vector<double> >  AllSkyLensing::DlSplines::DlSplines::mFactor3;

AllSkyLensing::DlSplines::DlSplines( int lmax, double theta, Anchor *a )
{
    dl00 = dl( 0, 0, theta, lmax, a );
    dl11 = dl( 1, 1, theta, lmax, a );
    dl1m1 = dl( 1, -1, theta, lmax, a );
    dl20 = dl( 2, 0, theta, lmax, a );
    dl2m2 = dl( 2, -2, theta, lmax, a );
    dl22 = dl( 2, 2, theta, lmax, a );
    dl31 = dl( 3, 1, theta, lmax, a );
    dl3m1 = dl( 3, -1, theta, lmax, a );
    dl3m3 = dl( 3, -3, theta, lmax, a );
    dl40 = dl( 4, 0, theta, lmax, a );
    dl4m2 = dl( 4, -2, theta, lmax, a );
    dl4m4 = dl( 4, -4, theta, lmax, a );
    dl02 = dl20;
    dlm24 = dl4m2;
}


void AllSkyLensing::DlSplines::addToVector( voip& v, int i, int j ) { v.push_back( poi( i, j ) );}

void AllSkyLensing::DlSplines::initializeFactorTables(unsigned int lmax)
{
  vector<poi> v;
  v.reserve( 12 );
  addToVector( v, 0, 0 ); addToVector( v, 2, 0 );
  addToVector( v, 4, 0 ); addToVector( v, 1, 1 );
  addToVector( v, 1, -1 ); addToVector( v, 2, 2 );
  addToVector( v, 2, -2 ); addToVector( v, 3, -1 );
  addToVector( v, 3, 1 ); addToVector( v, 3, -3 );
  addToVector( v, 4, -2 ); addToVector( v, 4, -4 );

  vector<double> l2;        // l^2
  vector<double> l2lm1;     // l*(2*l+1)
  vector<double> twolm1olm1;  // (2*l-1)/(l-1)
  l2.reserve( lmax );
  l2lm1.reserve( lmax );
  twolm1olm1.reserve( lmax );
  for ( double l = 0; l <= lmax; ++l )
  {
    l2.push_back( l*l );
    const double twolm1 = 2*l-1;
    l2lm1.push_back( l * twolm1 );
    twolm1olm1.push_back( (l>1)? twolm1/(l-1) : 0 );
  }


  vector<poi>::const_iterator it = v.begin();
  vector<poi>::const_iterator endit = v.end();
  for ( ; it != endit; ++it )
  {
    vector<double>& currentfac1 = mFactor1[*it];
    vector<double>& currentfac2 = mFactor2[*it];
    vector<double>& currentfac3 = mFactor3[*it];

    currentfac1.reserve( lmax );
    currentfac2.reserve( lmax );
    currentfac3.reserve( lmax );

    const double m = ( *it ).first; const double n = ( *it ).second;
    const double mn = m*n;
    const double mm = m*m;
    const double nn = n*n;
    for ( double l = 0; l <= m; ++l )
    {
      currentfac1.push_back( 0 );
      currentfac2.push_back( 0 );
      currentfac3.push_back( 0 );
    }
    double lastsqrtfactor = 0;
    for ( double l = m+1; l <= lmax; ++l )
    {
      const double sqrtfactor = sqrt( l2[ (int) l ]-mm )* sqrt( l2[ (int) l ]-nn );
      currentfac1.push_back( l2lm1[ (int) l ]/ sqrtfactor );
      currentfac2.push_back( mn * twolm1olm1[ (int) l ]/ sqrtfactor );
      currentfac3.push_back( l * lastsqrtfactor /( ( l - 1. ) * sqrtfactor ) );
      lastsqrtfactor = sqrtfactor;
    }
  }
}

// Blanco, Florez and Bermejo - Evaluation of the Rotation Matrices in the Basis of Real Spherical Harmonics
// eg: http://www1.elsevier.com/homepage/saa/eccc3/paper48/eccc3.html
const Spline * AllSkyLensing::DlSplines::dl( int m, int n, double alpha, int lmax, Anchor* a )
{
  if ( lmax < 1000 ) throw Bad_Error( "lmax too low in AllSkyLensing" );
  Spline * s = new Spline( lmax+1, "dlSpline", a );

  const double cosa = cos( alpha );
  double l = 0;
  for ( ; l < m; ++l )
    s->setY( ( int )l, 0 );

  l = m;
  s->setY(( int ) l,sqrt( factorial[ 2*( unsigned int )l ].value/( factorial[  ( unsigned int )l+n ].value*factorial[ ( unsigned int )l-n ].value ) )*pow( cos( alpha/2. ), l+n )*pow( -sin( alpha/2. ), l-n ) );

  l=m+1;
  s->setY(( int ) l, ( l*cosa-n )*sqrt( factorial[  2*( unsigned int )l-1 ].value/( factorial[  ( unsigned int )l+n ].value*factorial[ ( unsigned int )l-n ].value ) )*pow( cos( alpha/2. ), l-1+n )*pow( -sin( alpha/2. ), l-1-n ) );
  ++l;

  const pair<int, int> mnpair( m, n );
  double * fac1it = &mFactor1[ mnpair ][ 0 ];
  double * fac2it = &mFactor2[ mnpair ][ 0 ];
  double * fac3it = &mFactor3[ mnpair ][ 0 ];

  int il = ( int )l;

  for ( int i = 0; i < il; ++i )
  {
    ++fac1it; ++fac2it; ++fac3it;
  }

  for ( ; il <= lmax; ++il )
  {
    double value = ( *fac1it * cosa - *fac2it )*s->y( il-1 );
    value -= *fac3it * s->y( il-2 );
    s->setY( il, value );
    ++fac1it; ++fac2it; ++fac3it;
  }
  return s;
}

AllSkyLensing::Ximn::Ximn( const double theta, const double sigma2, const double c_gl2,
                           const int lstart, const int lend  )
{
  const int maxSize = (int) (lend / SUM_STEP_SIZE) + 1;
  X_000.reserve( maxSize );
  X_220.reserve( maxSize );
  X_022.reserve( maxSize );
  X_121.reserve( maxSize );
  X_132.reserve( maxSize );
  X_242.reserve( maxSize );
  Xprime_000.reserve( maxSize );
  Xprime_022.reserve( maxSize );
  for ( double l = lstart; l <= lend; l += SUM_STEP_SIZE )
  {
    X_000.push_back( exp( -l*( l+1. )*sigma2*0.25 ) );
    X_220.push_back( 0.25*sqrt(( l+2. )*( l-1. )*l*( l+1 ) )* exp( -( l*( l+1. )-2. )*sigma2*0.25 ) );
    X_022.push_back( exp( -0.25*( l*( l+1. )-4. )*sigma2 ) );
    X_121.push_back( -0.5*sqrt( ( l+2. )*( l-1. ) )*exp( -0.25*( l*( l+1 )-8./3. )*sigma2 ) );
    X_132.push_back( -0.5*sqrt( ( l+3. )*( l-2. ) )*exp( -0.25*( l*( l+1 )-20./3. )*sigma2 ) );
    X_242.push_back( .25*sqrt( ( l+4. )*( l+3. )*( l-2 )*( l-3 ) )*exp( -0.25*( l*( l+1 )-10. )*sigma2 ) );
    Xprime_000.push_back( -0.25*l*(l+1.)*X_000.back() );
    Xprime_022.push_back( -(l*( l+1.)-4.)*0.25*X_022.back() );
  }
}

AllSkyLensing::LensingDifference::LensingDifference( const int size )
{
  ts = new Spline( size, "tsDifference" );
  cs = new Spline( size, "csDifference" );
  e = new Spline( size, "esDifference" );
  b = new Spline( size, "bsDifference" );

  for ( int i = ts->first(); i <= ts->last(); ++i )
  {
    ts->setY( i, 0 );
    cs->setY( i, 0 );
    e->setY( i, 0 );
    b->setY( i, 0 );
  }
}

AllSkyLensing::LensingDifference::~LensingDifference()
{
  delete ts; delete cs;
  delete e; delete b;
}
