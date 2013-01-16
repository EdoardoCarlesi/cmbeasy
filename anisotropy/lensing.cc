#include "lensing.h"

#include "controlpanel.h"
#include "allskylensing.h"

Lensing::Lensing( Cosmos& cosmos, const ControlPanel& cp, CL& cl, const CmbCalc& cmbCalc )
                : mCosmos( cosmos ), mControlPanel( cp ), mCl( cl ), mCmbCalc( cmbCalc ),
                  mBesselFileName( string() )
{
}

void Lensing::doLensing( const LensingMethod method )
{
  CL* lensed = lensedCls( method );
  mCl = *lensed;
}

CL* Lensing::lensedCls( const LensingMethod method )
{
  if ( methodUsed( method ) == AllSky )
  {
    AllSkyLensing allSkyLensing( mCosmos, mControlPanel, mCl, mCmbCalc );
    return allSkyLensing.lensedCls();
  }
  else if ( methodUsed( method ) == FlatSky )
  {
    throw Bad_Error("FlatSky lensing is no longer supported, use AllSky instead");
    /*
    CmbCalc& cmbc = const_cast<CmbCalc&>( mCmbCalc );
    ControlPanel& cp = const_cast<ControlPanel&>( mControlPanel );
    CL* cl = new CL( mCl );
    // only lense cl, don't write anything to disc
    FlatSkyLensing( mCosmos, cp, *cl, cmbc.lval(), string(), mBesselFileName );
    return cl;
    */
  }
  return 0; //never happens
}

void Lensing::dumpLensedCls( const string& fileName,  const LensingMethod method )
{
  if ( methodUsed( method ) == AllSky )
  {
    AllSkyLensing allSkyLensing( mCosmos, mControlPanel, mCl, mCmbCalc );
    allSkyLensing.dumpLensedCls( fileName );
  }
  else if ( methodUsed( method ) == FlatSky )
  {
    throw Bad_Error("FlatSky lensing is no longer supported, use AllSky instead");
    /*
    CmbCalc& cmbc = const_cast<CmbCalc&>( mCmbCalc );
    ControlPanel& cp = const_cast<ControlPanel&>( mControlPanel );
    CL* cl = new CL( mCl );
    FlatSkyLensing( mCosmos, cp, *cl, cmbc.lval(), fileName, mBesselFileName );
    */
  }
}

Lensing::LensingMethod Lensing::methodUsed( LensingMethod m )
{
  if ( m != Automatic )
    return m;

  if ( const_cast<ControlPanel&>( mControlPanel ).isAllSkyLensing() )
  {
    return AllSky;
  }
  else
  {
    if ( mBesselFileName.empty() )
       throw Bad_Error( "Using FlatSkyLensing, but mBesselFileName is empty. Aborting. " );
    return FlatSky;
  }

  throw( Bad_Error( "Unable to determine lensing method" ) ); return Automatic; // never happens
}
