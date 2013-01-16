/****************************************************************
 *   This file is part of the C++ port of the Fortran likelihood
 *   code for the WMAP 7yr release provided by the WMAP team
 *   at http://lambda.gsfc.nasa.gov/ .
 *
 *   The code was ported by Georg Robbers for easier
 *   interfacing with cmbeasy (http://www.cmbeasy.org).
 *   Bugs in this port should be reported to the
 *   cmbeasy authors (bugs@cmbeasy.org).
 ****************************************************************/

#include "WMAP_7yr_util.h"

#include "port_helper.h"

#include <string>
#include <iostream>
#include <iomanip>
#include <ctime>

namespace  wmap_util
{

  using namespace std;

  bool  bwmap_likelihood_ok, bwmap_likelihood_good;
  const bool  bwmap_warnings_are_errors = false;

  namespace   // unnamed, simulates "private"
  {
     const int  n_errs_max = 10;
     string  err_msg[n_errs_max];
     int  n_errs, err_code[n_errs_max];

     const int  n_warn_max = 10;
     string  warn_msg[n_warn_max];
     int  n_warn, warn_code[n_warn_max];

     struct Timer
     {
         string desc;
         clock_t  start_count, end_count;
         int  checkpoint_count;
         bool  stub;
     };

     int  tlun;
     int  n_timers = 0;
     const int max_n_timers = 6;
     Timer timers[max_n_timers];
  }

//X  don't need this in c++:
//X   void  get_free_lun( lun )
//X   {
//X 	int, intent(out)  lun;
//X
//X 	int, save  last_lun = 19;
//X 	bool  usingd;
//X
//X 	lun = last_lun;
//X 	while ( true ) {
//X 		inquire( unit=lun, opened=usingd );
//X 		if (  ! usingd ) break;
//X 		lun = lun + 1;
//X 	}
//X
//X 	last_lun = lun;
//X   }


bool  wmap_likelihood_ok() { return bwmap_likelihood_ok; }
bool  wmap_likelihood_good() { return bwmap_likelihood_good; }



void  wmap_likelihood_error_init( )
{
  bwmap_likelihood_ok = true;
  bwmap_likelihood_good = true;
  n_errs = 0;
  n_warn = 0;
}


void  wmap_likelihood_error( string msg, int code )
{
  bwmap_likelihood_ok =  false ;
  bwmap_likelihood_good =  false ;

  if ( n_errs < n_errs_max )
  {
    n_errs = n_errs + 1;
    err_msg[n_errs] = msg;
    err_code[n_errs] = code;
  }
  else
  {
    cout <<  "*** error log full";
    exit( -1 );
  }
}

void  wmap_likelihood_warning( string msg, int code )
{
  if ( bwmap_warnings_are_errors )
    bwmap_likelihood_ok =  false;

  bwmap_likelihood_good =  false;

  if ( n_warn < n_warn_max )
  {
    n_warn = n_warn + 1;
    warn_msg[n_warn] = msg;
    warn_code[n_warn] = code;
  }
  else
  {
    cout <<  "*** warning log full";
    exit( -1 );
  }
}

void  wmap_likelihood_error_report( )
{
  int  i;

  cout <<  "------------------------------------------------------------" << endl;
  cout <<  "WMAP likelihood evaluation report:" << endl;

  if ( bwmap_likelihood_ok )
  {
    cout <<  "no errors" << endl;
  }
  else
  {
    cout <<  "number of errors = " << n_errs << endl;
    for (i = 1; i <= n_errs; ++i )
    {
      cout <<  "" << endl;
      cout <<  "error #" << i << endl;
      cout <<  err_msg[i] << endl;
      cout <<  "error code = " << err_code[i] << endl;
    }
  }

  if ( n_warn > 0 )
  {
    cout <<  "number of warnings = " << n_warn << endl ;
    for (i = 1; i <= n_warn; ++i )
    {
      cout <<  "" << endl ;
      cout <<  "warning #" << i << "::" << endl ;
      cout <<  warn_msg[i] << endl ;
      cout <<  "warning code = " << warn_code[i] << endl ;
    }
  }
  cout <<  "------------------------------------------------------------" << endl;

}

void wmap_timing_start(std::string desc)
{
    cout << setw(8) << setprecision(4);
    int  k;
    REAL elapsed_sec;

    //	if ( n_timers == 0 ) then
    //		call get_free_lun( tlun )
    //		open(tlun,file=ofn_timing,action='write',status='unknown')
    //	end if

    if ( n_timers >= max_n_timers )  {
        cout <<  "*** too many timing levels" << endl;
        exit(-1);
    }

    if ( n_timers > 0  &&  timers[n_timers].stub )  {
        timers[n_timers].stub = false;
        cout << "        .......   ";
        cout << timers[n_timers].desc << endl;
    }

    n_timers = n_timers + 1;
    k = n_timers;

    timers[k].desc = desc;
    timers[k].start_count = clock();
    timers[k].checkpoint_count = -1;
    timers[k].stub = true;
}


void  wmap_timing_checkpoint(std::string desc)
{
    clock_t start_count, end_count; 

    if ( timers[n_timers].checkpoint_count == -1 )  {
        start_count = timers[n_timers].start_count;
        if ( timers[n_timers].stub ) {
            timers[n_timers].stub = false;
            cout << "        ....... ";
            cout << timers[n_timers].desc << endl;
        }
    } else {
        start_count = timers[n_timers].checkpoint_count;
    }
    end_count = clock();

    REAL elapsed_sec =(float)(end_count-start_count)/CLOCKS_PER_SEC;
    cout << "        " << n_timers << "  " << " - ";
    cout << elapsed_sec << " " << desc << endl;

    timers[n_timers].checkpoint_count = end_count;
} // 


void wmap_timing_end( )
{
    clock_t /*start_count,*/ end_count;
    REAL elapsed_sec;

    if ( n_timers == 0 )  {
        cout <<  "*** n_timers == 0 in wmap_timing_end" << endl ;
        exit(-1);
    }

    end_count = clock();
    elapsed_sec = (float)(end_count-timers[n_timers].start_count)/CLOCKS_PER_SEC;

    cout << "        " << n_timers << "  ";
    cout << elapsed_sec << " - " << timers[n_timers].desc << endl;
    n_timers = n_timers - 1;

    //	if ( n_timers == 0 ) then
    //		close(tlun)
    //	end if
} // 


} // end namespace

