#ifndef READ_ARCHIVE_MAP_H
#define READ_ARCHIVE_MAP_H
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
// =============================================================================
// Interface block for the routine(s) in 'read_archive_map.f90'
// =============================================================================
#include "port_helper.h"

void  Read_Archive_Map( std::string FileName, real4_1d& Temp, real4_1d& NObs, long& NPix, int& Status);

#endif // READ_ARCHIVE_MAP_H
