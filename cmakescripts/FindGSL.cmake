##
## Try to find gnu scientific library GSL
## (see http://www.gnu.org/software/gsl/)
## Once run this will define:
##
## GSL_FOUND       = system has GSL lib
##
## GSL_LIBRARIES   = full path to the libraries
##    on Unix/Linux with additional linker flags from "gsl-config --libs"
##
## CMAKE_GSL_CXX_FLAGS  = Unix compiler flags for GSL, essentially
## "`gsl-config --cxxflags`"
##
## GSL_INCLUDE_DIR      = where to find headers
##
## GSL_LINK_DIRECTORIES = link directories, useful for rpath on Unix
## GSL_EXE_LINKER_FLAGS = rpath on Unix
##
## Felix Woelk 07/2004
## minor corrections Jan Woetzel
##
## modified by Georg Robbers to prefer the version in $GSLCMBDIR 08/2006
##
## www.mip.informatik.uni-kiel.de
## --------------------------------

IF(WIN32)
  MESSAGE(SEND_ERROR "FindGSL.cmake: gnu scientific library GSL not (yet)
supported on WIN32")

ELSE(WIN32)
  IF(UNIX)
    SET(GSL_CONFIG_PREFER_PATH $ENV{GSLCMBDIR}/bin)
    FIND_PROGRAM(GSL_CONFIG  gsl-config
                 ${GSL_CONFIG_PREFER_PATH}
                 DOC "Path to the gsl-config program"
                 NO_DEFAULT_PATH
                )
    IF(NOT GSL_CONFIG)
      FIND_PROGRAM(GSL_CONFIG  gsl-config)
    ENDIF(NOT GSL_CONFIG)
    # MESSAGE("DBG GSL_CONFIG ${GSL_CONFIG}")

    IF (GSL_CONFIG)
      # set CXXFLAGS to be fed into CXX_FLAGS by the user:
      # GR: changed from SET to use OUTPUT_VARIABLE and EXECUTE_PROCESS
      #     so that the flags get set at configure time, not build time
      #     (which doesn't work with Xcode)
      # SET(GSL_CXX_FLAGS "`${GSL_CONFIG} --cflags`")
      EXECUTE_PROCESS(COMMAND ${GSL_CONFIG} --cflags OUTPUT_VARIABLE GSL_CXX_FLAGS)

      # set INCLUDE_DIRS to prefix+include
      EXEC_PROGRAM(${GSL_CONFIG}
        ARGS --prefix
        OUTPUT_VARIABLE GSL_PREFIX)
      SET(GSL_INCLUDE_DIR ${GSL_PREFIX}/include CACHE STRING INTERNAL)

      # set link libraries and link flags
      # GR: see above
      # SET(GSL_LIBRARIES "`${GSL_CONFIG} --libs`")
      EXECUTE_PROCESS(COMMAND ${GSL_CONFIG} --libs OUTPUT_VARIABLE GSL_LIBRARIES
                                                   OUTPUT_STRIP_TRAILING_WHITESPACE)
      #      MESSAGE("DBG  GSL_LINK_DIRECTORIES=${GSL_LIBRARIES}")

      ## extract link dirs for rpath
      EXEC_PROGRAM(${GSL_CONFIG}
        ARGS --libs
        OUTPUT_VARIABLE GSL_CONFIG_LIBS )

      ## split off the link dirs (for rpath)
      ## use regular expression to match wildcard equivalent
      ## "-L*<endchar>"
      ## with <endchar> is a space or a semicolon
      STRING(REGEX MATCHALL "[-][L]([^ ;])+"
        GSL_LINK_DIRECTORIES_WITH_PREFIX
        "${GSL_CONFIG_LIBS}" )

        #      MESSAGE("DBG GSL_LINK_DIRECTORIES_WITH_PREFIX=${GSL_LINK_DIRECTORIES_WITH_PREFIX}")

      ## remove prefix -L because we need the pure directory for LINK_DIRECTORIES

      IF (GSL_LINK_DIRECTORIES_WITH_PREFIX)
        STRING(REGEX REPLACE "[-][L]" "" GSL_LINK_DIRECTORIES ${GSL_LINK_DIRECTORIES_WITH_PREFIX} )
      ENDIF (GSL_LINK_DIRECTORIES_WITH_PREFIX)
      SET(GSL_EXE_LINKER_FLAGS "-Wl,-rpath,${GSL_LINK_DIRECTORIES}" CACHE
STRING INTERNAL)
      #      MESSAGE("DBG  GSL_LINK_DIRECTORIES=${GSL_LINK_DIRECTORIES}")
      #      MESSAGE("DBG  GSL_EXE_LINKER_FLAGS=${GSL_EXE_LINKER_FLAGS}")

      #      ADD_DEFINITIONS("-DHAVE_GSL")
      #      SET(GSL_DEFINITIONS "-DHAVE_GSL")
      MARK_AS_ADVANCED(
        GSL_CXX_FLAGS
        GSL_INCLUDE_DIR
        GSL_LIBRARIES
        GSL_LINK_DIRECTORIES
        GSL_DEFINITIONS
        GSL_EXE_LINKER_FLAGS
      )
      MESSAGE(STATUS "Using GSL from ${GSL_PREFIX}")

    ELSE(GSL_CONFIG)
      MESSAGE("FindGSL.cmake: gsl-config not found. Please set it manually. GSL_CONFIG=${GSL_CONFIG}")
    ENDIF(GSL_CONFIG)

  ENDIF(UNIX)
ENDIF(WIN32)


IF(GSL_LIBRARIES)
  IF(GSL_INCLUDE_DIR OR GSL_CXX_FLAGS)

    SET(GSL_FOUND 1)

  ENDIF(GSL_INCLUDE_DIR OR GSL_CXX_FLAGS)
ENDIF(GSL_LIBRARIES)




