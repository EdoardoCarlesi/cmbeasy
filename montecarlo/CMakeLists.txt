# FIND_PACKAGE(MPImod QUIET)

FIND_PACKAGE(MPI)

IF(NOT MPI_FOUND)
MESSAGE(STATUS "No MPI C++ compiler found, the monte-carlo drivers in the cmbeasy/montecarlo subdir will not be built")
ENDIF(NOT MPI_FOUND)

IF(MPI_FOUND)
MESSAGE(STATUS "MPI C++ compiler found")
#SET(MPI_INCLUDE_PATH /home/carlesi/Programs/openmpi/include/)
#SET(MPI_LIBRARIES /home/carlesi/Programs/openmpi/lib/libmpi_cxx.so)
#SET(MPI_LIBRARIES)
#SET(MPI_INCLUDE_PATH /usr/include/mpi/)
#SET(MPI_LIBRARIES /usr/lib/libmpi_cxx.so)
SET(MPI_COMPILE_FLAGS)
SET(MPI_LINK_FLAGS -Wl --export-dynamic)

INCLUDE_DIRECTORIES(${MPI_INCLUDE_PATH} ${PROJECT_SOURCE_DIR}/montecarlo)
#INCLUDE_DIRECTORIES(/usr/local/include/ ${PROJECT_SOURCE_DIR}/montecarlo)

  #  SET(KMEANSDIR kmlocal-1.7.1/src)
  #SET(KMEANSFILES ${KMEANSDIR}/KM_ANN.cpp ${KMEANSDIR}/KMeans.cpp
  #                ${KMEANSDIR}/KMterm.cpp ${KMEANSDIR}/KMrand.cpp
  #                ${KMEANSDIR}/KCutil.cpp ${KMEANSDIR}/KCtree.cpp
  #                ${KMEANSDIR}/KMdata.cpp ${KMEANSDIR}/KMcenters.cpp
  #                ${KMEANSDIR}/KMfilterCenters.cpp ${KMEANSDIR}/KMlocal.cpp
  #   )

  # new style montecarlo driver(s)
  SET(MCSOURCES mcmain.cc mcmodel.cc mcrunner.cc mclikelihoodcalculator.cc mcsettings.cc
                mctaskinfo.cc mcutils.cc
                ${MATHOBJECT}/multigaussian.cc
                fastslowstepper.cc nestedsampler.cc
                ${KMEANSFILES}
     )

  MACRO(MCDRIVER source)
    PARSE_ARGUMENTS(THIS_MCDRIVER "EXTRA_SOURCES" "" ${ARGN})
    GET_FILENAME_COMPONENT(name ${source} NAME_WE)
    ADD_EXECUTABLE(${name} ${source} ${THIS_MCDRIVER_EXTRA_SOURCES} ${MCSOURCES} ${BASIC_WMAPSOURCES} ${VDEFILES})
    TARGET_LINK_LIBRARIES(${name} ${MPI_LIBRARIES} ${GSL_LIBRARIES} ${WMAPLIB})
    SET_TARGET_PROPERTIES(${name} PROPERTIES
      COMPILE_FLAGS "${MPI_COMPILE_FLAGS} -DMONTECARLO"
      LINK_FLAGS    "${MPI_LINK_FLAGS}")
  ENDMACRO(MCDRIVER)


  ################################################################
  # developement stuff, not needed for release
  SET(CMBEASY_MONTECARLO_MAINTAINER_MODE TRUE)
  INCLUDE(${CMBEASY_MAINTAINER_DIR}/CMakeLists.txt OPTIONAL)
  SET(CMBEASY_MONTECARLO_MAINTAINER_MODE)
  ################################################################

  # old style monolithic mc driver
  ADD_EXECUTABLE(mc_general mc_general.cc ${MATHOBJECT}/multigaussian.cc ${BASIC_WMAPSOURCES})
  INSTALL(TARGETS mc_general
          RUNTIME DESTINATION bin
         )
  TARGET_LINK_LIBRARIES(mc_general ${GSL_LIBRARIES} ${WMAPLIB}
                                   ${MPI_LIBRARIES}
                       )
  SET_TARGET_PROPERTIES(mc_general PROPERTIES
                                   COMPILE_FLAGS "${MPI_COMPILE_FLAGS}"
                                   LINK_FLAGS    "${MPI_LINK_FLAGS}"
                       )

  MCDRIVER(mcsimpleexample.cc)
  MCDRIVER(fastslowexample.cc)
  MCDRIVER(dummymodel.cc)
  MCDRIVER(mc_vde.cc)
  MCDRIVER(mc_decaying_dedm.cc)
  MCDRIVER(mc_uncoupled_quintessence.cc)
  MCDRIVER(mc_dedm_corasaniti.cc)
  MCDRIVER(mc_cn_cmb.cc)

ENDIF(MPI_FOUND)

