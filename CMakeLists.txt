PROJECT(cmbeasy)

cmake_minimum_required(VERSION 2.6.4)

SUBDIRS(montecarlo)


SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin
           CACHE PATH "Single output directory for building all executables.")

SET(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmakescripts ${CMAKE_MODULE_PATH})
INCLUDE(${PROJECT_SOURCE_DIR}/cmakescripts/extracommands.inc)

# make sure these refer to absolute paths
SET(MATHOBJECT ${PROJECT_SOURCE_DIR}/mathobject)
SET(ANISOTROPY ${PROJECT_SOURCE_DIR}/anisotropy)
SET(MASSFUNCTION ${PROJECT_SOURCE_DIR}/massfunction)
SET(COSMOS ${PROJECT_SOURCE_DIR}/cosmos)
SET(VDE ${PROJECT_SOURCE_DIR}/vde)
SET(CHAIN ${PROJECT_SOURCE_DIR}/chain)
SET(ANALYZE ${PROJECT_SOURCE_DIR}/analyze)
SET(WMAP ${ANALYZE}/wmap7)
SET(GUI ${PROJECT_SOURCE_DIR}/gui)
SET(COUPLED ${PROJECT_SOURCE_DIR}/coupled)

FIND_PACKAGE(Doxygen)
FIND_PACKAGE(GSL)

#find stuff even on 64-bit SUSE systems
SET(CMAKE_LIBRARY_PATH ${CMAKE_LIBRARY_PATH} /usr/lib64/qt4/lib64)
SET(CMAKE_INCLUDE_PATH ${CMAKE_INCLUDE_PATH} /usr/lib64/qt4/include)

SET(QT_MIN_VERSION 4.2.0)
FIND_PACKAGE(Qt4)
IF(QT_FOUND)
  MESSAGE(STATUS "Found Qt4 headers: ${QT_INCLUDE_DIR}")
  MESSAGE(STATUS "Found Qt4 libraries: ${QT_LIBRARY_DIR}")
ELSE(QT_FOUND)
  MESSAGE(STATUS "Qt4.2.0 or higher not found.")
  SET(QT_INCLUDE_DIR "")
  SET(QT_DEFINITIONS)
  SET(QT_LIBRARIES)
ENDIF(QT_FOUND)

OPTION(BUILD_GUI "Build the graphical user interface (needs the Qt library)" ON)
IF(NOT QT_FOUND AND BUILD_GUI)
  MESSAGE("Warning: The Qt library could not be found.
       Without Qt, the graphical user interface will not be built.")
  SET(BUILD_GUI 0)
ENDIF(NOT QT_FOUND AND BUILD_GUI)

OPTION(WMAP7 "Build the c++ WMAP7yr likelihood code" ON)

FIND_PACKAGE(Cfitsio)
IF(NOT CFITSIO_FOUND AND WMAP7)
  MESSAGE("Warning: The cfitsio library could not be found.
           Without cfitsio, the support for WMAP 7-year likelihoods will not be built.")
  SET(WMAP7 OFF)
ENDIF(NOT CFITSIO_FOUND AND WMAP7)

IF(WMAP7)
  SET(WMAPINC ${WMAP} ${CFITSIO_INCLUDE_DIR})
  SET(WMAPLIB ${GSL_LIBRARIES} ${CFITSIO_LIBRARY})
  SET(WMAPFILES ${WMAP}/read_archive_map.cc ${WMAP}/read_fits.cc ${WMAP}/WMAP_7yr_options.cc
                ${WMAP}/WMAP_7yr_util.cc ${WMAP}/WMAP_7yr_tt_pixlike.cc
                ${WMAP}/WMAP_7yr_tt_beam_and_ptsrc_chisq.cc
                ${WMAP}/WMAP_7yr_teeebb_pixlike.cc ${WMAP}/WMAP_7yr_tetbeebbeb_pixlike.cc
                ${WMAP}/WMAP_7yr_gibbs.cc ${WMAP} ${WMAP}/br_mod_dist.cc ${WMAP}/WMAP_7yr_likelihood.cc)
  SET(WMAPDEFINE -DOPTIMIZE -DFASTERTT)
  #SET(WMAPDEFINE -DOPTIMIZE -DFASTERTT -DUSE_LOWELL_TBEB -DUSE_HIGHELL_TB -DTIMING)
ELSE (WMAP7)
  SET(WMAPINC)
  SET(WMAPLIB)
  SET(WMAPFILES)
  SET(WMAPDEFINE -DNOWMAP7)
ENDIF(WMAP7)

SET(VDEFILES ${VDE}/vdecosmos.cc ${VDE}/vectorde.cc ${VDE}/grid.cc)

SET(COUPLEDFILES ${COUPLED}/massiveneutrinosnew.cc 
		 ${COUPLED}/coupledexp.cc ${COUPLED}/coupledquintcosmos.cc
#		 ${COUPLED}/coupledleaping.cc 
		${COUPLED}/perturbationtracker.cc
		 ${COUPLED}/exponentialcoupling.cc
		 ${COUPLED}/decayingdedmcosmos.cc
		 ${COUPLED}/cninvariant.cc 
		 ${COUPLED}/invariant.cc
		 ${COUPLED}/cncosmos.cc
		 ${COUPLED}/speedycoupledinvariant.cc
		 )

SET(MARGINALIZEFILES ${MARGINALIZE}/mychainshop.cc)

SET(COSMOSFILES ${COSMOS}/cosmos.cc ${COSMOS}/recfast.cc ${COSMOS}/quintessence.cc
                ${COSMOS}/exponential.cc ${COSMOS}/ratra.cc  ${COSMOS}/arthur.cc ${COSMOS}/quintcosmos.cc
                ${COSMOS}/perturbation.cc ${COSMOS}/synchronous.cc ${COSMOS}/perturbationfactory.cc
                ${COSMOS}/quintsynchronous.cc  ${COSMOS}/arbitrary.cc ${COSMOS}/recfastalpha.cc
                ${COSMOS}/crossoverfield.cc ${COSMOS}/crossover.cc ${COSMOS}/speedyinvariant.cc
                ${COSMOS}/massiveneutrinos.cc ${COSMOS}/speedydeinvariant.cc ${COSMOS}/distcosmos.cc
                ${COSMOS}/distcelestine.cc ${COSMOS}/arbitraryomega.cc ${COSMOS}/omegastep.cc
   )


SET(MATHOBJECTFILES ${MATHOBJECT}/miscmath.cc ${MATHOBJECT}/spline.cc
                    ${MATHOBJECT}/anchor.cc ${MATHOBJECT}/bessel.cc
                    ${MATHOBJECT}/data.cc  ${MATHOBJECT}/snedata.cc
                    ${MATHOBJECT}/splinetools.cc
                    ${MATHOBJECT}/rollingaverage.cc)

SET(ANALYZEFILES  ${ANALYZE}/LymanAlpha/PatMcDonald/mynrutils.cc
                  ${ANALYZE}/LymanAlpha/PatMcDonald/GlassInterpolator.cc
                  ${ANALYZE}/LymanAlpha/PatMcDonald/ObjectInSpace.cc
                  ${ANALYZE}/LymanAlpha/PatMcDonald/splint.cc
                  ${ANALYZE}/LymanAlpha/PatMcDonald/mysvdlsq.cc
                  ${ANALYZE}/LymanAlpha/PatMcDonald/pythag.cc
                  ${ANALYZE}/LymanAlpha/PatMcDonald/svbksb.cc
                  ${ANALYZE}/LymanAlpha/PatMcDonald/LyaFchi2Interpolator.cc
                  ${ANALYZE}/LymanAlpha/PatMcDonald/myutils.cc
                  ${ANALYZE}/LymanAlpha/PatMcDonald/spline.cc
                  ${ANALYZE}/LymanAlpha/PatMcDonald/svdcmp.cc
                  ${ANALYZE}/LymanAlpha/PatMcDonald/oldspline.cc
                  ${ANALYZE}/LymanAlpha/PatMcDonald/oldsplint.cc
                  ${ANALYZE}/sdss/sdss_lrg/sdsslrg.cc
                  ${ANALYZE}/sdss/sdss_lrg10/sdsslrg10.cc
                  ${ANALYZE}/cbi2/cbi2.cc
                  ${ANALYZE}/newdatcl/newdatclchi2.cc
                  ${ANALYZE}/sne/union/union.cc
                  ${ANALYZE}/sne/union2/union2.cc
                  ${ANALYZE}/sne/union2.1/union21.cc
                  ${ANALYZE}/sne/constitution/constitution.cc
                  ${ANALYZE}/sne/legacy/legacy.cc)

SET(ANISOTROPYFILES ${ANISOTROPY}/cmbcalc.cc ${ANISOTROPY}/scalarintegrator.cc
                    ${ANISOTROPY}/tensorintegrator.cc  ${ANISOTROPY}/integrator.cc
                    ${ANISOTROPY}/lensing.cc ${ANISOTROPY}/analyzethis.cc ${ANISOTROPY}/model.cc
                    ${ANISOTROPY}/DataManager.cc ${ANISOTROPY}/scalartrapezoidintegrator.cc
                    ${ANISOTROPY}/tensortrapezoidintegrator.cc ${ANISOTROPY}/allskylensing.cc
                    ${ANISOTROPY}/newdatreader.cc ${ANISOTROPY}/cl.cc)

SET(MASSFUNCTIONFILES ${MASSFUNCTION}/mass_function.cc)


SET(CHAINFILES ${CHAIN}/chainshop.cc ${CHAIN}/lowlevelplot.cc ${CHAIN}/postscriptplot.cc)

IF(BUILD_GUI)
  SET(GUIFILES ${GUI}/cmbeasywidget.cc ${GUI}/cmbmainwindow.cc ${GUI}/plotwidget.cc ${GUI}/tipdialog.cc
               ${GUI}/helpwindow.cc ${GUI}/configurationdialog.cc ${GUI}/q4buttongroup.cc
               ${GUI}/plotcontrol.cc
               ${GUI}/parameterdock.cc
               ${GUI}/modelinspectordock.cc
               ${GUI}/automaticrunner.cc
     )

  SET(GUIMOCHEADERS ${GUI}/cmbeasywidget.h ${GUI}/cmbmainwindow.h ${GUI}/plotwidget.h ${GUI}/tipdialog.h
               ${GUI}/helpwindow.h ${GUI}/configurationdialog.h ${GUI}/q4buttongroup.h
               ${GUI}/plotcontrol.h
               ${GUI}/parameterspinbox.h
               ${GUI}/parameterdock.h
               ${GUI}/modelinspectordock.h
               ${GUI}/automaticrunner.h
     )

  QT4_WRAP_CPP(GUIFILES ${GUIMOCHEADERS})

  SET(DESIGNERFILES ${GUI}/design_cmbeasywidget.ui ${GUI}/design_tipdialog.ui
                    ${GUI}/design_configurationdialog.ui
                    ${GUI}/design_plotcontrol.ui
                    ${GUI}/design_parameterdock.ui
                    ${GUI}/design_modelinspector.ui
     )

  QT4_WRAP_UI(UI_HEADER_FILES ${DESIGNERFILES})
  LIST(APPEND GUIFILES ${UI_HEADER_FILES})

ENDIF(BUILD_GUI)

########## Includes and Flags ######################################
ADD_DEFINITIONS(-DPRERELEASE ${WMAPDEFINE})

INCLUDE_DIRECTORIES(${MATHOBJECT} ${ANISOTROPY} ${COSMOS} ${GUI} ${GSL_INCLUDE_DIR} ${CHAIN}
                    ${WMAPINC} ${QT_INCLUDE_DIR} ${QT_QT_INCLUDE_DIR} ${QT_QTCORE_INCLUDE_DIR}
                    ${QT_QTGUI_INCLUDE_DIR}
                    ${ANALYZE}/LymanAlpha/PatMcDonald/
                    ${ANALYZE}/sdss/sdss_lrg 
		    ${ANALYZE}/sdss/sdss_lrg10
			${ANALYZE}/cbi2 ${ANALYZE}/newdatcl ${QT_INCLUDE_DIR}
                    ${CMAKE_BINARY_DIR}
		    ${COUPLED}
		    ${VDE} ${MASSFUNCTION}
	            ${MARGINALIZE}
		    )

IF(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE Release CACHE STRING
      "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
      FORCE)
ENDIF(NOT CMAKE_BUILD_TYPE)

SET( CMAKE_CXX_FLAGS_RELEASE "-O2 -DNDEBUG -fno-unit-at-a-time" CACHE STRING
    "Flags used by the C++ compiler during release builds."
   )

IF(CMAKE_COMPILER_IS_GNUCXX)
  SET( CMAKE_CXX_FLAGS_DEBUG "-O0 -g -Wall -ansi -pedantic -Wno-long-long" CACHE STRING
       "Flags used by the C++ compiler during debug builds."
     )
ENDIF(CMAKE_COMPILER_IS_GNUCXX)


IF(APPLE)
 MESSAGE("when building on mac os x panther, edit CMakeLists.txt and enable the MACOSX_PANTHER define")
 #ADD_DEFINITIONS(-DMACOSX_PANTHER)
ENDIF(APPLE)

OPTION(THREADS "Use multi-threading (needs OpenMP support from the compiler)" ON)
IF(THREADS)
  FIND_PACKAGE(OpenMP)
  IF(OPENMP_FOUND)
     MESSAGE("Compile with OpenMP/threads")
     SET(CMAKE_CXX_FLAGS "${OpenMP_CXX_FLAGS} ${CMAKE_CXX_FLAGS}")
     SET(CMAKE_C_FLAGS "${OpenMP_C_FLAGS} ${CMAKE_C_FLAGS}")
  ELSE(OPENMP_FOUND)
     MESSAGE("Compiler does not support OpenMP -> no threads")
  ENDIF(OPENMP_FOUND)
ENDIF(THREADS)


################################################################
#convenience stuff
SET(BASIC_SOURCES ${MATHOBJECTFILES} ${ANISOTROPYFILES} ${COSMOSFILES} ${ANALYZEFILES} ${COUPLEDFILES} ${MASSFUNCTIONFILES})
SET(ALL_SOURCES ${WMAPFILES} ${COSMOSFILES} ${MATHOBJECTFILES} ${ANISOTROPYFILES} ${CHAINFILES} ${ANALYZEFILES} ${COUPLEDFILES} ${MASSFUNCTIONFILES})
SET(BASIC_COUPLED_SOURCES ${MATHOBJECTFILES} ${ANISOTROPYFILES} ${COSMOSFILES}
${ANALYZEFILES} ${WMAPFILES} ${COUPLEDFILES} ${MASSFUNCTIONFILES})
SET(BASIC_WMAPSOURCES ${MATHOBJECTFILES} ${ANISOTROPYFILES} ${COSMOSFILES} ${ANALYZEFILES} ${COUPLEDFILES} ${MASSFUNCTIONFILES} ${WMAPFILES})
SET(BASIC_LIBS ${GSL_LIBRARIES} ${WMAPLIB})

################################################################
## maintainer specific stuff, not needed for release
SET(CMBEASY_MAINTAINER_DIR $ENV{CMBEASY_MAINTAINER_DIR} CACHE PATH "path to additional cmbeasy source files that are not part of the official distribution")

SET(CMBEASY_MAINTAINER_MODE TRUE)
INCLUDE(${CMBEASY_MAINTAINER_DIR}/CMakeLists.txt OPTIONAL)
SET(CMBEASY_MAINTAINER_MODE)

################################################################

#################################################

ADD_EXECUTABLE(cmb ${ANISOTROPY}/driver.cc ${BASIC_WMAPSOURCES})
TARGET_LINK_LIBRARIES(cmb ${BASIC_LIBS})

#################################################

ADD_EXECUTABLE(test_mass_function ${ANISOTROPY}/test_mass_function.cc ${BASIC_WMAPSOURCES})
TARGET_LINK_LIBRARIES(test_mass_function ${BASIC_LIBS})


#################################################

ADD_EXECUTABLE(neutrinodriver ${COUPLED}/neutrinodriver.cc ${BASIC_COUPLED_SOURCES})
TARGET_LINK_LIBRARIES(neutrinodriver ${BASIC_LIBS})

#################################################

ADD_EXECUTABLE(coupled_corasaniti ${COUPLED}/coupled_corasaniti.cc ${BASIC_COUPLED_SOURCES})
TARGET_LINK_LIBRARIES(coupled_corasaniti ${BASIC_LIBS})

#################################################

ADD_EXECUTABLE(decayingdedm_driver ${COUPLED}/decayingdedm_driver.cc ${BASIC_COUPLED_SOURCES})
TARGET_LINK_LIBRARIES(decayingdedm_driver ${BASIC_LIBS})

#################################################

ADD_EXECUTABLE(vde ${VDE}/driver_vde.cc ${VDEFILES} ${BASIC_WMAPSOURCES})
TARGET_LINK_LIBRARIES(vde ${BASIC_LIBS})

#################################################

ADD_EXECUTABLE(xcmb ${ANISOTROPY}/xdriver.cc ${BASIC_WMAPSOURCES})
TARGET_LINK_LIBRARIES(xcmb ${BASIC_LIBS})

#################################################
#
ADD_EXECUTABLE(distill ${ANISOTROPY}/distill.cc ${BASIC_SOURCES} ${CHAINFILES})
TARGET_LINK_LIBRARIES(distill ${GSL_LIBRARIES})
SET_TARGET_PROPERTIES(distill PROPERTIES COMPILE_FLAGS -DNOWMAP7)
#
###################################################
##
ADD_EXECUTABLE(testsneroutines ${ANISOTROPY}/testsneroutines.cc ${BASIC_SOURCES})
SET_TARGET_PROPERTIES(testsneroutines PROPERTIES COMPILE_FLAGS -DNOWMAP7)
TARGET_LINK_LIBRARIES(testsneroutines ${BASIC_LIBS})
##
###################################################
##
IF(WMAP7)
  ADD_EXECUTABLE(testwmap ${WMAP}/test.cc ${BASIC_WMAPSOURCES})
  TARGET_LINK_LIBRARIES(testwmap ${WMAPLIB})
ENDIF(WMAP7)
##
###################################################

ADD_LIBRARY(cmbeasy-static STATIC #EXCLUDE_FROM_ALL
                   ${BASIC_WMAPSOURCES})
ADD_LIBRARY(cmbeasy-shared SHARED #EXCLUDE_FROM_ALL
  ${BASIC_WMAPSOURCES})


###################################################
##
IF(BUILD_GUI)
  ADD_DEFINITIONS(${QT_DEFINITIONS})
  ADD_EXECUTABLE(cmbeasy ${GUI}/cmbeasy.cc ${ALL_SOURCES} ${GUIFILES})
  #  SET_TARGET_PROPERTIES(cmbeasy PROPERTIES COMPILE_FLAGS "${QT_DEFINITIONS}")
  TARGET_LINK_LIBRARIES(cmbeasy ${WMAPLIB} ${GSL_LIBRARIES} ${QT_QTCORE_LIBRARY} ${QT_QTGUI_LIBRARY})
ENDIF(BUILD_GUI)

###################################################

ADD_EXECUTABLE(jlgen ${ANISOTROPY}/jlgen.cc ${BASIC_WMAPSOURCES})
TARGET_LINK_LIBRARIES(jlgen ${BASIC_LIBS})

ADD_EXECUTABLE(jlens ${ANISOTROPY}/jlens.cc ${MATHOBJECT}/bessel.cc ${MATHOBJECT}/spline.cc
                     ${MATHOBJECT}/miscmath.cc ${MATHOBJECT}/anchor.cc)


INSTALL(TARGETS cmb xcmb distill testsneroutines jlgen jlens
        RUNTIME DESTINATION bin)

IF(WMAP7)
  INSTALL(TARGETS testwmap
          RUNTIME DESTINATION bin)
ENDIF(WMAP7)

IF(BUILD_GUI)
  INSTALL(TARGETS cmbeasy
          RUNTIME DESTINATION bin)
ENDIF(BUILD_GUI)

CONFIGURE_FILE(${PROJECT_SOURCE_DIR}/Doxyfile.in
                ${PROJECT_BINARY_DIR}/Doxyfile
                @ONLY IMMEDIATE)
ADD_CUSTOM_TARGET(doxy ${DOXYGEN} ${PROJECT_BINARY_DIR}/Doxyfile)

MESSAGE("**********************************************************") 
MESSAGE("Executables will be built in ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}") 
MESSAGE("Use 'make xcmb' to build the example driver,\n 'make cmbeasy' build the graphical user"
        "interface, and 'make mcsimpleexample' for the example montecarlo driver")
MESSAGE("**********************************************************") 

###################  Package targets ###################################

SET(CPACK_PACKAGE_DESCRIPTION_SUMMARY "cmbeasy")
#SET(CPACK_PACKAGE_VENDOR "")
#SET(CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_CURRENT_SOURCE_DIR}/ReadMe.txt")
SET(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/LICENSE_CMBEASY")
SET(CPACK_PACKAGE_VERSION_MAJOR "4")
SET(CPACK_PACKAGE_VERSION_MINOR "4")
SET(CPACK_PACKAGE_VERSION_PATCH "2pre7")
#SET(CPACK_PACKAGE_INSTALL_DIRECTORY "")
#SET(CPACK_STRIP_FILES "bin/MyExecutable")
SET(CPACK_SOURCE_STRIP_FILES "")
SET(CPACK_SOURCE_IGNORE_FILES ".svn" ".git" "resources/jlgen.dat")
SET(CPACK_SOURCE_GENERATOR TBZ2)
SET(CPACK_GENERATOR TBZ2)

INCLUDE(CPack)

