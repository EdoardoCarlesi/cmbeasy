SET(CFITSIO_DIR "" CACHE PATH "Path to the cfitsio installation")

FIND_PATH(CFITSIO_INCLUDE_DIR fitsio.h
                                 ${CFITSIO_DIR}/include
                                 /usr/include/
                                 /usr/local/include/)

FIND_LIBRARY(CFITSIO_LIBRARY NAMES cfitsio
                                PATHS ${CFITSIO_DIR}/lib
                                     /usr/lib /usr/local/lib)
IF (CFITSIO_INCLUDE_DIR AND CFITSIO_LIBRARY)

   SET(CFITSIO_FOUND TRUE)
ENDIF (CFITSIO_INCLUDE_DIR AND CFITSIO_LIBRARY)

IF (CFITSIO_FOUND)
   #Set HPPATH first to the lib path
   GET_FILENAME_COMPONENT(CFPATH ${CFITSIO_LIBRARY} PATH)
   #and then its parent dir
   GET_FILENAME_COMPONENT(CFPATH ${CFPATH} PATH)
   SET(CFITSIO_DIR ${CFPATH} CACHE PATH "Path to the cfitsio installation" FORCE)
   IF (NOT Cfitsio_FIND_QUIETLY)
      MESSAGE(STATUS "Found cfitsio: ${CFITSIO_LIBRARY}")
   ENDIF (NOT Cfitsio_FIND_QUIETLY)
ELSE (CFITSIO_FOUND)
   IF (Cfitsio_FIND_REQUIRED)
      MESSAGE(FATAL_ERROR "Could not find cfitsio")
   ENDIF (Cfitsio_FIND_REQUIRED)
ENDIF (CFITSIO_FOUND)

MARK_AS_ADVANCED(CFITSIO_INCLUDE_DIR CFITSIO_LIBRARY_DIR CFITSIO_LIBRARY)
