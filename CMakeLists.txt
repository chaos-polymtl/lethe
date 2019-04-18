# top level CMakeLists.txt
CMAKE_MINIMUM_REQUIRED(VERSION 2.8.12)
SET_PROPERTY(GLOBAL PROPERTY TARGET_SUPPORTS_SHARED_LIBS TRUE)

FIND_PACKAGE(deal.II 9.0.0 QUIET
  HINTS ${deal.II_DIR} ${DEAL_II_DIR} ../ ../../ $ENV{DEAL_II_DIR}
  )
IF(NOT ${deal.II_FOUND})
  MESSAGE(FATAL_ERROR "\n"
    "*** Could not locate a (sufficiently recent) version of deal.II. ***\n\n"
    "You may want to either pass a flag -DDEAL_II_DIR=/path/to/deal.II to cmake\n"
    "or set an environment variable \"DEAL_II_DIR\" that contains this path."
    )
ENDIF()

DEAL_II_INITIALIZE_CACHED_VARIABLES()

Project(lethe)
ENABLE_TESTING()
ADD_SUBDIRECTORY(src)
ADD_SUBDIRECTORY(applications)
ADD_SUBDIRECTORY(tests/taylorcouette)
ADD_SUBDIRECTORY(tests/taylorcouette_unstructured)
ADD_SUBDIRECTORY(tests/mms2d)
ADD_SUBDIRECTORY(tests/mms2d_unstructured)
ADD_SUBDIRECTORY(tests/mms3d)




