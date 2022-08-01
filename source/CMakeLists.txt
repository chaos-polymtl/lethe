MESSAGE(STATUS "Setting up Lethe library")
FIND_PACKAGE(deal.II 9.2.0 QUIET
  HINTS ${deal.II_DIR} ${DEAL_II_DIR} ../ ../../ $ENV{DEAL_II_DIR}
  )
IF(NOT ${deal.II_FOUND})
  MESSAGE(FATAL_ERROR "\n"
    "*** Could not locate a (sufficiently recent) version of deal.II. ***\n\n"
    "You may want to either pass a flag -DDEAL_II_DIR=/path/to/deal.II to cmake\n"
    "or set an environment variable \"DEAL_II_DIR\" that contains this path."
    )
ENDIF()

INCLUDE_DIRECTORIES(
  ${CMAKE_SOURCE_DIR}/include/
  SYSTEM
  )

ADD_SUBDIRECTORY(core)
ADD_SUBDIRECTORY(solvers)
ADD_SUBDIRECTORY(dem)
ADD_SUBDIRECTORY(rpt)
ADD_SUBDIRECTORY(fem-dem)

