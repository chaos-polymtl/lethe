# top level CMakeLists.txt
CMAKE_MINIMUM_REQUIRED(VERSION 2.8.12)
SET_PROPERTY(GLOBAL PROPERTY TARGET_SUPPORTS_SHARED_LIBS TRUE)

add_compile_options(-pedantic)
add_compile_options(-pedantic-errors)
add_compile_options(-Wextra)
add_compile_options(-Wunused-variable)
add_compile_options(-Wdisabled-optimization)
#add_compile_options(-Werror)
#add_compile_options(-Weffc++)

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

IF(EXISTS ${CMAKE_SOURCE_DIR}/bundled/CMakeLists.txt)
  SET(DEAL_II_HAVE_BUNDLED_DIRECTORY TRUE)
ENDIF()

IF(EXISTS ${CMAKE_SOURCE_DIR}/doc/CMakeLists.txt)
  SET(DEAL_II_HAVE_DOC_DIRECTORY TRUE)
ENDIF()

IF(EXISTS ${CMAKE_SOURCE_DIR}/tests/CMakeLists.txt)
  SET(DEAL_II_HAVE_TESTS_DIRECTORY TRUE)
ENDIF()

set(PROTOTYPES "FALSE" CACHE FORCE BOOL)


Project(lethe)
SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
ENABLE_TESTING()
ADD_SUBDIRECTORY(source)
ADD_SUBDIRECTORY(applications)
ADD_SUBDIRECTORY(applications_tests)
ADD_SUBDIRECTORY(tests)
IF(PROTOTYPES)
    ADD_SUBDIRECTORY(prototypes)
ENDIF()





