cmake_minimum_required(VERSION 3.8...3.21)
if(CMAKE_VERSION VERSION_LESS 3.12) # To support the version range
  cmake_policy(VERSION ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION})
endif()

SET_PROPERTY(GLOBAL PROPERTY TARGET_SUPPORTS_SHARED_LIBS TRUE)

set(CMAKE_CXX_STANDARD 17)
add_compile_options(-pedantic)
add_compile_options(-pedantic-errors)
add_compile_options(-Wextra)
add_compile_options(-Wunused-variable)
add_compile_options(-Wdisabled-optimization)
add_compile_options(-Wall)


# Disabled by default, used to be necessary on the Cedar cluster of Compute Canada
#add_compile_options(-fext-numeric-literals)

#add_compile_options(-Werror)
#add_compile_options(-Weffc++)

option(WARNINGERROR "Consider compilation warnings as errors" OFF)
IF(WARNINGERROR)
    add_compile_options(-Werror)
ENDIF()

find_package(deal.II 9.4.0 QUIET REQUIRED
  HINTS "${deal.II_DIR}" "${DEAL_II_DIR}" "$ENV{DEAL_II_DIR}")
deal_ii_initialize_cached_variables()

#set(PROTOTYPES "FALSE" CACHE BOOL)
option(PROTOTYPES "Enable prototype applications" OFF)
option(MARCH_NATIVE "Enable the -march=native flag" OFF)


Project(lethe)

include(GNUInstallDirs)
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
list(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES
  "${CMAKE_INSTALL_FULL_LIBDIR}" _libdir_index)
if(_libdir_index LESS 0)
  set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_FULL_LIBDIR}")
endif()

ENABLE_TESTING()
ADD_SUBDIRECTORY(source)
ADD_SUBDIRECTORY(applications)
ADD_SUBDIRECTORY(applications_tests)
ADD_SUBDIRECTORY(tests)
IF(PROTOTYPES)
    MESSAGE("Adding prototypes")
    ADD_SUBDIRECTORY(prototypes)
 ENDIF()
 IF(MARCH_NATIVE)
    MESSAGE("  Enabling -march=native flag. This may cause some tests to break.")
    add_compile_options(-march=native)
 ENDIF()


set(missing_dealii_features)

foreach(feat
    DEAL_II_WITH_METIS
    # MPI is required for p4est and Trilinos, but leave it here as
    # additional documentation.
    DEAL_II_WITH_MPI
    DEAL_II_WITH_P4EST
    DEAL_II_WITH_TRILINOS)
  message(STATUS "Checking required ${feat}: " ${${feat}})
  if(NOT DEFINED ${feat} OR NOT ${${feat}})
    list(APPEND missing_dealii_features ${feat})
  endif()
endforeach()

foreach(feat
    DEAL_II_WITH_OPENCASCADE
    DEAL_II_WITH_SUNDIALS)
  message(STATUS "Checking optional ${feat}: " ${${feat}})
endforeach()

list(LENGTH missing_dealii_features num_missing_dealii_features)
if(num_missing_dealii_features GREATER 0)
  message(FATAL_ERROR
    "Lethe requires the following features that deal.II was compiled without: "
    "${missing_dealii_features}")
endif()

unset(num_missing_dealii_features)
unset(missing_dealii_features)

find_program(TXR_EXECUTABLE
  NAMES txr
  DOC "Path to the TXR executable.")
if(TXR_EXECUTABLE STREQUAL "TXR_EXECUTABLE-NOTFOUND")
  message(STATUS "Could not find TXR: omitting target \"update-golden\"")
else()
  message(STATUS "Found TXR: ${TXR_EXECUTABLE}")
  add_custom_target(update-golden
    COMMAND ${TXR_EXECUTABLE}
            ${CMAKE_SOURCE_DIR}/contrib/utilities/update-golden.tl
            ${CMAKE_BINARY_DIR} ${CMAKE_SOURCE_DIR})
endif()
