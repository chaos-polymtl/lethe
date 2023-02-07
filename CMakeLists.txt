cmake_minimum_required(VERSION 3.8...3.21)
if(CMAKE_VERSION VERSION_LESS 3.12) # To support the version range
  cmake_policy(VERSION ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION})
endif()
# TARGET_SUPPORTS_SHARED_LIBS is set by PROJECT, but we have to call
# FIND_PACKAGE before that, and deal.II's configuration calls
# ADD_LIBRARY with the SHARED keyword, so setting
# TARGET_SUPPORTS_SHARED_LIBS beforehand may be necessary in some
# environments.
set_property(GLOBAL PROPERTY TARGET_SUPPORTS_SHARED_LIBS TRUE)

# Assume that deal.II_DIR won't be provided on the command line in favor
# of DEAL_II_DIR and can thus be used to detect whether CMake was
# already run.
if(DEFINED deal.II_DIR)
  set(_had_deal_ii TRUE)
endif()
find_package(deal.II 9.4.0 QUIET REQUIRED
  HINTS "${deal.II_DIR}" "${DEAL_II_DIR}" "$ENV{DEAL_II_DIR}")
deal_ii_initialize_cached_variables()

project(lethe VERSION 0.1 LANGUAGES CXX)

# Check for deal.II features after the PROJECT call — but still as close
# as possible to the FIND_PACKAGE call — so that the messages appear
# after the CXX compiler–related messages, which is where FIND_PACKAGE
# is usually invoked (i.e., after PROJECT).
if(NOT _had_deal_ii)
  set(_missing_deal_ii_features)
  foreach(_feat IN ITEMS
      DEAL_II_WITH_METIS
      # MPI is required for p4est and Trilinos, but leave it here as
      # additional documentation.
      DEAL_II_WITH_MPI
      DEAL_II_WITH_P4EST
      DEAL_II_WITH_TRILINOS)
    message(STATUS "Checking required ${_feat}: " ${${_feat}})
    if(NOT DEFINED ${_feat} OR NOT ${${_feat}})
      list(APPEND _missing_deal_ii_features ${_feat})
    endif()
  endforeach()

  foreach(_feat IN ITEMS
      DEAL_II_WITH_OPENCASCADE
      DEAL_II_WITH_SUNDIALS)
    message(STATUS "Checking optional ${_feat}: " ${${_feat}})
  endforeach()

  list(LENGTH _missing_deal_ii_features _num_missing_deal_ii_features)
  if(_num_missing_deal_ii_features GREATER 0)
    message(FATAL_ERROR
      "Lethe requires the following features that deal.II was compiled without: "
      "${_missing_deal_ii_features}")
  endif()

  unset(_num_missing_deal_ii_features)
  unset(_missing_deal_ii_features)
endif()
unset(_had_deal_ii)

if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
  set(CMAKE_CXX_STANDARD 17 CACHE STRING "The C++ standard to use")
  set(CMAKE_CXX_STANDARD_REQUIRED ON)
  set(CMAKE_CXX_EXTENSIONS OFF)

  add_compile_options(
    -Wall -Wextra -Wpedantic # Actually already in DEAL_II_CXX_FLAGS
    -Wunused-variable -Wdisabled-optimization)

  if(CMAKE_VERSION VERSION_LESS 3.24)
    option(CMAKE_COMPILE_WARNING_AS_ERROR "Treat warnings on compile as errors")
    if(CMAKE_COMPILE_WARNING_AS_ERROR)
      add_compile_options(-Werror)
    endif()
  endif()
  if(CMAKE_COMPILE_WARNING_AS_ERROR)
    add_compile_options(-pedantic-errors) # Probably not added by CMake
  endif()

  # Tests may break when enabled.
  option(ENABLE_NATIVE_OPTIMIZATION "Enable native optimization (-march=native)")
  if(ENABLE_NATIVE_OPTIMIZATION)
    add_compile_options(-march=native)
  endif()

  option(BUILD_PROTOTYPES "Build the prototype applications")
  mark_as_advanced(BUILD_PROTOTYPES)

  include(GNUInstallDirs)
  set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
  list(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES
    "${CMAKE_INSTALL_FULL_LIBDIR}" _libdir_index)
  if(_libdir_index LESS 0)
    set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_FULL_LIBDIR}")
  endif()

  option(BUILD_TESTING "Build the tests" ON)
  if(BUILD_TESTING)
    enable_testing()

    if(DEFINED TXR_EXECUTABLE)
      set(_had_txr_executable TRUE)
    endif()
    find_program(TXR_EXECUTABLE NAMES txr
      DOC "Path to the TXR executable.")
    mark_as_advanced(TXR_EXECUTABLE)
    # Print the status information only if TXR_EXECUTABLE was neither
    # provided on the command line nor found in the cache.
    # Ideally, users would be informed whenever TXR_EXECUTABLE changes
    # (which is what happens with PERL_EXECUTABLE), but there is no way
    # to access a variable's cached value if it is also provided on the
    # command line (except by defining additional variables).
    if(NOT _had_txr_executable)
      if(TXR_EXECUTABLE)
        message(STATUS "Found TXR: ${TXR_EXECUTABLE}")
      else()
        message(STATUS "Could not find TXR: omitting target \"update-golden\"")
      endif()
    endif()
    unset(_had_txr_executable)
  endif()
endif()

add_subdirectory(source)
add_subdirectory(applications)

if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
  if(BUILD_PROTOTYPES)
    add_subdirectory(prototypes)
  endif()

  if(BUILD_TESTING)
    add_subdirectory(tests)
    add_subdirectory(applications_tests)

    if(TXR_EXECUTABLE)
      add_custom_target(update-golden
        COMMAND "${TXR_EXECUTABLE}"
                "${CMAKE_SOURCE_DIR}/contrib/utilities/update-golden.tl"
                "${CMAKE_BINARY_DIR}" "${CMAKE_SOURCE_DIR}")
    endif()
  endif()
endif()
