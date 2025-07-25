# These functions are taken from Aspect.
# A function that checks if a test depends on another test.
# This is encoded in .prm files through lines of the form
#    '# DEPENDS-ON: testname'
# The result is returned in a variable _depends_on in the
# caller's scope. Having test dependencies is helpful if
# one test requires another test to finish first, for example
# because the result of the first test is used by the later test.
function(get_depends_on _filename)
    file(STRINGS ${_filename} _input_lines
            REGEX "DEPENDS-ON:")
    if("${_input_lines}" STREQUAL "")
        set(_depends_on "" PARENT_SCOPE)
    else()
        # go over the (possibly multiple) lines with DEPENDS-ON markers and choose the last
        foreach(_input_line ${_input_lines})
            set(_last_line ${_input_line})
        endforeach()
        string(REGEX REPLACE "^ *# *DEPENDS-ON: *(.*) *$" "\\1"
                _depends_on ${_last_line})
        set(_depends_on "${_depends_on}" PARENT_SCOPE)
    endif()
endfunction()


# A function that extracts from a file (presumably a .output file)
# the number of MPI processes this test is to be invoked as.
# This is encoded in the name of the output file as mpirun=N
# The result is returned in a variable
function(get_mpi_count _filename)
    string(REGEX MATCH ".mpirun=[0-9]+" _match "${_filename}")
    if(_match)
        set(_mpirun_string "${_match}" PARENT_SCOPE)
    else()
        set(_mpirun_string "" PARENT_SCOPE)
    endif()
endfunction()


