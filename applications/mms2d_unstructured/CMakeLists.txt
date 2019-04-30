# Usually, you will not need to modify anything beyond this point...
#Method of manufactured solution 2D

# Set the name of the project and target:
SET(TARGET "mms2d_unstructured")

INCLUDE_DIRECTORIES(
  lethe
  ${CMAKE_SOURCE_DIR}/src/
  )

ADD_EXECUTABLE(mms2d_unstructured mms2d_unstructured.cc)
DEAL_II_SETUP_TARGET(mms2d_unstructured)
TARGET_LINK_LIBRARIES(mms2d_unstructured lethe)

install(TARGETS mms2d_unstructured RUNTIME DESTINATION bin)
