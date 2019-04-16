# Usually, you will not need to modify anything beyond this point...
#Method of manufactured solution 2D

# Set the name of the project and target:
SET(TARGET "taylorcouette_unstructured")

INCLUDE_DIRECTORIES(
  lethe
  ${CMAKE_SOURCE_DIR}/src/
  )
ADD_EXECUTABLE(taylorcouette_unstructured taylorcouette_unstructured.cc)
DEAL_II_SETUP_TARGET(taylorcouette_unstructured)
TARGET_LINK_LIBRARIES(taylorcouette_unstructured lethe)


#ENABLE_TESTING()
#ADD_SUBDIRECTORY(tests)
