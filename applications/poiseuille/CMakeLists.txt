# Usually, you will not need to modify anything beyond this point...
#Method of manufactured solution 2D

# Set the name of the project and target:
SET(TARGET "poiseuille")

INCLUDE_DIRECTORIES(
  lethe
  ${CMAKE_SOURCE_DIR}/src/
  )
ADD_EXECUTABLE(poiseuille poiseuille.cc)
DEAL_II_SETUP_TARGET(poiseuille)
TARGET_LINK_LIBRARIES(poiseuille lethe)

install(TARGETS poiseuille RUNTIME DESTINATION bin)

#ENABLE_TESTING()
#ADD_SUBDIRECTORY(tests)
