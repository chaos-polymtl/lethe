# Set the name of the project and target:
SET(TARGET "chorin_navier_stokes")

ADD_EXECUTABLE(${TARGET} chorin_navier_stokes.cc)
DEAL_II_SETUP_TARGET(${TARGET})
TARGET_LINK_LIBRARIES(${TARGET}) 
