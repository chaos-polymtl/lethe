# top level CMakeLists.txt

ENABLE_TESTING()
ADD_SUBDIRECTORY(initial_conditions)
ADD_SUBDIRECTORY(mms2d)
ADD_SUBDIRECTORY(mms3d)
ADD_SUBDIRECTORY(mms2d_unstructured)
ADD_SUBDIRECTORY(poiseuille)
ADD_SUBDIRECTORY(cylinder)
ADD_SUBDIRECTORY(taylorcouette)
ADD_SUBDIRECTORY(taylorcouette_unstructured)
ADD_SUBDIRECTORY(tg_vortex)
