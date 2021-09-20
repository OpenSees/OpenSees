
opensees_elements(OPS_Element_beam3d
  SOURCES
     beam3d01.cpp
     beam3d02.cpp
     beam3d01.h
     beam3d02.h
)

target_include_directories(OPS_Element_beam3d PUBLIC ${CMAKE_CURRENT_LIST_DIR})

