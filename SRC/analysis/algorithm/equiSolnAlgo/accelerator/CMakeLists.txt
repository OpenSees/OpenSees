#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
target_sources(OPS_Analysis
    PRIVATE
      Accelerator.cpp
      MillerAccelerator.cpp 
      RaphsonAccelerator.cpp 
      PeriodicAccelerator.cpp 
      KrylovAccelerator.cpp 
      KrylovAccelerator2.cpp 
      DifferenceAccelerator.cpp 
      DifferenceAccelerator2.cpp
      SecantAccelerator1.cpp 
      SecantAccelerator2.cpp 
      SecantAccelerator3.cpp
    PUBLIC
      Accelerator.h
      MillerAccelerator.h 
      RaphsonAccelerator.h 
      PeriodicAccelerator.h 
      KrylovAccelerator.h 
      KrylovAccelerator2.h 
      DifferenceAccelerator.h 
      DifferenceAccelerator2.h
      SecantAccelerator1.h 
      SecantAccelerator2.h 
      SecantAccelerator3.h
)
# Build files are not properly generated for Windows
# when libraries contain both C++ and Fortran files
add_library(OPS_Analysis_Naccel OBJECT naccel.f)
target_link_libraries(OPS_Analysis OPS_Analysis_Naccel)

target_include_directories(OPS_Analysis PUBLIC ${CMAKE_CURRENT_LIST_DIR})
