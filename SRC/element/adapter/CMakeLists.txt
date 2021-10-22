#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
add_library(OPS_Element_adapter OBJECT)

target_sources(OPS_Element_adapter
    PRIVATE
        ActuatorCorot.cpp
        Actuator.cpp
        Adapter.cpp
        #TclActuatorCommand.cpp
        #TclActuatorCorotCommand.cpp
        #TclAdapterCommand.cpp
    PUBLIC
        ActuatorCorot.h
        Actuator.h
        Adapter.h
)
target_include_directories(OPS_Element_adapter PUBLIC ${CMAKE_CURRENT_LIST_DIR})

