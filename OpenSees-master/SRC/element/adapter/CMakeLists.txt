#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
#add_library(OPS_Element_adapter OBJECT)

target_sources(OPS_Element
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
target_include_directories(OPS_Element PUBLIC ${CMAKE_CURRENT_LIST_DIR})

