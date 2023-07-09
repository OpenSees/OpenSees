#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================

#opensees_elements(OPS_Element_feap
target_sources(OPS_Element
  PRIVATE
    fElement.cpp
    fElmt02.cpp
    fElmt05.cpp
  PUBLIC
    fElement.h
    fElmt02.h
    fElmt05.h
)

target_sources(OPS_ElementFortran
    PRIVATE
        dummyFeap.f
        elmt01.f
        elmt02.f
        elmt03.f
        elmt04.f
        elmt05.f
        elmt06.f
        elmt11.f
        fillCommon.f
        getCommon.f
)


#opensees_add_cxx_flag(OPS_Element GNU -fpermissive)

target_include_directories(OPS_Element PUBLIC ${CMAKE_CURRENT_LIST_DIR})


