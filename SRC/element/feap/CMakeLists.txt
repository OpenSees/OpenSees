
opensees_elements(OPS_Element_feap
  SOURCES
    fElement.cpp
    fElmt02.cpp
    fElmt05.cpp
    fElement.h
    fElmt02.h
    fElmt05.h
)

add_library(OPS_Element_feap_f)
target_sources(OPS_Element_feap_f
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
        #TclFeapElementCommand.cpp
)

opensees_add_cxx_flag(OPS_Element_feap GNU -fpermissive)
target_link_libraries(OPS_Element_feap PRIVATE OPS_Element_feap_f)

target_include_directories(OPS_Element_feap PUBLIC ${CMAKE_CURRENT_LIST_DIR})


