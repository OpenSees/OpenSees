add_library(OPS_Element_feap OBJECT)


target_sources(OPS_Element_feap
    PRIVATE
        dummyFeap.f
        elmt01.f
        elmt02.f
        elmt03.f
        elmt04.f
        elmt05.f
        elmt06.f
        elmt11.f
        fElement.cpp
        fElmt02.cpp
        fElmt05.cpp
        fillCommon.f
        getCommon.f
        #TclFeapElementCommand.cpp
    PUBLIC
        fElement.h
        fElmt02.h
        fElmt05.h
)

opensees_add_cxx_flag(OPS_Element_feap GNU -fpermissive)

target_include_directories(OPS_Element_feap PUBLIC ${CMAKE_CURRENT_LIST_DIR})

