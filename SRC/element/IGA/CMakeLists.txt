target_sources(OPS_Element
    PRIVATE
    IGAKLShell.cpp
    IGAKLShell_BendingStrip.cpp
    PUBLIC
    IGAKLShell.h
    IGAKLShell_BendingStrip.h
)

target_include_directories(OPS_Element PUBLIC ${CMAKE_CURRENT_LIST_DIR})
