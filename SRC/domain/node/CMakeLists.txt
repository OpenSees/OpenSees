
target_sources(G3_Domain
    PRIVATE
    Node.cpp
    NodalLoad.cpp
    PUBLIC
    Node.h
    NodalLoad.h
)

target_include_directories(OPS_Domain PUBLIC ${CMAKE_CURRENT_LIST_DIR})
