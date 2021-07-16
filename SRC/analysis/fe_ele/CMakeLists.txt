target_sources(G3_Analysis
    PRIVATE
      FE_Element.cpp
    PUBLIC
      FE_Element.h
)

#target_include_directories(G3_Analysis PUBLIC ${CMAKE_CURRENT_LIST_DIR})

add_subdirectory(penalty)
add_subdirectory(lagrange)
add_subdirectory(transformation)
