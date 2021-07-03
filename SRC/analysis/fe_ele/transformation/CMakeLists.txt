target_sources(analysis
    PRIVATE
      TransformationFE.cpp
    PUBLIC
      TransformationFE.h
)

target_include_directories(analysis PUBLIC $(CMAKE_CURRENT_LIST_DIR))
