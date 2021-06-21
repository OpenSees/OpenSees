target_sources(LG3_Analysis
    PRIVATE
      DOF_Group.cpp
      LagrangeDOF_Group.cpp
      TransformationDOF_Group.cpp
    PUBLIC
      DOF_Group.h
      LagrangeDOF_Group.h
      TransformationDOF_Group.h
)

target_include_directories(LG3_Analysis PUBLIC $(CMAKE_CURRENT_LIST_DIR))
