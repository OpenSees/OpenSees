target_sources(G3_Analysis
    PRIVATE
      ConstraintHandler.cpp 
      PlainHandler.cpp
      PenaltyConstraintHandler.cpp 
      LagrangeConstraintHandler.cpp
      TransformationConstraintHandler.cpp 
    PUBLIC
      ConstraintHandler.h 
      PlainHandler.h
      PenaltyConstraintHandler.h 
      LagrangeConstraintHandler.h
      TransformationConstraintHandler.h 
)

target_include_directories(G3_Analysis PUBLIC ${CMAKE_CURRENT_LIST_DIR})
