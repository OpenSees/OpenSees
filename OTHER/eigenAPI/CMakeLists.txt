#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
# target_sources(EIGENAPI
#     PUBLIC
#       EigenAPI.h
#       converters.h
#       operations.h
#       typedefs.h
# )

target_include_directories(EIGENAPI PUBLIC ${CMAKE_CURRENT_LIST_DIR})
target_include_directories(EIGENAPI PUBLIC ${CMAKE_CURRENT_LIST_DIR}/eigen)
