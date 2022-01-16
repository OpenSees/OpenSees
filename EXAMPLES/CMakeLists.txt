#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================

# Verification tests
#-------------------------------------------------------------------------
add_subdirectory(verification)

# ModelBuilder tests
#-------------------------------------------------------------------------
add_subdirectory(PlaneFrame)
add_subdirectory(Example1)

# Serialization tests
#-------------------------------------------------------------------------
add_executable(database_test Serialization/database.cpp)
target_link_libraries(database_test G3 G3_API)
add_test(SerializationTest database_test COMMAND database_test ./tmp/database/test1)


