#add_test(OPS_Serialize)
enable_testing()
add_executable(database_test Serialization/database.cpp)
target_link_libraries(database_test G3 G3_API)


# Add tests
add_test(SerializationTest database_test COMMAND database_test ./tmp/database/test1)

