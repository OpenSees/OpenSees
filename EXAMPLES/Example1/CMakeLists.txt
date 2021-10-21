
add_executable(trussExample main.cpp)
target_link_libraries(trussExample OPS_Element_truss G3_API G3)
add_test(TrussTest trussExample COMMAND trussExample)

