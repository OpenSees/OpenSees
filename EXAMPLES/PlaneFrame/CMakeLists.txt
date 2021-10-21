
add_executable(planeFrame main.cpp)
target_link_libraries(planeFrame  OPS_ModelBuilder_PlaneFrame G3_API G3)
add_test(PlaneFrameTest planeFrame COMMAND planeFrame)
