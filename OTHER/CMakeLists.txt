cmake_minimum_required(VERSION 3.16)
project(OpenSeesNumPack)

add_subdirectory(AMD/)
add_subdirectory(UMFPACK/)
add_subdirectory(ARPACK/)
add_subdirectory(SuperLU_5.1.1/)
add_subdirectory(CSPARSE/)

