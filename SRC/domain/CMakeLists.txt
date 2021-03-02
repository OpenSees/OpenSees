
add_library(domain "")

target_include_directories(domain PUBLIC $(CMAKE_CURRENT_LIST_DIR))

add_subdirectory(component)
add_subdirectory(domain)
add_subdirectory(load)	
add_subdirectory(node)	
add_subdirectory(subdomain)	
add_subdirectory(constraints)
add_subdirectory(pattern)
add_subdirectory(groundMotion)
add_subdirectory(region)
add_subdirectory(partitioner)

