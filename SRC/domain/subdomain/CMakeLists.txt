
target_sources(libg3
    PRIVATE
	Subdomain.cpp
	SubdomainNodIter.cpp 
	ActorSubdomain.cpp
    PUBLIC
	Subdomain.h
	SubdomainNodIter.h 
	ActorSubdomain.h
)

target_sources(domain
    PRIVATE
	ShadowSubdomain.cpp 
    PUBLIC
	ShadowSubdomain.h 
)

target_include_directories(domain PUBLIC $(CMAKE_CURRENT_LIST_DIR))
