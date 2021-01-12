
target_sources(domain
    PRIVATE
	Subdomain.cpp
	SubdomainNodIter.cpp 
	ShadowSubdomain.cpp 
	ActorSubdomain.cpp
    PUBLIC
	Subdomain.h
	SubdomainNodIter.h 
	ShadowSubdomain.h 
	ActorSubdomain.h
)

target_include_directories(domain PUBLIC $(CMAKE_CURRENT_LIST_DIR))
