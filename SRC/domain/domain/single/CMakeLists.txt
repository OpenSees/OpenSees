
target_sources(domain
    PRIVATE
	SingleDomEleIter.cpp 
	SingleDomNodIter.cpp
	SingleDomSP_Iter.cpp 
	SingleDomMP_Iter.cpp 
	SingleDomAllSP_Iter.cpp
	SingleDomParamIter.cpp  
	SingleDomPC_Iter.cpp
    PUBLIC
	SingleDomEleIter.h 
	SingleDomNodIter.h
	SingleDomSP_Iter.h 
	SingleDomMP_Iter.h 
	SingleDomAllSP_Iter.h
	SingleDomParamIter.h  
	SingleDomPC_Iter.h
)

target_include_directories(domain PUBLIC $(CMAKE_CURRENT_LIST_DIR))
