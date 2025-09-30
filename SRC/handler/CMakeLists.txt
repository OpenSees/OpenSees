#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
#add_library(handler OBJECT)

target_sources(OPS_Handler
    PRIVATE
	    BinaryFileStream.cpp
		ChannelStream.cpp
		DatabaseStream.cpp
		DataFileStream.cpp
		DataFileStreamAdd.cpp
		DummyStream.cpp
		FileStream.cpp
        OPS_Stream.cpp
        StandardStream.cpp
        TCP_Stream.cpp
        XmlFileStream.cpp
    PUBLIC
	    BinaryFileStream.h
		ChannelStream.h
		DatabaseStream.h
		DataFileStream.h
		DataFileStreamAdd.h
		DummyStream.h
		FileStream.h
        OPS_Stream.h
        StandardStream.h
        TCP_Stream.h
        XmlFileStream.h
)

target_include_directories(OPS_Handler PUBLIC ${CMAKE_CURRENT_LIST_DIR})
