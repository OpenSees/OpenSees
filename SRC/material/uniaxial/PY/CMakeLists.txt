#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
target_sources(OPS_Material
  PRIVATE
    PyLiq1.cpp
    PySimple1.cpp
    PySimple1Gen.cpp
    PySimple2.cpp
    PySimple3.cpp
	  QzLiq1.cpp
    QzSimple1.cpp
    QzSimple2.cpp
    ShallowFoundationGen.cpp
    TzLiq1.cpp
    TzSimple1.cpp
    TzSimple1Gen.cpp
    TzSimple2.cpp
  PUBLIC
    PyLiq1.h
    PySimple1.h
    PySimple1Gen.h
    PySimple2.h
    PySimple3.h
	  QzLiq1.h
    QzSimple1.h
    QzSimple2.h
    ShallowFoundationGen.h
    TzLiq1.h
    TzSimple1.h
    TzSimple1Gen.h
    TzSimple2.h
#=======
#	PyLiq1.cpp
#	PySimple1.cpp
#	PySimple1Gen.cpp
#	PySimple2.cpp
#	PySimple3.cpp
#	QzLiq1.cpp
#	QzSimple1.cpp
#	QzSimple2.cpp
#	ShallowFoundationGen.cpp
#	TzLiq1.cpp
#	TzSimple1.cpp
#	TzSimple1Gen.cpp
#	TzSimple2.cpp
#    PUBLIC
#	PyLiq1.h
#	PySimple1.h
#	PySimple1Gen.h
#	PySimple2.h
#	PySimple3.h
#	QzLiq1.h
#	QzSimple1.h
#	QzSimple2.h
#	ShallowFoundationGen.h
#	TzLiq1.h
#	TzSimple1.h
#	TzSimple1Gen.h
#	TzSimple2.h
)

target_include_directories(OPS_Material PUBLIC ${CMAKE_CURRENT_LIST_DIR})
