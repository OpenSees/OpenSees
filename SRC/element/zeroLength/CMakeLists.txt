target_sources(element
    PRIVATE
	ZeroLength.cpp
	ZeroLengthSection.cpp
	ZeroLengthContact2D.cpp
	ZeroLengthContact3D.cpp
	ZeroLengthND.cpp
	ZeroLengthContactNTS2D.cpp
	ZeroLengthInterface2D.cpp
	ZeroLengthRocking.cpp
	CoupledZeroLength.cpp
	ZeroLengthImpact3D.cpp
	ZeroLengthVG_HG.cpp
    PUBLIC
	ZeroLength.h
	ZeroLengthSection.h
	ZeroLengthContact2D.h
	ZeroLengthContact3D.h
	ZeroLengthND.h
	ZeroLengthContactNTS2D.h
	ZeroLengthInterface2D.h
	ZeroLengthRocking.h
	CoupledZeroLength.h
	ZeroLengthImpact3D.h
	ZeroLengthVG_HG.h
)

target_include_directories(element PUBLIC $(CMAKE_CURRENT_LIST_DIR))
