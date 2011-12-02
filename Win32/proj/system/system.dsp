# Microsoft Developer Studio Project File - Name="system" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Static Library" 0x0104

CFG=system - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "system.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "system.mak" CFG="system - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "system - Win32 Release" (based on "Win32 (x86) Static Library")
!MESSAGE "system - Win32 Debug" (based on "Win32 (x86) Static Library")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
F90=df.exe
RSC=rc.exe

!IF  "$(CFG)" == "system - Win32 Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "Release"
# PROP BASE Intermediate_Dir "Release"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "..\..\lib\release"
# PROP Intermediate_Dir "..\..\obj\system\release"
# PROP Target_Dir ""
# ADD BASE CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_MBCS" /D "_LIB" /YX /FD /c
# ADD CPP /nologo /MLd /W3 /GX /O2 /I "..\..\..\src\system_of_eqn\linearSOE\diagonal" /I "..\..\..\src\system_of_eqn\linearSOE" /I "..\..\..\src\system_of_eqn\linearSOE\sparseSYM" /I "..\..\..\src\analysis\integrator" /I "..\..\..\src\analysis\fe_ele" /I "..\..\..\src\analysis\dof_grp" /I "..\..\..\src\system_of_eqn\eigenSOE" /I "..\..\..\src\handler" /I "..\..\..\symSparse" /I "..\..\..\src\analysis\model\simple" /I "..\..\..\src\system_of_eqn\linearSOE\umfGEN" /I "..\..\..\src\system_of_eqn\linearSOE\fullGEN" /I "..\..\..\src\system_of_eqn\linearSOE\sparseGEN" /I "..\..\..\src\system_of_eqn\linearSOE\bandSPD" /I "..\..\..\src\system_of_eqn\linearSOE\bandGEN" /I "..\..\..\src\domain\domain" /I "..\..\..\src\analysis\model" /I "..\..\..\src\actor\objectBroker" /I "..\..\..\src\actor\channel" /I "..\..\..\src\tagged" /I "..\..\..\src\graph\graph" /I "..\..\..\src\system_of_eqn\linearSOE\profileSPD" /I "..\..\..\src" /I "..\..\..\src\matrix" /I "..\..\..\src\actor\actor" /I "..\..\..\src\system_of_eqn" /I "..\..\..\src\nDarray" /I "..\..\..\src\system_of_eqn\linearSOE\itpack" /I "..\..\..\src\system_of_eqn\linearSOE\cg" /I "..\..\..\other\superlu_3.0\src" /I "c:\Program Files\tcl" /I "c:\Program Files\tcl\include" /D "WIN32" /D "NDEBUG" /D "_LIB" /D "_MBCS" /D "_TCL84" /FD /c
# SUBTRACT CPP /YX
# ADD BASE RSC /l 0x409 /d "NDEBUG"
# ADD RSC /l 0x409 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo

!ELSEIF  "$(CFG)" == "system - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "Debug"
# PROP BASE Intermediate_Dir "Debug"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "..\..\lib\debug"
# PROP Intermediate_Dir "..\..\obj\system\debug"
# PROP Target_Dir ""
# ADD BASE CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_MBCS" /D "_LIB" /YX /FD /GZ /c
# ADD CPP /nologo /W3 /GX /ZI /Od /I "..\..\..\src\system_of_eqn\linearSOE\diagonal" /I "..\..\..\src\system_of_eqn\linearSOE" /I "..\..\..\src\system_of_eqn\linearSOE\sparseSYM" /I "..\..\..\src\analysis\integrator" /I "..\..\..\src\analysis\fe_ele" /I "..\..\..\src\analysis\dof_grp" /I "..\..\..\src\system_of_eqn\eigenSOE" /I "..\..\..\src\handler" /I "..\..\..\symSparse" /I "..\..\..\src\analysis\model\simple" /I "..\..\..\src\system_of_eqn\linearSOE\umfGEN" /I "..\..\..\src\system_of_eqn\linearSOE\fullGEN" /I "..\..\..\src\system_of_eqn\linearSOE\sparseGEN" /I "..\..\..\src\system_of_eqn\linearSOE\bandSPD" /I "..\..\..\src\system_of_eqn\linearSOE\bandGEN" /I "..\..\..\src\domain\domain" /I "..\..\..\src\analysis\model" /I "..\..\..\src\actor\objectBroker" /I "..\..\..\src\actor\channel" /I "..\..\..\src\tagged" /I "..\..\..\src\graph\graph" /I "..\..\..\src\system_of_eqn\linearSOE\profileSPD" /I "..\..\..\src" /I "..\..\..\src\matrix" /I "..\..\..\src\actor\actor" /I "..\..\..\src\system_of_eqn" /I "..\..\..\src\nDarray" /I "..\..\..\src\system_of_eqn\linearSOE\itpack" /I "..\..\..\src\system_of_eqn\linearSOE\cg" /I "..\..\..\other\superlu_3.0\src" /I "c:\Program Files\tcl" /I "c:\Program Files\tcl\include" /D "WIN32" /D "_DEBUG" /D "_LIB" /D "_WIN32" /D "_MBCS" /D "_TCL84" /FR /FD /GZ /c
# SUBTRACT CPP /YX
# ADD BASE RSC /l 0x409 /d "_DEBUG"
# ADD RSC /l 0x409 /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo

!ENDIF 

# Begin Target

# Name "system - Win32 Release"
# Name "system - Win32 Debug"
# Begin Group "Source Files"

# PROP Default_Filter "cpp;c;cxx;rc;def;r;odl;idl;hpj;bat"
# Begin Source File

SOURCE=..\..\..\SRC\system_of_eqn\linearSOE\DomainSolver.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\system_of_eqn\eigenSOE\EigenSOE.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\system_of_eqn\eigenSOE\EigenSolver.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\system_of_eqn\linearSOE\LinearSOE.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\system_of_eqn\linearSOE\LinearSOESolver.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\system_of_eqn\Solver.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\system_of_eqn\SystemOfEqn.cpp
# End Source File
# End Group
# Begin Group "Header Files"

# PROP Default_Filter "h;hpp;hxx;hm;inl"
# Begin Source File

SOURCE=..\..\..\SRC\system_of_eqn\linearSOE\DomainSolver.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\system_of_eqn\eigenSOE\EigenSOE.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\system_of_eqn\eigenSOE\EigenSolver.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\system_of_eqn\linearSOE\LinearSOE.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\system_of_eqn\linearSOE\LinearSOESolver.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\system_of_eqn\Solver.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\system_of_eqn\SystemOfEqn.h
# End Source File
# End Group
# Begin Group "profileSPD"

# PROP Default_Filter ""
# Begin Source File

SOURCE=..\..\..\SRC\system_of_eqn\linearSOE\profileSPD\ProfileSPDLinDirectSolver.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\system_of_eqn\linearSOE\profileSPD\ProfileSPDLinDirectSolver.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\system_of_eqn\linearSOE\profileSPD\ProfileSPDLinSOE.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\system_of_eqn\linearSOE\profileSPD\ProfileSPDLinSOE.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\system_of_eqn\linearSOE\profileSPD\ProfileSPDLinSolver.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\system_of_eqn\linearSOE\profileSPD\ProfileSPDLinSolver.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\system_of_eqn\linearSOE\profileSPD\ProfileSPDLinSubstrSolver.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\system_of_eqn\linearSOE\profileSPD\ProfileSPDLinSubstrSolver.h
# End Source File
# End Group
# Begin Group "symSparse"

# PROP Default_Filter ""
# Begin Source File

SOURCE=..\..\..\SRC\system_of_eqn\linearSOE\sparseSYM\FeStructs.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\system_of_eqn\linearSOE\sparseSYM\globalVars.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\system_of_eqn\linearSOE\sparseSYM\grcm.c
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\system_of_eqn\linearSOE\sparseSYM\nest.c
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\system_of_eqn\linearSOE\sparseSYM\newordr.c
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\system_of_eqn\linearSOE\sparseSYM\nmat.c
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\system_of_eqn\linearSOE\sparseSYM\nnsim.c
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\system_of_eqn\linearSOE\sparseSYM\symbolic.c
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\system_of_eqn\linearSOE\sparseSYM\SymSparseLinSOE.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\system_of_eqn\linearSOE\sparseSYM\SymSparseLinSOE.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\system_of_eqn\linearSOE\sparseSYM\SymSparseLinSolver.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\system_of_eqn\linearSOE\sparseSYM\SymSparseLinSolver.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\system_of_eqn\linearSOE\sparseSYM\tim.c
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\system_of_eqn\linearSOE\sparseSYM\tim.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\system_of_eqn\linearSOE\sparseSYM\utility.c
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\system_of_eqn\linearSOE\sparseSYM\utility.h
# End Source File
# End Group
# Begin Group "bandGEN"

# PROP Default_Filter ""
# Begin Source File

SOURCE=..\..\..\SRC\system_of_eqn\linearSOE\bandGEN\BandGenLinLapackSolver.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\system_of_eqn\linearSOE\bandGEN\BandGenLinLapackSolver.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\system_of_eqn\linearSOE\bandGEN\BandGenLinSOE.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\system_of_eqn\linearSOE\bandGEN\BandGenLinSOE.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\system_of_eqn\linearSOE\bandGEN\BandGenLinSolver.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\system_of_eqn\linearSOE\bandGEN\BandGenLinSolver.h
# End Source File
# End Group
# Begin Group "bandSPD"

# PROP Default_Filter ""
# Begin Source File

SOURCE=..\..\..\SRC\system_of_eqn\linearSOE\bandSPD\BandSPDLinLapackSolver.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\system_of_eqn\linearSOE\bandSPD\BandSPDLinLapackSolver.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\system_of_eqn\linearSOE\bandSPD\BandSPDLinSOE.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\system_of_eqn\linearSOE\bandSPD\BandSPDLinSOE.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\system_of_eqn\linearSOE\bandSPD\BandSPDLinSolver.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\system_of_eqn\linearSOE\bandSPD\BandSPDLinSolver.h
# End Source File
# End Group
# Begin Group "fullGEN"

# PROP Default_Filter ""
# Begin Source File

SOURCE=..\..\..\SRC\system_of_eqn\linearSOE\fullGEN\FullGenLinLapackSolver.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\system_of_eqn\linearSOE\fullGEN\FullGenLinLapackSolver.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\system_of_eqn\linearSOE\fullGEN\FullGenLinSOE.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\system_of_eqn\linearSOE\fullGEN\FullGenLinSOE.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\system_of_eqn\linearSOE\fullGEN\FullGenLinSolver.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\system_of_eqn\linearSOE\fullGEN\FullGenLinSolver.h
# End Source File
# End Group
# Begin Group "umfpack"

# PROP Default_Filter ".cpp; .h"
# Begin Source File

SOURCE=..\..\..\SRC\system_of_eqn\linearSOE\umfGEN\UmfpackGenLinSOE.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\system_of_eqn\linearSOE\umfGEN\UmfpackGenLinSOE.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\system_of_eqn\linearSOE\umfGEN\UmfpackGenLinSolver.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\system_of_eqn\linearSOE\umfGEN\UmfpackGenLinSolver.h
# End Source File
# End Group
# Begin Group "superLU"

# PROP Default_Filter ".cpp; .h"
# Begin Source File

SOURCE=..\..\..\SRC\system_of_eqn\linearSOE\sparseGEN\SparseGenColLinSOE.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\system_of_eqn\linearSOE\sparseGEN\SparseGenColLinSOE.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\system_of_eqn\linearSOE\sparseGEN\SparseGenColLinSolver.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\system_of_eqn\linearSOE\sparseGEN\SparseGenColLinSolver.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\system_of_eqn\linearSOE\sparseGEN\SuperLU.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\system_of_eqn\linearSOE\sparseGEN\SuperLU.h
# End Source File
# End Group
# Begin Group "arpack"

# PROP Default_Filter ""
# Begin Source File

SOURCE=..\..\..\SRC\system_of_eqn\eigenSOE\BandArpackSOE.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\system_of_eqn\eigenSOE\BandArpackSOE.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\system_of_eqn\eigenSOE\BandArpackSolver.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\system_of_eqn\eigenSOE\BandArpackSolver.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\system_of_eqn\eigenSOE\SymArpackSOE.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\system_of_eqn\eigenSOE\SymArpackSOE.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\system_of_eqn\eigenSOE\SymArpackSolver.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\system_of_eqn\eigenSOE\SymArpackSolver.h
# End Source File
# End Group
# Begin Group "symBandEigen"

# PROP Default_Filter ""
# Begin Source File

SOURCE=..\..\..\SRC\system_of_eqn\eigenSOE\SymBandEigenSOE.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\system_of_eqn\eigenSOE\SymBandEigenSOE.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\system_of_eqn\eigenSOE\SymBandEigenSolver.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\system_of_eqn\eigenSOE\SymBandEigenSolver.h
# End Source File
# End Group
# Begin Group "diagonal"

# PROP Default_Filter ""
# Begin Source File

SOURCE=..\..\..\SRC\system_of_eqn\linearSOE\diagonal\DiagonalDirectSolver.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\system_of_eqn\linearSOE\diagonal\DiagonalDirectSolver.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\system_of_eqn\linearSOE\diagonal\DiagonalSOE.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\system_of_eqn\linearSOE\diagonal\DiagonalSOE.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\system_of_eqn\linearSOE\diagonal\DiagonalSolver.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\system_of_eqn\linearSOE\diagonal\DiagonalSolver.h
# End Source File
# End Group
# End Target
# End Project
