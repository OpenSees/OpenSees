# Microsoft Developer Studio Project File - Name="quickMain" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Console Application" 0x0103

CFG=quickMain - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "quickMain.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "quickMain.mak" CFG="quickMain - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "quickMain - Win32 Release" (based on "Win32 (x86) Console Application")
!MESSAGE "quickMain - Win32 Debug" (based on "Win32 (x86) Console Application")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
F90=df.exe
RSC=rc.exe

!IF  "$(CFG)" == "quickMain - Win32 Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "Release"
# PROP BASE Intermediate_Dir "Release"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "..\..\bin"
# PROP Intermediate_Dir "..\..\obj\quickMain\release"
# PROP Ignore_Export_Lib 0
# PROP Target_Dir ""
# ADD BASE CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /c
# ADD CPP /nologo /MLd /W3 /GX /O2 /I "..\..\..\src\ndarray" /I "..\..\..\src\domain\groundMotion" /I "..\..\..\src\element\nonlinearBeamColumn\quadrule" /I "..\..\..\src\material\section\fiber" /I "..\..\..\src\element\fourNodeQuad" /I "..\..\..\src\element\nonlinearBeamColumn\element" /I "..\..\..\src\element\zeroLength" /I "..\..\..\src\element\beamWithHinges" /I "..\..\..\src\recorder\response" /I "..\..\..\src\coordTransformation" /I "..\..\..\src\material\section\repres\section" /I "..\..\..\src\analysis\algorithm\equiSolnAlgo" /I "..\..\..\src\system_of_eqn\eigenSOE" /I "..\..\..\src\analysis\algorithm\eigenAlgo" /I "..\..\..\src\material\nD" /I "..\..\..\src\material\uniaxial" /I "..\..\..\src\handler" /I "..\..\..\src\tcl" /I "..\..\..\src\actor\objectBroker" /I "..\..\..\src\system_of_eqn\linearSOE\umfGEN" /I "..\..\..\src\system_of_eqn\linearSOE\fullGEN" /I "..\..\..\src\system_of_eqn\linearSOE\sparseGEN" /I "..\..\..\src\system_of_eqn\linearSOE\bandSPD" /I "..\..\..\src\system_of_eqn\linearSOE\bandGEN" /I "..\..\..\src\element\nonlinearBeamColumn\tcl\repres\section" /I "..\..\..\src\matrix" /I "..\..\..\src\recorder" /I "..\..\..\src\graph\numberer" /I "..\..\..\src\material\section" /I "..\..\..\src\graph\graph" /I "..\..\..\src\element\beam2d" /I "..\..\..\src\element\beam3d" /I "..\..\..\src\system_of_eqn" /I "..\..\..\src\system_of_eqn\linearSOE" /I "..\..\..\src\system_of_eqn\linearSOE\profileSPD" /I "..\..\..\src\system_of_eqn\linearSOE\sparseSYM" /I "..\..\..\src\domain\pattern" /I "..\..\..\src\analysis\analysis" /I "..\..\..\src\analysis\integrator" /I "..\..\..\src\analysis\numberer" /I "..\..\..\src\analysis\handler" /I "..\..\..\src\renderer" /I "..\..\..\src\material" /I "..\..\..\src\analysis\algorithm" /I "..\..\..\src\convergenceTest" /I "..\..\..\src\analysis\model\simple" /I "..\..\..\src\domain\load" /I "..\..\..\src\analysis\model" /I "..\..\..\src\element\truss" /I "..\..\..\src\actor\channel" /I "..\..\..\src\utility" /I "..\..\..\src\actor\actor" /I "..\..\..\src\modelbuilder" /I "..\..\..\src\modelbuilder\tcl" /I "..\..\..\src\domain\constraints" /I "..\..\..\src\domain\component" /I "..\..\..\src\element" /I "..\..\..\src\domain\node" /I "..\..\..\src\domain\domain" /I "..\..\..\src\tagged\storage" /I "..\..\..\src" /I "..\..\..\src\tagged" /I "..\..\..\src\material\nD\template3Dep" /I "c:\Program Files\tcl" /I "c:\Program Files\tcl\include" /D "_WIN32" /D "WIN32" /D "NDEBUG" /D "_CONSOLE" /D "_MBCS" /D "_TCL84" /FD /c
# SUBTRACT CPP /YX
# ADD BASE RSC /l 0x409 /d "NDEBUG"
# ADD RSC /l 0x409 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:console /machine:I386
# ADD LINK32 reliability.lib database.lib renderer.lib blas.lib lapack.lib feap.lib arpack.lib umfpack.lib openSeesFortran.lib actor.lib analysis.lib cblas.lib convergence.lib domain.lib element.lib graph.lib material.lib matrix.lib modelbuilder.lib recorder.lib superLU.lib system.lib tagged.lib utility.lib tcl84.lib tk84.lib kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib OpenGL32.lib glu32.lib GlAux.lib /nologo /subsystem:console /machine:I386 /nodefaultlib:"libc.lib" /libpath:"c:\Program Files\tcl\lib" /libpath:"c:\msdev\lib" /libpath:"..\..\lib" /libpath:"..\..\lib\release" /FORCE:MULTIPLE
# SUBTRACT LINK32 /pdb:none

!ELSEIF  "$(CFG)" == "quickMain - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "Debug"
# PROP BASE Intermediate_Dir "Debug"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "..\..\bin"
# PROP Intermediate_Dir "..\..\obj\quickMain\debug"
# PROP Ignore_Export_Lib 0
# PROP Target_Dir ""
# ADD BASE CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /GZ /c
# ADD CPP /nologo /W3 /GX /ZI /Od /I "..\..\src\domain\groundMotion" /I "..\..\src\domain\load" /I "..\..\src\domain\pattern" /I "..\..\src\element" /I "..\..\src\element\truss" /I "..\..\src\graph\graph" /I "..\..\src\graph\numberer" /I "..\..\src\system_of_eqn" /I "..\..\src\system_of_eqn\linearSOE" /I "..\..\src\system_of_eqn\linearSOE\sparseGen" /I "..\..\src\analysis\analysis" /I "..\..\src\convergenceTest" /I "..\..\src\analysis\integrator" /I "..\..\src\analysis\numberer" /I "..\..\src\domain\constraints" /I "..\..\src\analysis\handler" /I "..\..\src\analysis\algorithm" /I "..\..\src\analysis\model\simple" /I "..\..\src\matrix" /I "..\..\src" /I "..\..\src\actor\actor" /I "..\..\src\tagged" /I "..\..\src\domain\component" /I "..\..\src\analysis\algorithm\equiSolnAlgo" /I "..\..\src\analysis\model" /I "..\..\src\domain\domain" /I "..\..\src\domain\node" /I "..\..\src\modelbuilder" /I "..\..\src\utility" /I "..\..\..\src\ndarray" /I "..\..\..\src\domain\groundMotion" /I "..\..\..\src\element\nonlinearBeamColumn\quadrule" /I "..\..\..\src\material\section\fiber" /I "..\..\..\src\element\fourNodeQuad" /I "..\..\..\src\element\nonlinearBeamColumn\element" /I "..\..\..\src\element\zeroLength" /I "..\..\..\src\element\beamWithHinges" /I "..\..\..\src\recorder\response" /I "..\..\..\src\coordTransformation" /I "..\..\..\src\material\section\repres\section" /I "..\..\..\src\analysis\algorithm\equiSolnAlgo" /I "..\..\..\src\system_of_eqn\eigenSOE" /I "..\..\..\src\analysis\algorithm\eigenAlgo" /I "..\..\..\src\material\nD" /I "..\..\..\src\material\uniaxial" /I "..\..\..\src\handler" /I "..\..\..\src\tcl" /I "..\..\..\src\actor\objectBroker" /I "..\..\..\src\system_of_eqn\linearSOE\umfGEN" /I "..\..\..\src\system_of_eqn\linearSOE\fullGEN" /I "..\..\..\src\system_of_eqn\linearSOE\sparseGEN" /I "..\..\..\src\system_of_eqn\linearSOE\bandSPD" /I "..\..\..\src\system_of_eqn\linearSOE\bandGEN" /I "..\..\..\src\element\nonlinearBeamColumn\tcl\repres\section" /I "..\..\..\src\matrix" /I "..\..\..\src\recorder" /I "..\..\..\src\graph\numberer" /I "..\..\..\src\material\section" /I "..\..\..\src\graph\graph" /I "..\..\..\src\element\beam2d" /I "..\..\..\src\element\beam3d" /I "..\..\..\src\system_of_eqn" /I "..\..\..\src\system_of_eqn\linearSOE" /I "..\..\..\src\system_of_eqn\linearSOE\profileSPD" /I "..\..\..\src\system_of_eqn\linearSOE\sparseSYM" /I "..\..\..\src\domain\pattern" /I "..\..\..\src\analysis\analysis" /I "..\..\..\src\analysis\integrator" /I "..\..\..\src\analysis\numberer" /I "..\..\..\src\analysis\handler" /I "..\..\..\src\renderer" /I "..\..\..\src\material" /I "..\..\..\src\analysis\algorithm" /I "..\..\..\src\convergenceTest" /I "..\..\..\src\analysis\model\simple" /I "..\..\..\src\domain\load" /I "..\..\..\src\analysis\model" /I "..\..\..\src\element\truss" /I "..\..\..\src\actor\channel" /I "..\..\..\src\utility" /I "..\..\..\src\actor\actor" /I "..\..\..\src\modelbuilder" /I "..\..\..\src\modelbuilder\tcl" /I "..\..\..\src\domain\constraints" /I "..\..\..\src\domain\component" /I "..\..\..\src\element" /I "..\..\..\src\domain\node" /I "..\..\..\src\domain\domain" /I "..\..\..\src\tagged\storage" /I "..\..\..\src" /I "..\..\..\src\tagged" /I "..\..\..\src\material\nD\template3Dep" /I "c:\Program Files\tcl" /I "c:\Program Files\tcl\include" /D "_WIN32" /D "WIN32" /D "_DEBUG" /D "_CONSOLE" /D "_MBCS" /D "_TCL84" /FR /FD /GZ /c
# SUBTRACT CPP /YX
# ADD BASE RSC /l 0x409 /d "_DEBUG"
# ADD RSC /l 0x409 /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:console /debug /machine:I386 /pdbtype:sept
# ADD LINK32 reliability.lib actor.lib analysis.lib cblas.lib convergence.lib domain.lib element.lib graph.lib material.lib matrix.lib modelbuilder.lib recorder.lib superLU.lib system.lib tagged.lib utility.lib database.lib renderer.lib blas.lib lapack.lib feap.lib arpack.lib umfpack.lib openSeesFortran.lib tcl83.lib tk83.lib kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib OpenGL32.lib glu32.lib GlAux.lib /nologo /subsystem:console /profile /debug /machine:I386 /nodefaultlib:"libc.lib" /libpath:"c:\msdev\lib" /libpath:"c:\Program Files\tcl\lib" /libpath:"..\..\lib" /libpath:"..\..\lib\debug" /FORCE:MULTIPLE

!ENDIF 

# Begin Target

# Name "quickMain - Win32 Release"
# Name "quickMain - Win32 Debug"
# Begin Group "Source Files"

# PROP Default_Filter "cpp;c;cxx;rc;def;r;odl;idl;hpj;bat"
# Begin Source File

SOURCE=..\..\..\EXAMPLES\Example1\main.cpp
# End Source File
# End Group
# Begin Group "Header Files"

# PROP Default_Filter "h;hpp;hxx;hm;inl"
# End Group
# End Target
# End Project
