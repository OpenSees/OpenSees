# Microsoft Developer Studio Project File - Name="recorder" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Static Library" 0x0104

CFG=recorder - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "recorder.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "recorder.mak" CFG="recorder - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "recorder - Win32 Release" (based on "Win32 (x86) Static Library")
!MESSAGE "recorder - Win32 Debug" (based on "Win32 (x86) Static Library")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
RSC=rc.exe

!IF  "$(CFG)" == "recorder - Win32 Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "Release"
# PROP BASE Intermediate_Dir "Release"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "..\..\lib\release"
# PROP Intermediate_Dir "..\..\obj\recorder\release"
# PROP Target_Dir ""
# ADD BASE CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_MBCS" /D "_LIB" /YX /FD /c
# ADD CPP /nologo /MLd /W3 /GX /O2 /I "..\..\..\src\modelbuilder" /I "..\..\..\src\modelbuilder\tcl" /I "..\..\..\src\damage" /I "..\..\..\src\domain\pattern" /I "..\..\..\src\system_of_eqn" /I "..\..\..\src\system_of_eqn\linearSOE" /I "..\..\..\src\analysis\integrator" /I "..\..\..\src\material\section\fiber" /I "..\..\..\src\material" /I "..\..\..\src\recorder\response" /I "..\..\..\src\analysis\algorithm" /I "..\..\..\src\analysis\algorithm\equiSolnAlgo" /I "..\..\..\src\handler" /I "..\..\..\src\renderer" /I "..\..\..\src\tcl" /I "..\..\..\src\actor\channel" /I "..\..\..\src\database" /I "..\..\..\src\tagged" /I "..\..\..\src\actor\actor" /I "..\..\..\src\domain\component" /I "..\..\..\src\element" /I "..\..\..\src\matrix" /I "..\..\..\src\domain\domain" /I "..\..\..\src" /I "..\..\..\src\recorder" /I "..\..\..\src\domain\node" /I "..\..\..\src\nDarray" /I "..\..\..\src\domain\region" /I "..\..\..\src\domain\constraints" /I "..\..\..\src\element\updatedlagrangianbeamcolumn" /I "..\..\..\src\material\yieldSurface\evolution" /I "..\..\..\src\material\yieldsurface\yieldSurfaceBC" /I "c:\Program Files\tcl" /D "NDEBUG" /D "WIN32" /D "_LIB" /D "_WGL" /D "_MBCS" /D "_TCL84" /FD /c
# SUBTRACT CPP /YX
# ADD BASE RSC /l 0x409 /d "NDEBUG"
# ADD RSC /l 0x409 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo

!ELSEIF  "$(CFG)" == "recorder - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "Debug"
# PROP BASE Intermediate_Dir "Debug"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "..\..\lib\debug"
# PROP Intermediate_Dir "..\..\obj\recorder\debug"
# PROP Target_Dir ""
# ADD BASE CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_MBCS" /D "_LIB" /YX /FD /GZ /c
# ADD CPP /nologo /W3 /GX /ZI /Od /I "..\..\..\src\modelbuilder" /I "..\..\..\src\modelbuilder\tcl" /I "..\..\..\src\damage" /I "..\..\..\src\domain\pattern" /I "..\..\..\src\system_of_eqn" /I "..\..\..\src\system_of_eqn\linearSOE" /I "..\..\..\src\analysis\integrator" /I "..\..\..\src\material\section\fiber" /I "..\..\..\src\material" /I "..\..\..\src\recorder\response" /I "..\..\..\src\analysis\algorithm" /I "..\..\..\src\analysis\algorithm\equiSolnAlgo" /I "..\..\..\src\handler" /I "..\..\..\src\renderer" /I "..\..\..\src\tcl" /I "..\..\..\src\actor\channel" /I "..\..\..\src\database" /I "..\..\..\src\tagged" /I "..\..\..\src\actor\actor" /I "..\..\..\src\domain\component" /I "..\..\..\src\element" /I "..\..\..\src\matrix" /I "..\..\..\src\domain\domain" /I "..\..\..\src" /I "..\..\..\src\recorder" /I "..\..\..\src\domain\node" /I "..\..\..\src\nDarray" /I "..\..\..\src\domain\region" /I "..\..\..\src\domain\constraints" /I "..\..\..\src\element\updatedlagrangianbeamcolumn" /I "..\..\..\src\material\yieldSurface\evolution" /I "..\..\..\src\material\yieldsurface\yieldSurfaceBC" /I "c:\Program Files\tcl" /D "_DEBUG" /D "WIN32" /D "_LIB" /D "_WGL" /D "_MBCS" /D "_TCL84" /FR /FD /GZ /c
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

# Name "recorder - Win32 Release"
# Name "recorder - Win32 Debug"
# Begin Group "Source Files"

# PROP Default_Filter "cpp;c;cxx;rc;def;r;odl;idl;hpj;bat"
# Begin Source File

SOURCE=..\..\..\SRC\recorder\AlgorithmIncrements.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\recorder\DamageRecorder.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\recorder\DatastoreRecorder.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\recorder\DriftRecorder.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\recorder\ElementRecorder.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\recorder\EnvelopeElementRecorder.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\recorder\EnvelopeNodeRecorder.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\recorder\FilePlotter.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\Src\recorder\GSA_Recorder.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\recorder\MaxNodeDispRecorder.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\recorder\NodeRecorder.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\recorder\PatternRecorder.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\recorder\TclRecorderCommands.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\recorder\YsVisual.cpp
# End Source File
# End Group
# Begin Group "Header Files"

# PROP Default_Filter "h;hpp;hxx;hm;inl"
# Begin Source File

SOURCE=..\..\..\SRC\recorder\AlgorithmIncrements.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\recorder\DamageRecorder.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\recorder\DatastoreRecorder.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\recorder\DriftRecorder.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\recorder\ElementRecorder.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\recorder\EnvelopeElementRecorder.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\recorder\EnvelopeNodeRecorder.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\recorder\FilePlotter.h
# End Source File
# Begin Source File

SOURCE=..\..\..\Src\recorder\GSA_Recorder.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\recorder\MaxNodeDispRecorder.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\recorder\NodeRecorder.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\recorder\PatternRecorder.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\recorder\Recorder.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\recorder\YsVisual.h
# End Source File
# End Group
# Begin Group "response"

# PROP Default_Filter ""
# Begin Source File

SOURCE=..\..\..\SRC\recorder\response\ElementResponse.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\recorder\response\ElementResponse.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\recorder\response\FiberResponse.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\recorder\response\FiberResponse.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\recorder\response\MaterialResponse.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\recorder\response\MaterialResponse.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\recorder\response\Response.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\recorder\response\Response.h
# End Source File
# End Group
# End Target
# End Project
