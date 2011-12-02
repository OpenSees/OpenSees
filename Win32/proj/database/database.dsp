# Microsoft Developer Studio Project File - Name="database" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Static Library" 0x0104

CFG=database - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "database.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "database.mak" CFG="database - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "database - Win32 Release" (based on "Win32 (x86) Static Library")
!MESSAGE "database - Win32 Debug" (based on "Win32 (x86) Static Library")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
RSC=rc.exe

!IF  "$(CFG)" == "database - Win32 Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "Release"
# PROP BASE Intermediate_Dir "Release"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "..\..\lib\release"
# PROP Intermediate_Dir "..\..\obj\database\release"
# PROP Target_Dir ""
# ADD BASE CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_MBCS" /D "_LIB" /YX /FD /c
# ADD CPP /nologo /MLd /W3 /GX /O2 /I "..\..\..\src\element\fourNodeQuad" /I "equiSolnAlgo" /I "..\..\..\src\analysis\integrator" /I "..\..\..\src\analysis\algorithm\equiSolnAlgo" /I "..\..\..\src\analysis\algorithm" /I "..\..\..\src\handler" /I "..\..\..\src\database" /I "..\..\..\src\actor\channel" /I ".\..\..\src\domain\groundMotion" /I "..\..\..\src\analysis\fe_ele" /I "..\..\..\src\utility" /I "..\..\..\src\domain\subdomain" /I "..\..\..\src\analysis\analysis" /I "..\..\..\src\recorder" /I "..\..\..\src\graph\graph" /I "..\..\..\src\modelbuilder" /I "..\..\..\src\domain\domain\single" /I "..\..\..\src\modelbuilder\tcl" /I "..\..\..\src\tagged\storage" /I "..\..\..\src\domain\constraints" /I "..\..\..\src\domain\pattern" /I "..\..\..\src\matrix" /I "..\..\..\src\analysis\dof_grp" /I "..\..\..\src\domain\domain" /I "..\..\..\src" /I "..\..\..\src\actor\actor" /I "..\..\..\src\actor\objectBroker" /I "..\..\..\src\tagged" /I "..\..\..\src\domain\component" /I "..\..\..\src\domain\load" /I "..\..\..\src\domain\node" /I "..\..\..\src\element" /I "..\..\..\src\renderer" /I "..\..\..\src\nDarray" /I "c:\Program Files\tcl" /D "WIN32" /D "NDEBUG" /D "_LIB" /D "_MBCS" /D "_TCL84" /FD /c
# SUBTRACT CPP /YX
# ADD BASE RSC /l 0x409 /d "NDEBUG"
# ADD RSC /l 0x409 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo

!ELSEIF  "$(CFG)" == "database - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "Debug"
# PROP BASE Intermediate_Dir "Debug"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "..\..\lib\debug"
# PROP Intermediate_Dir "..\..\obj\database\debug"
# PROP Target_Dir ""
# ADD BASE CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_MBCS" /D "_LIB" /YX /FD /GZ /c
# ADD CPP /nologo /W3 /Gm /GX /ZI /Od /I "..\..\..\src\analysis\integrator" /I "..\..\..\src\analysis\algorithm\equiSolnAlgo" /I "..\..\..\src\analysis\algorithm" /I "..\..\..\src\handler" /I "..\..\..\src\database" /I "..\..\..\src\actor\channel" /I ".\..\..\src\domain\groundMotion" /I "..\..\..\src\analysis\fe_ele" /I "..\..\..\src\utility" /I "..\..\..\src\domain\subdomain" /I "..\..\..\src\analysis\analysis" /I "..\..\..\src\recorder" /I "..\..\..\src\graph\graph" /I "..\..\..\src\modelbuilder" /I "..\..\..\src\domain\domain\single" /I "..\..\..\src\modelbuilder\tcl" /I "..\..\..\src\tagged\storage" /I "..\..\..\src\domain\constraints" /I "..\..\..\src\domain\pattern" /I "..\..\..\src\matrix" /I "..\..\..\src\analysis\dof_grp" /I "..\..\..\src\domain\domain" /I "..\..\..\src" /I "..\..\..\src\actor\actor" /I "..\..\..\src\actor\objectBroker" /I "..\..\..\src\tagged" /I "..\..\..\src\domain\component" /I "..\..\..\src\domain\load" /I "..\..\..\src\domain\node" /I "..\..\..\src\element" /I "..\..\..\src\renderer" /I "..\..\..\src\nDarray" /I "c:\Program Files\tcl" /D "WIN32" /D "_DEBUG" /D "_LIB" /D "_MBCS" /D "_TCL84" /FR /FD /GZ /c
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

# Name "database - Win32 Release"
# Name "database - Win32 Debug"
# Begin Group "Source Files"

# PROP Default_Filter "cpp;c;cxx;rc;def;r;odl;idl;hpj;bat"
# Begin Source File

SOURCE=..\..\..\SRC\database\FE_Datastore.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\database\FileDatastore.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\database\NEESData.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\database\TclDatabaseCommands.cpp
# End Source File
# End Group
# Begin Group "Header Files"

# PROP Default_Filter "h;hpp;hxx;hm;inl"
# Begin Source File

SOURCE=..\..\..\SRC\database\FE_Datastore.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\database\FileDatastore.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\database\NEESData.h
# End Source File
# End Group
# End Target
# End Project
