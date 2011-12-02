# Microsoft Developer Studio Project File - Name="handler" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Static Library" 0x0104

CFG=handler - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "handler.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "handler.mak" CFG="handler - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "handler - Win32 Release" (based on "Win32 (x86) Static Library")
!MESSAGE "handler - Win32 Debug" (based on "Win32 (x86) Static Library")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
F90=df.exe
RSC=rc.exe

!IF  "$(CFG)" == "handler - Win32 Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "Release"
# PROP BASE Intermediate_Dir "Release"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "..\..\lib\release"
# PROP Intermediate_Dir "..\..\obj\handler\release"
# PROP Target_Dir ""
# ADD BASE CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_MBCS" /D "_LIB" /YX /FD /c
# ADD CPP /nologo /MLd /W3 /GX /O2 /I "..\..\..\src\actor\channel" /I "..\..\..\src\objectBroker" /I "..\..\..\src\actor\message" /I "..\..\..\src\matrix" /I "..\..\..\src\database" /I "..\..\..\src\ndarray" /I "..\..\..\src" /I "..\..\..\src\actor\actor" /I "..\..\..\src\tagged" /I "..\..\..\src\handler" /I "c:\Program Files\tcl" /I "c:\Program Files\tcl\include" /D "WIN32" /D "NDEBUG" /D "_LIB" /D "_MBCS" /D "_TCL84" /FD /c
# SUBTRACT CPP /YX
# ADD BASE RSC /l 0x409 /d "NDEBUG"
# ADD RSC /l 0x409 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo

!ELSEIF  "$(CFG)" == "handler - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "Debug"
# PROP BASE Intermediate_Dir "Debug"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "..\..\lib\debug"
# PROP Intermediate_Dir "..\..\obj\handler\debug"
# PROP Target_Dir ""
# ADD BASE CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_MBCS" /D "_LIB" /YX /FD /GZ /c
# ADD CPP /nologo /W3 /Gm /GX /ZI /Od /I "..\..\..\src\actor\channel" /I "..\..\..\src\objectBroker" /I "..\..\..\src\actor\message" /I "..\..\..\src\matrix" /I "..\..\..\src\database" /I "..\..\..\src\ndarray" /I "..\..\..\src" /I "..\..\..\src\actor\actor" /I "..\..\..\src\tagged" /I "..\..\..\src\handler" /I "c:\Program Files\tcl" /I "c:\Program Files\tcl\include" /D "WIN32" /D "_DEBUG" /D "_LIB" /D "_MBCS" /D "_TCL84" /FR /FD /GZ /c
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

# Name "handler - Win32 Release"
# Name "handler - Win32 Debug"
# Begin Group "Source Files"

# PROP Default_Filter "cpp;c;cxx;rc;def;r;odl;idl;hpj;bat"
# Begin Source File

SOURCE=..\..\..\SRC\handler\DataOutputDatabaseHandler.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\handler\DataOutputFileHandler.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\handler\DataOutputHandler.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\handler\DataOutputStreamHandler.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\Src\handler\FileStream.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\Src\handler\OPS_Stream.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\Src\handler\StandardStream.cpp
# End Source File
# End Group
# Begin Group "Header Files"

# PROP Default_Filter "h;hpp;hxx;hm;inl"
# Begin Source File

SOURCE=..\..\..\Src\handler\FileStream.h
# End Source File
# Begin Source File

SOURCE=..\..\..\Src\handler\OPS_Stream.h
# End Source File
# Begin Source File

SOURCE=..\..\..\Src\handler\StandardStream.h
# End Source File
# End Group
# End Target
# End Project
