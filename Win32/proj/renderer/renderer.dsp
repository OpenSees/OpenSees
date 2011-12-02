# Microsoft Developer Studio Project File - Name="renderer" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Static Library" 0x0104

CFG=renderer - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "renderer.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "renderer.mak" CFG="renderer - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "renderer - Win32 Release" (based on "Win32 (x86) Static Library")
!MESSAGE "renderer - Win32 Debug" (based on "Win32 (x86) Static Library")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
RSC=rc.exe

!IF  "$(CFG)" == "renderer - Win32 Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "Release"
# PROP BASE Intermediate_Dir "Release"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "..\..\lib"
# PROP Intermediate_Dir "..\..\obj\renderer"
# PROP Target_Dir ""
# ADD BASE CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_MBCS" /D "_LIB" /YX /FD /c
# ADD CPP /nologo /MLd /W3 /GX /O2 /I "..\..\..\src\tcl" /I "..\..\..\src\handler" /I "..\..\..\src\domain\node" /I "..\..\..\src\actor\actor" /I "..\..\..\src\tagged" /I "..\..\..\src\domain\component" /I "..\..\..\src\renderer" /I "..\..\..\src\matrix" /I "..\..\..\src\utility" /I "..\..\..\src\domain\domain" /I "..\..\..\src" /I "..\..\..\src\element" /D "WIN32" /D "NDEBUG" /D "_MBCS" /D "_LIB" /FD /c
# SUBTRACT CPP /YX
# ADD BASE RSC /l 0x409 /d "NDEBUG"
# ADD RSC /l 0x409 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo

!ELSEIF  "$(CFG)" == "renderer - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "Debug"
# PROP BASE Intermediate_Dir "Debug"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "..\..\lib"
# PROP Intermediate_Dir "..\..\obj\render"
# PROP Target_Dir ""
# ADD BASE CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_MBCS" /D "_LIB" /YX /FD /GZ /c
# ADD CPP /nologo /W3 /Gm /GX /ZI /Od /I "..\..\..\src\tcl" /I "..\..\..\src\handler" /I "..\..\..\src\domain\node" /I "..\..\..\src\actor\actor" /I "..\..\..\src\tagged" /I "..\..\..\src\domain\component" /I "..\..\..\src\renderer" /I "..\..\..\src\matrix" /I "..\..\..\src\utility" /I "..\..\..\src\domain\domain" /I "..\..\..\src" /I "..\..\..\src\element" /D "_WIN32" /D "_DEBUG" /D "_MBCS" /D "_LIB" /FR /FD /GZ /c
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

# Name "renderer - Win32 Release"
# Name "renderer - Win32 Debug"
# Begin Group "Source Files"

# PROP Default_Filter "cpp;c;cxx;rc;def;r;odl;idl;hpj;bat"
# Begin Source File

SOURCE=..\..\..\SRC\renderer\Clipping.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\renderer\db.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\renderer\Device.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\renderer\DofColorMap.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\renderer\gMatrix.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\renderer\OpenGlRenderer.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\renderer\PlainMap.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\renderer\Projection.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\renderer\Renderer.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\renderer\Scan.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\renderer\View.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\renderer\Viewport.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\renderer\WindowDevice.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\renderer\WindowRenderer.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\renderer\X11Renderer.cpp
# End Source File
# End Group
# Begin Group "Header Files"

# PROP Default_Filter "h;hpp;hxx;hm;inl"
# Begin Source File

SOURCE=..\..\..\SRC\renderer\Clipping.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\renderer\ColorMap.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\renderer\container.H
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\renderer\db.H
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\renderer\Device.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\renderer\DofColorMap.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\renderer\gMatrix.H
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\renderer\OpenGLRenderer.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\renderer\PlainMap.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\renderer\Projection.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\renderer\Renderer.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\renderer\Scan.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\renderer\View.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\renderer\Viewport.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\renderer\WindowDevice.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\renderer\WindowRenderer.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\renderer\X11Renderer.h
# End Source File
# End Group
# End Target
# End Project
