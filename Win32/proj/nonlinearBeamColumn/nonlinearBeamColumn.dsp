# Microsoft Developer Studio Project File - Name="nonlinearBeamColumn" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Static Library" 0x0104

CFG=nonlinearBeamColumn - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "nonlinearBeamColumn.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "nonlinearBeamColumn.mak" CFG="nonlinearBeamColumn - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "nonlinearBeamColumn - Win32 Release" (based on "Win32 (x86) Static Library")
!MESSAGE "nonlinearBeamColumn - Win32 Debug" (based on "Win32 (x86) Static Library")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
RSC=rc.exe

!IF  "$(CFG)" == "nonlinearBeamColumn - Win32 Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "Release"
# PROP BASE Intermediate_Dir "Release"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "..\..\lib"
# PROP Intermediate_Dir "..\..\obj\nonlinearBeamColumn"
# PROP Target_Dir ""
# ADD BASE CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_MBCS" /D "_LIB" /YX /FD /c
# ADD CPP /nologo /MLd /W3 /GX /O2 /I "c:\Program Files\tcl\include" /I "..\..\..\src\coordTransformation" /I "..\..\..\src\handler" /I "..\..\..\src\material\uniaxial" /I "..\..\..\src\numberer" /I "..\..\..\src\material\section" /I "..\..\..\src\modelbuilder" /I "..\..\..\src\matrix" /I "..\..\..\src\renderer" /I "..\..\..\src\modelbuilder\tcl" /I "..\..\..\src\actor\objectBroker" /I "..\..\..\src\tagged\storage" /I "..\..\..\src\domain\node" /I "..\..\..\src\element\nonlinearBeamColumn\tcl\repres\section" /I "..\..\..\src\element\nonlinearBeamColumn\tcl\repres\patch" /I "..\..\..\src\element\nonlinearBeamColumn\tcl\repres\reinfBar" /I "..\..\..\src\element\nonlinearBeamColumn\tcl\repres\reinfLayer" /I "..\..\..\src\element\nonlinearBeamColumn\tcl\repres\cell" /I "..\..\..\src\domain\domain" /I "..\..\..\src\element\nonlinearBeamColumn" /I "..\..\..\src\element\nonlinearBEamCOlumn\quadrule" /I "..\..\..\src\element\nonlinearBeamColumn\matrixutil" /I "..\..\..\src" /I "..\..\..\src\element" /I "..\..\..\src\element\nonlinearBeamColumn\fiber" /I "..\..\..\src\domain\component" /I "..\..\..\src\tagged" /I "..\..\..\src\actor\actor" /I "..\..\..\src\actor\channel" /I "..\..\..\src\material" /I "..\..\..\src\element\nonlinearBeamColumn\element" /D "WIN32" /D "NDEBUG" /D "_MBCS" /D "_LIB" /FD /c
# SUBTRACT CPP /YX
# ADD BASE RSC /l 0x409 /d "NDEBUG"
# ADD RSC /l 0x409 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo

!ELSEIF  "$(CFG)" == "nonlinearBeamColumn - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "Debug"
# PROP BASE Intermediate_Dir "Debug"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "..\..\lib"
# PROP Intermediate_Dir "..\..\obj\nonlinearBeamColumn"
# PROP Target_Dir ""
# ADD BASE CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_MBCS" /D "_LIB" /YX /FD /GZ /c
# ADD CPP /nologo /W3 /GX /ZI /Od /I "c:\Program Files\tcl\include" /I "..\..\..\src\coordTransformation" /I "..\..\..\src\handler" /I "..\..\..\src\material\uniaxial" /I "..\..\..\src\numberer" /I "..\..\..\src\material\section" /I "..\..\..\src\modelbuilder" /I "..\..\..\src\matrix" /I "..\..\..\src\renderer" /I "..\..\..\src\modelbuilder\tcl" /I "..\..\..\src\actor\objectBroker" /I "..\..\..\src\tagged\storage" /I "..\..\..\src\domain\node" /I "..\..\..\src\element\nonlinearBeamColumn\tcl\repres\section" /I "..\..\..\src\element\nonlinearBeamColumn\tcl\repres\patch" /I "..\..\..\src\element\nonlinearBeamColumn\tcl\repres\reinfBar" /I "..\..\..\src\element\nonlinearBeamColumn\tcl\repres\reinfLayer" /I "..\..\..\src\element\nonlinearBeamColumn\tcl\repres\cell" /I "..\..\..\src\domain\domain" /I "..\..\..\src\element\nonlinearBeamColumn" /I "..\..\..\src\element\nonlinearBEamCOlumn\quadrule" /I "..\..\..\src\element\nonlinearBeamColumn\matrixutil" /I "..\..\..\src" /I "..\..\..\src\element" /I "..\..\..\src\element\nonlinearBeamColumn\fiber" /I "..\..\..\src\domain\component" /I "..\..\..\src\tagged" /I "..\..\..\src\actor\actor" /I "..\..\..\src\actor\channel" /I "..\..\..\src\material" /I "..\..\..\src\element\nonlinearBeamColumn\element" /D "WIN32" /D "_DEBUG" /D "_MBCS" /D "_LIB" /FR /FD /GZ /c
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

# Name "nonlinearBeamColumn - Win32 Release"
# Name "nonlinearBeamColumn - Win32 Debug"
# Begin Group "element"

# PROP Default_Filter ""
# Begin Source File

SOURCE=..\..\..\SRC\element\nonlinearBeamColumn\element\NLBeamColumn2d.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\element\nonlinearBeamColumn\element\NLBeamColumn2d.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\element\nonlinearBeamColumn\element\NLBeamColumn3d.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\element\nonlinearBeamColumn\element\NLBeamColumn3d.h
# End Source File
# End Group
# Begin Group "matrixutil"

# PROP Default_Filter ""
# Begin Source File

SOURCE=..\..\..\SRC\element\nonlinearBeamColumn\matrixutil\MatrixUtil.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\element\nonlinearBeamColumn\matrixutil\MatrixUtil.h
# End Source File
# End Group
# Begin Group "quadrule"

# PROP Default_Filter ""
# Begin Source File

SOURCE=..\..\..\SRC\element\nonlinearBeamColumn\quadrule\GaussLobattoQuadRule1d01.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\element\nonlinearBeamColumn\quadrule\GaussLobattoQuadRule1d01.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\element\nonlinearBeamColumn\quadrule\GaussQuadRule1d.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\element\nonlinearBeamColumn\quadrule\GaussQuadRule1d.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\element\nonlinearBeamColumn\quadrule\GaussQuadRule1d01.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\element\nonlinearBeamColumn\quadrule\GaussQuadRule1d01.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\element\nonlinearBeamColumn\quadrule\QuadRule.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\element\nonlinearBeamColumn\quadrule\QuadRule.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\element\nonlinearBeamColumn\quadrule\QuadRule1d.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\element\nonlinearBeamColumn\quadrule\QuadRule1d.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\element\nonlinearBeamColumn\quadrule\QuadRule1d01.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\element\nonlinearBeamColumn\quadrule\QuadRule1d01.h
# End Source File
# End Group
# Begin Source File

SOURCE=..\..\..\SRC\element\nonlinearBeamColumn\tcl\TclElmtBuilder.cpp
# End Source File
# End Target
# End Project
