# Microsoft Developer Studio Project File - Name="matrix" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Static Library" 0x0104

CFG=matrix - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "matrix.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "matrix.mak" CFG="matrix - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "matrix - Win32 Release" (based on "Win32 (x86) Static Library")
!MESSAGE "matrix - Win32 Debug" (based on "Win32 (x86) Static Library")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
F90=df.exe
RSC=rc.exe

!IF  "$(CFG)" == "matrix - Win32 Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "Release"
# PROP BASE Intermediate_Dir "Release"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "..\..\lib\release"
# PROP Intermediate_Dir "..\..\obj\matrix\release"
# PROP Target_Dir ""
# ADD BASE CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_MBCS" /D "_LIB" /YX /FD /c
# ADD CPP /nologo /MLd /W3 /GX /O2 /I "..\..\..\src\matrix" /I "..\..\..\src\handler" /I "..\..\..\src" /I "..\..\..\src\nDarray" /I "c:\Program Files\tcl" /I "c:\Program Files\tcl\include" /D "WIN32" /D "NDEBUG" /D "_LIB" /D "_MBCS" /D "_TCL84" /FD /c
# SUBTRACT CPP /YX
# ADD BASE RSC /l 0x409 /d "NDEBUG"
# ADD RSC /l 0x409 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo

!ELSEIF  "$(CFG)" == "matrix - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "Debug"
# PROP BASE Intermediate_Dir "Debug"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "..\..\lib\debug"
# PROP Intermediate_Dir "..\..\obj\matrix\debug"
# PROP Target_Dir ""
# ADD BASE CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_MBCS" /D "_LIB" /YX /FD /GZ /c
# ADD CPP /nologo /W3 /GX /ZI /Od /I "..\..\..\src\matrix" /I "..\..\..\src\handler" /I "..\..\..\src" /I "..\..\..\src\nDarray" /I "c:\Program Files\tcl" /I "c:\Program Files\tcl\include" /D "WIN32" /D "_DEBUG" /D "_LIB" /D "_MBCS" /D "_TCL84" /FR /FD /GZ /c
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

# Name "matrix - Win32 Release"
# Name "matrix - Win32 Debug"
# Begin Group "Source Files"

# PROP Default_Filter "cpp;c;cxx;rc;def;r;odl;idl;hpj;bat"
# Begin Source File

SOURCE=..\..\..\Src\nDarray\basics.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\Src\nDarray\BJmatrix.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\Src\nDarray\BJtensor.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\Src\nDarray\BJvector.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\matrix\ID.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\matrix\Matrix.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\Src\nDarray\nDarray.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\Src\nDarray\straint.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\Src\nDarray\stresst.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\matrix\Vector.cpp
# End Source File
# End Group
# Begin Group "Header Files"

# PROP Default_Filter "h;hpp;hxx;hm;inl"
# Begin Source File

SOURCE=..\..\..\Src\nDarray\basics.h
# End Source File
# Begin Source File

SOURCE=..\..\..\Src\nDarray\BJmatrix.h
# End Source File
# Begin Source File

SOURCE=..\..\..\Src\nDarray\BJtensor.h
# End Source File
# Begin Source File

SOURCE=..\..\..\Src\nDarray\BJvector.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\matrix\ID.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\matrix\Matrix.h
# End Source File
# Begin Source File

SOURCE=..\..\..\Src\nDarray\nDarray.h
# End Source File
# Begin Source File

SOURCE=..\..\..\Src\nDarray\straint.h
# End Source File
# Begin Source File

SOURCE=..\..\..\Src\nDarray\stresst.h
# End Source File
# Begin Source File

SOURCE=..\..\..\Src\nDarray\Tensor.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\matrix\Vector.h
# End Source File
# End Group
# End Target
# End Project
