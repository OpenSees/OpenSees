# Microsoft Developer Studio Project File - Name="graph" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Static Library" 0x0104

CFG=graph - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "graph.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "graph.mak" CFG="graph - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "graph - Win32 Release" (based on "Win32 (x86) Static Library")
!MESSAGE "graph - Win32 Debug" (based on "Win32 (x86) Static Library")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
RSC=rc.exe

!IF  "$(CFG)" == "graph - Win32 Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "Release"
# PROP BASE Intermediate_Dir "Release"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "..\..\lib\release"
# PROP Intermediate_Dir "..\..\obj\graph\release"
# PROP Target_Dir ""
F90=df.exe
# ADD BASE CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_MBCS" /D "_LIB" /YX /FD /c
# ADD CPP /nologo /MLd /W3 /GX /O2 /I "..\..\..\src\handler" /I "..\..\..\src\actor\channel" /I "..\..\..\src\graph\numberer" /I "..\..\..\src\actor\objectBroker" /I "..\..\..\src\analysis\dof_grp" /I "..\..\..\src\analysis\fe_ele" /I "..\..\..\src\actor\actor" /I "..\..\..\src\analysis\model\simple" /I "..\..\..\src\analysis\model" /I "..\..\..\src\tagged\storage" /I "..\..\..\src\matrix" /I "..\..\..\src\graph\graph" /I "..\..\..\src" /I "..\..\..\src\tagged" /I "..\..\..\src\nDarray" /I "c:\Program Files\tcl" /I "c:\Program Files\tcl\include" /D "WIN32" /D "NDEBUG" /D "_LIB" /D "_MBCS" /D "_TCL84" /D "_VC6" /FD /c
# SUBTRACT CPP /YX
# ADD BASE RSC /l 0x409 /d "NDEBUG"
# ADD RSC /l 0x409 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo

!ELSEIF  "$(CFG)" == "graph - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "Debug"
# PROP BASE Intermediate_Dir "Debug"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "..\..\lib\debug"
# PROP Intermediate_Dir "..\..\obj\graph\debug"
# PROP Target_Dir ""
F90=df.exe
# ADD BASE CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_MBCS" /D "_LIB" /YX /FD /GZ /c
# ADD CPP /nologo /W3 /GX /ZI /Od /I "..\..\..\src\handler" /I "..\..\..\src\actor\channel" /I "..\..\..\src\graph\numberer" /I "..\..\..\src\actor\objectBroker" /I "..\..\..\src\analysis\dof_grp" /I "..\..\..\src\analysis\fe_ele" /I "..\..\..\src\actor\actor" /I "..\..\..\src\analysis\model\simple" /I "..\..\..\src\analysis\model" /I "..\..\..\src\tagged\storage" /I "..\..\..\src\matrix" /I "..\..\..\src\graph\graph" /I "..\..\..\src" /I "..\..\..\src\tagged" /I "..\..\..\src\nDarray" /I "c:\Program Files\tcl" /I "c:\Program Files\tcl\include" /D "WIN32" /D "_DEBUG" /D "_LIB" /D "_MBCS" /D "_TCL84" /D "_VC6" /FR /FD /GZ /c
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

# Name "graph - Win32 Release"
# Name "graph - Win32 Debug"
# Begin Group "Source Files"

# PROP Default_Filter "cpp;c;cxx;rc;def;r;odl;idl;hpj;bat"
# Begin Source File

SOURCE=..\..\..\SRC\graph\graph\ArrayGraph.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\graph\graph\ArrayVertexIter.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\graph\graph\DOF_Graph.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\graph\graph\DOF_GroupGraph.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\graph\graph\Graph.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\graph\numberer\GraphNumberer.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\graph\numberer\MyRCM.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\graph\numberer\RCM.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\graph\numberer\SimpleNumberer.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\graph\graph\Vertex.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\graph\graph\VertexIter.cpp
# End Source File
# End Group
# Begin Group "Header Files"

# PROP Default_Filter "h;hpp;hxx;hm;inl"
# Begin Source File

SOURCE=..\..\..\SRC\graph\graph\ArrayGraph.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\graph\graph\ArrayVertexIter.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\graph\graph\DOF_Graph.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\graph\graph\DOF_GroupGraph.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\graph\graph\Graph.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\graph\numberer\GraphNumberer.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\graph\numberer\MyRCM.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\graph\numberer\RCM.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\graph\numberer\SimpleNumberer.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\graph\graph\Vertex.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\graph\graph\VertexIter.h
# End Source File
# End Group
# End Target
# End Project
