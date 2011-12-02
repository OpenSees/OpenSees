# Microsoft Developer Studio Project File - Name="actor" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Static Library" 0x0104

CFG=actor - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "actor.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "actor.mak" CFG="actor - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "actor - Win32 Release" (based on "Win32 (x86) Static Library")
!MESSAGE "actor - Win32 Debug" (based on "Win32 (x86) Static Library")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
RSC=rc.exe

!IF  "$(CFG)" == "actor - Win32 Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "Release"
# PROP BASE Intermediate_Dir "Release"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "..\..\lib\release"
# PROP Intermediate_Dir "..\..\obj\actor\release"
# PROP Target_Dir ""
# ADD BASE CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_MBCS" /D "_LIB" /YX /FD /c
# ADD CPP /nologo /MLd /W3 /GX /O2 /I "..\..\..\src\element\elasticBeamColumn" /I "..\..\..\src\element\shell" /I "..\..\..\src\element\dispBeamColumn" /I "..\..\..\src\element\beamWithHinges" /I "..\..\..\src\material\section\fiber" /I "..\..\..\src\element\fourNodeQuad" /I "..\..\..\src\material\backbone" /I "..\..\..\src\material\state" /I "..\..\..\src\material\state\stiffness" /I "..\..\..\src\material\state\deformation" /I "..\..\..\src\material\state\strength" /I "..\..\..\src\material\nD" /I "..\..\..\src\coordTransformation" /I "..\..\..\src\domain\groundMotion" /I "..\..\..\src\system_of_eqn\linearSOE\profileSPD" /I "..\..\..\src\system_of_eqn\linearSOE" /I "..\..\..\src\system_of_eqn\linearSOE\sparseSYM" /I "..\..\..\src\analysis\integrator" /I "..\..\..\src\analysis\fe_ele" /I "..\..\..\src\analysis\dof_grp" /I "..\..\..\src\system_of_eqn\eigenSOE" /I "..\..\..\src\handler" /I "..\..\..\symSparse" /I "..\..\..\src\analysis\model\simple" /I "..\..\..\src\system_of_eqn\linearSOE\umfGEN" /I "..\..\..\src\system_of_eqn\linearSOE\fullGEN" /I "..\..\..\src\system_of_eqn\linearSOE\sparseGEN" /I "..\..\..\src\system_of_eqn\linearSOE\bandSPD" /I "..\..\..\src\system_of_eqn\linearSOE\bandGEN" /I "..\..\..\src\analysis\algorithm\domaindecompAlgo" /I "..\..\..\src\analysis\algorithm\eigenAlgo" /I "..\..\..\src\domain\domain\single" /I "..\..\..\src\convergenceTest" /I "..\..\..\src\analysis\analysis" /I "..\..\..\src\recorder" /I "..\..\..\src\analysis\algorithm\equiSolnAlgo" /I "..\..\..\src\analysis\algorithm" /I "..\..\..\src\system_of_eqn" /I "..\..\..\src\graph\graph" /I "..\..\..\src\graph\numberer" /I "..\..\..\src\analysis\numberer" /I "..\..\..\src\analysis\fe_ele\transformation" /I "..\..\..\src\analysis\fe_ele\lagrange" /I "..\..\..\src\analysis\fe_ele\penalty" /I "..\..\..\src\actor\objectBroker" /I "..\..\..\src\actor\channel" /I "..\..\..\src\utility" /I "..\..\..\src\domain\subdomain" /I "..\..\..\src\domain\constraints" /I "..\..\..\src\tagged" /I "..\..\..\src\domain\component" /I "..\..\..\src\domain\node" /I "..\..\..\src\element" /I "..\..\..\src\matrix" /I "..\..\..\src\domain\domain" /I "..\..\..\src\analysis\model" /I "..\..\..\src" /I "..\..\..\src\actor\actor" /I "..\..\..\src\analysis\handler" /I "..\..\..\src\material\section" /I "..\..\..\src\modelbuilder" /I "..\..\..\src\renderer" /I "..\..\..\src\modelbuilder\tcl" /I "..\..\..\src\tagged\storage" /I "..\..\..\src\element\nonlinearBeamColumn\tcl\repres\section" /I "..\..\..\src\element\nonlinearBeamColumn\tcl\repres\patch" /I "..\..\..\src\element\nonlinearBeamColumn\tcl\repres\reinfBar" /I "..\..\..\src\element\nonlinearBeamColumn\tcl\repres\reinfLayer" /I "..\..\..\src\element\nonlinearBeamColumn\tcl\repres\cell" /I "..\..\..\src\element\nonlinearBeamColumn" /I "..\..\..\src\element\nonlinearBEamCOlumn\quadrule" /I "..\..\..\src\element\nonlinearBeamColumn\matrixutil" /I "..\..\..\src\element\nonlinearBeamColumn\fiber" /I "..\..\..\src\material" /I "..\..\..\src\element\nonlinearBeamColumn\element" /I "..\..\..\src\domain\load" /I "..\..\..\src\domain\pattern" /I "..\..\..\src\element\zeroLength" /I "..\..\..\src\element\feap" /I "..\..\..\src\element\truss" /I "..\..\..\src\element\beam3d" /I "..\..\..\src\element\beam2d" /I "..\..\..\src\material\uniaxial" /I "..\..\..\src\actor\address" /I "..\..\..\src\actor\message" /I "..\..\..\src\nDarray" /I "..\..\..\src\material\uniaxial\fedeas" /I "..\..\..\src\material\uniaxial\drain" /I "..\..\..\src\element\8nbrick" /I "..\..\..\src\element\brick" /I "..\..\..\src\material\uniaxial\py" /I "..\..\..\src\material\nd\soil" /I "..\..\..\src\element\joint" /D "WIN32" /D "NDEBUG" /D "_MBCS" /D "_LIB" /D "_COROTATIONAL" /FD /c
# SUBTRACT CPP /YX
# ADD BASE RSC /l 0x409 /d "NDEBUG"
# ADD RSC /l 0x409 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo

!ELSEIF  "$(CFG)" == "actor - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "Debug"
# PROP BASE Intermediate_Dir "Debug"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "..\..\lib\debug"
# PROP Intermediate_Dir "..\..\obj\actor\debug"
# PROP Target_Dir ""
# ADD BASE CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_MBCS" /D "_LIB" /YX /FD /GZ /c
# ADD CPP /nologo /W3 /GX /ZI /Od /I "..\..\..\src\element\elasticBeamColumn" /I "..\..\..\src\element\shell" /I "..\..\..\src\element\dispBeamColumn" /I "..\..\..\src\element\beamWithHinges" /I "..\..\..\src\material\section\fiber" /I "..\..\..\src\element\fourNodeQuad" /I "..\..\..\src\material\backbone" /I "..\..\..\src\material\state" /I "..\..\..\src\material\state\stiffness" /I "..\..\..\src\material\state\deformation" /I "..\..\..\src\material\state\strength" /I "..\..\..\src\material\nD" /I "..\..\..\src\coordTransformation" /I "..\..\..\src\domain\groundMotion" /I "..\..\..\src\system_of_eqn\linearSOE\profileSPD" /I "..\..\..\src\system_of_eqn\linearSOE" /I "..\..\..\src\system_of_eqn\linearSOE\sparseSYM" /I "..\..\..\src\analysis\integrator" /I "..\..\..\src\analysis\fe_ele" /I "..\..\..\src\analysis\dof_grp" /I "..\..\..\src\system_of_eqn\eigenSOE" /I "..\..\..\src\handler" /I "..\..\..\symSparse" /I "..\..\..\src\analysis\model\simple" /I "..\..\..\src\system_of_eqn\linearSOE\umfGEN" /I "..\..\..\src\system_of_eqn\linearSOE\fullGEN" /I "..\..\..\src\system_of_eqn\linearSOE\sparseGEN" /I "..\..\..\src\system_of_eqn\linearSOE\bandSPD" /I "..\..\..\src\system_of_eqn\linearSOE\bandGEN" /I "..\..\..\src\analysis\algorithm\domaindecompAlgo" /I "..\..\..\src\analysis\algorithm\eigenAlgo" /I "..\..\..\src\domain\domain\single" /I "..\..\..\src\convergenceTest" /I "..\..\..\src\analysis\analysis" /I "..\..\..\src\recorder" /I "..\..\..\src\analysis\algorithm\equiSolnAlgo" /I "..\..\..\src\analysis\algorithm" /I "..\..\..\src\system_of_eqn" /I "..\..\..\src\graph\graph" /I "..\..\..\src\graph\numberer" /I "..\..\..\src\analysis\numberer" /I "..\..\..\src\analysis\fe_ele\transformation" /I "..\..\..\src\analysis\fe_ele\lagrange" /I "..\..\..\src\analysis\fe_ele\penalty" /I "..\..\..\src\actor\objectBroker" /I "..\..\..\src\actor\channel" /I "..\..\..\src\utility" /I "..\..\..\src\domain\subdomain" /I "..\..\..\src\domain\constraints" /I "..\..\..\src\tagged" /I "..\..\..\src\domain\component" /I "..\..\..\src\domain\node" /I "..\..\..\src\element" /I "..\..\..\src\matrix" /I "..\..\..\src\domain\domain" /I "..\..\..\src\analysis\model" /I "..\..\..\src" /I "..\..\..\src\actor\actor" /I "..\..\..\src\analysis\handler" /I "..\..\..\src\material\section" /I "..\..\..\src\modelbuilder" /I "..\..\..\src\renderer" /I "..\..\..\src\modelbuilder\tcl" /I "..\..\..\src\tagged\storage" /I "..\..\..\src\element\nonlinearBeamColumn\tcl\repres\section" /I "..\..\..\src\element\nonlinearBeamColumn\tcl\repres\patch" /I "..\..\..\src\element\nonlinearBeamColumn\tcl\repres\reinfBar" /I "..\..\..\src\element\nonlinearBeamColumn\tcl\repres\reinfLayer" /I "..\..\..\src\element\nonlinearBeamColumn\tcl\repres\cell" /I "..\..\..\src\element\nonlinearBeamColumn" /I "..\..\..\src\element\nonlinearBEamCOlumn\quadrule" /I "..\..\..\src\element\nonlinearBeamColumn\matrixutil" /I "..\..\..\src\element\nonlinearBeamColumn\fiber" /I "..\..\..\src\material" /I "..\..\..\src\element\nonlinearBeamColumn\element" /I "..\..\..\src\domain\load" /I "..\..\..\src\domain\pattern" /I "..\..\..\src\element\zeroLength" /I "..\..\..\src\element\feap" /I "..\..\..\src\element\truss" /I "..\..\..\src\element\beam3d" /I "..\..\..\src\element\beam2d" /I "..\..\..\src\material\uniaxial" /I "..\..\..\src\actor\address" /I "..\..\..\src\actor\message" /I "..\..\..\src\nDarray" /I "..\..\..\src\material\uniaxial\fedeas" /I "..\..\..\src\material\uniaxial\drain" /I "..\..\..\src\element\8nbrick" /I "..\..\..\src\element\brick" /I "..\..\..\src\material\uniaxial\py" /I "..\..\..\src\material\nd\soil" /I "..\..\..\src\element\joint" /D "_DEBUG" /D "_MBCS" /D "_LIB" /D "_COROTATIONAL" /FR /FD /c
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

# Name "actor - Win32 Release"
# Name "actor - Win32 Debug"
# Begin Source File

SOURCE=..\..\..\SRC\actor\channel\Channel.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\actor\channel\Channel.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\actor\address\ChannelAddress.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\actor\address\ChannelAddress.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\actor\objectBroker\FEM_ObjectBroker.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\actor\objectBroker\FEM_ObjectBroker.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\actor\message\Message.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\actor\message\Message.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\actor\actor\MovableObject.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\actor\actor\MovableObject.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\actor\objectBroker\ObjectBroker.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\actor\objectBroker\ObjectBroker.h
# End Source File
# End Target
# End Project
