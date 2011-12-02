# Microsoft Developer Studio Project File - Name="domain" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Static Library" 0x0104

CFG=domain - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "domain.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "domain.mak" CFG="domain - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "domain - Win32 Release" (based on "Win32 (x86) Static Library")
!MESSAGE "domain - Win32 Debug" (based on "Win32 (x86) Static Library")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
RSC=rc.exe

!IF  "$(CFG)" == "domain - Win32 Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "Release"
# PROP BASE Intermediate_Dir "Release"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "..\..\lib"
# PROP Intermediate_Dir "..\..\obj\domain"
# PROP Target_Dir ""
# ADD BASE CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_MBCS" /D "_LIB" /YX /FD /c
# ADD CPP /nologo /MLd /W3 /GX /O2 /I "c:\Program Files\tcl\include" /I "..\..\..\src\handler" /I "..\..\..\src\domain\groundMotion" /I "..\..\..\src\analysis\fe_ele" /I "..\..\..\src\utility" /I "..\..\..\src\domain\subdomain" /I "..\..\..\src\database" /I "..\..\..\src\analysis\analysis" /I "..\..\..\src\recorder" /I "..\..\..\src\graph\graph" /I "..\..\..\src\modelbuilder" /I "..\..\..\src\domain\domain\single" /I "..\..\..\src\modelbuilder\tcl" /I "..\..\..\src\tagged\storage" /I "..\..\..\src\domain\constraints" /I "..\..\..\src\domain\pattern" /I "..\..\..\src\matrix" /I "..\..\..\src\analysis\dof_grp" /I "..\..\..\src\domain\domain" /I "..\..\..\src\actor\channel" /I "..\..\..\src" /I "..\..\..\src\actor\actor" /I "..\..\..\src\actor\objectBroker" /I "..\..\..\src\tagged" /I "..\..\..\src\domain\component" /I "..\..\..\src\domain\load" /I "..\..\..\src\domain\node" /I "..\..\..\src\element" /I "..\..\..\src\renderer" /D "WIN32" /D "NDEBUG" /D "_MBCS" /D "_LIB" /FD /c
# SUBTRACT CPP /Fr /YX
# ADD BASE RSC /l 0x409 /d "NDEBUG"
# ADD RSC /l 0x409 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo

!ELSEIF  "$(CFG)" == "domain - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "Debug"
# PROP BASE Intermediate_Dir "Debug"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "..\..\lib"
# PROP Intermediate_Dir "..\..\obj\domain"
# PROP Target_Dir ""
# ADD BASE CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_MBCS" /D "_LIB" /YX /FD /GZ /c
# ADD CPP /nologo /W3 /GX /ZI /Od /I "c:\Program Files\tcl\include" /I "..\..\..\src\handler" /I "..\..\..\src\domain\groundMotion" /I "..\..\..\src\analysis\fe_ele" /I "..\..\..\src\utility" /I "..\..\..\src\domain\subdomain" /I "..\..\..\src\database" /I "..\..\..\src\analysis\analysis" /I "..\..\..\src\recorder" /I "..\..\..\src\graph\graph" /I "..\..\..\src\modelbuilder" /I "..\..\..\src\domain\domain\single" /I "..\..\..\src\modelbuilder\tcl" /I "..\..\..\src\tagged\storage" /I "..\..\..\src\domain\constraints" /I "..\..\..\src\domain\pattern" /I "..\..\..\src\matrix" /I "..\..\..\src\analysis\dof_grp" /I "..\..\..\src\domain\domain" /I "..\..\..\src\actor\channel" /I "..\..\..\src" /I "..\..\..\src\actor\actor" /I "..\..\..\src\actor\objectBroker" /I "..\..\..\src\tagged" /I "..\..\..\src\domain\component" /I "..\..\..\src\domain\load" /I "..\..\..\src\domain\node" /I "..\..\..\src\element" /I "..\..\..\src\renderer" /D "WIN32" /D "_DEBUG" /D "_MBCS" /D "_LIB" /FR /FD /GZ /c
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

# Name "domain - Win32 Release"
# Name "domain - Win32 Debug"
# Begin Group "node"

# PROP Default_Filter ""
# Begin Source File

SOURCE=..\..\..\SRC\domain\node\NodalLoad.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\domain\node\NodalLoad.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\domain\node\Node.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\domain\node\Node.h
# End Source File
# End Group
# Begin Group "domain"

# PROP Default_Filter ""
# Begin Source File

SOURCE=..\..\..\SRC\domain\domain\Domain.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\domain\domain\Domain.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\domain\domain\ElementIter.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\domain\domain\MP_ConstraintIter.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\domain\domain\NodeIter.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\domain\domain\single\SingleDomAllSP_Iter.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\domain\domain\single\SingleDomAllSP_Iter.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\domain\domain\single\SingleDomEleIter.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\domain\domain\single\SingleDomEleIter.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\domain\domain\single\SingleDomMP_Iter.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\domain\domain\single\SingleDomMP_Iter.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\domain\domain\single\SingleDomNodIter.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\domain\domain\single\SingleDomNodIter.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\domain\domain\single\SingleDomSP_Iter.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\domain\domain\single\SingleDomSP_Iter.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\domain\domain\SP_ConstraintIter.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\domain\subdomain\Subdomain.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\domain\subdomain\Subdomain.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\domain\domain\SubdomainIter.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\domain\subdomain\SubdomainNodIter.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\domain\subdomain\SubdomainNodIter.h
# End Source File
# End Group
# Begin Group "component"

# PROP Default_Filter ""
# Begin Source File

SOURCE=..\..\..\SRC\domain\component\DomainComponent.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\domain\component\DomainComponent.h
# End Source File
# End Group
# Begin Group "load"

# PROP Default_Filter ""
# Begin Source File

SOURCE=..\..\..\SRC\domain\load\ElementalLoadIter.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\domain\load\ElementalLoadIter.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\domain\load\Load.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\domain\load\Load.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\domain\load\NodalLoadIter.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\domain\load\NodalLoadIter.h
# End Source File
# End Group
# Begin Group "constraints"

# PROP Default_Filter ""
# Begin Source File

SOURCE=..\..\..\SRC\domain\constraints\MP_Constraint.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\domain\constraints\MP_Constraint.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\domain\constraints\RigidBeam.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\domain\constraints\RigidBeam.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\domain\constraints\RigidDiaphragm.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\domain\constraints\RigidDiaphragm.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\domain\constraints\RigidRod.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\domain\constraints\RigidRod.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\domain\constraints\SP_Constraint.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\domain\constraints\SP_Constraint.h
# End Source File
# End Group
# Begin Group "pattern"

# PROP Default_Filter ""
# Begin Source File

SOURCE=..\..\..\SRC\domain\pattern\ConstantSeries.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\domain\pattern\ConstantSeries.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\domain\pattern\EarthquakePattern.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\domain\pattern\EarthquakePattern.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\domain\pattern\LinearSeries.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\domain\pattern\LinearSeries.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\domain\pattern\LoadPattern.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\domain\pattern\LoadPattern.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\domain\pattern\LoadPatternIter.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\domain\pattern\LoadPatternIter.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\domain\pattern\PathSeries.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\domain\pattern\PathSeries.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\domain\pattern\PathTimeSeries.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\domain\pattern\PathTimeSeries.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\domain\pattern\RectangularSeries.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\domain\pattern\RectangularSeries.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\domain\pattern\TclPatternCommand.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\domain\pattern\TimeSeries.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\domain\pattern\TimeSeries.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\domain\pattern\TrigSeries.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\domain\pattern\TrigSeries.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\domain\pattern\UniformExcitation.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\domain\pattern\UniformExcitation.h
# End Source File
# End Group
# Begin Group "groundMotion"

# PROP Default_Filter ""
# Begin Source File

SOURCE=..\..\..\SRC\domain\groundMotion\GroundMotion.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\domain\groundMotion\GroundMotion.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\domain\groundMotion\GroundMotionIntegrator.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\domain\groundMotion\GroundMotionIntegrator.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\domain\groundMotion\GroundMotionRecord.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\domain\groundMotion\GroundMotionRecord.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\domain\groundMotion\TrapezoidalGroundMotionIntegrator.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\domain\groundMotion\TrapezoidalGroundMotionIntegrator.h
# End Source File
# End Group
# End Target
# End Project
