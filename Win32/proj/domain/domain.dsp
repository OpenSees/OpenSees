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
# PROP Output_Dir "..\..\lib\release"
# PROP Intermediate_Dir "..\..\obj\domain\release"
# PROP Target_Dir ""
# ADD BASE CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_MBCS" /D "_LIB" /YX /FD /c
# ADD CPP /nologo /MLd /W3 /GX /O2 /I "..\..\..\src\reliability\domain\component" /I "..\..\..\src\reliability\domain\filter" /I "..\..\..\src\handler" /I "..\..\..\src\domain\groundMotion" /I "..\..\..\src\analysis\fe_ele" /I "..\..\..\src\utility" /I "..\..\..\src\domain\subdomain" /I "..\..\..\src\database" /I "..\..\..\src\analysis\analysis" /I "..\..\..\src\recorder" /I "..\..\..\src\graph\graph" /I "..\..\..\src\modelbuilder" /I "..\..\..\src\domain\domain\single" /I "..\..\..\src\modelbuilder\tcl" /I "..\..\..\src\tagged\storage" /I "..\..\..\src\domain\constraints" /I "..\..\..\src\domain\pattern" /I "..\..\..\src\matrix" /I "..\..\..\src\analysis\dof_grp" /I "..\..\..\src\domain\domain" /I "..\..\..\src\actor\channel" /I "..\..\..\src" /I "..\..\..\src\actor\actor" /I "..\..\..\src\actor\objectBroker" /I "..\..\..\src\tagged" /I "..\..\..\src\domain\component" /I "..\..\..\src\domain\load" /I "..\..\..\src\domain\node" /I "..\..\..\src\element" /I "..\..\..\src\renderer" /I "..\..\..\src\nDarray" /I "..\..\..\src\domain\region" /I "..\..\..\src\reliability\domain\modulatingFunction" /I "..\..\..\src\reliability\domain\components" /I "..\..\..\src\reliability\domain\spectrum" /I "..\..\..\src\reliability\analysis\randomNumber" /I "..\..\..\src\analysis\algorithm\equiSolnAlgo" /I "..\..\..\src\analysis\algorithm" /I "..\..\..\src\analysis\integrator" /I "..\..\..\src\system_of_eqn\linearSOE" /I "..\..\..\src\system_of_eqn" /I "c:\Program Files\tcl" /D "NDEBUG" /D "WIN32" /D "_LIB" /D "_RELIABILITY" /D "_MBCS" /D "_TCL84" /FD /c
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
# PROP Output_Dir "..\..\lib\debug"
# PROP Intermediate_Dir "..\..\obj\domain\debug"
# PROP Target_Dir ""
# ADD BASE CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_MBCS" /D "_LIB" /YX /FD /GZ /c
# ADD CPP /nologo /W3 /GX /ZI /Od /I "..\..\..\src\reliability\domain\filter" /I "..\..\..\src\handler" /I "..\..\..\src\domain\groundMotion" /I "..\..\..\src\analysis\fe_ele" /I "..\..\..\src\utility" /I "..\..\..\src\domain\subdomain" /I "..\..\..\src\database" /I "..\..\..\src\analysis\analysis" /I "..\..\..\src\recorder" /I "..\..\..\src\graph\graph" /I "..\..\..\src\modelbuilder" /I "..\..\..\src\domain\domain\single" /I "..\..\..\src\modelbuilder\tcl" /I "..\..\..\src\tagged\storage" /I "..\..\..\src\domain\constraints" /I "..\..\..\src\domain\pattern" /I "..\..\..\src\matrix" /I "..\..\..\src\analysis\dof_grp" /I "..\..\..\src\domain\domain" /I "..\..\..\src\actor\channel" /I "..\..\..\src" /I "..\..\..\src\actor\actor" /I "..\..\..\src\actor\objectBroker" /I "..\..\..\src\tagged" /I "..\..\..\src\domain\component" /I "..\..\..\src\domain\load" /I "..\..\..\src\domain\node" /I "..\..\..\src\element" /I "..\..\..\src\renderer" /I "..\..\..\src\nDarray" /I "..\..\..\src\domain\region" /I "..\..\..\src\reliability\domain\modulatingFunction" /I "..\..\..\src\reliability\domain\components" /I "..\..\..\src\reliability\domain\spectrum" /I "..\..\..\src\reliability\analysis\randomNumber" /I "..\..\..\src\analysis\algorithm\equiSolnAlgo" /I "..\..\..\src\analysis\algorithm" /I "..\..\..\src\analysis\integrator" /I "..\..\..\src\system_of_eqn\linearSOE" /I "..\..\..\src\system_of_eqn" /I "c:\Program Files\tcl" /D "_DEBUG" /D "WIN32" /D "_LIB" /D "_RELIABILITY" /D "_MBCS" /D "_TCL84" /FR /FD /GZ /c
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
# Begin Group "beam"

# PROP Default_Filter ""
# Begin Source File

SOURCE=..\..\..\Src\domain\load\Beam2dPointLoad.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\Src\domain\load\Beam2dPointLoad.h
# End Source File
# Begin Source File

SOURCE=..\..\..\Src\domain\load\Beam2dTempLoad.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\Src\domain\load\Beam2dTempLoad.h
# End Source File
# Begin Source File

SOURCE=..\..\..\Src\domain\load\Beam2dUniformLoad.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\Src\domain\load\Beam2dUniformLoad.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\domain\load\Beam3dPointLoad.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\domain\load\Beam3dPointLoad.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\domain\load\Beam3dUniformLoad.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\domain\load\Beam3dUniformLoad.h
# End Source File
# End Group
# Begin Group "brick"

# PROP Default_Filter ""
# Begin Source File

SOURCE=..\..\..\Src\domain\load\BrickSelfWeight.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\Src\domain\load\BrickSelfWeight.h
# End Source File
# End Group
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

SOURCE=..\..\..\SRC\domain\constraints\ImposedMotionSP.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\domain\constraints\ImposedMotionSP.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\domain\constraints\ImposedMotionSP1.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\domain\constraints\ImposedMotionSP1.h
# End Source File
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

SOURCE=..\..\..\SRC\domain\pattern\EarthquakePattern.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\domain\pattern\EarthquakePattern.h
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

SOURCE=..\..\..\SRC\domain\pattern\MultiSupportPattern.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\domain\pattern\MultiSupportPattern.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\domain\pattern\PBowlLoading.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\domain\pattern\PBowlLoading.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\domain\pattern\TclPatternCommand.cpp
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

SOURCE=..\..\..\SRC\domain\groundMotion\GroundMotionRecord.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\domain\groundMotion\GroundMotionRecord.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\domain\groundMotion\InterpolatedGroundMotion.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\domain\groundMotion\InterpolatedGroundMotion.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\domain\groundMotion\TclGroundMotionCommand.cpp
# End Source File
# End Group
# Begin Group "timeSeries"

# PROP Default_Filter ""
# Begin Source File

SOURCE=..\..\..\SRC\domain\pattern\ConstantSeries.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\domain\pattern\ConstantSeries.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\domain\pattern\DiscretizedRandomProcessSeries.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\domain\pattern\DiscretizedRandomProcessSeries.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\domain\pattern\LinearSeries.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\domain\pattern\LinearSeries.h
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

SOURCE=..\..\..\SRC\domain\pattern\SimulatedRandomProcessSeries.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\domain\pattern\SimulatedRandomProcessSeries.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\domain\pattern\TclSeriesCommand.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\domain\pattern\TclSeriesIntegratorCommand.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\domain\pattern\TimeSeries.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\domain\pattern\TimeSeries.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\domain\pattern\TimeSeriesIntegrator.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\domain\pattern\TimeSeriesIntegrator.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\domain\pattern\TrapezoidalTimeSeriesIntegrator.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\domain\pattern\TrapezoidalTimeSeriesIntegrator.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\domain\pattern\TrigSeries.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\domain\pattern\TrigSeries.h
# End Source File
# End Group
# Begin Group "region"

# PROP Default_Filter ""
# Begin Source File

SOURCE=..\..\..\SRC\domain\region\MeshRegion.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\domain\region\MeshRegion.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\domain\region\TclRegionCommands.cpp
# End Source File
# End Group
# End Target
# End Project
