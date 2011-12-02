# Microsoft Developer Studio Project File - Name="analysis" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Static Library" 0x0104

CFG=analysis - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "analysis.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "analysis.mak" CFG="analysis - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "analysis - Win32 Release" (based on "Win32 (x86) Static Library")
!MESSAGE "analysis - Win32 Debug" (based on "Win32 (x86) Static Library")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
RSC=rc.exe

!IF  "$(CFG)" == "analysis - Win32 Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "Release"
# PROP BASE Intermediate_Dir "Release"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "..\..\lib"
# PROP Intermediate_Dir "..\..\obj\analysis"
# PROP Target_Dir ""
# ADD BASE CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_MBCS" /D "_LIB" /YX /FD /c
# ADD CPP /nologo /MLd /W3 /GX /O2 /I "..\..\..\src\analysis\algorithm\domainDecompAlgo" /I "..\..\..\src\analysis\algorithm\eigenAlgo" /I "..\..\..\src\system_of_eqn\eigenSOE" /I "..\..\..\src\handler" /I "..\..\..\src\domain\domain\single" /I "..\..\..\src\convergenceTest" /I "..\..\..\src\analysis\analysis" /I "..\..\..\src\recorder" /I "..\..\..\src\analysis\integrator" /I "..\..\..\src\analysis\algorithm\equiSolnAlgo" /I "..\..\..\src\analysis\algorithm" /I "..\..\..\src\system_of_eqn" /I "..\..\..\src\system_of_eqn\linearSOE" /I "..\..\..\src\graph\graph" /I "..\..\..\src\graph\numberer" /I "..\..\..\src\analysis\numberer" /I "..\..\..\src\analysis\fe_ele\transformation" /I "..\..\..\src\analysis\fe_ele\lagrange" /I "..\..\..\src\analysis\fe_ele\penalty" /I "..\..\..\src\actor\objectBroker" /I "..\..\..\src\actor\channel" /I "..\..\..\src\utility" /I "..\..\..\src\domain\subdomain" /I "..\..\..\src\domain\constraints" /I "..\..\..\src\tagged" /I "..\..\..\src\domain\component" /I "..\..\..\src\domain\node" /I "..\..\..\src\element" /I "..\..\..\src\analysis\dof_grp" /I "..\..\..\src\matrix" /I "..\..\..\src\analysis\model\simple" /I "..\..\..\src\domain\domain" /I "..\..\..\src\analysis\model" /I "..\..\..\src" /I "..\..\..\src\actor\actor" /I "..\..\..\src\analysis\handler" /I "..\..\..\src\analysis\fe_ele" /D "WIN32" /D "NDEBUG" /D "_MBCS" /D "_LIB" /FD /c
# SUBTRACT CPP /YX
# ADD BASE RSC /l 0x409 /d "NDEBUG"
# ADD RSC /l 0x409 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo

!ELSEIF  "$(CFG)" == "analysis - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "Debug"
# PROP BASE Intermediate_Dir "Debug"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "..\..\lib"
# PROP Intermediate_Dir "..\..\obj\analysis"
# PROP Target_Dir ""
# ADD BASE CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_MBCS" /D "_LIB" /YX /FD /GZ /c
# ADD CPP /nologo /W3 /GX /ZI /Od /I "..\..\..\src\analysis\algorithm\domainDecompAlgo" /I "..\..\..\src\analysis\algorithm\eigenAlgo" /I "..\..\..\src\system_of_eqn\eigenSOE" /I "..\..\..\src\handler" /I "..\..\..\src\domain\domain\single" /I "..\..\..\src\convergenceTest" /I "..\..\..\src\analysis\analysis" /I "..\..\..\src\recorder" /I "..\..\..\src\analysis\integrator" /I "..\..\..\src\analysis\algorithm\equiSolnAlgo" /I "..\..\..\src\analysis\algorithm" /I "..\..\..\src\system_of_eqn" /I "..\..\..\src\system_of_eqn\linearSOE" /I "..\..\..\src\graph\graph" /I "..\..\..\src\graph\numberer" /I "..\..\..\src\analysis\numberer" /I "..\..\..\src\analysis\fe_ele\transformation" /I "..\..\..\src\analysis\fe_ele\lagrange" /I "..\..\..\src\analysis\fe_ele\penalty" /I "..\..\..\src\actor\objectBroker" /I "..\..\..\src\actor\channel" /I "..\..\..\src\utility" /I "..\..\..\src\domain\subdomain" /I "..\..\..\src\domain\constraints" /I "..\..\..\src\tagged" /I "..\..\..\src\domain\component" /I "..\..\..\src\domain\node" /I "..\..\..\src\element" /I "..\..\..\src\analysis\dof_grp" /I "..\..\..\src\matrix" /I "..\..\..\src\analysis\model\simple" /I "..\..\..\src\domain\domain" /I "..\..\..\src\analysis\model" /I "..\..\..\src" /I "..\..\..\src\actor\actor" /I "..\..\..\src\analysis\handler" /I "..\..\..\src\analysis\fe_ele" /D "WIN32" /D "_DEBUG" /D "_MBCS" /D "_LIB" /FR /FD /GZ /c
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

# Name "analysis - Win32 Release"
# Name "analysis - Win32 Debug"
# Begin Group "analysis"

# PROP Default_Filter ".cpp; .h"
# Begin Source File

SOURCE=..\..\..\SRC\analysis\analysis\Analysis.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\analysis\analysis\Analysis.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\analysis\analysis\DirectIntegrationAnalysis.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\analysis\analysis\DirectIntegrationAnalysis.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\analysis\analysis\DomainDecompositionAnalysis.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\analysis\analysis\DomainDecompositionAnalysis.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\analysis\analysis\EigenAnalysis.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\analysis\analysis\EigenAnalysis.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\analysis\analysis\StaticAnalysis.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\analysis\analysis\StaticAnalysis.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\analysis\analysis\TransientAnalysis.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\analysis\analysis\TransientAnalysis.h
# End Source File
# End Group
# Begin Group "model"

# PROP Default_Filter ".cpp; .h"
# Begin Source File

SOURCE=..\..\..\SRC\analysis\model\AnalysisModel.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\analysis\model\AnalysisModel.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\analysis\model\DOF_GrpIter.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\analysis\model\FE_EleIter.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\analysis\model\simple\SimpleDOF_Iter.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\analysis\model\simple\SimpleDOF_Iter.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\analysis\model\simple\SimpleFE_Iter.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\analysis\model\simple\SimpleFE_Iter.h
# End Source File
# End Group
# Begin Group "fe_ele"

# PROP Default_Filter ".cpp; .h"
# Begin Source File

SOURCE=..\..\..\SRC\analysis\fe_ele\FE_Element.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\analysis\fe_ele\FE_Element.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\analysis\fe_ele\lagrange\LagrangeMP_FE.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\analysis\fe_ele\lagrange\LagrangeMP_FE.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\analysis\fe_ele\lagrange\LagrangeSP_FE.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\analysis\fe_ele\lagrange\LagrangeSP_FE.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\analysis\fe_ele\penalty\PenaltyMP_FE.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\analysis\fe_ele\penalty\PenaltyMP_FE.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\analysis\fe_ele\penalty\PenaltySP_FE.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\analysis\fe_ele\penalty\PenaltySP_FE.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\analysis\fe_ele\transformation\TransformationFE.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\analysis\fe_ele\transformation\TransformationFE.h
# End Source File
# End Group
# Begin Group "dof_grp"

# PROP Default_Filter ".cpp; .h"
# Begin Source File

SOURCE=..\..\..\SRC\analysis\dof_grp\DOF_Group.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\analysis\dof_grp\DOF_Group.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\analysis\dof_grp\LagrangeDOF_Group.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\analysis\dof_grp\LagrangeDOF_Group.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\analysis\dof_grp\TransformationDOF_Group.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\analysis\dof_grp\TransformationDOF_Group.h
# End Source File
# End Group
# Begin Group "algorithm"

# PROP Default_Filter ".cpp; .h"
# Begin Source File

SOURCE=..\..\..\SRC\analysis\algorithm\domainDecompAlgo\DomainDecompAlgo.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\analysis\algorithm\domainDecompAlgo\DomainDecompAlgo.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\analysis\algorithm\eigenAlgo\EigenAlgorithm.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\analysis\algorithm\eigenAlgo\EigenAlgorithm.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\analysis\algorithm\equiSolnAlgo\EquiSolnAlgo.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\analysis\algorithm\equiSolnAlgo\EquiSolnAlgo.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\analysis\algorithm\eigenAlgo\FrequencyAlgo.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\analysis\algorithm\eigenAlgo\FrequencyAlgo.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\analysis\algorithm\equiSolnAlgo\Linear.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\analysis\algorithm\equiSolnAlgo\Linear.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\analysis\algorithm\equiSolnAlgo\ModifiedNewton.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\analysis\algorithm\equiSolnAlgo\ModifiedNewton.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\analysis\algorithm\equiSolnAlgo\NewtonRaphson.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\analysis\algorithm\equiSolnAlgo\NewtonRaphson.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\analysis\algorithm\SolutionAlgorithm.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\analysis\algorithm\SolutionAlgorithm.h
# End Source File
# End Group
# Begin Group "integrator"

# PROP Default_Filter ".cpp; .h"
# Begin Source File

SOURCE=..\..\..\SRC\analysis\integrator\ArcLength.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\analysis\integrator\ArcLength.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\analysis\integrator\ArcLength1.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\analysis\integrator\ArcLength1.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\analysis\integrator\CentralDifference.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\analysis\integrator\CentralDifference.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\analysis\integrator\DisplacementControl.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\analysis\integrator\DisplacementControl.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\analysis\integrator\EigenIntegrator.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\analysis\integrator\EigenIntegrator.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\analysis\integrator\HHT.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\analysis\integrator\HHT.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\analysis\integrator\HHT1.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\analysis\integrator\HHT1.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\analysis\integrator\IncrementalIntegrator.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\analysis\integrator\IncrementalIntegrator.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\analysis\integrator\Integrator.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\analysis\integrator\Integrator.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\analysis\integrator\LoadControl.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\analysis\integrator\LoadControl.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\analysis\integrator\LoadPath.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\analysis\integrator\LoadPath.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\analysis\integrator\MinUnbalDispNorm.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\analysis\integrator\MinUnbalDispNorm.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\analysis\integrator\Newmark.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\analysis\integrator\Newmark.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\analysis\integrator\Newmark1.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\analysis\integrator\Newmark1.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\analysis\integrator\StaticIntegrator.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\analysis\integrator\StaticIntegrator.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\analysis\integrator\TransientIntegrator.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\analysis\integrator\TransientIntegrator.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\analysis\integrator\WilsonTheta.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\analysis\integrator\WilsonTheta.h
# End Source File
# End Group
# Begin Group "numberer"

# PROP Default_Filter ".cpp; .h"
# Begin Source File

SOURCE=..\..\..\SRC\analysis\numberer\DOF_Numberer.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\analysis\numberer\DOF_Numberer.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\analysis\numberer\PlainNumberer.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\analysis\numberer\PlainNumberer.h
# End Source File
# End Group
# Begin Group "handler"

# PROP Default_Filter ".cpp; .h"
# Begin Source File

SOURCE=..\..\..\SRC\analysis\handler\ConstraintHandler.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\analysis\handler\ConstraintHandler.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\analysis\handler\LagrangeConstraintHandler.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\analysis\handler\LagrangeConstraintHandler.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\analysis\handler\PenaltyConstraintHandler.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\analysis\handler\PenaltyConstraintHandler.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\analysis\handler\PlainHandler.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\analysis\handler\PlainHandler.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\analysis\handler\TransformationConstraintHandler.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\analysis\handler\TransformationConstraintHandler.h
# End Source File
# End Group
# End Target
# End Project
