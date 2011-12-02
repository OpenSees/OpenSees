# Microsoft Developer Studio Project File - Name="element" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Static Library" 0x0104

CFG=element - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "element.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "element.mak" CFG="element - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "element - Win32 Release" (based on "Win32 (x86) Static Library")
!MESSAGE "element - Win32 Debug" (based on "Win32 (x86) Static Library")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
RSC=rc.exe

!IF  "$(CFG)" == "element - Win32 Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "Release"
# PROP BASE Intermediate_Dir "Release"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "..\..\lib\release"
# PROP Intermediate_Dir "..\..\obj\element\release"
# PROP Target_Dir ""
# ADD BASE CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_MBCS" /D "_LIB" /YX /FD /c
# ADD CPP /nologo /MLd /W3 /GX /O2 /I "c:\Program Files\tcl\include" /I "..\..\..\src\element\upU" /I "..\..\..\src\element\dispBeamColumn" /I "..\..\..\src\element\brick" /I "..\..\..\src\element\shell" /I "..\..\..\src\element\8nbrick" /I "..\..\..\src\recorder\response" /I "..\..\..\src\element\nonlinearBeamColumn\quadRule" /I "..\..\..\src\material\nD" /I "..\..\..\src\element\fourNodeQuad" /I "..\..\..\src\element\damper" /I "..\..\..\src\coordTransformation" /I "..\..\..\src\element\beamWithHinges" /I "..\..\..\src\element\nonlinearBeamColumn\matrixutil" /I "..\..\..\src\element\zeroLength" /I "..\..\..\src\modelbuilder" /I "..\..\..\src\modelbuilder\tcl" /I "..\..\..\src\element\feap" /I "..\..\..\src\handler" /I "..\..\..\src\element" /I "..\..\..\src\element\truss" /I "..\..\..\src\material\section" /I "..\..\..\src\element\beam3d" /I "..\..\..\src\element\beam2d" /I "..\..\..\src\material" /I "..\..\..\src\material\uniaxial" /I "..\..\..\src\actor\objectBroker" /I "..\..\..\src\matrix" /I "..\..\..\src\domain\load" /I "..\..\..\src\renderer" /I "..\..\..\src\actor\channel" /I "..\..\..\src\domain\node" /I "..\..\..\src\actor\actor" /I "..\..\..\src\tagged" /I "..\..\..\src\domain\component" /I "..\..\..\src" /I "..\..\..\src\domain\domain" /I "..\..\..\src\material\nd\template3dep" /I "..\..\..\src\nDarray" /I "..\..\..\src\element\20nbrick" /I "..\..\..\src\element\elasticBeamColumn" /I "..\..\..\src\element\joint" /I "..\..\..\src\domain\constraints" /I "..\..\..\src\element\updatedLagrangianBeamColumn" /I "..\..\..\src\material\yieldSurface\yieldSurfaceBC" /I "..\..\..\src\material\yieldSurface\evolution" /I "..\..\..\src\element\forceBeamColumn" /I "..\..\..\src\element\nonlinearBeamColumn\element" /I "..\..\..\src\element\nonlinearBeamColumn\matrixUtil" /D "NDEBUG" /D "WIN32" /D "_MBCS" /D "_LIB" /D "_WGL" /D "_COROTATIONAL" /FD /c
# SUBTRACT CPP /YX
# ADD BASE RSC /l 0x409 /d "NDEBUG"
# ADD RSC /l 0x409 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo

!ELSEIF  "$(CFG)" == "element - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "Debug"
# PROP BASE Intermediate_Dir "Debug"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "..\..\lib\debug"
# PROP Intermediate_Dir "..\..\obj\element\debug"
# PROP Target_Dir ""
# ADD BASE CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_MBCS" /D "_LIB" /YX /FD /GZ /c
# ADD CPP /nologo /W3 /GX /ZI /Od /I "..\..\..\src\element\material\nd\template3dep" /I "c:\Program Files\tcl\include" /I "..\..\..\src\element\upU" /I "..\..\..\src\element\dispBeamColumn" /I "..\..\..\src\element\brick" /I "..\..\..\src\element\shell" /I "..\..\..\src\element\8nbrick" /I "..\..\..\src\recorder\response" /I "..\..\..\src\element\nonlinearBeamColumn\quadRule" /I "..\..\..\src\material\nD" /I "..\..\..\src\element\fourNodeQuad" /I "..\..\..\src\element\damper" /I "..\..\..\src\coordTransformation" /I "..\..\..\src\element\beamWithHinges" /I "..\..\..\src\element\nonlinearBeamColumn\matrixutil" /I "..\..\..\src\element\zeroLength" /I "..\..\..\src\modelbuilder" /I "..\..\..\src\modelbuilder\tcl" /I "..\..\..\src\element\feap" /I "..\..\..\src\handler" /I "..\..\..\src\element" /I "..\..\..\src\element\truss" /I "..\..\..\src\material\section" /I "..\..\..\src\element\beam3d" /I "..\..\..\src\element\beam2d" /I "..\..\..\src\material" /I "..\..\..\src\material\uniaxial" /I "..\..\..\src\actor\objectBroker" /I "..\..\..\src\matrix" /I "..\..\..\src\domain\load" /I "..\..\..\src\renderer" /I "..\..\..\src\actor\channel" /I "..\..\..\src\domain\node" /I "..\..\..\src\actor\actor" /I "..\..\..\src\tagged" /I "..\..\..\src\domain\component" /I "..\..\..\src" /I "..\..\..\src\domain\domain" /I "..\..\..\src\material\nd\template3dep" /I "..\..\..\src\nDarray" /I "..\..\..\src\element\20nbrick" /I "..\..\..\src\element\elasticBeamColumn" /I "..\..\..\src\element\joint" /I "..\..\..\src\domain\constraints" /I "..\..\..\src\element\updatedLagrangianBeamColumn" /I "..\..\..\src\material\yieldSurface\yieldSurfaceBC" /I "..\..\..\src\material\yieldSurface\evolution" /I "..\..\..\src\element\forceBeamColumn" /I "..\..\..\src\element\nonlinearBeamColumn\element" /I "..\..\..\src\element\nonlinearBeamColumn\matrixUtil" /D "_DEBUG" /D "WIN32" /D "_MBCS" /D "_LIB" /D "_WGL" /D "_COROTATIONAL" /FR /FD /GZ /c
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

# Name "element - Win32 Release"
# Name "element - Win32 Debug"
# Begin Group "Source Files"

# PROP Default_Filter "cpp;c;cxx;rc;def;r;odl;idl;hpj;bat"
# Begin Source File

SOURCE=..\..\..\SRC\element\Element.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\element\ElementalLoad.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\element\Information.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\element\TclElementCommands.cpp
# End Source File
# End Group
# Begin Group "Header Files"

# PROP Default_Filter "h;hpp;hxx;hm;inl"
# Begin Source File

SOURCE=..\..\..\SRC\element\Element.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\element\ElementalLoad.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\element\Information.h
# End Source File
# End Group
# Begin Group "truss"

# PROP Default_Filter ""
# Begin Source File

SOURCE=..\..\..\SRC\element\truss\CorotTruss.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\element\truss\CorotTruss.h
# End Source File
# Begin Source File

SOURCE=..\..\..\Src\element\truss\CorotTrussSection.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\Src\element\truss\CorotTrussSection.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\element\truss\TclTrussCommand.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\element\truss\Truss.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\element\truss\Truss.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\element\truss\TrussSection.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\element\truss\TrussSection.h
# End Source File
# End Group
# Begin Group "feap"

# PROP Default_Filter ""
# Begin Source File

SOURCE=..\..\..\SRC\element\feap\fElement.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\element\feap\fElement.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\element\feap\fElmt02.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\element\feap\fElmt02.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\element\feap\TclFeapElementCommand.cpp
# End Source File
# End Group
# Begin Group "zeroLength"

# PROP Default_Filter ""
# Begin Source File

SOURCE=..\..\..\SRC\element\zeroLength\TclZeroLength.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\element\zeroLength\ZeroLength.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\element\zeroLength\ZeroLength.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\element\zeroLength\ZeroLengthSection.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\element\zeroLength\ZeroLengthSection.h
# End Source File
# End Group
# Begin Group "beamWithHinges"

# PROP Default_Filter ""
# Begin Source File

SOURCE=..\..\..\SRC\element\beamWithHinges\BeamWithHinges2d.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\element\beamWithHinges\BeamWithHinges2d.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\element\beamWithHinges\BeamWithHinges3d.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\element\beamWithHinges\BeamWithHinges3d.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\element\beamWithHinges\TclBeamWithHingesBuilder.cpp
# End Source File
# End Group
# Begin Group "crdTransf"

# PROP Default_Filter ""
# Begin Source File

SOURCE=..\..\..\SRC\coordTransformation\CrdTransf.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\coordTransformation\CrdTransf.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\coordTransformation\CrdTransf2d.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\coordTransformation\CrdTransf2d.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\coordTransformation\CrdTransf3d.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\coordTransformation\CrdTransf3d.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\coordTransformation\LinearCrdTransf2d.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\coordTransformation\LinearCrdTransf2d.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\coordTransformation\LinearCrdTransf3d.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\coordTransformation\LinearCrdTransf3d.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\coordTransformation\PDeltaCrdTransf2d.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\coordTransformation\PDeltaCrdTransf2d.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\coordTransformation\PDeltaCrdTransf3d.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\coordTransformation\PDeltaCrdTransf3d.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\coordTransformation\TclGeomTransfCommand.cpp
# End Source File
# End Group
# Begin Group "quad"

# PROP Default_Filter ""
# Begin Source File

SOURCE=..\..\..\SRC\element\fourNodeQuad\ConstantPressureVolumeQuad.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\element\fourNodeQuad\ConstantPressureVolumeQuad.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\element\fourNodeQuad\EnhancedQuad.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\element\fourNodeQuad\EnhancedQuad.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\element\fourNodeQuad\FourNodeQuad.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\element\fourNodeQuad\FourNodeQuad.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\element\fourNodeQuad\NineNodeMixedQuad.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\element\fourNodeQuad\NineNodeMixedQuad.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\element\fourNodeQuad\TclFourNodeQuadCommand.cpp
# End Source File
# End Group
# Begin Group "elasticBeamColumn"

# PROP Default_Filter ""
# Begin Source File

SOURCE=..\..\..\SRC\element\elasticBeamColumn\ElasticBeam2d.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\element\elasticBeamColumn\ElasticBeam2d.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\element\elasticBeamColumn\ElasticBeam3d.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\element\elasticBeamColumn\ElasticBeam3d.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\element\elasticBeamColumn\TclElasticBeamCommand.cpp
# End Source File
# End Group
# Begin Group "eightNodeBrick"

# PROP Default_Filter ""
# Begin Source File

SOURCE=..\..\..\SRC\element\brick\BbarBrick.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\element\brick\BbarBrick.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\element\brick\Brick.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\element\brick\Brick.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\element\8nbrick\EightNodeBrick.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\element\8nbrick\EightNodeBrick.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\element\brick\shp3d.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\element\brick\TclBrickCommand.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\element\8nbrick\TclEightNodeBrickCommand.cpp
# End Source File
# End Group
# Begin Group "shell"

# PROP Default_Filter ""
# Begin Source File

SOURCE=..\..\..\SRC\element\shell\R3vectors.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\element\shell\R3vectors.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\element\shell\ShellMITC4.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\element\shell\ShellMITC4.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\element\shell\TclShellCommand.cpp
# End Source File
# End Group
# Begin Group "dispBeamColumn"

# PROP Default_Filter ""
# Begin Source File

SOURCE=..\..\..\SRC\element\dispBeamColumn\DispBeamColumn2d.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\element\dispBeamColumn\DispBeamColumn2d.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\element\dispBeamColumn\DispBeamColumn3d.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\element\dispBeamColumn\DispBeamColumn3d.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\element\dispBeamColumn\TclDispBeamColumnCommand.cpp
# End Source File
# End Group
# Begin Group "twentyNodeBrick"

# PROP Default_Filter ""
# Begin Source File

SOURCE=..\..\..\SRC\element\20nbrick\TclTwentyNodeBrickCommand.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\element\20nbrick\TwentyNodeBrick.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\element\20nbrick\TwentyNodeBrick.h
# End Source File
# End Group
# Begin Group "upU"

# PROP Default_Filter ""
# Begin Source File

SOURCE=..\..\..\Src\element\upU\EightNodeBrick_u_p_U.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\Src\element\upU\EightNodeBrick_u_p_U.h
# End Source File
# Begin Source File

SOURCE=..\..\..\Src\element\upU\FourNodeQuadUP.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\Src\element\upU\FourNodeQuadUP.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\element\upU\TclEightNodeBrick_u_p_U.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\Src\element\upU\TclFourNodeQuadUPCommand.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\element\upU\TclTwentyNodeBrick_u_p_U.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\Src\element\upU\TwentyNodeBrick_u_p_U.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\Src\element\upU\TwentyNodeBrick_u_p_U.h
# End Source File
# End Group
# Begin Group "joint"

# PROP Default_Filter ""
# Begin Source File

SOURCE=..\..\..\Src\element\joint\Joint2D.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\Src\element\joint\Joint2D.h
# End Source File
# Begin Source File

SOURCE=..\..\..\Src\element\joint\MP_Joint2D.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\Src\element\joint\MP_Joint2D.h
# End Source File
# Begin Source File

SOURCE=..\..\..\Src\element\joint\TclJoint2dCommand.cpp
# End Source File
# End Group
# Begin Group "ulBeamColumn"

# PROP Default_Filter ""
# Begin Source File

SOURCE=..\..\..\SRC\element\updatedLagrangianBeamColumn\BilinearCyclic.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\element\updatedLagrangianBeamColumn\BilinearCyclic.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\element\updatedLagrangianBeamColumn\CyclicModel.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\element\updatedLagrangianBeamColumn\CyclicModel.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\element\updatedLagrangianBeamColumn\Elastic2DGNL.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\element\updatedLagrangianBeamColumn\Elastic2DGNL.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\element\updatedLagrangianBeamColumn\Inelastic2DYS01.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\element\updatedLagrangianBeamColumn\Inelastic2DYS01.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\element\updatedLagrangianBeamColumn\Inelastic2DYS02.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\element\updatedLagrangianBeamColumn\Inelastic2DYS02.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\element\updatedLagrangianBeamColumn\Inelastic2DYS03.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\element\updatedLagrangianBeamColumn\Inelastic2DYS03.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\element\updatedLagrangianBeamColumn\InelasticYS2DGNL.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\element\updatedLagrangianBeamColumn\InelasticYS2DGNL.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\element\updatedLagrangianBeamColumn\LinearCyclic.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\element\updatedLagrangianBeamColumn\LinearCyclic.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\element\updatedLagrangianBeamColumn\QuadraticCyclic.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\element\updatedLagrangianBeamColumn\QuadraticCyclic.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\element\updatedLagrangianBeamColumn\TclCyclicModelCommands.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\element\updatedLagrangianBeamColumn\TclElement2dGNL.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\element\updatedLagrangianBeamColumn\TclElement2dYS.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\element\updatedLagrangianBeamColumn\UpdatedLagrangianBeam2D.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\element\updatedLagrangianBeamColumn\UpdatedLagrangianBeam2D.h
# End Source File
# End Group
# Begin Group "forceBeamColumn"

# PROP Default_Filter ""
# Begin Group "beamIntegration"

# PROP Default_Filter ""
# Begin Source File

SOURCE=..\..\..\SRC\element\forceBeamColumn\BeamIntegration.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\element\forceBeamColumn\BeamIntegration.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\element\forceBeamColumn\HingeMidpointBeamIntegration2d.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\element\forceBeamColumn\HingeMidpointBeamIntegration2d.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\element\forceBeamColumn\HingeMidpointBeamIntegration3d.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\element\forceBeamColumn\HingeMidpointBeamIntegration3d.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\element\forceBeamColumn\HingeRadauBeamIntegration2d.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\element\forceBeamColumn\HingeRadauBeamIntegration2d.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\element\forceBeamColumn\HingeRadauBeamIntegration3d.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\element\forceBeamColumn\HingeRadauBeamIntegration3d.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\element\forceBeamColumn\HingeRadauTwoBeamIntegration2d.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\element\forceBeamColumn\HingeRadauTwoBeamIntegration2d.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\element\forceBeamColumn\HingeRadauTwoBeamIntegration3d.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\element\forceBeamColumn\HingeRadauTwoBeamIntegration3d.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\element\forceBeamColumn\LobattoBeamIntegration.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\element\forceBeamColumn\LobattoBeamIntegration.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\element\forceBeamColumn\UserDefinedBeamIntegration.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\element\forceBeamColumn\UserDefinedBeamIntegration.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\element\forceBeamColumn\UserDefinedHingeIntegration2d.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\element\forceBeamColumn\UserDefinedHingeIntegration2d.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\element\forceBeamColumn\UserDefinedHingeIntegration3d.cpp
# End Source File
# End Group
# Begin Source File

SOURCE=..\..\..\SRC\element\forceBeamColumn\ForceBeamColumn2d.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\element\forceBeamColumn\ForceBeamColumn2d.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\element\forceBeamColumn\ForceBeamColumn3d.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\element\forceBeamColumn\ForceBeamColumn3d.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\element\forceBeamColumn\TclForceBeamColumnCommand.cpp
# End Source File
# End Group
# Begin Group "nonlinearBeamColumn"

# PROP Default_Filter ""
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
# Begin Source File

SOURCE=..\..\..\SRC\element\nonlinearBeamColumn\element\TclNLBeamColumnCommand.cpp
# End Source File
# End Group
# Begin Source File

SOURCE=..\..\..\SRC\element\NewElement.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\element\NewElement.h
# End Source File
# End Target
# End Project
