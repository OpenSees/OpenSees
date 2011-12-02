# Microsoft Developer Studio Project File - Name="material" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Static Library" 0x0104

CFG=material - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "material.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "material.mak" CFG="material - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "material - Win32 Release" (based on "Win32 (x86) Static Library")
!MESSAGE "material - Win32 Debug" (based on "Win32 (x86) Static Library")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
RSC=rc.exe

!IF  "$(CFG)" == "material - Win32 Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "Release"
# PROP BASE Intermediate_Dir "Release"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "..\..\lib\release"
# PROP Intermediate_Dir "..\..\obj\material\release"
# PROP Target_Dir ""
# ADD BASE CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_MBCS" /D "_LIB" /YX /FD /c
# ADD CPP /nologo /MLd /W3 /GX /O2 /I "c:\Program Files\tcl\include" /I "..\..\..\src\material\uniaxial\fedeas" /I "..\..\..\src\material\uniaxial\drain" /I "..\..\..\src\domain\domain" /I "..\..\..\src\renderer" /I "..\..\..\src\material\nD\soil" /I "..\..\..\src\material\nD\template3dep" /I "..\..\..\src\recorder\response" /I "..\..\..\src\material\backbone" /I "..\..\..\src\material\state" /I "..\..\..\src\material\state\strength" /I "..\..\..\src\material\state\deformation" /I "..\..\..\src\material\state\stiffness" /I "..\..\..\src\material\section\repres\section" /I "..\..\..\src\material\section\repres\cell" /I "..\..\..\src\material\section\repres\patch" /I "..\..\..\src\material\section\repres\reinfBar" /I "..\..\..\src\material\section\repres\reinfLayer" /I "..\..\..\src\material\section\fiber" /I "..\..\..\src\element\nonlinearBeamColumn\fiber" /I "..\..\..\src\element\nonlinearBeamColumn\matrixutil" /I "..\..\..\src\material\section" /I "..\..\..\src\handler" /I "..\..\..\src\material\uniaxial" /I "..\..\..\src\material\nD" /I "..\..\..\src\element" /I "..\..\..\src\actor\channel" /I "..\..\..\src\actor\objectBroker" /I "..\..\..\src\matrix" /I "..\..\..\src" /I "..\..\..\src\actor\actor" /I "..\..\..\src\tagged" /I "..\..\..\src\modelbuilder" /I "..\..\..\src\domain\component" /I "..\..\..\src\material" /I "..\..\..\src\modelbuilder\tcl" /I "..\..\..\src\material\nd\template3dep" /I "..\..\..\src\nDarray" /I "..\..\..\src\material\uniaxial\py" /I "..\..\..\src\material\uniaxial\snap" /I "..\..\..\src\material\yieldSurface\yieldSurfaceBC" /I "..\..\..\src\material\yieldSurface\evolution" /I "..\..\..\src\material\yieldSurface\plasticHardeningMaterial" /I "..\..\..\src\material\section\yieldSurface" /I "..\..\..\src\material\nd\feap" /D "WIN32" /D "NDEBUG" /D "_MBCS" /D "_LIB" /FD /c
# SUBTRACT CPP /YX
# ADD BASE RSC /l 0x409 /d "NDEBUG"
# ADD RSC /l 0x409 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo

!ELSEIF  "$(CFG)" == "material - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "Debug"
# PROP BASE Intermediate_Dir "Debug"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "..\..\lib\debug"
# PROP Intermediate_Dir "..\..\obj\material\debug"
# PROP Target_Dir ""
# ADD BASE CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_MBCS" /D "_LIB" /YX /FD /GZ /c
# ADD CPP /nologo /W3 /GX /ZI /Od /I "c:\Program Files\tcl\include" /I "..\..\..\src\material\uniaxial\fedeas" /I "..\..\..\src\material\uniaxial\drain" /I "..\..\..\src\domain\domain" /I "..\..\..\src\renderer" /I "..\..\..\src\material\nD\soil" /I "..\..\..\src\material\nD\template3dep" /I "..\..\..\src\recorder\response" /I "..\..\..\src\material\backbone" /I "..\..\..\src\material\state" /I "..\..\..\src\material\state\strength" /I "..\..\..\src\material\state\deformation" /I "..\..\..\src\material\state\stiffness" /I "..\..\..\src\material\section\repres\section" /I "..\..\..\src\material\section\repres\cell" /I "..\..\..\src\material\section\repres\patch" /I "..\..\..\src\material\section\repres\reinfBar" /I "..\..\..\src\material\section\repres\reinfLayer" /I "..\..\..\src\material\section\fiber" /I "..\..\..\src\element\nonlinearBeamColumn\fiber" /I "..\..\..\src\element\nonlinearBeamColumn\matrixutil" /I "..\..\..\src\material\section" /I "..\..\..\src\handler" /I "..\..\..\src\material\uniaxial" /I "..\..\..\src\material\nD" /I "..\..\..\src\element" /I "..\..\..\src\actor\channel" /I "..\..\..\src\actor\objectBroker" /I "..\..\..\src\matrix" /I "..\..\..\src" /I "..\..\..\src\actor\actor" /I "..\..\..\src\tagged" /I "..\..\..\src\modelbuilder" /I "..\..\..\src\domain\component" /I "..\..\..\src\material" /I "..\..\..\src\modelbuilder\tcl" /I "..\..\..\src\material\nd\template3dep" /I "..\..\..\src\nDarray" /I "..\..\..\src\material\uniaxial\py" /I "..\..\..\src\material\uniaxial\snap" /I "..\..\..\src\material\yieldSurface\yieldSurfaceBC" /I "..\..\..\src\material\yieldSurface\evolution" /I "..\..\..\src\material\yieldSurface\plasticHardeningMaterial" /I "..\..\..\src\material\section\yieldSurface" /I "..\..\..\src\material\nd\feap" /D "WIN32" /D "_DEBUG" /D "_MBCS" /D "_LIB" /D "MHS_INSURE" /FR /FD /GZ /c
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

# Name "material - Win32 Release"
# Name "material - Win32 Debug"
# Begin Group "uniaxial"

# PROP Default_Filter ""
# Begin Group "fedeas"

# PROP Default_Filter ""
# Begin Source File

SOURCE=..\..\..\SRC\material\uniaxial\fedeas\FedeasBond1Material.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\uniaxial\fedeas\FedeasBond1Material.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\uniaxial\fedeas\FedeasBond2Material.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\uniaxial\fedeas\FedeasBond2Material.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\uniaxial\fedeas\FedeasConcr1Material.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\uniaxial\fedeas\FedeasConcr1Material.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\uniaxial\fedeas\FedeasConcr2Material.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\uniaxial\fedeas\FedeasConcr2Material.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\uniaxial\fedeas\FedeasConcr3Material.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\uniaxial\fedeas\FedeasConcr3Material.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\uniaxial\fedeas\FedeasHardeningMaterial.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\uniaxial\fedeas\FedeasHardeningMaterial.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\uniaxial\fedeas\FedeasHyster1Material.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\uniaxial\fedeas\FedeasHyster1Material.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\uniaxial\fedeas\FedeasHyster2Material.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\uniaxial\fedeas\FedeasHyster2Material.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\uniaxial\fedeas\FedeasSteel1Material.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\uniaxial\fedeas\FedeasSteel1Material.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\uniaxial\fedeas\FedeasSteel2Material.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\uniaxial\fedeas\FedeasSteel2Material.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\uniaxial\fedeas\TclFedeasMaterialCommand.cpp
# End Source File
# End Group
# Begin Group "drain"

# PROP Default_Filter ""
# Begin Source File

SOURCE=..\..\..\SRC\material\uniaxial\drain\DrainBilinearMaterial.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\uniaxial\drain\DrainBilinearMaterial.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\uniaxial\drain\DrainClough1Material.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\uniaxial\drain\DrainClough1Material.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\uniaxial\drain\DrainClough2Material.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\uniaxial\drain\DrainClough2Material.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\uniaxial\drain\DrainHardeningMaterial.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\uniaxial\drain\DrainHardeningMaterial.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\uniaxial\drain\DrainPinch1Material.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\uniaxial\drain\DrainPinch1Material.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\uniaxial\drain\TclDrainMaterialCommand.cpp
# End Source File
# End Group
# Begin Group "py_tz_qz"

# PROP Default_Filter ""
# Begin Source File

SOURCE=..\..\..\Src\material\uniaxial\Py\PySimple1.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\Src\material\uniaxial\Py\PySimple1.h
# End Source File
# Begin Source File

SOURCE=..\..\..\Src\material\uniaxial\Py\QzSimple1.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\Src\material\uniaxial\Py\QzSimple1.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\uniaxial\PY\TclPyTzQzMaterialCommand.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\Src\material\uniaxial\Py\TzSimple1.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\Src\material\uniaxial\Py\TzSimple1.h
# End Source File
# End Group
# Begin Group "snap"

# PROP Default_Filter ""
# Begin Source File

SOURCE=..\..\..\SRC\material\uniaxial\snap\Pinching.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\uniaxial\snap\Pinching.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\uniaxial\snap\TclSnapMaterialCommand.cpp
# End Source File
# End Group
# Begin Source File

SOURCE=..\..\..\SRC\material\uniaxial\BoucWenMaterial.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\uniaxial\BoucWenMaterial.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\uniaxial\CableMaterial.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\uniaxial\CableMaterial.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\uniaxial\Concrete01.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\uniaxial\Concrete01.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\uniaxial\DrainMaterial.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\uniaxial\DrainMaterial.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\uniaxial\ElasticMaterial.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\uniaxial\ElasticMaterial.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\uniaxial\ElasticPPMaterial.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\uniaxial\ElasticPPMaterial.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\uniaxial\ENTMaterial.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\uniaxial\ENTMaterial.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\uniaxial\EPPGapMaterial.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\uniaxial\EPPGapMaterial.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\uniaxial\FedeasMaterial.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\uniaxial\FedeasMaterial.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\uniaxial\HardeningMaterial.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\uniaxial\HardeningMaterial.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\uniaxial\HystereticMaterial.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\uniaxial\HystereticMaterial.h
# End Source File
# Begin Source File

SOURCE=..\..\..\Src\material\uniaxial\MinMaxMaterial.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\Src\material\uniaxial\MinMaxMaterial.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\uniaxial\NewUniaxialMaterial.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\uniaxial\NewUniaxialMaterial.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\uniaxial\ParallelMaterial.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\uniaxial\ParallelMaterial.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\uniaxial\PathIndependentMaterial.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\uniaxial\PathIndependentMaterial.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\uniaxial\SeriesMaterial.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\uniaxial\SeriesMaterial.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\uniaxial\Steel01.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\uniaxial\Steel01.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\uniaxial\TclModelBuilderUniaxialMaterialCommand.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\uniaxial\UniaxialMaterial.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\uniaxial\UniaxialMaterial.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\uniaxial\ViscousMaterial.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\uniaxial\ViscousMaterial.h
# End Source File
# End Group
# Begin Group "nD"

# PROP Default_Filter ""
# Begin Group "elasticIsotropic"

# PROP Default_Filter ""
# Begin Source File

SOURCE=..\..\..\Src\material\nD\ElasticCrossAnisotropic.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\Src\material\nD\ElasticCrossAnisotropic.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\ElasticIsotropic3D.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\ElasticIsotropic3D.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\ElasticIsotropicAxiSymm.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\ElasticIsotropicAxiSymm.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\ElasticIsotropicBeamFiber.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\ElasticIsotropicBeamFiber.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\ElasticIsotropicMaterial.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\ElasticIsotropicMaterial.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\ElasticIsotropicPlaneStrain2D.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\ElasticIsotropicPlaneStrain2D.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\ElasticIsotropicPlaneStress2D.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\ElasticIsotropicPlaneStress2D.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\ElasticIsotropicPlateFiber.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\ElasticIsotropicPlateFiber.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\PressureDependentElastic3D.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\PressureDependentElastic3D.h
# End Source File
# End Group
# Begin Group "j2Plasticity"

# PROP Default_Filter ""
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\J2AxiSymm.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\J2AxiSymm.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\J2PlaneStrain.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\J2PlaneStrain.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\J2PlaneStress.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\J2PlaneStress.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\J2Plasticity.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\J2Plasticity.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\J2PlateFiber.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\J2PlateFiber.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\J2ThreeDimensional.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\J2ThreeDimensional.h
# End Source File
# End Group
# Begin Group "soilModels"

# PROP Default_Filter ""
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\soil\FluidSolidPorousMaterial.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\soil\FluidSolidPorousMaterial.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\soil\MultiYieldSurface.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\soil\MultiYieldSurface.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\soil\PressureDependMultiYield.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\soil\PressureDependMultiYield.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\soil\PressureIndependMultiYield.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\soil\PressureIndependMultiYield.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\soil\T2Vector.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\soil\T2Vector.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\soil\TclUpdateMaterialStageCommand.cpp
# End Source File
# End Group
# Begin Group "template3dep"

# PROP Default_Filter ""
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\Template3Dep\CAM_PS.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\Template3Dep\CAM_PS.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\Template3Dep\CAM_YS.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\Template3Dep\CAM_YS.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\Template3Dep\DP_PS.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\Template3Dep\DP_PS.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\Template3Dep\DP_YS.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\Template3Dep\DP_YS.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\Template3Dep\DP_YS01.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\Template3Dep\DP_YS01.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\Template3Dep\EL_LEeq.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\Template3Dep\EL_LEeq.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\Template3Dep\EL_LEij.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\Template3Dep\EL_LEij.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\Template3Dep\EL_NLEeq.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\Template3Dep\EL_NLEeq.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\Template3Dep\EL_NLEij.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\Template3Dep\EL_NLEij.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\Template3Dep\EL_NLEijMD.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\Template3Dep\EL_NLEijMD.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\Template3Dep\EL_NLEp.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\Template3Dep\EL_NLEp.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\Template3Dep\EL_S.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\Template3Dep\EL_S.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\Template3Dep\EL_T.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\Template3Dep\EL_T.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\Template3Dep\EPState.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\Template3Dep\EPState.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\Template3Dep\MatPoint3D.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\Template3Dep\MatPoint3D.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\Template3Dep\MD_PS.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\Template3Dep\MD_PS.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\Template3Dep\MD_PS01.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\Template3Dep\MD_PS01.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\Template3Dep\MD_YS.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\Template3Dep\MD_YS.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\Template3Dep\RMC01.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\Template3Dep\RMC01.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\Template3Dep\RMC01_PS.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\Template3Dep\RMC01_PS.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\Template3Dep\RMC01_YS.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\Template3Dep\RMC01_YS.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\Template3Dep\TclTemplate3DepCommand.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\Template3Dep\Template3Dep.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\Template3Dep\Template3Dep.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\Template3Dep\VM_PS.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\Template3Dep\VM_PS.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\Template3Dep\VM_YS.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\Template3Dep\VM_YS.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\Template3Dep\YS.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\Template3Dep\YS.h
# End Source File
# End Group
# Begin Group "wrappers"

# PROP Default_Filter ""
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\BeamFiberMaterial.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\BeamFiberMaterial.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\PlaneStressMaterial.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\PlaneStressMaterial.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\PlateFiberMaterial.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\PlateFiberMaterial.h
# End Source File
# End Group
# Begin Group "feap"

# PROP Default_Filter ""
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\FeapMaterial.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\FeapMaterial.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\feap\FeapMaterial01.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\feap\FeapMaterial01.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\feap\FeapMaterial02.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\feap\FeapMaterial02.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\feap\FeapMaterial03.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\feap\FeapMaterial03.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\feap\TclFeapMaterialCommand.cpp
# End Source File
# End Group
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\NDMaterial.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\NDMaterial.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\TclModelBuilderNDMaterialCommand.cpp
# End Source File
# End Group
# Begin Group "section"

# PROP Default_Filter ""
# Begin Group "repres"

# PROP Default_Filter ""
# Begin Group "cell"

# PROP Default_Filter ".cpp; .h"
# Begin Source File

SOURCE=..\..\..\SRC\material\section\repres\cell\Cell.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\section\repres\cell\Cell.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\section\repres\cell\QuadCell.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\section\repres\cell\QuadCell.h
# End Source File
# End Group
# Begin Group "patch"

# PROP Default_Filter ""
# Begin Source File

SOURCE=..\..\..\SRC\material\section\repres\patch\CircPatch.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\section\repres\patch\CircPatch.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\section\repres\patch\Patch.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\section\repres\patch\Patch.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\section\repres\patch\QuadPatch.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\section\repres\patch\QuadPatch.h
# End Source File
# End Group
# Begin Group "reinfBar"

# PROP Default_Filter ".cpp; .h"
# Begin Source File

SOURCE=..\..\..\SRC\material\section\repres\reinfBar\ReinfBar.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\section\repres\reinfBar\ReinfBar.h
# End Source File
# End Group
# Begin Group "reinfLayer"

# PROP Default_Filter ".cpp; .h"
# Begin Source File

SOURCE=..\..\..\SRC\material\section\repres\reinfLayer\CircReinfLayer.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\section\repres\reinfLayer\CircReinfLayer.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\section\repres\reinfLayer\ReinfLayer.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\section\repres\reinfLayer\ReinfLayer.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\section\repres\reinfLayer\StraightReinfLayer.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\section\repres\reinfLayer\StraightReinfLayer.h
# End Source File
# End Group
# Begin Group "sect"

# PROP Default_Filter ".cpp; .h"
# Begin Source File

SOURCE=..\..\..\SRC\material\section\repres\section\FiberSectionRepr.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\section\repres\section\FiberSectionRepr.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\section\repres\section\SectionRepres.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\section\repres\section\SectionRepres.h
# End Source File
# End Group
# End Group
# Begin Group "fiber"

# PROP Default_Filter ""
# Begin Source File

SOURCE=..\..\..\SRC\material\section\fiber\Fiber.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\section\fiber\Fiber.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\section\fiber\UniaxialFiber2d.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\section\fiber\UniaxialFiber2d.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\section\fiber\UniaxialFiber3d.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\section\fiber\UniaxialFiber3d.h
# End Source File
# End Group
# Begin Group "ysSection"

# PROP Default_Filter ""
# Begin Source File

SOURCE=..\..\..\SRC\material\section\yieldSurface\TclModelBuilderYS_SectionCommand.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\section\yieldSurface\YieldSurfaceSection2d.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\section\yieldSurface\YieldSurfaceSection2d.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\section\yieldSurface\YS_Section2D01.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\section\yieldSurface\YS_Section2D01.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\section\yieldSurface\YS_Section2D02.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\section\yieldSurface\YS_Section2D02.h
# End Source File
# End Group
# Begin Source File

SOURCE=..\..\..\Src\material\section\Bidirectional.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\Src\material\section\Bidirectional.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\section\ElasticMembranePlateSection.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\section\ElasticMembranePlateSection.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\section\ElasticPlateSection.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\section\ElasticPlateSection.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\section\ElasticSection2d.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\section\ElasticSection2d.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\section\ElasticSection3d.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\section\ElasticSection3d.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\section\FiberSection2d.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\section\FiberSection2d.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\section\FiberSection3d.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\section\FiberSection3d.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\section\FiberSectionGJ.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\section\FiberSectionGJ.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\section\GenericSection1d.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\section\GenericSection1d.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\section\MembranePlateFiberSection.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\section\MembranePlateFiberSection.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\section\SectionAggregator.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\section\SectionAggregator.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\section\SectionForceDeformation.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\section\SectionForceDeformation.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\section\TclModelBuilderSectionCommand.cpp
# End Source File
# End Group
# Begin Group "yieldSurface"

# PROP Default_Filter ""
# Begin Group "yieldSurfaceBC"

# PROP Default_Filter ""
# Begin Source File

SOURCE=..\..\..\SRC\material\yieldSurface\yieldSurfaceBC\Attalla2D.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\yieldSurface\yieldSurfaceBC\Attalla2D.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\yieldSurface\yieldSurfaceBC\ElTawil2D.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\yieldSurface\yieldSurfaceBC\ElTawil2D.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\yieldSurface\yieldSurfaceBC\ElTawil2DUnSym.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\yieldSurface\yieldSurfaceBC\ElTawil2DUnSym.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\yieldSurface\yieldSurfaceBC\Hajjar2D.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\yieldSurface\yieldSurfaceBC\Hajjar2D.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\yieldSurface\yieldSurfaceBC\NullYS2D.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\yieldSurface\yieldSurfaceBC\NullYS2D.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\yieldSurface\yieldSurfaceBC\Orbison2D.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\yieldSurface\yieldSurfaceBC\Orbison2D.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\yieldSurface\yieldSurfaceBC\TclModelBuilderYieldSurfaceBCCommand.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\yieldSurface\yieldSurfaceBC\YieldSurface_BC.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\yieldSurface\yieldSurfaceBC\YieldSurface_BC.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\yieldSurface\yieldSurfaceBC\YieldSurface_BC2D.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\yieldSurface\yieldSurfaceBC\YieldSurface_BC2D.h
# End Source File
# End Group
# Begin Group "evolution"

# PROP Default_Filter ""
# Begin Source File

SOURCE=..\..\..\SRC\material\yieldSurface\evolution\BkStressLimSurface2D.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\yieldSurface\evolution\BkStressLimSurface2D.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\yieldSurface\evolution\BoundingSurface2D.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\yieldSurface\evolution\BoundingSurface2D.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\yieldSurface\evolution\CombinedIsoKin2D01.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\yieldSurface\evolution\CombinedIsoKin2D01.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\yieldSurface\evolution\CombinedIsoKin2D02.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\yieldSurface\evolution\CombinedIsoKin2D02.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\yieldSurface\evolution\Isotropic2D01.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\yieldSurface\evolution\Isotropic2D01.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\yieldSurface\evolution\Kinematic2D01.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\yieldSurface\evolution\Kinematic2D01.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\yieldSurface\evolution\Kinematic2D02.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\yieldSurface\evolution\Kinematic2D02.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\yieldSurface\evolution\NullEvolution.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\yieldSurface\evolution\NullEvolution.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\yieldSurface\evolution\PeakOriented2D01.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\yieldSurface\evolution\PeakOriented2D01.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\yieldSurface\evolution\PeakOriented2D02.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\yieldSurface\evolution\PeakOriented2D02.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\yieldSurface\evolution\PlasticHardening2D.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\yieldSurface\evolution\PlasticHardening2D.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\yieldSurface\evolution\TclModelBuilderYS_EvolutionCommand.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\yieldSurface\evolution\YS_Evolution.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\yieldSurface\evolution\YS_Evolution.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\yieldSurface\evolution\YS_Evolution2D.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\yieldSurface\evolution\YS_Evolution2D.h
# End Source File
# End Group
# Begin Group "plasticHardening"

# PROP Default_Filter ""
# Begin Source File

SOURCE=..\..\..\SRC\material\yieldSurface\plasticHardeningMaterial\ExponReducing.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\yieldSurface\plasticHardeningMaterial\ExponReducing.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\yieldSurface\plasticHardeningMaterial\MultiLinearKp.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\yieldSurface\plasticHardeningMaterial\MultiLinearKp.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\yieldSurface\plasticHardeningMaterial\NullPlasticMaterial.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\yieldSurface\plasticHardeningMaterial\NullPlasticMaterial.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\yieldSurface\plasticHardeningMaterial\PlasticHardeningMaterial.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\yieldSurface\plasticHardeningMaterial\PlasticHardeningMaterial.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\yieldSurface\plasticHardeningMaterial\TclModelBuilderYSPlasticMaterialCommand.cpp
# End Source File
# End Group
# End Group
# Begin Source File

SOURCE=..\..\..\SRC\material\Material.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\Material.h
# End Source File
# End Target
# End Project
