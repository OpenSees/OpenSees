# Microsoft Developer Studio Project File - Name="openSees" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Console Application" 0x0103

CFG=openSees - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "openSees.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "openSees.mak" CFG="openSees - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "openSees - Win32 Release" (based on "Win32 (x86) Console Application")
!MESSAGE "openSees - Win32 Debug" (based on "Win32 (x86) Console Application")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
RSC=rc.exe

!IF  "$(CFG)" == "openSees - Win32 Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "Release"
# PROP BASE Intermediate_Dir "Release"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "..\..\bin"
# PROP Intermediate_Dir "..\..\obj\openSees\release"
# PROP Ignore_Export_Lib 0
# PROP Target_Dir ""
# ADD BASE CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /c
# ADD CPP /nologo /MLd /W3 /GX /O2 /I "..\..\..\src\damage" /I "..\..\..\src\material\uniaxial\py" /I "..\..\..\src\reliability\fesensitivity" /I "..\..\..\src\reliability\tcl" /I "..\..\..\src\coordTransformation" /I "..\..\..\src\material\section\repres\section" /I "..\..\..\src\analysis\algorithm\equiSolnAlgo" /I "..\..\..\src\system_of_eqn\eigenSOE" /I "..\..\..\src\analysis\algorithm\eigenAlgo" /I "..\..\..\src\material\nD" /I "..\..\..\src\material\uniaxial" /I "..\..\..\src\handler" /I "..\..\..\src\tcl" /I "..\..\..\src\actor\objectBroker" /I "..\..\..\src\system_of_eqn\linearSOE\umfGEN" /I "..\..\..\src\system_of_eqn\linearSOE\fullGEN" /I "..\..\..\src\system_of_eqn\linearSOE\sparseGEN" /I "..\..\..\src\system_of_eqn\linearSOE\bandSPD" /I "..\..\..\src\system_of_eqn\linearSOE\bandGEN" /I "..\..\..\src\element\nonlinearBeamColumn\tcl\repres\section" /I "..\..\..\src\matrix" /I "..\..\..\src\recorder" /I "..\..\..\src\graph\numberer" /I "..\..\..\src\material\section" /I "..\..\..\src\graph\graph" /I "..\..\..\src\element\beam2d" /I "..\..\..\src\element\beam3d" /I "..\..\..\src\system_of_eqn" /I "..\..\..\src\system_of_eqn\linearSOE" /I "..\..\..\src\system_of_eqn\linearSOE\profileSPD" /I "..\..\..\src\system_of_eqn\linearSOE\sparseSYM" /I "..\..\..\src\domain\pattern" /I "..\..\..\src\analysis\analysis" /I "..\..\..\src\analysis\integrator" /I "..\..\..\src\analysis\numberer" /I "..\..\..\src\analysis\handler" /I "..\..\..\src\renderer" /I "..\..\..\src\material" /I "..\..\..\src\analysis\algorithm" /I "..\..\..\src\convergenceTest" /I "..\..\..\src\analysis\model\simple" /I "..\..\..\src\domain\load" /I "..\..\..\src\analysis\model" /I "..\..\..\src\element\truss" /I "..\..\..\src\actor\channel" /I "..\..\..\src\utility" /I "..\..\..\src\actor\actor" /I "..\..\..\src\modelbuilder" /I "..\..\..\src\modelbuilder\tcl" /I "..\..\..\src\domain\constraints" /I "..\..\..\src\domain\component" /I "..\..\..\src\element" /I "..\..\..\src\domain\node" /I "..\..\..\src\domain\domain" /I "..\..\..\src\tagged\storage" /I "..\..\..\src" /I "..\..\..\src\tagged" /I "..\..\..\src\reliability\domain" /I "..\..\..\src\reliability\domain\components" /I "..\..\..\src\reliability\domain\distributions" /I "..\..\..\src\reliability\analysis" /I "..\..\..\src\reliability\analysis\analysis" /I "..\..\..\src\reliability\analysis\curvature" /I "..\..\..\src\reliability\analysis\designPoint" /I "..\..\..\src\reliability\analysis\direction" /I "..\..\..\src\reliability\analysis\gFunction" /I "..\..\..\src\reliability\analysis\misc" /I "..\..\..\src\reliability\analysis\randomNumber" /I "..\..\..\src\reliability\analysis\sensitivity" /I "..\..\..\src\reliability\analysis\stepSize" /I "..\..\..\src\reliability\analysis\transformation" /I "..\..\..\src\nDarray" /I "..\..\..\src\system_of_eqn\linearSOE\cg" /I "..\..\..\src\system_of_eqn\linearSOE\itpack" /I "..\..\..\other\superlu_mt" /I "..\..\..\src\database" /I "..\..\..\src\element\updatedLagrangianBeamColumn" /I "..\..\..\src\material\yieldSurface\yieldSurfaceBC" /I "..\..\..\src\material\yieldSurface\evolution" /I "..\..\..\src\material\yieldSurface\plasticHardeningMaterial" /I "..\..\..\src\reliability\domain\modulatingFunction" /I "..\..\..\src\reliability\domain\spectrum" /I "..\..\..\src\reliability\domain\filter" /I "..\..\..\src\reliability\analysis\hessianApproximation" /I "..\..\..\src\reliability\analysis\convergenceCheck" /I "..\..\..\src\reliability\analysis\meritFunction" /I "..\..\..\src\reliability\analysis\rootFinding" /I "c:\Program Files\tcl" /D "NDEBUG" /D "_WGL" /D "_RELIABILITY" /D "_WIN32" /D "_FORTRAN" /D "WIN32" /D "_CONSOLE" /D "BUILD_tcl" /D "_MBCS" /D "_TCL84" /FD /c
# SUBTRACT CPP /YX
# ADD BASE RSC /l 0x409 /d "NDEBUG"
# ADD RSC /l 0x409 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlopenSees2.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib kernel32.lib user32.lib gdi32.lib winspool.lib comdlopenSees2.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:console /machine:I386
# ADD LINK32 OpenGL32.lib glu32.lib GlAux.lib damage.lib corotational.lib fedeas.lib drain.lib reliability.lib database.lib renderer.lib blas.lib lapack.lib feap.lib arpack.lib umfpack.lib openSeesFortran.lib actor.lib analysis.lib cblas.lib convergence.lib domain.lib element.lib graph.lib material.lib matrix.lib modelbuilder.lib recorder.lib superLU.lib system.lib tagged.lib utility.lib tcl84.lib tk84.lib kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:console /pdb:none /machine:I386 /nodefaultlib:"libc.lib" /libpath:"..\..\lib" /libpath:"..\..\lib\release" /FORCE:MULTIPLE

!ELSEIF  "$(CFG)" == "openSees - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "Debug"
# PROP BASE Intermediate_Dir "Debug"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "..\..\bin"
# PROP Intermediate_Dir "..\..\obj\openSees\debug"
# PROP Ignore_Export_Lib 0
# PROP Target_Dir ""
# ADD BASE CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /GZ /c
# ADD CPP /nologo /W3 /GX /ZI /Od /I "..\..\..\src\damage" /I "..\..\..\src\material\uniaxial\py" /I "..\..\..\src\reliability\fesensitivity" /I "..\..\..\src\reliability\tcl" /I "..\..\..\src\coordTransformation" /I "..\..\..\src\material\section\repres\section" /I "..\..\..\src\analysis\algorithm\equiSolnAlgo" /I "..\..\..\src\system_of_eqn\eigenSOE" /I "..\..\..\src\analysis\algorithm\eigenAlgo" /I "..\..\..\src\material\nD" /I "..\..\..\src\material\uniaxial" /I "..\..\..\src\handler" /I "..\..\..\src\tcl" /I "..\..\..\src\actor\objectBroker" /I "..\..\..\src\system_of_eqn\linearSOE\umfGEN" /I "..\..\..\src\system_of_eqn\linearSOE\fullGEN" /I "..\..\..\src\system_of_eqn\linearSOE\sparseGEN" /I "..\..\..\src\system_of_eqn\linearSOE\bandSPD" /I "..\..\..\src\system_of_eqn\linearSOE\bandGEN" /I "..\..\..\src\element\nonlinearBeamColumn\tcl\repres\section" /I "..\..\..\src\matrix" /I "..\..\..\src\recorder" /I "..\..\..\src\graph\numberer" /I "..\..\..\src\material\section" /I "..\..\..\src\graph\graph" /I "..\..\..\src\element\beam2d" /I "..\..\..\src\element\beam3d" /I "..\..\..\src\system_of_eqn" /I "..\..\..\src\system_of_eqn\linearSOE" /I "..\..\..\src\system_of_eqn\linearSOE\profileSPD" /I "..\..\..\src\system_of_eqn\linearSOE\sparseSYM" /I "..\..\..\src\domain\pattern" /I "..\..\..\src\analysis\analysis" /I "..\..\..\src\analysis\integrator" /I "..\..\..\src\analysis\numberer" /I "..\..\..\src\analysis\handler" /I "..\..\..\src\renderer" /I "..\..\..\src\material" /I "..\..\..\src\analysis\algorithm" /I "..\..\..\src\convergenceTest" /I "..\..\..\src\analysis\model\simple" /I "..\..\..\src\domain\load" /I "..\..\..\src\analysis\model" /I "..\..\..\src\element\truss" /I "..\..\..\src\actor\channel" /I "..\..\..\src\utility" /I "..\..\..\src\actor\actor" /I "..\..\..\src\modelbuilder" /I "..\..\..\src\modelbuilder\tcl" /I "..\..\..\src\domain\constraints" /I "..\..\..\src\domain\component" /I "..\..\..\src\element" /I "..\..\..\src\domain\node" /I "..\..\..\src\domain\domain" /I "..\..\..\src\tagged\storage" /I "..\..\..\src" /I "..\..\..\src\tagged" /I "..\..\..\src\reliability\domain" /I "..\..\..\src\reliability\domain\components" /I "..\..\..\src\reliability\domain\distributions" /I "..\..\..\src\reliability\analysis" /I "..\..\..\src\reliability\analysis\analysis" /I "..\..\..\src\reliability\analysis\curvature" /I "..\..\..\src\reliability\analysis\designPoint" /I "..\..\..\src\reliability\analysis\direction" /I "..\..\..\src\reliability\analysis\gFunction" /I "..\..\..\src\reliability\analysis\misc" /I "..\..\..\src\reliability\analysis\randomNumber" /I "..\..\..\src\reliability\analysis\sensitivity" /I "..\..\..\src\reliability\analysis\stepSize" /I "..\..\..\src\reliability\analysis\transformation" /I "..\..\..\src\nDarray" /I "..\..\..\src\system_of_eqn\linearSOE\cg" /I "..\..\..\src\system_of_eqn\linearSOE\itpack" /I "..\..\..\other\superlu_mt" /I "..\..\..\src\database" /I "..\..\..\src\element\updatedLagrangianBeamColumn" /I "..\..\..\src\material\yieldSurface\yieldSurfaceBC" /I "..\..\..\src\material\yieldSurface\evolution" /I "..\..\..\src\material\yieldSurface\plasticHardeningMaterial" /I "..\..\..\src\reliability\domain\modulatingFunction" /I "..\..\..\src\reliability\domain\spectrum" /I "..\..\..\src\reliability\domain\filter" /I "..\..\..\src\reliability\analysis\hessianApproximation" /I "..\..\..\src\reliability\analysis\convergenceCheck" /I "..\..\..\src\reliability\analysis\meritFunction" /I "..\..\..\src\reliability\analysis\rootFinding" /I "c:\Program Files\tcl" /D "_DEBUG" /D "_WGL" /D "_RELIABILITY" /D "_WIN32" /D "_FORTRAN" /D "WIN32" /D "_CONSOLE" /D "BUILD_tcl" /D "_MBCS" /D "_TCL84" /FR /FD /GZ
# ADD BASE RSC /l 0x409 /d "_DEBUG"
# ADD RSC /l 0x409 /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlopenSees2.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib kernel32.lib user32.lib gdi32.lib winspool.lib comdlopenSees2.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:console /debug /machine:I386 /pdbtype:sept
# ADD LINK32 OpenGL32.lib glu32.lib GlAux.lib damage.lib corotational.lib fedeas.lib drain.lib reliability.lib database.lib renderer.lib blas.lib lapack.lib feap.lib arpack.lib umfpack.lib openSeesFortran.lib actor.lib analysis.lib cblas.lib convergence.lib domain.lib element.lib graph.lib material.lib matrix.lib modelbuilder.lib recorder.lib superLU.lib system.lib tagged.lib utility.lib tcl84.lib tk84.lib kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:console /debug /machine:I386 /nodefaultlib:"libc.lib" /libpath:"..\..\lib" /libpath:"..\..\lib\debug" /FORCE:MULTIPLE
# SUBTRACT LINK32 /profile

!ENDIF 

# Begin Target

# Name "openSees - Win32 Release"
# Name "openSees - Win32 Debug"
# Begin Group "Source Files"

# PROP Default_Filter "cpp;c;cxx;rc;def;r;odl;idl;hpj;bat"
# Begin Source File

SOURCE=..\..\..\SRC\tcl\commands.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\modelbuilder\tcl\myCommands.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\tcl\tclAppInit.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\tcl\TclFeViewer.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\tcl\tclMain.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\modelbuilder\tcl\TclModelBuilder.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\Src\modelbuilder\tcl\TclUniaxialMaterialTester.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\tcl\TclVideoPlayer.cpp
# End Source File
# End Group
# Begin Group "Header Files"

# PROP Default_Filter "h;hpp;hxx;hm;inl"
# Begin Source File

SOURCE=..\..\..\SRC\classTags.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\tcl\commands.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\FiniteDeformation\FiniteDeformationElastic3D.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\OPS_Globals.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\tcl\TclFeViewer.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\modelbuilder\tcl\TclModelBuilder.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\tcl\TclVideoPlayer.h
# End Source File
# End Group
# End Target
# End Project
