# Microsoft Developer Studio Generated NMAKE File, Based on g3.dsp
!IF "$(CFG)" == ""
CFG=g3 - Win32 Debug
!MESSAGE No configuration specified. Defaulting to g3 - Win32 Debug.
!ENDIF 

!IF "$(CFG)" != "g3 - Win32 Release" && "$(CFG)" != "g3 - Win32 Debug"
!MESSAGE Invalid configuration "$(CFG)" specified.
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "g3.mak" CFG="g3 - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "g3 - Win32 Release" (based on "Win32 (x86) Console Application")
!MESSAGE "g3 - Win32 Debug" (based on "Win32 (x86) Console Application")
!MESSAGE 
!ERROR An invalid configuration is specified.
!ENDIF 

!IF "$(OS)" == "Windows_NT"
NULL=
!ELSE 
NULL=nul
!ENDIF 

!IF  "$(CFG)" == "g3 - Win32 Release"

OUTDIR=.\Release
INTDIR=.\Release
# Begin Custom Macros
OutDir=.\Release
# End Custom Macros

!IF "$(RECURSE)" == "0" 

ALL : "$(OUTDIR)\g3.exe"

!ELSE 

ALL : "utility - Win32 Release" "tagged - Win32 Release" "system - Win32 Release" "superLU - Win32 Release" "recorder - Win32 Release" "nonlinearBeamColumn - Win32 Release" "modelbuilder - Win32 Release" "matrix - Win32 Release" "material - Win32 Release" "graph - Win32 Release" "element - Win32 Release" "domain - Win32 Release" "convergence - Win32 Release" "cblas - Win32 Release" "analysis - Win32 Release" "actor - Win32 Release" "$(OUTDIR)\g3.exe"

!ENDIF 

!IF "$(RECURSE)" == "1" 
CLEAN :"actor - Win32 ReleaseCLEAN" "analysis - Win32 ReleaseCLEAN" "cblas - Win32 ReleaseCLEAN" "convergence - Win32 ReleaseCLEAN" "domain - Win32 ReleaseCLEAN" "element - Win32 ReleaseCLEAN" "graph - Win32 ReleaseCLEAN" "material - Win32 ReleaseCLEAN" "matrix - Win32 ReleaseCLEAN" "modelbuilder - Win32 ReleaseCLEAN" "nonlinearBeamColumn - Win32 ReleaseCLEAN" "recorder - Win32 ReleaseCLEAN" "superLU - Win32 ReleaseCLEAN" "system - Win32 ReleaseCLEAN" "tagged - Win32 ReleaseCLEAN" "utility - Win32 ReleaseCLEAN" 
!ELSE 
CLEAN :
!ENDIF 
	-@erase "$(INTDIR)\commands.obj"
	-@erase "$(INTDIR)\myCommands.obj"
	-@erase "$(INTDIR)\tclAppInit.obj"
	-@erase "$(INTDIR)\tclMain.obj"
	-@erase "$(INTDIR)\TclModelBuilder.obj"
	-@erase "$(INTDIR)\vc60.idb"
	-@erase "$(OUTDIR)\g3.exe"

"$(OUTDIR)" :
    if not exist "$(OUTDIR)/$(NULL)" mkdir "$(OUTDIR)"

CPP=cl.exe
CPP_PROJ=/nologo /ML /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_CONSOLE" /D "_MBCS" /Fp"$(INTDIR)\g3.pch" /YX /Fo"$(INTDIR)\\" /Fd"$(INTDIR)\\" /FD /c 

.c{$(INTDIR)}.obj::
   $(CPP) @<<
   $(CPP_PROJ) $< 
<<

.cpp{$(INTDIR)}.obj::
   $(CPP) @<<
   $(CPP_PROJ) $< 
<<

.cxx{$(INTDIR)}.obj::
   $(CPP) @<<
   $(CPP_PROJ) $< 
<<

.c{$(INTDIR)}.sbr::
   $(CPP) @<<
   $(CPP_PROJ) $< 
<<

.cpp{$(INTDIR)}.sbr::
   $(CPP) @<<
   $(CPP_PROJ) $< 
<<

.cxx{$(INTDIR)}.sbr::
   $(CPP) @<<
   $(CPP_PROJ) $< 
<<

RSC=rc.exe
BSC32=bscmake.exe
BSC32_FLAGS=/nologo /o"$(OUTDIR)\g3.bsc" 
BSC32_SBRS= \
	
LINK32=link.exe
LINK32_FLAGS=kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:console /incremental:no /pdb:"$(OUTDIR)\g3.pdb" /machine:I386 /out:"$(OUTDIR)\g3.exe" 
LINK32_OBJS= \
	"$(INTDIR)\commands.obj" \
	"$(INTDIR)\myCommands.obj" \
	"$(INTDIR)\tclAppInit.obj" \
	"$(INTDIR)\tclMain.obj" \
	"$(INTDIR)\TclModelBuilder.obj" \
	"..\actor\Release\actor.lib" \
	"..\analysis\Release\analysis.lib" \
	"..\cblas\Release\cblas.lib" \
	"..\convergence\Release\convergence.lib" \
	"..\domain\Release\domain.lib" \
	"..\element\Release\element.lib" \
	"..\graph\Release\graph.lib" \
	"..\material\Release\material.lib" \
	"..\matrix\Release\matrix.lib" \
	"..\modelbuilder\Release\modelbuilder.lib" \
	"..\nonlinearBeamColumn\Release\nonlinearBeamColumn.lib" \
	"..\recorder\Release\recorder.lib" \
	"..\superLU\Release\superLU.lib" \
	"..\system\Release\system.lib" \
	"..\tagged\Release\tagged.lib" \
	"..\utility\Release\utility.lib"

"$(OUTDIR)\g3.exe" : "$(OUTDIR)" $(DEF_FILE) $(LINK32_OBJS)
    $(LINK32) @<<
  $(LINK32_FLAGS) $(LINK32_OBJS)
<<

!ELSEIF  "$(CFG)" == "g3 - Win32 Debug"

OUTDIR=.\..\..\bin
INTDIR=.\..\..\obj\g3
# Begin Custom Macros
OutDir=.\..\..\bin
# End Custom Macros

!IF "$(RECURSE)" == "0" 

ALL : "$(OUTDIR)\g3.exe"

!ELSE 

ALL : "utility - Win32 Debug" "tagged - Win32 Debug" "system - Win32 Debug" "superLU - Win32 Debug" "recorder - Win32 Debug" "nonlinearBeamColumn - Win32 Debug" "modelbuilder - Win32 Debug" "matrix - Win32 Debug" "material - Win32 Debug" "graph - Win32 Debug" "element - Win32 Debug" "domain - Win32 Debug" "convergence - Win32 Debug" "cblas - Win32 Debug" "analysis - Win32 Debug" "actor - Win32 Debug" "$(OUTDIR)\g3.exe"

!ENDIF 

!IF "$(RECURSE)" == "1" 
CLEAN :"actor - Win32 DebugCLEAN" "analysis - Win32 DebugCLEAN" "cblas - Win32 DebugCLEAN" "convergence - Win32 DebugCLEAN" "domain - Win32 DebugCLEAN" "element - Win32 DebugCLEAN" "graph - Win32 DebugCLEAN" "material - Win32 DebugCLEAN" "matrix - Win32 DebugCLEAN" "modelbuilder - Win32 DebugCLEAN" "nonlinearBeamColumn - Win32 DebugCLEAN" "recorder - Win32 DebugCLEAN" "superLU - Win32 DebugCLEAN" "system - Win32 DebugCLEAN" "tagged - Win32 DebugCLEAN" "utility - Win32 DebugCLEAN" 
!ELSE 
CLEAN :
!ENDIF 
	-@erase "$(INTDIR)\commands.obj"
	-@erase "$(INTDIR)\myCommands.obj"
	-@erase "$(INTDIR)\tclAppInit.obj"
	-@erase "$(INTDIR)\tclMain.obj"
	-@erase "$(INTDIR)\TclModelBuilder.obj"
	-@erase "$(INTDIR)\vc60.idb"
	-@erase "$(OUTDIR)\g3.exe"

"$(OUTDIR)" :
    if not exist "$(OUTDIR)/$(NULL)" mkdir "$(OUTDIR)"

"$(INTDIR)" :
    if not exist "$(INTDIR)/$(NULL)" mkdir "$(INTDIR)"

CPP=cl.exe
CPP_PROJ=/nologo /MLd /W3 /GX /Od /I "D:\Program Files\tcl\include" /I "..\..\src\tcl" /I "..\..\src\actor\objectBroker" /I "..\..\src\system_of_eqn\linearSOE\umfGEN" /I "..\..\src\system_of_eqn\linearSOE\fullGEN" /I "..\..\src\system_of_eqn\linearSOE\sparseGEN" /I "..\..\src\system_of_eqn\linearSOE\bandSPD" /I "..\..\src\system_of_eqn\linearSOE\bandGEN" /I "..\..\src\element\nonlinearBeamColumn\tcl\repres\section" /I "..\..\src\matrix" /I "..\..\src\recorder" /I "..\..\src\graph\numberer" /I "..\..\src\element\nonlinearBeamColumn\section" /I "..\..\src\graph\graph" /I "..\..\src\element\beam2d" /I "..\..\src\element\beam3d" /I "..\..\src\system_of_eqn" /I "..\..\src\system_of_eqn\linearSOE" /I "..\..\src\system_of_eqn\linearSOE\profileSPD" /I "..\..\src\system_of_eqn\linearSOE\sparseSYM" /I "..\..\src\domain\pattern" /I "..\..\src\analysis\analysis" /I "..\..\src\analysis\integrator" /I "..\..\src\analysis\numberer" /I "..\..\src\analysis\handler" /I "..\..\src\renderer" /I "..\..\src\material" /I "..\..\src\analysis\algorithm\equiSolnAlgo" /I "..\..\src\analysis\algorithm" /I "..\..\src\convergenceTest" /I "..\..\src\analysis\model\simple" /I "..\..\src\domain\load" /I\
 "..\..\src\analysis\model" /I "..\..\src\element\truss" /I "..\..\src\actor\channel" /I "..\..\src\utility" /I "..\..\src\actor\actor" /I "..\..\src\modelbuilder" /I "..\..\src\modelbuilder\tcl" /I "..\..\src\domain\constraints" /I "..\..\src\domain\component" /I "..\..\src\element" /I "..\..\src\domain\node" /I "..\..\src\domain\domain" /I "..\..\src\tagged\storage" /I "..\..\src" /I "..\..\src\tagged" /D "_FORTRAN" /D "WIN32" /D "_DEBUG" /D "_CONSOLE" /D "_MBCS" /D "BUILD_tcl" /Fp"$(INTDIR)\g3.pch" /YX /Fo"$(INTDIR)\\" /Fd"$(INTDIR)\\" /FD /I /GZ Files\Tcl\include "   

.c{$(INTDIR)}.obj::
   $(CPP) @<<
   $(CPP_PROJ) $< 
<<

.cpp{$(INTDIR)}.obj::
   $(CPP) @<<
   $(CPP_PROJ) $< 
<<

.cxx{$(INTDIR)}.obj::
   $(CPP) @<<
   $(CPP_PROJ) $< 
<<

.c{$(INTDIR)}.sbr::
   $(CPP) @<<
   $(CPP_PROJ) $< 
<<

.cpp{$(INTDIR)}.sbr::
   $(CPP) @<<
   $(CPP_PROJ) $< 
<<

.cxx{$(INTDIR)}.sbr::
   $(CPP) @<<
   $(CPP_PROJ) $< 
<<

RSC=rc.exe
BSC32=bscmake.exe
BSC32_FLAGS=/nologo /o"$(OUTDIR)\g3.bsc" 
BSC32_SBRS= \
	
LINK32=link.exe
LINK32_FLAGS=g3Fortran.lib actor.lib analysis.lib cblas.lib convergence.lib domain.lib element.lib graph.lib material.lib matrix.lib modelbuilder.lib nonlinearBeamColumn.lib recorder.lib superLU.lib system.lib tagged.lib utility.lib tcl82.lib kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:console /incremental:no /pdb:"$(OUTDIR)\g3.pdb" /machine:I386 /nodefaultlib:"libc.lib" /out:"$(OUTDIR)\g3.exe" /pdbtype:sept /libpath:"f:\Program Files\tcl\lib" /libpath:"c:\msdev\lib" /libpath:"..\..\lib" /FORCE:MULTIPLE 
LINK32_OBJS= \
	"$(INTDIR)\commands.obj" \
	"$(INTDIR)\myCommands.obj" \
	"$(INTDIR)\tclAppInit.obj" \
	"$(INTDIR)\tclMain.obj" \
	"$(INTDIR)\TclModelBuilder.obj" \
	"..\..\lib\actor.lib" \
	"..\..\lib\analysis.lib" \
	"..\..\lib\cblas.lib" \
	"..\..\lib\convergence.lib" \
	"..\..\lib\domain.lib" \
	"..\..\lib\element.lib" \
	"..\..\lib\graph.lib" \
	"..\..\lib\material.lib" \
	"..\..\lib\matrix.lib" \
	"..\..\lib\modelbuilder.lib" \
	"..\..\lib\nonlinearBeamColumn.lib" \
	"..\..\lib\recorder.lib" \
	"..\..\lib\superLU.lib" \
	"..\..\lib\system.lib" \
	"..\..\lib\tagged.lib" \
	"..\..\lib\utility.lib"

"$(OUTDIR)\g3.exe" : "$(OUTDIR)" $(DEF_FILE) $(LINK32_OBJS)
    $(LINK32) @<<
  $(LINK32_FLAGS) $(LINK32_OBJS)
<<

!ENDIF 


!IF "$(NO_EXTERNAL_DEPS)" != "1"
!IF EXISTS("g3.dep")
!INCLUDE "g3.dep"
!ELSE 
!MESSAGE Warning: cannot find "g3.dep"
!ENDIF 
!ENDIF 


!IF "$(CFG)" == "g3 - Win32 Release" || "$(CFG)" == "g3 - Win32 Debug"
SOURCE=..\..\Src\tcl\commands.cpp

"$(INTDIR)\commands.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\Src\modelbuilder\tcl\myCommands.cpp

"$(INTDIR)\myCommands.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\Src\tcl\tclAppInit.cpp

"$(INTDIR)\tclAppInit.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\Src\tcl\tclMain.cpp

"$(INTDIR)\tclMain.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\Src\modelbuilder\tcl\TclModelBuilder.cpp

"$(INTDIR)\TclModelBuilder.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!IF  "$(CFG)" == "g3 - Win32 Release"

"actor - Win32 Release" : 
   cd "\g3\proj\actor"
   $(MAKE) /$(MAKEFLAGS) /F .\actor.mak CFG="actor - Win32 Release" 
   cd "..\g3"

"actor - Win32 ReleaseCLEAN" : 
   cd "\g3\proj\actor"
   $(MAKE) /$(MAKEFLAGS) /F .\actor.mak CFG="actor - Win32 Release" RECURSE=1 CLEAN 
   cd "..\g3"

!ELSEIF  "$(CFG)" == "g3 - Win32 Debug"

"actor - Win32 Debug" : 
   cd "\g3\proj\actor"
   $(MAKE) /$(MAKEFLAGS) /F .\actor.mak CFG="actor - Win32 Debug" 
   cd "..\g3"

"actor - Win32 DebugCLEAN" : 
   cd "\g3\proj\actor"
   $(MAKE) /$(MAKEFLAGS) /F .\actor.mak CFG="actor - Win32 Debug" RECURSE=1 CLEAN 
   cd "..\g3"

!ENDIF 

!IF  "$(CFG)" == "g3 - Win32 Release"

"analysis - Win32 Release" : 
   cd "\g3\proj\analysis"
   $(MAKE) /$(MAKEFLAGS) /F .\analysis.mak CFG="analysis - Win32 Release" 
   cd "..\g3"

"analysis - Win32 ReleaseCLEAN" : 
   cd "\g3\proj\analysis"
   $(MAKE) /$(MAKEFLAGS) /F .\analysis.mak CFG="analysis - Win32 Release" RECURSE=1 CLEAN 
   cd "..\g3"

!ELSEIF  "$(CFG)" == "g3 - Win32 Debug"

"analysis - Win32 Debug" : 
   cd "\g3\proj\analysis"
   $(MAKE) /$(MAKEFLAGS) /F .\analysis.mak CFG="analysis - Win32 Debug" 
   cd "..\g3"

"analysis - Win32 DebugCLEAN" : 
   cd "\g3\proj\analysis"
   $(MAKE) /$(MAKEFLAGS) /F .\analysis.mak CFG="analysis - Win32 Debug" RECURSE=1 CLEAN 
   cd "..\g3"

!ENDIF 

!IF  "$(CFG)" == "g3 - Win32 Release"

"cblas - Win32 Release" : 
   cd "\g3\proj\cblas"
   $(MAKE) /$(MAKEFLAGS) /F .\cblas.mak CFG="cblas - Win32 Release" 
   cd "..\g3"

"cblas - Win32 ReleaseCLEAN" : 
   cd "\g3\proj\cblas"
   $(MAKE) /$(MAKEFLAGS) /F .\cblas.mak CFG="cblas - Win32 Release" RECURSE=1 CLEAN 
   cd "..\g3"

!ELSEIF  "$(CFG)" == "g3 - Win32 Debug"

"cblas - Win32 Debug" : 
   cd "\g3\proj\cblas"
   $(MAKE) /$(MAKEFLAGS) /F .\cblas.mak CFG="cblas - Win32 Debug" 
   cd "..\g3"

"cblas - Win32 DebugCLEAN" : 
   cd "\g3\proj\cblas"
   $(MAKE) /$(MAKEFLAGS) /F .\cblas.mak CFG="cblas - Win32 Debug" RECURSE=1 CLEAN 
   cd "..\g3"

!ENDIF 

!IF  "$(CFG)" == "g3 - Win32 Release"

"convergence - Win32 Release" : 
   cd "\g3\proj\convergence"
   $(MAKE) /$(MAKEFLAGS) /F .\convergence.mak CFG="convergence - Win32 Release" 
   cd "..\g3"

"convergence - Win32 ReleaseCLEAN" : 
   cd "\g3\proj\convergence"
   $(MAKE) /$(MAKEFLAGS) /F .\convergence.mak CFG="convergence - Win32 Release" RECURSE=1 CLEAN 
   cd "..\g3"

!ELSEIF  "$(CFG)" == "g3 - Win32 Debug"

"convergence - Win32 Debug" : 
   cd "\g3\proj\convergence"
   $(MAKE) /$(MAKEFLAGS) /F .\convergence.mak CFG="convergence - Win32 Debug" 
   cd "..\g3"

"convergence - Win32 DebugCLEAN" : 
   cd "\g3\proj\convergence"
   $(MAKE) /$(MAKEFLAGS) /F .\convergence.mak CFG="convergence - Win32 Debug" RECURSE=1 CLEAN 
   cd "..\g3"

!ENDIF 

!IF  "$(CFG)" == "g3 - Win32 Release"

"domain - Win32 Release" : 
   cd "\g3\proj\domain"
   $(MAKE) /$(MAKEFLAGS) /F .\domain.mak CFG="domain - Win32 Release" 
   cd "..\g3"

"domain - Win32 ReleaseCLEAN" : 
   cd "\g3\proj\domain"
   $(MAKE) /$(MAKEFLAGS) /F .\domain.mak CFG="domain - Win32 Release" RECURSE=1 CLEAN 
   cd "..\g3"

!ELSEIF  "$(CFG)" == "g3 - Win32 Debug"

"domain - Win32 Debug" : 
   cd "\g3\proj\domain"
   $(MAKE) /$(MAKEFLAGS) /F .\domain.mak CFG="domain - Win32 Debug" 
   cd "..\g3"

"domain - Win32 DebugCLEAN" : 
   cd "\g3\proj\domain"
   $(MAKE) /$(MAKEFLAGS) /F .\domain.mak CFG="domain - Win32 Debug" RECURSE=1 CLEAN 
   cd "..\g3"

!ENDIF 

!IF  "$(CFG)" == "g3 - Win32 Release"

"element - Win32 Release" : 
   cd "\g3\proj\element"
   $(MAKE) /$(MAKEFLAGS) /F .\element.mak CFG="element - Win32 Release" 
   cd "..\g3"

"element - Win32 ReleaseCLEAN" : 
   cd "\g3\proj\element"
   $(MAKE) /$(MAKEFLAGS) /F .\element.mak CFG="element - Win32 Release" RECURSE=1 CLEAN 
   cd "..\g3"

!ELSEIF  "$(CFG)" == "g3 - Win32 Debug"

"element - Win32 Debug" : 
   cd "\g3\proj\element"
   $(MAKE) /$(MAKEFLAGS) /F .\element.mak CFG="element - Win32 Debug" 
   cd "..\g3"

"element - Win32 DebugCLEAN" : 
   cd "\g3\proj\element"
   $(MAKE) /$(MAKEFLAGS) /F .\element.mak CFG="element - Win32 Debug" RECURSE=1 CLEAN 
   cd "..\g3"

!ENDIF 

!IF  "$(CFG)" == "g3 - Win32 Release"

"graph - Win32 Release" : 
   cd "\g3\proj\graph"
   $(MAKE) /$(MAKEFLAGS) /F .\graph.mak CFG="graph - Win32 Release" 
   cd "..\g3"

"graph - Win32 ReleaseCLEAN" : 
   cd "\g3\proj\graph"
   $(MAKE) /$(MAKEFLAGS) /F .\graph.mak CFG="graph - Win32 Release" RECURSE=1 CLEAN 
   cd "..\g3"

!ELSEIF  "$(CFG)" == "g3 - Win32 Debug"

"graph - Win32 Debug" : 
   cd "\g3\proj\graph"
   $(MAKE) /$(MAKEFLAGS) /F .\graph.mak CFG="graph - Win32 Debug" 
   cd "..\g3"

"graph - Win32 DebugCLEAN" : 
   cd "\g3\proj\graph"
   $(MAKE) /$(MAKEFLAGS) /F .\graph.mak CFG="graph - Win32 Debug" RECURSE=1 CLEAN 
   cd "..\g3"

!ENDIF 

!IF  "$(CFG)" == "g3 - Win32 Release"

"material - Win32 Release" : 
   cd "\g3\proj\material"
   $(MAKE) /$(MAKEFLAGS) /F .\material.mak CFG="material - Win32 Release" 
   cd "..\g3"

"material - Win32 ReleaseCLEAN" : 
   cd "\g3\proj\material"
   $(MAKE) /$(MAKEFLAGS) /F .\material.mak CFG="material - Win32 Release" RECURSE=1 CLEAN 
   cd "..\g3"

!ELSEIF  "$(CFG)" == "g3 - Win32 Debug"

"material - Win32 Debug" : 
   cd "\g3\proj\material"
   $(MAKE) /$(MAKEFLAGS) /F .\material.mak CFG="material - Win32 Debug" 
   cd "..\g3"

"material - Win32 DebugCLEAN" : 
   cd "\g3\proj\material"
   $(MAKE) /$(MAKEFLAGS) /F .\material.mak CFG="material - Win32 Debug" RECURSE=1 CLEAN 
   cd "..\g3"

!ENDIF 

!IF  "$(CFG)" == "g3 - Win32 Release"

"matrix - Win32 Release" : 
   cd "\g3\proj\matrix"
   $(MAKE) /$(MAKEFLAGS) /F .\matrix.mak CFG="matrix - Win32 Release" 
   cd "..\g3"

"matrix - Win32 ReleaseCLEAN" : 
   cd "\g3\proj\matrix"
   $(MAKE) /$(MAKEFLAGS) /F .\matrix.mak CFG="matrix - Win32 Release" RECURSE=1 CLEAN 
   cd "..\g3"

!ELSEIF  "$(CFG)" == "g3 - Win32 Debug"

"matrix - Win32 Debug" : 
   cd "\g3\proj\matrix"
   $(MAKE) /$(MAKEFLAGS) /F .\matrix.mak CFG="matrix - Win32 Debug" 
   cd "..\g3"

"matrix - Win32 DebugCLEAN" : 
   cd "\g3\proj\matrix"
   $(MAKE) /$(MAKEFLAGS) /F .\matrix.mak CFG="matrix - Win32 Debug" RECURSE=1 CLEAN 
   cd "..\g3"

!ENDIF 

!IF  "$(CFG)" == "g3 - Win32 Release"

"modelbuilder - Win32 Release" : 
   cd "\g3\proj\modelbuilder"
   $(MAKE) /$(MAKEFLAGS) /F .\modelbuilder.mak CFG="modelbuilder - Win32 Release" 
   cd "..\g3"

"modelbuilder - Win32 ReleaseCLEAN" : 
   cd "\g3\proj\modelbuilder"
   $(MAKE) /$(MAKEFLAGS) /F .\modelbuilder.mak CFG="modelbuilder - Win32 Release" RECURSE=1 CLEAN 
   cd "..\g3"

!ELSEIF  "$(CFG)" == "g3 - Win32 Debug"

"modelbuilder - Win32 Debug" : 
   cd "\g3\proj\modelbuilder"
   $(MAKE) /$(MAKEFLAGS) /F .\modelbuilder.mak CFG="modelbuilder - Win32 Debug" 
   cd "..\g3"

"modelbuilder - Win32 DebugCLEAN" : 
   cd "\g3\proj\modelbuilder"
   $(MAKE) /$(MAKEFLAGS) /F .\modelbuilder.mak CFG="modelbuilder - Win32 Debug" RECURSE=1 CLEAN 
   cd "..\g3"

!ENDIF 

!IF  "$(CFG)" == "g3 - Win32 Release"

"nonlinearBeamColumn - Win32 Release" : 
   cd "\g3\proj\nonlinearBeamColumn"
   $(MAKE) /$(MAKEFLAGS) /F .\nonlinearBeamColumn.mak CFG="nonlinearBeamColumn - Win32 Release" 
   cd "..\g3"

"nonlinearBeamColumn - Win32 ReleaseCLEAN" : 
   cd "\g3\proj\nonlinearBeamColumn"
   $(MAKE) /$(MAKEFLAGS) /F .\nonlinearBeamColumn.mak CFG="nonlinearBeamColumn - Win32 Release" RECURSE=1 CLEAN 
   cd "..\g3"

!ELSEIF  "$(CFG)" == "g3 - Win32 Debug"

"nonlinearBeamColumn - Win32 Debug" : 
   cd "\g3\proj\nonlinearBeamColumn"
   $(MAKE) /$(MAKEFLAGS) /F .\nonlinearBeamColumn.mak CFG="nonlinearBeamColumn - Win32 Debug" 
   cd "..\g3"

"nonlinearBeamColumn - Win32 DebugCLEAN" : 
   cd "\g3\proj\nonlinearBeamColumn"
   $(MAKE) /$(MAKEFLAGS) /F .\nonlinearBeamColumn.mak CFG="nonlinearBeamColumn - Win32 Debug" RECURSE=1 CLEAN 
   cd "..\g3"

!ENDIF 

!IF  "$(CFG)" == "g3 - Win32 Release"

"recorder - Win32 Release" : 
   cd "\g3\proj\recorder"
   $(MAKE) /$(MAKEFLAGS) /F .\recorder.mak CFG="recorder - Win32 Release" 
   cd "..\g3"

"recorder - Win32 ReleaseCLEAN" : 
   cd "\g3\proj\recorder"
   $(MAKE) /$(MAKEFLAGS) /F .\recorder.mak CFG="recorder - Win32 Release" RECURSE=1 CLEAN 
   cd "..\g3"

!ELSEIF  "$(CFG)" == "g3 - Win32 Debug"

"recorder - Win32 Debug" : 
   cd "\g3\proj\recorder"
   $(MAKE) /$(MAKEFLAGS) /F .\recorder.mak CFG="recorder - Win32 Debug" 
   cd "..\g3"

"recorder - Win32 DebugCLEAN" : 
   cd "\g3\proj\recorder"
   $(MAKE) /$(MAKEFLAGS) /F .\recorder.mak CFG="recorder - Win32 Debug" RECURSE=1 CLEAN 
   cd "..\g3"

!ENDIF 

!IF  "$(CFG)" == "g3 - Win32 Release"

"superLU - Win32 Release" : 
   cd "\g3\proj\superLU"
   $(MAKE) /$(MAKEFLAGS) /F .\superLU.mak CFG="superLU - Win32 Release" 
   cd "..\g3"

"superLU - Win32 ReleaseCLEAN" : 
   cd "\g3\proj\superLU"
   $(MAKE) /$(MAKEFLAGS) /F .\superLU.mak CFG="superLU - Win32 Release" RECURSE=1 CLEAN 
   cd "..\g3"

!ELSEIF  "$(CFG)" == "g3 - Win32 Debug"

"superLU - Win32 Debug" : 
   cd "\g3\proj\superLU"
   $(MAKE) /$(MAKEFLAGS) /F .\superLU.mak CFG="superLU - Win32 Debug" 
   cd "..\g3"

"superLU - Win32 DebugCLEAN" : 
   cd "\g3\proj\superLU"
   $(MAKE) /$(MAKEFLAGS) /F .\superLU.mak CFG="superLU - Win32 Debug" RECURSE=1 CLEAN 
   cd "..\g3"

!ENDIF 

!IF  "$(CFG)" == "g3 - Win32 Release"

"system - Win32 Release" : 
   cd "\g3\proj\system"
   $(MAKE) /$(MAKEFLAGS) /F .\system.mak CFG="system - Win32 Release" 
   cd "..\g3"

"system - Win32 ReleaseCLEAN" : 
   cd "\g3\proj\system"
   $(MAKE) /$(MAKEFLAGS) /F .\system.mak CFG="system - Win32 Release" RECURSE=1 CLEAN 
   cd "..\g3"

!ELSEIF  "$(CFG)" == "g3 - Win32 Debug"

"system - Win32 Debug" : 
   cd "\g3\proj\system"
   $(MAKE) /$(MAKEFLAGS) /F .\system.mak CFG="system - Win32 Debug" 
   cd "..\g3"

"system - Win32 DebugCLEAN" : 
   cd "\g3\proj\system"
   $(MAKE) /$(MAKEFLAGS) /F .\system.mak CFG="system - Win32 Debug" RECURSE=1 CLEAN 
   cd "..\g3"

!ENDIF 

!IF  "$(CFG)" == "g3 - Win32 Release"

"tagged - Win32 Release" : 
   cd "\g3\proj\tagged"
   $(MAKE) /$(MAKEFLAGS) /F .\tagged.mak CFG="tagged - Win32 Release" 
   cd "..\g3"

"tagged - Win32 ReleaseCLEAN" : 
   cd "\g3\proj\tagged"
   $(MAKE) /$(MAKEFLAGS) /F .\tagged.mak CFG="tagged - Win32 Release" RECURSE=1 CLEAN 
   cd "..\g3"

!ELSEIF  "$(CFG)" == "g3 - Win32 Debug"

"tagged - Win32 Debug" : 
   cd "\g3\proj\tagged"
   $(MAKE) /$(MAKEFLAGS) /F .\tagged.mak CFG="tagged - Win32 Debug" 
   cd "..\g3"

"tagged - Win32 DebugCLEAN" : 
   cd "\g3\proj\tagged"
   $(MAKE) /$(MAKEFLAGS) /F .\tagged.mak CFG="tagged - Win32 Debug" RECURSE=1 CLEAN 
   cd "..\g3"

!ENDIF 

!IF  "$(CFG)" == "g3 - Win32 Release"

"utility - Win32 Release" : 
   cd "\g3\proj\utility"
   $(MAKE) /$(MAKEFLAGS) /F .\utility.mak CFG="utility - Win32 Release" 
   cd "..\g3"

"utility - Win32 ReleaseCLEAN" : 
   cd "\g3\proj\utility"
   $(MAKE) /$(MAKEFLAGS) /F .\utility.mak CFG="utility - Win32 Release" RECURSE=1 CLEAN 
   cd "..\g3"

!ELSEIF  "$(CFG)" == "g3 - Win32 Debug"

"utility - Win32 Debug" : 
   cd "\g3\proj\utility"
   $(MAKE) /$(MAKEFLAGS) /F .\utility.mak CFG="utility - Win32 Debug" 
   cd "..\g3"

"utility - Win32 DebugCLEAN" : 
   cd "\g3\proj\utility"
   $(MAKE) /$(MAKEFLAGS) /F .\utility.mak CFG="utility - Win32 Debug" RECURSE=1 CLEAN 
   cd "..\g3"

!ENDIF 


!ENDIF 

