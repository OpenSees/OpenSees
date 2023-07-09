# Microsoft Developer Studio Generated NMAKE File, Based on quickMain.dsp
!IF "$(CFG)" == ""
CFG=quickMain - Win32 Debug
!MESSAGE No configuration specified. Defaulting to quickMain - Win32 Debug.
!ENDIF 

!IF "$(CFG)" != "quickMain - Win32 Release" && "$(CFG)" != "quickMain - Win32 Debug"
!MESSAGE Invalid configuration "$(CFG)" specified.
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "quickMain.mak" CFG="quickMain - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "quickMain - Win32 Release" (based on "Win32 (x86) Console Application")
!MESSAGE "quickMain - Win32 Debug" (based on "Win32 (x86) Console Application")
!MESSAGE 
!ERROR An invalid configuration is specified.
!ENDIF 

!IF "$(OS)" == "Windows_NT"
NULL=
!ELSE 
NULL=nul
!ENDIF 

!IF  "$(CFG)" == "quickMain - Win32 Release"

OUTDIR=.\Release
INTDIR=.\Release
# Begin Custom Macros
OutDir=.\Release
# End Custom Macros

!IF "$(RECURSE)" == "0" 

ALL : "$(OUTDIR)\quickMain.exe"

!ELSE 

ALL : "utility - Win32 Release" "tagged - Win32 Release" "system - Win32 Release" "superLU - Win32 Release" "recorder - Win32 Release" "nonlinearBeamColumn - Win32 Release" "modelbuilder - Win32 Release" "matrix - Win32 Release" "material - Win32 Release" "graph - Win32 Release" "element - Win32 Release" "domain - Win32 Release" "convergence - Win32 Release" "cblas - Win32 Release" "analysis - Win32 Release" "actor - Win32 Release" "$(OUTDIR)\quickMain.exe"

!ENDIF 

!IF "$(RECURSE)" == "1" 
CLEAN :"actor - Win32 ReleaseCLEAN" "analysis - Win32 ReleaseCLEAN" "cblas - Win32 ReleaseCLEAN" "convergence - Win32 ReleaseCLEAN" "domain - Win32 ReleaseCLEAN" "element - Win32 ReleaseCLEAN" "graph - Win32 ReleaseCLEAN" "material - Win32 ReleaseCLEAN" "matrix - Win32 ReleaseCLEAN" "modelbuilder - Win32 ReleaseCLEAN" "nonlinearBeamColumn - Win32 ReleaseCLEAN" "recorder - Win32 ReleaseCLEAN" "superLU - Win32 ReleaseCLEAN" "system - Win32 ReleaseCLEAN" "tagged - Win32 ReleaseCLEAN" "utility - Win32 ReleaseCLEAN" 
!ELSE 
CLEAN :
!ENDIF 
	-@erase "$(INTDIR)\myMain.obj"
	-@erase "$(INTDIR)\vc60.idb"
	-@erase "$(OUTDIR)\quickMain.exe"

"$(OUTDIR)" :
    if not exist "$(OUTDIR)/$(NULL)" mkdir "$(OUTDIR)"

CPP=cl.exe
CPP_PROJ=/nologo /ML /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_CONSOLE" /D "_MBCS" /Fp"$(INTDIR)\quickMain.pch" /YX /Fo"$(INTDIR)\\" /Fd"$(INTDIR)\\" /FD /c 

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
BSC32_FLAGS=/nologo /o"$(OUTDIR)\quickMain.bsc" 
BSC32_SBRS= \
	
LINK32=link.exe
LINK32_FLAGS=kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:console /incremental:no /pdb:"$(OUTDIR)\quickMain.pdb" /machine:I386 /out:"$(OUTDIR)\quickMain.exe" 
LINK32_OBJS= \
	"$(INTDIR)\myMain.obj" \
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

"$(OUTDIR)\quickMain.exe" : "$(OUTDIR)" $(DEF_FILE) $(LINK32_OBJS)
    $(LINK32) @<<
  $(LINK32_FLAGS) $(LINK32_OBJS)
<<

!ELSEIF  "$(CFG)" == "quickMain - Win32 Debug"

OUTDIR=.\..\..\bin
INTDIR=.\..\..\obj\quickMain
# Begin Custom Macros
OutDir=.\..\..\bin
# End Custom Macros

!IF "$(RECURSE)" == "0" 

ALL : "$(OUTDIR)\quickMain.exe"

!ELSE 

ALL : "utility - Win32 Debug" "tagged - Win32 Debug" "system - Win32 Debug" "superLU - Win32 Debug" "recorder - Win32 Debug" "nonlinearBeamColumn - Win32 Debug" "modelbuilder - Win32 Debug" "matrix - Win32 Debug" "material - Win32 Debug" "graph - Win32 Debug" "element - Win32 Debug" "domain - Win32 Debug" "convergence - Win32 Debug" "cblas - Win32 Debug" "analysis - Win32 Debug" "actor - Win32 Debug" "$(OUTDIR)\quickMain.exe"

!ENDIF 

!IF "$(RECURSE)" == "1" 
CLEAN :"actor - Win32 DebugCLEAN" "analysis - Win32 DebugCLEAN" "cblas - Win32 DebugCLEAN" "convergence - Win32 DebugCLEAN" "domain - Win32 DebugCLEAN" "element - Win32 DebugCLEAN" "graph - Win32 DebugCLEAN" "material - Win32 DebugCLEAN" "matrix - Win32 DebugCLEAN" "modelbuilder - Win32 DebugCLEAN" "nonlinearBeamColumn - Win32 DebugCLEAN" "recorder - Win32 DebugCLEAN" "superLU - Win32 DebugCLEAN" "system - Win32 DebugCLEAN" "tagged - Win32 DebugCLEAN" "utility - Win32 DebugCLEAN" 
!ELSE 
CLEAN :
!ENDIF 
	-@erase "$(INTDIR)\myMain.obj"
	-@erase "$(INTDIR)\vc60.idb"
	-@erase "$(OUTDIR)\quickMain.exe"

"$(OUTDIR)" :
    if not exist "$(OUTDIR)/$(NULL)" mkdir "$(OUTDIR)"

"$(INTDIR)" :
    if not exist "$(INTDIR)/$(NULL)" mkdir "$(INTDIR)"

CPP=cl.exe
CPP_PROJ=/nologo /MLd /W3 /GX /Od /I "..\..\src\domain\groundMotion" /I "..\..\src\domain\load" /I "..\..\src\domain\pattern" /I "..\..\src\element" /I "..\..\src\element\truss" /I "..\..\src\graph\graph" /I "..\..\src\graph\numberer" /I "..\..\src\system_of_eqn" /I "..\..\src\system_of_eqn\linearSOE" /I "..\..\src\system_of_eqn\linearSOE\sparseGen" /I "..\..\src\analysis\analysis" /I "..\..\src\convergenceTest" /I "..\..\src\analysis\integrator" /I "..\..\src\analysis\numberer" /I "..\..\src\domain\constraints" /I "..\..\src\analysis\handler" /I "..\..\src\analysis\algorithm" /I "..\..\src\analysis\model\simple" /I "..\..\src\matrix" /I "..\..\src" /I "..\..\src\actor\actor" /I "..\..\src\tagged" /I "..\..\src\domain\component" /I "..\..\src\analysis\algorithm\equiSolnAlgo" /I "..\..\src\analysis\model" /I "..\..\src\domain\domain" /I "..\..\src\domain\node" /I "..\..\src\modelbuilder" /I "..\..\src\utility" /D "WIN32" /D "_DEBUG" /D "_CONSOLE" /D "_MBCS" /Fp"$(INTDIR)\quickMain.pch" /YX /Fo"$(INTDIR)\\" /Fd"$(INTDIR)\\" /FD /GZ /c 

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
BSC32_FLAGS=/nologo /o"$(OUTDIR)\quickMain.bsc" 
BSC32_SBRS= \
	
LINK32=link.exe
LINK32_FLAGS=actor.lib analysis.lib cblas.lib convergence.lib domain.lib element.lib graph.lib material.lib matrix.lib modelbuilder.lib nonlinearBeamColumn.lib recorder.lib superLU.lib system.lib tagged.lib utility.lib kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:console /incremental:no /pdb:"$(OUTDIR)\quickMain.pdb" /machine:I386 /out:"$(OUTDIR)\quickMain.exe" /pdbtype:sept /libpath:"c:\msdev\lib" /libpath:"..\..\lib" 
LINK32_OBJS= \
	"$(INTDIR)\myMain.obj" \
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

"$(OUTDIR)\quickMain.exe" : "$(OUTDIR)" $(DEF_FILE) $(LINK32_OBJS)
    $(LINK32) @<<
  $(LINK32_FLAGS) $(LINK32_OBJS)
<<

!ENDIF 


!IF "$(NO_EXTERNAL_DEPS)" != "1"
!IF EXISTS("quickMain.dep")
!INCLUDE "quickMain.dep"
!ELSE 
!MESSAGE Warning: cannot find "quickMain.dep"
!ENDIF 
!ENDIF 


!IF "$(CFG)" == "quickMain - Win32 Release" || "$(CFG)" == "quickMain - Win32 Debug"
SOURCE=..\..\mhs\myMain.cpp

"$(INTDIR)\myMain.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!IF  "$(CFG)" == "quickMain - Win32 Release"

"actor - Win32 Release" : 
   cd "\g3\proj\actor"
   $(MAKE) /$(MAKEFLAGS) /F .\actor.mak CFG="actor - Win32 Release" 
   cd "..\quickMain"

"actor - Win32 ReleaseCLEAN" : 
   cd "\g3\proj\actor"
   $(MAKE) /$(MAKEFLAGS) /F .\actor.mak CFG="actor - Win32 Release" RECURSE=1 CLEAN 
   cd "..\quickMain"

!ELSEIF  "$(CFG)" == "quickMain - Win32 Debug"

"actor - Win32 Debug" : 
   cd "\g3\proj\actor"
   $(MAKE) /$(MAKEFLAGS) /F .\actor.mak CFG="actor - Win32 Debug" 
   cd "..\quickMain"

"actor - Win32 DebugCLEAN" : 
   cd "\g3\proj\actor"
   $(MAKE) /$(MAKEFLAGS) /F .\actor.mak CFG="actor - Win32 Debug" RECURSE=1 CLEAN 
   cd "..\quickMain"

!ENDIF 

!IF  "$(CFG)" == "quickMain - Win32 Release"

"analysis - Win32 Release" : 
   cd "\g3\proj\analysis"
   $(MAKE) /$(MAKEFLAGS) /F .\analysis.mak CFG="analysis - Win32 Release" 
   cd "..\quickMain"

"analysis - Win32 ReleaseCLEAN" : 
   cd "\g3\proj\analysis"
   $(MAKE) /$(MAKEFLAGS) /F .\analysis.mak CFG="analysis - Win32 Release" RECURSE=1 CLEAN 
   cd "..\quickMain"

!ELSEIF  "$(CFG)" == "quickMain - Win32 Debug"

"analysis - Win32 Debug" : 
   cd "\g3\proj\analysis"
   $(MAKE) /$(MAKEFLAGS) /F .\analysis.mak CFG="analysis - Win32 Debug" 
   cd "..\quickMain"

"analysis - Win32 DebugCLEAN" : 
   cd "\g3\proj\analysis"
   $(MAKE) /$(MAKEFLAGS) /F .\analysis.mak CFG="analysis - Win32 Debug" RECURSE=1 CLEAN 
   cd "..\quickMain"

!ENDIF 

!IF  "$(CFG)" == "quickMain - Win32 Release"

"cblas - Win32 Release" : 
   cd "\g3\proj\cblas"
   $(MAKE) /$(MAKEFLAGS) /F .\cblas.mak CFG="cblas - Win32 Release" 
   cd "..\quickMain"

"cblas - Win32 ReleaseCLEAN" : 
   cd "\g3\proj\cblas"
   $(MAKE) /$(MAKEFLAGS) /F .\cblas.mak CFG="cblas - Win32 Release" RECURSE=1 CLEAN 
   cd "..\quickMain"

!ELSEIF  "$(CFG)" == "quickMain - Win32 Debug"

"cblas - Win32 Debug" : 
   cd "\g3\proj\cblas"
   $(MAKE) /$(MAKEFLAGS) /F .\cblas.mak CFG="cblas - Win32 Debug" 
   cd "..\quickMain"

"cblas - Win32 DebugCLEAN" : 
   cd "\g3\proj\cblas"
   $(MAKE) /$(MAKEFLAGS) /F .\cblas.mak CFG="cblas - Win32 Debug" RECURSE=1 CLEAN 
   cd "..\quickMain"

!ENDIF 

!IF  "$(CFG)" == "quickMain - Win32 Release"

"convergence - Win32 Release" : 
   cd "\g3\proj\convergence"
   $(MAKE) /$(MAKEFLAGS) /F .\convergence.mak CFG="convergence - Win32 Release" 
   cd "..\quickMain"

"convergence - Win32 ReleaseCLEAN" : 
   cd "\g3\proj\convergence"
   $(MAKE) /$(MAKEFLAGS) /F .\convergence.mak CFG="convergence - Win32 Release" RECURSE=1 CLEAN 
   cd "..\quickMain"

!ELSEIF  "$(CFG)" == "quickMain - Win32 Debug"

"convergence - Win32 Debug" : 
   cd "\g3\proj\convergence"
   $(MAKE) /$(MAKEFLAGS) /F .\convergence.mak CFG="convergence - Win32 Debug" 
   cd "..\quickMain"

"convergence - Win32 DebugCLEAN" : 
   cd "\g3\proj\convergence"
   $(MAKE) /$(MAKEFLAGS) /F .\convergence.mak CFG="convergence - Win32 Debug" RECURSE=1 CLEAN 
   cd "..\quickMain"

!ENDIF 

!IF  "$(CFG)" == "quickMain - Win32 Release"

"domain - Win32 Release" : 
   cd "\g3\proj\domain"
   $(MAKE) /$(MAKEFLAGS) /F .\domain.mak CFG="domain - Win32 Release" 
   cd "..\quickMain"

"domain - Win32 ReleaseCLEAN" : 
   cd "\g3\proj\domain"
   $(MAKE) /$(MAKEFLAGS) /F .\domain.mak CFG="domain - Win32 Release" RECURSE=1 CLEAN 
   cd "..\quickMain"

!ELSEIF  "$(CFG)" == "quickMain - Win32 Debug"

"domain - Win32 Debug" : 
   cd "\g3\proj\domain"
   $(MAKE) /$(MAKEFLAGS) /F .\domain.mak CFG="domain - Win32 Debug" 
   cd "..\quickMain"

"domain - Win32 DebugCLEAN" : 
   cd "\g3\proj\domain"
   $(MAKE) /$(MAKEFLAGS) /F .\domain.mak CFG="domain - Win32 Debug" RECURSE=1 CLEAN 
   cd "..\quickMain"

!ENDIF 

!IF  "$(CFG)" == "quickMain - Win32 Release"

"element - Win32 Release" : 
   cd "\g3\proj\element"
   $(MAKE) /$(MAKEFLAGS) /F .\element.mak CFG="element - Win32 Release" 
   cd "..\quickMain"

"element - Win32 ReleaseCLEAN" : 
   cd "\g3\proj\element"
   $(MAKE) /$(MAKEFLAGS) /F .\element.mak CFG="element - Win32 Release" RECURSE=1 CLEAN 
   cd "..\quickMain"

!ELSEIF  "$(CFG)" == "quickMain - Win32 Debug"

"element - Win32 Debug" : 
   cd "\g3\proj\element"
   $(MAKE) /$(MAKEFLAGS) /F .\element.mak CFG="element - Win32 Debug" 
   cd "..\quickMain"

"element - Win32 DebugCLEAN" : 
   cd "\g3\proj\element"
   $(MAKE) /$(MAKEFLAGS) /F .\element.mak CFG="element - Win32 Debug" RECURSE=1 CLEAN 
   cd "..\quickMain"

!ENDIF 

!IF  "$(CFG)" == "quickMain - Win32 Release"

"graph - Win32 Release" : 
   cd "\g3\proj\graph"
   $(MAKE) /$(MAKEFLAGS) /F .\graph.mak CFG="graph - Win32 Release" 
   cd "..\quickMain"

"graph - Win32 ReleaseCLEAN" : 
   cd "\g3\proj\graph"
   $(MAKE) /$(MAKEFLAGS) /F .\graph.mak CFG="graph - Win32 Release" RECURSE=1 CLEAN 
   cd "..\quickMain"

!ELSEIF  "$(CFG)" == "quickMain - Win32 Debug"

"graph - Win32 Debug" : 
   cd "\g3\proj\graph"
   $(MAKE) /$(MAKEFLAGS) /F .\graph.mak CFG="graph - Win32 Debug" 
   cd "..\quickMain"

"graph - Win32 DebugCLEAN" : 
   cd "\g3\proj\graph"
   $(MAKE) /$(MAKEFLAGS) /F .\graph.mak CFG="graph - Win32 Debug" RECURSE=1 CLEAN 
   cd "..\quickMain"

!ENDIF 

!IF  "$(CFG)" == "quickMain - Win32 Release"

"material - Win32 Release" : 
   cd "\g3\proj\material"
   $(MAKE) /$(MAKEFLAGS) /F .\material.mak CFG="material - Win32 Release" 
   cd "..\quickMain"

"material - Win32 ReleaseCLEAN" : 
   cd "\g3\proj\material"
   $(MAKE) /$(MAKEFLAGS) /F .\material.mak CFG="material - Win32 Release" RECURSE=1 CLEAN 
   cd "..\quickMain"

!ELSEIF  "$(CFG)" == "quickMain - Win32 Debug"

"material - Win32 Debug" : 
   cd "\g3\proj\material"
   $(MAKE) /$(MAKEFLAGS) /F .\material.mak CFG="material - Win32 Debug" 
   cd "..\quickMain"

"material - Win32 DebugCLEAN" : 
   cd "\g3\proj\material"
   $(MAKE) /$(MAKEFLAGS) /F .\material.mak CFG="material - Win32 Debug" RECURSE=1 CLEAN 
   cd "..\quickMain"

!ENDIF 

!IF  "$(CFG)" == "quickMain - Win32 Release"

"matrix - Win32 Release" : 
   cd "\g3\proj\matrix"
   $(MAKE) /$(MAKEFLAGS) /F .\matrix.mak CFG="matrix - Win32 Release" 
   cd "..\quickMain"

"matrix - Win32 ReleaseCLEAN" : 
   cd "\g3\proj\matrix"
   $(MAKE) /$(MAKEFLAGS) /F .\matrix.mak CFG="matrix - Win32 Release" RECURSE=1 CLEAN 
   cd "..\quickMain"

!ELSEIF  "$(CFG)" == "quickMain - Win32 Debug"

"matrix - Win32 Debug" : 
   cd "\g3\proj\matrix"
   $(MAKE) /$(MAKEFLAGS) /F .\matrix.mak CFG="matrix - Win32 Debug" 
   cd "..\quickMain"

"matrix - Win32 DebugCLEAN" : 
   cd "\g3\proj\matrix"
   $(MAKE) /$(MAKEFLAGS) /F .\matrix.mak CFG="matrix - Win32 Debug" RECURSE=1 CLEAN 
   cd "..\quickMain"

!ENDIF 

!IF  "$(CFG)" == "quickMain - Win32 Release"

"modelbuilder - Win32 Release" : 
   cd "\g3\proj\modelbuilder"
   $(MAKE) /$(MAKEFLAGS) /F .\modelbuilder.mak CFG="modelbuilder - Win32 Release" 
   cd "..\quickMain"

"modelbuilder - Win32 ReleaseCLEAN" : 
   cd "\g3\proj\modelbuilder"
   $(MAKE) /$(MAKEFLAGS) /F .\modelbuilder.mak CFG="modelbuilder - Win32 Release" RECURSE=1 CLEAN 
   cd "..\quickMain"

!ELSEIF  "$(CFG)" == "quickMain - Win32 Debug"

"modelbuilder - Win32 Debug" : 
   cd "\g3\proj\modelbuilder"
   $(MAKE) /$(MAKEFLAGS) /F .\modelbuilder.mak CFG="modelbuilder - Win32 Debug" 
   cd "..\quickMain"

"modelbuilder - Win32 DebugCLEAN" : 
   cd "\g3\proj\modelbuilder"
   $(MAKE) /$(MAKEFLAGS) /F .\modelbuilder.mak CFG="modelbuilder - Win32 Debug" RECURSE=1 CLEAN 
   cd "..\quickMain"

!ENDIF 

!IF  "$(CFG)" == "quickMain - Win32 Release"

"nonlinearBeamColumn - Win32 Release" : 
   cd "\g3\proj\nonlinearBeamColumn"
   $(MAKE) /$(MAKEFLAGS) /F .\nonlinearBeamColumn.mak CFG="nonlinearBeamColumn - Win32 Release" 
   cd "..\quickMain"

"nonlinearBeamColumn - Win32 ReleaseCLEAN" : 
   cd "\g3\proj\nonlinearBeamColumn"
   $(MAKE) /$(MAKEFLAGS) /F .\nonlinearBeamColumn.mak CFG="nonlinearBeamColumn - Win32 Release" RECURSE=1 CLEAN 
   cd "..\quickMain"

!ELSEIF  "$(CFG)" == "quickMain - Win32 Debug"

"nonlinearBeamColumn - Win32 Debug" : 
   cd "\g3\proj\nonlinearBeamColumn"
   $(MAKE) /$(MAKEFLAGS) /F .\nonlinearBeamColumn.mak CFG="nonlinearBeamColumn - Win32 Debug" 
   cd "..\quickMain"

"nonlinearBeamColumn - Win32 DebugCLEAN" : 
   cd "\g3\proj\nonlinearBeamColumn"
   $(MAKE) /$(MAKEFLAGS) /F .\nonlinearBeamColumn.mak CFG="nonlinearBeamColumn - Win32 Debug" RECURSE=1 CLEAN 
   cd "..\quickMain"

!ENDIF 

!IF  "$(CFG)" == "quickMain - Win32 Release"

"recorder - Win32 Release" : 
   cd "\g3\proj\recorder"
   $(MAKE) /$(MAKEFLAGS) /F .\recorder.mak CFG="recorder - Win32 Release" 
   cd "..\quickMain"

"recorder - Win32 ReleaseCLEAN" : 
   cd "\g3\proj\recorder"
   $(MAKE) /$(MAKEFLAGS) /F .\recorder.mak CFG="recorder - Win32 Release" RECURSE=1 CLEAN 
   cd "..\quickMain"

!ELSEIF  "$(CFG)" == "quickMain - Win32 Debug"

"recorder - Win32 Debug" : 
   cd "\g3\proj\recorder"
   $(MAKE) /$(MAKEFLAGS) /F .\recorder.mak CFG="recorder - Win32 Debug" 
   cd "..\quickMain"

"recorder - Win32 DebugCLEAN" : 
   cd "\g3\proj\recorder"
   $(MAKE) /$(MAKEFLAGS) /F .\recorder.mak CFG="recorder - Win32 Debug" RECURSE=1 CLEAN 
   cd "..\quickMain"

!ENDIF 

!IF  "$(CFG)" == "quickMain - Win32 Release"

"superLU - Win32 Release" : 
   cd "\g3\proj\superLU"
   $(MAKE) /$(MAKEFLAGS) /F .\superLU.mak CFG="superLU - Win32 Release" 
   cd "..\quickMain"

"superLU - Win32 ReleaseCLEAN" : 
   cd "\g3\proj\superLU"
   $(MAKE) /$(MAKEFLAGS) /F .\superLU.mak CFG="superLU - Win32 Release" RECURSE=1 CLEAN 
   cd "..\quickMain"

!ELSEIF  "$(CFG)" == "quickMain - Win32 Debug"

"superLU - Win32 Debug" : 
   cd "\g3\proj\superLU"
   $(MAKE) /$(MAKEFLAGS) /F .\superLU.mak CFG="superLU - Win32 Debug" 
   cd "..\quickMain"

"superLU - Win32 DebugCLEAN" : 
   cd "\g3\proj\superLU"
   $(MAKE) /$(MAKEFLAGS) /F .\superLU.mak CFG="superLU - Win32 Debug" RECURSE=1 CLEAN 
   cd "..\quickMain"

!ENDIF 

!IF  "$(CFG)" == "quickMain - Win32 Release"

"system - Win32 Release" : 
   cd "\g3\proj\system"
   $(MAKE) /$(MAKEFLAGS) /F .\system.mak CFG="system - Win32 Release" 
   cd "..\quickMain"

"system - Win32 ReleaseCLEAN" : 
   cd "\g3\proj\system"
   $(MAKE) /$(MAKEFLAGS) /F .\system.mak CFG="system - Win32 Release" RECURSE=1 CLEAN 
   cd "..\quickMain"

!ELSEIF  "$(CFG)" == "quickMain - Win32 Debug"

"system - Win32 Debug" : 
   cd "\g3\proj\system"
   $(MAKE) /$(MAKEFLAGS) /F .\system.mak CFG="system - Win32 Debug" 
   cd "..\quickMain"

"system - Win32 DebugCLEAN" : 
   cd "\g3\proj\system"
   $(MAKE) /$(MAKEFLAGS) /F .\system.mak CFG="system - Win32 Debug" RECURSE=1 CLEAN 
   cd "..\quickMain"

!ENDIF 

!IF  "$(CFG)" == "quickMain - Win32 Release"

"tagged - Win32 Release" : 
   cd "\g3\proj\tagged"
   $(MAKE) /$(MAKEFLAGS) /F .\tagged.mak CFG="tagged - Win32 Release" 
   cd "..\quickMain"

"tagged - Win32 ReleaseCLEAN" : 
   cd "\g3\proj\tagged"
   $(MAKE) /$(MAKEFLAGS) /F .\tagged.mak CFG="tagged - Win32 Release" RECURSE=1 CLEAN 
   cd "..\quickMain"

!ELSEIF  "$(CFG)" == "quickMain - Win32 Debug"

"tagged - Win32 Debug" : 
   cd "\g3\proj\tagged"
   $(MAKE) /$(MAKEFLAGS) /F .\tagged.mak CFG="tagged - Win32 Debug" 
   cd "..\quickMain"

"tagged - Win32 DebugCLEAN" : 
   cd "\g3\proj\tagged"
   $(MAKE) /$(MAKEFLAGS) /F .\tagged.mak CFG="tagged - Win32 Debug" RECURSE=1 CLEAN 
   cd "..\quickMain"

!ENDIF 

!IF  "$(CFG)" == "quickMain - Win32 Release"

"utility - Win32 Release" : 
   cd "\g3\proj\utility"
   $(MAKE) /$(MAKEFLAGS) /F .\utility.mak CFG="utility - Win32 Release" 
   cd "..\quickMain"

"utility - Win32 ReleaseCLEAN" : 
   cd "\g3\proj\utility"
   $(MAKE) /$(MAKEFLAGS) /F .\utility.mak CFG="utility - Win32 Release" RECURSE=1 CLEAN 
   cd "..\quickMain"

!ELSEIF  "$(CFG)" == "quickMain - Win32 Debug"

"utility - Win32 Debug" : 
   cd "\g3\proj\utility"
   $(MAKE) /$(MAKEFLAGS) /F .\utility.mak CFG="utility - Win32 Debug" 
   cd "..\quickMain"

"utility - Win32 DebugCLEAN" : 
   cd "\g3\proj\utility"
   $(MAKE) /$(MAKEFLAGS) /F .\utility.mak CFG="utility - Win32 Debug" RECURSE=1 CLEAN 
   cd "..\quickMain"

!ENDIF 


!ENDIF 

