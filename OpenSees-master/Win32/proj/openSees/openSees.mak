# Microsoft Developer Studio Generated NMAKE File, Based on openSees.dsp
!IF "$(CFG)" == ""
CFG=openSees - Win32 Debug
!MESSAGE No configuration specified. Defaulting to openSees - Win32 Debug.
!ENDIF 

!IF "$(CFG)" != "openSees - Win32 Release" && "$(CFG)" != "openSees - Win32 Debug"
!MESSAGE Invalid configuration "$(CFG)" specified.
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
!ERROR An invalid configuration is specified.
!ENDIF 

!IF "$(OS)" == "Windows_NT"
NULL=
!ELSE 
NULL=nul
!ENDIF 

!IF  "$(CFG)" == "openSees - Win32 Release"

OUTDIR=.\..\..\bin
INTDIR=.\..\..\obj\openSees
# Begin Custom Macros
OutDir=.\..\..\bin
# End Custom Macros

!IF "$(RECURSE)" == "0" 

ALL : "$(OUTDIR)\openSees.exe"

!ELSE 

ALL : "utility - Win32 Release" "tagged - Win32 Release" "system - Win32 Release" "superLU - Win32 Release" "string - Win32 Release" "renderer - Win32 Release" "recorder - Win32 Release" "nonlinearBeamColumn - Win32 Release" "modelbuilder - Win32 Release" "matrix - Win32 Release" "material - Win32 Release" "handler - Win32 Release" "graph - Win32 Release" "element - Win32 Release" "domain - Win32 Release" "database - Win32 Release" "convergence - Win32 Release" "cblas - Win32 Release" "analysis - Win32 Release" "actor - Win32 Release" "$(OUTDIR)\openSees.exe"

!ENDIF 

!IF "$(RECURSE)" == "1" 
CLEAN :"actor - Win32 ReleaseCLEAN" "analysis - Win32 ReleaseCLEAN" "cblas - Win32 ReleaseCLEAN" "convergence - Win32 ReleaseCLEAN" "database - Win32 ReleaseCLEAN" "domain - Win32 ReleaseCLEAN" "element - Win32 ReleaseCLEAN" "graph - Win32 ReleaseCLEAN" "handler - Win32 ReleaseCLEAN" "material - Win32 ReleaseCLEAN" "matrix - Win32 ReleaseCLEAN" "modelbuilder - Win32 ReleaseCLEAN" "nonlinearBeamColumn - Win32 ReleaseCLEAN" "recorder - Win32 ReleaseCLEAN" "renderer - Win32 ReleaseCLEAN" "string - Win32 ReleaseCLEAN" "superLU - Win32 ReleaseCLEAN" "system - Win32 ReleaseCLEAN" "tagged - Win32 ReleaseCLEAN" "utility - Win32 ReleaseCLEAN" 
!ELSE 
CLEAN :
!ENDIF 
	-@erase "$(INTDIR)\commands.obj"
	-@erase "$(INTDIR)\myCommands.obj"
	-@erase "$(INTDIR)\tclAppInit.obj"
	-@erase "$(INTDIR)\TclFeViewer.obj"
	-@erase "$(INTDIR)\tclMain.obj"
	-@erase "$(INTDIR)\TclModelBuilder.obj"
	-@erase "$(INTDIR)\TclVideoPlayer.obj"
	-@erase "$(INTDIR)\vc60.idb"
	-@erase "$(OUTDIR)\openSees.exe"

"$(OUTDIR)" :
    if not exist "$(OUTDIR)/$(NULL)" mkdir "$(OUTDIR)"

"$(INTDIR)" :
    if not exist "$(INTDIR)/$(NULL)" mkdir "$(INTDIR)"

CPP=cl.exe
CPP_PROJ=/nologo /MLd /W3 /GX /O2 /I "C:\Program Files\tcl\include" /I "..\..\..\src\material\backbone" /I "..\..\..\src\material\state" /I "..\..\..\src\material\state\stiffness" /I "..\..\..\src\material\state\deformation" /I "..\..\..\src\material\state\strength" /I "..\..\..\src\coordTransformation" /I "..\..\..\src\material\section\repres\section" /I "..\..\..\src\analysis\algorithm\equiSolnAlgo" /I "..\..\..\src\system_of_eqn\eigenSOE" /I "..\..\..\src\analysis\algorithm\eigenAlgo" /I "..\..\..\src\material\nD" /I "..\..\..\src\material\uniaxial" /I "..\..\..\src\handler" /I "..\..\..\src\tcl" /I "..\..\..\src\actor\objectBroker" /I "..\..\..\src\system_of_eqn\linearSOE\umfGEN" /I "..\..\..\src\system_of_eqn\linearSOE\fullGEN" /I "..\..\..\src\system_of_eqn\linearSOE\sparseGEN" /I "..\..\..\src\system_of_eqn\linearSOE\bandSPD" /I "..\..\..\src\system_of_eqn\linearSOE\bandGEN" /I "..\..\..\src\element\nonlinearBeamColumn\tcl\repres\section" /I "..\..\..\src\matrix" /I "..\..\..\src\recorder" /I "..\..\..\src\graph\numberer" /I "..\..\..\src\material\section" /I "..\..\..\src\graph\graph" /I "..\..\..\src\element\beam2d" /I "..\..\..\src\element\beam3d" /I\
 "..\..\..\src\system_of_eqn" /I "..\..\..\src\system_of_eqn\linearSOE" /I "..\..\..\src\system_of_eqn\linearSOE\profileSPD" /I "..\..\..\src\system_of_eqn\linearSOE\sparseSYM" /I "..\..\..\src\domain\pattern" /I "..\..\..\src\analysis\analysis" /I "..\..\..\src\analysis\integrator" /I "..\..\..\src\analysis\numberer" /I "..\..\..\src\analysis\handler" /I "..\..\..\src\renderer" /I "..\..\..\src\material" /I "..\..\..\src\analysis\algorithm" /I "..\..\..\src\convergenceTest" /I "..\..\..\src\analysis\model\simple" /I "..\..\..\src\domain\load" /I "..\..\..\src\analysis\model" /I "..\..\..\src\element\truss" /I "..\..\..\src\actor\channel" /I "..\..\..\src\utility" /I "..\..\..\src\actor\actor" /I "..\..\..\src\modelbuilder" /I "..\..\..\src\modelbuilder\tcl" /I "..\..\..\src\domain\constraints" /I "..\..\..\src\domain\component" /I "..\..\..\src\element" /I "..\..\..\src\domain\node" /I "..\..\..\src\domain\domain" /I "..\..\..\src\tagged\storage" /I "..\..\..\src" /I "..\..\..\src\tagged" /D "WIN32" /D "NDEBUG" /D "_CONSOLE" /D "_MBCS" /D "_WIN32" /D "_FORTRAN" /D "BUILD_tcl" /Fo"$(INTDIR)\\" /Fd"$(INTDIR)\\" /FD /c 

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
BSC32_FLAGS=/nologo /o"$(OUTDIR)\openSees.bsc" 
BSC32_SBRS= \
	
LINK32=link.exe
LINK32_FLAGS=database.lib renderer.lib blas.lib lapack.lib feap.lib arpack.lib umfpack.lib openSeesFortran.lib actor.lib analysis.lib cblas.lib renderer.lib convergence.lib domain.lib element.lib graph.lib material.lib matrix.lib modelbuilder.lib nonlinearBeamColumn.lib recorder.lib superLU.lib system.lib tagged.lib utility.lib tcl82.lib kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib OpenGL32.lib glu32.lib GlAux.lib tcl82.lib /nologo /subsystem:console /pdb:none /machine:I386 /nodefaultlib:"libc.lib" /out:"$(OUTDIR)\openSees.exe" /libpath:"c:\Program Files\tcl\lib" /libpath:"c:\msdev\lib" /libpath:"..\..\lib" /FORCE:MULTIPLE 
LINK32_OBJS= \
	"$(INTDIR)\commands.obj" \
	"$(INTDIR)\myCommands.obj" \
	"$(INTDIR)\tclAppInit.obj" \
	"$(INTDIR)\TclFeViewer.obj" \
	"$(INTDIR)\tclMain.obj" \
	"$(INTDIR)\TclModelBuilder.obj" \
	"$(INTDIR)\TclVideoPlayer.obj" \
	"..\..\lib\actor.lib" \
	"..\..\lib\analysis.lib" \
	"..\..\lib\cblas.lib" \
	"..\..\lib\convergence.lib" \
	"..\..\lib\database.lib" \
	"..\..\lib\domain.lib" \
	"..\..\lib\element.lib" \
	"..\..\lib\graph.lib" \
	"..\..\lib\handler.lib" \
	"..\..\lib\material.lib" \
	"..\..\lib\matrix.lib" \
	"..\..\lib\modelbuilder.lib" \
	"..\..\lib\nonlinearBeamColumn.lib" \
	"..\..\lib\recorder.lib" \
	"..\..\lib\renderer.lib" \
	"..\..\lib\string.lib" \
	"..\..\lib\superLU.lib" \
	"..\..\lib\system.lib" \
	"..\..\lib\tagged.lib" \
	"..\..\lib\utility.lib"

"$(OUTDIR)\openSees.exe" : "$(OUTDIR)" $(DEF_FILE) $(LINK32_OBJS)
    $(LINK32) @<<
  $(LINK32_FLAGS) $(LINK32_OBJS)
<<

!ELSEIF  "$(CFG)" == "openSees - Win32 Debug"

OUTDIR=.\..\..\bin
INTDIR=.\..\..\obj\openSees
# Begin Custom Macros
OutDir=.\..\..\bin
# End Custom Macros

!IF "$(RECURSE)" == "0" 

ALL : "$(OUTDIR)\openSees.exe" "$(OUTDIR)\openSees.bsc"

!ELSE 

ALL : "utility - Win32 Debug" "tagged - Win32 Debug" "system - Win32 Debug" "superLU - Win32 Debug" "string - Win32 Debug" "renderer - Win32 Debug" "recorder - Win32 Debug" "nonlinearBeamColumn - Win32 Debug" "modelbuilder - Win32 Debug" "matrix - Win32 Debug" "material - Win32 Debug" "handler - Win32 Debug" "graph - Win32 Debug" "element - Win32 Debug" "domain - Win32 Debug" "database - Win32 Debug" "convergence - Win32 Debug" "cblas - Win32 Debug" "analysis - Win32 Debug" "actor - Win32 Debug" "$(OUTDIR)\openSees.exe" "$(OUTDIR)\openSees.bsc"

!ENDIF 

!IF "$(RECURSE)" == "1" 
CLEAN :"actor - Win32 DebugCLEAN" "analysis - Win32 DebugCLEAN" "cblas - Win32 DebugCLEAN" "convergence - Win32 DebugCLEAN" "database - Win32 DebugCLEAN" "domain - Win32 DebugCLEAN" "element - Win32 DebugCLEAN" "graph - Win32 DebugCLEAN" "handler - Win32 DebugCLEAN" "material - Win32 DebugCLEAN" "matrix - Win32 DebugCLEAN" "modelbuilder - Win32 DebugCLEAN" "nonlinearBeamColumn - Win32 DebugCLEAN" "recorder - Win32 DebugCLEAN" "renderer - Win32 DebugCLEAN" "string - Win32 DebugCLEAN" "superLU - Win32 DebugCLEAN" "system - Win32 DebugCLEAN" "tagged - Win32 DebugCLEAN" "utility - Win32 DebugCLEAN" 
!ELSE 
CLEAN :
!ENDIF 
	-@erase "$(INTDIR)\commands.obj"
	-@erase "$(INTDIR)\commands.sbr"
	-@erase "$(INTDIR)\myCommands.obj"
	-@erase "$(INTDIR)\myCommands.sbr"
	-@erase "$(INTDIR)\tclAppInit.obj"
	-@erase "$(INTDIR)\tclAppInit.sbr"
	-@erase "$(INTDIR)\TclFeViewer.obj"
	-@erase "$(INTDIR)\TclFeViewer.sbr"
	-@erase "$(INTDIR)\tclMain.obj"
	-@erase "$(INTDIR)\tclMain.sbr"
	-@erase "$(INTDIR)\TclModelBuilder.obj"
	-@erase "$(INTDIR)\TclModelBuilder.sbr"
	-@erase "$(INTDIR)\TclVideoPlayer.obj"
	-@erase "$(INTDIR)\TclVideoPlayer.sbr"
	-@erase "$(INTDIR)\vc60.idb"
	-@erase "$(INTDIR)\vc60.pdb"
	-@erase "$(OUTDIR)\openSees.bsc"
	-@erase "$(OUTDIR)\openSees.exe"

"$(OUTDIR)" :
    if not exist "$(OUTDIR)/$(NULL)" mkdir "$(OUTDIR)"

"$(INTDIR)" :
    if not exist "$(INTDIR)/$(NULL)" mkdir "$(INTDIR)"

CPP=cl.exe
CPP_PROJ=/nologo /MLd /W3 /GX /ZI /Od /I "..\..\..\src\material\state\stiffness" /I "..\..\..\src\material\state\deformation" /I "..\..\..\src\material\state\strength" /I "C:\Program Files\tcl\include" /I "..\..\..\src\material\backbone" /I "..\..\..\src\material\state" /I "..\..\..\src\coordTransformation" /I "..\..\..\src\material\section\repres\section" /I "..\..\..\src\analysis\algorithm\equiSolnAlgo" /I "..\..\..\src\system_of_eqn\eigenSOE" /I "..\..\..\src\analysis\algorithm\eigenAlgo" /I "..\..\..\src\material\nD" /I "..\..\..\src\material\uniaxial" /I "..\..\..\src\handler" /I "..\..\..\src\tcl" /I "..\..\..\src\actor\objectBroker" /I "..\..\..\src\system_of_eqn\linearSOE\umfGEN" /I "..\..\..\src\system_of_eqn\linearSOE\fullGEN" /I "..\..\..\src\system_of_eqn\linearSOE\sparseGEN" /I "..\..\..\src\system_of_eqn\linearSOE\bandSPD" /I "..\..\..\src\system_of_eqn\linearSOE\bandGEN" /I "..\..\..\src\element\nonlinearBeamColumn\tcl\repres\section" /I "..\..\..\src\matrix" /I "..\..\..\src\recorder" /I "..\..\..\src\graph\numberer" /I "..\..\..\src\material\section" /I "..\..\..\src\graph\graph" /I "..\..\..\src\element\beam2d" /I "..\..\..\src\element\beam3d" /I\
 "..\..\..\src\system_of_eqn" /I "..\..\..\src\system_of_eqn\linearSOE" /I "..\..\..\src\system_of_eqn\linearSOE\profileSPD" /I "..\..\..\src\system_of_eqn\linearSOE\sparseSYM" /I "..\..\..\src\domain\pattern" /I "..\..\..\src\analysis\analysis" /I "..\..\..\src\analysis\integrator" /I "..\..\..\src\analysis\numberer" /I "..\..\..\src\analysis\handler" /I "..\..\..\src\renderer" /I "..\..\..\src\material" /I "..\..\..\src\analysis\algorithm" /I "..\..\..\src\convergenceTest" /I "..\..\..\src\analysis\model\simple" /I "..\..\..\src\domain\load" /I "..\..\..\src\analysis\model" /I "..\..\..\src\element\truss" /I "..\..\..\src\actor\channel" /I "..\..\..\src\utility" /I "..\..\..\src\actor\actor" /I "..\..\..\src\modelbuilder" /I "..\..\..\src\modelbuilder\tcl" /I "..\..\..\src\domain\constraints" /I "..\..\..\src\domain\component" /I "..\..\..\src\element" /I "..\..\..\src\domain\node" /I "..\..\..\src\domain\domain" /I "..\..\..\src\tagged\storage" /I "..\..\..\src" /I "..\..\..\src\tagged" /D "_WIN32" /D "_FORTRAN" /D "WIN32" /D "_DEBUG" /D "_CONSOLE" /D "_MBCS" /D "BUILD_tcl" /FR"$(INTDIR)\\" /Fo"$(INTDIR)\\" /Fd"$(INTDIR)\\" /FD /I /GZ 

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
BSC32_FLAGS=/nologo /o"$(OUTDIR)\openSees.bsc" 
BSC32_SBRS= \
	"$(INTDIR)\commands.sbr" \
	"$(INTDIR)\myCommands.sbr" \
	"$(INTDIR)\tclAppInit.sbr" \
	"$(INTDIR)\TclFeViewer.sbr" \
	"$(INTDIR)\tclMain.sbr" \
	"$(INTDIR)\TclModelBuilder.sbr" \
	"$(INTDIR)\TclVideoPlayer.sbr"

"$(OUTDIR)\openSees.bsc" : "$(OUTDIR)" $(BSC32_SBRS)
    $(BSC32) @<<
  $(BSC32_FLAGS) $(BSC32_SBRS)
<<

LINK32=link.exe
LINK32_FLAGS=database.lib renderer.lib blas.lib lapack.lib feap.lib arpack.lib umfpack.lib openSeesFortran.lib actor.lib analysis.lib cblas.lib renderer.lib convergence.lib domain.lib element.lib graph.lib material.lib matrix.lib modelbuilder.lib nonlinearBeamColumn.lib recorder.lib superLU.lib system.lib tagged.lib utility.lib tcl82.lib kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib OpenGL32.lib glu32.lib GlAux.lib tcl82.lib /nologo /subsystem:console /pdb:none /debug /machine:I386 /nodefaultlib:"libc.lib" /out:"$(OUTDIR)\openSees.exe" /libpath:"c:\Program Files\tcl\lib" /libpath:"c:\msdev\lib" /libpath:"..\..\lib" /FORCE:MULTIPLE 
LINK32_OBJS= \
	"$(INTDIR)\commands.obj" \
	"$(INTDIR)\myCommands.obj" \
	"$(INTDIR)\tclAppInit.obj" \
	"$(INTDIR)\TclFeViewer.obj" \
	"$(INTDIR)\tclMain.obj" \
	"$(INTDIR)\TclModelBuilder.obj" \
	"$(INTDIR)\TclVideoPlayer.obj" \
	"..\..\lib\actor.lib" \
	"..\..\lib\analysis.lib" \
	"..\..\lib\cblas.lib" \
	"..\..\lib\convergence.lib" \
	"..\..\lib\database.lib" \
	"..\..\lib\domain.lib" \
	"..\..\lib\element.lib" \
	"..\..\lib\graph.lib" \
	"..\..\lib\handler.lib" \
	"..\..\lib\material.lib" \
	"..\..\lib\matrix.lib" \
	"..\..\lib\modelbuilder.lib" \
	"..\..\lib\nonlinearBeamColumn.lib" \
	"..\..\lib\recorder.lib" \
	"..\..\lib\renderer.lib" \
	"..\..\lib\string.lib" \
	"..\..\lib\superLU.lib" \
	"..\..\lib\system.lib" \
	"..\..\lib\tagged.lib" \
	"..\..\lib\utility.lib"

"$(OUTDIR)\openSees.exe" : "$(OUTDIR)" $(DEF_FILE) $(LINK32_OBJS)
    $(LINK32) @<<
  $(LINK32_FLAGS) $(LINK32_OBJS)
<<

!ENDIF 


!IF "$(NO_EXTERNAL_DEPS)" != "1"
!IF EXISTS("openSees.dep")
!INCLUDE "openSees.dep"
!ELSE 
!MESSAGE Warning: cannot find "openSees.dep"
!ENDIF 
!ENDIF 


!IF "$(CFG)" == "openSees - Win32 Release" || "$(CFG)" == "openSees - Win32 Debug"
SOURCE=..\..\..\SRC\tcl\commands.cpp

!IF  "$(CFG)" == "openSees - Win32 Release"


"$(INTDIR)\commands.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "openSees - Win32 Debug"


"$(INTDIR)\commands.obj"	"$(INTDIR)\commands.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\..\SRC\modelbuilder\tcl\myCommands.cpp

!IF  "$(CFG)" == "openSees - Win32 Release"


"$(INTDIR)\myCommands.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "openSees - Win32 Debug"


"$(INTDIR)\myCommands.obj"	"$(INTDIR)\myCommands.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\..\SRC\tcl\tclAppInit.cpp

!IF  "$(CFG)" == "openSees - Win32 Release"


"$(INTDIR)\tclAppInit.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "openSees - Win32 Debug"


"$(INTDIR)\tclAppInit.obj"	"$(INTDIR)\tclAppInit.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\..\SRC\tcl\TclFeViewer.cpp

!IF  "$(CFG)" == "openSees - Win32 Release"


"$(INTDIR)\TclFeViewer.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "openSees - Win32 Debug"


"$(INTDIR)\TclFeViewer.obj"	"$(INTDIR)\TclFeViewer.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\..\SRC\tcl\tclMain.cpp

!IF  "$(CFG)" == "openSees - Win32 Release"


"$(INTDIR)\tclMain.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "openSees - Win32 Debug"


"$(INTDIR)\tclMain.obj"	"$(INTDIR)\tclMain.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\..\SRC\modelbuilder\tcl\TclModelBuilder.cpp

!IF  "$(CFG)" == "openSees - Win32 Release"


"$(INTDIR)\TclModelBuilder.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "openSees - Win32 Debug"


"$(INTDIR)\TclModelBuilder.obj"	"$(INTDIR)\TclModelBuilder.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\..\SRC\tcl\TclVideoPlayer.cpp

!IF  "$(CFG)" == "openSees - Win32 Release"


"$(INTDIR)\TclVideoPlayer.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "openSees - Win32 Debug"


"$(INTDIR)\TclVideoPlayer.obj"	"$(INTDIR)\TclVideoPlayer.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

!IF  "$(CFG)" == "openSees - Win32 Release"

"actor - Win32 Release" : 
   cd "\OpenSees\Win32\proj\actor"
   $(MAKE) /$(MAKEFLAGS) /F .\actor.mak CFG="actor - Win32 Release" 
   cd "..\openSees"

"actor - Win32 ReleaseCLEAN" : 
   cd "\OpenSees\Win32\proj\actor"
   $(MAKE) /$(MAKEFLAGS) /F .\actor.mak CFG="actor - Win32 Release" RECURSE=1 CLEAN 
   cd "..\openSees"

!ELSEIF  "$(CFG)" == "openSees - Win32 Debug"

"actor - Win32 Debug" : 
   cd "\OpenSees\Win32\proj\actor"
   $(MAKE) /$(MAKEFLAGS) /F .\actor.mak CFG="actor - Win32 Debug" 
   cd "..\openSees"

"actor - Win32 DebugCLEAN" : 
   cd "\OpenSees\Win32\proj\actor"
   $(MAKE) /$(MAKEFLAGS) /F .\actor.mak CFG="actor - Win32 Debug" RECURSE=1 CLEAN 
   cd "..\openSees"

!ENDIF 

!IF  "$(CFG)" == "openSees - Win32 Release"

"analysis - Win32 Release" : 
   cd "\OpenSees\Win32\proj\analysis"
   $(MAKE) /$(MAKEFLAGS) /F .\analysis.mak CFG="analysis - Win32 Release" 
   cd "..\openSees"

"analysis - Win32 ReleaseCLEAN" : 
   cd "\OpenSees\Win32\proj\analysis"
   $(MAKE) /$(MAKEFLAGS) /F .\analysis.mak CFG="analysis - Win32 Release" RECURSE=1 CLEAN 
   cd "..\openSees"

!ELSEIF  "$(CFG)" == "openSees - Win32 Debug"

"analysis - Win32 Debug" : 
   cd "\OpenSees\Win32\proj\analysis"
   $(MAKE) /$(MAKEFLAGS) /F .\analysis.mak CFG="analysis - Win32 Debug" 
   cd "..\openSees"

"analysis - Win32 DebugCLEAN" : 
   cd "\OpenSees\Win32\proj\analysis"
   $(MAKE) /$(MAKEFLAGS) /F .\analysis.mak CFG="analysis - Win32 Debug" RECURSE=1 CLEAN 
   cd "..\openSees"

!ENDIF 

!IF  "$(CFG)" == "openSees - Win32 Release"

"cblas - Win32 Release" : 
   cd "\OpenSees\Win32\proj\cblas"
   $(MAKE) /$(MAKEFLAGS) /F .\cblas.mak CFG="cblas - Win32 Release" 
   cd "..\openSees"

"cblas - Win32 ReleaseCLEAN" : 
   cd "\OpenSees\Win32\proj\cblas"
   $(MAKE) /$(MAKEFLAGS) /F .\cblas.mak CFG="cblas - Win32 Release" RECURSE=1 CLEAN 
   cd "..\openSees"

!ELSEIF  "$(CFG)" == "openSees - Win32 Debug"

"cblas - Win32 Debug" : 
   cd "\OpenSees\Win32\proj\cblas"
   $(MAKE) /$(MAKEFLAGS) /F .\cblas.mak CFG="cblas - Win32 Debug" 
   cd "..\openSees"

"cblas - Win32 DebugCLEAN" : 
   cd "\OpenSees\Win32\proj\cblas"
   $(MAKE) /$(MAKEFLAGS) /F .\cblas.mak CFG="cblas - Win32 Debug" RECURSE=1 CLEAN 
   cd "..\openSees"

!ENDIF 

!IF  "$(CFG)" == "openSees - Win32 Release"

"convergence - Win32 Release" : 
   cd "\OpenSees\Win32\proj\convergence"
   $(MAKE) /$(MAKEFLAGS) /F .\convergence.mak CFG="convergence - Win32 Release" 
   cd "..\openSees"

"convergence - Win32 ReleaseCLEAN" : 
   cd "\OpenSees\Win32\proj\convergence"
   $(MAKE) /$(MAKEFLAGS) /F .\convergence.mak CFG="convergence - Win32 Release" RECURSE=1 CLEAN 
   cd "..\openSees"

!ELSEIF  "$(CFG)" == "openSees - Win32 Debug"

"convergence - Win32 Debug" : 
   cd "\OpenSees\Win32\proj\convergence"
   $(MAKE) /$(MAKEFLAGS) /F .\convergence.mak CFG="convergence - Win32 Debug" 
   cd "..\openSees"

"convergence - Win32 DebugCLEAN" : 
   cd "\OpenSees\Win32\proj\convergence"
   $(MAKE) /$(MAKEFLAGS) /F .\convergence.mak CFG="convergence - Win32 Debug" RECURSE=1 CLEAN 
   cd "..\openSees"

!ENDIF 

!IF  "$(CFG)" == "openSees - Win32 Release"

"database - Win32 Release" : 
   cd "\OpenSees\Win32\proj\database"
   $(MAKE) /$(MAKEFLAGS) /F .\database.mak CFG="database - Win32 Release" 
   cd "..\openSees"

"database - Win32 ReleaseCLEAN" : 
   cd "\OpenSees\Win32\proj\database"
   $(MAKE) /$(MAKEFLAGS) /F .\database.mak CFG="database - Win32 Release" RECURSE=1 CLEAN 
   cd "..\openSees"

!ELSEIF  "$(CFG)" == "openSees - Win32 Debug"

"database - Win32 Debug" : 
   cd "\OpenSees\Win32\proj\database"
   $(MAKE) /$(MAKEFLAGS) /F .\database.mak CFG="database - Win32 Debug" 
   cd "..\openSees"

"database - Win32 DebugCLEAN" : 
   cd "\OpenSees\Win32\proj\database"
   $(MAKE) /$(MAKEFLAGS) /F .\database.mak CFG="database - Win32 Debug" RECURSE=1 CLEAN 
   cd "..\openSees"

!ENDIF 

!IF  "$(CFG)" == "openSees - Win32 Release"

"domain - Win32 Release" : 
   cd "\OpenSees\Win32\proj\domain"
   $(MAKE) /$(MAKEFLAGS) /F .\domain.mak CFG="domain - Win32 Release" 
   cd "..\openSees"

"domain - Win32 ReleaseCLEAN" : 
   cd "\OpenSees\Win32\proj\domain"
   $(MAKE) /$(MAKEFLAGS) /F .\domain.mak CFG="domain - Win32 Release" RECURSE=1 CLEAN 
   cd "..\openSees"

!ELSEIF  "$(CFG)" == "openSees - Win32 Debug"

"domain - Win32 Debug" : 
   cd "\OpenSees\Win32\proj\domain"
   $(MAKE) /$(MAKEFLAGS) /F .\domain.mak CFG="domain - Win32 Debug" 
   cd "..\openSees"

"domain - Win32 DebugCLEAN" : 
   cd "\OpenSees\Win32\proj\domain"
   $(MAKE) /$(MAKEFLAGS) /F .\domain.mak CFG="domain - Win32 Debug" RECURSE=1 CLEAN 
   cd "..\openSees"

!ENDIF 

!IF  "$(CFG)" == "openSees - Win32 Release"

"element - Win32 Release" : 
   cd "\OpenSees\Win32\proj\element"
   $(MAKE) /$(MAKEFLAGS) /F .\element.mak CFG="element - Win32 Release" 
   cd "..\openSees"

"element - Win32 ReleaseCLEAN" : 
   cd "\OpenSees\Win32\proj\element"
   $(MAKE) /$(MAKEFLAGS) /F .\element.mak CFG="element - Win32 Release" RECURSE=1 CLEAN 
   cd "..\openSees"

!ELSEIF  "$(CFG)" == "openSees - Win32 Debug"

"element - Win32 Debug" : 
   cd "\OpenSees\Win32\proj\element"
   $(MAKE) /$(MAKEFLAGS) /F .\element.mak CFG="element - Win32 Debug" 
   cd "..\openSees"

"element - Win32 DebugCLEAN" : 
   cd "\OpenSees\Win32\proj\element"
   $(MAKE) /$(MAKEFLAGS) /F .\element.mak CFG="element - Win32 Debug" RECURSE=1 CLEAN 
   cd "..\openSees"

!ENDIF 

!IF  "$(CFG)" == "openSees - Win32 Release"

"graph - Win32 Release" : 
   cd "\OpenSees\Win32\proj\graph"
   $(MAKE) /$(MAKEFLAGS) /F .\graph.mak CFG="graph - Win32 Release" 
   cd "..\openSees"

"graph - Win32 ReleaseCLEAN" : 
   cd "\OpenSees\Win32\proj\graph"
   $(MAKE) /$(MAKEFLAGS) /F .\graph.mak CFG="graph - Win32 Release" RECURSE=1 CLEAN 
   cd "..\openSees"

!ELSEIF  "$(CFG)" == "openSees - Win32 Debug"

"graph - Win32 Debug" : 
   cd "\OpenSees\Win32\proj\graph"
   $(MAKE) /$(MAKEFLAGS) /F .\graph.mak CFG="graph - Win32 Debug" 
   cd "..\openSees"

"graph - Win32 DebugCLEAN" : 
   cd "\OpenSees\Win32\proj\graph"
   $(MAKE) /$(MAKEFLAGS) /F .\graph.mak CFG="graph - Win32 Debug" RECURSE=1 CLEAN 
   cd "..\openSees"

!ENDIF 

!IF  "$(CFG)" == "openSees - Win32 Release"

"handler - Win32 Release" : 
   cd "\OpenSees\Win32\proj\handler"
   $(MAKE) /$(MAKEFLAGS) /F .\handler.mak CFG="handler - Win32 Release" 
   cd "..\openSees"

"handler - Win32 ReleaseCLEAN" : 
   cd "\OpenSees\Win32\proj\handler"
   $(MAKE) /$(MAKEFLAGS) /F .\handler.mak CFG="handler - Win32 Release" RECURSE=1 CLEAN 
   cd "..\openSees"

!ELSEIF  "$(CFG)" == "openSees - Win32 Debug"

"handler - Win32 Debug" : 
   cd "\OpenSees\Win32\proj\handler"
   $(MAKE) /$(MAKEFLAGS) /F .\handler.mak CFG="handler - Win32 Debug" 
   cd "..\openSees"

"handler - Win32 DebugCLEAN" : 
   cd "\OpenSees\Win32\proj\handler"
   $(MAKE) /$(MAKEFLAGS) /F .\handler.mak CFG="handler - Win32 Debug" RECURSE=1 CLEAN 
   cd "..\openSees"

!ENDIF 

!IF  "$(CFG)" == "openSees - Win32 Release"

"material - Win32 Release" : 
   cd "\OpenSees\Win32\proj\material"
   $(MAKE) /$(MAKEFLAGS) /F .\material.mak CFG="material - Win32 Release" 
   cd "..\openSees"

"material - Win32 ReleaseCLEAN" : 
   cd "\OpenSees\Win32\proj\material"
   $(MAKE) /$(MAKEFLAGS) /F .\material.mak CFG="material - Win32 Release" RECURSE=1 CLEAN 
   cd "..\openSees"

!ELSEIF  "$(CFG)" == "openSees - Win32 Debug"

"material - Win32 Debug" : 
   cd "\OpenSees\Win32\proj\material"
   $(MAKE) /$(MAKEFLAGS) /F .\material.mak CFG="material - Win32 Debug" 
   cd "..\openSees"

"material - Win32 DebugCLEAN" : 
   cd "\OpenSees\Win32\proj\material"
   $(MAKE) /$(MAKEFLAGS) /F .\material.mak CFG="material - Win32 Debug" RECURSE=1 CLEAN 
   cd "..\openSees"

!ENDIF 

!IF  "$(CFG)" == "openSees - Win32 Release"

"matrix - Win32 Release" : 
   cd "\OpenSees\Win32\proj\matrix"
   $(MAKE) /$(MAKEFLAGS) /F .\matrix.mak CFG="matrix - Win32 Release" 
   cd "..\openSees"

"matrix - Win32 ReleaseCLEAN" : 
   cd "\OpenSees\Win32\proj\matrix"
   $(MAKE) /$(MAKEFLAGS) /F .\matrix.mak CFG="matrix - Win32 Release" RECURSE=1 CLEAN 
   cd "..\openSees"

!ELSEIF  "$(CFG)" == "openSees - Win32 Debug"

"matrix - Win32 Debug" : 
   cd "\OpenSees\Win32\proj\matrix"
   $(MAKE) /$(MAKEFLAGS) /F .\matrix.mak CFG="matrix - Win32 Debug" 
   cd "..\openSees"

"matrix - Win32 DebugCLEAN" : 
   cd "\OpenSees\Win32\proj\matrix"
   $(MAKE) /$(MAKEFLAGS) /F .\matrix.mak CFG="matrix - Win32 Debug" RECURSE=1 CLEAN 
   cd "..\openSees"

!ENDIF 

!IF  "$(CFG)" == "openSees - Win32 Release"

"modelbuilder - Win32 Release" : 
   cd "\OpenSees\Win32\proj\modelbuilder"
   $(MAKE) /$(MAKEFLAGS) /F .\modelbuilder.mak CFG="modelbuilder - Win32 Release" 
   cd "..\openSees"

"modelbuilder - Win32 ReleaseCLEAN" : 
   cd "\OpenSees\Win32\proj\modelbuilder"
   $(MAKE) /$(MAKEFLAGS) /F .\modelbuilder.mak CFG="modelbuilder - Win32 Release" RECURSE=1 CLEAN 
   cd "..\openSees"

!ELSEIF  "$(CFG)" == "openSees - Win32 Debug"

"modelbuilder - Win32 Debug" : 
   cd "\OpenSees\Win32\proj\modelbuilder"
   $(MAKE) /$(MAKEFLAGS) /F .\modelbuilder.mak CFG="modelbuilder - Win32 Debug" 
   cd "..\openSees"

"modelbuilder - Win32 DebugCLEAN" : 
   cd "\OpenSees\Win32\proj\modelbuilder"
   $(MAKE) /$(MAKEFLAGS) /F .\modelbuilder.mak CFG="modelbuilder - Win32 Debug" RECURSE=1 CLEAN 
   cd "..\openSees"

!ENDIF 

!IF  "$(CFG)" == "openSees - Win32 Release"

"nonlinearBeamColumn - Win32 Release" : 
   cd "\OpenSees\Win32\proj\nonlinearBeamColumn"
   $(MAKE) /$(MAKEFLAGS) /F .\nonlinearBeamColumn.mak CFG="nonlinearBeamColumn - Win32 Release" 
   cd "..\openSees"

"nonlinearBeamColumn - Win32 ReleaseCLEAN" : 
   cd "\OpenSees\Win32\proj\nonlinearBeamColumn"
   $(MAKE) /$(MAKEFLAGS) /F .\nonlinearBeamColumn.mak CFG="nonlinearBeamColumn - Win32 Release" RECURSE=1 CLEAN 
   cd "..\openSees"

!ELSEIF  "$(CFG)" == "openSees - Win32 Debug"

"nonlinearBeamColumn - Win32 Debug" : 
   cd "\OpenSees\Win32\proj\nonlinearBeamColumn"
   $(MAKE) /$(MAKEFLAGS) /F .\nonlinearBeamColumn.mak CFG="nonlinearBeamColumn - Win32 Debug" 
   cd "..\openSees"

"nonlinearBeamColumn - Win32 DebugCLEAN" : 
   cd "\OpenSees\Win32\proj\nonlinearBeamColumn"
   $(MAKE) /$(MAKEFLAGS) /F .\nonlinearBeamColumn.mak CFG="nonlinearBeamColumn - Win32 Debug" RECURSE=1 CLEAN 
   cd "..\openSees"

!ENDIF 

!IF  "$(CFG)" == "openSees - Win32 Release"

"recorder - Win32 Release" : 
   cd "\OpenSees\Win32\proj\recorder"
   $(MAKE) /$(MAKEFLAGS) /F .\recorder.mak CFG="recorder - Win32 Release" 
   cd "..\openSees"

"recorder - Win32 ReleaseCLEAN" : 
   cd "\OpenSees\Win32\proj\recorder"
   $(MAKE) /$(MAKEFLAGS) /F .\recorder.mak CFG="recorder - Win32 Release" RECURSE=1 CLEAN 
   cd "..\openSees"

!ELSEIF  "$(CFG)" == "openSees - Win32 Debug"

"recorder - Win32 Debug" : 
   cd "\OpenSees\Win32\proj\recorder"
   $(MAKE) /$(MAKEFLAGS) /F .\recorder.mak CFG="recorder - Win32 Debug" 
   cd "..\openSees"

"recorder - Win32 DebugCLEAN" : 
   cd "\OpenSees\Win32\proj\recorder"
   $(MAKE) /$(MAKEFLAGS) /F .\recorder.mak CFG="recorder - Win32 Debug" RECURSE=1 CLEAN 
   cd "..\openSees"

!ENDIF 

!IF  "$(CFG)" == "openSees - Win32 Release"

"renderer - Win32 Release" : 
   cd "\OpenSees\Win32\proj\renderer"
   $(MAKE) /$(MAKEFLAGS) /F .\renderer.mak CFG="renderer - Win32 Release" 
   cd "..\openSees"

"renderer - Win32 ReleaseCLEAN" : 
   cd "\OpenSees\Win32\proj\renderer"
   $(MAKE) /$(MAKEFLAGS) /F .\renderer.mak CFG="renderer - Win32 Release" RECURSE=1 CLEAN 
   cd "..\openSees"

!ELSEIF  "$(CFG)" == "openSees - Win32 Debug"

"renderer - Win32 Debug" : 
   cd "\OpenSees\Win32\proj\renderer"
   $(MAKE) /$(MAKEFLAGS) /F .\renderer.mak CFG="renderer - Win32 Debug" 
   cd "..\openSees"

"renderer - Win32 DebugCLEAN" : 
   cd "\OpenSees\Win32\proj\renderer"
   $(MAKE) /$(MAKEFLAGS) /F .\renderer.mak CFG="renderer - Win32 Debug" RECURSE=1 CLEAN 
   cd "..\openSees"

!ENDIF 

!IF  "$(CFG)" == "openSees - Win32 Release"

"string - Win32 Release" : 
   cd "\OpenSees\Win32\proj\string"
   $(MAKE) /$(MAKEFLAGS) /F .\string.mak CFG="string - Win32 Release" 
   cd "..\openSees"

"string - Win32 ReleaseCLEAN" : 
   cd "\OpenSees\Win32\proj\string"
   $(MAKE) /$(MAKEFLAGS) /F .\string.mak CFG="string - Win32 Release" RECURSE=1 CLEAN 
   cd "..\openSees"

!ELSEIF  "$(CFG)" == "openSees - Win32 Debug"

"string - Win32 Debug" : 
   cd "\OpenSees\Win32\proj\string"
   $(MAKE) /$(MAKEFLAGS) /F .\string.mak CFG="string - Win32 Debug" 
   cd "..\openSees"

"string - Win32 DebugCLEAN" : 
   cd "\OpenSees\Win32\proj\string"
   $(MAKE) /$(MAKEFLAGS) /F .\string.mak CFG="string - Win32 Debug" RECURSE=1 CLEAN 
   cd "..\openSees"

!ENDIF 

!IF  "$(CFG)" == "openSees - Win32 Release"

"superLU - Win32 Release" : 
   cd "\OpenSees\Win32\proj\superLU"
   $(MAKE) /$(MAKEFLAGS) /F .\superLU.mak CFG="superLU - Win32 Release" 
   cd "..\openSees"

"superLU - Win32 ReleaseCLEAN" : 
   cd "\OpenSees\Win32\proj\superLU"
   $(MAKE) /$(MAKEFLAGS) /F .\superLU.mak CFG="superLU - Win32 Release" RECURSE=1 CLEAN 
   cd "..\openSees"

!ELSEIF  "$(CFG)" == "openSees - Win32 Debug"

"superLU - Win32 Debug" : 
   cd "\OpenSees\Win32\proj\superLU"
   $(MAKE) /$(MAKEFLAGS) /F .\superLU.mak CFG="superLU - Win32 Debug" 
   cd "..\openSees"

"superLU - Win32 DebugCLEAN" : 
   cd "\OpenSees\Win32\proj\superLU"
   $(MAKE) /$(MAKEFLAGS) /F .\superLU.mak CFG="superLU - Win32 Debug" RECURSE=1 CLEAN 
   cd "..\openSees"

!ENDIF 

!IF  "$(CFG)" == "openSees - Win32 Release"

"system - Win32 Release" : 
   cd "\OpenSees\Win32\proj\system"
   $(MAKE) /$(MAKEFLAGS) /F .\system.mak CFG="system - Win32 Release" 
   cd "..\openSees"

"system - Win32 ReleaseCLEAN" : 
   cd "\OpenSees\Win32\proj\system"
   $(MAKE) /$(MAKEFLAGS) /F .\system.mak CFG="system - Win32 Release" RECURSE=1 CLEAN 
   cd "..\openSees"

!ELSEIF  "$(CFG)" == "openSees - Win32 Debug"

"system - Win32 Debug" : 
   cd "\OpenSees\Win32\proj\system"
   $(MAKE) /$(MAKEFLAGS) /F .\system.mak CFG="system - Win32 Debug" 
   cd "..\openSees"

"system - Win32 DebugCLEAN" : 
   cd "\OpenSees\Win32\proj\system"
   $(MAKE) /$(MAKEFLAGS) /F .\system.mak CFG="system - Win32 Debug" RECURSE=1 CLEAN 
   cd "..\openSees"

!ENDIF 

!IF  "$(CFG)" == "openSees - Win32 Release"

"tagged - Win32 Release" : 
   cd "\OpenSees\Win32\proj\tagged"
   $(MAKE) /$(MAKEFLAGS) /F .\tagged.mak CFG="tagged - Win32 Release" 
   cd "..\openSees"

"tagged - Win32 ReleaseCLEAN" : 
   cd "\OpenSees\Win32\proj\tagged"
   $(MAKE) /$(MAKEFLAGS) /F .\tagged.mak CFG="tagged - Win32 Release" RECURSE=1 CLEAN 
   cd "..\openSees"

!ELSEIF  "$(CFG)" == "openSees - Win32 Debug"

"tagged - Win32 Debug" : 
   cd "\OpenSees\Win32\proj\tagged"
   $(MAKE) /$(MAKEFLAGS) /F .\tagged.mak CFG="tagged - Win32 Debug" 
   cd "..\openSees"

"tagged - Win32 DebugCLEAN" : 
   cd "\OpenSees\Win32\proj\tagged"
   $(MAKE) /$(MAKEFLAGS) /F .\tagged.mak CFG="tagged - Win32 Debug" RECURSE=1 CLEAN 
   cd "..\openSees"

!ENDIF 

!IF  "$(CFG)" == "openSees - Win32 Release"

"utility - Win32 Release" : 
   cd "\OpenSees\Win32\proj\utility"
   $(MAKE) /$(MAKEFLAGS) /F .\utility.mak CFG="utility - Win32 Release" 
   cd "..\openSees"

"utility - Win32 ReleaseCLEAN" : 
   cd "\OpenSees\Win32\proj\utility"
   $(MAKE) /$(MAKEFLAGS) /F .\utility.mak CFG="utility - Win32 Release" RECURSE=1 CLEAN 
   cd "..\openSees"

!ELSEIF  "$(CFG)" == "openSees - Win32 Debug"

"utility - Win32 Debug" : 
   cd "\OpenSees\Win32\proj\utility"
   $(MAKE) /$(MAKEFLAGS) /F .\utility.mak CFG="utility - Win32 Debug" 
   cd "..\openSees"

"utility - Win32 DebugCLEAN" : 
   cd "\OpenSees\Win32\proj\utility"
   $(MAKE) /$(MAKEFLAGS) /F .\utility.mak CFG="utility - Win32 Debug" RECURSE=1 CLEAN 
   cd "..\openSees"

!ENDIF 


!ENDIF 

