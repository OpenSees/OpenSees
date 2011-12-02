# Microsoft Developer Studio Generated NMAKE File, Based on system.dsp
!IF "$(CFG)" == ""
CFG=system - Win32 Debug
!MESSAGE No configuration specified. Defaulting to system - Win32 Debug.
!ENDIF 

!IF "$(CFG)" != "system - Win32 Release" && "$(CFG)" != "system - Win32 Debug"
!MESSAGE Invalid configuration "$(CFG)" specified.
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "system.mak" CFG="system - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "system - Win32 Release" (based on "Win32 (x86) Static Library")
!MESSAGE "system - Win32 Debug" (based on "Win32 (x86) Static Library")
!MESSAGE 
!ERROR An invalid configuration is specified.
!ENDIF 

!IF "$(OS)" == "Windows_NT"
NULL=
!ELSE 
NULL=nul
!ENDIF 

!IF  "$(CFG)" == "system - Win32 Release"

OUTDIR=.\Release
INTDIR=.\Release
# Begin Custom Macros
OutDir=.\Release
# End Custom Macros

ALL : "$(OUTDIR)\system.lib"


CLEAN :
	-@erase "$(INTDIR)\BandGenLinLapackSolver.obj"
	-@erase "$(INTDIR)\BandGenLinSOE.obj"
	-@erase "$(INTDIR)\BandGenLinSolver.obj"
	-@erase "$(INTDIR)\BandSPDLinLapackSolver.obj"
	-@erase "$(INTDIR)\BandSPDLinSOE.obj"
	-@erase "$(INTDIR)\BandSPDLinSolver.obj"
	-@erase "$(INTDIR)\FullGenLinSOE.obj"
	-@erase "$(INTDIR)\FullGenLinSolver.obj"
	-@erase "$(INTDIR)\LinearSOE.obj"
	-@erase "$(INTDIR)\LinearSOESolver.obj"
	-@erase "$(INTDIR)\ProfileSPDLinDirectSolver.obj"
	-@erase "$(INTDIR)\ProfileSPDLinSOE.obj"
	-@erase "$(INTDIR)\ProfileSPDLinSolver.obj"
	-@erase "$(INTDIR)\Solver.obj"
	-@erase "$(INTDIR)\SparseGenColLinSOE.obj"
	-@erase "$(INTDIR)\SparseGenColLinSolver.obj"
	-@erase "$(INTDIR)\SuperLU.obj"
	-@erase "$(INTDIR)\SymSparseLinSOE.obj"
	-@erase "$(INTDIR)\SymSparseLinSolver.obj"
	-@erase "$(INTDIR)\SystemOfEqn.obj"
	-@erase "$(INTDIR)\UmfpackGenLinSOE.obj"
	-@erase "$(INTDIR)\UmfpackGenLinSolver.obj"
	-@erase "$(INTDIR)\vc60.idb"
	-@erase "$(OUTDIR)\system.lib"

"$(OUTDIR)" :
    if not exist "$(OUTDIR)/$(NULL)" mkdir "$(OUTDIR)"

CPP=cl.exe
CPP_PROJ=/nologo /ML /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_MBCS" /D "_LIB" /Fp"$(INTDIR)\system.pch" /YX /Fo"$(INTDIR)\\" /Fd"$(INTDIR)\\" /FD /c 

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
BSC32_FLAGS=/nologo /o"$(OUTDIR)\system.bsc" 
BSC32_SBRS= \
	
LIB32=link.exe -lib
LIB32_FLAGS=/nologo /out:"$(OUTDIR)\system.lib" 
LIB32_OBJS= \
	"$(INTDIR)\LinearSOE.obj" \
	"$(INTDIR)\LinearSOESolver.obj" \
	"$(INTDIR)\Solver.obj" \
	"$(INTDIR)\SystemOfEqn.obj" \
	"$(INTDIR)\ProfileSPDLinDirectSolver.obj" \
	"$(INTDIR)\ProfileSPDLinSOE.obj" \
	"$(INTDIR)\ProfileSPDLinSolver.obj" \
	"$(INTDIR)\SymSparseLinSOE.obj" \
	"$(INTDIR)\SymSparseLinSolver.obj" \
	"$(INTDIR)\BandGenLinLapackSolver.obj" \
	"$(INTDIR)\BandGenLinSOE.obj" \
	"$(INTDIR)\BandGenLinSolver.obj" \
	"$(INTDIR)\BandSPDLinLapackSolver.obj" \
	"$(INTDIR)\BandSPDLinSOE.obj" \
	"$(INTDIR)\BandSPDLinSolver.obj" \
	"$(INTDIR)\FullGenLinSOE.obj" \
	"$(INTDIR)\FullGenLinSolver.obj" \
	"$(INTDIR)\UmfpackGenLinSOE.obj" \
	"$(INTDIR)\UmfpackGenLinSolver.obj" \
	"$(INTDIR)\SparseGenColLinSOE.obj" \
	"$(INTDIR)\SparseGenColLinSolver.obj" \
	"$(INTDIR)\SuperLU.obj"

"$(OUTDIR)\system.lib" : "$(OUTDIR)" $(DEF_FILE) $(LIB32_OBJS)
    $(LIB32) @<<
  $(LIB32_FLAGS) $(DEF_FLAGS) $(LIB32_OBJS)
<<

!ELSEIF  "$(CFG)" == "system - Win32 Debug"

OUTDIR=.\..\..\lib
INTDIR=.\..\..\obj\system
# Begin Custom Macros
OutDir=.\..\..\lib
# End Custom Macros

ALL : "$(OUTDIR)\system.lib"


CLEAN :
	-@erase "$(INTDIR)\BandGenLinLapackSolver.obj"
	-@erase "$(INTDIR)\BandGenLinSOE.obj"
	-@erase "$(INTDIR)\BandGenLinSolver.obj"
	-@erase "$(INTDIR)\BandSPDLinLapackSolver.obj"
	-@erase "$(INTDIR)\BandSPDLinSOE.obj"
	-@erase "$(INTDIR)\BandSPDLinSolver.obj"
	-@erase "$(INTDIR)\FullGenLinSOE.obj"
	-@erase "$(INTDIR)\FullGenLinSolver.obj"
	-@erase "$(INTDIR)\LinearSOE.obj"
	-@erase "$(INTDIR)\LinearSOESolver.obj"
	-@erase "$(INTDIR)\ProfileSPDLinDirectSolver.obj"
	-@erase "$(INTDIR)\ProfileSPDLinSOE.obj"
	-@erase "$(INTDIR)\ProfileSPDLinSolver.obj"
	-@erase "$(INTDIR)\Solver.obj"
	-@erase "$(INTDIR)\SparseGenColLinSOE.obj"
	-@erase "$(INTDIR)\SparseGenColLinSolver.obj"
	-@erase "$(INTDIR)\SuperLU.obj"
	-@erase "$(INTDIR)\SymSparseLinSOE.obj"
	-@erase "$(INTDIR)\SymSparseLinSolver.obj"
	-@erase "$(INTDIR)\SystemOfEqn.obj"
	-@erase "$(INTDIR)\UmfpackGenLinSOE.obj"
	-@erase "$(INTDIR)\UmfpackGenLinSolver.obj"
	-@erase "$(INTDIR)\vc60.idb"
	-@erase "$(OUTDIR)\system.lib"

"$(OUTDIR)" :
    if not exist "$(OUTDIR)/$(NULL)" mkdir "$(OUTDIR)"

"$(INTDIR)" :
    if not exist "$(INTDIR)/$(NULL)" mkdir "$(INTDIR)"

CPP=cl.exe
CPP_PROJ=/nologo /MLd /W3 /GX /Od /I "..\..\symSparse" /I "..\..\src\analysis\model\simple" /I "..\..\src\system_of_eqn\linearSOE\umfGEN" /I "..\..\src\system_of_eqn\linearSOE\fullGEN" /I "..\..\src\system_of_eqn\linearSOE\sparseGEN" /I "..\..\src\system_of_eqn\linearSOE\bandSPD" /I "..\..\src\system_of_eqn\linearSOE\bandGEN" /I "..\..\src\domain\domain" /I "..\..\src\analysis\model" /I "..\..\src\actor\objectBroker" /I "..\..\src\actor\channel" /I "..\..\src\tagged" /I "..\..\src\graph\graph" /I "..\..\src\system_of_eqn\linearSOE\profileSPD" /I "..\..\src" /I "..\..\src\matrix" /I "..\..\src\actor\actor" /I "..\..\src\system_of_eqn" /I "..\..\src\system_of_eqn\linearSOE" /D "WIN32" /D "_DEBUG" /D "_MBCS" /D "_LIB" /D "_WIN32" /Fp"$(INTDIR)\system.pch" /YX /Fo"$(INTDIR)\\" /Fd"$(INTDIR)\\" /FD /GZ /c 

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
BSC32_FLAGS=/nologo /o"$(OUTDIR)\system.bsc" 
BSC32_SBRS= \
	
LIB32=link.exe -lib
LIB32_FLAGS=/nologo /out:"$(OUTDIR)\system.lib" 
LIB32_OBJS= \
	"$(INTDIR)\LinearSOE.obj" \
	"$(INTDIR)\LinearSOESolver.obj" \
	"$(INTDIR)\Solver.obj" \
	"$(INTDIR)\SystemOfEqn.obj" \
	"$(INTDIR)\ProfileSPDLinDirectSolver.obj" \
	"$(INTDIR)\ProfileSPDLinSOE.obj" \
	"$(INTDIR)\ProfileSPDLinSolver.obj" \
	"$(INTDIR)\SymSparseLinSOE.obj" \
	"$(INTDIR)\SymSparseLinSolver.obj" \
	"$(INTDIR)\BandGenLinLapackSolver.obj" \
	"$(INTDIR)\BandGenLinSOE.obj" \
	"$(INTDIR)\BandGenLinSolver.obj" \
	"$(INTDIR)\BandSPDLinLapackSolver.obj" \
	"$(INTDIR)\BandSPDLinSOE.obj" \
	"$(INTDIR)\BandSPDLinSolver.obj" \
	"$(INTDIR)\FullGenLinSOE.obj" \
	"$(INTDIR)\FullGenLinSolver.obj" \
	"$(INTDIR)\UmfpackGenLinSOE.obj" \
	"$(INTDIR)\UmfpackGenLinSolver.obj" \
	"$(INTDIR)\SparseGenColLinSOE.obj" \
	"$(INTDIR)\SparseGenColLinSolver.obj" \
	"$(INTDIR)\SuperLU.obj"

"$(OUTDIR)\system.lib" : "$(OUTDIR)" $(DEF_FILE) $(LIB32_OBJS)
    $(LIB32) @<<
  $(LIB32_FLAGS) $(DEF_FLAGS) $(LIB32_OBJS)
<<

!ENDIF 


!IF "$(NO_EXTERNAL_DEPS)" != "1"
!IF EXISTS("system.dep")
!INCLUDE "system.dep"
!ELSE 
!MESSAGE Warning: cannot find "system.dep"
!ENDIF 
!ENDIF 


!IF "$(CFG)" == "system - Win32 Release" || "$(CFG)" == "system - Win32 Debug"
SOURCE=..\..\Src\system_of_eqn\linearSOE\LinearSOE.cpp

"$(INTDIR)\LinearSOE.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\Src\system_of_eqn\linearSOE\LinearSOESolver.cpp

"$(INTDIR)\LinearSOESolver.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\Src\system_of_eqn\Solver.cpp

"$(INTDIR)\Solver.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\Src\system_of_eqn\SystemOfEqn.cpp

"$(INTDIR)\SystemOfEqn.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\Src\system_of_eqn\linearSOE\profileSPD\ProfileSPDLinDirectSolver.cpp

"$(INTDIR)\ProfileSPDLinDirectSolver.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\Src\system_of_eqn\linearSOE\profileSPD\ProfileSPDLinSOE.cpp

"$(INTDIR)\ProfileSPDLinSOE.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\Src\system_of_eqn\linearSOE\profileSPD\ProfileSPDLinSolver.cpp

"$(INTDIR)\ProfileSPDLinSolver.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\Src\system_of_eqn\linearSOE\sparseSYM\SymSparseLinSOE.cpp

"$(INTDIR)\SymSparseLinSOE.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\Src\system_of_eqn\linearSOE\sparseSYM\SymSparseLinSolver.cpp

"$(INTDIR)\SymSparseLinSolver.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\src\system_of_eqn\linearSOE\bandGEN\BandGenLinLapackSolver.cpp

"$(INTDIR)\BandGenLinLapackSolver.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\Src\system_of_eqn\linearSOE\bandGEN\BandGenLinSOE.cpp

"$(INTDIR)\BandGenLinSOE.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\Src\system_of_eqn\linearSOE\bandGEN\BandGenLinSolver.cpp

"$(INTDIR)\BandGenLinSolver.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\src\system_of_eqn\linearSOE\bandSPD\BandSPDLinLapackSolver.cpp

"$(INTDIR)\BandSPDLinLapackSolver.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\Src\system_of_eqn\linearSOE\bandSPD\BandSPDLinSOE.cpp

"$(INTDIR)\BandSPDLinSOE.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\Src\system_of_eqn\linearSOE\bandSPD\BandSPDLinSolver.cpp

"$(INTDIR)\BandSPDLinSolver.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\Src\system_of_eqn\linearSOE\fullGEN\FullGenLinSOE.cpp

"$(INTDIR)\FullGenLinSOE.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\Src\system_of_eqn\linearSOE\fullGEN\FullGenLinSolver.cpp

"$(INTDIR)\FullGenLinSolver.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\Src\system_of_eqn\linearSOE\umfGEN\UmfpackGenLinSOE.cpp

"$(INTDIR)\UmfpackGenLinSOE.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\Src\system_of_eqn\linearSOE\umfGEN\UmfpackGenLinSolver.cpp

"$(INTDIR)\UmfpackGenLinSolver.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\Src\system_of_eqn\linearSOE\sparseGEN\SparseGenColLinSOE.cpp

"$(INTDIR)\SparseGenColLinSOE.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\Src\system_of_eqn\linearSOE\sparseGEN\SparseGenColLinSolver.cpp

"$(INTDIR)\SparseGenColLinSolver.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\Src\system_of_eqn\linearSOE\sparseGEN\SuperLU.cpp

"$(INTDIR)\SuperLU.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)



!ENDIF 

