# Microsoft Developer Studio Generated NMAKE File, Based on convergence.dsp
!IF "$(CFG)" == ""
CFG=convergence - Win32 Debug
!MESSAGE No configuration specified. Defaulting to convergence - Win32 Debug.
!ENDIF 

!IF "$(CFG)" != "convergence - Win32 Release" && "$(CFG)" != "convergence - Win32 Debug"
!MESSAGE Invalid configuration "$(CFG)" specified.
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "convergence.mak" CFG="convergence - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "convergence - Win32 Release" (based on "Win32 (x86) Static Library")
!MESSAGE "convergence - Win32 Debug" (based on "Win32 (x86) Static Library")
!MESSAGE 
!ERROR An invalid configuration is specified.
!ENDIF 

!IF "$(OS)" == "Windows_NT"
NULL=
!ELSE 
NULL=nul
!ENDIF 

!IF  "$(CFG)" == "convergence - Win32 Release"

OUTDIR=.\Release
INTDIR=.\Release
# Begin Custom Macros
OutDir=.\Release
# End Custom Macros

ALL : "$(OUTDIR)\convergence.lib"


CLEAN :
	-@erase "$(INTDIR)\ConvergenceTest.obj"
	-@erase "$(INTDIR)\CTestEnergyIncr.obj"
	-@erase "$(INTDIR)\CTestNormDispIncr.obj"
	-@erase "$(INTDIR)\CTestNormUnbalance.obj"
	-@erase "$(INTDIR)\vc60.idb"
	-@erase "$(OUTDIR)\convergence.lib"

"$(OUTDIR)" :
    if not exist "$(OUTDIR)/$(NULL)" mkdir "$(OUTDIR)"

CPP=cl.exe
CPP_PROJ=/nologo /ML /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_MBCS" /D "_LIB" /Fp"$(INTDIR)\convergence.pch" /YX /Fo"$(INTDIR)\\" /Fd"$(INTDIR)\\" /FD /c 

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
BSC32_FLAGS=/nologo /o"$(OUTDIR)\convergence.bsc" 
BSC32_SBRS= \
	
LIB32=link.exe -lib
LIB32_FLAGS=/nologo /out:"$(OUTDIR)\convergence.lib" 
LIB32_OBJS= \
	"$(INTDIR)\ConvergenceTest.obj" \
	"$(INTDIR)\CTestEnergyIncr.obj" \
	"$(INTDIR)\CTestNormDispIncr.obj" \
	"$(INTDIR)\CTestNormUnbalance.obj"

"$(OUTDIR)\convergence.lib" : "$(OUTDIR)" $(DEF_FILE) $(LIB32_OBJS)
    $(LIB32) @<<
  $(LIB32_FLAGS) $(DEF_FLAGS) $(LIB32_OBJS)
<<

!ELSEIF  "$(CFG)" == "convergence - Win32 Debug"

OUTDIR=.\..\..\lib
INTDIR=.\..\..\obj\convergence
# Begin Custom Macros
OutDir=.\..\..\lib
# End Custom Macros

ALL : "$(OUTDIR)\convergence.lib"


CLEAN :
	-@erase "$(INTDIR)\ConvergenceTest.obj"
	-@erase "$(INTDIR)\CTestEnergyIncr.obj"
	-@erase "$(INTDIR)\CTestNormDispIncr.obj"
	-@erase "$(INTDIR)\CTestNormUnbalance.obj"
	-@erase "$(INTDIR)\vc60.idb"
	-@erase "$(OUTDIR)\convergence.lib"

"$(OUTDIR)" :
    if not exist "$(OUTDIR)/$(NULL)" mkdir "$(OUTDIR)"

"$(INTDIR)" :
    if not exist "$(INTDIR)/$(NULL)" mkdir "$(INTDIR)"

CPP=cl.exe
CPP_PROJ=/nologo /MLd /W3 /GX /Od /I "..\..\src\actor\objectBroker" /I "..\..\src\system_of_eqn" /I "..\..\src\system_of_eqn\linearSOE" /I "..\..\src\analysis\algorithm" /I "..\..\src\analysis\algorithm\equiSolnAlgo" /I "..\..\src\actor\channel" /I "..\..\src\matrix" /I "..\..\src" /I "..\..\src\actor\actor" /I "..\..\src\convergenceTest" /D "WIN32" /D "_DEBUG" /D "_MBCS" /D "_LIB" /Fp"$(INTDIR)\convergence.pch" /YX /Fo"$(INTDIR)\\" /Fd"$(INTDIR)\\" /FD /GZ /c 

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
BSC32_FLAGS=/nologo /o"$(OUTDIR)\convergence.bsc" 
BSC32_SBRS= \
	
LIB32=link.exe -lib
LIB32_FLAGS=/nologo /out:"$(OUTDIR)\convergence.lib" 
LIB32_OBJS= \
	"$(INTDIR)\ConvergenceTest.obj" \
	"$(INTDIR)\CTestEnergyIncr.obj" \
	"$(INTDIR)\CTestNormDispIncr.obj" \
	"$(INTDIR)\CTestNormUnbalance.obj"

"$(OUTDIR)\convergence.lib" : "$(OUTDIR)" $(DEF_FILE) $(LIB32_OBJS)
    $(LIB32) @<<
  $(LIB32_FLAGS) $(DEF_FLAGS) $(LIB32_OBJS)
<<

!ENDIF 


!IF "$(NO_EXTERNAL_DEPS)" != "1"
!IF EXISTS("convergence.dep")
!INCLUDE "convergence.dep"
!ELSE 
!MESSAGE Warning: cannot find "convergence.dep"
!ENDIF 
!ENDIF 


!IF "$(CFG)" == "convergence - Win32 Release" || "$(CFG)" == "convergence - Win32 Debug"
SOURCE=..\..\Src\convergenceTest\ConvergenceTest.cpp

"$(INTDIR)\ConvergenceTest.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\Src\convergenceTest\CTestEnergyIncr.cpp

"$(INTDIR)\CTestEnergyIncr.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\Src\convergenceTest\CTestNormDispIncr.cpp

"$(INTDIR)\CTestNormDispIncr.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\Src\convergenceTest\CTestNormUnbalance.cpp

"$(INTDIR)\CTestNormUnbalance.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)



!ENDIF 

