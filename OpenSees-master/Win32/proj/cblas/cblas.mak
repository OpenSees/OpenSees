# Microsoft Developer Studio Generated NMAKE File, Based on cblas.dsp
!IF "$(CFG)" == ""
CFG=cblas - Win32 Debug
!MESSAGE No configuration specified. Defaulting to cblas - Win32 Debug.
!ENDIF 

!IF "$(CFG)" != "cblas - Win32 Release" && "$(CFG)" != "cblas - Win32 Debug"
!MESSAGE Invalid configuration "$(CFG)" specified.
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "cblas.mak" CFG="cblas - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "cblas - Win32 Release" (based on "Win32 (x86) Static Library")
!MESSAGE "cblas - Win32 Debug" (based on "Win32 (x86) Static Library")
!MESSAGE 
!ERROR An invalid configuration is specified.
!ENDIF 

!IF "$(OS)" == "Windows_NT"
NULL=
!ELSE 
NULL=nul
!ENDIF 

!IF  "$(CFG)" == "cblas - Win32 Release"

OUTDIR=.\Release
INTDIR=.\Release
# Begin Custom Macros
OutDir=.\Release
# End Custom Macros

ALL : "$(OUTDIR)\cblas.lib"


CLEAN :
	-@erase "$(INTDIR)\dasum.obj"
	-@erase "$(INTDIR)\daxpy.obj"
	-@erase "$(INTDIR)\dcabs1.obj"
	-@erase "$(INTDIR)\dcopy.obj"
	-@erase "$(INTDIR)\ddot.obj"
	-@erase "$(INTDIR)\dgemv.obj"
	-@erase "$(INTDIR)\dger.obj"
	-@erase "$(INTDIR)\dmyblas2.obj"
	-@erase "$(INTDIR)\dnrm2.obj"
	-@erase "$(INTDIR)\drot.obj"
	-@erase "$(INTDIR)\dscal.obj"
	-@erase "$(INTDIR)\dsymv.obj"
	-@erase "$(INTDIR)\dsyr2.obj"
	-@erase "$(INTDIR)\dtrsv.obj"
	-@erase "$(INTDIR)\dzasum.obj"
	-@erase "$(INTDIR)\dznrm2.obj"
	-@erase "$(INTDIR)\idamax.obj"
	-@erase "$(INTDIR)\vc60.idb"
	-@erase "$(OUTDIR)\cblas.lib"

"$(OUTDIR)" :
    if not exist "$(OUTDIR)/$(NULL)" mkdir "$(OUTDIR)"

CPP=cl.exe
CPP_PROJ=/nologo /ML /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_MBCS" /D "_LIB" /Fp"$(INTDIR)\cblas.pch" /YX /Fo"$(INTDIR)\\" /Fd"$(INTDIR)\\" /FD /c 

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
BSC32_FLAGS=/nologo /o"$(OUTDIR)\cblas.bsc" 
BSC32_SBRS= \
	
LIB32=link.exe -lib
LIB32_FLAGS=/nologo /out:"$(OUTDIR)\cblas.lib" 
LIB32_OBJS= \
	"$(INTDIR)\dasum.obj" \
	"$(INTDIR)\daxpy.obj" \
	"$(INTDIR)\dcabs1.obj" \
	"$(INTDIR)\dcopy.obj" \
	"$(INTDIR)\ddot.obj" \
	"$(INTDIR)\dgemv.obj" \
	"$(INTDIR)\dger.obj" \
	"$(INTDIR)\dmyblas2.obj" \
	"$(INTDIR)\dnrm2.obj" \
	"$(INTDIR)\drot.obj" \
	"$(INTDIR)\dscal.obj" \
	"$(INTDIR)\dsymv.obj" \
	"$(INTDIR)\dsyr2.obj" \
	"$(INTDIR)\dtrsv.obj" \
	"$(INTDIR)\dzasum.obj" \
	"$(INTDIR)\dznrm2.obj" \
	"$(INTDIR)\idamax.obj"

"$(OUTDIR)\cblas.lib" : "$(OUTDIR)" $(DEF_FILE) $(LIB32_OBJS)
    $(LIB32) @<<
  $(LIB32_FLAGS) $(DEF_FLAGS) $(LIB32_OBJS)
<<

!ELSEIF  "$(CFG)" == "cblas - Win32 Debug"

OUTDIR=.\..\..\lib
INTDIR=.\..\..\obj\cblas
# Begin Custom Macros
OutDir=.\..\..\lib
# End Custom Macros

ALL : "$(OUTDIR)\cblas.lib"


CLEAN :
	-@erase "$(INTDIR)\dasum.obj"
	-@erase "$(INTDIR)\daxpy.obj"
	-@erase "$(INTDIR)\dcabs1.obj"
	-@erase "$(INTDIR)\dcopy.obj"
	-@erase "$(INTDIR)\ddot.obj"
	-@erase "$(INTDIR)\dgemv.obj"
	-@erase "$(INTDIR)\dger.obj"
	-@erase "$(INTDIR)\dmyblas2.obj"
	-@erase "$(INTDIR)\dnrm2.obj"
	-@erase "$(INTDIR)\drot.obj"
	-@erase "$(INTDIR)\dscal.obj"
	-@erase "$(INTDIR)\dsymv.obj"
	-@erase "$(INTDIR)\dsyr2.obj"
	-@erase "$(INTDIR)\dtrsv.obj"
	-@erase "$(INTDIR)\dzasum.obj"
	-@erase "$(INTDIR)\dznrm2.obj"
	-@erase "$(INTDIR)\idamax.obj"
	-@erase "$(INTDIR)\vc60.idb"
	-@erase "$(OUTDIR)\cblas.lib"

"$(OUTDIR)" :
    if not exist "$(OUTDIR)/$(NULL)" mkdir "$(OUTDIR)"

"$(INTDIR)" :
    if not exist "$(INTDIR)/$(NULL)" mkdir "$(INTDIR)"

CPP=cl.exe
CPP_PROJ=/nologo /MLd /W3 /GX /Od /I "..\..\cblas" /D "WIN32" /D "_DEBUG" /D "_MBCS" /D "_LIB" /Fp"$(INTDIR)\cblas.pch" /YX /Fo"$(INTDIR)\\" /Fd"$(INTDIR)\\" /FD /GZ /c 

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
BSC32_FLAGS=/nologo /o"$(OUTDIR)\cblas.bsc" 
BSC32_SBRS= \
	
LIB32=link.exe -lib
LIB32_FLAGS=/nologo /out:"$(OUTDIR)\cblas.lib" 
LIB32_OBJS= \
	"$(INTDIR)\dasum.obj" \
	"$(INTDIR)\daxpy.obj" \
	"$(INTDIR)\dcabs1.obj" \
	"$(INTDIR)\dcopy.obj" \
	"$(INTDIR)\ddot.obj" \
	"$(INTDIR)\dgemv.obj" \
	"$(INTDIR)\dger.obj" \
	"$(INTDIR)\dmyblas2.obj" \
	"$(INTDIR)\dnrm2.obj" \
	"$(INTDIR)\drot.obj" \
	"$(INTDIR)\dscal.obj" \
	"$(INTDIR)\dsymv.obj" \
	"$(INTDIR)\dsyr2.obj" \
	"$(INTDIR)\dtrsv.obj" \
	"$(INTDIR)\dzasum.obj" \
	"$(INTDIR)\dznrm2.obj" \
	"$(INTDIR)\idamax.obj"

"$(OUTDIR)\cblas.lib" : "$(OUTDIR)" $(DEF_FILE) $(LIB32_OBJS)
    $(LIB32) @<<
  $(LIB32_FLAGS) $(DEF_FLAGS) $(LIB32_OBJS)
<<

!ENDIF 


!IF "$(NO_EXTERNAL_DEPS)" != "1"
!IF EXISTS("cblas.dep")
!INCLUDE "cblas.dep"
!ELSE 
!MESSAGE Warning: cannot find "cblas.dep"
!ENDIF 
!ENDIF 


!IF "$(CFG)" == "cblas - Win32 Release" || "$(CFG)" == "cblas - Win32 Debug"
SOURCE=..\..\CBLAS\dasum.c

"$(INTDIR)\dasum.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\CBLAS\daxpy.c

"$(INTDIR)\daxpy.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\CBLAS\dcabs1.c

"$(INTDIR)\dcabs1.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\CBLAS\dcopy.c

"$(INTDIR)\dcopy.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\CBLAS\ddot.c

"$(INTDIR)\ddot.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\CBLAS\dgemv.c

"$(INTDIR)\dgemv.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\CBLAS\dger.c

"$(INTDIR)\dger.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\CBLAS\dmyblas2.c

"$(INTDIR)\dmyblas2.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\CBLAS\dnrm2.c

"$(INTDIR)\dnrm2.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\CBLAS\drot.c

"$(INTDIR)\drot.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\CBLAS\dscal.c

"$(INTDIR)\dscal.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\CBLAS\dsymv.c

"$(INTDIR)\dsymv.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\CBLAS\dsyr2.c

"$(INTDIR)\dsyr2.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\CBLAS\dtrsv.c

"$(INTDIR)\dtrsv.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\CBLAS\dzasum.c

"$(INTDIR)\dzasum.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\CBLAS\dznrm2.c

"$(INTDIR)\dznrm2.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\CBLAS\idamax.c

"$(INTDIR)\idamax.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)



!ENDIF 

