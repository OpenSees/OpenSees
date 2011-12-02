# Microsoft Developer Studio Generated NMAKE File, Based on modelbuilder.dsp
!IF "$(CFG)" == ""
CFG=modelbuilder - Win32 Debug
!MESSAGE No configuration specified. Defaulting to modelbuilder - Win32 Debug.
!ENDIF 

!IF "$(CFG)" != "modelbuilder - Win32 Release" && "$(CFG)" != "modelbuilder - Win32 Debug"
!MESSAGE Invalid configuration "$(CFG)" specified.
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "modelbuilder.mak" CFG="modelbuilder - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "modelbuilder - Win32 Release" (based on "Win32 (x86) Static Library")
!MESSAGE "modelbuilder - Win32 Debug" (based on "Win32 (x86) Static Library")
!MESSAGE 
!ERROR An invalid configuration is specified.
!ENDIF 

!IF "$(OS)" == "Windows_NT"
NULL=
!ELSE 
NULL=nul
!ENDIF 

!IF  "$(CFG)" == "modelbuilder - Win32 Release"

OUTDIR=.\Release
INTDIR=.\Release
# Begin Custom Macros
OutDir=.\Release
# End Custom Macros

ALL : "$(OUTDIR)\modelbuilder.lib"


CLEAN :
	-@erase "$(INTDIR)\ModelBuilder.obj"
	-@erase "$(INTDIR)\PlaneFrame.obj"
	-@erase "$(INTDIR)\Quick2dFrame.obj"
	-@erase "$(INTDIR)\Quick3dFrame.obj"
	-@erase "$(INTDIR)\vc60.idb"
	-@erase "$(OUTDIR)\modelbuilder.lib"

"$(OUTDIR)" :
    if not exist "$(OUTDIR)/$(NULL)" mkdir "$(OUTDIR)"

CPP=cl.exe
CPP_PROJ=/nologo /ML /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_MBCS" /D "_LIB" /Fp"$(INTDIR)\modelbuilder.pch" /YX /Fo"$(INTDIR)\\" /Fd"$(INTDIR)\\" /FD /c 

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
BSC32_FLAGS=/nologo /o"$(OUTDIR)\modelbuilder.bsc" 
BSC32_SBRS= \
	
LIB32=link.exe -lib
LIB32_FLAGS=/nologo /out:"$(OUTDIR)\modelbuilder.lib" 
LIB32_OBJS= \
	"$(INTDIR)\ModelBuilder.obj" \
	"$(INTDIR)\PlaneFrame.obj" \
	"$(INTDIR)\Quick2dFrame.obj" \
	"$(INTDIR)\Quick3dFrame.obj"

"$(OUTDIR)\modelbuilder.lib" : "$(OUTDIR)" $(DEF_FILE) $(LIB32_OBJS)
    $(LIB32) @<<
  $(LIB32_FLAGS) $(DEF_FLAGS) $(LIB32_OBJS)
<<

!ELSEIF  "$(CFG)" == "modelbuilder - Win32 Debug"

OUTDIR=.\..\..\lib
INTDIR=.\..\..\obj\modelbuilder
# Begin Custom Macros
OutDir=.\..\..\lib
# End Custom Macros

ALL : "$(OUTDIR)\modelbuilder.lib"


CLEAN :
	-@erase "$(INTDIR)\ModelBuilder.obj"
	-@erase "$(INTDIR)\PlaneFrame.obj"
	-@erase "$(INTDIR)\Quick2dFrame.obj"
	-@erase "$(INTDIR)\Quick3dFrame.obj"
	-@erase "$(INTDIR)\vc60.idb"
	-@erase "$(OUTDIR)\modelbuilder.lib"

"$(OUTDIR)" :
    if not exist "$(OUTDIR)/$(NULL)" mkdir "$(OUTDIR)"

"$(INTDIR)" :
    if not exist "$(INTDIR)/$(NULL)" mkdir "$(INTDIR)"

CPP=cl.exe
CPP_PROJ=/nologo /MLd /W3 /GX /Od /I "..\..\src\domain\pattern" /I "..\..\src\domain\load" /I "..\..\src\actor\channel" /I "..\..\src\domain\node" /I "..\..\src\domain\constraints" /I "..\..\src\tagged" /I "..\..\src\actor\actor" /I "..\..\src\domain\component" /I "..\..\src\element" /I "..\..\src\element\beam3d" /I "..\..\src" /I "..\..\src\matrix" /I "..\..\src\domain\domain" /I "..\..\src\modelbuilder" /I "..\..\src\element\beam2d" /D "WIN32" /D "_DEBUG" /D "_MBCS" /D "_LIB" /Fp"$(INTDIR)\modelbuilder.pch" /YX /Fo"$(INTDIR)\\" /Fd"$(INTDIR)\\" /FD /GZ /c 

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
BSC32_FLAGS=/nologo /o"$(OUTDIR)\modelbuilder.bsc" 
BSC32_SBRS= \
	
LIB32=link.exe -lib
LIB32_FLAGS=/nologo /out:"$(OUTDIR)\modelbuilder.lib" 
LIB32_OBJS= \
	"$(INTDIR)\ModelBuilder.obj" \
	"$(INTDIR)\PlaneFrame.obj" \
	"$(INTDIR)\Quick2dFrame.obj" \
	"$(INTDIR)\Quick3dFrame.obj"

"$(OUTDIR)\modelbuilder.lib" : "$(OUTDIR)" $(DEF_FILE) $(LIB32_OBJS)
    $(LIB32) @<<
  $(LIB32_FLAGS) $(DEF_FLAGS) $(LIB32_OBJS)
<<

!ENDIF 


!IF "$(NO_EXTERNAL_DEPS)" != "1"
!IF EXISTS("modelbuilder.dep")
!INCLUDE "modelbuilder.dep"
!ELSE 
!MESSAGE Warning: cannot find "modelbuilder.dep"
!ENDIF 
!ENDIF 


!IF "$(CFG)" == "modelbuilder - Win32 Release" || "$(CFG)" == "modelbuilder - Win32 Debug"
SOURCE=..\..\Src\modelbuilder\ModelBuilder.cpp

"$(INTDIR)\ModelBuilder.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\Src\modelbuilder\PlaneFrame.cpp

"$(INTDIR)\PlaneFrame.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\Src\modelbuilder\Quick2dFrame.cpp

"$(INTDIR)\Quick2dFrame.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\Src\modelbuilder\Quick3dFrame.cpp

"$(INTDIR)\Quick3dFrame.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)



!ENDIF 

