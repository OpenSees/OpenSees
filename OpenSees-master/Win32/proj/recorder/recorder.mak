# Microsoft Developer Studio Generated NMAKE File, Based on recorder.dsp
!IF "$(CFG)" == ""
CFG=recorder - Win32 Debug
!MESSAGE No configuration specified. Defaulting to recorder - Win32 Debug.
!ENDIF 

!IF "$(CFG)" != "recorder - Win32 Release" && "$(CFG)" != "recorder - Win32 Debug"
!MESSAGE Invalid configuration "$(CFG)" specified.
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "recorder.mak" CFG="recorder - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "recorder - Win32 Release" (based on "Win32 (x86) Static Library")
!MESSAGE "recorder - Win32 Debug" (based on "Win32 (x86) Static Library")
!MESSAGE 
!ERROR An invalid configuration is specified.
!ENDIF 

!IF "$(OS)" == "Windows_NT"
NULL=
!ELSE 
NULL=nul
!ENDIF 

!IF  "$(CFG)" == "recorder - Win32 Release"

OUTDIR=.\Release
INTDIR=.\Release
# Begin Custom Macros
OutDir=.\Release
# End Custom Macros

ALL : "$(OUTDIR)\recorder.lib"


CLEAN :
	-@erase "$(INTDIR)\DatastoreRecorder.obj"
	-@erase "$(INTDIR)\ElementRecorder.obj"
	-@erase "$(INTDIR)\FileNodeDispRecorder.obj"
	-@erase "$(INTDIR)\MaxNodeDispRecorder.obj"
	-@erase "$(INTDIR)\NodeIncrDispRecorder.obj"
	-@erase "$(INTDIR)\TclRecorderCommands.obj"
	-@erase "$(INTDIR)\vc60.idb"
	-@erase "$(OUTDIR)\recorder.lib"

"$(OUTDIR)" :
    if not exist "$(OUTDIR)/$(NULL)" mkdir "$(OUTDIR)"

CPP=cl.exe
CPP_PROJ=/nologo /ML /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_MBCS" /D "_LIB" /Fp"$(INTDIR)\recorder.pch" /YX /Fo"$(INTDIR)\\" /Fd"$(INTDIR)\\" /FD /c 

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
BSC32_FLAGS=/nologo /o"$(OUTDIR)\recorder.bsc" 
BSC32_SBRS= \
	
LIB32=link.exe -lib
LIB32_FLAGS=/nologo /out:"$(OUTDIR)\recorder.lib" 
LIB32_OBJS= \
	"$(INTDIR)\DatastoreRecorder.obj" \
	"$(INTDIR)\ElementRecorder.obj" \
	"$(INTDIR)\FileNodeDispRecorder.obj" \
	"$(INTDIR)\MaxNodeDispRecorder.obj" \
	"$(INTDIR)\NodeIncrDispRecorder.obj" \
	"$(INTDIR)\TclRecorderCommands.obj"

"$(OUTDIR)\recorder.lib" : "$(OUTDIR)" $(DEF_FILE) $(LIB32_OBJS)
    $(LIB32) @<<
  $(LIB32_FLAGS) $(DEF_FLAGS) $(LIB32_OBJS)
<<

!ELSEIF  "$(CFG)" == "recorder - Win32 Debug"

OUTDIR=.\..\..\lib
INTDIR=.\..\..\obj\recorder
# Begin Custom Macros
OutDir=.\..\..\lib
# End Custom Macros

ALL : "$(OUTDIR)\recorder.lib"


CLEAN :
	-@erase "$(INTDIR)\DatastoreRecorder.obj"
	-@erase "$(INTDIR)\ElementRecorder.obj"
	-@erase "$(INTDIR)\FileNodeDispRecorder.obj"
	-@erase "$(INTDIR)\MaxNodeDispRecorder.obj"
	-@erase "$(INTDIR)\NodeIncrDispRecorder.obj"
	-@erase "$(INTDIR)\TclRecorderCommands.obj"
	-@erase "$(INTDIR)\vc60.idb"
	-@erase "$(OUTDIR)\recorder.lib"

"$(OUTDIR)" :
    if not exist "$(OUTDIR)/$(NULL)" mkdir "$(OUTDIR)"

"$(INTDIR)" :
    if not exist "$(INTDIR)/$(NULL)" mkdir "$(INTDIR)"

CPP=cl.exe
CPP_PROJ=/nologo /MLd /W3 /GX /Od /I "c:\Program Files\tcl\include" /I "..\..\src\actor\channel" /I "..\..\src\modelbuilder" /I "..\..\src\database" /I "..\..\src\tagged" /I "..\..\src\actor\actor" /I "..\..\src\domain\component" /I "..\..\src\element" /I "..\..\src\matrix" /I "..\..\src\domain\domain" /I "..\..\src" /I "..\..\src\recorder" /I "..\..\src\domain\node" /D "WIN32" /D "_DEBUG" /D "_MBCS" /D "_LIB" /Fp"$(INTDIR)\recorder.pch" /YX /Fo"$(INTDIR)\\" /Fd"$(INTDIR)\\" /FD /GZ /c 

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
BSC32_FLAGS=/nologo /o"$(OUTDIR)\recorder.bsc" 
BSC32_SBRS= \
	
LIB32=link.exe -lib
LIB32_FLAGS=/nologo /out:"$(OUTDIR)\recorder.lib" 
LIB32_OBJS= \
	"$(INTDIR)\DatastoreRecorder.obj" \
	"$(INTDIR)\ElementRecorder.obj" \
	"$(INTDIR)\FileNodeDispRecorder.obj" \
	"$(INTDIR)\MaxNodeDispRecorder.obj" \
	"$(INTDIR)\NodeIncrDispRecorder.obj" \
	"$(INTDIR)\TclRecorderCommands.obj"

"$(OUTDIR)\recorder.lib" : "$(OUTDIR)" $(DEF_FILE) $(LIB32_OBJS)
    $(LIB32) @<<
  $(LIB32_FLAGS) $(DEF_FLAGS) $(LIB32_OBJS)
<<

!ENDIF 


!IF "$(NO_EXTERNAL_DEPS)" != "1"
!IF EXISTS("recorder.dep")
!INCLUDE "recorder.dep"
!ELSE 
!MESSAGE Warning: cannot find "recorder.dep"
!ENDIF 
!ENDIF 


!IF "$(CFG)" == "recorder - Win32 Release" || "$(CFG)" == "recorder - Win32 Debug"
SOURCE=..\..\Src\recorder\DatastoreRecorder.cpp

"$(INTDIR)\DatastoreRecorder.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\Src\recorder\ElementRecorder.cpp

"$(INTDIR)\ElementRecorder.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\Src\recorder\FileNodeDispRecorder.cpp

"$(INTDIR)\FileNodeDispRecorder.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\Src\recorder\MaxNodeDispRecorder.cpp

"$(INTDIR)\MaxNodeDispRecorder.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\Src\recorder\NodeIncrDispRecorder.cpp

"$(INTDIR)\NodeIncrDispRecorder.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\Src\recorder\TclRecorderCommands.cpp

"$(INTDIR)\TclRecorderCommands.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)



!ENDIF 

