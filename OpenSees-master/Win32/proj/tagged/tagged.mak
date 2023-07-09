# Microsoft Developer Studio Generated NMAKE File, Based on tagged.dsp
!IF "$(CFG)" == ""
CFG=tagged - Win32 Debug
!MESSAGE No configuration specified. Defaulting to tagged - Win32 Debug.
!ENDIF 

!IF "$(CFG)" != "tagged - Win32 Release" && "$(CFG)" != "tagged - Win32 Debug"
!MESSAGE Invalid configuration "$(CFG)" specified.
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "tagged.mak" CFG="tagged - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "tagged - Win32 Release" (based on "Win32 (x86) Static Library")
!MESSAGE "tagged - Win32 Debug" (based on "Win32 (x86) Static Library")
!MESSAGE 
!ERROR An invalid configuration is specified.
!ENDIF 

!IF "$(OS)" == "Windows_NT"
NULL=
!ELSE 
NULL=nul
!ENDIF 

!IF  "$(CFG)" == "tagged - Win32 Release"

OUTDIR=.\Release
INTDIR=.\Release
# Begin Custom Macros
OutDir=.\Release
# End Custom Macros

ALL : "$(OUTDIR)\tagged.lib"


CLEAN :
	-@erase "$(INTDIR)\ArrayOfTaggedObjects.obj"
	-@erase "$(INTDIR)\ArrayOfTaggedObjectsIter.obj"
	-@erase "$(INTDIR)\TaggedObject.obj"
	-@erase "$(INTDIR)\vc60.idb"
	-@erase "$(OUTDIR)\tagged.lib"

"$(OUTDIR)" :
    if not exist "$(OUTDIR)/$(NULL)" mkdir "$(OUTDIR)"

CPP=cl.exe
CPP_PROJ=/nologo /ML /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_MBCS" /D "_LIB" /Fp"$(INTDIR)\tagged.pch" /YX /Fo"$(INTDIR)\\" /Fd"$(INTDIR)\\" /FD /c 

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
BSC32_FLAGS=/nologo /o"$(OUTDIR)\tagged.bsc" 
BSC32_SBRS= \
	
LIB32=link.exe -lib
LIB32_FLAGS=/nologo /out:"$(OUTDIR)\tagged.lib" 
LIB32_OBJS= \
	"$(INTDIR)\TaggedObject.obj" \
	"$(INTDIR)\ArrayOfTaggedObjects.obj" \
	"$(INTDIR)\ArrayOfTaggedObjectsIter.obj"

"$(OUTDIR)\tagged.lib" : "$(OUTDIR)" $(DEF_FILE) $(LIB32_OBJS)
    $(LIB32) @<<
  $(LIB32_FLAGS) $(DEF_FLAGS) $(LIB32_OBJS)
<<

!ELSEIF  "$(CFG)" == "tagged - Win32 Debug"

OUTDIR=.\..\..\lib
INTDIR=.\..\..\obj\tagged
# Begin Custom Macros
OutDir=.\..\..\lib
# End Custom Macros

ALL : "$(OUTDIR)\tagged.lib"


CLEAN :
	-@erase "$(INTDIR)\ArrayOfTaggedObjects.obj"
	-@erase "$(INTDIR)\ArrayOfTaggedObjectsIter.obj"
	-@erase "$(INTDIR)\TaggedObject.obj"
	-@erase "$(INTDIR)\vc60.idb"
	-@erase "$(OUTDIR)\tagged.lib"

"$(OUTDIR)" :
    if not exist "$(OUTDIR)/$(NULL)" mkdir "$(OUTDIR)"

"$(INTDIR)" :
    if not exist "$(INTDIR)/$(NULL)" mkdir "$(INTDIR)"

CPP=cl.exe
CPP_PROJ=/nologo /MLd /W3 /GX /Od /I "..\..\src\matrix" /I "..\..\src\domain\domain" /I "..\..\src\tagged\storage" /I "..\..\src" /I "..\..\src\tagged" /D "WIN32" /D "_DEBUG" /D "_MBCS" /D "_LIB" /Fp"$(INTDIR)\tagged.pch" /YX /Fo"$(INTDIR)\\" /Fd"$(INTDIR)\\" /FD /GZ /c 

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
BSC32_FLAGS=/nologo /o"$(OUTDIR)\tagged.bsc" 
BSC32_SBRS= \
	
LIB32=link.exe -lib
LIB32_FLAGS=/nologo /out:"$(OUTDIR)\tagged.lib" 
LIB32_OBJS= \
	"$(INTDIR)\TaggedObject.obj" \
	"$(INTDIR)\ArrayOfTaggedObjects.obj" \
	"$(INTDIR)\ArrayOfTaggedObjectsIter.obj"

"$(OUTDIR)\tagged.lib" : "$(OUTDIR)" $(DEF_FILE) $(LIB32_OBJS)
    $(LIB32) @<<
  $(LIB32_FLAGS) $(DEF_FLAGS) $(LIB32_OBJS)
<<

!ENDIF 


!IF "$(NO_EXTERNAL_DEPS)" != "1"
!IF EXISTS("tagged.dep")
!INCLUDE "tagged.dep"
!ELSE 
!MESSAGE Warning: cannot find "tagged.dep"
!ENDIF 
!ENDIF 


!IF "$(CFG)" == "tagged - Win32 Release" || "$(CFG)" == "tagged - Win32 Debug"
SOURCE=..\..\Src\tagged\TaggedObject.cpp

"$(INTDIR)\TaggedObject.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\Src\tagged\storage\ArrayOfTaggedObjects.cpp

"$(INTDIR)\ArrayOfTaggedObjects.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\Src\tagged\storage\ArrayOfTaggedObjectsIter.cpp

"$(INTDIR)\ArrayOfTaggedObjectsIter.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)



!ENDIF 

