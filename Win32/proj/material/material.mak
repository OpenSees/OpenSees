# Microsoft Developer Studio Generated NMAKE File, Based on material.dsp
!IF "$(CFG)" == ""
CFG=material - Win32 Debug
!MESSAGE No configuration specified. Defaulting to material - Win32 Debug.
!ENDIF 

!IF "$(CFG)" != "material - Win32 Release" && "$(CFG)" != "material - Win32 Debug"
!MESSAGE Invalid configuration "$(CFG)" specified.
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "material.mak" CFG="material - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "material - Win32 Release" (based on "Win32 (x86) Static Library")
!MESSAGE "material - Win32 Debug" (based on "Win32 (x86) Static Library")
!MESSAGE 
!ERROR An invalid configuration is specified.
!ENDIF 

!IF "$(OS)" == "Windows_NT"
NULL=
!ELSE 
NULL=nul
!ENDIF 

!IF  "$(CFG)" == "material - Win32 Release"

OUTDIR=.\Release
INTDIR=.\Release
# Begin Custom Macros
OutDir=.\Release
# End Custom Macros

ALL : "$(OUTDIR)\material.lib"


CLEAN :
	-@erase "$(INTDIR)\Concrete01.obj"
	-@erase "$(INTDIR)\ElasticMaterialModel.obj"
	-@erase "$(INTDIR)\ElasticPPMaterialModel.obj"
	-@erase "$(INTDIR)\MaterialModel.obj"
	-@erase "$(INTDIR)\ParallelMaterialModel.obj"
	-@erase "$(INTDIR)\Steel01.obj"
	-@erase "$(INTDIR)\TclModelBuilderMaterialCommand.obj"
	-@erase "$(INTDIR)\UniaxialMaterialModel.obj"
	-@erase "$(INTDIR)\vc60.idb"
	-@erase "$(OUTDIR)\material.lib"

"$(OUTDIR)" :
    if not exist "$(OUTDIR)/$(NULL)" mkdir "$(OUTDIR)"

CPP=cl.exe
CPP_PROJ=/nologo /ML /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_MBCS" /D "_LIB" /Fp"$(INTDIR)\material.pch" /YX /Fo"$(INTDIR)\\" /Fd"$(INTDIR)\\" /FD /c 

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
BSC32_FLAGS=/nologo /o"$(OUTDIR)\material.bsc" 
BSC32_SBRS= \
	
LIB32=link.exe -lib
LIB32_FLAGS=/nologo /out:"$(OUTDIR)\material.lib" 
LIB32_OBJS= \
	"$(INTDIR)\Concrete01.obj" \
	"$(INTDIR)\ElasticMaterialModel.obj" \
	"$(INTDIR)\ElasticPPMaterialModel.obj" \
	"$(INTDIR)\MaterialModel.obj" \
	"$(INTDIR)\ParallelMaterialModel.obj" \
	"$(INTDIR)\Steel01.obj" \
	"$(INTDIR)\TclModelBuilderMaterialCommand.obj" \
	"$(INTDIR)\UniaxialMaterialModel.obj"

"$(OUTDIR)\material.lib" : "$(OUTDIR)" $(DEF_FILE) $(LIB32_OBJS)
    $(LIB32) @<<
  $(LIB32_FLAGS) $(DEF_FLAGS) $(LIB32_OBJS)
<<

!ELSEIF  "$(CFG)" == "material - Win32 Debug"

OUTDIR=.\..\..\lib
INTDIR=.\..\..\obj\material
# Begin Custom Macros
OutDir=.\..\..\lib
# End Custom Macros

ALL : "$(OUTDIR)\material.lib"


CLEAN :
	-@erase "$(INTDIR)\Concrete01.obj"
	-@erase "$(INTDIR)\ElasticMaterialModel.obj"
	-@erase "$(INTDIR)\ElasticPPMaterialModel.obj"
	-@erase "$(INTDIR)\MaterialModel.obj"
	-@erase "$(INTDIR)\ParallelMaterialModel.obj"
	-@erase "$(INTDIR)\Steel01.obj"
	-@erase "$(INTDIR)\TclModelBuilderMaterialCommand.obj"
	-@erase "$(INTDIR)\UniaxialMaterialModel.obj"
	-@erase "$(INTDIR)\vc60.idb"
	-@erase "$(OUTDIR)\material.lib"

"$(OUTDIR)" :
    if not exist "$(OUTDIR)/$(NULL)" mkdir "$(OUTDIR)"

"$(INTDIR)" :
    if not exist "$(INTDIR)/$(NULL)" mkdir "$(INTDIR)"

CPP=cl.exe
CPP_PROJ=/nologo /MLd /W3 /GX /Od /I "c:\Program Files\tcl\include" /I "..\..\src\element" /I "..\..\src\actor\channel" /I "..\..\src\actor\objectBroker" /I "..\..\src\matrix" /I "..\..\src" /I "..\..\src\actor\actor" /I "..\..\src\tagged" /I "..\..\src\modelbuilder" /I "..\..\src\domain\component" /I "..\..\src\material" /I "..\..\src\modelbuilder\tcl" /D "WIN32" /D "_DEBUG" /D "_MBCS" /D "_LIB" /Fp"$(INTDIR)\material.pch" /YX /Fo"$(INTDIR)\\" /Fd"$(INTDIR)\\" /FD /GZ /c 

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
BSC32_FLAGS=/nologo /o"$(OUTDIR)\material.bsc" 
BSC32_SBRS= \
	
LIB32=link.exe -lib
LIB32_FLAGS=/nologo /out:"$(OUTDIR)\material.lib" 
LIB32_OBJS= \
	"$(INTDIR)\Concrete01.obj" \
	"$(INTDIR)\ElasticMaterialModel.obj" \
	"$(INTDIR)\ElasticPPMaterialModel.obj" \
	"$(INTDIR)\MaterialModel.obj" \
	"$(INTDIR)\ParallelMaterialModel.obj" \
	"$(INTDIR)\Steel01.obj" \
	"$(INTDIR)\TclModelBuilderMaterialCommand.obj" \
	"$(INTDIR)\UniaxialMaterialModel.obj"

"$(OUTDIR)\material.lib" : "$(OUTDIR)" $(DEF_FILE) $(LIB32_OBJS)
    $(LIB32) @<<
  $(LIB32_FLAGS) $(DEF_FLAGS) $(LIB32_OBJS)
<<

!ENDIF 


!IF "$(NO_EXTERNAL_DEPS)" != "1"
!IF EXISTS("material.dep")
!INCLUDE "material.dep"
!ELSE 
!MESSAGE Warning: cannot find "material.dep"
!ENDIF 
!ENDIF 


!IF "$(CFG)" == "material - Win32 Release" || "$(CFG)" == "material - Win32 Debug"
SOURCE=..\..\Src\material\Concrete01.cpp

"$(INTDIR)\Concrete01.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\Src\material\ElasticMaterialModel.cpp

"$(INTDIR)\ElasticMaterialModel.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\Src\material\ElasticPPMaterialModel.cpp

"$(INTDIR)\ElasticPPMaterialModel.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\Src\material\MaterialModel.cpp

"$(INTDIR)\MaterialModel.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\Src\material\ParallelMaterialModel.cpp

"$(INTDIR)\ParallelMaterialModel.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\Src\material\Steel01.cpp

"$(INTDIR)\Steel01.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\Src\material\TclModelBuilderMaterialCommand.cpp

"$(INTDIR)\TclModelBuilderMaterialCommand.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\Src\material\UniaxialMaterialModel.cpp

"$(INTDIR)\UniaxialMaterialModel.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)



!ENDIF 

