# Microsoft Developer Studio Generated NMAKE File, Based on element.dsp
!IF "$(CFG)" == ""
CFG=element - Win32 Debug
!MESSAGE No configuration specified. Defaulting to element - Win32 Debug.
!ENDIF 

!IF "$(CFG)" != "element - Win32 Release" && "$(CFG)" != "element - Win32 Debug"
!MESSAGE Invalid configuration "$(CFG)" specified.
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "element.mak" CFG="element - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "element - Win32 Release" (based on "Win32 (x86) Static Library")
!MESSAGE "element - Win32 Debug" (based on "Win32 (x86) Static Library")
!MESSAGE 
!ERROR An invalid configuration is specified.
!ENDIF 

!IF "$(OS)" == "Windows_NT"
NULL=
!ELSE 
NULL=nul
!ENDIF 

!IF  "$(CFG)" == "element - Win32 Release"

OUTDIR=.\Release
INTDIR=.\Release
# Begin Custom Macros
OutDir=.\Release
# End Custom Macros

ALL : "$(OUTDIR)\element.lib"


CLEAN :
	-@erase "$(INTDIR)\beam2d02.obj"
	-@erase "$(INTDIR)\beam2d03.obj"
	-@erase "$(INTDIR)\beam2d04.obj"
	-@erase "$(INTDIR)\beam3d01.obj"
	-@erase "$(INTDIR)\beam3d02.obj"
	-@erase "$(INTDIR)\ElasticBeam2d.obj"
	-@erase "$(INTDIR)\ElasticBeam3d.obj"
	-@erase "$(INTDIR)\Element.obj"
	-@erase "$(INTDIR)\ElementalLoad.obj"
	-@erase "$(INTDIR)\Information.obj"
	-@erase "$(INTDIR)\Truss.obj"
	-@erase "$(INTDIR)\TrussSection.obj"
	-@erase "$(INTDIR)\vc60.idb"
	-@erase "$(OUTDIR)\element.lib"

"$(OUTDIR)" :
    if not exist "$(OUTDIR)/$(NULL)" mkdir "$(OUTDIR)"

CPP=cl.exe
CPP_PROJ=/nologo /ML /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_MBCS" /D "_LIB" /Fp"$(INTDIR)\element.pch" /YX /Fo"$(INTDIR)\\" /Fd"$(INTDIR)\\" /FD /c 

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
BSC32_FLAGS=/nologo /o"$(OUTDIR)\element.bsc" 
BSC32_SBRS= \
	
LIB32=link.exe -lib
LIB32_FLAGS=/nologo /out:"$(OUTDIR)\element.lib" 
LIB32_OBJS= \
	"$(INTDIR)\Element.obj" \
	"$(INTDIR)\ElementalLoad.obj" \
	"$(INTDIR)\Information.obj" \
	"$(INTDIR)\Truss.obj" \
	"$(INTDIR)\TrussSection.obj" \
	"$(INTDIR)\beam2d02.obj" \
	"$(INTDIR)\beam2d03.obj" \
	"$(INTDIR)\beam2d04.obj" \
	"$(INTDIR)\ElasticBeam2d.obj" \
	"$(INTDIR)\beam3d01.obj" \
	"$(INTDIR)\beam3d02.obj" \
	"$(INTDIR)\ElasticBeam3d.obj"

"$(OUTDIR)\element.lib" : "$(OUTDIR)" $(DEF_FILE) $(LIB32_OBJS)
    $(LIB32) @<<
  $(LIB32_FLAGS) $(DEF_FLAGS) $(LIB32_OBJS)
<<

!ELSEIF  "$(CFG)" == "element - Win32 Debug"

OUTDIR=.\..\..\lib
INTDIR=.\..\..\obj\element
# Begin Custom Macros
OutDir=.\..\..\lib
# End Custom Macros

ALL : "$(OUTDIR)\element.lib"


CLEAN :
	-@erase "$(INTDIR)\beam2d02.obj"
	-@erase "$(INTDIR)\beam2d03.obj"
	-@erase "$(INTDIR)\beam2d04.obj"
	-@erase "$(INTDIR)\beam3d01.obj"
	-@erase "$(INTDIR)\beam3d02.obj"
	-@erase "$(INTDIR)\ElasticBeam2d.obj"
	-@erase "$(INTDIR)\ElasticBeam3d.obj"
	-@erase "$(INTDIR)\Element.obj"
	-@erase "$(INTDIR)\ElementalLoad.obj"
	-@erase "$(INTDIR)\Information.obj"
	-@erase "$(INTDIR)\Truss.obj"
	-@erase "$(INTDIR)\TrussSection.obj"
	-@erase "$(INTDIR)\vc60.idb"
	-@erase "$(OUTDIR)\element.lib"

"$(OUTDIR)" :
    if not exist "$(OUTDIR)/$(NULL)" mkdir "$(OUTDIR)"

"$(INTDIR)" :
    if not exist "$(INTDIR)/$(NULL)" mkdir "$(INTDIR)"

CPP=cl.exe
CPP_PROJ=/nologo /MLd /W3 /GX /Od /I "..\..\src\element" /I "..\..\src\element\truss" /I "..\..\src\element\nonlinearBeamColumn\section" /I "..\..\src\element\beam3d" /I "..\..\src\element\beam2d" /I "..\..\src\material" /I "..\..\src\actor\objectBroker" /I "..\..\src\matrix" /I "..\..\src\domain\load" /I "..\..\src\renderer" /I "..\..\src\actor\channel" /I "..\..\src\domain\node" /I "..\..\src\actor\actor" /I "..\..\src\tagged" /I "..\..\src\domain\component" /I "..\..\src" /I "..\..\src\domain\domain" /D "WIN32" /D "_DEBUG" /D "_MBCS" /D "_LIB" /Fp"$(INTDIR)\element.pch" /YX /Fo"$(INTDIR)\\" /Fd"$(INTDIR)\\" /FD /GZ /c 

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
BSC32_FLAGS=/nologo /o"$(OUTDIR)\element.bsc" 
BSC32_SBRS= \
	
LIB32=link.exe -lib
LIB32_FLAGS=/nologo /out:"$(OUTDIR)\element.lib" 
LIB32_OBJS= \
	"$(INTDIR)\Element.obj" \
	"$(INTDIR)\ElementalLoad.obj" \
	"$(INTDIR)\Information.obj" \
	"$(INTDIR)\Truss.obj" \
	"$(INTDIR)\TrussSection.obj" \
	"$(INTDIR)\beam2d02.obj" \
	"$(INTDIR)\beam2d03.obj" \
	"$(INTDIR)\beam2d04.obj" \
	"$(INTDIR)\ElasticBeam2d.obj" \
	"$(INTDIR)\beam3d01.obj" \
	"$(INTDIR)\beam3d02.obj" \
	"$(INTDIR)\ElasticBeam3d.obj"

"$(OUTDIR)\element.lib" : "$(OUTDIR)" $(DEF_FILE) $(LIB32_OBJS)
    $(LIB32) @<<
  $(LIB32_FLAGS) $(DEF_FLAGS) $(LIB32_OBJS)
<<

!ENDIF 


!IF "$(NO_EXTERNAL_DEPS)" != "1"
!IF EXISTS("element.dep")
!INCLUDE "element.dep"
!ELSE 
!MESSAGE Warning: cannot find "element.dep"
!ENDIF 
!ENDIF 


!IF "$(CFG)" == "element - Win32 Release" || "$(CFG)" == "element - Win32 Debug"
SOURCE=..\..\Src\element\Element.cpp

"$(INTDIR)\Element.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\Src\element\ElementalLoad.cpp

"$(INTDIR)\ElementalLoad.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\src\element\Information.cpp

"$(INTDIR)\Information.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\Src\element\truss\Truss.cpp

"$(INTDIR)\Truss.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\Src\element\truss\TrussSection.cpp

"$(INTDIR)\TrussSection.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\Src\element\beam2d\beam2d02.cpp

"$(INTDIR)\beam2d02.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\Src\element\beam2d\beam2d03.cpp

"$(INTDIR)\beam2d03.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\Src\element\beam2d\beam2d04.cpp

"$(INTDIR)\beam2d04.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\Src\element\beam2d\ElasticBeam2d.cpp

"$(INTDIR)\ElasticBeam2d.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\Src\element\beam3d\beam3d01.cpp

"$(INTDIR)\beam3d01.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\Src\element\beam3d\beam3d02.cpp

"$(INTDIR)\beam3d02.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\Src\element\beam3d\ElasticBeam3d.cpp

"$(INTDIR)\ElasticBeam3d.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)



!ENDIF 

