# Microsoft Developer Studio Generated NMAKE File, Based on nonlinearBeamColumn.dsp
!IF "$(CFG)" == ""
CFG=nonlinearBeamColumn - Win32 Debug
!MESSAGE No configuration specified. Defaulting to nonlinearBeamColumn - Win32 Debug.
!ENDIF 

!IF "$(CFG)" != "nonlinearBeamColumn - Win32 Release" && "$(CFG)" != "nonlinearBeamColumn - Win32 Debug"
!MESSAGE Invalid configuration "$(CFG)" specified.
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "nonlinearBeamColumn.mak" CFG="nonlinearBeamColumn - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "nonlinearBeamColumn - Win32 Release" (based on "Win32 (x86) Static Library")
!MESSAGE "nonlinearBeamColumn - Win32 Debug" (based on "Win32 (x86) Static Library")
!MESSAGE 
!ERROR An invalid configuration is specified.
!ENDIF 

!IF "$(OS)" == "Windows_NT"
NULL=
!ELSE 
NULL=nul
!ENDIF 

!IF  "$(CFG)" == "nonlinearBeamColumn - Win32 Release"

OUTDIR=.\Release
INTDIR=.\Release
# Begin Custom Macros
OutDir=.\Release
# End Custom Macros

ALL : "$(OUTDIR)\nonlinearBeamColumn.lib"


CLEAN :
	-@erase "$(INTDIR)\Cell.obj"
	-@erase "$(INTDIR)\CircPatch.obj"
	-@erase "$(INTDIR)\CircReinfLayer.obj"
	-@erase "$(INTDIR)\Fiber.obj"
	-@erase "$(INTDIR)\FiberSection2d.obj"
	-@erase "$(INTDIR)\FiberSection3d.obj"
	-@erase "$(INTDIR)\FiberSectionRepr.obj"
	-@erase "$(INTDIR)\GaussLobattoQuadRule1d01.obj"
	-@erase "$(INTDIR)\GaussQuadRule1d.obj"
	-@erase "$(INTDIR)\GaussQuadRule1d01.obj"
	-@erase "$(INTDIR)\MatrixUtil.obj"
	-@erase "$(INTDIR)\mySection.obj"
	-@erase "$(INTDIR)\NLBeamColumn2d.obj"
	-@erase "$(INTDIR)\NLBeamColumn3d.obj"
	-@erase "$(INTDIR)\Patch.obj"
	-@erase "$(INTDIR)\QuadCell.obj"
	-@erase "$(INTDIR)\QuadPatch.obj"
	-@erase "$(INTDIR)\QuadRule.obj"
	-@erase "$(INTDIR)\QuadRule1d.obj"
	-@erase "$(INTDIR)\QuadRule1d01.obj"
	-@erase "$(INTDIR)\ReinfBar.obj"
	-@erase "$(INTDIR)\ReinfLayer.obj"
	-@erase "$(INTDIR)\Section2d.obj"
	-@erase "$(INTDIR)\Section3d.obj"
	-@erase "$(INTDIR)\SectionRepres.obj"
	-@erase "$(INTDIR)\StraightReinfLayer.obj"
	-@erase "$(INTDIR)\TclElmtBuilder.obj"
	-@erase "$(INTDIR)\UniaxialFiber.obj"
	-@erase "$(INTDIR)\UniaxialFiber2d.obj"
	-@erase "$(INTDIR)\UniaxialFiber3d.obj"
	-@erase "$(INTDIR)\vc60.idb"
	-@erase "$(OUTDIR)\nonlinearBeamColumn.lib"

"$(OUTDIR)" :
    if not exist "$(OUTDIR)/$(NULL)" mkdir "$(OUTDIR)"

CPP=cl.exe
CPP_PROJ=/nologo /ML /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_MBCS" /D "_LIB" /Fp"$(INTDIR)\nonlinearBeamColumn.pch" /YX /Fo"$(INTDIR)\\" /Fd"$(INTDIR)\\" /FD /c 

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
BSC32_FLAGS=/nologo /o"$(OUTDIR)\nonlinearBeamColumn.bsc" 
BSC32_SBRS= \
	
LIB32=link.exe -lib
LIB32_FLAGS=/nologo /out:"$(OUTDIR)\nonlinearBeamColumn.lib" 
LIB32_OBJS= \
	"$(INTDIR)\NLBeamColumn2d.obj" \
	"$(INTDIR)\NLBeamColumn3d.obj" \
	"$(INTDIR)\Fiber.obj" \
	"$(INTDIR)\UniaxialFiber.obj" \
	"$(INTDIR)\UniaxialFiber2d.obj" \
	"$(INTDIR)\UniaxialFiber3d.obj" \
	"$(INTDIR)\FiberSection2d.obj" \
	"$(INTDIR)\FiberSection3d.obj" \
	"$(INTDIR)\mySection.obj" \
	"$(INTDIR)\Section2d.obj" \
	"$(INTDIR)\Section3d.obj" \
	"$(INTDIR)\MatrixUtil.obj" \
	"$(INTDIR)\GaussLobattoQuadRule1d01.obj" \
	"$(INTDIR)\GaussQuadRule1d.obj" \
	"$(INTDIR)\GaussQuadRule1d01.obj" \
	"$(INTDIR)\QuadRule.obj" \
	"$(INTDIR)\QuadRule1d.obj" \
	"$(INTDIR)\QuadRule1d01.obj" \
	"$(INTDIR)\Cell.obj" \
	"$(INTDIR)\QuadCell.obj" \
	"$(INTDIR)\CircPatch.obj" \
	"$(INTDIR)\Patch.obj" \
	"$(INTDIR)\QuadPatch.obj" \
	"$(INTDIR)\ReinfBar.obj" \
	"$(INTDIR)\CircReinfLayer.obj" \
	"$(INTDIR)\ReinfLayer.obj" \
	"$(INTDIR)\StraightReinfLayer.obj" \
	"$(INTDIR)\FiberSectionRepr.obj" \
	"$(INTDIR)\SectionRepres.obj" \
	"$(INTDIR)\TclElmtBuilder.obj"

"$(OUTDIR)\nonlinearBeamColumn.lib" : "$(OUTDIR)" $(DEF_FILE) $(LIB32_OBJS)
    $(LIB32) @<<
  $(LIB32_FLAGS) $(DEF_FLAGS) $(LIB32_OBJS)
<<

!ELSEIF  "$(CFG)" == "nonlinearBeamColumn - Win32 Debug"

OUTDIR=.\..\..\lib
INTDIR=.\..\..\obj\nonlinearBeamColumn
# Begin Custom Macros
OutDir=.\..\..\lib
# End Custom Macros

ALL : "$(OUTDIR)\nonlinearBeamColumn.lib"


CLEAN :
	-@erase "$(INTDIR)\Cell.obj"
	-@erase "$(INTDIR)\CircPatch.obj"
	-@erase "$(INTDIR)\CircReinfLayer.obj"
	-@erase "$(INTDIR)\Fiber.obj"
	-@erase "$(INTDIR)\FiberSection2d.obj"
	-@erase "$(INTDIR)\FiberSection3d.obj"
	-@erase "$(INTDIR)\FiberSectionRepr.obj"
	-@erase "$(INTDIR)\GaussLobattoQuadRule1d01.obj"
	-@erase "$(INTDIR)\GaussQuadRule1d.obj"
	-@erase "$(INTDIR)\GaussQuadRule1d01.obj"
	-@erase "$(INTDIR)\MatrixUtil.obj"
	-@erase "$(INTDIR)\mySection.obj"
	-@erase "$(INTDIR)\NLBeamColumn2d.obj"
	-@erase "$(INTDIR)\NLBeamColumn3d.obj"
	-@erase "$(INTDIR)\Patch.obj"
	-@erase "$(INTDIR)\QuadCell.obj"
	-@erase "$(INTDIR)\QuadPatch.obj"
	-@erase "$(INTDIR)\QuadRule.obj"
	-@erase "$(INTDIR)\QuadRule1d.obj"
	-@erase "$(INTDIR)\QuadRule1d01.obj"
	-@erase "$(INTDIR)\ReinfBar.obj"
	-@erase "$(INTDIR)\ReinfLayer.obj"
	-@erase "$(INTDIR)\Section2d.obj"
	-@erase "$(INTDIR)\Section3d.obj"
	-@erase "$(INTDIR)\SectionRepres.obj"
	-@erase "$(INTDIR)\StraightReinfLayer.obj"
	-@erase "$(INTDIR)\TclElmtBuilder.obj"
	-@erase "$(INTDIR)\UniaxialFiber.obj"
	-@erase "$(INTDIR)\UniaxialFiber2d.obj"
	-@erase "$(INTDIR)\UniaxialFiber3d.obj"
	-@erase "$(INTDIR)\vc60.idb"
	-@erase "$(OUTDIR)\nonlinearBeamColumn.lib"

"$(OUTDIR)" :
    if not exist "$(OUTDIR)/$(NULL)" mkdir "$(OUTDIR)"

"$(INTDIR)" :
    if not exist "$(INTDIR)/$(NULL)" mkdir "$(INTDIR)"

CPP=cl.exe
CPP_PROJ=/nologo /MLd /W3 /GX /Od /I "c:\Program Files\tcl\include" /I "..\..\src\numberer" /I "..\..\src\element\nonlinearBeamColumn\section" /I "..\..\src\modelbuilder" /I "..\..\src\matrix" /I "..\..\src\renderer" /I "..\..\src\modelbuilder\tcl" /I "..\..\src\actor\objectBroker" /I "..\..\src\tagged\storage" /I "..\..\src\domain\node" /I "..\..\src\element\nonlinearBeamColumn\tcl\repres\section" /I "..\..\src\element\nonlinearBeamColumn\tcl\repres\patch" /I "..\..\src\element\nonlinearBeamColumn\tcl\repres\reinfBar" /I "..\..\src\element\nonlinearBeamColumn\tcl\repres\reinfLayer" /I "..\..\src\element\nonlinearBeamColumn\tcl\repres\cell" /I "..\..\src\domain\domain" /I "..\..\src\element\nonlinearBeamColumn" /I "..\..\src\element\nonlinearBEamCOlumn\quadrule" /I "..\..\src\element\nonlinearBeamColumn\matrixutil" /I "..\..\src" /I "..\..\src\element" /I "..\..\src\element\nonlinearBeamColumn\fiber" /I "..\..\src\domain\component" /I "..\..\src\tagged" /I "..\..\src\actor\actor" /I "..\..\src\actor\channel" /I "..\..\src\material" /I "..\..\src\element\nonlinearBeamColumn\element" /D "WIN32" /D "_DEBUG" /D "_MBCS" /D "_LIB" /Fp"$(INTDIR)\nonlinearBeamColumn.pch" /YX\
 /Fo"$(INTDIR)\\" /Fd"$(INTDIR)\\" /FD /GZ /c 

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
BSC32_FLAGS=/nologo /o"$(OUTDIR)\nonlinearBeamColumn.bsc" 
BSC32_SBRS= \
	
LIB32=link.exe -lib
LIB32_FLAGS=/nologo /out:"$(OUTDIR)\nonlinearBeamColumn.lib" 
LIB32_OBJS= \
	"$(INTDIR)\NLBeamColumn2d.obj" \
	"$(INTDIR)\NLBeamColumn3d.obj" \
	"$(INTDIR)\Fiber.obj" \
	"$(INTDIR)\UniaxialFiber.obj" \
	"$(INTDIR)\UniaxialFiber2d.obj" \
	"$(INTDIR)\UniaxialFiber3d.obj" \
	"$(INTDIR)\FiberSection2d.obj" \
	"$(INTDIR)\FiberSection3d.obj" \
	"$(INTDIR)\mySection.obj" \
	"$(INTDIR)\Section2d.obj" \
	"$(INTDIR)\Section3d.obj" \
	"$(INTDIR)\MatrixUtil.obj" \
	"$(INTDIR)\GaussLobattoQuadRule1d01.obj" \
	"$(INTDIR)\GaussQuadRule1d.obj" \
	"$(INTDIR)\GaussQuadRule1d01.obj" \
	"$(INTDIR)\QuadRule.obj" \
	"$(INTDIR)\QuadRule1d.obj" \
	"$(INTDIR)\QuadRule1d01.obj" \
	"$(INTDIR)\Cell.obj" \
	"$(INTDIR)\QuadCell.obj" \
	"$(INTDIR)\CircPatch.obj" \
	"$(INTDIR)\Patch.obj" \
	"$(INTDIR)\QuadPatch.obj" \
	"$(INTDIR)\ReinfBar.obj" \
	"$(INTDIR)\CircReinfLayer.obj" \
	"$(INTDIR)\ReinfLayer.obj" \
	"$(INTDIR)\StraightReinfLayer.obj" \
	"$(INTDIR)\FiberSectionRepr.obj" \
	"$(INTDIR)\SectionRepres.obj" \
	"$(INTDIR)\TclElmtBuilder.obj"

"$(OUTDIR)\nonlinearBeamColumn.lib" : "$(OUTDIR)" $(DEF_FILE) $(LIB32_OBJS)
    $(LIB32) @<<
  $(LIB32_FLAGS) $(DEF_FLAGS) $(LIB32_OBJS)
<<

!ENDIF 


!IF "$(NO_EXTERNAL_DEPS)" != "1"
!IF EXISTS("nonlinearBeamColumn.dep")
!INCLUDE "nonlinearBeamColumn.dep"
!ELSE 
!MESSAGE Warning: cannot find "nonlinearBeamColumn.dep"
!ENDIF 
!ENDIF 


!IF "$(CFG)" == "nonlinearBeamColumn - Win32 Release" || "$(CFG)" == "nonlinearBeamColumn - Win32 Debug"
SOURCE=..\..\Src\element\nonlinearBeamColumn\element\NLBeamColumn2d.cpp

"$(INTDIR)\NLBeamColumn2d.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\Src\element\nonlinearBeamColumn\element\NLBeamColumn3d.cpp

"$(INTDIR)\NLBeamColumn3d.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\Src\element\nonlinearBeamColumn\fiber\Fiber.cpp

"$(INTDIR)\Fiber.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\Src\element\nonlinearBeamColumn\fiber\UniaxialFiber.cpp

"$(INTDIR)\UniaxialFiber.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\Src\element\nonlinearBeamColumn\fiber\UniaxialFiber2d.cpp

"$(INTDIR)\UniaxialFiber2d.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\Src\element\nonlinearBeamColumn\fiber\UniaxialFiber3d.cpp

"$(INTDIR)\UniaxialFiber3d.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\Src\element\nonlinearBeamColumn\section\FiberSection2d.cpp

"$(INTDIR)\FiberSection2d.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\Src\element\nonlinearBeamColumn\section\FiberSection3d.cpp

"$(INTDIR)\FiberSection3d.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\Src\element\nonlinearBeamColumn\section\mySection.cpp

"$(INTDIR)\mySection.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\Src\element\nonlinearBeamColumn\section\Section2d.cpp

"$(INTDIR)\Section2d.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\Src\element\nonlinearBeamColumn\section\Section3d.cpp

"$(INTDIR)\Section3d.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\Src\element\nonlinearBeamColumn\matrixutil\MatrixUtil.cpp

"$(INTDIR)\MatrixUtil.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\Src\element\nonlinearBeamColumn\quadrule\GaussLobattoQuadRule1d01.cpp

"$(INTDIR)\GaussLobattoQuadRule1d01.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\Src\element\nonlinearBeamColumn\quadrule\GaussQuadRule1d.cpp

"$(INTDIR)\GaussQuadRule1d.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\Src\element\nonlinearBeamColumn\quadrule\GaussQuadRule1d01.cpp

"$(INTDIR)\GaussQuadRule1d01.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\Src\element\nonlinearBeamColumn\quadrule\QuadRule.cpp

"$(INTDIR)\QuadRule.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\Src\element\nonlinearBeamColumn\quadrule\QuadRule1d.cpp

"$(INTDIR)\QuadRule1d.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\Src\element\nonlinearBeamColumn\quadrule\QuadRule1d01.cpp

"$(INTDIR)\QuadRule1d01.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\Src\element\nonlinearBeamColumn\tcl\repres\cell\Cell.cpp

"$(INTDIR)\Cell.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\Src\element\nonlinearBeamColumn\tcl\repres\cell\QuadCell.cpp

"$(INTDIR)\QuadCell.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\Src\element\nonlinearBeamColumn\tcl\repres\patch\CircPatch.cpp

"$(INTDIR)\CircPatch.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\Src\element\nonlinearBeamColumn\tcl\repres\patch\Patch.cpp

"$(INTDIR)\Patch.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\Src\element\nonlinearBeamColumn\tcl\repres\patch\QuadPatch.cpp

"$(INTDIR)\QuadPatch.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\Src\element\nonlinearBeamColumn\tcl\repres\reinfBar\ReinfBar.cpp

"$(INTDIR)\ReinfBar.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\Src\element\nonlinearBeamColumn\tcl\repres\reinfLayer\CircReinfLayer.cpp

"$(INTDIR)\CircReinfLayer.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\Src\element\nonlinearBeamColumn\tcl\repres\reinfLayer\ReinfLayer.cpp

"$(INTDIR)\ReinfLayer.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\Src\element\nonlinearBeamColumn\tcl\repres\reinfLayer\StraightReinfLayer.cpp

"$(INTDIR)\StraightReinfLayer.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\Src\element\nonlinearBeamColumn\tcl\repres\section\FiberSectionRepr.cpp

"$(INTDIR)\FiberSectionRepr.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\Src\element\nonlinearBeamColumn\tcl\repres\section\SectionRepres.cpp

"$(INTDIR)\SectionRepres.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\Src\element\nonlinearBeamColumn\tcl\TclElmtBuilder.cpp

"$(INTDIR)\TclElmtBuilder.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)



!ENDIF 

