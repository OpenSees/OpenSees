# Microsoft Developer Studio Generated NMAKE File, Based on domain.dsp
!IF "$(CFG)" == ""
CFG=domain - Win32 Debug
!MESSAGE No configuration specified. Defaulting to domain - Win32 Debug.
!ENDIF 

!IF "$(CFG)" != "domain - Win32 Release" && "$(CFG)" != "domain - Win32 Debug"
!MESSAGE Invalid configuration "$(CFG)" specified.
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "domain.mak" CFG="domain - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "domain - Win32 Release" (based on "Win32 (x86) Static Library")
!MESSAGE "domain - Win32 Debug" (based on "Win32 (x86) Static Library")
!MESSAGE 
!ERROR An invalid configuration is specified.
!ENDIF 

!IF "$(OS)" == "Windows_NT"
NULL=
!ELSE 
NULL=nul
!ENDIF 

!IF  "$(CFG)" == "domain - Win32 Release"

OUTDIR=.\Release
INTDIR=.\Release
# Begin Custom Macros
OutDir=.\Release
# End Custom Macros

ALL : "$(OUTDIR)\domain.lib"


CLEAN :
	-@erase "$(INTDIR)\Domain.obj"
	-@erase "$(INTDIR)\DomainComponent.obj"
	-@erase "$(INTDIR)\EarthquakePattern.obj"
	-@erase "$(INTDIR)\ElementalLoadIter.obj"
	-@erase "$(INTDIR)\LinearSeries.obj"
	-@erase "$(INTDIR)\Load.obj"
	-@erase "$(INTDIR)\LoadPattern.obj"
	-@erase "$(INTDIR)\LoadPatternIter.obj"
	-@erase "$(INTDIR)\MP_Constraint.obj"
	-@erase "$(INTDIR)\NodalLoad.obj"
	-@erase "$(INTDIR)\NodalLoadIter.obj"
	-@erase "$(INTDIR)\Node.obj"
	-@erase "$(INTDIR)\PathSeries.obj"
	-@erase "$(INTDIR)\PathTimeSeries.obj"
	-@erase "$(INTDIR)\RectangularSeries.obj"
	-@erase "$(INTDIR)\RigidFloor.obj"
	-@erase "$(INTDIR)\RigidLink.obj"
	-@erase "$(INTDIR)\SingleDomEleIter.obj"
	-@erase "$(INTDIR)\SingleDomMP_Iter.obj"
	-@erase "$(INTDIR)\SingleDomNodIter.obj"
	-@erase "$(INTDIR)\SingleDomSP_Iter.obj"
	-@erase "$(INTDIR)\SP_Constraint.obj"
	-@erase "$(INTDIR)\Subdomain.obj"
	-@erase "$(INTDIR)\SubdomainNodIter.obj"
	-@erase "$(INTDIR)\TclPatternCommand.obj"
	-@erase "$(INTDIR)\TimeSeries.obj"
	-@erase "$(INTDIR)\UniformExcitation.obj"
	-@erase "$(INTDIR)\vc60.idb"
	-@erase "$(OUTDIR)\domain.lib"

"$(OUTDIR)" :
    if not exist "$(OUTDIR)/$(NULL)" mkdir "$(OUTDIR)"

CPP=cl.exe
CPP_PROJ=/nologo /ML /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_MBCS" /D "_LIB" /Fp"$(INTDIR)\domain.pch" /YX /Fo"$(INTDIR)\\" /Fd"$(INTDIR)\\" /FD /c 

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
BSC32_FLAGS=/nologo /o"$(OUTDIR)\domain.bsc" 
BSC32_SBRS= \
	
LIB32=link.exe -lib
LIB32_FLAGS=/nologo /out:"$(OUTDIR)\domain.lib" 
LIB32_OBJS= \
	"$(INTDIR)\NodalLoad.obj" \
	"$(INTDIR)\Node.obj" \
	"$(INTDIR)\Domain.obj" \
	"$(INTDIR)\SingleDomEleIter.obj" \
	"$(INTDIR)\SingleDomMP_Iter.obj" \
	"$(INTDIR)\SingleDomNodIter.obj" \
	"$(INTDIR)\SingleDomSP_Iter.obj" \
	"$(INTDIR)\Subdomain.obj" \
	"$(INTDIR)\SubdomainNodIter.obj" \
	"$(INTDIR)\DomainComponent.obj" \
	"$(INTDIR)\ElementalLoadIter.obj" \
	"$(INTDIR)\Load.obj" \
	"$(INTDIR)\NodalLoadIter.obj" \
	"$(INTDIR)\MP_Constraint.obj" \
	"$(INTDIR)\RigidFloor.obj" \
	"$(INTDIR)\RigidLink.obj" \
	"$(INTDIR)\SP_Constraint.obj" \
	"$(INTDIR)\EarthquakePattern.obj" \
	"$(INTDIR)\LinearSeries.obj" \
	"$(INTDIR)\LoadPattern.obj" \
	"$(INTDIR)\LoadPatternIter.obj" \
	"$(INTDIR)\PathSeries.obj" \
	"$(INTDIR)\PathTimeSeries.obj" \
	"$(INTDIR)\RectangularSeries.obj" \
	"$(INTDIR)\TclPatternCommand.obj" \
	"$(INTDIR)\TimeSeries.obj" \
	"$(INTDIR)\UniformExcitation.obj"

"$(OUTDIR)\domain.lib" : "$(OUTDIR)" $(DEF_FILE) $(LIB32_OBJS)
    $(LIB32) @<<
  $(LIB32_FLAGS) $(DEF_FLAGS) $(LIB32_OBJS)
<<

!ELSEIF  "$(CFG)" == "domain - Win32 Debug"

OUTDIR=.\..\..\lib
INTDIR=.\..\..\obj\domain
# Begin Custom Macros
OutDir=.\..\..\lib
# End Custom Macros

ALL : "$(OUTDIR)\domain.lib"


CLEAN :
	-@erase "$(INTDIR)\Domain.obj"
	-@erase "$(INTDIR)\DomainComponent.obj"
	-@erase "$(INTDIR)\EarthquakePattern.obj"
	-@erase "$(INTDIR)\ElementalLoadIter.obj"
	-@erase "$(INTDIR)\LinearSeries.obj"
	-@erase "$(INTDIR)\Load.obj"
	-@erase "$(INTDIR)\LoadPattern.obj"
	-@erase "$(INTDIR)\LoadPatternIter.obj"
	-@erase "$(INTDIR)\MP_Constraint.obj"
	-@erase "$(INTDIR)\NodalLoad.obj"
	-@erase "$(INTDIR)\NodalLoadIter.obj"
	-@erase "$(INTDIR)\Node.obj"
	-@erase "$(INTDIR)\PathSeries.obj"
	-@erase "$(INTDIR)\PathTimeSeries.obj"
	-@erase "$(INTDIR)\RectangularSeries.obj"
	-@erase "$(INTDIR)\RigidFloor.obj"
	-@erase "$(INTDIR)\RigidLink.obj"
	-@erase "$(INTDIR)\SingleDomEleIter.obj"
	-@erase "$(INTDIR)\SingleDomMP_Iter.obj"
	-@erase "$(INTDIR)\SingleDomNodIter.obj"
	-@erase "$(INTDIR)\SingleDomSP_Iter.obj"
	-@erase "$(INTDIR)\SP_Constraint.obj"
	-@erase "$(INTDIR)\Subdomain.obj"
	-@erase "$(INTDIR)\SubdomainNodIter.obj"
	-@erase "$(INTDIR)\TclPatternCommand.obj"
	-@erase "$(INTDIR)\TimeSeries.obj"
	-@erase "$(INTDIR)\UniformExcitation.obj"
	-@erase "$(INTDIR)\vc60.idb"
	-@erase "$(OUTDIR)\domain.lib"

"$(OUTDIR)" :
    if not exist "$(OUTDIR)/$(NULL)" mkdir "$(OUTDIR)"

"$(INTDIR)" :
    if not exist "$(INTDIR)/$(NULL)" mkdir "$(INTDIR)"

CPP=cl.exe
CPP_PROJ=/nologo /MLd /W3 /GX /Od /I "c:\Program Files\tcl\include" /I "..\..\src\domain\groundMotion" /I "..\..\src\analysis\fe_ele" /I "..\..\src\utility" /I "..\..\src\domain\subdomain" /I "..\..\src\database" /I "..\..\src\analysis\analysis" /I "..\..\src\recorder" /I "..\..\src\graph\graph" /I "..\..\src\modelbuilder" /I "..\..\src\domain\domain\single" /I "..\..\src\modelbuilder\tcl" /I "..\..\src\tagged\storage" /I "..\..\src\domain\constraints" /I "..\..\src\domain\pattern" /I "..\..\src\matrix" /I "..\..\src\analysis\dof_grp" /I "..\..\src\domain\domain" /I "..\..\src\actor\channel" /I "..\..\src" /I "..\..\src\actor\actor" /I "..\..\src\actor\objectBroker" /I "..\..\src\tagged" /I "..\..\src\domain\component" /I "..\..\src\domain\load" /I "..\..\src\domain\node" /I "..\..\src\element" /I "..\..\src\renderer" /D "WIN32" /D "_DEBUG" /D "_MBCS" /D "_LIB" /Fp"$(INTDIR)\domain.pch" /YX /Fo"$(INTDIR)\\" /Fd"$(INTDIR)\\" /FD /GZ /c 

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
BSC32_FLAGS=/nologo /o"$(OUTDIR)\domain.bsc" 
BSC32_SBRS= \
	
LIB32=link.exe -lib
LIB32_FLAGS=/nologo /out:"$(OUTDIR)\domain.lib" 
LIB32_OBJS= \
	"$(INTDIR)\NodalLoad.obj" \
	"$(INTDIR)\Node.obj" \
	"$(INTDIR)\Domain.obj" \
	"$(INTDIR)\SingleDomEleIter.obj" \
	"$(INTDIR)\SingleDomMP_Iter.obj" \
	"$(INTDIR)\SingleDomNodIter.obj" \
	"$(INTDIR)\SingleDomSP_Iter.obj" \
	"$(INTDIR)\Subdomain.obj" \
	"$(INTDIR)\SubdomainNodIter.obj" \
	"$(INTDIR)\DomainComponent.obj" \
	"$(INTDIR)\ElementalLoadIter.obj" \
	"$(INTDIR)\Load.obj" \
	"$(INTDIR)\NodalLoadIter.obj" \
	"$(INTDIR)\MP_Constraint.obj" \
	"$(INTDIR)\RigidFloor.obj" \
	"$(INTDIR)\RigidLink.obj" \
	"$(INTDIR)\SP_Constraint.obj" \
	"$(INTDIR)\EarthquakePattern.obj" \
	"$(INTDIR)\LinearSeries.obj" \
	"$(INTDIR)\LoadPattern.obj" \
	"$(INTDIR)\LoadPatternIter.obj" \
	"$(INTDIR)\PathSeries.obj" \
	"$(INTDIR)\PathTimeSeries.obj" \
	"$(INTDIR)\RectangularSeries.obj" \
	"$(INTDIR)\TclPatternCommand.obj" \
	"$(INTDIR)\TimeSeries.obj" \
	"$(INTDIR)\UniformExcitation.obj"

"$(OUTDIR)\domain.lib" : "$(OUTDIR)" $(DEF_FILE) $(LIB32_OBJS)
    $(LIB32) @<<
  $(LIB32_FLAGS) $(DEF_FLAGS) $(LIB32_OBJS)
<<

!ENDIF 


!IF "$(NO_EXTERNAL_DEPS)" != "1"
!IF EXISTS("domain.dep")
!INCLUDE "domain.dep"
!ELSE 
!MESSAGE Warning: cannot find "domain.dep"
!ENDIF 
!ENDIF 


!IF "$(CFG)" == "domain - Win32 Release" || "$(CFG)" == "domain - Win32 Debug"
SOURCE=..\..\Src\domain\node\NodalLoad.cpp

"$(INTDIR)\NodalLoad.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\Src\domain\node\Node.cpp

"$(INTDIR)\Node.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\Src\domain\domain\Domain.cpp

"$(INTDIR)\Domain.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\Src\domain\domain\single\SingleDomEleIter.cpp

"$(INTDIR)\SingleDomEleIter.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\Src\domain\domain\single\SingleDomMP_Iter.cpp

"$(INTDIR)\SingleDomMP_Iter.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\Src\domain\domain\single\SingleDomNodIter.cpp

"$(INTDIR)\SingleDomNodIter.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\Src\domain\domain\single\SingleDomSP_Iter.cpp

"$(INTDIR)\SingleDomSP_Iter.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\Src\domain\subdomain\Subdomain.cpp

"$(INTDIR)\Subdomain.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\Src\domain\subdomain\SubdomainNodIter.cpp

"$(INTDIR)\SubdomainNodIter.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\Src\domain\component\DomainComponent.cpp

"$(INTDIR)\DomainComponent.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\Src\domain\load\ElementalLoadIter.cpp

"$(INTDIR)\ElementalLoadIter.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\Src\domain\load\Load.cpp

"$(INTDIR)\Load.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\Src\domain\load\NodalLoadIter.cpp

"$(INTDIR)\NodalLoadIter.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\Src\domain\constraints\MP_Constraint.cpp

"$(INTDIR)\MP_Constraint.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\Src\domain\constraints\RigidFloor.cpp

"$(INTDIR)\RigidFloor.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\Src\domain\constraints\RigidLink.cpp

"$(INTDIR)\RigidLink.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\Src\domain\constraints\SP_Constraint.cpp

"$(INTDIR)\SP_Constraint.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\Src\domain\pattern\EarthquakePattern.cpp

"$(INTDIR)\EarthquakePattern.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\Src\domain\pattern\LinearSeries.cpp

"$(INTDIR)\LinearSeries.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\Src\domain\pattern\LoadPattern.cpp

"$(INTDIR)\LoadPattern.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\Src\domain\pattern\LoadPatternIter.cpp

"$(INTDIR)\LoadPatternIter.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\Src\domain\pattern\PathSeries.cpp

"$(INTDIR)\PathSeries.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\Src\domain\pattern\PathTimeSeries.cpp

"$(INTDIR)\PathTimeSeries.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\Src\domain\pattern\RectangularSeries.cpp

"$(INTDIR)\RectangularSeries.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\Src\domain\pattern\TclPatternCommand.cpp

"$(INTDIR)\TclPatternCommand.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\Src\domain\pattern\TimeSeries.cpp

"$(INTDIR)\TimeSeries.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\Src\domain\pattern\UniformExcitation.cpp

"$(INTDIR)\UniformExcitation.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)



!ENDIF 

