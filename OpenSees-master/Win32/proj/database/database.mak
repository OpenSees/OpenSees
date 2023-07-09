# Microsoft Developer Studio Generated NMAKE File, Based on database.dsp
!IF "$(CFG)" == ""
CFG=database - Win32 Debug
!MESSAGE No configuration specified. Defaulting to database - Win32 Debug.
!ENDIF 

!IF "$(CFG)" != "database - Win32 Release" && "$(CFG)" != "database - Win32 Debug"
!MESSAGE Invalid configuration "$(CFG)" specified.
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "database.mak" CFG="database - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "database - Win32 Release" (based on "Win32 (x86) Static Library")
!MESSAGE "database - Win32 Debug" (based on "Win32 (x86) Static Library")
!MESSAGE 
!ERROR An invalid configuration is specified.
!ENDIF 

!IF "$(OS)" == "Windows_NT"
NULL=
!ELSE 
NULL=nul
!ENDIF 

!IF  "$(CFG)" == "database - Win32 Release"

OUTDIR=.\..\..\lib
INTDIR=.\..\..\obj\database
# Begin Custom Macros
OutDir=.\..\..\lib
# End Custom Macros

ALL : "$(OUTDIR)\database.lib"


CLEAN :
	-@erase "$(INTDIR)\FE_Datastore.obj"
	-@erase "$(INTDIR)\FileDatastore.obj"
	-@erase "$(INTDIR)\TclDatabaseCommands.obj"
	-@erase "$(INTDIR)\vc60.idb"
	-@erase "$(OUTDIR)\database.lib"

"$(OUTDIR)" :
    if not exist "$(OUTDIR)/$(NULL)" mkdir "$(OUTDIR)"

"$(INTDIR)" :
    if not exist "$(INTDIR)/$(NULL)" mkdir "$(INTDIR)"

CPP=cl.exe
CPP_PROJ=/nologo /MLd /W3 /GX /O2 /I "..\..\..\src\element\fourNodeQuad" /I "..\..\..\src\analysis\algorithm" /I "equiSolnAlgo" /I "..\..\..\src\handler" /I "c:\Program Files\tcl\include" /I "..\..\..\src\analysis\algorithm\equiSolnAlgo" /I "..\..\..\src\database" /I "..\..\..\src\actor\channel" /I ".\..\..\src\domain\groundMotion" /I "..\..\..\src\analysis\fe_ele" /I "..\..\..\src\utility" /I "..\..\..\src\domain\subdomain" /I "..\..\..\src\analysis\analysis" /I "..\..\..\src\recorder" /I "..\..\..\src\graph\graph" /I "..\..\..\src\modelbuilder" /I "..\..\..\src\domain\domain\single" /I "..\..\..\src\modelbuilder\tcl" /I "..\..\..\src\tagged\storage" /I "..\..\..\src\domain\constraints" /I "..\..\..\src\domain\pattern" /I "..\..\..\src\matrix" /I "..\..\..\src\analysis\dof_grp" /I "..\..\..\src\domain\domain" /I "..\..\..\src" /I "..\..\..\src\actor\actor" /I "..\..\..\src\actor\objectBroker" /I "..\..\..\src\tagged" /I "..\..\..\src\domain\component" /I "..\..\..\src\domain\load" /I "..\..\..\src\domain\node" /I "..\..\..\src\element" /I "..\..\..\src\renderer" /D "WIN32" /D "NDEBUG" /D "_MBCS" /D "_LIB" /Fo"$(INTDIR)\\" /Fd"$(INTDIR)\\" /FD /c 

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
BSC32_FLAGS=/nologo /o"$(OUTDIR)\database.bsc" 
BSC32_SBRS= \
	
LIB32=link.exe -lib
LIB32_FLAGS=/nologo /out:"$(OUTDIR)\database.lib" 
LIB32_OBJS= \
	"$(INTDIR)\FE_Datastore.obj" \
	"$(INTDIR)\FileDatastore.obj" \
	"$(INTDIR)\TclDatabaseCommands.obj"

"$(OUTDIR)\database.lib" : "$(OUTDIR)" $(DEF_FILE) $(LIB32_OBJS)
    $(LIB32) @<<
  $(LIB32_FLAGS) $(DEF_FLAGS) $(LIB32_OBJS)
<<

!ELSEIF  "$(CFG)" == "database - Win32 Debug"

OUTDIR=.\..\..\lib
INTDIR=.\..\..\obj\database
# Begin Custom Macros
OutDir=.\..\..\lib
# End Custom Macros

ALL : "$(OUTDIR)\database.lib" "$(OUTDIR)\database.bsc"


CLEAN :
	-@erase "$(INTDIR)\FE_Datastore.obj"
	-@erase "$(INTDIR)\FE_Datastore.sbr"
	-@erase "$(INTDIR)\FileDatastore.obj"
	-@erase "$(INTDIR)\FileDatastore.sbr"
	-@erase "$(INTDIR)\TclDatabaseCommands.obj"
	-@erase "$(INTDIR)\TclDatabaseCommands.sbr"
	-@erase "$(INTDIR)\vc60.idb"
	-@erase "$(INTDIR)\vc60.pdb"
	-@erase "$(OUTDIR)\database.bsc"
	-@erase "$(OUTDIR)\database.lib"

"$(OUTDIR)" :
    if not exist "$(OUTDIR)/$(NULL)" mkdir "$(OUTDIR)"

"$(INTDIR)" :
    if not exist "$(INTDIR)/$(NULL)" mkdir "$(INTDIR)"

CPP=cl.exe
CPP_PROJ=/nologo /MLd /W3 /Gm /GX /ZI /Od /I "..\..\..\src\analysis\algorithm" /I "equiSolnAlgo" /I "..\..\..\src\handler" /I "c:\Program Files\tcl\include" /I "..\..\..\src\analysis\algorithm\equiSolnAlgo" /I "..\..\..\src\database" /I "..\..\..\src\actor\channel" /I ".\..\..\src\domain\groundMotion" /I "..\..\..\src\analysis\fe_ele" /I "..\..\..\src\utility" /I "..\..\..\src\domain\subdomain" /I "..\..\..\src\analysis\analysis" /I "..\..\..\src\recorder" /I "..\..\..\src\graph\graph" /I "..\..\..\src\modelbuilder" /I "..\..\..\src\domain\domain\single" /I "..\..\..\src\modelbuilder\tcl" /I "..\..\..\src\tagged\storage" /I "..\..\..\src\domain\constraints" /I "..\..\..\src\domain\pattern" /I "..\..\..\src\matrix" /I "..\..\..\src\analysis\dof_grp" /I "..\..\..\src\domain\domain" /I "..\..\..\src" /I "..\..\..\src\actor\actor" /I "..\..\..\src\actor\objectBroker" /I "..\..\..\src\tagged" /I "..\..\..\src\domain\component" /I "..\..\..\src\domain\load" /I "..\..\..\src\domain\node" /I "..\..\..\src\element" /I "..\..\..\src\renderer" /D "WIN32" /D "_DEBUG" /D "_MBCS" /D "_LIB" /FR"$(INTDIR)\\" /Fo"$(INTDIR)\\" /Fd"$(INTDIR)\\" /FD /GZ /c 

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
BSC32_FLAGS=/nologo /o"$(OUTDIR)\database.bsc" 
BSC32_SBRS= \
	"$(INTDIR)\FE_Datastore.sbr" \
	"$(INTDIR)\FileDatastore.sbr" \
	"$(INTDIR)\TclDatabaseCommands.sbr"

"$(OUTDIR)\database.bsc" : "$(OUTDIR)" $(BSC32_SBRS)
    $(BSC32) @<<
  $(BSC32_FLAGS) $(BSC32_SBRS)
<<

LIB32=link.exe -lib
LIB32_FLAGS=/nologo /out:"$(OUTDIR)\database.lib" 
LIB32_OBJS= \
	"$(INTDIR)\FE_Datastore.obj" \
	"$(INTDIR)\FileDatastore.obj" \
	"$(INTDIR)\TclDatabaseCommands.obj"

"$(OUTDIR)\database.lib" : "$(OUTDIR)" $(DEF_FILE) $(LIB32_OBJS)
    $(LIB32) @<<
  $(LIB32_FLAGS) $(DEF_FLAGS) $(LIB32_OBJS)
<<

!ENDIF 


!IF "$(NO_EXTERNAL_DEPS)" != "1"
!IF EXISTS("database.dep")
!INCLUDE "database.dep"
!ELSE 
!MESSAGE Warning: cannot find "database.dep"
!ENDIF 
!ENDIF 


!IF "$(CFG)" == "database - Win32 Release" || "$(CFG)" == "database - Win32 Debug"
SOURCE=..\..\..\SRC\database\FE_Datastore.cpp

!IF  "$(CFG)" == "database - Win32 Release"


"$(INTDIR)\FE_Datastore.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "database - Win32 Debug"


"$(INTDIR)\FE_Datastore.obj"	"$(INTDIR)\FE_Datastore.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\..\SRC\database\FileDatastore.cpp

!IF  "$(CFG)" == "database - Win32 Release"


"$(INTDIR)\FileDatastore.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "database - Win32 Debug"


"$(INTDIR)\FileDatastore.obj"	"$(INTDIR)\FileDatastore.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\..\SRC\database\TclDatabaseCommands.cpp

!IF  "$(CFG)" == "database - Win32 Release"


"$(INTDIR)\TclDatabaseCommands.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "database - Win32 Debug"


"$(INTDIR)\TclDatabaseCommands.obj"	"$(INTDIR)\TclDatabaseCommands.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 


!ENDIF 

