# Microsoft Developer Studio Generated NMAKE File, Based on superLU.dsp
!IF "$(CFG)" == ""
CFG=superLU - Win32 Debug
!MESSAGE No configuration specified. Defaulting to superLU - Win32 Debug.
!ENDIF 

!IF "$(CFG)" != "superLU - Win32 Release" && "$(CFG)" != "superLU - Win32 Debug"
!MESSAGE Invalid configuration "$(CFG)" specified.
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "superLU.mak" CFG="superLU - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "superLU - Win32 Release" (based on "Win32 (x86) Static Library")
!MESSAGE "superLU - Win32 Debug" (based on "Win32 (x86) Static Library")
!MESSAGE 
!ERROR An invalid configuration is specified.
!ENDIF 

!IF "$(OS)" == "Windows_NT"
NULL=
!ELSE 
NULL=nul
!ENDIF 

!IF  "$(CFG)" == "superLU - Win32 Release"

OUTDIR=.\Release
INTDIR=.\Release
# Begin Custom Macros
OutDir=.\Release
# End Custom Macros

ALL : "$(OUTDIR)\superLU.lib"


CLEAN :
	-@erase "$(INTDIR)\dcolumn_bmod.obj"
	-@erase "$(INTDIR)\dcolumn_dfs.obj"
	-@erase "$(INTDIR)\dcopy_to_ucol.obj"
	-@erase "$(INTDIR)\dgscon.obj"
	-@erase "$(INTDIR)\dgsequ.obj"
	-@erase "$(INTDIR)\dgsrfs.obj"
	-@erase "$(INTDIR)\dgssv.obj"
	-@erase "$(INTDIR)\dgssvx.obj"
	-@erase "$(INTDIR)\dgstrf.obj"
	-@erase "$(INTDIR)\dgstrs.obj"
	-@erase "$(INTDIR)\dlacon.obj"
	-@erase "$(INTDIR)\dlamch.obj"
	-@erase "$(INTDIR)\dlangs.obj"
	-@erase "$(INTDIR)\dlaqgs.obj"
	-@erase "$(INTDIR)\dmemory.obj"
	-@erase "$(INTDIR)\dpanel_bmod.obj"
	-@erase "$(INTDIR)\dpanel_dfs.obj"
	-@erase "$(INTDIR)\dpivotgrowth.obj"
	-@erase "$(INTDIR)\dpivotL.obj"
	-@erase "$(INTDIR)\dpruneL.obj"
	-@erase "$(INTDIR)\dreadhb.obj"
	-@erase "$(INTDIR)\dsnode_bmod.obj"
	-@erase "$(INTDIR)\dsnode_dfs.obj"
	-@erase "$(INTDIR)\dsp_blas2.obj"
	-@erase "$(INTDIR)\dsp_blas3.obj"
	-@erase "$(INTDIR)\dutil.obj"
	-@erase "$(INTDIR)\get_perm_c.obj"
	-@erase "$(INTDIR)\lsame.obj"
	-@erase "$(INTDIR)\memory.obj"
	-@erase "$(INTDIR)\mmd.obj"
	-@erase "$(INTDIR)\relax_snode.obj"
	-@erase "$(INTDIR)\sp_coletree.obj"
	-@erase "$(INTDIR)\sp_ienv.obj"
	-@erase "$(INTDIR)\sp_preorder.obj"
	-@erase "$(INTDIR)\superlu_timer.obj"
	-@erase "$(INTDIR)\sutil.obj"
	-@erase "$(INTDIR)\util.obj"
	-@erase "$(INTDIR)\vc60.idb"
	-@erase "$(INTDIR)\xerbla.obj"
	-@erase "$(OUTDIR)\superLU.lib"

"$(OUTDIR)" :
    if not exist "$(OUTDIR)/$(NULL)" mkdir "$(OUTDIR)"

CPP=cl.exe
CPP_PROJ=/nologo /ML /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_MBCS" /D "_LIB" /Fp"$(INTDIR)\superLU.pch" /YX /Fo"$(INTDIR)\\" /Fd"$(INTDIR)\\" /FD /c 

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
BSC32_FLAGS=/nologo /o"$(OUTDIR)\superLU.bsc" 
BSC32_SBRS= \
	
LIB32=link.exe -lib
LIB32_FLAGS=/nologo /out:"$(OUTDIR)\superLU.lib" 
LIB32_OBJS= \
	"$(INTDIR)\dcolumn_bmod.obj" \
	"$(INTDIR)\dcolumn_dfs.obj" \
	"$(INTDIR)\dcopy_to_ucol.obj" \
	"$(INTDIR)\dgscon.obj" \
	"$(INTDIR)\dgsequ.obj" \
	"$(INTDIR)\dgsrfs.obj" \
	"$(INTDIR)\dgssv.obj" \
	"$(INTDIR)\dgssvx.obj" \
	"$(INTDIR)\dgstrf.obj" \
	"$(INTDIR)\dgstrs.obj" \
	"$(INTDIR)\dlacon.obj" \
	"$(INTDIR)\dlamch.obj" \
	"$(INTDIR)\dlangs.obj" \
	"$(INTDIR)\dlaqgs.obj" \
	"$(INTDIR)\dmemory.obj" \
	"$(INTDIR)\dpanel_bmod.obj" \
	"$(INTDIR)\dpanel_dfs.obj" \
	"$(INTDIR)\dpivotgrowth.obj" \
	"$(INTDIR)\dpivotL.obj" \
	"$(INTDIR)\dpruneL.obj" \
	"$(INTDIR)\dreadhb.obj" \
	"$(INTDIR)\dsnode_bmod.obj" \
	"$(INTDIR)\dsnode_dfs.obj" \
	"$(INTDIR)\dsp_blas2.obj" \
	"$(INTDIR)\dsp_blas3.obj" \
	"$(INTDIR)\dutil.obj" \
	"$(INTDIR)\get_perm_c.obj" \
	"$(INTDIR)\lsame.obj" \
	"$(INTDIR)\memory.obj" \
	"$(INTDIR)\mmd.obj" \
	"$(INTDIR)\relax_snode.obj" \
	"$(INTDIR)\sp_coletree.obj" \
	"$(INTDIR)\sp_ienv.obj" \
	"$(INTDIR)\sp_preorder.obj" \
	"$(INTDIR)\superlu_timer.obj" \
	"$(INTDIR)\sutil.obj" \
	"$(INTDIR)\util.obj" \
	"$(INTDIR)\xerbla.obj"

"$(OUTDIR)\superLU.lib" : "$(OUTDIR)" $(DEF_FILE) $(LIB32_OBJS)
    $(LIB32) @<<
  $(LIB32_FLAGS) $(DEF_FLAGS) $(LIB32_OBJS)
<<

!ELSEIF  "$(CFG)" == "superLU - Win32 Debug"

OUTDIR=.\..\..\lib
INTDIR=.\..\..\obj\superLU
# Begin Custom Macros
OutDir=.\..\..\lib
# End Custom Macros

ALL : "$(OUTDIR)\superLU.lib"


CLEAN :
	-@erase "$(INTDIR)\dcolumn_bmod.obj"
	-@erase "$(INTDIR)\dcolumn_dfs.obj"
	-@erase "$(INTDIR)\dcopy_to_ucol.obj"
	-@erase "$(INTDIR)\dgscon.obj"
	-@erase "$(INTDIR)\dgsequ.obj"
	-@erase "$(INTDIR)\dgsrfs.obj"
	-@erase "$(INTDIR)\dgssv.obj"
	-@erase "$(INTDIR)\dgssvx.obj"
	-@erase "$(INTDIR)\dgstrf.obj"
	-@erase "$(INTDIR)\dgstrs.obj"
	-@erase "$(INTDIR)\dlacon.obj"
	-@erase "$(INTDIR)\dlamch.obj"
	-@erase "$(INTDIR)\dlangs.obj"
	-@erase "$(INTDIR)\dlaqgs.obj"
	-@erase "$(INTDIR)\dmemory.obj"
	-@erase "$(INTDIR)\dpanel_bmod.obj"
	-@erase "$(INTDIR)\dpanel_dfs.obj"
	-@erase "$(INTDIR)\dpivotgrowth.obj"
	-@erase "$(INTDIR)\dpivotL.obj"
	-@erase "$(INTDIR)\dpruneL.obj"
	-@erase "$(INTDIR)\dreadhb.obj"
	-@erase "$(INTDIR)\dsnode_bmod.obj"
	-@erase "$(INTDIR)\dsnode_dfs.obj"
	-@erase "$(INTDIR)\dsp_blas2.obj"
	-@erase "$(INTDIR)\dsp_blas3.obj"
	-@erase "$(INTDIR)\dutil.obj"
	-@erase "$(INTDIR)\get_perm_c.obj"
	-@erase "$(INTDIR)\lsame.obj"
	-@erase "$(INTDIR)\memory.obj"
	-@erase "$(INTDIR)\mmd.obj"
	-@erase "$(INTDIR)\relax_snode.obj"
	-@erase "$(INTDIR)\sp_coletree.obj"
	-@erase "$(INTDIR)\sp_ienv.obj"
	-@erase "$(INTDIR)\sp_preorder.obj"
	-@erase "$(INTDIR)\superlu_timer.obj"
	-@erase "$(INTDIR)\sutil.obj"
	-@erase "$(INTDIR)\util.obj"
	-@erase "$(INTDIR)\vc60.idb"
	-@erase "$(INTDIR)\xerbla.obj"
	-@erase "$(OUTDIR)\superLU.lib"

"$(OUTDIR)" :
    if not exist "$(OUTDIR)/$(NULL)" mkdir "$(OUTDIR)"

"$(INTDIR)" :
    if not exist "$(INTDIR)/$(NULL)" mkdir "$(INTDIR)"

CPP=cl.exe
CPP_PROJ=/nologo /MLd /W3 /GX /Od /I "..\..\superLU" /D "WIN32" /D "_DEBUG" /D "_MBCS" /D "_LIB" /D "WINDOWS" /Fp"$(INTDIR)\superLU.pch" /YX /Fo"$(INTDIR)\\" /Fd"$(INTDIR)\\" /FD /GZ /c 

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
BSC32_FLAGS=/nologo /o"$(OUTDIR)\superLU.bsc" 
BSC32_SBRS= \
	
LIB32=link.exe -lib
LIB32_FLAGS=/nologo /out:"$(OUTDIR)\superLU.lib" 
LIB32_OBJS= \
	"$(INTDIR)\dcolumn_bmod.obj" \
	"$(INTDIR)\dcolumn_dfs.obj" \
	"$(INTDIR)\dcopy_to_ucol.obj" \
	"$(INTDIR)\dgscon.obj" \
	"$(INTDIR)\dgsequ.obj" \
	"$(INTDIR)\dgsrfs.obj" \
	"$(INTDIR)\dgssv.obj" \
	"$(INTDIR)\dgssvx.obj" \
	"$(INTDIR)\dgstrf.obj" \
	"$(INTDIR)\dgstrs.obj" \
	"$(INTDIR)\dlacon.obj" \
	"$(INTDIR)\dlamch.obj" \
	"$(INTDIR)\dlangs.obj" \
	"$(INTDIR)\dlaqgs.obj" \
	"$(INTDIR)\dmemory.obj" \
	"$(INTDIR)\dpanel_bmod.obj" \
	"$(INTDIR)\dpanel_dfs.obj" \
	"$(INTDIR)\dpivotgrowth.obj" \
	"$(INTDIR)\dpivotL.obj" \
	"$(INTDIR)\dpruneL.obj" \
	"$(INTDIR)\dreadhb.obj" \
	"$(INTDIR)\dsnode_bmod.obj" \
	"$(INTDIR)\dsnode_dfs.obj" \
	"$(INTDIR)\dsp_blas2.obj" \
	"$(INTDIR)\dsp_blas3.obj" \
	"$(INTDIR)\dutil.obj" \
	"$(INTDIR)\get_perm_c.obj" \
	"$(INTDIR)\lsame.obj" \
	"$(INTDIR)\memory.obj" \
	"$(INTDIR)\mmd.obj" \
	"$(INTDIR)\relax_snode.obj" \
	"$(INTDIR)\sp_coletree.obj" \
	"$(INTDIR)\sp_ienv.obj" \
	"$(INTDIR)\sp_preorder.obj" \
	"$(INTDIR)\superlu_timer.obj" \
	"$(INTDIR)\sutil.obj" \
	"$(INTDIR)\util.obj" \
	"$(INTDIR)\xerbla.obj"

"$(OUTDIR)\superLU.lib" : "$(OUTDIR)" $(DEF_FILE) $(LIB32_OBJS)
    $(LIB32) @<<
  $(LIB32_FLAGS) $(DEF_FLAGS) $(LIB32_OBJS)
<<

!ENDIF 


!IF "$(NO_EXTERNAL_DEPS)" != "1"
!IF EXISTS("superLU.dep")
!INCLUDE "superLU.dep"
!ELSE 
!MESSAGE Warning: cannot find "superLU.dep"
!ENDIF 
!ENDIF 


!IF "$(CFG)" == "superLU - Win32 Release" || "$(CFG)" == "superLU - Win32 Debug"
SOURCE=..\..\SuperLU\dcolumn_bmod.c

"$(INTDIR)\dcolumn_bmod.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\SuperLU\dcolumn_dfs.c

"$(INTDIR)\dcolumn_dfs.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\SuperLU\dcopy_to_ucol.c

"$(INTDIR)\dcopy_to_ucol.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\SuperLU\dgscon.c

"$(INTDIR)\dgscon.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\SuperLU\dgsequ.c

"$(INTDIR)\dgsequ.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\SuperLU\dgsrfs.c

"$(INTDIR)\dgsrfs.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\SuperLU\dgssv.c

"$(INTDIR)\dgssv.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\SuperLU\dgssvx.c

"$(INTDIR)\dgssvx.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\SuperLU\dgstrf.c

"$(INTDIR)\dgstrf.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\SuperLU\dgstrs.c

"$(INTDIR)\dgstrs.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\SuperLU\dlacon.c

"$(INTDIR)\dlacon.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\SuperLU\dlamch.c

"$(INTDIR)\dlamch.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\SuperLU\dlangs.c

"$(INTDIR)\dlangs.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\SuperLU\dlaqgs.c

"$(INTDIR)\dlaqgs.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\SuperLU\dmemory.c

"$(INTDIR)\dmemory.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\SuperLU\dpanel_bmod.c

"$(INTDIR)\dpanel_bmod.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\SuperLU\dpanel_dfs.c

"$(INTDIR)\dpanel_dfs.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\SuperLU\dpivotgrowth.c

"$(INTDIR)\dpivotgrowth.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\SuperLU\dpivotL.c

"$(INTDIR)\dpivotL.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\SuperLU\dpruneL.c

"$(INTDIR)\dpruneL.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\SuperLU\dreadhb.c

"$(INTDIR)\dreadhb.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\SuperLU\dsnode_bmod.c

"$(INTDIR)\dsnode_bmod.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\SuperLU\dsnode_dfs.c

"$(INTDIR)\dsnode_dfs.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\SuperLU\dsp_blas2.c

"$(INTDIR)\dsp_blas2.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\SuperLU\dsp_blas3.c

"$(INTDIR)\dsp_blas3.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\SuperLU\dutil.c

"$(INTDIR)\dutil.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\SuperLU\get_perm_c.c

"$(INTDIR)\get_perm_c.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\SuperLU\lsame.c

"$(INTDIR)\lsame.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\SuperLU\memory.c

"$(INTDIR)\memory.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\SuperLU\mmd.c

"$(INTDIR)\mmd.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\SuperLU\relax_snode.c

"$(INTDIR)\relax_snode.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\SuperLU\sp_coletree.c

"$(INTDIR)\sp_coletree.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\SuperLU\sp_ienv.c

"$(INTDIR)\sp_ienv.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\SuperLU\sp_preorder.c

"$(INTDIR)\sp_preorder.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\SuperLU\superlu_timer.c

"$(INTDIR)\superlu_timer.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\SuperLU\sutil.c

"$(INTDIR)\sutil.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\SuperLU\util.c

"$(INTDIR)\util.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\SuperLU\xerbla.c

"$(INTDIR)\xerbla.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)



!ENDIF 

