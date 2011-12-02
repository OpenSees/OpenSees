# Microsoft Developer Studio Generated NMAKE File, Format Version 4.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Static Library" 0x0104

!IF "$(CFG)" == ""
CFG=OpenSees - Win32 Debug
!MESSAGE No configuration specified.  Defaulting to OpenSees - Win32 Debug.
!ENDIF 

!IF "$(CFG)" != "g3Fortran - Win32 Release" && "$(CFG)" !=\
 "g3Fortran - Win32 Debug" && "$(CFG)" != "ARPACK - Win32 Release" && "$(CFG)"\
 != "ARPACK - Win32 Debug" && "$(CFG)" != "LAPACK - Win32 Release" && "$(CFG)"\
 != "LAPACK - Win32 Debug" && "$(CFG)" != "BLAS - Win32 Release" && "$(CFG)" !=\
 "BLAS - Win32 Debug" && "$(CFG)" != "UmfPack - Win32 Release" && "$(CFG)" !=\
 "UmfPack - Win32 Debug" && "$(CFG)" != "feap - Win32 Release" && "$(CFG)" !=\
 "feap - Win32 Debug" && "$(CFG)" != "OpenSees - Win32 Release" && "$(CFG)" !=\
 "OpenSees - Win32 Debug"
!MESSAGE Invalid configuration "$(CFG)" specified.
!MESSAGE You can specify a configuration when running NMAKE on this makefile
!MESSAGE by defining the macro CFG on the command line.  For example:
!MESSAGE 
!MESSAGE NMAKE /f "OpenSeesFortran.mak" CFG="OpenSees - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "g3Fortran - Win32 Release" (based on "Win32 (x86) Static Library")
!MESSAGE "g3Fortran - Win32 Debug" (based on "Win32 (x86) Static Library")
!MESSAGE "ARPACK - Win32 Release" (based on "Win32 (x86) Static Library")
!MESSAGE "ARPACK - Win32 Debug" (based on "Win32 (x86) Static Library")
!MESSAGE "LAPACK - Win32 Release" (based on "Win32 (x86) Static Library")
!MESSAGE "LAPACK - Win32 Debug" (based on "Win32 (x86) Static Library")
!MESSAGE "BLAS - Win32 Release" (based on "Win32 (x86) Static Library")
!MESSAGE "BLAS - Win32 Debug" (based on "Win32 (x86) Static Library")
!MESSAGE "UmfPack - Win32 Release" (based on "Win32 (x86) Static Library")
!MESSAGE "UmfPack - Win32 Debug" (based on "Win32 (x86) Static Library")
!MESSAGE "feap - Win32 Release" (based on "Win32 (x86) Static Library")
!MESSAGE "feap - Win32 Debug" (based on "Win32 (x86) Static Library")
!MESSAGE "OpenSees - Win32 Release" (based on "Win32 (x86) Static Library")
!MESSAGE "OpenSees - Win32 Debug" (based on "Win32 (x86) Static Library")
!MESSAGE 
!ERROR An invalid configuration is specified.
!ENDIF 

!IF "$(OS)" == "Windows_NT"
NULL=
!ELSE 
NULL=nul
!ENDIF 
################################################################################
# Begin Project
# PROP Target_Last_Scanned "UmfPack - Win32 Release"
F90=fl32.exe

!IF  "$(CFG)" == "g3Fortran - Win32 Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "Release"
# PROP BASE Intermediate_Dir "Release"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir ".\lib"
# PROP Intermediate_Dir ".\obj\fortran"
# PROP Target_Dir ""
OUTDIR=.\.\lib
INTDIR=.\.\obj\fortran

ALL : 

CLEAN : 
	-@erase 

"$(OUTDIR)" :
    if not exist "$(OUTDIR)/$(NULL)" mkdir "$(OUTDIR)"

"$(INTDIR)" :
    if not exist "$(INTDIR)/$(NULL)" mkdir "$(INTDIR)"

# ADD BASE F90 /Ox /I "Release/" /c /nologo
# ADD F90 /Ox /I ".\obj\fortran/" /c /nologo
F90_PROJ=/Ox /I ".\obj\fortran/" /c /nologo /Fo".\obj\fortran/" 
F90_OBJS=.\.\obj\fortran/
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
BSC32_FLAGS=/nologo /o"$(OUTDIR)/OpenSeesFortran.bsc" 
BSC32_SBRS=
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo
LIB32_FLAGS=/nologo /out:"$(OUTDIR)/OpenSeesFortran.lib" 
LIB32_OBJS=

!ELSEIF  "$(CFG)" == "g3Fortran - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "Debug"
# PROP BASE Intermediate_Dir "Debug"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir ".\lib"
# PROP Intermediate_Dir ".\obj\fortran"
# PROP Target_Dir ""
OUTDIR=.\.\lib
INTDIR=.\.\obj\fortran

ALL : 

CLEAN : 
	-@erase 

"$(OUTDIR)" :
    if not exist "$(OUTDIR)/$(NULL)" mkdir "$(OUTDIR)"

"$(INTDIR)" :
    if not exist "$(INTDIR)/$(NULL)" mkdir "$(INTDIR)"

# ADD BASE F90 /Z7 /I "Debug/" /c /nologo
# ADD F90 /Z7 /I ".\obj\fortran/" /c /nologo
F90_PROJ=/Z7 /I ".\obj\fortran/" /c /nologo /Fo".\obj\fortran/" 
F90_OBJS=.\.\obj\fortran/
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
BSC32_FLAGS=/nologo /o"$(OUTDIR)/OpenSeesFortran.bsc" 
BSC32_SBRS=
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo
LIB32_FLAGS=/nologo /out:"$(OUTDIR)/OpenSeesFortran.lib" 
LIB32_OBJS=

!ELSEIF  "$(CFG)" == "ARPACK - Win32 Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "ARPACK\Release"
# PROP BASE Intermediate_Dir "ARPACK\Release"
# PROP BASE Target_Dir "ARPACK"
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir ".\lib"
# PROP Intermediate_Dir ".\obj\fortran"
# PROP Target_Dir "ARPACK"
OUTDIR=.\.\lib
INTDIR=.\.\obj\fortran

ALL : "BLAS - Win32 Release" "$(OUTDIR)\ARPACK.lib"

CLEAN : 
	-@erase ".\lib\ARPACK.lib"
	-@erase ".\obj\fortran\iset.obj"
	-@erase ".\obj\fortran\dgetv0.obj"
	-@erase ".\obj\fortran\dsortc.obj"
	-@erase ".\obj\fortran\dnconv.obj"
	-@erase ".\obj\fortran\zvout.obj"
	-@erase ".\obj\fortran\smout.obj"
	-@erase ".\obj\fortran\icnteq.obj"
	-@erase ".\obj\fortran\cmout.obj"
	-@erase ".\obj\fortran\iswap.obj"
	-@erase ".\obj\fortran\dsaitr.obj"
	-@erase ".\obj\fortran\dstqrb.obj"
	-@erase ".\obj\fortran\dvout.obj"
	-@erase ".\obj\fortran\dsapps.obj"
	-@erase ".\obj\fortran\dsaup2.obj"
	-@erase ".\obj\fortran\dsaupd.obj"
	-@erase ".\obj\fortran\ivout.obj"
	-@erase ".\obj\fortran\dneupd.obj"
	-@erase ".\obj\fortran\dsesrt.obj"
	-@erase ".\obj\fortran\dstats.obj"
	-@erase ".\obj\fortran\svout.obj"
	-@erase ".\obj\fortran\cvout.obj"
	-@erase ".\obj\fortran\dseupd.obj"
	-@erase ".\obj\fortran\dsconv.obj"
	-@erase ".\obj\fortran\dstatn.obj"
	-@erase ".\obj\fortran\zmout.obj"
	-@erase ".\obj\fortran\dnaitr.obj"
	-@erase ".\obj\fortran\dseigt.obj"
	-@erase ".\obj\fortran\dsortr.obj"
	-@erase ".\obj\fortran\icopy.obj"
	-@erase ".\obj\fortran\dlaqrb.obj"
	-@erase ".\obj\fortran\dneigh.obj"
	-@erase ".\obj\fortran\dngets.obj"
	-@erase ".\obj\fortran\dnapps.obj"
	-@erase ".\obj\fortran\dnaup2.obj"
	-@erase ".\obj\fortran\dnaupd.obj"
	-@erase ".\obj\fortran\dmout.obj"
	-@erase ".\obj\fortran\second.obj"
	-@erase ".\obj\fortran\dsgets.obj"

"$(OUTDIR)" :
    if not exist "$(OUTDIR)/$(NULL)" mkdir "$(OUTDIR)"

"$(INTDIR)" :
    if not exist "$(INTDIR)/$(NULL)" mkdir "$(INTDIR)"

# ADD BASE F90 /Ox /I "ARPACK\Release/" /c /nologo
# ADD F90 /Ox /I ".\obj\fortran/" /c /nologo
F90_PROJ=/Ox /I ".\obj\fortran/" /c /nologo /Fo".\obj\fortran/" 
F90_OBJS=.\.\obj\fortran/
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
BSC32_FLAGS=/nologo /o"$(OUTDIR)/ARPACK.bsc" 
BSC32_SBRS=
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo
LIB32_FLAGS=/nologo /out:"$(OUTDIR)/ARPACK.lib" 
LIB32_OBJS= \
	"$(INTDIR)/iset.obj" \
	"$(INTDIR)/dgetv0.obj" \
	"$(INTDIR)/dsortc.obj" \
	"$(INTDIR)/dnconv.obj" \
	"$(INTDIR)/zvout.obj" \
	"$(INTDIR)/smout.obj" \
	"$(INTDIR)/icnteq.obj" \
	"$(INTDIR)/cmout.obj" \
	"$(INTDIR)/iswap.obj" \
	"$(INTDIR)/dsaitr.obj" \
	"$(INTDIR)/dstqrb.obj" \
	"$(INTDIR)/dvout.obj" \
	"$(INTDIR)/dsapps.obj" \
	"$(INTDIR)/dsaup2.obj" \
	"$(INTDIR)/dsaupd.obj" \
	"$(INTDIR)/ivout.obj" \
	"$(INTDIR)/dneupd.obj" \
	"$(INTDIR)/dsesrt.obj" \
	"$(INTDIR)/dstats.obj" \
	"$(INTDIR)/svout.obj" \
	"$(INTDIR)/cvout.obj" \
	"$(INTDIR)/dseupd.obj" \
	"$(INTDIR)/dsconv.obj" \
	"$(INTDIR)/dstatn.obj" \
	"$(INTDIR)/zmout.obj" \
	"$(INTDIR)/dnaitr.obj" \
	"$(INTDIR)/dseigt.obj" \
	"$(INTDIR)/dsortr.obj" \
	"$(INTDIR)/icopy.obj" \
	"$(INTDIR)/dlaqrb.obj" \
	"$(INTDIR)/dneigh.obj" \
	"$(INTDIR)/dngets.obj" \
	"$(INTDIR)/dnapps.obj" \
	"$(INTDIR)/dnaup2.obj" \
	"$(INTDIR)/dnaupd.obj" \
	"$(INTDIR)/dmout.obj" \
	"$(INTDIR)/second.obj" \
	"$(INTDIR)/dsgets.obj" \
	".\lib\BLAS.lib"

"$(OUTDIR)\ARPACK.lib" : "$(OUTDIR)" $(DEF_FILE) $(LIB32_OBJS)
    $(LIB32) @<<
  $(LIB32_FLAGS) $(DEF_FLAGS) $(LIB32_OBJS)
<<

!ELSEIF  "$(CFG)" == "ARPACK - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "ARPACK\Debug"
# PROP BASE Intermediate_Dir "ARPACK\Debug"
# PROP BASE Target_Dir "ARPACK"
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir ".\lib"
# PROP Intermediate_Dir ".\obj\fortran"
# PROP Target_Dir "ARPACK"
OUTDIR=.\.\lib
INTDIR=.\.\obj\fortran

ALL : "BLAS - Win32 Debug" "$(OUTDIR)\ARPACK.lib"

CLEAN : 
	-@erase ".\lib\ARPACK.lib"
	-@erase ".\obj\fortran\iset.obj"
	-@erase ".\obj\fortran\dgetv0.obj"
	-@erase ".\obj\fortran\dsortc.obj"
	-@erase ".\obj\fortran\dnconv.obj"
	-@erase ".\obj\fortran\zvout.obj"
	-@erase ".\obj\fortran\smout.obj"
	-@erase ".\obj\fortran\icnteq.obj"
	-@erase ".\obj\fortran\cmout.obj"
	-@erase ".\obj\fortran\iswap.obj"
	-@erase ".\obj\fortran\dsaitr.obj"
	-@erase ".\obj\fortran\dstqrb.obj"
	-@erase ".\obj\fortran\dvout.obj"
	-@erase ".\obj\fortran\dsapps.obj"
	-@erase ".\obj\fortran\dsaup2.obj"
	-@erase ".\obj\fortran\dsaupd.obj"
	-@erase ".\obj\fortran\ivout.obj"
	-@erase ".\obj\fortran\dneupd.obj"
	-@erase ".\obj\fortran\dsesrt.obj"
	-@erase ".\obj\fortran\dstats.obj"
	-@erase ".\obj\fortran\svout.obj"
	-@erase ".\obj\fortran\cvout.obj"
	-@erase ".\obj\fortran\dseupd.obj"
	-@erase ".\obj\fortran\dsconv.obj"
	-@erase ".\obj\fortran\dstatn.obj"
	-@erase ".\obj\fortran\zmout.obj"
	-@erase ".\obj\fortran\dnaitr.obj"
	-@erase ".\obj\fortran\dseigt.obj"
	-@erase ".\obj\fortran\dsortr.obj"
	-@erase ".\obj\fortran\icopy.obj"
	-@erase ".\obj\fortran\dlaqrb.obj"
	-@erase ".\obj\fortran\dneigh.obj"
	-@erase ".\obj\fortran\dngets.obj"
	-@erase ".\obj\fortran\dnapps.obj"
	-@erase ".\obj\fortran\dnaup2.obj"
	-@erase ".\obj\fortran\dnaupd.obj"
	-@erase ".\obj\fortran\dmout.obj"
	-@erase ".\obj\fortran\second.obj"
	-@erase ".\obj\fortran\dsgets.obj"

"$(OUTDIR)" :
    if not exist "$(OUTDIR)/$(NULL)" mkdir "$(OUTDIR)"

"$(INTDIR)" :
    if not exist "$(INTDIR)/$(NULL)" mkdir "$(INTDIR)"

# ADD BASE F90 /Z7 /I "ARPACK\Debug/" /c /nologo
# ADD F90 /Z7 /I ".\obj\fortran/" /c /nologo
F90_PROJ=/Z7 /I ".\obj\fortran/" /c /nologo /Fo".\obj\fortran/" 
F90_OBJS=.\.\obj\fortran/
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
BSC32_FLAGS=/nologo /o"$(OUTDIR)/ARPACK.bsc" 
BSC32_SBRS=
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo
LIB32_FLAGS=/nologo /out:"$(OUTDIR)/ARPACK.lib" 
LIB32_OBJS= \
	"$(INTDIR)/iset.obj" \
	"$(INTDIR)/dgetv0.obj" \
	"$(INTDIR)/dsortc.obj" \
	"$(INTDIR)/dnconv.obj" \
	"$(INTDIR)/zvout.obj" \
	"$(INTDIR)/smout.obj" \
	"$(INTDIR)/icnteq.obj" \
	"$(INTDIR)/cmout.obj" \
	"$(INTDIR)/iswap.obj" \
	"$(INTDIR)/dsaitr.obj" \
	"$(INTDIR)/dstqrb.obj" \
	"$(INTDIR)/dvout.obj" \
	"$(INTDIR)/dsapps.obj" \
	"$(INTDIR)/dsaup2.obj" \
	"$(INTDIR)/dsaupd.obj" \
	"$(INTDIR)/ivout.obj" \
	"$(INTDIR)/dneupd.obj" \
	"$(INTDIR)/dsesrt.obj" \
	"$(INTDIR)/dstats.obj" \
	"$(INTDIR)/svout.obj" \
	"$(INTDIR)/cvout.obj" \
	"$(INTDIR)/dseupd.obj" \
	"$(INTDIR)/dsconv.obj" \
	"$(INTDIR)/dstatn.obj" \
	"$(INTDIR)/zmout.obj" \
	"$(INTDIR)/dnaitr.obj" \
	"$(INTDIR)/dseigt.obj" \
	"$(INTDIR)/dsortr.obj" \
	"$(INTDIR)/icopy.obj" \
	"$(INTDIR)/dlaqrb.obj" \
	"$(INTDIR)/dneigh.obj" \
	"$(INTDIR)/dngets.obj" \
	"$(INTDIR)/dnapps.obj" \
	"$(INTDIR)/dnaup2.obj" \
	"$(INTDIR)/dnaupd.obj" \
	"$(INTDIR)/dmout.obj" \
	"$(INTDIR)/second.obj" \
	"$(INTDIR)/dsgets.obj" \
	".\lib\BLAS.lib"

"$(OUTDIR)\ARPACK.lib" : "$(OUTDIR)" $(DEF_FILE) $(LIB32_OBJS)
    $(LIB32) @<<
  $(LIB32_FLAGS) $(DEF_FLAGS) $(LIB32_OBJS)
<<

!ELSEIF  "$(CFG)" == "LAPACK - Win32 Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "LAPACK\Release"
# PROP BASE Intermediate_Dir "LAPACK\Release"
# PROP BASE Target_Dir "LAPACK"
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir ".\lib"
# PROP Intermediate_Dir ".\obj\fortran"
# PROP Target_Dir "LAPACK"
OUTDIR=.\.\lib
INTDIR=.\.\obj\fortran

ALL : "BLAS - Win32 Release" "$(OUTDIR)\LAPACK.lib"

CLEAN : 
	-@erase ".\lib\LAPACK.lib"
	-@erase ".\obj\fortran\dlabad.obj"
	-@erase ".\obj\fortran\dpbsv.obj"
	-@erase ".\obj\fortran\dlanv2.obj"
	-@erase ".\obj\fortran\dlaptm.obj"
	-@erase ".\obj\fortran\dgetf2.obj"
	-@erase ".\obj\fortran\dgttrf.obj"
	-@erase ".\obj\fortran\dpbtrf.obj"
	-@erase ".\obj\fortran\dlanhs.obj"
	-@erase ".\obj\fortran\dtrevc.obj"
	-@erase ".\obj\fortran\dlaran.obj"
	-@erase ".\obj\fortran\dlassq.obj"
	-@erase ".\obj\fortran\dgttrs.obj"
	-@erase ".\obj\fortran\dlamch.obj"
	-@erase ".\obj\fortran\ilaenv.obj"
	-@erase ".\obj\fortran\dpbtrs.obj"
	-@erase ".\obj\fortran\dladiv.obj"
	-@erase ".\obj\fortran\dlaln2.obj"
	-@erase ".\obj\fortran\izmax1.obj"
	-@erase ".\obj\fortran\dgesv.obj"
	-@erase ".\obj\fortran\dgetrf.obj"
	-@erase ".\obj\fortran\dlaswp.obj"
	-@erase ".\obj\fortran\dorm2r.obj"
	-@erase ".\obj\fortran\dlarfg.obj"
	-@erase ".\obj\fortran\dlanst.obj"
	-@erase ".\obj\fortran\dlasy2.obj"
	-@erase ".\obj\fortran\xerbla.obj"
	-@erase ".\obj\fortran\dlapy3.obj"
	-@erase ".\obj\fortran\dgbtf2.obj"
	-@erase ".\obj\fortran\dlartg.obj"
	-@erase ".\obj\fortran\dgetrs.obj"
	-@erase ".\obj\fortran\dzsum1.obj"
	-@erase ".\obj\fortran\dlange.obj"
	-@erase ".\obj\fortran\dlaset.obj"
	-@erase ".\obj\fortran\dlaexc.obj"
	-@erase ".\obj\fortran\xlaenv.obj"
	-@erase ".\obj\fortran\dpttrf.obj"
	-@erase ".\obj\fortran\dpotf2.obj"
	-@erase ".\obj\fortran\dgbsv.obj"
	-@erase ".\obj\fortran\dlacon.obj"
	-@erase ".\obj\fortran\dlarfx.obj"
	-@erase ".\obj\fortran\lsame.obj"
	-@erase ".\obj\fortran\dtrexc.obj"
	-@erase ".\obj\fortran\dgbtrf.obj"
	-@erase ".\obj\fortran\dlapy2.obj"
	-@erase ".\obj\fortran\lsamen.obj"
	-@erase ".\obj\fortran\dpttrs.obj"
	-@erase ".\obj\fortran\dsteqr.obj"
	-@erase ".\obj\fortran\dtrsyl.obj"
	-@erase ".\obj\fortran\dgbtrs.obj"
	-@erase ".\obj\fortran\dlarnd.obj"
	-@erase ".\obj\fortran\icmax1.obj"
	-@erase ".\obj\fortran\dlarnv.obj"
	-@erase ".\obj\fortran\dlaev2.obj"
	-@erase ".\obj\fortran\dlascl.obj"
	-@erase ".\obj\fortran\dlaruv.obj"
	-@erase ".\obj\fortran\dlasrt.obj"
	-@erase ".\obj\fortran\dlagtm.obj"
	-@erase ".\obj\fortran\dlarf.obj"
	-@erase ".\obj\fortran\dpbtf2.obj"
	-@erase ".\obj\fortran\dlacpy.obj"
	-@erase ".\obj\fortran\dlahqr.obj"
	-@erase ".\obj\fortran\dtrsen.obj"
	-@erase ".\obj\fortran\dlae2.obj"
	-@erase ".\obj\fortran\dlasr.obj"
	-@erase ".\obj\fortran\dgeqr2.obj"

"$(OUTDIR)" :
    if not exist "$(OUTDIR)/$(NULL)" mkdir "$(OUTDIR)"

"$(INTDIR)" :
    if not exist "$(INTDIR)/$(NULL)" mkdir "$(INTDIR)"

# ADD BASE F90 /Ox /I "LAPACK\Release/" /c /nologo
# ADD F90 /Ox /I ".\obj\fortran/" /c /nologo
F90_PROJ=/Ox /I ".\obj\fortran/" /c /nologo /Fo".\obj\fortran/" 
F90_OBJS=.\.\obj\fortran/
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
BSC32_FLAGS=/nologo /o"$(OUTDIR)/LAPACK.bsc" 
BSC32_SBRS=
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo
LIB32_FLAGS=/nologo /out:"$(OUTDIR)/LAPACK.lib" 
LIB32_OBJS= \
	"$(INTDIR)/dlabad.obj" \
	"$(INTDIR)/dpbsv.obj" \
	"$(INTDIR)/dlanv2.obj" \
	"$(INTDIR)/dlaptm.obj" \
	"$(INTDIR)/dgetf2.obj" \
	"$(INTDIR)/dgttrf.obj" \
	"$(INTDIR)/dpbtrf.obj" \
	"$(INTDIR)/dlanhs.obj" \
	"$(INTDIR)/dtrevc.obj" \
	"$(INTDIR)/dlaran.obj" \
	"$(INTDIR)/dlassq.obj" \
	"$(INTDIR)/dgttrs.obj" \
	"$(INTDIR)/dlamch.obj" \
	"$(INTDIR)/ilaenv.obj" \
	"$(INTDIR)/dpbtrs.obj" \
	"$(INTDIR)/dladiv.obj" \
	"$(INTDIR)/dlaln2.obj" \
	"$(INTDIR)/izmax1.obj" \
	"$(INTDIR)/dgesv.obj" \
	"$(INTDIR)/dgetrf.obj" \
	"$(INTDIR)/dlaswp.obj" \
	"$(INTDIR)/dorm2r.obj" \
	"$(INTDIR)/dlarfg.obj" \
	"$(INTDIR)/dlanst.obj" \
	"$(INTDIR)/dlasy2.obj" \
	"$(INTDIR)/xerbla.obj" \
	"$(INTDIR)/dlapy3.obj" \
	"$(INTDIR)/dgbtf2.obj" \
	"$(INTDIR)/dlartg.obj" \
	"$(INTDIR)/dgetrs.obj" \
	"$(INTDIR)/dzsum1.obj" \
	"$(INTDIR)/dlange.obj" \
	"$(INTDIR)/dlaset.obj" \
	"$(INTDIR)/dlaexc.obj" \
	"$(INTDIR)/xlaenv.obj" \
	"$(INTDIR)/dpttrf.obj" \
	"$(INTDIR)/dpotf2.obj" \
	"$(INTDIR)/dgbsv.obj" \
	"$(INTDIR)/dlacon.obj" \
	"$(INTDIR)/dlarfx.obj" \
	"$(INTDIR)/lsame.obj" \
	"$(INTDIR)/dtrexc.obj" \
	"$(INTDIR)/dgbtrf.obj" \
	"$(INTDIR)/dlapy2.obj" \
	"$(INTDIR)/lsamen.obj" \
	"$(INTDIR)/dpttrs.obj" \
	"$(INTDIR)/dsteqr.obj" \
	"$(INTDIR)/dtrsyl.obj" \
	"$(INTDIR)/dgbtrs.obj" \
	"$(INTDIR)/dlarnd.obj" \
	"$(INTDIR)/icmax1.obj" \
	"$(INTDIR)/dlarnv.obj" \
	"$(INTDIR)/dlaev2.obj" \
	"$(INTDIR)/dlascl.obj" \
	"$(INTDIR)/dlaruv.obj" \
	"$(INTDIR)/dlasrt.obj" \
	"$(INTDIR)/dlagtm.obj" \
	"$(INTDIR)/dlarf.obj" \
	"$(INTDIR)/dpbtf2.obj" \
	"$(INTDIR)/dlacpy.obj" \
	"$(INTDIR)/dlahqr.obj" \
	"$(INTDIR)/dtrsen.obj" \
	"$(INTDIR)/dlae2.obj" \
	"$(INTDIR)/dlasr.obj" \
	"$(INTDIR)/dgeqr2.obj" \
	".\lib\BLAS.lib"

"$(OUTDIR)\LAPACK.lib" : "$(OUTDIR)" $(DEF_FILE) $(LIB32_OBJS)
    $(LIB32) @<<
  $(LIB32_FLAGS) $(DEF_FLAGS) $(LIB32_OBJS)
<<

!ELSEIF  "$(CFG)" == "LAPACK - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "LAPACK\Debug"
# PROP BASE Intermediate_Dir "LAPACK\Debug"
# PROP BASE Target_Dir "LAPACK"
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir ".\lib"
# PROP Intermediate_Dir ".\obj\fortran"
# PROP Target_Dir "LAPACK"
OUTDIR=.\.\lib
INTDIR=.\.\obj\fortran

ALL : "BLAS - Win32 Debug" "$(OUTDIR)\LAPACK.lib"

CLEAN : 
	-@erase ".\lib\LAPACK.lib"
	-@erase ".\obj\fortran\dlabad.obj"
	-@erase ".\obj\fortran\dpbsv.obj"
	-@erase ".\obj\fortran\dlanv2.obj"
	-@erase ".\obj\fortran\dlaptm.obj"
	-@erase ".\obj\fortran\dgetf2.obj"
	-@erase ".\obj\fortran\dgttrf.obj"
	-@erase ".\obj\fortran\dpbtrf.obj"
	-@erase ".\obj\fortran\dlanhs.obj"
	-@erase ".\obj\fortran\dtrevc.obj"
	-@erase ".\obj\fortran\dlaran.obj"
	-@erase ".\obj\fortran\dlassq.obj"
	-@erase ".\obj\fortran\dgttrs.obj"
	-@erase ".\obj\fortran\dlamch.obj"
	-@erase ".\obj\fortran\ilaenv.obj"
	-@erase ".\obj\fortran\dpbtrs.obj"
	-@erase ".\obj\fortran\dladiv.obj"
	-@erase ".\obj\fortran\dlaln2.obj"
	-@erase ".\obj\fortran\izmax1.obj"
	-@erase ".\obj\fortran\dgesv.obj"
	-@erase ".\obj\fortran\dgetrf.obj"
	-@erase ".\obj\fortran\dlaswp.obj"
	-@erase ".\obj\fortran\dorm2r.obj"
	-@erase ".\obj\fortran\dlarfg.obj"
	-@erase ".\obj\fortran\dlanst.obj"
	-@erase ".\obj\fortran\dlasy2.obj"
	-@erase ".\obj\fortran\xerbla.obj"
	-@erase ".\obj\fortran\dlapy3.obj"
	-@erase ".\obj\fortran\dgbtf2.obj"
	-@erase ".\obj\fortran\dlartg.obj"
	-@erase ".\obj\fortran\dgetrs.obj"
	-@erase ".\obj\fortran\dzsum1.obj"
	-@erase ".\obj\fortran\dlange.obj"
	-@erase ".\obj\fortran\dlaset.obj"
	-@erase ".\obj\fortran\dlaexc.obj"
	-@erase ".\obj\fortran\xlaenv.obj"
	-@erase ".\obj\fortran\dpttrf.obj"
	-@erase ".\obj\fortran\dpotf2.obj"
	-@erase ".\obj\fortran\dgbsv.obj"
	-@erase ".\obj\fortran\dlacon.obj"
	-@erase ".\obj\fortran\dlarfx.obj"
	-@erase ".\obj\fortran\lsame.obj"
	-@erase ".\obj\fortran\dtrexc.obj"
	-@erase ".\obj\fortran\dgbtrf.obj"
	-@erase ".\obj\fortran\dlapy2.obj"
	-@erase ".\obj\fortran\lsamen.obj"
	-@erase ".\obj\fortran\dpttrs.obj"
	-@erase ".\obj\fortran\dsteqr.obj"
	-@erase ".\obj\fortran\dtrsyl.obj"
	-@erase ".\obj\fortran\dgbtrs.obj"
	-@erase ".\obj\fortran\dlarnd.obj"
	-@erase ".\obj\fortran\icmax1.obj"
	-@erase ".\obj\fortran\dlarnv.obj"
	-@erase ".\obj\fortran\dlaev2.obj"
	-@erase ".\obj\fortran\dlascl.obj"
	-@erase ".\obj\fortran\dlaruv.obj"
	-@erase ".\obj\fortran\dlasrt.obj"
	-@erase ".\obj\fortran\dlagtm.obj"
	-@erase ".\obj\fortran\dlarf.obj"
	-@erase ".\obj\fortran\dpbtf2.obj"
	-@erase ".\obj\fortran\dlacpy.obj"
	-@erase ".\obj\fortran\dlahqr.obj"
	-@erase ".\obj\fortran\dtrsen.obj"
	-@erase ".\obj\fortran\dlae2.obj"
	-@erase ".\obj\fortran\dlasr.obj"
	-@erase ".\obj\fortran\dgeqr2.obj"

"$(OUTDIR)" :
    if not exist "$(OUTDIR)/$(NULL)" mkdir "$(OUTDIR)"

"$(INTDIR)" :
    if not exist "$(INTDIR)/$(NULL)" mkdir "$(INTDIR)"

# ADD BASE F90 /Z7 /I "LAPACK\Debug/" /c /nologo
# ADD F90 /Z7 /I ".\obj\fortran/" /c /nologo
F90_PROJ=/Z7 /I ".\obj\fortran/" /c /nologo /Fo".\obj\fortran/" 
F90_OBJS=.\.\obj\fortran/
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
BSC32_FLAGS=/nologo /o"$(OUTDIR)/LAPACK.bsc" 
BSC32_SBRS=
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo
LIB32_FLAGS=/nologo /out:"$(OUTDIR)/LAPACK.lib" 
LIB32_OBJS= \
	"$(INTDIR)/dlabad.obj" \
	"$(INTDIR)/dpbsv.obj" \
	"$(INTDIR)/dlanv2.obj" \
	"$(INTDIR)/dlaptm.obj" \
	"$(INTDIR)/dgetf2.obj" \
	"$(INTDIR)/dgttrf.obj" \
	"$(INTDIR)/dpbtrf.obj" \
	"$(INTDIR)/dlanhs.obj" \
	"$(INTDIR)/dtrevc.obj" \
	"$(INTDIR)/dlaran.obj" \
	"$(INTDIR)/dlassq.obj" \
	"$(INTDIR)/dgttrs.obj" \
	"$(INTDIR)/dlamch.obj" \
	"$(INTDIR)/ilaenv.obj" \
	"$(INTDIR)/dpbtrs.obj" \
	"$(INTDIR)/dladiv.obj" \
	"$(INTDIR)/dlaln2.obj" \
	"$(INTDIR)/izmax1.obj" \
	"$(INTDIR)/dgesv.obj" \
	"$(INTDIR)/dgetrf.obj" \
	"$(INTDIR)/dlaswp.obj" \
	"$(INTDIR)/dorm2r.obj" \
	"$(INTDIR)/dlarfg.obj" \
	"$(INTDIR)/dlanst.obj" \
	"$(INTDIR)/dlasy2.obj" \
	"$(INTDIR)/xerbla.obj" \
	"$(INTDIR)/dlapy3.obj" \
	"$(INTDIR)/dgbtf2.obj" \
	"$(INTDIR)/dlartg.obj" \
	"$(INTDIR)/dgetrs.obj" \
	"$(INTDIR)/dzsum1.obj" \
	"$(INTDIR)/dlange.obj" \
	"$(INTDIR)/dlaset.obj" \
	"$(INTDIR)/dlaexc.obj" \
	"$(INTDIR)/xlaenv.obj" \
	"$(INTDIR)/dpttrf.obj" \
	"$(INTDIR)/dpotf2.obj" \
	"$(INTDIR)/dgbsv.obj" \
	"$(INTDIR)/dlacon.obj" \
	"$(INTDIR)/dlarfx.obj" \
	"$(INTDIR)/lsame.obj" \
	"$(INTDIR)/dtrexc.obj" \
	"$(INTDIR)/dgbtrf.obj" \
	"$(INTDIR)/dlapy2.obj" \
	"$(INTDIR)/lsamen.obj" \
	"$(INTDIR)/dpttrs.obj" \
	"$(INTDIR)/dsteqr.obj" \
	"$(INTDIR)/dtrsyl.obj" \
	"$(INTDIR)/dgbtrs.obj" \
	"$(INTDIR)/dlarnd.obj" \
	"$(INTDIR)/icmax1.obj" \
	"$(INTDIR)/dlarnv.obj" \
	"$(INTDIR)/dlaev2.obj" \
	"$(INTDIR)/dlascl.obj" \
	"$(INTDIR)/dlaruv.obj" \
	"$(INTDIR)/dlasrt.obj" \
	"$(INTDIR)/dlagtm.obj" \
	"$(INTDIR)/dlarf.obj" \
	"$(INTDIR)/dpbtf2.obj" \
	"$(INTDIR)/dlacpy.obj" \
	"$(INTDIR)/dlahqr.obj" \
	"$(INTDIR)/dtrsen.obj" \
	"$(INTDIR)/dlae2.obj" \
	"$(INTDIR)/dlasr.obj" \
	"$(INTDIR)/dgeqr2.obj" \
	".\lib\BLAS.lib"

"$(OUTDIR)\LAPACK.lib" : "$(OUTDIR)" $(DEF_FILE) $(LIB32_OBJS)
    $(LIB32) @<<
  $(LIB32_FLAGS) $(DEF_FLAGS) $(LIB32_OBJS)
<<

!ELSEIF  "$(CFG)" == "BLAS - Win32 Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "BLAS\Release"
# PROP BASE Intermediate_Dir "BLAS\Release"
# PROP BASE Target_Dir "BLAS"
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir ".\lib"
# PROP Intermediate_Dir ".\obj\fortran"
# PROP Target_Dir "BLAS"
OUTDIR=.\.\lib
INTDIR=.\.\obj\fortran

ALL : "$(OUTDIR)\BLAS.lib"

CLEAN : 
	-@erase ".\lib\BLAS.lib"
	-@erase ".\obj\fortran\dgemv.obj"
	-@erase ".\obj\fortran\dsymv.obj"
	-@erase ".\obj\fortran\dsyrk.obj"
	-@erase ".\obj\fortran\dsyr.obj"
	-@erase ".\obj\fortran\dtrsm.obj"
	-@erase ".\obj\fortran\dtrsv.obj"
	-@erase ".\obj\fortran\idamax.obj"
	-@erase ".\obj\fortran\dcopy.obj"
	-@erase ".\obj\fortran\dtbsv.obj"
	-@erase ".\obj\fortran\dgbmv.obj"
	-@erase ".\obj\fortran\dger.obj"
	-@erase ".\obj\fortran\ddot.obj"
	-@erase ".\obj\fortran\daxpy.obj"
	-@erase ".\obj\fortran\dsyr2.obj"
	-@erase ".\obj\fortran\dzasum.obj"
	-@erase ".\obj\fortran\dswap.obj"
	-@erase ".\obj\fortran\drot.obj"
	-@erase ".\obj\fortran\dscal.obj"
	-@erase ".\obj\fortran\dznrm2.obj"
	-@erase ".\obj\fortran\dtrmm.obj"
	-@erase ".\obj\fortran\dasum.obj"
	-@erase ".\obj\fortran\drotg.obj"
	-@erase ".\obj\fortran\isamax.obj"
	-@erase ".\obj\fortran\icamax.obj"
	-@erase ".\obj\fortran\dnrm2.obj"
	-@erase ".\obj\fortran\izamax.obj"
	-@erase ".\obj\fortran\dgemm.obj"

"$(OUTDIR)" :
    if not exist "$(OUTDIR)/$(NULL)" mkdir "$(OUTDIR)"

"$(INTDIR)" :
    if not exist "$(INTDIR)/$(NULL)" mkdir "$(INTDIR)"

# ADD BASE F90 /Ox /I "BLAS\Release/" /c /nologo
# ADD F90 /Ox /I ".\obj\fortran/" /c /nologo
F90_PROJ=/Ox /I ".\obj\fortran/" /c /nologo /Fo".\obj\fortran/" 
F90_OBJS=.\.\obj\fortran/
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
BSC32_FLAGS=/nologo /o"$(OUTDIR)/BLAS.bsc" 
BSC32_SBRS=
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo
LIB32_FLAGS=/nologo /out:"$(OUTDIR)/BLAS.lib" 
LIB32_OBJS= \
	"$(INTDIR)/dgemv.obj" \
	"$(INTDIR)/dsymv.obj" \
	"$(INTDIR)/dsyrk.obj" \
	"$(INTDIR)/dsyr.obj" \
	"$(INTDIR)/dtrsm.obj" \
	"$(INTDIR)/dtrsv.obj" \
	"$(INTDIR)/idamax.obj" \
	"$(INTDIR)/dcopy.obj" \
	"$(INTDIR)/dtbsv.obj" \
	"$(INTDIR)/dgbmv.obj" \
	"$(INTDIR)/dger.obj" \
	"$(INTDIR)/ddot.obj" \
	"$(INTDIR)/daxpy.obj" \
	"$(INTDIR)/dsyr2.obj" \
	"$(INTDIR)/dzasum.obj" \
	"$(INTDIR)/dswap.obj" \
	"$(INTDIR)/drot.obj" \
	"$(INTDIR)/dscal.obj" \
	"$(INTDIR)/dznrm2.obj" \
	"$(INTDIR)/dtrmm.obj" \
	"$(INTDIR)/dasum.obj" \
	"$(INTDIR)/drotg.obj" \
	"$(INTDIR)/isamax.obj" \
	"$(INTDIR)/icamax.obj" \
	"$(INTDIR)/dnrm2.obj" \
	"$(INTDIR)/izamax.obj" \
	"$(INTDIR)/dgemm.obj"

"$(OUTDIR)\BLAS.lib" : "$(OUTDIR)" $(DEF_FILE) $(LIB32_OBJS)
    $(LIB32) @<<
  $(LIB32_FLAGS) $(DEF_FLAGS) $(LIB32_OBJS)
<<

!ELSEIF  "$(CFG)" == "BLAS - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "BLAS\Debug"
# PROP BASE Intermediate_Dir "BLAS\Debug"
# PROP BASE Target_Dir "BLAS"
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir ".\lib"
# PROP Intermediate_Dir ".\obj\fortran"
# PROP Target_Dir "BLAS"
OUTDIR=.\.\lib
INTDIR=.\.\obj\fortran

ALL : "$(OUTDIR)\BLAS.lib"

CLEAN : 
	-@erase ".\lib\BLAS.lib"
	-@erase ".\obj\fortran\dgemv.obj"
	-@erase ".\obj\fortran\dsymv.obj"
	-@erase ".\obj\fortran\dsyrk.obj"
	-@erase ".\obj\fortran\dsyr.obj"
	-@erase ".\obj\fortran\dtrsm.obj"
	-@erase ".\obj\fortran\dtrsv.obj"
	-@erase ".\obj\fortran\idamax.obj"
	-@erase ".\obj\fortran\dcopy.obj"
	-@erase ".\obj\fortran\dtbsv.obj"
	-@erase ".\obj\fortran\dgbmv.obj"
	-@erase ".\obj\fortran\dger.obj"
	-@erase ".\obj\fortran\ddot.obj"
	-@erase ".\obj\fortran\daxpy.obj"
	-@erase ".\obj\fortran\dsyr2.obj"
	-@erase ".\obj\fortran\dzasum.obj"
	-@erase ".\obj\fortran\dswap.obj"
	-@erase ".\obj\fortran\drot.obj"
	-@erase ".\obj\fortran\dscal.obj"
	-@erase ".\obj\fortran\dznrm2.obj"
	-@erase ".\obj\fortran\dtrmm.obj"
	-@erase ".\obj\fortran\dasum.obj"
	-@erase ".\obj\fortran\drotg.obj"
	-@erase ".\obj\fortran\isamax.obj"
	-@erase ".\obj\fortran\icamax.obj"
	-@erase ".\obj\fortran\dnrm2.obj"
	-@erase ".\obj\fortran\izamax.obj"
	-@erase ".\obj\fortran\dgemm.obj"

"$(OUTDIR)" :
    if not exist "$(OUTDIR)/$(NULL)" mkdir "$(OUTDIR)"

"$(INTDIR)" :
    if not exist "$(INTDIR)/$(NULL)" mkdir "$(INTDIR)"

# ADD BASE F90 /Z7 /I "BLAS\Debug/" /c /nologo
# ADD F90 /Z7 /I ".\obj\fortran/" /c /nologo
F90_PROJ=/Z7 /I ".\obj\fortran/" /c /nologo /Fo".\obj\fortran/" 
F90_OBJS=.\.\obj\fortran/
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
BSC32_FLAGS=/nologo /o"$(OUTDIR)/BLAS.bsc" 
BSC32_SBRS=
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo
LIB32_FLAGS=/nologo /out:"$(OUTDIR)/BLAS.lib" 
LIB32_OBJS= \
	"$(INTDIR)/dgemv.obj" \
	"$(INTDIR)/dsymv.obj" \
	"$(INTDIR)/dsyrk.obj" \
	"$(INTDIR)/dsyr.obj" \
	"$(INTDIR)/dtrsm.obj" \
	"$(INTDIR)/dtrsv.obj" \
	"$(INTDIR)/idamax.obj" \
	"$(INTDIR)/dcopy.obj" \
	"$(INTDIR)/dtbsv.obj" \
	"$(INTDIR)/dgbmv.obj" \
	"$(INTDIR)/dger.obj" \
	"$(INTDIR)/ddot.obj" \
	"$(INTDIR)/daxpy.obj" \
	"$(INTDIR)/dsyr2.obj" \
	"$(INTDIR)/dzasum.obj" \
	"$(INTDIR)/dswap.obj" \
	"$(INTDIR)/drot.obj" \
	"$(INTDIR)/dscal.obj" \
	"$(INTDIR)/dznrm2.obj" \
	"$(INTDIR)/dtrmm.obj" \
	"$(INTDIR)/dasum.obj" \
	"$(INTDIR)/drotg.obj" \
	"$(INTDIR)/isamax.obj" \
	"$(INTDIR)/icamax.obj" \
	"$(INTDIR)/dnrm2.obj" \
	"$(INTDIR)/izamax.obj" \
	"$(INTDIR)/dgemm.obj"

"$(OUTDIR)\BLAS.lib" : "$(OUTDIR)" $(DEF_FILE) $(LIB32_OBJS)
    $(LIB32) @<<
  $(LIB32_FLAGS) $(DEF_FLAGS) $(LIB32_OBJS)
<<

!ELSEIF  "$(CFG)" == "UmfPack - Win32 Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "UmfPack\Release"
# PROP BASE Intermediate_Dir "UmfPack\Release"
# PROP BASE Target_Dir "UmfPack"
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir ".\lib"
# PROP Intermediate_Dir ".\obj\fortran"
# PROP Target_Dir "UmfPack"
OUTDIR=.\.\lib
INTDIR=.\.\obj\fortran

ALL : "BLAS - Win32 Release" "$(OUTDIR)\UmfPack.lib"

CLEAN : 
	-@erase ".\lib\UmfPack.lib"
	-@erase ".\obj\fortran\umd2p2.obj"
	-@erase ".\obj\fortran\umd2rf.obj"
	-@erase ".\obj\fortran\umd2er.obj"
	-@erase ".\obj\fortran\umd2f1.obj"
	-@erase ".\obj\fortran\umd21i.obj"
	-@erase ".\obj\fortran\umd2ut.obj"
	-@erase ".\obj\fortran\umd2ra.obj"
	-@erase ".\obj\fortran\umd2s2.obj"
	-@erase ".\obj\fortran\umd2fg.obj"
	-@erase ".\obj\fortran\umd2p1.obj"
	-@erase ".\obj\fortran\umd2f0.obj"
	-@erase ".\obj\fortran\umd2co.obj"
	-@erase ".\obj\fortran\umd2fb.obj"
	-@erase ".\obj\fortran\mc21b.obj"
	-@erase ".\obj\fortran\mc13e.obj"
	-@erase ".\obj\fortran\umd2sl.obj"
	-@erase ".\obj\fortran\umd2su.obj"
	-@erase ".\obj\fortran\umd2r2.obj"
	-@erase ".\obj\fortran\umd2fa.obj"
	-@erase ".\obj\fortran\umd2lt.obj"
	-@erase ".\obj\fortran\umd2so.obj"
	-@erase ".\obj\fortran\umd2in.obj"
	-@erase ".\obj\fortran\umd2rg.obj"
	-@erase ".\obj\fortran\umd2f2.obj"
	-@erase ".\obj\fortran\umd2of.obj"
	-@erase ".\obj\fortran\umd2r0.obj"

"$(OUTDIR)" :
    if not exist "$(OUTDIR)/$(NULL)" mkdir "$(OUTDIR)"

"$(INTDIR)" :
    if not exist "$(INTDIR)/$(NULL)" mkdir "$(INTDIR)"

# ADD BASE F90 /Ox /I "UmfPack\Release/" /c /nologo
# ADD F90 /Ox /I ".\obj\fortran/" /c /nologo
F90_PROJ=/Ox /I ".\obj\fortran/" /c /nologo /Fo".\obj\fortran/" 
F90_OBJS=.\.\obj\fortran/
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
BSC32_FLAGS=/nologo /o"$(OUTDIR)/UmfPack.bsc" 
BSC32_SBRS=
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo
LIB32_FLAGS=/nologo /out:"$(OUTDIR)/UmfPack.lib" 
LIB32_OBJS= \
	"$(INTDIR)/umd2p2.obj" \
	"$(INTDIR)/umd2rf.obj" \
	"$(INTDIR)/umd2er.obj" \
	"$(INTDIR)/umd2f1.obj" \
	"$(INTDIR)/umd21i.obj" \
	"$(INTDIR)/umd2ut.obj" \
	"$(INTDIR)/umd2ra.obj" \
	"$(INTDIR)/umd2s2.obj" \
	"$(INTDIR)/umd2fg.obj" \
	"$(INTDIR)/umd2p1.obj" \
	"$(INTDIR)/umd2f0.obj" \
	"$(INTDIR)/umd2co.obj" \
	"$(INTDIR)/umd2fb.obj" \
	"$(INTDIR)/mc21b.obj" \
	"$(INTDIR)/mc13e.obj" \
	"$(INTDIR)/umd2sl.obj" \
	"$(INTDIR)/umd2su.obj" \
	"$(INTDIR)/umd2r2.obj" \
	"$(INTDIR)/umd2fa.obj" \
	"$(INTDIR)/umd2lt.obj" \
	"$(INTDIR)/umd2so.obj" \
	"$(INTDIR)/umd2in.obj" \
	"$(INTDIR)/umd2rg.obj" \
	"$(INTDIR)/umd2f2.obj" \
	"$(INTDIR)/umd2of.obj" \
	"$(INTDIR)/umd2r0.obj" \
	".\lib\BLAS.lib"

"$(OUTDIR)\UmfPack.lib" : "$(OUTDIR)" $(DEF_FILE) $(LIB32_OBJS)
    $(LIB32) @<<
  $(LIB32_FLAGS) $(DEF_FLAGS) $(LIB32_OBJS)
<<

!ELSEIF  "$(CFG)" == "UmfPack - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "UmfPack\Debug"
# PROP BASE Intermediate_Dir "UmfPack\Debug"
# PROP BASE Target_Dir "UmfPack"
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir ".\lib"
# PROP Intermediate_Dir ".\obj\fortran"
# PROP Target_Dir "UmfPack"
OUTDIR=.\.\lib
INTDIR=.\.\obj\fortran

ALL : "BLAS - Win32 Debug" "$(OUTDIR)\UmfPack.lib"

CLEAN : 
	-@erase ".\lib\UmfPack.lib"
	-@erase ".\obj\fortran\umd2p2.obj"
	-@erase ".\obj\fortran\umd2rf.obj"
	-@erase ".\obj\fortran\umd2er.obj"
	-@erase ".\obj\fortran\umd2f1.obj"
	-@erase ".\obj\fortran\umd21i.obj"
	-@erase ".\obj\fortran\umd2ut.obj"
	-@erase ".\obj\fortran\umd2ra.obj"
	-@erase ".\obj\fortran\umd2s2.obj"
	-@erase ".\obj\fortran\umd2fg.obj"
	-@erase ".\obj\fortran\umd2p1.obj"
	-@erase ".\obj\fortran\umd2f0.obj"
	-@erase ".\obj\fortran\umd2co.obj"
	-@erase ".\obj\fortran\umd2fb.obj"
	-@erase ".\obj\fortran\mc21b.obj"
	-@erase ".\obj\fortran\mc13e.obj"
	-@erase ".\obj\fortran\umd2sl.obj"
	-@erase ".\obj\fortran\umd2su.obj"
	-@erase ".\obj\fortran\umd2r2.obj"
	-@erase ".\obj\fortran\umd2fa.obj"
	-@erase ".\obj\fortran\umd2lt.obj"
	-@erase ".\obj\fortran\umd2so.obj"
	-@erase ".\obj\fortran\umd2in.obj"
	-@erase ".\obj\fortran\umd2rg.obj"
	-@erase ".\obj\fortran\umd2f2.obj"
	-@erase ".\obj\fortran\umd2of.obj"
	-@erase ".\obj\fortran\umd2r0.obj"

"$(OUTDIR)" :
    if not exist "$(OUTDIR)/$(NULL)" mkdir "$(OUTDIR)"

"$(INTDIR)" :
    if not exist "$(INTDIR)/$(NULL)" mkdir "$(INTDIR)"

# ADD BASE F90 /Z7 /I "UmfPack\Debug/" /c /nologo
# ADD F90 /Z7 /I ".\obj\fortran/" /c /nologo
F90_PROJ=/Z7 /I ".\obj\fortran/" /c /nologo /Fo".\obj\fortran/" 
F90_OBJS=.\.\obj\fortran/
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
BSC32_FLAGS=/nologo /o"$(OUTDIR)/UmfPack.bsc" 
BSC32_SBRS=
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo
LIB32_FLAGS=/nologo /out:"$(OUTDIR)/UmfPack.lib" 
LIB32_OBJS= \
	"$(INTDIR)/umd2p2.obj" \
	"$(INTDIR)/umd2rf.obj" \
	"$(INTDIR)/umd2er.obj" \
	"$(INTDIR)/umd2f1.obj" \
	"$(INTDIR)/umd21i.obj" \
	"$(INTDIR)/umd2ut.obj" \
	"$(INTDIR)/umd2ra.obj" \
	"$(INTDIR)/umd2s2.obj" \
	"$(INTDIR)/umd2fg.obj" \
	"$(INTDIR)/umd2p1.obj" \
	"$(INTDIR)/umd2f0.obj" \
	"$(INTDIR)/umd2co.obj" \
	"$(INTDIR)/umd2fb.obj" \
	"$(INTDIR)/mc21b.obj" \
	"$(INTDIR)/mc13e.obj" \
	"$(INTDIR)/umd2sl.obj" \
	"$(INTDIR)/umd2su.obj" \
	"$(INTDIR)/umd2r2.obj" \
	"$(INTDIR)/umd2fa.obj" \
	"$(INTDIR)/umd2lt.obj" \
	"$(INTDIR)/umd2so.obj" \
	"$(INTDIR)/umd2in.obj" \
	"$(INTDIR)/umd2rg.obj" \
	"$(INTDIR)/umd2f2.obj" \
	"$(INTDIR)/umd2of.obj" \
	"$(INTDIR)/umd2r0.obj" \
	".\lib\BLAS.lib"

"$(OUTDIR)\UmfPack.lib" : "$(OUTDIR)" $(DEF_FILE) $(LIB32_OBJS)
    $(LIB32) @<<
  $(LIB32_FLAGS) $(DEF_FLAGS) $(LIB32_OBJS)
<<

!ELSEIF  "$(CFG)" == "feap - Win32 Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "feap\Release"
# PROP BASE Intermediate_Dir "feap\Release"
# PROP BASE Target_Dir "feap"
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir ".\lib"
# PROP Intermediate_Dir ".\obj\fortran"
# PROP Target_Dir "feap"
OUTDIR=.\.\lib
INTDIR=.\.\obj\fortran

ALL : "$(OUTDIR)\feap.lib"

CLEAN : 
	-@erase ".\lib\feap.lib"
	-@erase ".\obj\fortran\getCommon.obj"
	-@erase ".\obj\fortran\elmt06.obj"
	-@erase ".\obj\fortran\elmt01.obj"
	-@erase ".\obj\fortran\elmt05.obj"
	-@erase ".\obj\fortran\dummyFeap.obj"
	-@erase ".\obj\fortran\elmt04.obj"
	-@erase ".\obj\fortran\elmt03.obj"
	-@erase ".\obj\fortran\fillCommon.obj"
	-@erase ".\obj\fortran\elmt02.obj"

"$(OUTDIR)" :
    if not exist "$(OUTDIR)/$(NULL)" mkdir "$(OUTDIR)"

"$(INTDIR)" :
    if not exist "$(INTDIR)/$(NULL)" mkdir "$(INTDIR)"

# ADD BASE F90 /Ox /I "feap\Release/" /c /nologo
# ADD F90 /Ox /I ".\obj\fortran/" /c /nologo
F90_PROJ=/Ox /I ".\obj\fortran/" /c /nologo /Fo".\obj\fortran/" 
F90_OBJS=.\.\obj\fortran/
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
BSC32_FLAGS=/nologo /o"$(OUTDIR)/feap.bsc" 
BSC32_SBRS=
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo
LIB32_FLAGS=/nologo /out:"$(OUTDIR)/feap.lib" 
LIB32_OBJS= \
	"$(INTDIR)/getCommon.obj" \
	"$(INTDIR)/elmt06.obj" \
	"$(INTDIR)/elmt01.obj" \
	"$(INTDIR)/elmt05.obj" \
	"$(INTDIR)/dummyFeap.obj" \
	"$(INTDIR)/elmt04.obj" \
	"$(INTDIR)/elmt03.obj" \
	"$(INTDIR)/fillCommon.obj" \
	"$(INTDIR)/elmt02.obj"

"$(OUTDIR)\feap.lib" : "$(OUTDIR)" $(DEF_FILE) $(LIB32_OBJS)
    $(LIB32) @<<
  $(LIB32_FLAGS) $(DEF_FLAGS) $(LIB32_OBJS)
<<

!ELSEIF  "$(CFG)" == "feap - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "feap\Debug"
# PROP BASE Intermediate_Dir "feap\Debug"
# PROP BASE Target_Dir "feap"
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir ".\lib"
# PROP Intermediate_Dir ".\obj\fortran"
# PROP Target_Dir "feap"
OUTDIR=.\.\lib
INTDIR=.\.\obj\fortran

ALL : "$(OUTDIR)\feap.lib"

CLEAN : 
	-@erase ".\lib\feap.lib"
	-@erase ".\obj\fortran\getCommon.obj"
	-@erase ".\obj\fortran\elmt06.obj"
	-@erase ".\obj\fortran\elmt01.obj"
	-@erase ".\obj\fortran\elmt05.obj"
	-@erase ".\obj\fortran\dummyFeap.obj"
	-@erase ".\obj\fortran\elmt04.obj"
	-@erase ".\obj\fortran\elmt03.obj"
	-@erase ".\obj\fortran\fillCommon.obj"
	-@erase ".\obj\fortran\elmt02.obj"

"$(OUTDIR)" :
    if not exist "$(OUTDIR)/$(NULL)" mkdir "$(OUTDIR)"

"$(INTDIR)" :
    if not exist "$(INTDIR)/$(NULL)" mkdir "$(INTDIR)"

# ADD BASE F90 /Z7 /I "feap\Debug/" /c /nologo
# ADD F90 /Z7 /I ".\obj\fortran/" /c /nologo
F90_PROJ=/Z7 /I ".\obj\fortran/" /c /nologo /Fo".\obj\fortran/" 
F90_OBJS=.\.\obj\fortran/
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
BSC32_FLAGS=/nologo /o"$(OUTDIR)/feap.bsc" 
BSC32_SBRS=
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo
LIB32_FLAGS=/nologo /out:"$(OUTDIR)/feap.lib" 
LIB32_OBJS= \
	"$(INTDIR)/getCommon.obj" \
	"$(INTDIR)/elmt06.obj" \
	"$(INTDIR)/elmt01.obj" \
	"$(INTDIR)/elmt05.obj" \
	"$(INTDIR)/dummyFeap.obj" \
	"$(INTDIR)/elmt04.obj" \
	"$(INTDIR)/elmt03.obj" \
	"$(INTDIR)/fillCommon.obj" \
	"$(INTDIR)/elmt02.obj"

"$(OUTDIR)\feap.lib" : "$(OUTDIR)" $(DEF_FILE) $(LIB32_OBJS)
    $(LIB32) @<<
  $(LIB32_FLAGS) $(DEF_FLAGS) $(LIB32_OBJS)
<<

!ELSEIF  "$(CFG)" == "OpenSees - Win32 Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "OpenSees\Release"
# PROP BASE Intermediate_Dir "OpenSees\Release"
# PROP BASE Target_Dir "OpenSees"
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir ".\lib"
# PROP Intermediate_Dir ".\obj\fortran"
# PROP Target_Dir "OpenSees"
OUTDIR=.\.\lib
INTDIR=.\.\obj\fortran

ALL : "$(OUTDIR)\OpenSeesFortran.lib"

CLEAN : 
	-@erase ".\lib\OpenSeesFortran.lib"
	-@erase ".\obj\fortran\genmmd.obj"

"$(OUTDIR)" :
    if not exist "$(OUTDIR)/$(NULL)" mkdir "$(OUTDIR)"

"$(INTDIR)" :
    if not exist "$(INTDIR)/$(NULL)" mkdir "$(INTDIR)"

# ADD BASE F90 /Ox /I "OpenSees\Release/" /c /nologo
# ADD F90 /Ox /I ".\obj\fortran/" /c /nologo
F90_PROJ=/Ox /I ".\obj\fortran/" /c /nologo /Fo".\obj\fortran/" 
F90_OBJS=.\.\obj\fortran/
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
BSC32_FLAGS=/nologo /o"$(OUTDIR)/OpenSees.bsc" 
BSC32_SBRS=
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo /out:".\lib\OpenSeesFortran.lib"
LIB32_FLAGS=/nologo /out:"$(OUTDIR)/OpenSeesFortran.lib" 
LIB32_OBJS= \
	"$(INTDIR)/genmmd.obj"

"$(OUTDIR)\OpenSeesFortran.lib" : "$(OUTDIR)" $(DEF_FILE) $(LIB32_OBJS)
    $(LIB32) @<<
  $(LIB32_FLAGS) $(DEF_FLAGS) $(LIB32_OBJS)
<<

!ELSEIF  "$(CFG)" == "OpenSees - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "OpenSees\Debug"
# PROP BASE Intermediate_Dir "OpenSees\Debug"
# PROP BASE Target_Dir "OpenSees"
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir ".\lib"
# PROP Intermediate_Dir ".\obj\fortran"
# PROP Target_Dir "OpenSees"
OUTDIR=.\.\lib
INTDIR=.\.\obj\fortran

ALL : "$(OUTDIR)\OpenSees.lib"

CLEAN : 
	-@erase ".\lib\OpenSees.lib"
	-@erase ".\obj\fortran\genmmd.obj"

"$(OUTDIR)" :
    if not exist "$(OUTDIR)/$(NULL)" mkdir "$(OUTDIR)"

"$(INTDIR)" :
    if not exist "$(INTDIR)/$(NULL)" mkdir "$(INTDIR)"

# ADD BASE F90 /Z7 /I "OpenSees\Debug/" /c /nologo
# ADD F90 /Z7 /I ".\obj\fortran/" /c /nologo
F90_PROJ=/Z7 /I ".\obj\fortran/" /c /nologo /Fo".\obj\fortran/" 
F90_OBJS=.\.\obj\fortran/
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
BSC32_FLAGS=/nologo /o"$(OUTDIR)/OpenSees.bsc" 
BSC32_SBRS=
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo
LIB32_FLAGS=/nologo /out:"$(OUTDIR)/OpenSees.lib" 
LIB32_OBJS= \
	"$(INTDIR)/genmmd.obj"

"$(OUTDIR)\OpenSees.lib" : "$(OUTDIR)" $(DEF_FILE) $(LIB32_OBJS)
    $(LIB32) @<<
  $(LIB32_FLAGS) $(DEF_FLAGS) $(LIB32_OBJS)
<<

!ENDIF 

.for{$(F90_OBJS)}.obj:
   $(F90) $(F90_PROJ) $<  

.f{$(F90_OBJS)}.obj:
   $(F90) $(F90_PROJ) $<  

.f90{$(F90_OBJS)}.obj:
   $(F90) $(F90_PROJ) $<  

################################################################################
# Begin Target

# Name "g3Fortran - Win32 Release"
# Name "g3Fortran - Win32 Debug"

!IF  "$(CFG)" == "g3Fortran - Win32 Release"

!ELSEIF  "$(CFG)" == "g3Fortran - Win32 Debug"

!ENDIF 

# End Target
################################################################################
# Begin Target

# Name "ARPACK - Win32 Release"
# Name "ARPACK - Win32 Debug"

!IF  "$(CFG)" == "ARPACK - Win32 Release"

!ELSEIF  "$(CFG)" == "ARPACK - Win32 Debug"

!ENDIF 

################################################################################
# Begin Project Dependency

# Project_Dep_Name "BLAS"

!IF  "$(CFG)" == "ARPACK - Win32 Release"

"BLAS - Win32 Release" : 
   $(MAKE) /$(MAKEFLAGS) /F .\OpenSeesFortran.mak CFG="BLAS - Win32 Release" 

!ELSEIF  "$(CFG)" == "ARPACK - Win32 Debug"

"BLAS - Win32 Debug" : 
   $(MAKE) /$(MAKEFLAGS) /F .\OpenSeesFortran.mak CFG="BLAS - Win32 Debug" 

!ENDIF 

# End Project Dependency
################################################################################
# Begin Source File

SOURCE=\g3\OTHER\ARPACK\zvout.f

"$(INTDIR)\zvout.obj" : $(SOURCE) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


# End Source File
################################################################################
# Begin Source File

SOURCE=\g3\OTHER\ARPACK\zmout.f

"$(INTDIR)\zmout.obj" : $(SOURCE) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


# End Source File
################################################################################
# Begin Source File

SOURCE=\g3\OTHER\ARPACK\svout.f

"$(INTDIR)\svout.obj" : $(SOURCE) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


# End Source File
################################################################################
# Begin Source File

SOURCE=\g3\OTHER\ARPACK\smout.f

"$(INTDIR)\smout.obj" : $(SOURCE) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


# End Source File
################################################################################
# Begin Source File

SOURCE=\g3\OTHER\ARPACK\second.f

"$(INTDIR)\second.obj" : $(SOURCE) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


# End Source File
################################################################################
# Begin Source File

SOURCE=\g3\OTHER\ARPACK\ivout.f

"$(INTDIR)\ivout.obj" : $(SOURCE) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


# End Source File
################################################################################
# Begin Source File

SOURCE=\g3\OTHER\ARPACK\iswap.f

"$(INTDIR)\iswap.obj" : $(SOURCE) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


# End Source File
################################################################################
# Begin Source File

SOURCE=\g3\OTHER\ARPACK\iset.f

"$(INTDIR)\iset.obj" : $(SOURCE) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


# End Source File
################################################################################
# Begin Source File

SOURCE=\g3\OTHER\ARPACK\icopy.f

"$(INTDIR)\icopy.obj" : $(SOURCE) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


# End Source File
################################################################################
# Begin Source File

SOURCE=\g3\OTHER\ARPACK\icnteq.f

"$(INTDIR)\icnteq.obj" : $(SOURCE) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


# End Source File
################################################################################
# Begin Source File

SOURCE=\g3\OTHER\ARPACK\dvout.f

"$(INTDIR)\dvout.obj" : $(SOURCE) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


# End Source File
################################################################################
# Begin Source File

SOURCE=\g3\OTHER\ARPACK\dstqrb.f

"$(INTDIR)\dstqrb.obj" : $(SOURCE) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


# End Source File
################################################################################
# Begin Source File

SOURCE=\g3\OTHER\ARPACK\dstats.f
DEP_F90_DSTAT=\
	".\..\OTHER\ARPACK\stat.h"\
	

"$(INTDIR)\dstats.obj" : $(SOURCE) $(DEP_F90_DSTAT) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


# End Source File
################################################################################
# Begin Source File

SOURCE=\g3\OTHER\ARPACK\dstatn.f
DEP_F90_DSTATN=\
	".\..\OTHER\ARPACK\stat.h"\
	

"$(INTDIR)\dstatn.obj" : $(SOURCE) $(DEP_F90_DSTATN) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


# End Source File
################################################################################
# Begin Source File

SOURCE=\g3\OTHER\ARPACK\dsortr.f

"$(INTDIR)\dsortr.obj" : $(SOURCE) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


# End Source File
################################################################################
# Begin Source File

SOURCE=\g3\OTHER\ARPACK\dsortc.f

"$(INTDIR)\dsortc.obj" : $(SOURCE) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


# End Source File
################################################################################
# Begin Source File

SOURCE=\g3\OTHER\ARPACK\dsgets.f
DEP_F90_DSGET=\
	".\..\OTHER\ARPACK\debug.h"\
	".\..\OTHER\ARPACK\stat.h"\
	

"$(INTDIR)\dsgets.obj" : $(SOURCE) $(DEP_F90_DSGET) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


# End Source File
################################################################################
# Begin Source File

SOURCE=\g3\OTHER\ARPACK\dseupd.f
DEP_F90_DSEUP=\
	".\..\OTHER\ARPACK\debug.h"\
	".\..\OTHER\ARPACK\stat.h"\
	

"$(INTDIR)\dseupd.obj" : $(SOURCE) $(DEP_F90_DSEUP) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


# End Source File
################################################################################
# Begin Source File

SOURCE=\g3\OTHER\ARPACK\dsesrt.f

"$(INTDIR)\dsesrt.obj" : $(SOURCE) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


# End Source File
################################################################################
# Begin Source File

SOURCE=\g3\OTHER\ARPACK\dseigt.f
DEP_F90_DSEIG=\
	".\..\OTHER\ARPACK\debug.h"\
	".\..\OTHER\ARPACK\stat.h"\
	

"$(INTDIR)\dseigt.obj" : $(SOURCE) $(DEP_F90_DSEIG) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


# End Source File
################################################################################
# Begin Source File

SOURCE=\g3\OTHER\ARPACK\dsconv.f
DEP_F90_DSCON=\
	".\..\OTHER\ARPACK\debug.h"\
	".\..\OTHER\ARPACK\stat.h"\
	

"$(INTDIR)\dsconv.obj" : $(SOURCE) $(DEP_F90_DSCON) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


# End Source File
################################################################################
# Begin Source File

SOURCE=\g3\OTHER\ARPACK\dsaupd.f
DEP_F90_DSAUP=\
	".\..\OTHER\ARPACK\debug.h"\
	".\..\OTHER\ARPACK\stat.h"\
	

"$(INTDIR)\dsaupd.obj" : $(SOURCE) $(DEP_F90_DSAUP) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


# End Source File
################################################################################
# Begin Source File

SOURCE=\g3\OTHER\ARPACK\dsaup2.f
DEP_F90_DSAUP2=\
	".\..\OTHER\ARPACK\debug.h"\
	".\..\OTHER\ARPACK\stat.h"\
	

"$(INTDIR)\dsaup2.obj" : $(SOURCE) $(DEP_F90_DSAUP2) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


# End Source File
################################################################################
# Begin Source File

SOURCE=\g3\OTHER\ARPACK\dsapps.f
DEP_F90_DSAPP=\
	".\..\OTHER\ARPACK\debug.h"\
	".\..\OTHER\ARPACK\stat.h"\
	

"$(INTDIR)\dsapps.obj" : $(SOURCE) $(DEP_F90_DSAPP) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


# End Source File
################################################################################
# Begin Source File

SOURCE=\g3\OTHER\ARPACK\dsaitr.f
DEP_F90_DSAIT=\
	".\..\OTHER\ARPACK\debug.h"\
	".\..\OTHER\ARPACK\stat.h"\
	

"$(INTDIR)\dsaitr.obj" : $(SOURCE) $(DEP_F90_DSAIT) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


# End Source File
################################################################################
# Begin Source File

SOURCE=\g3\OTHER\ARPACK\dngets.f
DEP_F90_DNGET=\
	".\..\OTHER\ARPACK\debug.h"\
	".\..\OTHER\ARPACK\stat.h"\
	

"$(INTDIR)\dngets.obj" : $(SOURCE) $(DEP_F90_DNGET) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


# End Source File
################################################################################
# Begin Source File

SOURCE=\g3\OTHER\ARPACK\dneupd.f
DEP_F90_DNEUP=\
	".\..\OTHER\ARPACK\debug.h"\
	".\..\OTHER\ARPACK\stat.h"\
	

"$(INTDIR)\dneupd.obj" : $(SOURCE) $(DEP_F90_DNEUP) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


# End Source File
################################################################################
# Begin Source File

SOURCE=\g3\OTHER\ARPACK\dneigh.f
DEP_F90_DNEIG=\
	".\..\OTHER\ARPACK\debug.h"\
	".\..\OTHER\ARPACK\stat.h"\
	

"$(INTDIR)\dneigh.obj" : $(SOURCE) $(DEP_F90_DNEIG) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


# End Source File
################################################################################
# Begin Source File

SOURCE=\g3\OTHER\ARPACK\dnconv.f
DEP_F90_DNCON=\
	".\..\OTHER\ARPACK\debug.h"\
	".\..\OTHER\ARPACK\stat.h"\
	

"$(INTDIR)\dnconv.obj" : $(SOURCE) $(DEP_F90_DNCON) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


# End Source File
################################################################################
# Begin Source File

SOURCE=\g3\OTHER\ARPACK\dnaupd.f
DEP_F90_DNAUP=\
	".\..\OTHER\ARPACK\debug.h"\
	".\..\OTHER\ARPACK\stat.h"\
	

"$(INTDIR)\dnaupd.obj" : $(SOURCE) $(DEP_F90_DNAUP) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


# End Source File
################################################################################
# Begin Source File

SOURCE=\g3\OTHER\ARPACK\dnaup2.f
DEP_F90_DNAUP2=\
	".\..\OTHER\ARPACK\debug.h"\
	".\..\OTHER\ARPACK\stat.h"\
	

"$(INTDIR)\dnaup2.obj" : $(SOURCE) $(DEP_F90_DNAUP2) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


# End Source File
################################################################################
# Begin Source File

SOURCE=\g3\OTHER\ARPACK\dnapps.f
DEP_F90_DNAPP=\
	".\..\OTHER\ARPACK\debug.h"\
	".\..\OTHER\ARPACK\stat.h"\
	

"$(INTDIR)\dnapps.obj" : $(SOURCE) $(DEP_F90_DNAPP) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


# End Source File
################################################################################
# Begin Source File

SOURCE=\g3\OTHER\ARPACK\dnaitr.f
DEP_F90_DNAIT=\
	".\..\OTHER\ARPACK\debug.h"\
	".\..\OTHER\ARPACK\stat.h"\
	

"$(INTDIR)\dnaitr.obj" : $(SOURCE) $(DEP_F90_DNAIT) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


# End Source File
################################################################################
# Begin Source File

SOURCE=\g3\OTHER\ARPACK\dmout.f

"$(INTDIR)\dmout.obj" : $(SOURCE) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


# End Source File
################################################################################
# Begin Source File

SOURCE=\g3\OTHER\ARPACK\dlaqrb.f

"$(INTDIR)\dlaqrb.obj" : $(SOURCE) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


# End Source File
################################################################################
# Begin Source File

SOURCE=\g3\OTHER\ARPACK\dgetv0.f
DEP_F90_DGETV=\
	".\..\OTHER\ARPACK\debug.h"\
	".\..\OTHER\ARPACK\stat.h"\
	

"$(INTDIR)\dgetv0.obj" : $(SOURCE) $(DEP_F90_DGETV) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


# End Source File
################################################################################
# Begin Source File

SOURCE=\g3\OTHER\ARPACK\cvout.f

"$(INTDIR)\cvout.obj" : $(SOURCE) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


# End Source File
################################################################################
# Begin Source File

SOURCE=\g3\OTHER\ARPACK\cmout.f

"$(INTDIR)\cmout.obj" : $(SOURCE) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


# End Source File
# End Target
################################################################################
# Begin Target

# Name "LAPACK - Win32 Release"
# Name "LAPACK - Win32 Debug"

!IF  "$(CFG)" == "LAPACK - Win32 Release"

!ELSEIF  "$(CFG)" == "LAPACK - Win32 Debug"

!ENDIF 

################################################################################
# Begin Project Dependency

# Project_Dep_Name "BLAS"

!IF  "$(CFG)" == "LAPACK - Win32 Release"

"BLAS - Win32 Release" : 
   $(MAKE) /$(MAKEFLAGS) /F .\OpenSeesFortran.mak CFG="BLAS - Win32 Release" 

!ELSEIF  "$(CFG)" == "LAPACK - Win32 Debug"

"BLAS - Win32 Debug" : 
   $(MAKE) /$(MAKEFLAGS) /F .\OpenSeesFortran.mak CFG="BLAS - Win32 Debug" 

!ENDIF 

# End Project Dependency
################################################################################
# Begin Source File

SOURCE=\g3\OTHER\LAPACK\xlaenv.f

"$(INTDIR)\xlaenv.obj" : $(SOURCE) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


# End Source File
################################################################################
# Begin Source File

SOURCE=\g3\OTHER\LAPACK\xerbla.f

"$(INTDIR)\xerbla.obj" : $(SOURCE) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


# End Source File
################################################################################
# Begin Source File

SOURCE=\g3\OTHER\LAPACK\lsamen.f

"$(INTDIR)\lsamen.obj" : $(SOURCE) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


# End Source File
################################################################################
# Begin Source File

SOURCE=\g3\OTHER\LAPACK\lsame.f

"$(INTDIR)\lsame.obj" : $(SOURCE) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


# End Source File
################################################################################
# Begin Source File

SOURCE=\g3\OTHER\LAPACK\izmax1.f

"$(INTDIR)\izmax1.obj" : $(SOURCE) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


# End Source File
################################################################################
# Begin Source File

SOURCE=\g3\OTHER\LAPACK\ilaenv.f

"$(INTDIR)\ilaenv.obj" : $(SOURCE) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


# End Source File
################################################################################
# Begin Source File

SOURCE=\g3\OTHER\LAPACK\icmax1.f

"$(INTDIR)\icmax1.obj" : $(SOURCE) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


# End Source File
################################################################################
# Begin Source File

SOURCE=\g3\OTHER\LAPACK\dzsum1.f

"$(INTDIR)\dzsum1.obj" : $(SOURCE) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


# End Source File
################################################################################
# Begin Source File

SOURCE=\g3\OTHER\LAPACK\dtrsyl.f

"$(INTDIR)\dtrsyl.obj" : $(SOURCE) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


# End Source File
################################################################################
# Begin Source File

SOURCE=\g3\OTHER\LAPACK\dtrsen.f

"$(INTDIR)\dtrsen.obj" : $(SOURCE) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


# End Source File
################################################################################
# Begin Source File

SOURCE=\g3\OTHER\LAPACK\dtrexc.f

"$(INTDIR)\dtrexc.obj" : $(SOURCE) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


# End Source File
################################################################################
# Begin Source File

SOURCE=\g3\OTHER\LAPACK\dtrevc.f

"$(INTDIR)\dtrevc.obj" : $(SOURCE) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


# End Source File
################################################################################
# Begin Source File

SOURCE=\g3\OTHER\LAPACK\dsteqr.f

"$(INTDIR)\dsteqr.obj" : $(SOURCE) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


# End Source File
################################################################################
# Begin Source File

SOURCE=\g3\OTHER\LAPACK\dpttrs.f

"$(INTDIR)\dpttrs.obj" : $(SOURCE) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


# End Source File
################################################################################
# Begin Source File

SOURCE=\g3\OTHER\LAPACK\dpttrf.f

"$(INTDIR)\dpttrf.obj" : $(SOURCE) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


# End Source File
################################################################################
# Begin Source File

SOURCE=\g3\OTHER\LAPACK\dpotf2.f

"$(INTDIR)\dpotf2.obj" : $(SOURCE) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


# End Source File
################################################################################
# Begin Source File

SOURCE=\g3\OTHER\LAPACK\dpbtrs.f

"$(INTDIR)\dpbtrs.obj" : $(SOURCE) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


# End Source File
################################################################################
# Begin Source File

SOURCE=\g3\OTHER\LAPACK\dpbtrf.f

"$(INTDIR)\dpbtrf.obj" : $(SOURCE) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


# End Source File
################################################################################
# Begin Source File

SOURCE=\g3\OTHER\LAPACK\dpbtf2.f

"$(INTDIR)\dpbtf2.obj" : $(SOURCE) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


# End Source File
################################################################################
# Begin Source File

SOURCE=\g3\OTHER\LAPACK\dpbsv.f

"$(INTDIR)\dpbsv.obj" : $(SOURCE) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


# End Source File
################################################################################
# Begin Source File

SOURCE=\g3\OTHER\LAPACK\dorm2r.f

"$(INTDIR)\dorm2r.obj" : $(SOURCE) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


# End Source File
################################################################################
# Begin Source File

SOURCE=\g3\OTHER\LAPACK\dlasy2.f

"$(INTDIR)\dlasy2.obj" : $(SOURCE) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


# End Source File
################################################################################
# Begin Source File

SOURCE=\g3\OTHER\LAPACK\dlaswp.f

"$(INTDIR)\dlaswp.obj" : $(SOURCE) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


# End Source File
################################################################################
# Begin Source File

SOURCE=\g3\OTHER\LAPACK\dlassq.f

"$(INTDIR)\dlassq.obj" : $(SOURCE) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


# End Source File
################################################################################
# Begin Source File

SOURCE=\g3\OTHER\LAPACK\dlasrt.f

"$(INTDIR)\dlasrt.obj" : $(SOURCE) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


# End Source File
################################################################################
# Begin Source File

SOURCE=\g3\OTHER\LAPACK\dlasr.f

"$(INTDIR)\dlasr.obj" : $(SOURCE) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


# End Source File
################################################################################
# Begin Source File

SOURCE=\g3\OTHER\LAPACK\dlaset.f

"$(INTDIR)\dlaset.obj" : $(SOURCE) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


# End Source File
################################################################################
# Begin Source File

SOURCE=\g3\OTHER\LAPACK\dlascl.f

"$(INTDIR)\dlascl.obj" : $(SOURCE) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


# End Source File
################################################################################
# Begin Source File

SOURCE=\g3\OTHER\LAPACK\dlaruv.f

"$(INTDIR)\dlaruv.obj" : $(SOURCE) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


# End Source File
################################################################################
# Begin Source File

SOURCE=\g3\OTHER\LAPACK\dlartg.f

"$(INTDIR)\dlartg.obj" : $(SOURCE) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


# End Source File
################################################################################
# Begin Source File

SOURCE=\g3\OTHER\LAPACK\dlarnv.f

"$(INTDIR)\dlarnv.obj" : $(SOURCE) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


# End Source File
################################################################################
# Begin Source File

SOURCE=\g3\OTHER\LAPACK\dlarnd.f

"$(INTDIR)\dlarnd.obj" : $(SOURCE) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


# End Source File
################################################################################
# Begin Source File

SOURCE=\g3\OTHER\LAPACK\dlarfx.f

"$(INTDIR)\dlarfx.obj" : $(SOURCE) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


# End Source File
################################################################################
# Begin Source File

SOURCE=\g3\OTHER\LAPACK\dlarfg.f

"$(INTDIR)\dlarfg.obj" : $(SOURCE) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


# End Source File
################################################################################
# Begin Source File

SOURCE=\g3\OTHER\LAPACK\dlarf.f

"$(INTDIR)\dlarf.obj" : $(SOURCE) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


# End Source File
################################################################################
# Begin Source File

SOURCE=\g3\OTHER\LAPACK\dlaran.f

"$(INTDIR)\dlaran.obj" : $(SOURCE) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


# End Source File
################################################################################
# Begin Source File

SOURCE=\g3\OTHER\LAPACK\dlapy3.f

"$(INTDIR)\dlapy3.obj" : $(SOURCE) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


# End Source File
################################################################################
# Begin Source File

SOURCE=\g3\OTHER\LAPACK\dlapy2.f

"$(INTDIR)\dlapy2.obj" : $(SOURCE) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


# End Source File
################################################################################
# Begin Source File

SOURCE=\g3\OTHER\LAPACK\dlaptm.f

"$(INTDIR)\dlaptm.obj" : $(SOURCE) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


# End Source File
################################################################################
# Begin Source File

SOURCE=\g3\OTHER\LAPACK\dlanv2.f

"$(INTDIR)\dlanv2.obj" : $(SOURCE) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


# End Source File
################################################################################
# Begin Source File

SOURCE=\g3\OTHER\LAPACK\dlanst.f

"$(INTDIR)\dlanst.obj" : $(SOURCE) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


# End Source File
################################################################################
# Begin Source File

SOURCE=\g3\OTHER\LAPACK\dlanhs.f

"$(INTDIR)\dlanhs.obj" : $(SOURCE) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


# End Source File
################################################################################
# Begin Source File

SOURCE=\g3\OTHER\LAPACK\dlange.f

"$(INTDIR)\dlange.obj" : $(SOURCE) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


# End Source File
################################################################################
# Begin Source File

SOURCE=\g3\OTHER\LAPACK\dlamch.f

"$(INTDIR)\dlamch.obj" : $(SOURCE) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


# End Source File
################################################################################
# Begin Source File

SOURCE=\g3\OTHER\LAPACK\dlaln2.f

"$(INTDIR)\dlaln2.obj" : $(SOURCE) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


# End Source File
################################################################################
# Begin Source File

SOURCE=\g3\OTHER\LAPACK\dlahqr.f

"$(INTDIR)\dlahqr.obj" : $(SOURCE) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


# End Source File
################################################################################
# Begin Source File

SOURCE=\g3\OTHER\LAPACK\dlagtm.f

"$(INTDIR)\dlagtm.obj" : $(SOURCE) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


# End Source File
################################################################################
# Begin Source File

SOURCE=\g3\OTHER\LAPACK\dlaexc.f

"$(INTDIR)\dlaexc.obj" : $(SOURCE) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


# End Source File
################################################################################
# Begin Source File

SOURCE=\g3\OTHER\LAPACK\dlaev2.f

"$(INTDIR)\dlaev2.obj" : $(SOURCE) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


# End Source File
################################################################################
# Begin Source File

SOURCE=\g3\OTHER\LAPACK\dlae2.f

"$(INTDIR)\dlae2.obj" : $(SOURCE) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


# End Source File
################################################################################
# Begin Source File

SOURCE=\g3\OTHER\LAPACK\dladiv.f

"$(INTDIR)\dladiv.obj" : $(SOURCE) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


# End Source File
################################################################################
# Begin Source File

SOURCE=\g3\OTHER\LAPACK\dlacpy.f

"$(INTDIR)\dlacpy.obj" : $(SOURCE) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


# End Source File
################################################################################
# Begin Source File

SOURCE=\g3\OTHER\LAPACK\dlacon.f

"$(INTDIR)\dlacon.obj" : $(SOURCE) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


# End Source File
################################################################################
# Begin Source File

SOURCE=\g3\OTHER\LAPACK\dlabad.f

"$(INTDIR)\dlabad.obj" : $(SOURCE) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


# End Source File
################################################################################
# Begin Source File

SOURCE=\g3\OTHER\LAPACK\dgttrs.f

"$(INTDIR)\dgttrs.obj" : $(SOURCE) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


# End Source File
################################################################################
# Begin Source File

SOURCE=\g3\OTHER\LAPACK\dgttrf.f

"$(INTDIR)\dgttrf.obj" : $(SOURCE) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


# End Source File
################################################################################
# Begin Source File

SOURCE=\g3\OTHER\LAPACK\dgetrs.f

"$(INTDIR)\dgetrs.obj" : $(SOURCE) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


# End Source File
################################################################################
# Begin Source File

SOURCE=\g3\OTHER\LAPACK\dgetrf.f

"$(INTDIR)\dgetrf.obj" : $(SOURCE) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


# End Source File
################################################################################
# Begin Source File

SOURCE=\g3\OTHER\LAPACK\dgetf2.f

"$(INTDIR)\dgetf2.obj" : $(SOURCE) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


# End Source File
################################################################################
# Begin Source File

SOURCE=\g3\OTHER\LAPACK\dgesv.f

"$(INTDIR)\dgesv.obj" : $(SOURCE) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


# End Source File
################################################################################
# Begin Source File

SOURCE=\g3\OTHER\LAPACK\dgeqr2.f

"$(INTDIR)\dgeqr2.obj" : $(SOURCE) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


# End Source File
################################################################################
# Begin Source File

SOURCE=\g3\OTHER\LAPACK\dgbtrs.f

"$(INTDIR)\dgbtrs.obj" : $(SOURCE) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


# End Source File
################################################################################
# Begin Source File

SOURCE=\g3\OTHER\LAPACK\dgbtrf.f

"$(INTDIR)\dgbtrf.obj" : $(SOURCE) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


# End Source File
################################################################################
# Begin Source File

SOURCE=\g3\OTHER\LAPACK\dgbtf2.f

"$(INTDIR)\dgbtf2.obj" : $(SOURCE) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


# End Source File
################################################################################
# Begin Source File

SOURCE=\g3\OTHER\LAPACK\dgbsv.f

"$(INTDIR)\dgbsv.obj" : $(SOURCE) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


# End Source File
# End Target
################################################################################
# Begin Target

# Name "BLAS - Win32 Release"
# Name "BLAS - Win32 Debug"

!IF  "$(CFG)" == "BLAS - Win32 Release"

!ELSEIF  "$(CFG)" == "BLAS - Win32 Debug"

!ENDIF 

################################################################################
# Begin Source File

SOURCE=\g3\OTHER\BLAS\izamax.f

"$(INTDIR)\izamax.obj" : $(SOURCE) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


# End Source File
################################################################################
# Begin Source File

SOURCE=\g3\OTHER\BLAS\isamax.f

"$(INTDIR)\isamax.obj" : $(SOURCE) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


# End Source File
################################################################################
# Begin Source File

SOURCE=\g3\OTHER\BLAS\idamax.f

"$(INTDIR)\idamax.obj" : $(SOURCE) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


# End Source File
################################################################################
# Begin Source File

SOURCE=\g3\OTHER\BLAS\icamax.f

"$(INTDIR)\icamax.obj" : $(SOURCE) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


# End Source File
################################################################################
# Begin Source File

SOURCE=\g3\OTHER\BLAS\dznrm2.f

"$(INTDIR)\dznrm2.obj" : $(SOURCE) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


# End Source File
################################################################################
# Begin Source File

SOURCE=\g3\OTHER\BLAS\dzasum.f

"$(INTDIR)\dzasum.obj" : $(SOURCE) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


# End Source File
################################################################################
# Begin Source File

SOURCE=\g3\OTHER\BLAS\dtrsv.f

"$(INTDIR)\dtrsv.obj" : $(SOURCE) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


# End Source File
################################################################################
# Begin Source File

SOURCE=\g3\OTHER\BLAS\dtrsm.f

"$(INTDIR)\dtrsm.obj" : $(SOURCE) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


# End Source File
################################################################################
# Begin Source File

SOURCE=\g3\OTHER\BLAS\dtrmm.f

"$(INTDIR)\dtrmm.obj" : $(SOURCE) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


# End Source File
################################################################################
# Begin Source File

SOURCE=\g3\OTHER\BLAS\dtbsv.f

"$(INTDIR)\dtbsv.obj" : $(SOURCE) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


# End Source File
################################################################################
# Begin Source File

SOURCE=\g3\OTHER\BLAS\dsyrk.f

"$(INTDIR)\dsyrk.obj" : $(SOURCE) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


# End Source File
################################################################################
# Begin Source File

SOURCE=\g3\OTHER\BLAS\dsyr2.f

"$(INTDIR)\dsyr2.obj" : $(SOURCE) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


# End Source File
################################################################################
# Begin Source File

SOURCE=\g3\OTHER\BLAS\dsyr.f

"$(INTDIR)\dsyr.obj" : $(SOURCE) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


# End Source File
################################################################################
# Begin Source File

SOURCE=\g3\OTHER\BLAS\dsymv.f

"$(INTDIR)\dsymv.obj" : $(SOURCE) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


# End Source File
################################################################################
# Begin Source File

SOURCE=\g3\OTHER\BLAS\dswap.f

"$(INTDIR)\dswap.obj" : $(SOURCE) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


# End Source File
################################################################################
# Begin Source File

SOURCE=\g3\OTHER\BLAS\dscal.f

"$(INTDIR)\dscal.obj" : $(SOURCE) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


# End Source File
################################################################################
# Begin Source File

SOURCE=\g3\OTHER\BLAS\drotg.f

"$(INTDIR)\drotg.obj" : $(SOURCE) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


# End Source File
################################################################################
# Begin Source File

SOURCE=\g3\OTHER\BLAS\drot.f

"$(INTDIR)\drot.obj" : $(SOURCE) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


# End Source File
################################################################################
# Begin Source File

SOURCE=\g3\OTHER\BLAS\dnrm2.f

"$(INTDIR)\dnrm2.obj" : $(SOURCE) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


# End Source File
################################################################################
# Begin Source File

SOURCE=\g3\OTHER\BLAS\dger.f

"$(INTDIR)\dger.obj" : $(SOURCE) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


# End Source File
################################################################################
# Begin Source File

SOURCE=\g3\OTHER\BLAS\dgemv.f

"$(INTDIR)\dgemv.obj" : $(SOURCE) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


# End Source File
################################################################################
# Begin Source File

SOURCE=\g3\OTHER\BLAS\dgemm.f

"$(INTDIR)\dgemm.obj" : $(SOURCE) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


# End Source File
################################################################################
# Begin Source File

SOURCE=\g3\OTHER\BLAS\dgbmv.f

"$(INTDIR)\dgbmv.obj" : $(SOURCE) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


# End Source File
################################################################################
# Begin Source File

SOURCE=\g3\OTHER\BLAS\ddot.f

"$(INTDIR)\ddot.obj" : $(SOURCE) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


# End Source File
################################################################################
# Begin Source File

SOURCE=\g3\OTHER\BLAS\dcopy.f

"$(INTDIR)\dcopy.obj" : $(SOURCE) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


# End Source File
################################################################################
# Begin Source File

SOURCE=\g3\OTHER\BLAS\daxpy.f

"$(INTDIR)\daxpy.obj" : $(SOURCE) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


# End Source File
################################################################################
# Begin Source File

SOURCE=\g3\OTHER\BLAS\dasum.f

"$(INTDIR)\dasum.obj" : $(SOURCE) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


# End Source File
# End Target
################################################################################
# Begin Target

# Name "UmfPack - Win32 Release"
# Name "UmfPack - Win32 Debug"

!IF  "$(CFG)" == "UmfPack - Win32 Release"

!ELSEIF  "$(CFG)" == "UmfPack - Win32 Debug"

!ENDIF 

################################################################################
# Begin Project Dependency

# Project_Dep_Name "BLAS"

!IF  "$(CFG)" == "UmfPack - Win32 Release"

"BLAS - Win32 Release" : 
   $(MAKE) /$(MAKEFLAGS) /F .\OpenSeesFortran.mak CFG="BLAS - Win32 Release" 

!ELSEIF  "$(CFG)" == "UmfPack - Win32 Debug"

"BLAS - Win32 Debug" : 
   $(MAKE) /$(MAKEFLAGS) /F .\OpenSeesFortran.mak CFG="BLAS - Win32 Debug" 

!ENDIF 

# End Project Dependency
################################################################################
# Begin Source File

SOURCE=\g3\OTHER\UMFPACK\umd2ut.f

"$(INTDIR)\umd2ut.obj" : $(SOURCE) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


# End Source File
################################################################################
# Begin Source File

SOURCE=\g3\OTHER\UMFPACK\umd2su.f

"$(INTDIR)\umd2su.obj" : $(SOURCE) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


# End Source File
################################################################################
# Begin Source File

SOURCE=\g3\OTHER\UMFPACK\umd2so.f

"$(INTDIR)\umd2so.obj" : $(SOURCE) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


# End Source File
################################################################################
# Begin Source File

SOURCE=\g3\OTHER\UMFPACK\umd2sl.f

"$(INTDIR)\umd2sl.obj" : $(SOURCE) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


# End Source File
################################################################################
# Begin Source File

SOURCE=\g3\OTHER\UMFPACK\umd2s2.f

"$(INTDIR)\umd2s2.obj" : $(SOURCE) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


# End Source File
################################################################################
# Begin Source File

SOURCE=\g3\OTHER\UMFPACK\umd2rg.f

"$(INTDIR)\umd2rg.obj" : $(SOURCE) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


# End Source File
################################################################################
# Begin Source File

SOURCE=\g3\OTHER\UMFPACK\umd2rf.f

"$(INTDIR)\umd2rf.obj" : $(SOURCE) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


# End Source File
################################################################################
# Begin Source File

SOURCE=\g3\OTHER\UMFPACK\umd2ra.f

"$(INTDIR)\umd2ra.obj" : $(SOURCE) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


# End Source File
################################################################################
# Begin Source File

SOURCE=\g3\OTHER\UMFPACK\umd2r2.f

"$(INTDIR)\umd2r2.obj" : $(SOURCE) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


# End Source File
################################################################################
# Begin Source File

SOURCE=\g3\OTHER\UMFPACK\umd2r0.f

"$(INTDIR)\umd2r0.obj" : $(SOURCE) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


# End Source File
################################################################################
# Begin Source File

SOURCE=\g3\OTHER\UMFPACK\umd2p2.f

"$(INTDIR)\umd2p2.obj" : $(SOURCE) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


# End Source File
################################################################################
# Begin Source File

SOURCE=\g3\OTHER\UMFPACK\umd2p1.f

"$(INTDIR)\umd2p1.obj" : $(SOURCE) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


# End Source File
################################################################################
# Begin Source File

SOURCE=\g3\OTHER\UMFPACK\umd2of.f

"$(INTDIR)\umd2of.obj" : $(SOURCE) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


# End Source File
################################################################################
# Begin Source File

SOURCE=\g3\OTHER\UMFPACK\umd2lt.f

"$(INTDIR)\umd2lt.obj" : $(SOURCE) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


# End Source File
################################################################################
# Begin Source File

SOURCE=\g3\OTHER\UMFPACK\umd2in.f

"$(INTDIR)\umd2in.obj" : $(SOURCE) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


# End Source File
################################################################################
# Begin Source File

SOURCE=\g3\OTHER\UMFPACK\umd2fg.f

"$(INTDIR)\umd2fg.obj" : $(SOURCE) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


# End Source File
################################################################################
# Begin Source File

SOURCE=\g3\OTHER\UMFPACK\umd2fb.f

"$(INTDIR)\umd2fb.obj" : $(SOURCE) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


# End Source File
################################################################################
# Begin Source File

SOURCE=\g3\OTHER\UMFPACK\umd2fa.f

"$(INTDIR)\umd2fa.obj" : $(SOURCE) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


# End Source File
################################################################################
# Begin Source File

SOURCE=\g3\OTHER\UMFPACK\umd2f2.f

"$(INTDIR)\umd2f2.obj" : $(SOURCE) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


# End Source File
################################################################################
# Begin Source File

SOURCE=\g3\OTHER\UMFPACK\umd2f1.f

"$(INTDIR)\umd2f1.obj" : $(SOURCE) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


# End Source File
################################################################################
# Begin Source File

SOURCE=\g3\OTHER\UMFPACK\umd2f0.f

"$(INTDIR)\umd2f0.obj" : $(SOURCE) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


# End Source File
################################################################################
# Begin Source File

SOURCE=\g3\OTHER\UMFPACK\umd2er.f

"$(INTDIR)\umd2er.obj" : $(SOURCE) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


# End Source File
################################################################################
# Begin Source File

SOURCE=\g3\OTHER\UMFPACK\umd2co.f

"$(INTDIR)\umd2co.obj" : $(SOURCE) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


# End Source File
################################################################################
# Begin Source File

SOURCE=\g3\OTHER\UMFPACK\umd21i.f

"$(INTDIR)\umd21i.obj" : $(SOURCE) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


# End Source File
################################################################################
# Begin Source File

SOURCE=\g3\OTHER\UMFPACK\mc21b.f

"$(INTDIR)\mc21b.obj" : $(SOURCE) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


# End Source File
################################################################################
# Begin Source File

SOURCE=\g3\OTHER\UMFPACK\mc13e.f

"$(INTDIR)\mc13e.obj" : $(SOURCE) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


# End Source File
# End Target
################################################################################
# Begin Target

# Name "feap - Win32 Release"
# Name "feap - Win32 Debug"

!IF  "$(CFG)" == "feap - Win32 Release"

!ELSEIF  "$(CFG)" == "feap - Win32 Debug"

!ENDIF 

################################################################################
# Begin Source File

SOURCE=\g3\SRC\element\feap\getCommon.f

"$(INTDIR)\getCommon.obj" : $(SOURCE) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


# End Source File
################################################################################
# Begin Source File

SOURCE=\g3\SRC\element\feap\fillCommon.f

"$(INTDIR)\fillCommon.obj" : $(SOURCE) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


# End Source File
################################################################################
# Begin Source File

SOURCE=\g3\SRC\element\feap\elmt06.f

"$(INTDIR)\elmt06.obj" : $(SOURCE) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


# End Source File
################################################################################
# Begin Source File

SOURCE=\g3\SRC\element\feap\elmt05.f

"$(INTDIR)\elmt05.obj" : $(SOURCE) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


# End Source File
################################################################################
# Begin Source File

SOURCE=\g3\SRC\element\feap\elmt04.f

"$(INTDIR)\elmt04.obj" : $(SOURCE) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


# End Source File
################################################################################
# Begin Source File

SOURCE=\g3\SRC\element\feap\elmt03.f

"$(INTDIR)\elmt03.obj" : $(SOURCE) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


# End Source File
################################################################################
# Begin Source File

SOURCE=\g3\SRC\element\feap\elmt02.f

"$(INTDIR)\elmt02.obj" : $(SOURCE) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


# End Source File
################################################################################
# Begin Source File

SOURCE=\g3\SRC\element\feap\elmt01.f

"$(INTDIR)\elmt01.obj" : $(SOURCE) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


# End Source File
################################################################################
# Begin Source File

SOURCE=\g3\SRC\element\feap\dummyFeap.f

"$(INTDIR)\dummyFeap.obj" : $(SOURCE) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


# End Source File
# End Target
################################################################################
# Begin Target

# Name "OpenSees - Win32 Release"
# Name "OpenSees - Win32 Debug"

!IF  "$(CFG)" == "OpenSees - Win32 Release"

!ELSEIF  "$(CFG)" == "OpenSees - Win32 Debug"

!ENDIF 

################################################################################
# Begin Source File

SOURCE=\g3\SRC\system_of_eqn\linearSOE\sparseSYM\genmmd.f

"$(INTDIR)\genmmd.obj" : $(SOURCE) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


# End Source File
# End Target
# End Project
################################################################################
