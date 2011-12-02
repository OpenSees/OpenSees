# Microsoft Developer Studio Project File - Name="superLU" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Static Library" 0x0104

CFG=superLU - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "superLU.mak".
!MESSAGE 
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

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
F90=df.exe
RSC=rc.exe

!IF  "$(CFG)" == "superLU - Win32 Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "Release"
# PROP BASE Intermediate_Dir "Release"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "..\..\lib\release"
# PROP Intermediate_Dir "..\..\obj\superLU\release"
# PROP Target_Dir ""
# ADD BASE CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_MBCS" /D "_LIB" /YX /FD /c
# ADD CPP /nologo /MLd /W3 /GX /O2 /I "..\..\..\other\superLU_3.0\SRC" /I "c:\Program Files\tcl" /I "c:\Program Files\tcl\include" /D "NDEBUG" /D "NO_TIMER" /D "WIN32" /D "_LIB" /D "_MBCS" /D "_TCL84" /FD /c
# SUBTRACT CPP /YX
# ADD BASE RSC /l 0x409 /d "NDEBUG"
# ADD RSC /l 0x409 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo

!ELSEIF  "$(CFG)" == "superLU - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "Debug"
# PROP BASE Intermediate_Dir "Debug"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "..\..\lib\debug"
# PROP Intermediate_Dir "..\..\obj\superLU\debug"
# PROP Target_Dir ""
# ADD BASE CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_MBCS" /D "_LIB" /YX /FD /GZ /c
# ADD CPP /nologo /W3 /GX /ZI /Od /I "..\..\..\other\superLU_3.0\SRC" /I "c:\Program Files\tcl" /I "c:\Program Files\tcl\include" /D "_DEBUG" /D "WINDOWS" /D "WIN32" /D "_LIB" /D "_MBCS" /D "_TCL84" /D "NO_TIMER" /FR /FD /GZ /c
# SUBTRACT CPP /YX
# ADD BASE RSC /l 0x409 /d "_DEBUG"
# ADD RSC /l 0x409 /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo

!ENDIF 

# Begin Target

# Name "superLU - Win32 Release"
# Name "superLU - Win32 Debug"
# Begin Group "Source Files"

# PROP Default_Filter "cpp;c;cxx;rc;def;r;odl;idl;hpj;bat"
# Begin Source File

SOURCE=..\..\..\OTHER\SuperLU_3.0\SRC\colamd.c
# End Source File
# Begin Source File

SOURCE=..\..\..\OTHER\SuperLU_3.0\SRC\dcolumn_bmod.c
# End Source File
# Begin Source File

SOURCE=..\..\..\OTHER\SuperLU_3.0\SRC\dcolumn_dfs.c
# End Source File
# Begin Source File

SOURCE=..\..\..\OTHER\SuperLU_3.0\SRC\dcopy_to_ucol.c
# End Source File
# Begin Source File

SOURCE=..\..\..\OTHER\SuperLU_3.0\SRC\dGetDiagU.c
# End Source File
# Begin Source File

SOURCE=..\..\..\OTHER\SuperLU_3.0\SRC\dgscon.c
# End Source File
# Begin Source File

SOURCE=..\..\..\OTHER\SuperLU_3.0\SRC\dgsequ.c
# End Source File
# Begin Source File

SOURCE=..\..\..\OTHER\SuperLU_3.0\SRC\dgsrfs.c
# End Source File
# Begin Source File

SOURCE=..\..\..\OTHER\SuperLU_3.0\SRC\dgssv.c
# End Source File
# Begin Source File

SOURCE=..\..\..\OTHER\SuperLU_3.0\SRC\dgssvx.c
# End Source File
# Begin Source File

SOURCE=..\..\..\OTHER\SuperLU_3.0\SRC\dgstrf.c
# End Source File
# Begin Source File

SOURCE=..\..\..\OTHER\SuperLU_3.0\SRC\dgstrs.c
# End Source File
# Begin Source File

SOURCE=..\..\..\OTHER\SuperLU_3.0\SRC\dgstrsL.c
# End Source File
# Begin Source File

SOURCE=..\..\..\OTHER\SuperLU_3.0\SRC\dlacon.c
# End Source File
# Begin Source File

SOURCE=..\..\..\OTHER\SuperLU_3.0\SRC\dlamch.c
# End Source File
# Begin Source File

SOURCE=..\..\..\OTHER\SuperLU_3.0\SRC\dlangs.c
# End Source File
# Begin Source File

SOURCE=..\..\..\OTHER\SuperLU_3.0\SRC\dlaqgs.c
# End Source File
# Begin Source File

SOURCE=..\..\..\OTHER\SuperLU_3.0\SRC\dmemory.c
# End Source File
# Begin Source File

SOURCE=..\..\..\OTHER\SuperLU_3.0\SRC\dmyblas2.c
# End Source File
# Begin Source File

SOURCE=..\..\..\OTHER\SuperLU_3.0\SRC\dpanel_bmod.c
# End Source File
# Begin Source File

SOURCE=..\..\..\OTHER\SuperLU_3.0\SRC\dpanel_dfs.c
# End Source File
# Begin Source File

SOURCE=..\..\..\OTHER\SuperLU_3.0\SRC\dpivotgrowth.c
# End Source File
# Begin Source File

SOURCE=..\..\..\OTHER\SuperLU_3.0\SRC\dpivotL.c
# End Source File
# Begin Source File

SOURCE=..\..\..\OTHER\SuperLU_3.0\SRC\dpruneL.c
# End Source File
# Begin Source File

SOURCE=..\..\..\OTHER\SuperLU_3.0\SRC\dreadhb.c
# End Source File
# Begin Source File

SOURCE=..\..\..\OTHER\SuperLU_3.0\SRC\dsnode_bmod.c
# End Source File
# Begin Source File

SOURCE=..\..\..\OTHER\SuperLU_3.0\SRC\dsnode_dfs.c
# End Source File
# Begin Source File

SOURCE=..\..\..\OTHER\SuperLU_3.0\SRC\dsp_blas2.c
# End Source File
# Begin Source File

SOURCE=..\..\..\OTHER\SuperLU_3.0\SRC\dsp_blas3.c
# End Source File
# Begin Source File

SOURCE=..\..\..\OTHER\SuperLU_3.0\SRC\dutil.c
# End Source File
# Begin Source File

SOURCE=..\..\..\OTHER\SuperLU_3.0\SRC\get_perm_c.c
# End Source File
# Begin Source File

SOURCE=..\..\..\OTHER\SuperLU_3.0\SRC\heap_relax_snode.c
# End Source File
# Begin Source File

SOURCE=..\..\..\OTHER\SuperLU_3.0\SRC\lsame.c
# End Source File
# Begin Source File

SOURCE=..\..\..\OTHER\SuperLU_3.0\SRC\memory.c
# End Source File
# Begin Source File

SOURCE=..\..\..\OTHER\SuperLU_3.0\SRC\mmd.c
# End Source File
# Begin Source File

SOURCE=..\..\..\OTHER\SuperLU_3.0\SRC\relax_snode.c
# End Source File
# Begin Source File

SOURCE=..\..\..\OTHER\SuperLU_3.0\SRC\sp_coletree.c
# End Source File
# Begin Source File

SOURCE=..\..\..\OTHER\SuperLU_3.0\SRC\sp_ienv.c
# End Source File
# Begin Source File

SOURCE=..\..\..\OTHER\SuperLU_3.0\SRC\sp_preorder.c
# End Source File
# Begin Source File

SOURCE=..\..\..\OTHER\SuperLU_3.0\SRC\superlu_timer.c
# End Source File
# Begin Source File

SOURCE=..\..\..\OTHER\SuperLU_3.0\SRC\util.c
# End Source File
# Begin Source File

SOURCE=..\..\..\OTHER\SuperLU_3.0\SRC\xerbla.c
# End Source File
# End Group
# Begin Group "Header Files"

# PROP Default_Filter "h;hpp;hxx;hm;inl"
# End Group
# End Target
# End Project
