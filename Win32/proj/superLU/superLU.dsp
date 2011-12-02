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
# ADD CPP /nologo /MLd /W3 /GX /O2 /I "..\..\..\other\superLU" /D "WIN32" /D "NDEBUG" /D "_MBCS" /D "_LIB" /FD /c
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
# ADD CPP /nologo /W3 /GX /ZI /Od /I "..\..\..\other\superLU" /D "WIN32" /D "_DEBUG" /D "_MBCS" /D "_LIB" /D "WINDOWS" /FR /FD /GZ /c
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

SOURCE=..\..\..\OTHER\SuperLU\colamd.c
# End Source File
# Begin Source File

SOURCE=..\..\..\OTHER\SuperLU\dcolumn_bmod.c
# End Source File
# Begin Source File

SOURCE=..\..\..\OTHER\SuperLU\dcolumn_dfs.c
# End Source File
# Begin Source File

SOURCE=..\..\..\OTHER\SuperLU\dcopy_to_ucol.c
# End Source File
# Begin Source File

SOURCE=..\..\..\OTHER\SuperLU\dgscon.c
# End Source File
# Begin Source File

SOURCE=..\..\..\OTHER\SuperLU\dgsequ.c
# End Source File
# Begin Source File

SOURCE=..\..\..\OTHER\SuperLU\dgsrfs.c
# End Source File
# Begin Source File

SOURCE=..\..\..\OTHER\SuperLU\dgssv.c
# End Source File
# Begin Source File

SOURCE=..\..\..\OTHER\SuperLU\dgssvx.c
# End Source File
# Begin Source File

SOURCE=..\..\..\OTHER\SuperLU\dgstrf.c
# End Source File
# Begin Source File

SOURCE=..\..\..\OTHER\SuperLU\dgstrs.c
# End Source File
# Begin Source File

SOURCE=..\..\..\OTHER\SuperLU\dlacon.c
# End Source File
# Begin Source File

SOURCE=..\..\..\OTHER\SuperLU\dlamch.c
# End Source File
# Begin Source File

SOURCE=..\..\..\OTHER\SuperLU\dlangs.c
# End Source File
# Begin Source File

SOURCE=..\..\..\OTHER\SuperLU\dlaqgs.c
# End Source File
# Begin Source File

SOURCE=..\..\..\OTHER\SuperLU\dmemory.c
# End Source File
# Begin Source File

SOURCE=..\..\..\OTHER\SuperLU\dpanel_bmod.c
# End Source File
# Begin Source File

SOURCE=..\..\..\OTHER\SuperLU\dpanel_dfs.c
# End Source File
# Begin Source File

SOURCE=..\..\..\OTHER\SuperLU\dpivotgrowth.c
# End Source File
# Begin Source File

SOURCE=..\..\..\OTHER\SuperLU\dpivotL.c
# End Source File
# Begin Source File

SOURCE=..\..\..\OTHER\SuperLU\dpruneL.c
# End Source File
# Begin Source File

SOURCE=..\..\..\OTHER\SuperLU\dreadhb.c
# End Source File
# Begin Source File

SOURCE=..\..\..\OTHER\SuperLU\dsnode_bmod.c
# End Source File
# Begin Source File

SOURCE=..\..\..\OTHER\SuperLU\dsnode_dfs.c
# End Source File
# Begin Source File

SOURCE=..\..\..\OTHER\SuperLU\dsp_blas2.c
# End Source File
# Begin Source File

SOURCE=..\..\..\OTHER\SuperLU\dsp_blas3.c
# End Source File
# Begin Source File

SOURCE=..\..\..\OTHER\SuperLU\dutil.c
# End Source File
# Begin Source File

SOURCE=..\..\..\OTHER\SuperLU\get_perm_c.c
# End Source File
# Begin Source File

SOURCE=..\..\..\OTHER\SuperLU\lsame.c
# End Source File
# Begin Source File

SOURCE=..\..\..\OTHER\SuperLU\memory.c
# End Source File
# Begin Source File

SOURCE=..\..\..\OTHER\SuperLU\mmd.c
# End Source File
# Begin Source File

SOURCE=..\..\..\OTHER\SuperLU\relax_snode.c
# End Source File
# Begin Source File

SOURCE=..\..\..\OTHER\SuperLU\sp_coletree.c
# End Source File
# Begin Source File

SOURCE=..\..\..\OTHER\SuperLU\sp_ienv.c
# End Source File
# Begin Source File

SOURCE=..\..\..\OTHER\SuperLU\sp_preorder.c
# End Source File
# Begin Source File

SOURCE=..\..\..\OTHER\SuperLU\superlu_timer.c
# End Source File
# Begin Source File

SOURCE=..\..\..\OTHER\SuperLU\util.c
# End Source File
# Begin Source File

SOURCE=..\..\..\OTHER\SuperLU\xerbla.c
# End Source File
# End Group
# Begin Group "Header Files"

# PROP Default_Filter "h;hpp;hxx;hm;inl"
# Begin Source File

SOURCE=..\..\..\OTHER\SuperLU\Cnames.h
# End Source File
# Begin Source File

SOURCE=..\..\..\OTHER\SuperLU\colamd.h
# End Source File
# Begin Source File

SOURCE=..\..\..\OTHER\SuperLU\dsp_defs.h
# End Source File
# Begin Source File

SOURCE=..\..\..\OTHER\SuperLU\supermatrix.h
# End Source File
# Begin Source File

SOURCE=..\..\..\OTHER\SuperLU\util.h
# End Source File
# End Group
# End Target
# End Project
