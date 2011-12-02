# Microsoft Developer Studio Project File - Name="cblas" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Static Library" 0x0104

CFG=cblas - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "cblas.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "cblas.mak" CFG="cblas - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "cblas - Win32 Release" (based on "Win32 (x86) Static Library")
!MESSAGE "cblas - Win32 Debug" (based on "Win32 (x86) Static Library")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
RSC=rc.exe

!IF  "$(CFG)" == "cblas - Win32 Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "Release"
# PROP BASE Intermediate_Dir "Release"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "..\..\lib"
# PROP Intermediate_Dir "..\..\obj\cblas"
# PROP Target_Dir ""
# ADD BASE CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_MBCS" /D "_LIB" /YX /FD /c
# ADD CPP /nologo /MLd /W3 /GX /O2 /I "..\..\..\cblas" /D "WIN32" /D "NDEBUG" /D "_MBCS" /D "_LIB" /FD /c
# SUBTRACT CPP /YX
# ADD BASE RSC /l 0x409 /d "NDEBUG"
# ADD RSC /l 0x409 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo

!ELSEIF  "$(CFG)" == "cblas - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "Debug"
# PROP BASE Intermediate_Dir "Debug"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "..\..\lib"
# PROP Intermediate_Dir "..\..\obj\cblas"
# PROP Target_Dir ""
# ADD BASE CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_MBCS" /D "_LIB" /YX /FD /GZ /c
# ADD CPP /nologo /W3 /GX /ZI /Od /I "..\..\..\cblas" /D "WIN32" /D "_DEBUG" /D "_MBCS" /D "_LIB" /FR /FD /GZ /c
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

# Name "cblas - Win32 Release"
# Name "cblas - Win32 Debug"
# Begin Source File

SOURCE=..\..\..\OTHER\CBLAS\Cnames.h
# End Source File
# Begin Source File

SOURCE=..\..\..\OTHER\CBLAS\dasum.c
# End Source File
# Begin Source File

SOURCE=..\..\..\OTHER\CBLAS\daxpy.c
# End Source File
# Begin Source File

SOURCE=..\..\..\OTHER\CBLAS\dcabs1.c
# End Source File
# Begin Source File

SOURCE=..\..\..\OTHER\CBLAS\dcopy.c
# End Source File
# Begin Source File

SOURCE=..\..\..\OTHER\CBLAS\ddot.c
# End Source File
# Begin Source File

SOURCE=..\..\..\OTHER\CBLAS\dgemv.c
# End Source File
# Begin Source File

SOURCE=..\..\..\OTHER\CBLAS\dger.c
# End Source File
# Begin Source File

SOURCE=..\..\..\OTHER\CBLAS\dmyblas2.c
# End Source File
# Begin Source File

SOURCE=..\..\..\OTHER\CBLAS\dnrm2.c
# End Source File
# Begin Source File

SOURCE=..\..\..\OTHER\CBLAS\drot.c
# End Source File
# Begin Source File

SOURCE=..\..\..\OTHER\CBLAS\dscal.c
# End Source File
# Begin Source File

SOURCE=..\..\..\OTHER\CBLAS\dsymv.c
# End Source File
# Begin Source File

SOURCE=..\..\..\OTHER\CBLAS\dsyr2.c
# End Source File
# Begin Source File

SOURCE=..\..\..\OTHER\CBLAS\dtrsv.c
# End Source File
# Begin Source File

SOURCE=..\..\..\OTHER\CBLAS\dzasum.c
# End Source File
# Begin Source File

SOURCE=..\..\..\OTHER\CBLAS\dznrm2.c
# End Source File
# Begin Source File

SOURCE=..\..\..\OTHER\CBLAS\f2c.h
# End Source File
# Begin Source File

SOURCE=..\..\..\OTHER\CBLAS\idamax.c
# End Source File
# End Target
# End Project
