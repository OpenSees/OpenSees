# Microsoft Developer Studio Generated NMAKE File, Based on graph.dsp
!IF "$(CFG)" == ""
CFG=graph - Win32 Debug
!MESSAGE No configuration specified. Defaulting to graph - Win32 Debug.
!ENDIF 

!IF "$(CFG)" != "graph - Win32 Release" && "$(CFG)" != "graph - Win32 Debug"
!MESSAGE Invalid configuration "$(CFG)" specified.
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "graph.mak" CFG="graph - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "graph - Win32 Release" (based on "Win32 (x86) Static Library")
!MESSAGE "graph - Win32 Debug" (based on "Win32 (x86) Static Library")
!MESSAGE 
!ERROR An invalid configuration is specified.
!ENDIF 

!IF "$(OS)" == "Windows_NT"
NULL=
!ELSE 
NULL=nul
!ENDIF 

!IF  "$(CFG)" == "graph - Win32 Release"

OUTDIR=.\Release
INTDIR=.\Release
# Begin Custom Macros
OutDir=.\Release
# End Custom Macros

ALL : "$(OUTDIR)\graph.lib"


CLEAN :
	-@erase "$(INTDIR)\ArrayGraph.obj"
	-@erase "$(INTDIR)\ArrayVertexIter.obj"
	-@erase "$(INTDIR)\DOF_Graph.obj"
	-@erase "$(INTDIR)\DOF_GroupGraph.obj"
	-@erase "$(INTDIR)\Graph.obj"
	-@erase "$(INTDIR)\GraphNumberer.obj"
	-@erase "$(INTDIR)\RCM.obj"
	-@erase "$(INTDIR)\SimpleNumberer.obj"
	-@erase "$(INTDIR)\vc60.idb"
	-@erase "$(INTDIR)\Vertex.obj"
	-@erase "$(INTDIR)\VertexIter.obj"
	-@erase "$(OUTDIR)\graph.lib"

"$(OUTDIR)" :
    if not exist "$(OUTDIR)/$(NULL)" mkdir "$(OUTDIR)"

CPP=cl.exe
CPP_PROJ=/nologo /ML /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_MBCS" /D "_LIB" /Fp"$(INTDIR)\graph.pch" /YX /Fo"$(INTDIR)\\" /Fd"$(INTDIR)\\" /FD /c 

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
BSC32_FLAGS=/nologo /o"$(OUTDIR)\graph.bsc" 
BSC32_SBRS= \
	
LIB32=link.exe -lib
LIB32_FLAGS=/nologo /out:"$(OUTDIR)\graph.lib" 
LIB32_OBJS= \
	"$(INTDIR)\ArrayGraph.obj" \
	"$(INTDIR)\ArrayVertexIter.obj" \
	"$(INTDIR)\DOF_Graph.obj" \
	"$(INTDIR)\DOF_GroupGraph.obj" \
	"$(INTDIR)\Graph.obj" \
	"$(INTDIR)\GraphNumberer.obj" \
	"$(INTDIR)\RCM.obj" \
	"$(INTDIR)\SimpleNumberer.obj" \
	"$(INTDIR)\Vertex.obj" \
	"$(INTDIR)\VertexIter.obj"

"$(OUTDIR)\graph.lib" : "$(OUTDIR)" $(DEF_FILE) $(LIB32_OBJS)
    $(LIB32) @<<
  $(LIB32_FLAGS) $(DEF_FLAGS) $(LIB32_OBJS)
<<

!ELSEIF  "$(CFG)" == "graph - Win32 Debug"

OUTDIR=.\..\..\lib
INTDIR=.\..\..\obj\graph
# Begin Custom Macros
OutDir=.\..\..\lib
# End Custom Macros

ALL : "$(OUTDIR)\graph.lib"


CLEAN :
	-@erase "$(INTDIR)\ArrayGraph.obj"
	-@erase "$(INTDIR)\ArrayVertexIter.obj"
	-@erase "$(INTDIR)\DOF_Graph.obj"
	-@erase "$(INTDIR)\DOF_GroupGraph.obj"
	-@erase "$(INTDIR)\Graph.obj"
	-@erase "$(INTDIR)\GraphNumberer.obj"
	-@erase "$(INTDIR)\RCM.obj"
	-@erase "$(INTDIR)\SimpleNumberer.obj"
	-@erase "$(INTDIR)\vc60.idb"
	-@erase "$(INTDIR)\Vertex.obj"
	-@erase "$(INTDIR)\VertexIter.obj"
	-@erase "$(OUTDIR)\graph.lib"

"$(OUTDIR)" :
    if not exist "$(OUTDIR)/$(NULL)" mkdir "$(OUTDIR)"

"$(INTDIR)" :
    if not exist "$(INTDIR)/$(NULL)" mkdir "$(INTDIR)"

CPP=cl.exe
CPP_PROJ=/nologo /MLd /W3 /GX /Od /I "..\..\src\actor\channel" /I "..\..\src\graph\numberer" /I "..\..\src\actor\objectBroker" /I "..\..\src\analysis\dof_grp" /I "..\..\src\analysis\fe_ele" /I "..\..\src\actor\actor" /I "..\..\src\analysis\model\simple" /I "..\..\src\analysis\model" /I "..\..\src\tagged\storage" /I "..\..\src\matrix" /I "..\..\src\graph\graph" /I "..\..\src" /I "..\..\src\tagged" /D "WIN32" /D "_DEBUG" /D "_MBCS" /D "_LIB" /Fp"$(INTDIR)\graph.pch" /YX /Fo"$(INTDIR)\\" /Fd"$(INTDIR)\\" /FD /GZ /c 

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
BSC32_FLAGS=/nologo /o"$(OUTDIR)\graph.bsc" 
BSC32_SBRS= \
	
LIB32=link.exe -lib
LIB32_FLAGS=/nologo /out:"$(OUTDIR)\graph.lib" 
LIB32_OBJS= \
	"$(INTDIR)\ArrayGraph.obj" \
	"$(INTDIR)\ArrayVertexIter.obj" \
	"$(INTDIR)\DOF_Graph.obj" \
	"$(INTDIR)\DOF_GroupGraph.obj" \
	"$(INTDIR)\Graph.obj" \
	"$(INTDIR)\GraphNumberer.obj" \
	"$(INTDIR)\RCM.obj" \
	"$(INTDIR)\SimpleNumberer.obj" \
	"$(INTDIR)\Vertex.obj" \
	"$(INTDIR)\VertexIter.obj"

"$(OUTDIR)\graph.lib" : "$(OUTDIR)" $(DEF_FILE) $(LIB32_OBJS)
    $(LIB32) @<<
  $(LIB32_FLAGS) $(DEF_FLAGS) $(LIB32_OBJS)
<<

!ENDIF 


!IF "$(NO_EXTERNAL_DEPS)" != "1"
!IF EXISTS("graph.dep")
!INCLUDE "graph.dep"
!ELSE 
!MESSAGE Warning: cannot find "graph.dep"
!ENDIF 
!ENDIF 


!IF "$(CFG)" == "graph - Win32 Release" || "$(CFG)" == "graph - Win32 Debug"
SOURCE=..\..\Src\graph\graph\ArrayGraph.cpp

"$(INTDIR)\ArrayGraph.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\Src\graph\graph\ArrayVertexIter.cpp

"$(INTDIR)\ArrayVertexIter.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\Src\graph\graph\DOF_Graph.cpp

"$(INTDIR)\DOF_Graph.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\Src\graph\graph\DOF_GroupGraph.cpp

"$(INTDIR)\DOF_GroupGraph.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\Src\graph\graph\Graph.cpp

"$(INTDIR)\Graph.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\Src\graph\numberer\GraphNumberer.cpp

"$(INTDIR)\GraphNumberer.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\Src\graph\numberer\RCM.cpp

"$(INTDIR)\RCM.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\Src\graph\numberer\SimpleNumberer.cpp

"$(INTDIR)\SimpleNumberer.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\Src\graph\graph\Vertex.cpp

"$(INTDIR)\Vertex.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\Src\graph\graph\VertexIter.cpp

"$(INTDIR)\VertexIter.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)



!ENDIF 

