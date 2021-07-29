
are these files used for anything?

- `SRC/interpreter/tclMain.cpp`
- `tcl/winMain.cpp`

Is there any documentation of the `loadPackage` command? ()

Which is prefered;
- `api/elementAPI_TCL.cpp:theDomain`

Pointers to `class TclModelBuilder`
- `SRC/coordTransformation/TclGeomTransfCommand.cpp:static TclModelBuilder *theTclModelBuilder = 0;`
- `SRC/api/elementAPI_TCL.cpp:static TclModelBuilder* theModelBuilder = 0;`
- `SRC/element/nonlinearBeamColumn/tcl/TclElmtBuilder.cpp:static TclModelBuilder *theTclModelBuilder =0;`
- `SRC/modelbuilder/tcl/TclModelBuilder.cpp:static TclModelBuilder *theTclBuilder =0;`

- `modelbuilder/tcl/TclModelBuilder.cpp:theTclBuilder`
- `api/elementAPI_TCL.cpp:theModelBuilder`
