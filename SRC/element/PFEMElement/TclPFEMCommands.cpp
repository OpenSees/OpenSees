/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
** University of California, Berkeley, is strictly prohibited.  See   **
** file 'COPYRIGHT'  in main directory for information on usage and   **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
** Developed by:                                                      **
**   Frank McKenna (fmckenna@ce.berkeley.edu)                         **
**   Gregory L. Fenves (fenves@ce.berkeley.edu)                       **
**   Filip C. Filippou (filippou@ce.berkeley.edu)                     **
**                                                                    **
** ****************************************************************** */
                                                                        
// $Revision: 1.0 $
// $Date: 2012-11-15 20:13:23 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/PFEMElement/TclPFEMCommands.cpp,v $
                                                                        
// Written: Minjie
//
// Description: This file contains the function that is invoked
// by the interpreter when the comand 'PFEM' is invoked by the 
// user.
//
// What: "@(#) Region.h, revA"

#include <PFEMMesher2D.h>
#include <PFEMMesher3D.h>
#include <Domain.h>
#include <Node.h>
#include <Vector.h>
#include <ID.h>
#include <tcl.h>

static PFEMMesher2D theMesher2D;
static PFEMMesher3D theMesher3D;

int
TclModelBuilderPFEM2DCommand(ClientData clientData, Tcl_Interp *interp, int argc,   
                             TCL_Char **argv, Domain* theDomain)
{
    
    if(argc < 2) {
        opserr << "WARNING: at least 2 arguments -- PFEM2D subcommands: discretize, ";
        opserr << "doTriangulation, save \n";
        return TCL_ERROR;
    }

    if(strcmp(argv[1], "discretize") == 0) {
        if(argc < 3) {
            opserr << "WARNING: wrong num of args -- ";
            opserr << "PFEM2D discretize type\n";
            return TCL_ERROR;
        }
        if(strcmp(argv[2], "PSLG") == 0) {
            if(argc < 6) {
                opserr << "WARNING: wrong num of args -- ";
                opserr << "PFEM2D discretize PSLG startnodetag maxarea ndf ";
                opserr << "<-points (x1 y1 ...) -segments (p1 p2 ...) ";
                opserr << "-fix (constrValues) -mass (massValues) -holes (x1 y1 ...)>\n";
                return TCL_ERROR;
            }
            int startnodetag;
            if(Tcl_GetInt(interp, argv[3], &startnodetag) != TCL_OK) {
                opserr<<"WARNING: invalid startnodetag "<<argv[3]<<" -- PFEM2D discretize\n";
                return TCL_ERROR; 
            }
            double maxarea;
            if(Tcl_GetDouble(interp, argv[4], &maxarea) != TCL_OK) {
                opserr<<"WARNING: invalid maxarea "<<argv[4]<<" -- PFEM2D discretize\n";
                return TCL_ERROR; 
            }
            if(maxarea <= 0) {
                opserr<<"WARNING: nonpositive area ";
                opserr<<" -- PFEM2D discretize\n";
                return TCL_ERROR;
            }
            int ndf;
            if(Tcl_GetInt(interp, argv[5], &ndf) != TCL_OK) {
                opserr<<"WARNING: invalid ndf "<<argv[5]<<" -- PFEM2D discretize\n";
                return TCL_ERROR; 
            }
            if(ndf <= 0) {
                opserr<<"WARNING: nonpositive ndf ";
                opserr<<" -- PFEM2D discretize\n";
                return TCL_ERROR;
            }
            
            int loc = 6;
            Vector points(argc-loc), mass(argc-loc), holes(argc-loc);
            ID segments, fix;
            int state=-1, numpoints=0, numsegments=0, nummass=0, numfix=0, numholes=0;
            while(loc < argc) {
                if(strcmp(argv[loc], "-points") == 0) {
                    state = 1;
                } else if(strcmp(argv[loc], "-segments") == 0) {
                    state = 2;
                } else if(strcmp(argv[loc], "-mass") == 0) {
                    state = 3;
                } else if(strcmp(argv[loc], "-fix") == 0) {
                    state = 4;
                } else if(strcmp(argv[loc], "-holes") == 0) {
                    state = 5;
                } else {
                    if(state == 1) {
                        double x;
                        if(Tcl_GetDouble(interp, argv[loc], &x) != TCL_OK) {
                            opserr<<"WARNING: invalid coordinate "<<argv[loc];
                            opserr<<" -- PFEM2D discretize\n";
                            return TCL_ERROR; 
                        }
                        points(numpoints++) = x;
                    } else if(state==2) {
                        int tag;
                        if(Tcl_GetInt(interp, argv[loc], &tag) != TCL_OK) {
                            opserr<<"WARNING: invalid tag "<<argv[loc];
                            opserr<<" -- PFEM2D discretize\n";
                            return TCL_ERROR; 
                        }
                        segments[numsegments++] = tag;
                    } else if(state==3) {
                        double m;
                        if(Tcl_GetDouble(interp, argv[loc], &m) != TCL_OK) {
                            opserr<<"WARNING: invalid mass "<<argv[loc];
                            opserr<<" -- PFEM2D discretize\n";
                            return TCL_ERROR; 
                        }
                        mass(nummass++) = m;
                    } else if(state==4) {
                        int fixity;
                        if(Tcl_GetInt(interp, argv[loc], &fixity) != TCL_OK) {
                            opserr<<"WARNNG: invalid fixity "<<argv[loc];
                            opserr<<" -- PFEM2D discretize\n";
                            return TCL_ERROR;
                        }
                        fix[numfix++] = fixity;
                    } else if(state==5) {
                        double x;
                        if(Tcl_GetDouble(interp, argv[loc], &x) != TCL_OK) {
                            opserr<<"WARNING: invalid coordinate "<<argv[loc];
                            opserr<<" -- PFEM2D discretize\n";
                            return TCL_ERROR; 
                        }
                        holes(numholes++) = x;
                    }
                }
            
                loc++;
            }
            points.resize(numpoints);
            mass.resize(nummass);
            holes.resize(numholes);
            int res = theMesher2D.discretize(startnodetag,points,segments,holes,maxarea,ndf,
                                             fix,mass,theDomain);

            if(res < 0) {
                opserr << "WARNING: failed to discretize\n";
                return TCL_ERROR;
            }
            
            Tcl_SetObjResult(interp, Tcl_NewIntObj(res));

        } else if(strcmp(argv[2], "rectangle")==0 || strcmp(argv[2], "Rectangle")==0) {
            if(argc < 12) {
                opserr<<"WARNING: wrong num of args -- ";
                opserr<<"PFEM2D discretize Rectangle startnodetag x1 y1 hx hy angle nx ny ndf ";
                opserr<<"<-fix (constrValues) -mass (massValues)>\n";
                return TCL_ERROR;
            }
            double x1,y1,hx,hy,angle;
            int nx, ny, ndf, startnodetag;
            if(Tcl_GetInt(interp, argv[3], &startnodetag) != TCL_OK) {
                opserr<<"WARNING: invalid startnodetag "<<argv[3]<<" -- PFEM2D discretize Rectangle\n";
                return TCL_ERROR; 
            }
            if(Tcl_GetDouble(interp, argv[4], &x1) != TCL_OK) {
                opserr<<"WARNING: invalid x1 "<<argv[4]<<" -- PFEM2D discretize Rectangle\n";
                return TCL_ERROR; 
            }
            if(Tcl_GetDouble(interp, argv[5], &y1) != TCL_OK) {
                opserr<<"WARNING: invalid y1 "<<argv[5]<<" -- PFEM2D discretize Rectangle\n";
                return TCL_ERROR; 
            }
            if(Tcl_GetDouble(interp, argv[6], &hx) != TCL_OK) {
                opserr<<"WARNING: invalid hx "<<argv[6]<<" -- PFEM2D discretize Rectangle\n";
                return TCL_ERROR; 
            }
            if(Tcl_GetDouble(interp, argv[7], &hy) != TCL_OK) {
                opserr<<"WARNING: invalid hy "<<argv[7]<<" -- PFEM2D discretize Rectangle\n";
                return TCL_ERROR; 
            }
            if(Tcl_GetDouble(interp, argv[8], &angle) != TCL_OK) {
                opserr<<"WARNING: invalid angle "<<argv[8]<<" -- PFEM2D discretize Rectangle\n";
                return TCL_ERROR; 
            }
            if(Tcl_GetInt(interp, argv[9], &nx) != TCL_OK) {
                opserr<<"WARNING: invalid nx "<<argv[9]<<" -- PFEM2D discretize Rectangle\n";
                return TCL_ERROR; 
            }
            if(Tcl_GetInt(interp, argv[10], &ny) != TCL_OK) {
                opserr<<"WARNING: invalid ny "<<argv[10]<<" -- PFEM2D discretize Rectangle\n";
                return TCL_ERROR; 
            }
            if(Tcl_GetInt(interp, argv[11], &ndf) != TCL_OK) {
                opserr<<"WARNING: invalid ndf "<<argv[11]<<" -- PFEM2D discretize Rectangle\n";
                return TCL_ERROR; 
            }

            int loc = 12;
            Vector mass(argc-loc);
            ID fix;
            int state=-1, numfix=0, nummass=0;
            while(loc < argc) {
                if(strcmp(argv[loc], "-fix") == 0) {
                    state = 1;
                } else if(strcmp(argv[loc], "-mass") == 0) {
                    state = 2;
                } else {
                    if(state == 1) {
                        int fixity;
                        if(Tcl_GetInt(interp, argv[loc], &fixity) != TCL_OK) {
                            opserr<<"WARNING: invalid constraint value "<<argv[loc];
                            opserr<<" -- PFEM2D discretizeRectangle\n";
                            return TCL_ERROR; 
                        }
                        fix[numfix++] = fixity;
                    } else if(state==2) {
                        double m;
                        if(Tcl_GetDouble(interp, argv[loc], &m) != TCL_OK) {
                            opserr<<"WARNING: invalid mass "<<argv[loc];
                            opserr<<" -- PFEM2D discretizeRectangle\n";
                            return TCL_ERROR; 
                        }
                        mass(nummass++) = m;
                    }
                }
            
                loc++;
            }
            mass.resize(nummass);

            int res = theMesher2D.discretize(startnodetag,x1,y1,hx,hy,angle,nx,ny,ndf,fix,mass,theDomain);

            if(res < 0) {
                opserr << "WARNING: failed to discretize\n";
                return TCL_ERROR;
            }
 
            Tcl_SetObjResult(interp, Tcl_NewIntObj(res));

        } else if(strcmp(argv[2], "line")==0 || strcmp(argv[2], "Line")==0) {
            if(argc < 10) {
                opserr<<"WARNING: wrong num of args -- ";
                opserr<<"PFEM2D discretize line startnodetag x1 y1 h angle num ndf ";
                opserr<<"<-fix (constrValues) -mass (massValues)>\n";
                return TCL_ERROR;
            }
            double x1,y1,h,angle;
            int num, ndf, startnodetag;
            if(Tcl_GetInt(interp, argv[3], &startnodetag) != TCL_OK) {
                opserr<<"WARNING: invalid startnodetag "<<argv[3]<<" -- PFEM2D discretize Rectangle\n";
                return TCL_ERROR; 
            }
            if(Tcl_GetDouble(interp, argv[4], &x1) != TCL_OK) {
                opserr<<"WARNING: invalid x1 "<<argv[4]<<" -- PFEM2D discretize Rectangle\n";
                return TCL_ERROR; 
            }
            if(Tcl_GetDouble(interp, argv[5], &y1) != TCL_OK) {
                opserr<<"WARNING: invalid y1 "<<argv[5]<<" -- PFEM2D discretize Rectangle\n";
                return TCL_ERROR; 
            }
            if(Tcl_GetDouble(interp, argv[6], &h) != TCL_OK) {
                opserr<<"WARNING: invalid h "<<argv[6]<<" -- PFEM2D discretize Rectangle\n";
                return TCL_ERROR; 
            }
            if(Tcl_GetDouble(interp, argv[7], &angle) != TCL_OK) {
                opserr<<"WARNING: invalid angle "<<argv[7]<<" -- PFEM2D discretize Rectangle\n";
                return TCL_ERROR; 
            }
            if(Tcl_GetInt(interp, argv[8], &num) != TCL_OK) {
                opserr<<"WARNING: invalid num "<<argv[8]<<" -- PFEM2D discretize Rectangle\n";
                return TCL_ERROR; 
            }
            if(Tcl_GetInt(interp, argv[9], &ndf) != TCL_OK) {
                opserr<<"WARNING: invalid ndf "<<argv[9]<<" -- PFEM2D discretize Rectangle\n";
                return TCL_ERROR; 
            }

            int loc = 10;
            Vector mass(argc-loc);
            ID fix;
            int state=-1, numfix=0, nummass=0;
            while(loc < argc) {
                if(strcmp(argv[loc], "-fix") == 0) {
                    state = 1;
                } else if(strcmp(argv[loc], "-mass") == 0) {
                    state = 2;
                } else {
                    if(state == 1) {
                        int fixity;
                        if(Tcl_GetInt(interp, argv[loc], &fixity) != TCL_OK) {
                            opserr<<"WARNING: invalid constraint value "<<argv[loc];
                            opserr<<" -- PFEM2D discretizeRectangle\n";
                            return TCL_ERROR; 
                        }
                        fix[numfix++] = fixity;
                    } else if(state==2) {
                        double m;
                        if(Tcl_GetDouble(interp, argv[loc], &m) != TCL_OK) {
                            opserr<<"WARNING: invalid mass "<<argv[loc];
                            opserr<<" -- PFEM2D discretizeRectangle\n";
                            return TCL_ERROR; 
                        }
                        mass(nummass++) = m;
                    }
                }
            
                loc++;
            }
            mass.resize(nummass);

            int res = theMesher2D.discretize(startnodetag, x1,y1,h,angle,num,ndf,fix,mass,theDomain);

            if(res < 0) {
                opserr << "WARNING: failed to discretize\n";
                return TCL_ERROR;
            }

            Tcl_SetObjResult(interp, Tcl_NewIntObj(res));
        }

    } else if(strcmp(argv[1], "doTriangulation") == 0) {
        if(argc < 4) {
            opserr << "WARNING: wrong num of args -- ";
            opserr << "PFEM2D doTriangulation alpha -nodes (start end)  ";
            opserr << "<-addnodes (start end) -PFEM2DElement2D starteletag rho mu b1 b2>\n";
            return TCL_ERROR;
        }
        double alpha;
        if(Tcl_GetDouble(interp, argv[2], &alpha) != TCL_OK) {
            opserr<<"WARNING: invalid alpha "<<argv[2]<<" -- PFEM2D doTriangulation\n";
            return TCL_ERROR; 
        }

        int startnodetag=0, endnodetag=-1, startanodetag=0, endanodetag=-1, starteletag=0;
        double rho, mu, b1, b2;
        bool pfem = false;
        int loc = 3;
        while(loc < argc) {
            
            if(strcmp(argv[loc], "-nodes") == 0) {
                if(argc < loc+3) {
                    opserr << "WARNING: wrong num of args -- ";
                    opserr << "PFEM2D doTriangulation alpha -nodes (start end)  ";
                    opserr << "-addnodes (start end) -PFEMElement2D starteletag rho mu b1 b2 \n";
                    return TCL_ERROR;
                }
                if(Tcl_GetInt(interp, argv[loc+1], &startnodetag) != TCL_OK) {
                    opserr<<"WARNING: invalid startnodetag "<<argv[loc+1];
                    opserr<<" -- PFEM2D doTriangulation\n";
                    return TCL_ERROR; 
                }
                if(Tcl_GetInt(interp, argv[loc+2], &endnodetag) != TCL_OK) {
                    opserr<<"WARNING: invalid endnodetag "<<argv[loc+2];
                    opserr<<" -- PFEM2D doTriangulation\n";
                    return TCL_ERROR; 
                }
                loc += 2;
            } else if(strcmp(argv[loc], "-addnodes") == 0) {
                if(argc < loc+3) {
                    opserr << "WARNING: wrong num of args -- ";
                    opserr << "PFEM2D doTriangulation alpha -nodes (start end)  ";
                    opserr << "-addnodes (start end) -PFEMElement2D starteletag rho mu b1 b2 \n";
                    return TCL_ERROR;
                }
                if(Tcl_GetInt(interp, argv[loc+1], &startanodetag) != TCL_OK) {
                    opserr<<"WARNING: invalid startnodetag "<<argv[loc+1];
                    opserr<<" -- PFEM2D doTriangulation\n";
                    return TCL_ERROR; 
                }
                if(Tcl_GetInt(interp, argv[loc+2], &endanodetag) != TCL_OK) {
                    opserr<<"WARNING: invalid endnodetag "<<argv[loc+2];
                    opserr<<" -- PFEM2D doTriangulation\n";
                    return TCL_ERROR; 
                }
                loc += 2;
            } else if(strcmp(argv[loc], "-PFEMElement2D") == 0) {
                pfem = true;
                if(argc < loc+6) {
                    opserr << "WARNING: wrong num of args -- ";
                    opserr << "PFEM2D doTriangulation alpha -nodes (start end)  ";
                    opserr << "-addnodes (start end) -PFEMElement2D starteletag rho mu b1 b2\n";
                    return TCL_ERROR;
                }
                if(Tcl_GetInt(interp, argv[loc+1], &starteletag) != TCL_OK) {
                    opserr<<"WARNING: invalid starteletag "<<argv[loc+1];
                    opserr<<" -- PFEM2D doTriangulation\n";
                    return TCL_ERROR; 
                }
                if(Tcl_GetDouble(interp, argv[loc+2], &rho) != TCL_OK) {
                    opserr<<"WARNING: invalid node rho "<<argv[loc+2];
                    opserr<<" -- PFEM2D doTriangulation\n";
                    return TCL_ERROR; 
                }
                if(Tcl_GetDouble(interp, argv[loc+3], &mu) != TCL_OK) {
                    opserr<<"WARNING: invalid node mu "<<argv[loc+3];
                    opserr<<" -- PFEM2D doTriangulation\n";
                    return TCL_ERROR; 
                }
                if(Tcl_GetDouble(interp, argv[loc+4], &b1) != TCL_OK) {
                    opserr<<"WARNING: invalid node b1 "<<argv[loc+4];
                    opserr<<" -- PFEM2D doTriangulation\n";
                    return TCL_ERROR; 
                }
                if(Tcl_GetDouble(interp, argv[loc+5], &b2) != TCL_OK) {
                    opserr<<"WARNING: invalid node b2 "<<argv[loc+5];
                    opserr<<" -- PFEM2D doTriangulation\n";
                    return TCL_ERROR; 
                }
                loc += 5;
            }
            
            loc++;
        }

        // triangulation
        int res;
        ID eles;
        if(pfem) {
            res = theMesher2D.doTriangulation(starteletag,startnodetag,endnodetag,
                                              alpha,startanodetag,endanodetag,theDomain,
                                              rho,mu,b1,b2);
        } else {
            res = theMesher2D.doTriangulation(startnodetag,endnodetag,
                                              alpha,startanodetag,endanodetag,theDomain,eles);
        }

        if(res < 0) {
            opserr<<"WARNING: failed to do triangulation -- ";
            opserr<<" -- PFEM2D doTriangulation\n";
            return TCL_ERROR; 
        }

        if(pfem) {
            Tcl_SetObjResult(interp, Tcl_NewIntObj(res));
        } else {
            for(int i=0; i<eles.Size(); i++) {
                char buffer[100];
                sprintf(buffer, "%d ", eles(i));
                Tcl_AppendResult(interp, buffer, NULL);
            }
        }

    } else if(strcmp(argv[1], "save") == 0) {

        if(argc < 3) {
            opserr << "WARNING: wrong num of args -- ";
            opserr << "PFEM2D save filename <-structure (start end)> ";
            return TCL_ERROR;            
        }

        int startsnode=0, endsnode=-1;
        if(argc > 5) {
            if(strcmp(argv[3], "-structure") == 0) {
                if(Tcl_GetInt(interp, argv[4], &startsnode) != TCL_OK) {
                    opserr<<"WARNING: invalid startsnode "<<argv[4]<<" -- PFEM2D save\n";
                    return TCL_ERROR; 
                }
                if(Tcl_GetInt(interp, argv[5], &endsnode) != TCL_OK) {
                    opserr<<"WARNING: invalid endsnode "<<argv[5]<<" -- PFEM2D save\n";
                    return TCL_ERROR; 
                }
            }
        }
        theMesher2D.save(argv[2], startsnode, endsnode, theDomain);


    } else if(strcmp(argv[1], "setBoundary") == 0) {
        
        if(argc < 6) {
            opserr << "WARNING: wrong num of args -- ";
            opserr << "PFEM2D setBoundary x1 y1 x2 y2\n";
            return TCL_ERROR;            
        }
        double x1,y1,x2,y2;
        if(Tcl_GetDouble(interp, argv[2], &x1) != TCL_OK) {
            opserr<<"WARNING: invalid node x1 "<<argv[2];
            opserr<<" -- PFEM2D setBoundary\n";
            return TCL_ERROR; 
        }
        if(Tcl_GetDouble(interp, argv[3], &y1) != TCL_OK) {
            opserr<<"WARNING: invalid node y1 "<<argv[3];
            opserr<<" -- PFEM2D setBoundary\n";
            return TCL_ERROR; 
        }
        if(Tcl_GetDouble(interp, argv[4], &x2) != TCL_OK) {
            opserr<<"WARNING: invalid node x2 "<<argv[4];
            opserr<<" -- PFEM2D setBoundary\n";
            return TCL_ERROR; 
        }
        if(Tcl_GetDouble(interp, argv[5], &y2) != TCL_OK) {
            opserr<<"WARNING: invalid node y2 "<<argv[5];
            opserr<<" -- PFEM2D setBoundary\n";
            return TCL_ERROR; 
        }
        theMesher2D.setBoundary(x1,y1,x2,y2);

    } else if(strcmp(argv[1], "removeOutBoundNodes") == 0) {

        if(argc < 2) {
            opserr << "WARNING: wrong num of args -- ";
            opserr << "PFEM2D removeOutBoundNodes <start end>\n";
            return TCL_ERROR;            
        }

        int startnode=0, endnode=-1;
        if(argc > 3) {
            if(Tcl_GetInt(interp, argv[2], &startnode) != TCL_OK) {
                opserr<<"WARNING: invalid startnode "<<argv[2];
                opserr<<" -- PFEM2D removeOutBoundNodes\n";
                return TCL_ERROR; 
            }
            if(Tcl_GetInt(interp, argv[3], &endnode) != TCL_OK) {
                opserr<<"WARNING: invalid endnode "<<argv[3];
                opserr<<" -- PFEM2D removeOutBoundNodes\n";
                return TCL_ERROR; 
            }
            theMesher2D.removeOutBoundNodes(startnode, endnode, theDomain);
        }
    }

    return TCL_OK;
}


int
TclModelBuilderPFEM3DCommand(ClientData clientData, Tcl_Interp *interp, int argc,   
                             TCL_Char **argv, Domain* theDomain)
{
    
    if(argc < 2) {
        opserr << "WARNING: at least 2 arguments -- PFEM3D subcommands: discretize, ";
        opserr << "doTriangulation, save \n";
        return TCL_ERROR;
    }

    if(strcmp(argv[1], "discretize") == 0) {
        // if(argc < 3) {
        //     opserr << "WARNING: wrong num of args -- ";
        //     opserr << "PFEM3D discretize type\n";
        //     return TCL_ERROR;
        // }
        // if(strcmp(argv[2], "PSLG") == 0) {
        //     if(argc < 5) {
        //         opserr << "WARNING: wrong num of args -- ";
        //         opserr << "PFEM3D discretize PSLG maxarea ndf ";
        //         opserr << "<-points (x1 y1 ...) -segments (p1 p2 ...) ";
        //         opserr << "-fix (constrValues) -mass (massValues)>\n";
        //         return TCL_ERROR;
        //     }
        //     double maxarea;
        //     if(Tcl_GetDouble(interp, argv[3], &maxarea) != TCL_OK) {
        //         opserr<<"WARNING: invalid maxarea "<<argv[3]<<" -- PFEM3D discretize\n";
        //         return TCL_ERROR; 
        //     }
        //     if(maxarea <= 0) {
        //         opserr<<"WARNING: nonpositive area ";
        //         opserr<<" -- PFEM3D discretize\n";
        //         return TCL_ERROR;
        //     }
        //     int ndf;
        //     if(Tcl_GetInt(interp, argv[4], &ndf) != TCL_OK) {
        //         opserr<<"WARNING: invalid ndf "<<argv[4]<<" -- PFEM3D discretize\n";
        //         return TCL_ERROR; 
        //     }
        //     if(ndf <= 0) {
        //         opserr<<"WARNING: nonpositive ndf ";
        //         opserr<<" -- PFEM3D discretize\n";
        //         return TCL_ERROR;
        //     }
            
        //     int loc = 5;
        //     Vector points(argc-loc), mass(argc-loc);
        //     ID segments, fix;
        //     int state=-1, numpoints=0, numsegments=0, nummass=0, numfix=0;
        //     while(loc < argc) {
        //         if(strcmp(argv[loc], "-points") == 0) {
        //             state = 1;
        //         } else if(strcmp(argv[loc], "-segments") == 0) {
        //             state = 2;
        //         } else if(strcmp(argv[loc], "-mass") == 0) {
        //             state = 3;
        //         } else if(strcmp(argv[loc], "-fix") == 0) {
        //             state = 4;
        //         } else {
        //             if(state == 1) {
        //                 double x;
        //                 if(Tcl_GetDouble(interp, argv[loc], &x) != TCL_OK) {
        //                     opserr<<"WARNING: invalid coordinate "<<argv[loc];
        //                     opserr<<" -- PFEM3D discretize\n";
        //                     return TCL_ERROR; 
        //                 }
        //                 points(numpoints++) = x;
        //             } else if(state==2) {
        //                 int tag;
        //                 if(Tcl_GetInt(interp, argv[loc], &tag) != TCL_OK) {
        //                     opserr<<"WARNING: invalid tag "<<argv[loc];
        //                     opserr<<" -- PFEM3D discretize\n";
        //                     return TCL_ERROR; 
        //                 }
        //                 segments[numsegments++] = tag;
        //             } else if(state==3) {
        //                 double m;
        //                 if(Tcl_GetDouble(interp, argv[loc], &m) != TCL_OK) {
        //                     opserr<<"WARNING: invalid mass "<<argv[loc];
        //                     opserr<<" -- PFEM3D discretize\n";
        //                     return TCL_ERROR; 
        //                 }
        //                 mass(nummass++) = m;
        //             } else if(state==4) {
        //                 int fixity;
        //                 if(Tcl_GetInt(interp, argv[loc], &fixity) != TCL_OK) {
        //                     opserr<<"WARNNG: invalid fixity "<<argv[loc];
        //                     opserr<<" -- PFEM3D discretize\n";
        //                     return TCL_ERROR;
        //                 }
        //                 fix[numfix++] = fixity;
        //             }
        //         }
            
        //         loc++;
        //     }
        //     points.resize(numpoints);
        //     mass.resize(nummass);
        //     ID nodes;
        //     int res = theMesher3D.discretize(points,segments,maxarea,ndf,fix,mass,theDomain,nodes);
            
        //     if(res < 0) {
        //         opserr << "WARNING: failed to discretize\n";
        //         return TCL_ERROR;
        //     }
            
        //     for(int i=0; i<nodes.Size(); i++) {
        //         char buffer[100];
        //         sprintf(buffer, "%d ", nodes(i));
        //         Tcl_AppendResult(interp, buffer, NULL);
        //     }

        } else if(strcmp(argv[2], "rectangle")==0 || strcmp(argv[2], "Rectangle")==0) {
            // if(argc < 11) {
            //     opserr<<"WARNING: wrong num of args -- ";
            //     opserr<<"PFEM3D discretize Rectangle x1 y1 hx hy angle nx ny ndf ";
            //     opserr<<"<-fix (constrValues) -mass (massValues)>\n";
            //     return TCL_ERROR;
            // }
            // double x1,y1,hx,hy,angle;
            // int nx, ny, ndf;
            // if(Tcl_GetDouble(interp, argv[3], &x1) != TCL_OK) {
            //     opserr<<"WARNING: invalid x1 "<<argv[3]<<" -- PFEM3D discretize Rectangle\n";
            //     return TCL_ERROR; 
            // }
            // if(Tcl_GetDouble(interp, argv[4], &y1) != TCL_OK) {
            //     opserr<<"WARNING: invalid y1 "<<argv[4]<<" -- PFEM3D discretize Rectangle\n";
            //     return TCL_ERROR; 
            // }
            // if(Tcl_GetDouble(interp, argv[5], &hx) != TCL_OK) {
            //     opserr<<"WARNING: invalid hx "<<argv[5]<<" -- PFEM3D discretize Rectangle\n";
            //     return TCL_ERROR; 
            // }
            // if(Tcl_GetDouble(interp, argv[6], &hy) != TCL_OK) {
            //     opserr<<"WARNING: invalid hy "<<argv[6]<<" -- PFEM3D discretize Rectangle\n";
            //     return TCL_ERROR; 
            // }
            // if(Tcl_GetDouble(interp, argv[7], &angle) != TCL_OK) {
            //     opserr<<"WARNING: invalid angle "<<argv[7]<<" -- PFEM3D discretize Rectangle\n";
            //     return TCL_ERROR; 
            // }
            // if(Tcl_GetInt(interp, argv[8], &nx) != TCL_OK) {
            //     opserr<<"WARNING: invalid nx "<<argv[8]<<" -- PFEM3D discretize Rectangle\n";
            //     return TCL_ERROR; 
            // }
            // if(Tcl_GetInt(interp, argv[9], &ny) != TCL_OK) {
            //     opserr<<"WARNING: invalid ny "<<argv[9]<<" -- PFEM3D discretize Rectangle\n";
            //     return TCL_ERROR; 
            // }
            // if(Tcl_GetInt(interp, argv[10], &ndf) != TCL_OK) {
            //     opserr<<"WARNING: invalid ndf "<<argv[10]<<" -- PFEM3D discretize Rectangle\n";
            //     return TCL_ERROR; 
            // }

            // int loc = 11;
            // Vector mass(argc-loc);
            // ID fix;
            // int state=-1, numfix=0, nummass=0;
            // while(loc < argc) {
            //     if(strcmp(argv[loc], "-fix") == 0) {
            //         state = 1;
            //     } else if(strcmp(argv[loc], "-mass") == 0) {
            //         state = 2;
            //     } else {
            //         if(state == 1) {
            //             int fixity;
            //             if(Tcl_GetInt(interp, argv[loc], &fixity) != TCL_OK) {
            //                 opserr<<"WARNING: invalid constraint value "<<argv[loc];
            //                 opserr<<" -- PFEM3D discretizeRectangle\n";
            //                 return TCL_ERROR; 
            //             }
            //             fix[numfix++] = fixity;
            //         } else if(state==2) {
            //             double m;
            //             if(Tcl_GetDouble(interp, argv[loc], &m) != TCL_OK) {
            //                 opserr<<"WARNING: invalid mass "<<argv[loc];
            //                 opserr<<" -- PFEM3D discretizeRectangle\n";
            //                 return TCL_ERROR; 
            //             }
            //             mass(nummass++) = m;
            //         }
            //     }
            
            //     loc++;
            // }
            // mass.resize(nummass);

            // ID nodes(0, (nx+1)*(ny+1));
            // int res = theMesher3D.discretize(x1,y1,hx,hy,angle,nx,ny,ndf,fix,mass,theDomain,nodes);

            // if(res < 0) {
            //     opserr << "WARNING: failed to discretize\n";
            //     return TCL_ERROR;
            // }
 
            // for(int i=0; i<nodes.Size(); i++) {
            //     char buffer[100];
            //     sprintf(buffer, "%d ", nodes(i));
            //     Tcl_AppendResult(interp, buffer, NULL);
            // }


        } else if(strcmp(argv[2], "line")==0 || strcmp(argv[2], "Line")==0) {
            // if(argc < 9) {
        //         opserr<<"WARNING: wrong num of args -- ";
        //         opserr<<"PFEM3D discretize line x1 y1 h angle num ndf ";
        //         opserr<<"<-fix (constrValues) -mass (massValues)>\n";
        //         return TCL_ERROR;
        //     }
        //     double x1,y1,h,angle;
        //     int num, ndf;
        //     if(Tcl_GetDouble(interp, argv[3], &x1) != TCL_OK) {
        //         opserr<<"WARNING: invalid x1 "<<argv[3]<<" -- PFEM3D discretize Rectangle\n";
        //         return TCL_ERROR; 
        //     }
        //     if(Tcl_GetDouble(interp, argv[4], &y1) != TCL_OK) {
        //         opserr<<"WARNING: invalid y1 "<<argv[4]<<" -- PFEM3D discretize Rectangle\n";
        //         return TCL_ERROR; 
        //     }
        //     if(Tcl_GetDouble(interp, argv[5], &h) != TCL_OK) {
        //         opserr<<"WARNING: invalid h "<<argv[5]<<" -- PFEM3D discretize Rectangle\n";
        //         return TCL_ERROR; 
        //     }
        //     if(Tcl_GetDouble(interp, argv[6], &angle) != TCL_OK) {
        //         opserr<<"WARNING: invalid angle "<<argv[6]<<" -- PFEM3D discretize Rectangle\n";
        //         return TCL_ERROR; 
        //     }
        //     if(Tcl_GetInt(interp, argv[7], &num) != TCL_OK) {
        //         opserr<<"WARNING: invalid num "<<argv[7]<<" -- PFEM3D discretize Rectangle\n";
        //         return TCL_ERROR; 
        //     }
        //     if(Tcl_GetInt(interp, argv[8], &ndf) != TCL_OK) {
        //         opserr<<"WARNING: invalid ndf "<<argv[8]<<" -- PFEM3D discretize Rectangle\n";
        //         return TCL_ERROR; 
        //     }

        //     int loc = 9;
        //     Vector mass(argc-loc);
        //     ID fix;
        //     int state=-1, numfix=0, nummass=0;
        //     while(loc < argc) {
        //         if(strcmp(argv[loc], "-fix") == 0) {
        //             state = 1;
        //         } else if(strcmp(argv[loc], "-mass") == 0) {
        //             state = 2;
        //         } else {
        //             if(state == 1) {
        //                 int fixity;
        //                 if(Tcl_GetInt(interp, argv[loc], &fixity) != TCL_OK) {
        //                     opserr<<"WARNING: invalid constraint value "<<argv[loc];
        //                     opserr<<" -- PFEM3D discretizeRectangle\n";
        //                     return TCL_ERROR; 
        //                 }
        //                 fix[numfix++] = fixity;
        //             } else if(state==2) {
        //                 double m;
        //                 if(Tcl_GetDouble(interp, argv[loc], &m) != TCL_OK) {
        //                     opserr<<"WARNING: invalid mass "<<argv[loc];
        //                     opserr<<" -- PFEM3D discretizeRectangle\n";
        //                     return TCL_ERROR; 
        //                 }
        //                 mass(nummass++) = m;
        //             }
        //         }
            
        //         loc++;
        //     }
        //     mass.resize(nummass);

        //     ID nodes(0,num+1);
        //     int res = theMesher3D.discretize(x1,y1,h,angle,num,ndf,fix,mass,theDomain,nodes);

        //     if(res < 0) {
        //         opserr << "WARNING: failed to discretize\n";
        //         return TCL_ERROR;
        //     }
 
        //     for(int i=0; i<nodes.Size(); i++) {
        //         char buffer[100];
        //         sprintf(buffer, "%d ", nodes(i));
        //         Tcl_AppendResult(interp, buffer, NULL);
        //     }

        // }

    } else if(strcmp(argv[1], "doTriangulation") == 0) {
        if(argc < 4) {
            opserr << "WARNING: wrong num of args -- ";
            opserr << "PFEM3D doTriangulation alpha -nodes (start end)  ";
            opserr << "<-addnodes (start end) -PFEMElement3D startele rho mu b1 b2>\n";
            return TCL_ERROR;
        }
        double alpha;
        if(Tcl_GetDouble(interp, argv[2], &alpha) != TCL_OK) {
            opserr<<"WARNING: invalid alpha "<<argv[2]<<" -- PFEM3D doTriangulation\n";
            return TCL_ERROR; 
        }

        int startnode=0, endnode=-1, startanode=0, endanode=-1, startele=0;
        double rho, mu, b1, b2, b3;
        bool pfem = false;
        int loc = 3;
        while(loc < argc) {
            
            if(strcmp(argv[loc], "-nodes") == 0) {
                if(argc < loc+3) {
                    opserr << "WARNING: wrong num of args -- ";
                    opserr << "PFEM3D doTriangulation alpha -nodes (start end) ";
                    opserr << "-addnodes (start end) -PFEMElement3D startele rho mu b1 b2 b3\n";
                    return TCL_ERROR;
                }
                if(Tcl_GetInt(interp, argv[loc+1], &startnode) != TCL_OK) {
                    opserr<<"WARNING: invalid startnode "<<argv[loc+1];
                    opserr<<" -- PFEM3D doTriangulation\n";
                    return TCL_ERROR; 
                }
                if(Tcl_GetInt(interp, argv[loc+2], &endnode) != TCL_OK) {
                    opserr<<"WARNING: invalid endnode "<<argv[loc+2];
                    opserr<<" -- PFEM3D doTriangulation\n";
                    return TCL_ERROR; 
                }
                loc += 2;
            } else if(strcmp(argv[loc], "-addnodes") == 0) {
                if(argc < loc+3) {
                    opserr << "WARNING: wrong num of args -- ";
                    opserr << "PFEM3D doTriangulation alpha -nodes (start end) ";
                    opserr << "-addnodes (start end) -PFEMElement3D starteletag rho mu b1 b2 b3\n";
                    return TCL_ERROR;
                }
                if(Tcl_GetInt(interp, argv[loc+1], &startanode) != TCL_OK) {
                    opserr<<"WARNING: invalid startnodetag "<<argv[loc+1];
                    opserr<<" -- PFEM3D doTriangulation\n";
                    return TCL_ERROR; 
                }
                if(Tcl_GetInt(interp, argv[loc+2], &endanode) != TCL_OK) {
                    opserr<<"WARNING: invalid endnodetag "<<argv[loc+2];
                    opserr<<" -- PFEM3D doTriangulation\n";
                    return TCL_ERROR; 
                }
                loc += 2;
            } else if(strcmp(argv[loc], "-PFEMElement3D") == 0) {
                pfem = true;
                if(argc < loc+7) {
                    opserr << "WARNING: wrong num of args -- ";
                    opserr << "PFEM2D doTriangulation alpha -nodes (start end)  ";
                    opserr << "-addnodes (start end) -PFEMElement3D startele rho mu b1 b2 b3\n";
                    return TCL_ERROR;
                }
                if(Tcl_GetInt(interp, argv[loc+1], &startele) != TCL_OK) {
                    opserr<<"WARNING: invalid startele "<<argv[loc+1];
                    opserr<<" -- PFEM2D doTriangulation\n";
                    return TCL_ERROR; 
                }
                if(Tcl_GetDouble(interp, argv[loc+2], &rho) != TCL_OK) {
                    opserr<<"WARNING: invalid node rho "<<argv[loc+2];
                    opserr<<" -- PFEM2D doTriangulation\n";
                    return TCL_ERROR; 
                }
                if(Tcl_GetDouble(interp, argv[loc+3], &mu) != TCL_OK) {
                    opserr<<"WARNING: invalid node mu "<<argv[loc+3];
                    opserr<<" -- PFEM2D doTriangulation\n";
                    return TCL_ERROR; 
                }
                if(Tcl_GetDouble(interp, argv[loc+4], &b1) != TCL_OK) {
                    opserr<<"WARNING: invalid node b1 "<<argv[loc+4];
                    opserr<<" -- PFEM2D doTriangulation\n";
                    return TCL_ERROR; 
                }
                if(Tcl_GetDouble(interp, argv[loc+5], &b2) != TCL_OK) {
                    opserr<<"WARNING: invalid node b2 "<<argv[loc+5];
                    opserr<<" -- PFEM2D doTriangulation\n";
                    return TCL_ERROR; 
                }
                if(Tcl_GetDouble(interp, argv[loc+6], &b3) != TCL_OK) {
                    opserr<<"WARNING: invalid node b3 "<<argv[loc+6];
                    opserr<<" -- PFEM2D doTriangulation\n";
                    return TCL_ERROR; 
                }
                loc += 6;
            }
            
            loc++;
        }

        // triangulation
        ID eles;
        int res = 0;
        if(pfem) {
            // res = theMesher3D.doTriangulation(nodes,alpha,addNodes,theDomain,
            //                                   eles,rho,mu,b1,b2);
        } else {
            res = theMesher3D.doTriangulation(startnode,endnode,alpha,
                                              startanode,endanode,theDomain,eles);
        }
        if(res < 0) {
            opserr<<"WARNING: failed to do triangulation -- ";
            opserr<<" -- PFEM3D doTriangulation\n";
            return TCL_ERROR; 
        }

        if(pfem) {
            Tcl_SetObjResult(interp, Tcl_NewIntObj(res));
        } else {
            for(int i=0; i<eles.Size(); i++) {
                char buffer[100];
                sprintf(buffer, "%d ", eles(i));
                Tcl_AppendResult(interp, buffer, NULL);
            }
        }

    } else if(strcmp(argv[1], "save") == 0) {

        // if(argc < 3) {
        //     opserr << "WARNING: wrong num of args -- ";
        //     opserr << "PFEM3D save filename <-structure (nd1 nd2 ...)> ";
        //     return TCL_ERROR;            
        // }

        // if(argc > 3) {
        //     ID nodes(0,argc-3);
        //     int numnodes = 0;
        //     if(strcmp(argv[3], "-structure") == 0) {
        //         for(int i=4; i<argc; i++) {
        //             int tag;
        //             if(Tcl_GetInt(interp, argv[i], &tag) != TCL_OK) {
        //                 opserr<<"WARNING: invalid tag "<<argv[i]<<" -- PFEM3D save\n";
        //                 return TCL_ERROR; 
        //             }
        //             nodes[numnodes++] = tag;
        //         }
        //     }
        //     theMesher3D.save(argv[2], nodes, theDomain);
        // } else {
        //     ID nodes;
        //     theMesher3D.save(argv[2], nodes, theDomain);
        // }

    } else if(strcmp(argv[1], "setBoundary") == 0) {
        
        // if(argc < 6) {
        //     opserr << "WARNING: wrong num of args -- ";
        //     opserr << "PFEM3D setBoundary x1 y1 x2 y2\n";
        //     return TCL_ERROR;            
        // }
        // double x1,y1,x2,y2;
        // if(Tcl_GetDouble(interp, argv[2], &x1) != TCL_OK) {
        //     opserr<<"WARNING: invalid node x1 "<<argv[2];
        //     opserr<<" -- PFEM3D setBoundary\n";
        //     return TCL_ERROR; 
        // }
        // if(Tcl_GetDouble(interp, argv[3], &y1) != TCL_OK) {
        //     opserr<<"WARNING: invalid node y1 "<<argv[3];
        //     opserr<<" -- PFEM3D setBoundary\n";
        //     return TCL_ERROR; 
        // }
        // if(Tcl_GetDouble(interp, argv[4], &x2) != TCL_OK) {
        //     opserr<<"WARNING: invalid node x2 "<<argv[4];
        //     opserr<<" -- PFEM3D setBoundary\n";
        //     return TCL_ERROR; 
        // }
        // if(Tcl_GetDouble(interp, argv[5], &y2) != TCL_OK) {
        //     opserr<<"WARNING: invalid node y2 "<<argv[5];
        //     opserr<<" -- PFEM3D setBoundary\n";
        //     return TCL_ERROR; 
        // }
        // theMesher3D.setBoundary(x1,y1,x2,y2);

    } else if(strcmp(argv[1], "removeOutBoundNodes") == 0) {

        // if(argc < 2) {
        //     opserr << "WARNING: wrong num of args -- ";
        //     opserr << "PFEM3D removeOutBoundNodes <nd1 ...>\n";
        //     return TCL_ERROR;            
        // }

        // if(argc > 2) {
        //     ID nodes(argc-2);
        //     for(int i=2; i<argc; i++) {
        //         int tag = 0;
        //         if(Tcl_GetInt(interp, argv[i], &tag) != TCL_OK) {
        //             opserr<<"WARNING: invalid node tag "<<argv[i];
        //             opserr<<" -- PFEM3D removeOutBoundNodes\n";
        //             return TCL_ERROR; 
        //         }
        //         nodes(i-2) = tag;
        //     }

        //     ID nodes2(0, argc-2);
        //     theMesher3D.removeOutBoundNodes(nodes, nodes2, theDomain);

        //     for(int i=0; i<nodes2.Size(); i++) {
        //         char buffer[100];
        //         sprintf(buffer, "%d ", nodes2(i));
        //         Tcl_AppendResult(interp, buffer, NULL);
        //     }
        // }

    }

    return TCL_OK;
}
