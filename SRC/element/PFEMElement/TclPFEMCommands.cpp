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
                                                                        
// $Revision: 1.1 $
// $Date: 2012-11-15 20:13:23 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/PFEMElement/TclPFEMCommands.cpp,v $
                                                                        
// Written: Minjie
//
// Description: This file contains the function that is invoked
// by the interpreter when the command 'PFEM2D' is invoked by the 
// user.
//

#include <PFEMMesher2D.h>
#include <PFEMMesher3D.h>
#include <Domain.h>
#include <Node.h>
#include <NodeIter.h>
#include <Vector.h>
#include <ID.h>
#include <tcl.h>
#include <string>
#include <PFEMElement2DBubble.h>
#include <MeshRegion.h>

static PFEMMesher2D theMesher2D;
static PFEMMesher3D theMesher3D;

int
TclModelBuilderPFEM2DCommand(ClientData clientData, Tcl_Interp *interp, int argc,   
                             TCL_Char **argv, Domain* theDomain)
{
    
    if(argc < 2) {
        opserr << "WARNING: at least 2 arguments -- PFEM2D subcommands\n";
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
                opserr << "<-points {x1 y1 ...} -segments {p1 p2 ...} -vel {vx vy} ";
                opserr << "-fix {constrValues} -mass {massValues} -holes {x1 y1 ...} >\n";
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
            Vector points, mass, holes, vel, segments, fix;
            while(loc < argc) {

                Vector* vecPtr = 0;
                if(strcmp(argv[loc], "-points") == 0) {
                    vecPtr = &points;
                    loc++;

                } else if(strcmp(argv[loc], "-segments") == 0) {
                    vecPtr = &segments;
                    loc++;

                } else if(strcmp(argv[loc], "-mass") == 0) {
                    vecPtr = &mass;
                    loc++;

                } else if(strcmp(argv[loc], "-fix") == 0) {
                    vecPtr = &fix;
                    loc++;
  
                } else if(strcmp(argv[loc], "-holes") == 0) {
                    vecPtr = &holes;
                    loc++;

                } else if(strcmp(argv[loc], "-vel") == 0) {
                    vecPtr = &vel;
                    loc++;
                }

                if(vecPtr != 0) {
                    int num = 0;
                    const char** argvPtr = 0;
                    if(Tcl_SplitList(interp, argv[loc], &num, &argvPtr)!=TCL_OK) {
                        opserr<<"WARNING: failed to read the list "<<argv[loc];
                        opserr<<" -- PFEM2D discretize\n";
                        return TCL_ERROR;
                    }
                    if(num > 0) vecPtr->resize(num);
                    for(int i=0; i<num; i++) {
                        double x;
                        if(Tcl_GetDouble(interp, argvPtr[i], &x) != TCL_OK) {
                            opserr<<"WARNING: invalid input "<<argvPtr[i];
                            opserr<<" -- PFEM2D discretize\n";
                            return TCL_ERROR; 
                        }
                        (*vecPtr)(i) = x;
                    }

                    if(argvPtr != 0) Tcl_Free((char *) argvPtr);
                }

                loc++;
            }

            int endnodetag = startnodetag;
            int res = theMesher2D.discretize(startnodetag,points,segments,holes,maxarea,ndf,
                                             fix,vel,mass,theDomain, endnodetag);

            if(res < 0) {
                opserr << "WARNING: failed to discretize\n";
                return TCL_ERROR;
            }
            
            Tcl_SetObjResult(interp, Tcl_NewIntObj(endnodetag));

        } else if(strcmp(argv[2], "particles") == 0 || strcmp(argv[2], "Particles") == 0) {
            if(argc < 5) {
                opserr << "WARNING: wrong num of args -- ";
                opserr << "PFEM2D discretize particles startnodetag ndf ";
                opserr << "<-points {x1 y1 ...} -vel {vx vy} ";
                opserr << "-fix {constrValues} -mass {massValues} >\n";
                return TCL_ERROR;
            }
            int startnodetag;
            if(Tcl_GetInt(interp, argv[3], &startnodetag) != TCL_OK) {
                opserr<<"WARNING: invalid startnodetag "<<argv[3]<<" -- PFEM2D discretize\n";
                return TCL_ERROR; 
            }
            int ndf;
            if(Tcl_GetInt(interp, argv[4], &ndf) != TCL_OK) {
                opserr<<"WARNING: invalid ndf "<<argv[4]<<" -- PFEM2D discretize\n";
                return TCL_ERROR; 
            }
            if(ndf <= 0) {
                opserr<<"WARNING: nonpositive ndf ";
                opserr<<" -- PFEM2D discretize\n";
                return TCL_ERROR;
            }

            int loc = 5;
            Vector points, mass, vel, fix;
            while(loc < argc) {

                Vector* vecPtr = 0;
                if(strcmp(argv[loc], "-points") == 0) {
                    vecPtr = &points;
                    loc++;

                } else if(strcmp(argv[loc], "-mass") == 0) {
                    vecPtr = &mass;
                    loc++;

                } else if(strcmp(argv[loc], "-fix") == 0) {
                    vecPtr = &fix;
                    loc++;

                } else if(strcmp(argv[loc], "-vel") == 0) {
                    vecPtr = &vel;
                    loc++;
                }

                if(vecPtr != 0) {
                    int num = 0;
                    const char** argvPtr = 0;
                    if(Tcl_SplitList(interp, argv[loc], &num, &argvPtr)!=TCL_OK) {
                        opserr<<"WARNING: failed to read the list "<<argv[loc];
                        opserr<<" -- PFEM2D discretize\n";
                        return TCL_ERROR;
                    }
                    if(num > 0) vecPtr->resize(num);
                    for(int i=0; i<num; i++) {
                        double x;
                        if(Tcl_GetDouble(interp, argvPtr[i], &x) != TCL_OK) {
                            opserr<<"WARNING: invalid input "<<argvPtr[i];
                            opserr<<" -- PFEM2D discretize\n";
                            return TCL_ERROR; 
                        }
                        (*vecPtr)(i) = x;
                    }

                    if(argvPtr != 0) Tcl_Free((char *) argvPtr);
                }

                loc++;
            }

            int endnodetag = startnodetag;
            int res = theMesher2D.discretize(startnodetag,points,ndf,
                                             fix,vel,mass,theDomain, 
                                             endnodetag);

            if(res < 0) {
                opserr << "WARNING: failed to discretize\n";
                return TCL_ERROR;
            }
            
            Tcl_SetObjResult(interp, Tcl_NewIntObj(endnodetag));

        } else if(strcmp(argv[2], "rectangle")==0 || strcmp(argv[2], "Rectangle")==0) {
            if(argc < 12) {
                opserr<<"WARNING: wrong num of args -- ";
                opserr<<"PFEM2D discretize Rectangle startnodetag x1 y1 hx hy angle nx ny ndf ";
                opserr<<"<-fix {constrValues} -mass {massValues} -vel {vx vy}";
                opserr<<" -boundary {boundaryValues}>\n";
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
            Vector mass, vel, fix, boundary(4);
            boundary += 1;
            while(loc < argc) {

                Vector* vecPtr = 0;

                if(strcmp(argv[loc], "-mass") == 0) {
                    vecPtr = &mass;
                    loc++;

                } else if(strcmp(argv[loc], "-fix") == 0) {
                    vecPtr = &fix;
                    loc++;
  
                } else if(strcmp(argv[loc], "-vel") == 0) {
                    vecPtr = &vel;
                    loc++;
                } else if(strcmp(argv[loc], "-boundary") == 0) {
                    vecPtr = &boundary;
                    loc++;
                }

                if(vecPtr != 0) {
                    int num = 0;
                    const char** argvPtr = 0;
                    if(Tcl_SplitList(interp, argv[loc], &num, &argvPtr)!=TCL_OK) {
                        opserr<<"WARNING: failed to read the list "<<argv[loc];
                        opserr<<" -- PFEM2D discretize\n";
                        return TCL_ERROR;
                    }
                    if(num > 0) vecPtr->resize(num);
                    for(int i=0; i<num; i++) {
                        double x;
                        if(Tcl_GetDouble(interp, argvPtr[i], &x) != TCL_OK) {
                            opserr<<"WARNING: invalid input "<<argvPtr[i];
                            opserr<<" -- PFEM2D discretize\n";
                            return TCL_ERROR; 
                        }
                        (*vecPtr)(i) = x;
                    }

                    if(argvPtr != 0) Tcl_Free((char *) argvPtr);
                }

                loc++;
            }

            int endnodetag = startnodetag;
            int res = theMesher2D.discretize(startnodetag,x1,y1,hx,hy,angle,nx,ny,ndf,fix,vel,mass,
                                             boundary,theDomain,endnodetag);

            if(res < 0) {
                opserr << "WARNING: failed to discretize\n";
                return TCL_ERROR;
            }
 
            Tcl_SetObjResult(interp, Tcl_NewIntObj(endnodetag));

        } else if(strcmp(argv[2], "line")==0 || strcmp(argv[2], "Line")==0) {
            if(argc < 10) {
                opserr<<"WARNING: wrong num of args -- ";
                opserr<<"PFEM2D discretize line startnodetag x1 y1 h angle num ndf ";
                opserr<<"-fix {constrValues} -mass {massValues} -vel {vx vy}";
                opserr<<"-boundary {boundaryValues}\n";
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
            Vector mass, vel, fix, boundary(2);
            boundary += 1;
            while(loc < argc) {

                Vector* vecPtr = 0;

                if(strcmp(argv[loc], "-mass") == 0) {
                    vecPtr = &mass;
                    loc++;

                } else if(strcmp(argv[loc], "-fix") == 0) {
                    vecPtr = &fix;
                    loc++;
  
                } else if(strcmp(argv[loc], "-vel") == 0) {
                    vecPtr = &vel;
                    loc++;
                } else if(strcmp(argv[loc], "-boundary") == 0) {
                    vecPtr = &boundary;
                    loc++;
                }

                if(vecPtr != 0) {
                    int num = 0;
                    const char** argvPtr = 0;
                    if(Tcl_SplitList(interp, argv[loc], &num, &argvPtr)!=TCL_OK) {
                        opserr<<"WARNING: failed to read the list "<<argv[loc];
                        opserr<<" -- PFEM2D discretize\n";
                        return TCL_ERROR;
                    }
                    if(num > 0) vecPtr->resize(num);
                    for(int i=0; i<num; i++) {
                        double x;
                        if(Tcl_GetDouble(interp, argvPtr[i], &x) != TCL_OK) {
                            opserr<<"WARNING: invalid input "<<argvPtr[i];
                            opserr<<" -- PFEM2D discretize\n";
                            return TCL_ERROR; 
                        }
                        (*vecPtr)(i) = x;
                    }

                    if(argvPtr != 0) Tcl_Free((char *) argvPtr);
                }

                loc++;
            }

            int endnodetag = startnodetag;
            int res = theMesher2D.discretize(startnodetag, x1,y1,h,angle,num,ndf,fix,vel,mass,
                                             boundary,theDomain,endnodetag);

            if(res < 0) {
                opserr << "WARNING: failed to discretize\n";
                return TCL_ERROR;
            }

            Tcl_SetObjResult(interp, Tcl_NewIntObj(endnodetag));

        } else if(strcmp(argv[2], "triangle")==0 || strcmp(argv[2], "Triangle")==0) {

            if(argc < 13) {
                opserr<<"WARNING: wrong num of args -- ";
                opserr<<"PFEM2D discretize triangle startnodetag x1 y1 x2 y2 x3 y3 n1 n2 ndf ";
                opserr<<"<-fix {constrValues} -mass {massValues} -vel {vx vy} ";
                opserr<<"-boundary {boundaryValues}>\n";
                return TCL_ERROR;
            }
            double x1,y1,x2,y2,x3,y3;
            int n1,n2,ndf,startnode;
            if(Tcl_GetInt(interp, argv[3], &startnode) != TCL_OK) {
                opserr<<"WARNING: invalid startnod "<<argv[3]<<" -- PFEM2D discretize triangle\n";
                return TCL_ERROR; 
            }
            if(Tcl_GetDouble(interp, argv[4], &x1) != TCL_OK) {
                opserr<<"WARNING: invalid x1 "<<argv[4]<<" -- PFEM2D discretize triangle\n";
                return TCL_ERROR; 
            }
            if(Tcl_GetDouble(interp, argv[5], &y1) != TCL_OK) {
                opserr<<"WARNING: invalid y1 "<<argv[5]<<" -- PFEM2D discretize triangle\n";
                return TCL_ERROR; 
            }
            if(Tcl_GetDouble(interp, argv[6], &x2) != TCL_OK) {
                opserr<<"WARNING: invalid x2 "<<argv[6]<<" -- PFEM2D discretize triangle\n";
                return TCL_ERROR; 
            }
            if(Tcl_GetDouble(interp, argv[7], &y2) != TCL_OK) {
                opserr<<"WARNING: invalid y2 "<<argv[7]<<" -- PFEM2D discretize triangle\n";
                return TCL_ERROR; 
            }
            if(Tcl_GetDouble(interp, argv[8], &x3) != TCL_OK) {
                opserr<<"WARNING: invalid x3 "<<argv[8]<<" -- PFEM2D discretize triangle\n";
                return TCL_ERROR; 
            }
            if(Tcl_GetDouble(interp, argv[9], &y3) != TCL_OK) {
                opserr<<"WARNING: invalid y3 "<<argv[9]<<" -- PFEM2D discretize triangle\n";
                return TCL_ERROR; 
            }
            if(Tcl_GetInt(interp, argv[10], &n1) != TCL_OK) {
                opserr<<"WARNING: invalid n1 "<<argv[10]<<" -- PFEM2D discretize triangle\n";
                return TCL_ERROR; 
            }
            if(Tcl_GetInt(interp, argv[11], &n2) != TCL_OK) {
                opserr<<"WARNING: invalid n2 "<<argv[11]<<" -- PFEM2D discretize triangle\n";
                return TCL_ERROR; 
            }
            if(Tcl_GetInt(interp, argv[12], &ndf) != TCL_OK) {
                opserr<<"WARNING: invalid ndf "<<argv[12]<<" -- PFEM2D discretize triangle\n";
                return TCL_ERROR; 
            }

            int loc = 13;
            Vector mass, vel, fix, boundary(3);
            boundary += 1;
            while(loc < argc) {

                Vector* vecPtr = 0;

                if(strcmp(argv[loc], "-mass") == 0) {
                    vecPtr = &mass;
                    loc++;

                } else if(strcmp(argv[loc], "-fix") == 0) {
                    vecPtr = &fix;
                    loc++;
  
                } else if(strcmp(argv[loc], "-vel") == 0) {
                    vecPtr = &vel;
                    loc++;
                } else if(strcmp(argv[loc], "-boundary") == 0) {
                    vecPtr = &boundary;
                    loc++;
                }

                if(vecPtr != 0) {
                    int num = 0;
                    const char** argvPtr = 0;
                    if(Tcl_SplitList(interp, argv[loc], &num, &argvPtr)!=TCL_OK) {
                        opserr<<"WARNING: failed to read the list "<<argv[loc];
                        opserr<<" -- PFEM2D discretize\n";
                        return TCL_ERROR;
                    }
                    if(num > 0) vecPtr->resize(num);
                    for(int i=0; i<num; i++) {
                        double x;
                        if(Tcl_GetDouble(interp, argvPtr[i], &x) != TCL_OK) {
                            opserr<<"WARNING: invalid input "<<argvPtr[i];
                            opserr<<" -- PFEM2D discretize\n";
                            return TCL_ERROR; 
                        }
                        (*vecPtr)(i) = x;
                    }

                    if(argvPtr != 0) Tcl_Free((char *) argvPtr);
                }

                loc++;
            }

            int endnode = startnode;

            int res = theMesher2D.discretize(startnode,x1,y1,x2,y2,x3,y3,n1,n2,ndf,fix,vel,mass,
                                             boundary,theDomain,endnode);
            if(res < 0) {
                opserr << "WARNING: failed to discretize\n";
                return TCL_ERROR;
            }

            Tcl_SetObjResult(interp, Tcl_NewIntObj(endnode));

        } else if(strcmp(argv[2], "circle")==0 || strcmp(argv[2], "Circle")==0) {

            if(argc < 11) {
                opserr<<"WARNING: wrong num of args -- ";
                opserr<<"PFEM2D discretize circle startnode xc yc r1 r2 nc nr ndf ";
                opserr<<"<-fix {constrValues} -mass {massValues} -vel {vx vy} ";
                opserr<<"-angleRange {start end} "; 
                opserr<<"-boundary {boundaryValues}>\n";
                return TCL_ERROR;
            }

            double xc, yc, r1, r2;
            int nc,nr,ndf,startnode;
            if(Tcl_GetInt(interp, argv[3], &startnode) != TCL_OK) {
                opserr<<"WARNING: invalid startnod "<<argv[3]<<" -- PFEM2D discretize circle\n";
                return TCL_ERROR; 
            }
            if(Tcl_GetDouble(interp, argv[4], &xc) != TCL_OK) {
                opserr<<"WARNING: invalid xc "<<argv[4]<<" -- PFEM2D discretize circle\n";
                return TCL_ERROR; 
            }
            if(Tcl_GetDouble(interp, argv[5], &yc) != TCL_OK) {
                opserr<<"WARNING: invalid yc "<<argv[5]<<" -- PFEM2D discretize circle\n";
                return TCL_ERROR; 
            }
            if(Tcl_GetDouble(interp, argv[6], &r1) != TCL_OK) {
                opserr<<"WARNING: invalid r1 "<<argv[6]<<" -- PFEM2D discretize circle\n";
                return TCL_ERROR; 
            }
            if(Tcl_GetDouble(interp, argv[7], &r2) != TCL_OK) {
                opserr<<"WARNING: invalid r2 "<<argv[7]<<" -- PFEM2D discretize circle\n";
                return TCL_ERROR; 
            }
            if(Tcl_GetInt(interp, argv[8], &nc) != TCL_OK) {
                opserr<<"WARNING: invalid nc "<<argv[8]<<" -- PFEM2D discretize circle\n";
                return TCL_ERROR; 
            }
            if(Tcl_GetInt(interp, argv[9], &nr) != TCL_OK) {
                opserr<<"WARNING: invalid nr "<<argv[9]<<" -- PFEM2D discretize circle\n";
                return TCL_ERROR; 
            }
            if(Tcl_GetInt(interp, argv[10], &ndf) != TCL_OK) {
                opserr<<"WARNING: invalid ndf "<<argv[10]<<" -- PFEM2D discretize circle\n";
                return TCL_ERROR; 
            }
            int loc = 11;
            Vector mass, vel, fix, boundary(2), angle(2);
            boundary += 1;
            while(loc < argc) {

                Vector* vecPtr = 0;

                if(strcmp(argv[loc], "-mass") == 0) {
                    vecPtr = &mass;
                    loc++;

                } else if(strcmp(argv[loc], "-fix") == 0) {
                    vecPtr = &fix;
                    loc++;
  
                } else if(strcmp(argv[loc], "-vel") == 0) {
                    vecPtr = &vel;
                    loc++;
                } else if(strcmp(argv[loc], "-boundary") == 0) {
                    vecPtr = &boundary;
                    loc++;
                } else if(strcmp(argv[loc], "-angleRange") == 0) {
                    vecPtr = &angle;
                    loc++;
                }

                if(vecPtr != 0) {
                    int num = 0;
                    const char** argvPtr = 0;
                    if(Tcl_SplitList(interp, argv[loc], &num, &argvPtr)!=TCL_OK) {
                        opserr<<"WARNING: failed to read the list "<<argv[loc];
                        opserr<<" -- PFEM2D discretize\n";
                        return TCL_ERROR;
                    }
                    if(num > 0) vecPtr->resize(num);
                    for(int i=0; i<num; i++) {
                        double x;
                        if(Tcl_GetDouble(interp, argvPtr[i], &x) != TCL_OK) {
                            opserr<<"WARNING: invalid input "<<argvPtr[i];
                            opserr<<" -- PFEM2D discretize\n";
                            return TCL_ERROR; 
                        }
                        (*vecPtr)(i) = x;
                    }

                    if(argvPtr != 0) Tcl_Free((char *) argvPtr);
                }

                loc++;
            }

            int endnode = startnode;

            int res = theMesher2D.discretize(startnode,xc,yc,r1,r2,nc,nr,ndf,fix,vel,mass,
                                             boundary,angle,theDomain,endnode);
            if(res < 0) {
                opserr << "WARNING: failed to discretize\n";
                return TCL_ERROR;
            }

            Tcl_SetObjResult(interp, Tcl_NewIntObj(endnode));


        } else if(strcmp(argv[2], "frame")==0 || strcmp(argv[2], "Frame")==0) {

            if(argc < 11) {
                opserr<<"WARNING: wrong num of args -- ";
                opserr<<"PFEM2D discretize frame startnode x y hcol hbeam ndf ";
                opserr<<"{spanlist} {hlist} -fix {constrValues} ";
                opserr<<"-colmass {masslist} -beammass {masslist} ";
                opserr<<"-roofmass {masslist}";
                opserr<<"-colvel {vellist} -beamvel {vellist}>\n";
                return TCL_ERROR;
            }
            int startnode, ndf;
            double x,y,hcol, hbeam;
            if(Tcl_GetInt(interp, argv[3], &startnode) != TCL_OK) {
                opserr<<"WARNING: invalid startnod "<<argv[3]<<" -- PFEM2D discretize frame\n";
                return TCL_ERROR; 
            }
            if(Tcl_GetDouble(interp, argv[4], &x) != TCL_OK) {
                opserr<<"WARNING: invalid x "<<argv[4]<<" -- PFEM2D discretize frame\n";
                return TCL_ERROR; 
            }
            if(Tcl_GetDouble(interp, argv[5], &y) != TCL_OK) {
                opserr<<"WARNING: invalid y "<<argv[5]<<" -- PFEM2D discretize frame\n";
                return TCL_ERROR; 
            }
            if(Tcl_GetDouble(interp, argv[6], &hcol) != TCL_OK) {
                opserr<<"WARNING: invalid ncol "<<argv[6]<<" -- PFEM2D discretize frame\n";
                return TCL_ERROR; 
            }
            if(Tcl_GetDouble(interp, argv[7], &hbeam) != TCL_OK) {
                opserr<<"WARNING: invalid hbeam "<<argv[7]<<" -- PFEM2D discretize frame\n";
                return TCL_ERROR; 
            }
            if(Tcl_GetInt(interp, argv[8], &ndf) != TCL_OK) {
                opserr<<"WARNING: invalid nbeam "<<argv[8]<<" -- PFEM2D discretize frame\n";
                return TCL_ERROR; 
            }

            int loc = 9;
            Vector colmass, beammass, colvel, beamvel, roofmass;
            Vector fix, spans, heights;
            while(loc < argc) {

                Vector* vecPtr = 0;

                if(loc == 9) {
                    vecPtr = &spans;

                } else if(loc == 10) {
                    vecPtr = &heights;

                } else if(strcmp(argv[loc], "-colmass") == 0) {
                    vecPtr = &colmass;
                    loc++;

                } else if(strcmp(argv[loc], "-beammass") == 0) {
                    vecPtr = &beammass;
                    loc++;

                } else if(strcmp(argv[loc], "-fix") == 0) {
                    vecPtr = &fix;
                    loc++;
  
                } else if(strcmp(argv[loc], "-colvel") == 0) {
                    vecPtr = &colvel;
                    loc++;

                } else if(strcmp(argv[loc], "-beamvel") == 0) {
                    vecPtr = &beamvel;
                    loc++;
                } else if(strcmp(argv[loc], "-roofmass") == 0) {
                    vecPtr = &roofmass;
                    loc++;
                }

                if(vecPtr != 0) {
                    int num = 0;
                    const char** argvPtr = 0;
                    if(Tcl_SplitList(interp, argv[loc], &num, &argvPtr)!=TCL_OK) {
                        opserr<<"WARNING: failed to read the list "<<argv[loc];
                        opserr<<" -- PFEM2D discretize\n";
                        return TCL_ERROR;
                    }
                    if(num > 0) vecPtr->resize(num);
                    for(int i=0; i<num; i++) {
                        double x;
                        if(Tcl_GetDouble(interp, argvPtr[i], &x) != TCL_OK) {
                            opserr<<"WARNING: invalid input "<<argvPtr[i];
                            opserr<<" -- PFEM2D discretize\n";
                            return TCL_ERROR; 
                        }
                        (*vecPtr)(i) = x;
                    }

                    if(argvPtr != 0) Tcl_Free((char *) argvPtr);
                }

                loc++;
            }

            ID nodelist;
            int res = theMesher2D.discretize(startnode,x,y,hcol,hbeam,ndf,
                                             spans,heights,fix,colvel,
                                             beamvel,colmass,beammass,
                                             roofmass,theDomain,nodelist);
            if(res < 0) {
                opserr << "WARNING: failed to discretize\n";
                return TCL_ERROR;
            }

            for(int i=0; i<nodelist.Size(); i++) {
                char buffer[100];
                sprintf(buffer, "%d ", nodelist(i));
                Tcl_AppendResult(interp, buffer, NULL);
            }

        } else {
            opserr<<"WARNING: discretization type "<<argv[1]<<" is unknown ";
            opserr<<"-- PFEM2D discretize\n";
            return TCL_ERROR;
        }

    } else if(strcmp(argv[1], "doTriangulation") == 0) {
        if(argc < 4) {
            opserr << "WARNING: wrong num of args -- ";
            opserr << "PFEM2D doTriangulation alpha -groups {regTag1 ...}  ";
            opserr << "-addGroups {regTag2 ...} ";
            opserr << "-PFEMElement2D {starteletag rho mu b1 b2 <thk kappa>} ";
            opserr << "-PFEMElement2DCompressible {starteletag rho mu b1 b2 <thk kappa>} ";
            opserr << "-PFEMElement2DBubble {starteletag rho mu b1 b2 <thk kappa>} ";
            opserr << "-Tri31 {starteletag thk type matTag <pressure rho b1 b2>}\n ";
            return TCL_ERROR;
        }
        double alpha;
        if(Tcl_GetDouble(interp, argv[2], &alpha) != TCL_OK) {
            opserr<<"WARNING: invalid alpha "<<argv[2]<<" -- PFEM2D doTriangulation\n";
            return TCL_ERROR; 
        }

        int loc = 3;
        Vector params;
        int eletype = 0;
        std::string type;
        ID nodes, addnodes;
        while(loc < argc) {

            Vector* vecPtr = 0;
            ID* idPtr = 0;

            if(strcmp(argv[loc], "-groups") == 0) {
                idPtr = &nodes;
                loc++;

            } else if(strcmp(argv[loc], "-addgroups") == 0) {
                idPtr = &addnodes;
                loc++;
  
            } else if(strcmp(argv[loc], "-PFEMElement2D") == 0) {
                vecPtr = &params;
                eletype = 1;
                loc++;

            } else if(strcmp(argv[loc], "-Tri31") == 0 || strcmp(argv[loc], "-tri31") == 0) {
                vecPtr = &params;
                eletype = 2;
                loc++;
            } else if(strcmp(argv[loc], "-PFEMElement2DCompressible") == 0) {
                vecPtr = &params;
                eletype = 3;
                loc++;
            } else if(strcmp(argv[loc], "-PFEMElement2DBubble") == 0) {
                vecPtr = &params;
                eletype = 4;
                loc++;
            }

            if(vecPtr!=0 || idPtr!=0) {
                int num = 0;
                const char** argvPtr = 0;
                if(Tcl_SplitList(interp, argv[loc], &num, &argvPtr)!=TCL_OK) {
                    opserr<<"WARNING: failed to read the list "<<argv[loc];
                    opserr<<" -- PFEM2D doTriangulation\n";
                    return TCL_ERROR;
                }
                if(num > 0) {
                    if(vecPtr!=0) vecPtr->resize(num);
                    if(idPtr!=0) idPtr->resize(num);
                }
                if(vecPtr != 0) {
                    for(int i=0; i<num; i++) {
                        if(eletype==2 && i==2) {
                            type = argvPtr[i];
                            continue;
                        }
                        double x;
                        if(Tcl_GetDouble(interp, argvPtr[i], &x) != TCL_OK) {
                            opserr<<"WARNING: invalid input "<<argvPtr[i];
                            opserr<<" -- PFEM2D doTriangulation\n";
                            return TCL_ERROR; 
                        }
                        (*vecPtr)(i) = x;
                    }
                } else if(idPtr != 0) {
                    for(int i=0; i<num; i++) {
                        int x;
                        if(Tcl_GetInt(interp, argvPtr[i], &x) != TCL_OK) {
                            opserr<<"WARNING: invalid input "<<argvPtr[i];
                            opserr<<" -- PFEM2D doTriangulation\n";
                            return TCL_ERROR; 
                        }
                        (*idPtr)(i) = x;
                    }
                }

                if(argvPtr != 0) Tcl_Free((char *) argvPtr);
            }

            loc++;
        }

        // triangulation
        int res = -1, endele;
        ID eles;
        //opserr<<"params = "<<params;
        if(eletype>=1 && eletype<=7 && eletype!=2) {
            double rho,mu,b1,b2,thk=1.0,kappa=2.2e9;
            int startele;
            if(params.Size() >= 5) {
                startele = (int)params(0);
                rho = params(1);
                mu = params(2);
                b1 = params(3);
                b2 = params(4);
            } else {
                opserr<<"WARNING: not sufficient args ";
                opserr<<" -- PFEM2D doTriangulation\n";
                return TCL_ERROR;
            }
            if(params.Size() > 5) {
                thk = params(5);
            }
            if(params.Size() > 6) {
                kappa = params(6);
            }
            
            res = theMesher2D.doTriangulation(startele,alpha,nodes,
                                              addnodes,theDomain,
                                              rho,mu,b1,b2,thk,kappa,
                                              eletype,endele);

        } else if(eletype == 2) {
            double thk=0, p=0,rho=0,b1=0,b2=0;
            int startele=0, matTag=0;
            if(params.Size() >= 4) {
                startele = (int)params(0);
                thk = params(1);
                matTag = (int)params(3);
            } else {
                opserr<<"WARNING: not sufficient args ";
                opserr<<" -- PFEM2D doTriangulation\n";
                return TCL_ERROR;
            }
            if(params.Size() > 4) {
                p = params(4);
            }
            if(params.Size() > 5) {
                rho = params(5);
            }
            if(params.Size() > 6) {
                b1 = params(6);
            }
            if(params.Size() > 7) {
                b2 = params(7);
            }
            
            if(params.Size() >= 4) {
                res = theMesher2D.doTriangulation(startele,alpha,nodes,addnodes,
                                                  theDomain,
                                                  thk,type.c_str(),matTag,p,
                                                  rho,b1,b2,endele);
            }

            if(res < 0) {
                opserr<<"WARNING: failed to do triangulation -- ";
                opserr<<" -- PFEM2D doTriangulation\n";
                return TCL_ERROR; 
            }

        } else {
            res = theMesher2D.doTriangulation(alpha,nodes,addnodes,theDomain,eles);
            if(res < 0) {
                opserr<<"WARNING: failed to do triangulation -- ";
                opserr<<" -- PFEM2D doTriangulation\n";
                return TCL_ERROR; 
            }

            for(int i=0; i<eles.Size(); i++) {
                char buffer[100];
                sprintf(buffer, "%d ", eles(i));
                Tcl_AppendResult(interp, buffer, NULL);
            }
        }

        

        if(eletype>=1 && eletype<=4) {
            Tcl_SetObjResult(interp, Tcl_NewIntObj(endele));
        }

    } else if(strcmp(argv[1], "save") == 0) {

        if(argc < 3) {
            opserr << "WARNING: wrong num of args -- ";
            opserr << "PFEM2D save filename <maxelenodes>";
            return TCL_ERROR;            
        }
        int maxelenodes = 3;
        if(argc > 3) {
            if(Tcl_GetInt(interp, argv[3], &maxelenodes) != TCL_OK) {
                opserr<<"WANRING: invalid integer "<<argv[3];
                opserr<<" -- PFEM2D save\n";
                return TCL_ERROR; 
            }
        }

        theMesher2D.save(argv[2], theDomain, maxelenodes);

    } else if(strcmp(argv[1], "removeOutBoundNodes") == 0) {

        if(argc < 2) {
            opserr << "WARNING: wrong num of args -- ";
            opserr << "PFEM2D removeOutBoundNodes x1 y1 x2 y2 {groupList}\n";
            return TCL_ERROR;            
        }

        double x1,y1,x2,y2;
        if(Tcl_GetDouble(interp, argv[2], &x1) != TCL_OK) {
            opserr<<"WARNING: invalid node x1 "<<argv[2];
            opserr<<" -- PFEM2D removeOutBoundNodes\n";
            return TCL_ERROR; 
        }
        if(Tcl_GetDouble(interp, argv[3], &y1) != TCL_OK) {
            opserr<<"WARNING: invalid node y1 "<<argv[3];
            opserr<<" -- PFEM2D removeOutBoundNodes\n";
            return TCL_ERROR; 
        }
        if(Tcl_GetDouble(interp, argv[4], &x2) != TCL_OK) {
            opserr<<"WARNING: invalid node x2 "<<argv[4];
            opserr<<" -- PFEM2D removeOutBoundNodes\n";
            return TCL_ERROR; 
        }
        if(Tcl_GetDouble(interp, argv[5], &y2) != TCL_OK) {
            opserr<<"WARNING: invalid node y2 "<<argv[5];
            opserr<<" -- PFEM2D removeOutBoundNodes\n";
            return TCL_ERROR; 
        }

        ID groups;
        ID* idPtr = &groups;
        int num = 0, loc = 6;
        const char** argvPtr = 0;
        if(Tcl_SplitList(interp, argv[loc], &num, &argvPtr)!=TCL_OK) {
            opserr<<"WARNING: failed to read the list "<<argv[loc];
            opserr<<" -- PFEM2D removeOutBoundNodes\n";
            return TCL_ERROR;
        }
        if(num > 0) idPtr->resize(num);
        for(int i=0; i<num; i++) {
            int tag;
            if(Tcl_GetInt(interp, argvPtr[i], &tag) != TCL_OK) {
                opserr<<"WARNING: invalid input "<<argvPtr[i];
                opserr<<" -- PFEM2D removeOutBoundNodes\n";
                return TCL_ERROR; 
            }
            (*idPtr)(i) = tag;
        }
        
        if(argvPtr != 0) Tcl_Free((char *) argvPtr);
        theMesher2D.removeOutBoundNodes(groups, x1,y1,x2,y2, theDomain);

    } else if(strcmp(argv[1], "calculateForces") == 0) {

        if(argc < 2) {
            opserr<<"WARNING: wrong num of args -- ";
            opserr<<"PFEM2D calculateForces -boundary {$nd1 ...} ";
            opserr<<"-basenode $nd -dragdir {$dx $dy} -liftdir {$dx $dy}"; 
            return TCL_ERROR;            
        }
        int loc = 2;
        Vector drag, lift;
        ID boundary;
        int basenode;
        while(loc < argc) {
            Vector* vecPtr = 0;
            ID* idPtr = 0;

            if(strcmp(argv[loc], "-boundary") == 0) {
                idPtr = &boundary;
                loc++;

            } else if(strcmp(argv[loc], "-dragdir") == 0) {
                vecPtr = &drag;
                loc++;
  
            } else if(strcmp(argv[loc], "-liftdir") == 0) {
                vecPtr = &lift;
                loc++;
            } else if(strcmp(argv[loc], "-basenode") == 0) {
                loc++;
                if(Tcl_GetInt(interp, argv[loc], &basenode) != TCL_OK) {
                    opserr<<"WARNING: invalid input "<<argv[loc];
                    opserr<<" -- PFEM2D calculateForces\n";
                    return TCL_ERROR; 
                }
                loc++;
            }

            if(vecPtr!=0 || idPtr!=0) {
                int num = 0;
                const char** argvPtr = 0;
                if(Tcl_SplitList(interp, argv[loc], &num, &argvPtr)!=TCL_OK) {
                    opserr<<"WARNING: failed to read the list "<<argv[loc];
                    opserr<<" -- PFEM2D calculateForces\n";
                    return TCL_ERROR;
                }
                if(num > 0) {
                    if(vecPtr!=0) vecPtr->resize(num);
                    if(idPtr!=0) idPtr->resize(num);
                }
                if(vecPtr != 0) {
                    for(int i=0; i<num; i++) {
                        double x;
                        if(Tcl_GetDouble(interp, argvPtr[i], &x) != TCL_OK) {
                            opserr<<"WARNING: invalid input "<<argvPtr[i];
                            opserr<<" -- PFEM2D calculateForces\n";
                            return TCL_ERROR; 
                        }
                        (*vecPtr)(i) = x;
                    }
                } else if(idPtr != 0) {
                    for(int i=0; i<num; i++) {
                        int x;
                        if(Tcl_GetInt(interp, argvPtr[i], &x) != TCL_OK) {
                            opserr<<"WARNING: invalid input "<<argvPtr[i];
                            opserr<<" -- PFEM2D calculateForce\n";
                            return TCL_ERROR; 
                        }
                        (*idPtr)[i] = x;
                    }
                }

                if(argvPtr != 0) Tcl_Free((char *) argvPtr);
            }

            loc++;
        }

        Vector forces = theMesher2D.calculateForces(boundary,basenode,drag,lift,theDomain);

        for(int i=0; i<forces.Size(); i++) {
            char buffer[60];
            sprintf(buffer, "%35.20f ", forces(i));
            Tcl_AppendResult(interp, buffer, NULL);
        }
    } else if(strcmp(argv[1], "pc") == 0 || strcmp(argv[1], "PC") == 0) {
        if(argc < 5) {
            opserr<<"WARNING: wrong num of args -- ";
            opserr<<"PFEM2D pc endNodeList startpressurenode ndf\n";
            return TCL_ERROR;
        }

        ID nodes;
        int num = 0;
        const char** argvPtr = 0;
        if(Tcl_SplitList(interp, argv[2], &num, &argvPtr)!=TCL_OK) {
            opserr<<"WARNING: failed to read the list "<<argv[2];
            opserr<<" -- PFEM2D pc\n";
            return TCL_ERROR;
        }
        if(num > 0) nodes.resize(num);
        for(int i=0; i<num; i++) {
            int tag;
            if(Tcl_GetInt(interp, argvPtr[i], &tag) != TCL_OK) {
                opserr<<"WARNING: invalid input "<<argvPtr[i];
                opserr<<" -- PFEM2D pc\n";
                return TCL_ERROR; 
            }
            nodes(i) = tag;
        }

        if(argvPtr != 0) Tcl_Free((char *) argvPtr);

        int startpnode;
        if(Tcl_GetInt(interp, argv[3], &startpnode) != TCL_OK) {
            opserr<<"WARNING: invalid startpnode "<<argv[3];
            opserr<<" -- PFEM2D pc\n";
            return TCL_ERROR; 
        }

        int ndf;
        if(Tcl_GetInt(interp, argv[4], &ndf) != TCL_OK) {
            opserr<<"WARNING: invalid ndf "<<argv[4];
            opserr<<" -- PFEM2D pc\n";
            return TCL_ERROR; 
        }

        int endpnode;
        int res = theMesher2D.addPC(nodes,ndf,startpnode,theDomain,endpnode);

        if(res < 0) {
            opserr << "WARNING: failed to add PC\n";
            return TCL_ERROR;
        }
        
        Tcl_SetObjResult(interp, Tcl_NewIntObj(endpnode));

    } else if(strcmp(argv[1], "group") == 0 || strcmp(argv[1], "Group") == 0) {
        if(argc < 3) {
            opserr<<"WARNING: wrong num of args -- ";
            opserr<<"PFEM2D group regTag nodelist <-append -remove -series>\n";
            return TCL_ERROR;
        }

        int type = -1;
        if(Tcl_GetInt(interp, argv[2], &type) != TCL_OK) {
            opserr<<"WARNING: invalid type "<<argv[2];
            opserr<<" -- PFEM2D group \n";
            return TCL_ERROR;
        }

        ID nodes;
        int num = 0;
        const char** argvPtr = 0;
        if(Tcl_SplitList(interp, argv[3], &num, &argvPtr)!=TCL_OK) {
            opserr<<"WARNING: failed to read the list "<<argv[3];
            opserr<<" -- PFEM2D "<<argv[1]<<"\n";
            return TCL_ERROR;
        }
        if(num > 0) nodes.resize(num);
        for(int i=0; i<num; i++) {
            int tag;
            if(Tcl_GetInt(interp, argvPtr[i], &tag) != TCL_OK) {
                opserr<<"WARNING: invalid input "<<argvPtr[i];
                opserr<<" -- PFEM2D "<<argv[1]<<"\n";
                return TCL_ERROR; 
            }
            nodes(i) = tag;
        }

        if(argvPtr != 0) Tcl_Free((char *) argvPtr);

        int action = 1;        
        bool series = false;
        for(int i=4; i<argc; i++) {
            if(strcmp(argv[i], "-append")==0 || strcmp(argv[i], "-Append")==0) {
                action = 1;
            } else if(strcmp(argv[i], "-remove")==0 || strcmp(argv[i], "-Remove")==0) {
                action = 2;
            } else if(strcmp(argv[i], "-series")==0 || strcmp(argv[i], "-Series")==0) {
                series = true;
            }
        }

        theMesher2D.setNodes(nodes,type,series,action,theDomain);

    } else if(strcmp(argv[1], "updateNode") == 0 || strcmp(argv[1], "updatenode") == 0) {        
        if(argc < 6) {
            opserr<<"Warning: wrong num of args -- ";
            opserr<<"PFEM2D updateNode /coord/disp/vel/accel tag dof val\n";
            return TCL_ERROR;
        }
        int tag, dof, type;
        double value;
        if(strcmp(argv[2],"disp")==0 || strcmp(argv[2],"Disp")==0) {
            type = 1;
        } else if(strcmp(argv[2],"vel")==0 || strcmp(argv[2],"Vel")==0) {
            type = 2;
        } else if(strcmp(argv[2],"accel")==0 || strcmp(argv[2],"Accel")==0) {
            type = 3;
        } else {
            opserr<<"WARNING: unknown type for updating -- ";
            opserr<<"PFEM2D updateNode coord/disp/vel/accel tag dof val\n";
            return -1;
        }
        if(Tcl_GetInt(interp, argv[3], &tag) != TCL_OK) {
            opserr<<"WARNING: invalid input "<<argv[3];
            opserr<<" -- PFEM2D updateNode\n";
            return TCL_ERROR; 
        }
        if(Tcl_GetInt(interp, argv[4], &dof) != TCL_OK) {
            opserr<<"WARNING: invalid input "<<argv[4];
            opserr<<" -- PFEM2D updateNode\n";
            return TCL_ERROR; 
        }
        if(Tcl_GetDouble(interp, argv[5], &value) != TCL_OK) {
            opserr<<"WARNING: invalid input "<<argv[5];
            opserr<<" -- PFEM2D updateNode\n";
            return TCL_ERROR; 
        }
        theMesher2D.updateNode(tag,dof,value,type,theDomain);

    } else if(strcmp(argv[1], "region") == 0 || 
              strcmp(argv[1], "Region") == 0) {

        if(strcmp(argv[2], "PSLG") == 0) {
            if(argc < 6) {
                opserr << "WARNING: wrong num of args -- ";
                opserr << "PFEM2D region PSLG regTag maxarea ndf ";
                opserr << "<-points {x1 y1 ...} -segments {p1 p2 ...} ";
                opserr << "-vel {vx vy} ";
                opserr << "-fix {constrValues} -mass {massValues} ";
                opserr << "-holes {x1 y1 ...} >\n";
                return TCL_ERROR;
            }
            int regTag;
            if(Tcl_GetInt(interp, argv[3], &regTag) != TCL_OK) {
                opserr<<"WARNING: invalid region tag "<<argv[3];
                opserr<<" -- PFEM2D region PSLG\n";
                return TCL_ERROR; 
            }
            double maxarea;
            if(Tcl_GetDouble(interp, argv[4], &maxarea) != TCL_OK) {
                opserr<<"WARNING: invalid maxarea "<<argv[4];
                opserr<<" -- PFEM2D region PSLG\n";
                return TCL_ERROR; 
            }
            if(maxarea <= 0) {
                opserr<<"WARNING: nonpositive area ";
                opserr<<" -- PFEM2D region PSLG\n";
                return TCL_ERROR;
            }
            int ndf;
            if(Tcl_GetInt(interp, argv[5], &ndf) != TCL_OK) {
                opserr<<"WARNING: invalid ndf "<<argv[5];
                opserr<<" -- PFEM2D region PSLG\n";
                return TCL_ERROR; 
            }
            if(ndf <= 0) {
                opserr<<"WARNING: nonpositive ndf ";
                opserr<<" -- PFEM2D region PSLG\n";
                return TCL_ERROR;
            }

            int loc = 6;
            Vector points, mass, holes, vel, segments, fix;
            while(loc < argc) {

                Vector* vecPtr = 0;
                if(strcmp(argv[loc], "-points") == 0) {
                    vecPtr = &points;
                    loc++;

                } else if(strcmp(argv[loc], "-segments") == 0) {
                    vecPtr = &segments;
                    loc++;

                } else if(strcmp(argv[loc], "-mass") == 0) {
                    vecPtr = &mass;
                    loc++;

                } else if(strcmp(argv[loc], "-fix") == 0) {
                    vecPtr = &fix;
                    loc++;
  
                } else if(strcmp(argv[loc], "-holes") == 0) {
                    vecPtr = &holes;
                    loc++;

                } else if(strcmp(argv[loc], "-vel") == 0) {
                    vecPtr = &vel;
                    loc++;
                }

                if(vecPtr != 0) {
                    int num = 0;
                    const char** argvPtr = 0;
                    if(Tcl_SplitList(interp, argv[loc], &num, &argvPtr)
                       !=TCL_OK) {
                        opserr<<"WARNING: failed to read the list "<<argv[loc];
                        opserr<<" -- PFEM2D region PSLG\n";
                        return TCL_ERROR;
                    }
                    if(num > 0) vecPtr->resize(num);
                    for(int i=0; i<num; i++) {
                        double x;
                        if(Tcl_GetDouble(interp, argvPtr[i], &x) != TCL_OK) {
                            opserr<<"WARNING: invalid input "<<argvPtr[i];
                            opserr<<" -- PFEM2D region PSLG\n";
                            return TCL_ERROR; 
                        }
                        (*vecPtr)(i) = x;
                    }

                    if(argvPtr != 0) Tcl_Free((char *) argvPtr);
                }

                loc++;
            }

            int startnodetag = theMesher2D.findNodeTag(theDomain);
            int endnodetag = startnodetag;
            int res = theMesher2D.discretize(startnodetag,points,segments,
                                             holes,maxarea,ndf,
                                             fix,vel,mass,theDomain, 
                                             endnodetag);
            if(res < 0) {
                opserr << "WARNING: failed to discretize";
                opserr<<" -- PFEM2D region PSLG\n";
                return TCL_ERROR;
            }

            // add region
            ID regionnodes(2);
            regionnodes(0) = startnodetag;
            regionnodes(1) = endnodetag;
            bool series = false;
            int action = 1; // append
            theMesher2D.setNodes(regionnodes, regTag, series, action, 
                                 theDomain);

        } else if(strcmp(argv[2], "particles") == 0 || 
                  strcmp(argv[2], "Particles") == 0) {
            if(argc < 5) {
                opserr << "WARNING: wrong num of args -- ";
                opserr << "PFEM2D region particles regTag ndf ";
                opserr << "<-points {x1 y1 ...} -vel {vx vy} ";
                opserr << "-fix {constrValues} -mass {massValues} >\n";
                return TCL_ERROR;
            }
            int regTag;
            if(Tcl_GetInt(interp, argv[3], &regTag) != TCL_OK) {
                opserr<<"WARNING: invalid region tag "<<argv[3]<<" -- PFEM2D region\n";
                return TCL_ERROR; 
            }
            int ndf;
            if(Tcl_GetInt(interp, argv[4], &ndf) != TCL_OK) {
                opserr<<"WARNING: invalid ndf "<<argv[4]<<" -- PFEM2D region\n";
                return TCL_ERROR; 
            }
            if(ndf <= 0) {
                opserr<<"WARNING: nonpositive ndf ";
                opserr<<" -- PFEM2D region\n";
                return TCL_ERROR;
            }

            int loc = 5;
            Vector points, mass, vel, fix;
            while(loc < argc) {

                Vector* vecPtr = 0;
                if(strcmp(argv[loc], "-points") == 0) {
                    vecPtr = &points;
                    loc++;

                } else if(strcmp(argv[loc], "-mass") == 0) {
                    vecPtr = &mass;
                    loc++;

                } else if(strcmp(argv[loc], "-fix") == 0) {
                    vecPtr = &fix;
                    loc++;

                } else if(strcmp(argv[loc], "-vel") == 0) {
                    vecPtr = &vel;
                    loc++;
                }

                if(vecPtr != 0) {
                    int num = 0;
                    const char** argvPtr = 0;
                    if(Tcl_SplitList(interp, argv[loc], &num, &argvPtr)!=TCL_OK) {
                        opserr<<"WARNING: failed to read the list "<<argv[loc];
                        opserr<<" -- PFEM2D region\n";
                        return TCL_ERROR;
                    }
                    if(num > 0) vecPtr->resize(num);
                    for(int i=0; i<num; i++) {
                        double x;
                        if(Tcl_GetDouble(interp, argvPtr[i], &x) != TCL_OK) {
                            opserr<<"WARNING: invalid input "<<argvPtr[i];
                            opserr<<" -- PFEM2D region\n";
                            return TCL_ERROR; 
                        }
                        (*vecPtr)(i) = x;
                    }

                    if(argvPtr != 0) Tcl_Free((char *) argvPtr);
                }

                loc++;
            }
            
            int startnodetag = theMesher2D.findNodeTag(theDomain);
            int endnodetag = startnodetag;
            int res = theMesher2D.discretize(startnodetag,points,ndf,
                                             fix,vel,mass,theDomain, 
                                             endnodetag);
            
            if(res < 0) {
                opserr << "WARNING: failed to discretize\n";
                return TCL_ERROR;
            }

            // add region
            ID regionnodes(2);
            regionnodes(0) = startnodetag;
            regionnodes(1) = endnodetag;
            bool series = false;
            int action = 1; // append
            theMesher2D.setNodes(regionnodes, regTag, series, action, 
                                 theDomain);

        } else if(strcmp(argv[2], "rectangle")==0 || strcmp(argv[2], "Rectangle")==0) {
            if(argc < 12) {
                opserr<<"WARNING: wrong num of args -- ";
                opserr<<"PFEM2D region Rectangle regTag x1 y1 hx hy angle nx ny ndf ";
                opserr<<"<-fix {constrValues} -mass {massValues} -vel {vx vy}";
                opserr<<" -boundary {boundaryValues}>\n";
                return TCL_ERROR;
            }
            double x1,y1,hx,hy,angle;
            int nx, ny, ndf, regTag;
            if(Tcl_GetInt(interp, argv[3], &regTag) != TCL_OK) {
                opserr<<"WARNING: invalid region tag "<<argv[3]<<" -- PFEM2D region Rectangle\n";
                return TCL_ERROR; 
            }
            if(Tcl_GetDouble(interp, argv[4], &x1) != TCL_OK) {
                opserr<<"WARNING: invalid x1 "<<argv[4]<<" -- PFEM2D region Rectangle\n";
                return TCL_ERROR; 
            }
            if(Tcl_GetDouble(interp, argv[5], &y1) != TCL_OK) {
                opserr<<"WARNING: invalid y1 "<<argv[5]<<" -- PFEM2D region Rectangle\n";
                return TCL_ERROR; 
            }
            if(Tcl_GetDouble(interp, argv[6], &hx) != TCL_OK) {
                opserr<<"WARNING: invalid hx "<<argv[6]<<" -- PFEM2D region Rectangle\n";
                return TCL_ERROR; 
            }
            if(Tcl_GetDouble(interp, argv[7], &hy) != TCL_OK) {
                opserr<<"WARNING: invalid hy "<<argv[7]<<" -- PFEM2D region Rectangle\n";
                return TCL_ERROR; 
            }
            if(Tcl_GetDouble(interp, argv[8], &angle) != TCL_OK) {
                opserr<<"WARNING: invalid angle "<<argv[8]<<" -- PFEM2D region Rectangle\n";
                return TCL_ERROR; 
            }
            if(Tcl_GetInt(interp, argv[9], &nx) != TCL_OK) {
                opserr<<"WARNING: invalid nx "<<argv[9]<<" -- PFEM2D region Rectangle\n";
                return TCL_ERROR; 
            }
            if(Tcl_GetInt(interp, argv[10], &ny) != TCL_OK) {
                opserr<<"WARNING: invalid ny "<<argv[10]<<" -- PFEM2D region Rectangle\n";
                return TCL_ERROR; 
            }
            if(Tcl_GetInt(interp, argv[11], &ndf) != TCL_OK) {
                opserr<<"WARNING: invalid ndf "<<argv[11]<<" -- PFEM2D region Rectangle\n";
                return TCL_ERROR; 
            }

            int loc = 12;
            Vector mass, vel, fix, boundary(4);
            boundary += 1;
            while(loc < argc) {

                Vector* vecPtr = 0;

                if(strcmp(argv[loc], "-mass") == 0) {
                    vecPtr = &mass;
                    loc++;

                } else if(strcmp(argv[loc], "-fix") == 0) {
                    vecPtr = &fix;
                    loc++;
  
                } else if(strcmp(argv[loc], "-vel") == 0) {
                    vecPtr = &vel;
                    loc++;
                } else if(strcmp(argv[loc], "-boundary") == 0) {
                    vecPtr = &boundary;
                    loc++;
                }

                if(vecPtr != 0) {
                    int num = 0;
                    const char** argvPtr = 0;
                    if(Tcl_SplitList(interp, argv[loc], &num, &argvPtr)!=TCL_OK) {
                        opserr<<"WARNING: failed to read the list "<<argv[loc];
                        opserr<<" -- PFEM2D region\n";
                        return TCL_ERROR;
                    }
                    if(num > 0) vecPtr->resize(num);
                    for(int i=0; i<num; i++) {
                        double x;
                        if(Tcl_GetDouble(interp, argvPtr[i], &x) != TCL_OK) {
                            opserr<<"WARNING: invalid input "<<argvPtr[i];
                            opserr<<" -- PFEM2D region\n";
                            return TCL_ERROR; 
                        }
                        (*vecPtr)(i) = x;
                    }

                    if(argvPtr != 0) Tcl_Free((char *) argvPtr);
                }

                loc++;
            }

            int startnodetag = theMesher2D.findNodeTag(theDomain);
            int endnodetag = startnodetag;
            int res = theMesher2D.discretize(startnodetag,x1,y1,hx,hy,
                                             angle,nx,ny,ndf,fix,vel,mass,
                                             boundary,theDomain,endnodetag);

            if(res < 0) {
                opserr << "WARNING: failed to discretize\n";
                return TCL_ERROR;
            }
 
            // add region
            ID regionnodes(2);
            regionnodes(0) = startnodetag;
            regionnodes(1) = endnodetag;
            bool series = false;
            int action = 1; // append
            theMesher2D.setNodes(regionnodes, regTag, series, action, 
                                 theDomain);

        } else if(strcmp(argv[2], "line")==0 || strcmp(argv[2], "Line")==0) {
            if(argc < 10) {
                opserr<<"WARNING: wrong num of args -- ";
                opserr<<"PFEM2D region line regTag x1 y1 h angle num ndf ";
                opserr<<"-fix {constrValues} -mass {massValues} -vel {vx vy}";
                opserr<<"-boundary {boundaryValues}\n";
                return TCL_ERROR;
            }
            double x1,y1,h,angle;
            int num, ndf, regTag;
            if(Tcl_GetInt(interp, argv[3], &regTag) != TCL_OK) {
                opserr<<"WARNING: invalid region tag "<<argv[3]<<" -- PFEM2D region Rectangle\n";
                return TCL_ERROR; 
            }
            if(Tcl_GetDouble(interp, argv[4], &x1) != TCL_OK) {
                opserr<<"WARNING: invalid x1 "<<argv[4]<<" -- PFEM2D region Rectangle\n";
                return TCL_ERROR; 
            }
            if(Tcl_GetDouble(interp, argv[5], &y1) != TCL_OK) {
                opserr<<"WARNING: invalid y1 "<<argv[5]<<" -- PFEM2D region Rectangle\n";
                return TCL_ERROR; 
            }
            if(Tcl_GetDouble(interp, argv[6], &h) != TCL_OK) {
                opserr<<"WARNING: invalid h "<<argv[6]<<" -- PFEM2D region Rectangle\n";
                return TCL_ERROR; 
            }
            if(Tcl_GetDouble(interp, argv[7], &angle) != TCL_OK) {
                opserr<<"WARNING: invalid angle "<<argv[7]<<" -- PFEM2D region Rectangle\n";
                return TCL_ERROR; 
            }
            if(Tcl_GetInt(interp, argv[8], &num) != TCL_OK) {
                opserr<<"WARNING: invalid num "<<argv[8]<<" -- PFEM2D region Rectangle\n";
                return TCL_ERROR; 
            }
            if(Tcl_GetInt(interp, argv[9], &ndf) != TCL_OK) {
                opserr<<"WARNING: invalid ndf "<<argv[9]<<" -- PFEM2D region Rectangle\n";
                return TCL_ERROR; 
            }

            int loc = 10;
            Vector mass, vel, fix, boundary(2);
            boundary += 1;
            while(loc < argc) {

                Vector* vecPtr = 0;

                if(strcmp(argv[loc], "-mass") == 0) {
                    vecPtr = &mass;
                    loc++;

                } else if(strcmp(argv[loc], "-fix") == 0) {
                    vecPtr = &fix;
                    loc++;
  
                } else if(strcmp(argv[loc], "-vel") == 0) {
                    vecPtr = &vel;
                    loc++;
                } else if(strcmp(argv[loc], "-boundary") == 0) {
                    vecPtr = &boundary;
                    loc++;
                }

                if(vecPtr != 0) {
                    int num = 0;
                    const char** argvPtr = 0;
                    if(Tcl_SplitList(interp, argv[loc], &num, &argvPtr)!=TCL_OK) {
                        opserr<<"WARNING: failed to read the list "<<argv[loc];
                        opserr<<" -- PFEM2D region\n";
                        return TCL_ERROR;
                    }
                    if(num > 0) vecPtr->resize(num);
                    for(int i=0; i<num; i++) {
                        double x;
                        if(Tcl_GetDouble(interp, argvPtr[i], &x) != TCL_OK) {
                            opserr<<"WARNING: invalid input "<<argvPtr[i];
                            opserr<<" -- PFEM2D region\n";
                            return TCL_ERROR; 
                        }
                        (*vecPtr)(i) = x;
                    }

                    if(argvPtr != 0) Tcl_Free((char *) argvPtr);
                }

                loc++;
            }

            int startnodetag = theMesher2D.findNodeTag(theDomain);
            int endnodetag = startnodetag;
            int res = theMesher2D.discretize(startnodetag, x1,y1,h,angle,num,
                                             ndf,fix,vel,mass,
                                             boundary,theDomain,endnodetag);

            if(res < 0) {
                opserr << "WARNING: failed to discretize\n";
                return TCL_ERROR;
            }

            // add region
            ID regionnodes(2);
            regionnodes(0) = startnodetag;
            regionnodes(1) = endnodetag;
            bool series = false;
            int action = 1; // append
            theMesher2D.setNodes(regionnodes, regTag, series, action, 
                                 theDomain);

        } else if(strcmp(argv[2], "triangle")==0 || strcmp(argv[2], "Triangle")==0) {

            if(argc < 13) {
                opserr<<"WARNING: wrong num of args -- ";
                opserr<<"PFEM2D region triangle regTag x1 y1 x2 y2 x3 y3 n1 n2 ndf ";
                opserr<<"<-fix {constrValues} -mass {massValues} -vel {vx vy} ";
                opserr<<"-boundary {boundaryValues}>\n";
                return TCL_ERROR;
            }
            double x1,y1,x2,y2,x3,y3;
            int n1,n2,ndf,regTag;
            if(Tcl_GetInt(interp, argv[3], &regTag) != TCL_OK) {
                opserr<<"WARNING: invalid region tag "<<argv[3]<<" -- PFEM2D region triangle\n";
                return TCL_ERROR; 
            }
            if(Tcl_GetDouble(interp, argv[4], &x1) != TCL_OK) {
                opserr<<"WARNING: invalid x1 "<<argv[4]<<" -- PFEM2D region triangle\n";
                return TCL_ERROR; 
            }
            if(Tcl_GetDouble(interp, argv[5], &y1) != TCL_OK) {
                opserr<<"WARNING: invalid y1 "<<argv[5]<<" -- PFEM2D region triangle\n";
                return TCL_ERROR; 
            }
            if(Tcl_GetDouble(interp, argv[6], &x2) != TCL_OK) {
                opserr<<"WARNING: invalid x2 "<<argv[6]<<" -- PFEM2D region triangle\n";
                return TCL_ERROR; 
            }
            if(Tcl_GetDouble(interp, argv[7], &y2) != TCL_OK) {
                opserr<<"WARNING: invalid y2 "<<argv[7]<<" -- PFEM2D region triangle\n";
                return TCL_ERROR; 
            }
            if(Tcl_GetDouble(interp, argv[8], &x3) != TCL_OK) {
                opserr<<"WARNING: invalid x3 "<<argv[8]<<" -- PFEM2D region triangle\n";
                return TCL_ERROR; 
            }
            if(Tcl_GetDouble(interp, argv[9], &y3) != TCL_OK) {
                opserr<<"WARNING: invalid y3 "<<argv[9]<<" -- PFEM2D region triangle\n";
                return TCL_ERROR; 
            }
            if(Tcl_GetInt(interp, argv[10], &n1) != TCL_OK) {
                opserr<<"WARNING: invalid n1 "<<argv[10]<<" -- PFEM2D region triangle\n";
                return TCL_ERROR; 
            }
            if(Tcl_GetInt(interp, argv[11], &n2) != TCL_OK) {
                opserr<<"WARNING: invalid n2 "<<argv[11]<<" -- PFEM2D region triangle\n";
                return TCL_ERROR; 
            }
            if(Tcl_GetInt(interp, argv[12], &ndf) != TCL_OK) {
                opserr<<"WARNING: invalid ndf "<<argv[12]<<" -- PFEM2D region triangle\n";
                return TCL_ERROR; 
            }

            int loc = 13;
            Vector mass, vel, fix, boundary(3);
            boundary += 1;
            while(loc < argc) {

                Vector* vecPtr = 0;

                if(strcmp(argv[loc], "-mass") == 0) {
                    vecPtr = &mass;
                    loc++;

                } else if(strcmp(argv[loc], "-fix") == 0) {
                    vecPtr = &fix;
                    loc++;
  
                } else if(strcmp(argv[loc], "-vel") == 0) {
                    vecPtr = &vel;
                    loc++;
                } else if(strcmp(argv[loc], "-boundary") == 0) {
                    vecPtr = &boundary;
                    loc++;
                }

                if(vecPtr != 0) {
                    int num = 0;
                    const char** argvPtr = 0;
                    if(Tcl_SplitList(interp, argv[loc], &num, &argvPtr)!=TCL_OK) {
                        opserr<<"WARNING: failed to read the list "<<argv[loc];
                        opserr<<" -- PFEM2D region\n";
                        return TCL_ERROR;
                    }
                    if(num > 0) vecPtr->resize(num);
                    for(int i=0; i<num; i++) {
                        double x;
                        if(Tcl_GetDouble(interp, argvPtr[i], &x) != TCL_OK) {
                            opserr<<"WARNING: invalid input "<<argvPtr[i];
                            opserr<<" -- PFEM2D region\n";
                            return TCL_ERROR; 
                        }
                        (*vecPtr)(i) = x;
                    }

                    if(argvPtr != 0) Tcl_Free((char *) argvPtr);
                }

                loc++;
            }

            int startnode = theMesher2D.findNodeTag(theDomain);
            int endnode = startnode;

            int res = theMesher2D.discretize(startnode,x1,y1,x2,y2,x3,y3,n1,n2,
                                             ndf,fix,vel,mass,
                                             boundary,theDomain,endnode);
            if(res < 0) {
                opserr << "WARNING: failed to discretize\n";
                return TCL_ERROR;
            }

            // add region
            ID regionnodes(2);
            regionnodes(0) = startnode;
            regionnodes(1) = endnode;
            bool series = false;
            int action = 1; // append
            theMesher2D.setNodes(regionnodes, regTag, series, action, 
                                 theDomain);

        } else if(strcmp(argv[2], "circle")==0 || strcmp(argv[2], "Circle")==0) {

            if(argc < 11) {
                opserr<<"WARNING: wrong num of args -- ";
                opserr<<"PFEM2D region circle regTag xc yc r1 r2 nc nr ndf ";
                opserr<<"<-fix {constrValues} -mass {massValues} -vel {vx vy} ";
                opserr<<"-angleRange {start end} "; 
                opserr<<"-boundary {boundaryValues}>\n";
                return TCL_ERROR;
            }

            double xc, yc, r1, r2;
            int nc,nr,ndf,regTag;
            if(Tcl_GetInt(interp, argv[3], &regTag) != TCL_OK) {
                opserr<<"WARNING: invalid region tag "<<argv[3]<<" -- PFEM2D region circle\n";
                return TCL_ERROR; 
            }
            if(Tcl_GetDouble(interp, argv[4], &xc) != TCL_OK) {
                opserr<<"WARNING: invalid xc "<<argv[4]<<" -- PFEM2D region circle\n";
                return TCL_ERROR; 
            }
            if(Tcl_GetDouble(interp, argv[5], &yc) != TCL_OK) {
                opserr<<"WARNING: invalid yc "<<argv[5]<<" -- PFEM2D region circle\n";
                return TCL_ERROR; 
            }
            if(Tcl_GetDouble(interp, argv[6], &r1) != TCL_OK) {
                opserr<<"WARNING: invalid r1 "<<argv[6]<<" -- PFEM2D region circle\n";
                return TCL_ERROR; 
            }
            if(Tcl_GetDouble(interp, argv[7], &r2) != TCL_OK) {
                opserr<<"WARNING: invalid r2 "<<argv[7]<<" -- PFEM2D region circle\n";
                return TCL_ERROR; 
            }
            if(Tcl_GetInt(interp, argv[8], &nc) != TCL_OK) {
                opserr<<"WARNING: invalid nc "<<argv[8]<<" -- PFEM2D region circle\n";
                return TCL_ERROR; 
            }
            if(Tcl_GetInt(interp, argv[9], &nr) != TCL_OK) {
                opserr<<"WARNING: invalid nr "<<argv[9]<<" -- PFEM2D region circle\n";
                return TCL_ERROR; 
            }
            if(Tcl_GetInt(interp, argv[10], &ndf) != TCL_OK) {
                opserr<<"WARNING: invalid ndf "<<argv[10]<<" -- PFEM2D region circle\n";
                return TCL_ERROR; 
            }
            int loc = 11;
            Vector mass, vel, fix, boundary(2), angle(2);
            boundary += 1;
            while(loc < argc) {

                Vector* vecPtr = 0;

                if(strcmp(argv[loc], "-mass") == 0) {
                    vecPtr = &mass;
                    loc++;

                } else if(strcmp(argv[loc], "-fix") == 0) {
                    vecPtr = &fix;
                    loc++;
  
                } else if(strcmp(argv[loc], "-vel") == 0) {
                    vecPtr = &vel;
                    loc++;
                } else if(strcmp(argv[loc], "-boundary") == 0) {
                    vecPtr = &boundary;
                    loc++;
                } else if(strcmp(argv[loc], "-angleRange") == 0) {
                    vecPtr = &angle;
                    loc++;
                }

                if(vecPtr != 0) {
                    int num = 0;
                    const char** argvPtr = 0;
                    if(Tcl_SplitList(interp, argv[loc], &num, &argvPtr)!=TCL_OK) {
                        opserr<<"WARNING: failed to read the list "<<argv[loc];
                        opserr<<" -- PFEM2D region\n";
                        return TCL_ERROR;
                    }
                    if(num > 0) vecPtr->resize(num);
                    for(int i=0; i<num; i++) {
                        double x;
                        if(Tcl_GetDouble(interp, argvPtr[i], &x) != TCL_OK) {
                            opserr<<"WARNING: invalid input "<<argvPtr[i];
                            opserr<<" -- PFEM2D region\n";
                            return TCL_ERROR; 
                        }
                        (*vecPtr)(i) = x;
                    }

                    if(argvPtr != 0) Tcl_Free((char *) argvPtr);
                }

                loc++;
            }

            int startnode = theMesher2D.findNodeTag(theDomain);
            int endnode = startnode;

            int res = theMesher2D.discretize(startnode,xc,yc,r1,r2,nc,nr,ndf,fix,vel,mass,
                                             boundary,angle,theDomain,endnode);
            if(res < 0) {
                opserr << "WARNING: failed to discretize\n";
                return TCL_ERROR;
            }


            // add region
            ID regionnodes(2);
            regionnodes(0) = startnode;
            regionnodes(1) = endnode;
            bool series = false;
            int action = 1; // append
            theMesher2D.setNodes(regionnodes, regTag, series, action, 
                                 theDomain);

        } else if(strcmp(argv[2], "frame")==0 || strcmp(argv[2], "Frame")==0) {

            if(argc < 11) {
                opserr<<"WARNING: wrong num of args -- ";
                opserr<<"PFEM2D region frame regTag x y hcol hbeam ndf ";
                opserr<<"{spanlist} {hlist} -fix {constrValues} ";
                opserr<<"-colmass {masslist} -beammass {masslist} ";
                opserr<<"-roofmass {masslist}";
                opserr<<"-colvel {vellist} -beamvel {vellist}>\n";
                return TCL_ERROR;
            }
            int regTag, ndf;
            double x,y,hcol, hbeam;
            if(Tcl_GetInt(interp, argv[3], &regTag) != TCL_OK) {
                opserr<<"WARNING: invalid region tag "<<argv[3]<<" -- PFEM2D region frame\n";
                return TCL_ERROR; 
            }
            if(Tcl_GetDouble(interp, argv[4], &x) != TCL_OK) {
                opserr<<"WARNING: invalid x "<<argv[4]<<" -- PFEM2D region frame\n";
                return TCL_ERROR; 
            }
            if(Tcl_GetDouble(interp, argv[5], &y) != TCL_OK) {
                opserr<<"WARNING: invalid y "<<argv[5]<<" -- PFEM2D region frame\n";
                return TCL_ERROR; 
            }
            if(Tcl_GetDouble(interp, argv[6], &hcol) != TCL_OK) {
                opserr<<"WARNING: invalid ncol "<<argv[6]<<" -- PFEM2D region frame\n";
                return TCL_ERROR; 
            }
            if(Tcl_GetDouble(interp, argv[7], &hbeam) != TCL_OK) {
                opserr<<"WARNING: invalid hbeam "<<argv[7]<<" -- PFEM2D region frame\n";
                return TCL_ERROR; 
            }
            if(Tcl_GetInt(interp, argv[8], &ndf) != TCL_OK) {
                opserr<<"WARNING: invalid nbeam "<<argv[8]<<" -- PFEM2D region frame\n";
                return TCL_ERROR; 
            }

            int loc = 9;
            Vector colmass, beammass, colvel, beamvel, roofmass;
            Vector fix, spans, heights;
            while(loc < argc) {

                Vector* vecPtr = 0;

                if(loc == 9) {
                    vecPtr = &spans;

                } else if(loc == 10) {
                    vecPtr = &heights;

                } else if(strcmp(argv[loc], "-colmass") == 0) {
                    vecPtr = &colmass;
                    loc++;

                } else if(strcmp(argv[loc], "-beammass") == 0) {
                    vecPtr = &beammass;
                    loc++;

                } else if(strcmp(argv[loc], "-fix") == 0) {
                    vecPtr = &fix;
                    loc++;
  
                } else if(strcmp(argv[loc], "-colvel") == 0) {
                    vecPtr = &colvel;
                    loc++;

                } else if(strcmp(argv[loc], "-beamvel") == 0) {
                    vecPtr = &beamvel;
                    loc++;
                } else if(strcmp(argv[loc], "-roofmass") == 0) {
                    vecPtr = &roofmass;
                    loc++;
                }

                if(vecPtr != 0) {
                    int num = 0;
                    const char** argvPtr = 0;
                    if(Tcl_SplitList(interp, argv[loc], &num, &argvPtr)!=TCL_OK) {
                        opserr<<"WARNING: failed to read the list "<<argv[loc];
                        opserr<<" -- PFEM2D region\n";
                        return TCL_ERROR;
                    }
                    if(num > 0) vecPtr->resize(num);
                    for(int i=0; i<num; i++) {
                        double x;
                        if(Tcl_GetDouble(interp, argvPtr[i], &x) != TCL_OK) {
                            opserr<<"WARNING: invalid input "<<argvPtr[i];
                            opserr<<" -- PFEM2D region\n";
                            return TCL_ERROR; 
                        }
                        (*vecPtr)(i) = x;
                    }

                    if(argvPtr != 0) Tcl_Free((char *) argvPtr);
                }

                loc++;
            }

            ID nodelist;
            int startnode = theMesher2D.findNodeTag(theDomain);
            int res = theMesher2D.discretize(startnode,x,y,hcol,hbeam,ndf,
                                             spans,heights,fix,colvel,
                                             beamvel,colmass,beammass,
                                             roofmass,theDomain,nodelist);
            if(res < 0) {
                opserr << "WARNING: failed to discretize\n";
                return TCL_ERROR;
            }

            // add region
            ID regionnodes(2);
            regionnodes(0) = startnode;
            regionnodes(1) = nodelist(nodelist.Size()-1);
            bool series = false;
            int action = 1; // append
            theMesher2D.setNodes(regionnodes, regTag, series, action, 
                                 theDomain);

        } else if(strcmp(argv[2],"mesh")==0 ||
                  strcmp(argv[2],"Mesh")==0) {
            if(argc < 5) {
                opserr << "WARNING: wrong num of args -- ";
                opserr << "PFEM2D region mesh alpha eleRegTag";
                opserr << "-nodeRegions {regTag1 ...}  ";
                opserr << "-addNodeRegions {regTag2 ...} ";
                opserr << "-PFEMElement2D {rho mu b1 b2 <thk kappa>} ";
                opserr << "-PFEMElement2DCompressible {rho mu b1 b2 <thk kappa>} ";
                opserr << "-PFEMElement2DBubble {rho mu b1 b2 <thk kappa>} ";
                opserr << "-Tri31 {thk type matTag <pressure rho b1 b2>}\n ";
                return TCL_ERROR;
            }
            double alpha;
            if(Tcl_GetDouble(interp, argv[3], &alpha) != TCL_OK) {
                opserr<<"WARNING: invalid alpha "<<argv[3]<<" -- PFEM2D region mesh\n";
                return TCL_ERROR; 
            }
            int eleRegTag;
            if(Tcl_GetInt(interp, argv[4], &eleRegTag) != TCL_OK) {
                opserr<<"WARNING: invalid element region tag "<<argv[4]<<" -- PFEM2D region mesh\n";
                return TCL_ERROR; 
            }

            int loc = 5;
            Vector params;
            int eletype = 0;
            std::string type;
            ID nodes, addnodes;
            while(loc < argc) {

                Vector* vecPtr = 0;
                ID* idPtr = 0;

                if(strcmp(argv[loc], "-nodeRegions") == 0) {
                    idPtr = &nodes;
                    loc++;

                } else if(strcmp(argv[loc], "-addNodeRegions") == 0) {
                    idPtr = &addnodes;
                    loc++;

                } else if(strcmp(argv[loc], "-PFEMElement2D") == 0) {
                    vecPtr = &params;
                    eletype = 1;
                    loc++;

                } else if(strcmp(argv[loc], "-Tri31") == 0 || 
                          strcmp(argv[loc], "-tri31") == 0) {
                    vecPtr = &params;
                    eletype = 2;
                    loc++;
                } else if(strcmp(argv[loc], 
                                 "-PFEMElement2DCompressible") == 0) {
                    vecPtr = &params;
                    eletype = 3;
                    loc++;
                } else if(strcmp(argv[loc], "-PFEMElement2DBubble") == 0) {
                    vecPtr = &params;
                    eletype = 4;
                    loc++;
                }

                if(vecPtr!=0 || idPtr!=0) {
                    int num = 0;
                    const char** argvPtr = 0;
                    if(Tcl_SplitList(interp, argv[loc], &num, &argvPtr)!=TCL_OK) {
                        opserr<<"WARNING: failed to read the list "<<argv[loc];
                        opserr<<" -- PFEM2D region mesh\n";
                        return TCL_ERROR;
                    }
                    if(num > 0) {
                        if(vecPtr!=0) vecPtr->resize(num);
                        if(idPtr!=0) idPtr->resize(num);
                    }
                    if(vecPtr != 0) {
                        for(int i=0; i<num; i++) {
                            if(eletype==2 && i==2) {
                                type = argvPtr[i];
                                continue;
                            }
                            double x;
                            if(Tcl_GetDouble(interp, argvPtr[i], &x) != TCL_OK) {
                                opserr<<"WARNING: invalid input "<<argvPtr[i];
                                opserr<<" -- PFEM2D region mesh\n";
                                return TCL_ERROR; 
                            }
                            (*vecPtr)(i) = x;
                        }
                    } else if(idPtr != 0) {
                        for(int i=0; i<num; i++) {
                            int x;
                            if(Tcl_GetInt(interp, argvPtr[i], &x) != TCL_OK) {
                                opserr<<"WARNING: invalid input "<<argvPtr[i];
                                opserr<<" -- PFEM2D region mesh\n";
                                return TCL_ERROR; 
                            }
                            (*idPtr)(i) = x;
                        }
                    }

                    if(argvPtr != 0) Tcl_Free((char *) argvPtr);
                }

                loc++;
            }

            // triangulation
            int res = -1;
            ID eles;
            //opserr<<"params = "<<params;
            if(eletype>=1 && eletype!=2 && eletype<=7) {
                double rho,mu,b1,b2,thk=1.0,kappa=2.2e9;
                if(params.Size() > 3) {
                    rho = params(0);
                    mu = params(1);
                    b1 = params(2);
                    b2 = params(3);
                } else {
                    opserr<<"WARNING: not sufficient args ";
                    opserr<<" -- PFEM2D region mesh\n";
                    return TCL_ERROR;
                }
                if(params.Size() > 4) {
                    thk = params(4);
                }
                if(params.Size() > 5) {
                    kappa = params(5);
                }

                // remove all elements in eleReg
                theMesher2D.removeElements(eleRegTag,theDomain);

                // triangulation
                int startele = theMesher2D.findEleTag(theDomain);
                int endele = startele;
            
                res = theMesher2D.doTriangulation(startele,alpha,nodes,
                                                  addnodes,theDomain,
                                                  rho,mu,b1,b2,thk,kappa,
                                                  eletype,endele);
            
                if(res < 0) {
                    opserr<<"WARNING: failed to mesh -- ";
                    opserr<<" -- PFEM2D region mesh\n";
                    return TCL_ERROR; 
                }

                // add to element region
                ID regioneles(2);
                regioneles(0) = startele;
                regioneles(1) = endele;
                bool series = false;
                int action = 1; // append
                theMesher2D.setElements(regioneles,eleRegTag,series,action,
                                        theDomain);

            } else if(eletype == 2) {
                double thk=0, p=0,rho=0,b1=0,b2=0;
                int matTag=0;
                if(params.Size() > 2) {
                    thk = params(0);
                    matTag = (int)params(2);
                } else {
                    opserr<<"WARNING: not sufficient args ";
                    opserr<<" -- PFEM2D region mesh\n";
                    return TCL_ERROR;
                }
                if(params.Size() > 3) {
                    p = params(3);
                }
                if(params.Size() > 4) {
                    rho = params(4);
                }
                if(params.Size() > 5) {
                    b1 = params(5);
                }
                if(params.Size() > 6) {
                    b2 = params(6);
                }

                // remove all elements in eleReg
                theMesher2D.removeElements(eleRegTag,theDomain);
            
                // triangulation
                int startele = theMesher2D.findEleTag(theDomain);
                int endele = startele;
                res = theMesher2D.doTriangulation(startele,alpha,nodes,addnodes,
                                                  theDomain,
                                                  thk,type.c_str(),matTag,p,rho,
                                                  b1,b2,endele);

                if(res < 0) {
                    opserr<<"WARNING: failed to do triangulation -- ";
                    opserr<<" -- PFEM2D region mesh\n";
                    return TCL_ERROR; 
                }

                // add to element region
                ID regioneles(2);
                regioneles(0) = startele;
                regioneles(1) = endele;
                bool series = false;
                int action = 1; // append
                theMesher2D.setElements(regioneles,eleRegTag,series,action,
                                        theDomain);

            } else {
                opserr<<"WARNING: no element type is set -- ";
                opserr<<" -- PFEM2D region mesh\n";
                return TCL_ERROR; 
            }

        } else if(strcmp(argv[2],"node")==0||strcmp(argv[2],"Node")==0) {
            if(argc < 4) {
                opserr<<"WARNING: wrong no of args -- ";
                opserr<<"PFEM2D region node op \n";
                return TCL_ERROR;
            }
            int action=0,regTag=0;
            bool series=false;
            ID nodelist;
            if(strcmp(argv[3],"add")==0||strcmp(argv[3],"Add")==0) {
                action = 1;
                series = true;
            } else if(strcmp(argv[3],"remove")==0||
                      strcmp(argv[3],"Remove")==0) {
                action = 2;
                series = true;
            } else if(strcmp(argv[3],"set")==0||
                      strcmp(argv[3],"Set")==0) {
                action = 0;
                series = true;
            } else if(strcmp(argv[3],"addRange")==0) {
                action = 1;
                series = false;
            } else if(strcmp(argv[3],"removeRange")==0) {
                action = 2;
                series = false;
            } else if(strcmp(argv[3],"setRange")==0) {
                action = 0;
                series = false;
            } else if(strcmp(argv[3],"get")==0||
                      strcmp(argv[3],"Get")==0) {

                if(Tcl_GetInt(interp, argv[4], &regTag) != TCL_OK) {
                    opserr<<"WARNING: invalid region tag "<<argv[4];
                    opserr<<" -- PFEM2D region node "<<argv[3]<<"\n";
                    return TCL_ERROR;
                }

                MeshRegion* region = theDomain->getRegion(regTag);
                if(region == 0) return TCL_OK;
                const ID& regNodes = region->getNodes();
                char buffer[20];
                for(int i=0; i<regNodes.Size(); i++) {
                    sprintf(buffer, "%d ", regNodes(i));
                    Tcl_AppendResult(interp,buffer,NULL);
                }

                return TCL_OK;
                
            } else {
                opserr<<"WARNING: unknown operation on region -- ";
                opserr<<"PFEM2D region node "<<argv[3]<<"\n";
                return TCL_ERROR;
            }

            if(argc < 6) {
                opserr<<"WARNING: wrong no of args -- ";
                opserr<<"PFEM2D region node op regTag nodelist \n";
                return TCL_ERROR;
            }

            if(Tcl_GetInt(interp, argv[4], &regTag) != TCL_OK) {
                opserr<<"WARNING: invalid region tag "<<argv[4];
                opserr<<" -- PFEM2D region node "<<argv[3]<<"\n";
                return TCL_ERROR;
            }

            int num = 0;
            const char** argvPtr = 0;
            if(Tcl_SplitList(interp, argv[5], &num, &argvPtr)!=TCL_OK) {
                opserr<<"WARNING: failed to read the list "<<argv[5];
                opserr<<" -- PFEM2D region node "<<argv[3]<<"\n";
                return TCL_ERROR;
            }
            if(num > 0) nodelist.resize(num);
            for(int i=0; i<num; i++) {
                int tag;
                if(Tcl_GetInt(interp, argvPtr[i], &tag) != TCL_OK) {
                    opserr<<"WARNING: invalid input "<<argvPtr[i];
                    opserr<<" -- PFEM2D region node "<<argv[3]<<"\n";
                    return TCL_ERROR; 
                }
                nodelist(i) = tag;
            }

            if(argvPtr != 0) Tcl_Free((char *) argvPtr);

            theMesher2D.setNodes(nodelist,regTag,series,action,theDomain);
            
        } else if(strcmp(argv[2],"ele")==0||strcmp(argv[2],"Ele")==0) {

            if(argc < 4) {
                opserr<<"WARNING: wrong no of args -- ";
                opserr<<"PFEM2D region ele op \n";
                return TCL_ERROR;
            }
            int action=0,regTag=0;
            bool series=false;
            ID elelist;
            if(strcmp(argv[3],"add")==0||strcmp(argv[3],"Add")==0) {
                action = 1;
                series = true;
            } else if(strcmp(argv[3],"remove")==0||
                      strcmp(argv[3],"Remove")==0) {
                action = 2;
                series = true;
            } else if(strcmp(argv[3],"set")==0||
                      strcmp(argv[3],"Set")==0) {
                action = 0;
                series = true;
            } else if(strcmp(argv[3],"addRange")==0) {
                action = 1;
                series = false;
            } else if(strcmp(argv[3],"removeRange")==0) {
                action = 2;
                series = false;
            } else if(strcmp(argv[3],"setRange")==0) {
                action = 0;
                series = false;
            } else if(strcmp(argv[3],"get")==0||
                      strcmp(argv[3],"Get")==0) {

                if(Tcl_GetInt(interp, argv[4], &regTag) != TCL_OK) {
                    opserr<<"WARNING: invalid region tag "<<argv[4];
                    opserr<<" -- PFEM2D region ele "<<argv[3]<<"\n";
                    return TCL_ERROR;
                }

                MeshRegion* region = theDomain->getRegion(regTag);
                if(region == 0) return TCL_OK;
                const ID& regEles = region->getElements();
                char buffer[20];
                for(int i=0; i<regEles.Size(); i++) {
                    sprintf(buffer, "%d ", regEles(i));
                    Tcl_AppendResult(interp,buffer,NULL);
                }

                return TCL_OK;
                
            } else {
                opserr<<"WARNING: unknown operation on region -- ";
                opserr<<"PFEM2D region ele "<<argv[3]<<"\n";
                return TCL_ERROR;
            }

            if(argc < 6) {
                opserr<<"WARNING: wrong no of args -- ";
                opserr<<"PFEM2D region ele op regTag elelist \n";
                return TCL_ERROR;
            }

            if(Tcl_GetInt(interp, argv[4], &regTag) != TCL_OK) {
                opserr<<"WARNING: invalid region tag "<<argv[4];
                opserr<<" -- PFEM2D region ele "<<argv[3]<<"\n";
                return TCL_ERROR;
            }

            int num = 0;
            const char** argvPtr = 0;
            if(Tcl_SplitList(interp, argv[5], &num, &argvPtr)!=TCL_OK) {
                opserr<<"WARNING: failed to read the list "<<argv[5];
                opserr<<" -- PFEM2D region ele "<<argv[3]<<"\n";
                return TCL_ERROR;
            }
            if(num > 0) elelist.resize(num);
            for(int i=0; i<num; i++) {
                int tag;
                if(Tcl_GetInt(interp, argvPtr[i], &tag) != TCL_OK) {
                    opserr<<"WARNING: invalid input "<<argvPtr[i];
                    opserr<<" -- PFEM2D region node "<<argv[3]<<"\n";
                    return TCL_ERROR; 
                }
                elelist(i) = tag;
            }

            if(argvPtr != 0) Tcl_Free((char *) argvPtr);

            theMesher2D.setElements(elelist,regTag,series,action,theDomain);

        } else if(strcmp(argv[2],"save")==0) {

            if(argc < 6) {
                opserr<<"WARNING: wrong no of args -- ";
                opserr<<"PFEM2D region save nodeRegions eleRegions ";
                opserr<<"filename <maxelenodes> \n";
                return TCL_ERROR;
            }
            ID nodeRegs,eleRegs;
            int num = 0;
            const char** argvPtr = 0;
            if(Tcl_SplitList(interp, argv[3], &num, &argvPtr)!=TCL_OK) {
                opserr<<"WARNING: failed to read the list "<<argv[3];
                opserr<<" -- PFEM2D region save \n";
                return TCL_ERROR;
            }
            if(num > 0) nodeRegs.resize(num);
            for(int i=0; i<num; i++) {
                int tag;
                if(Tcl_GetInt(interp, argvPtr[i], &tag) != TCL_OK) {
                    opserr<<"WARNING: invalid input "<<argvPtr[i];
                    opserr<<" -- PFEM2D region save \n";
                    return TCL_ERROR; 
                }
                nodeRegs(i) = tag;
            }
            if(argvPtr != 0) Tcl_Free((char *) argvPtr);
            if(Tcl_SplitList(interp, argv[4], &num, &argvPtr)!=TCL_OK) {
                opserr<<"WARNING: failed to read the list "<<argv[4];
                opserr<<" -- PFEM2D region save \n";
                return TCL_ERROR;
            }
            if(num > 0) eleRegs.resize(num);
            for(int i=0; i<num; i++) {
                int tag;
                if(Tcl_GetInt(interp, argvPtr[i], &tag) != TCL_OK) {
                    opserr<<"WARNING: invalid input "<<argvPtr[i];
                    opserr<<" -- PFEM2D region save \n";
                    return TCL_ERROR; 
                }
                eleRegs(i) = tag;
            }
            if(argvPtr != 0) Tcl_Free((char *) argvPtr);

            int max = 3;
            if(argc > 6) {
                if(Tcl_GetInt(interp, argv[6], &max) != TCL_OK) {
                    opserr<<"WARNING: invalid integer "<<argv[6];
                    opserr<<" -- PFEM2D region save \n";
                    return TCL_ERROR;
                }
            }

            theMesher2D.save(argv[5],nodeRegs,eleRegs,theDomain,max);

        } else {
            opserr<<"WARNING: unknown subcommand -- PFEM2D region ";
            opserr<<argv[2]<<"\n";
            return TCL_ERROR;            
        }

    } else {
        opserr<<"WARNING: unknown command -- PFEM2D "<<argv[1]<<"\n";
        return TCL_ERROR;
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
        if(argc < 3) {
            opserr << "WARNING: wrong num of args -- ";
            opserr << "PFEM3D discretize type\n";
            return TCL_ERROR;
        }
        if(strcmp(argv[2], "PLC") == 0) {

            // discretize
            if(argc < 6) {
                opserr<<"WARNING: wrong num of args -- ";
                opserr<<"PFEM3D discretize PLC startnode maxvol ndf ";
                opserr<<"<-points (x1 y1 ...) ";
                opserr<<"-facet (p nd1 nd2 ... h x1 y1 z1 ...) ";
                opserr<<"-fix (constrValues) -mass (massValues) -holes (x1 y1 ...)>\n";
                return TCL_ERROR;
            }
            int startnode;
            if(Tcl_GetInt(interp, argv[3], &startnode) != TCL_OK) {
                opserr<<"WARNING: invalid startnode "<<argv[3]<<" -- PFEM3D discretize\n";
                return TCL_ERROR; 
            }
            double maxvol;
            if(Tcl_GetDouble(interp, argv[4], &maxvol) != TCL_OK) {
                opserr<<"WARNING: invalid maxvol "<<argv[4]<<" -- PFEM3D discretize\n";
                return TCL_ERROR; 
            }
            if(maxvol <= 0) {
                opserr<<"WARNING: nonpositive area ";
                opserr<<" -- PFEM3D discretize\n";
                return TCL_ERROR;
            }
            int ndf;
            if(Tcl_GetInt(interp, argv[5], &ndf) != TCL_OK) {
                opserr<<"WARNING: invalid ndf "<<argv[5]<<" -- PFEM3D discretize\n";
                return TCL_ERROR; 
            }
            if(ndf <= 0) {
                opserr<<"WARNING: nonpositive ndf ";
                opserr<<" -- PFEM3D discretize\n";
                return TCL_ERROR;
            }

            PFEMMesher3D::PointVec points, holes;
            std::vector<double> mass;
            std::vector<int> fix;
            PFEMMesher3D::FacetVec facets;

            // scan for switches
            std::vector<int> states;
            std::vector<int> locs;
            for(int loc=6; loc<argc; loc++) {
                if(strcmp(argv[loc], "-points") == 0) {
                    states.push_back(1);
                    locs.push_back(loc);
                } else if(strcmp(argv[loc], "-facet") == 0) {
                    states.push_back(2);
                    locs.push_back(loc);
                } else if(strcmp(argv[loc], "-mass") == 0) {
                    states.push_back(3);
                    locs.push_back(loc);
                } else if(strcmp(argv[loc], "-fix") == 0) {
                    states.push_back(4);
                    locs.push_back(loc);
                } else if(strcmp(argv[loc], "-holes") == 0) {
                    states.push_back(5);
                    locs.push_back(loc);
                }
            }
            locs.push_back(argc);
            for(int i=0; i<(int)states.size(); i++) {
                int state = states[i];
                int l0 = locs[i]+1;
                int l1 = locs[i+1];
                if(state == 1) {         // points
                    double coord[3];
                    int count = 0;
                    for(int loc=l0; loc<l1; loc++) {
                        if(Tcl_GetDouble(interp, argv[loc], &coord[count++]) != TCL_OK) {
                            opserr<<"WARNING: invalid coordinate "<<argv[loc];
                            opserr<<" -- PFEM3D discretize\n";
                            return TCL_ERROR; 
                        }
                        if(count == 3) {
                            points.push_back(PFEMMesher3D::Point(coord[0],coord[1],coord[2]));
                            count = 0;
                        }
                    }
                } else if(state == 2) {  // facet

                    // scan for 'h' and 'p'
                    std::vector<int> fstates;
                    std::vector<int> flocs;
                    for(int loc=l0; loc<l1; loc++) {
                        if(strcmp(argv[loc], "h") == 0) {
                            fstates.push_back(1);
                            flocs.push_back(loc);
                        } else if(strcmp(argv[loc], "p") == 0) {
                            fstates.push_back(2);
                            flocs.push_back(loc);
                        }
                    }
                    flocs.push_back(l1);

                    if(!fstates.empty()) {
                        facets.push_back(PFEMMesher3D::Facet());
                    }

                    for(int j=0; j<(int)fstates.size(); j++)
                    {
                        int fstate = fstates[j];
                        int fl0 = flocs[j]+1;
                        int fl1 = flocs[j+1];
                        if(fstate == 1) {  // facet hole
                            double coord[3];
                            int count = 0;
                            for(int floc=fl0; floc<fl1; floc++) {
                                if(Tcl_GetDouble(interp, argv[floc], &coord[count++]) != TCL_OK) {
                                    opserr<<"WARNING: invalid coordinate "<<argv[floc];
                                    opserr<<" -- PFEM3D discretize\n";
                                    return TCL_ERROR; 
                                }
                                if(count == 3) {
                                    facets.back().addhole(coord[0],coord[1],coord[2]);
                                    count = 0;
                                }
                            }
                        } else if(fstate == 2) {  // facet polygon
                            PFEMMesher3D::Polygon polygon;
                            for(int floc=fl0; floc<fl1; floc++) {
                                int tag;
                                if(Tcl_GetInt(interp, argv[floc], &tag) != TCL_OK) {
                                    opserr<<"WARNING: invalid point "<<argv[floc];
                                    opserr<<" -- PFEM3D discretize\n";
                                    return TCL_ERROR; 
                                }
                                polygon.push_back(tag);
                            }
                            facets.back().addpolygon(polygon);
                        }
                    }

                } else if(state == 3) {  // mass
                    for(int loc=l0; loc<l1; loc++) {
                        double m;
                        if(Tcl_GetDouble(interp, argv[loc], &m) != TCL_OK) {
                            opserr<<"WARNING: invalid mass "<<argv[loc];
                            opserr<<" -- PFEM3D discretize\n";
                            return TCL_ERROR; 
                        }
                        mass.push_back(m);
                    }
                } else if(state == 4) {  // fix
                    for(int loc=l0; loc<l1; loc++) {
                        int fixity;
                        if(Tcl_GetInt(interp, argv[loc], &fixity) != TCL_OK) {
                            opserr<<"WARNNG: invalid fixity "<<argv[loc];
                            opserr<<" -- PFEM3D discretize\n";
                            return TCL_ERROR;
                        }
                        fix.push_back(fixity);
                    }
                } else if(state == 5) {  // holes
                    double coord[3];
                    int count = 0;
                    for(int loc=l0; loc<l1; loc++) {
                        if(Tcl_GetDouble(interp, argv[loc], &coord[count++]) != TCL_OK) {
                            opserr<<"WARNING: invalid coordinate "<<argv[loc];
                            opserr<<" -- PFEM3D discretize\n";
                            return TCL_ERROR; 
                        }
                        if(count == 3) {
                            holes.push_back(PFEMMesher3D::Point(coord[0],coord[1],coord[2]));
                            count = 0;
                        }
                    }
                }
            }

            int endnode = startnode;
            int res = theMesher3D.discretize(startnode,points,facets,holes,maxvol,ndf,
                                             fix,mass,theDomain,endnode);

            if(res < 0) {
                opserr << "WARNING: failed to discretize\n";
                return TCL_ERROR;
            }
            
            Tcl_SetObjResult(interp, Tcl_NewIntObj(endnode));
            
        }

    } else if(strcmp(argv[1], "doTriangulation") == 0) {

        // triangulation
        if(argc < 4) {
            opserr << "WARNING: wrong num of args -- ";
            opserr << "PFEM3D doTriangulation alpha -volthresh volthres ";
            opserr << "-nodes (start1 end1 start2 end2 ...)  ";
            opserr << "<-PFEMElement3D starteletag rho mu b1 b2 b3>\n";
            return TCL_ERROR;
        }
        double alpha;
        if(Tcl_GetDouble(interp, argv[2], &alpha) != TCL_OK) {
            opserr<<"WARNING: invalid alpha "<<argv[2]<<" -- PFEM3D doTriangulation\n";
            return TCL_ERROR; 
        }

        // scan for switches
        std::vector<int> states;
        std::vector<int> locs;
        for(int loc=3; loc<argc; loc++) {
            if(strcmp(argv[loc], "-nodes") == 0) {
                states.push_back(1);
                locs.push_back(loc);
            } else if(strcmp(argv[loc], "-addnodes") == 0) {
                states.push_back(2);
                locs.push_back(loc);
            } else if(strcmp(argv[loc], "-PFEMElement3D") == 0) {
                states.push_back(3);
                locs.push_back(loc);
            } else if(strcmp(argv[loc], "-volthresh") == 0) {
                states.push_back(4);
                locs.push_back(loc);
            }
        }
        locs.push_back(argc);

        std::vector<int> nodes, addnodes;
        std::vector<double> pfemparams(5);
        int starteletag;
        bool pfem = false;
        double volthresh = 1e-6;
        for(int i=0; i<(int)states.size(); i++) {
            int state = states[i];
            int l0 = locs[i]+1;
            int l1 = locs[i+1];
            if(state == 1) {          // nodes
                if(l1-l0<2) {
                    opserr<<"WARNING: wrong no of arguments -- PFEM3D doTriangulation\n";
                    return TCL_ERROR;
                }
                for(int loc=l0; loc<l1; loc++) {
                    int tag;
                    if(Tcl_GetInt(interp, argv[loc], &tag) != TCL_OK) {
                        opserr<<"WARNING: invalid node tag "<<argv[loc];
                        opserr<<" -- PFEM3D doTriangulation\n";
                        return TCL_ERROR; 
                    }
                    nodes.push_back(tag);
                }
            } else if(state == 2) {   // addnodes
                if(l1-l0<2) {
                    opserr<<"WARNING: wrong no of arguments -- PFEM3D doTriangulation\n";
                    return TCL_ERROR;
                }
                for(int loc=l0; loc<l1; loc++) {
                    int tag;
                    if(Tcl_GetInt(interp, argv[loc], &tag) != TCL_OK) {
                        opserr<<"WARNING: invalid node tag "<<argv[loc];
                        opserr<<" -- PFEM3D doTriangulation\n";
                        return TCL_ERROR; 
                    }
                    addnodes.push_back(tag);
                }
            } else if(state == 3) {   // PFEMElement3D
                if(l1-l0<6) {
                    opserr<<"WARNING: wrong no of arguments -- PFEM3D doTriangulation\n";
                    return TCL_ERROR;
                }
                pfem = true;
                if(Tcl_GetInt(interp, argv[l0], &starteletag) != TCL_OK) {
                    opserr<<"WARNING: invalid element tag "<<argv[l0];
                    opserr<<" -- PFEM3D doTriangulation\n";
                    return TCL_ERROR; 
                }
                for(int loc=l0+1; loc<l0+6; loc++) {
                    if(Tcl_GetDouble(interp, argv[loc], &pfemparams[loc-l0-1]) != TCL_OK) {
                        opserr<<"WARNING: invalid parameter "<<argv[loc];
                        opserr<<" -- PFEM3D doTriangulation\n";
                        return TCL_ERROR; 
                    }
                }
            } else if(state == 4) {
                if(l1-l0 < 1) {
                    opserr<<"WARNING: wrong no of arguments -- PFEM3D doTriangulation\n";
                    return TCL_ERROR;
                }
                if(Tcl_GetDouble(interp, argv[l0], &volthresh) != TCL_OK) {
                    opserr<<"WARNING: invalid volthresh "<<argv[l0];
                    opserr<<" -- PFEM3D doTriangulation\n";
                    return TCL_ERROR; 
                }
            }
        }

        // triangulation
        int res;
        std::vector<int> eles;
        if(pfem) {
            res = theMesher3D.doTriangulation(starteletag,nodes,
                                              alpha,volthresh,addnodes,theDomain,
                                              pfemparams[0],pfemparams[1],
                                              pfemparams[2],pfemparams[3],pfemparams[4]);
        } else {
            res = theMesher3D.doTriangulation(nodes,alpha,volthresh,
                                              addnodes,theDomain,eles);
            if(res < 0) {
                opserr<<"WARNING: failed to do triangulation -- ";
                opserr<<" -- PFEM3D doTriangulation\n";
                return TCL_ERROR; 
            }
        }

        

        if(pfem) {
            Tcl_SetObjResult(interp, Tcl_NewIntObj(res));
        } else {
            for(int i=0; i<(int)eles.size(); i++) {
                char buffer[100];
                sprintf(buffer, "%d ", eles[i]);
                Tcl_AppendResult(interp, buffer, NULL);
            }
        }

    } else if(strcmp(argv[1], "save") == 0) {

        if(argc < 4) {
            opserr << "WARNING: wrong num of args -- ";
            opserr << "PFEM3D save filename step <-structure (start1 end1 ...)> ";
            return TCL_ERROR;            
        }

        int step;
        if(Tcl_GetInt(interp, argv[3], &step) != TCL_OK) {
            opserr<<"WARNING: invalid step "<<argv[3]<<" -- PFEM3D save\n";
            return TCL_ERROR; 
        }
        
        if(argc>6 && strcmp(argv[4], "-structure")==0) {
            int numnodes = argc-5;
            ID snodes(numnodes);
            for(int i=0; i<numnodes; i++) {
                if(Tcl_GetInt(interp, argv[5+i], &snodes(i)) != TCL_OK) {
                    opserr<<"WARNING: invalid strcutural node "<<argv[4+i]<<" -- PFEM3D save\n";
                    return TCL_ERROR; 
                }
            }
            theMesher3D.save(argv[2], snodes, step, theDomain);
        
        } else {
            ID snodes;
            theMesher3D.save(argv[2], snodes, step, theDomain);
        }

        // if(argc < 4) {
        //     opserr << "WARNING: wrong num of args -- ";
        //     opserr << "PFEM3D save filename step";
        //     return TCL_ERROR;            
        // }

        // int step = 0;
        // if(Tcl_GetInt(interp, argv[3], &step) != TCL_OK) {
        //     opserr<<"WARNING: invalid step "<<argv[3]<<" -- PFEM3D save\n";
        //     return TCL_ERROR; 
        // }

        // theMesher3D.save(argv[2], step, theDomain);

    } else if(strcmp(argv[1], "setBoundary") == 0) {
        
        if(argc < 8) {
            opserr << "WARNING: wrong num of args -- ";
            opserr << "PFEM3D setBoundary x1 y1 z1 x2 y2 z2\n";
            return TCL_ERROR;            
        }
        double bound[6];
        for(int i=0; i<6; i++) {
            if(Tcl_GetDouble(interp, argv[i+2], &bound[i]) != TCL_OK) {
                opserr<<"WARNING: invalid coordinate "<<argv[i+2];
                opserr<<" -- PFEM3D setBoundary\n";
                return TCL_ERROR; 
            }
        }
        
        theMesher3D.setBoundary(bound[0],bound[1],bound[2],bound[3],bound[4],bound[5]);

    } else if(strcmp(argv[1], "removeOutBoundNodes") == 0) {

        if(argc < 2) {
            opserr << "WARNING: wrong num of args -- ";
            opserr << "PFEM3D removeOutBoundNodes <-nodes startnode end ...>\n";
            return TCL_ERROR;            
        }

        if(argc>4 && strcmp(argv[3], "-nodes")==0) {
            int numnodes = argc-3;
            ID nodes(numnodes);
            for(int i=0; i<numnodes; i++) {
                if(Tcl_GetInt(interp, argv[3+i], &nodes(i)) != TCL_OK) {
                    opserr<<"WARNING: invalid node "<<argv[3+i]<<" -- PFEM3D removeOutBoundNodes\n";
                    return TCL_ERROR; 
                }
            }
            theMesher3D.removeOutBoundNodes(nodes, theDomain);
        
        }
    }

    return TCL_OK;
}

