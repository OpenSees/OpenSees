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
                opserr << "-fix (constrValues) -mass (massValues) -holes (x1 y1 ...) >\n";
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
            Vector points(argc-loc), mass(argc-loc), holes(argc-loc), vel(argc-loc);
            ID segments, fix;
            int state=-1, numpoints=0, numsegments=0, nummass=0, numfix=0, numvel=0,numholes=0;
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
                } else if(strcmp(argv[loc], "-vel") == 0) {
                    state = 6;
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
                    }  else if(state==6) {
                        double v;
                        if(Tcl_GetDouble(interp, argv[loc], &v) != TCL_OK) {
                            opserr<<"WARNING: invalid vel "<<argv[loc];
                            opserr<<" -- PFEM2D discretize\n";
                            return TCL_ERROR; 
                        }
                        vel(numvel++) = v;
                    }
                }
            
                loc++;
            }
            points.resize(numpoints);
            mass.resize(nummass);
            holes.resize(numholes);
            vel.resize(numvel);
            int endnodetag = startnodetag;
            int res = theMesher2D.discretize(startnodetag,points,segments,holes,maxarea,ndf,
                                             fix,vel,mass,theDomain, endnodetag);

            if(res < 0) {
                opserr << "WARNING: failed to discretize\n";
                return TCL_ERROR;
            }
            
            Tcl_SetObjResult(interp, Tcl_NewIntObj(endnodetag));

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
            Vector mass(argc-loc), vel(argc-loc);
            ID fix;
            int state=-1, numfix=0, nummass=0, numvel=0;
            while(loc < argc) {
                if(strcmp(argv[loc], "-fix") == 0) {
                    state = 1;
                } else if(strcmp(argv[loc], "-mass") == 0) {
                    state = 2;
                } else if(strcmp(argv[loc], "-vel") == 0) {
                    state = 3;
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
                    } else if(state==3) {
                        double v;
                        if(Tcl_GetDouble(interp, argv[loc], &v) != TCL_OK) {
                            opserr<<"WARNING: invalid vel "<<argv[loc];
                            opserr<<" -- PFEM2D discretizeRectangle\n";
                            return TCL_ERROR; 
                        }
                        vel(numvel++) = v;
                    }
                }
            
                loc++;
            }
            mass.resize(nummass);
            vel.resize(numvel);
            int endnodetag = startnodetag;
            int res = theMesher2D.discretize(startnodetag,x1,y1,hx,hy,angle,nx,ny,ndf,fix,vel,mass,
                                             theDomain,endnodetag);

            if(res < 0) {
                opserr << "WARNING: failed to discretize\n";
                return TCL_ERROR;
            }
 
            Tcl_SetObjResult(interp, Tcl_NewIntObj(endnodetag));

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
            Vector mass(argc-loc), vel(argc-loc);
            ID fix;
            int state=-1, numfix=0, nummass=0, numvel=0;
            while(loc < argc) {
                if(strcmp(argv[loc], "-fix") == 0) {
                    state = 1;
                } else if(strcmp(argv[loc], "-mass") == 0) {
                    state = 2;
                } else if(strcmp(argv[loc], "-vel") == 0) {
                    state = 3;
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
                    } else if(state==3) {
                        double v;
                        if(Tcl_GetDouble(interp, argv[loc], &v) != TCL_OK) {
                            opserr<<"WARNING: invalid vel "<<argv[loc];
                            opserr<<" -- PFEM2D discretizeRectangle\n";
                            return TCL_ERROR; 
                        }
                        vel(numvel++) = v;
                    }
                }
            
                loc++;
            }
            mass.resize(nummass);
            vel.resize(numvel);
            int endnodetag = startnodetag;
            int res = theMesher2D.discretize(startnodetag, x1,y1,h,angle,num,ndf,fix,vel,mass,
                                             theDomain,endnodetag);

            if(res < 0) {
                opserr << "WARNING: failed to discretize\n";
                return TCL_ERROR;
            }

            Tcl_SetObjResult(interp, Tcl_NewIntObj(endnodetag));

        } else if(strcmp(argv[2], "triangle")==0 || strcmp(argv[2], "Triangle")==0) {

            if(argc < 13) {
                opserr<<"WARNING: wrong num of args -- ";
                opserr<<"PFEM2D discretize triangle startnodetag x1 y1 x2 y2 x3 y3 n1 n2 ndf ";
                opserr<<"<-fix (constrValues) -mass (massValues)>\n";
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
            Vector mass(argc-loc), vel(argc-loc);
            ID fix;
            int state=-1, numfix=0, nummass=0, numvel=0;
            while(loc < argc) {
                if(strcmp(argv[loc], "-fix") == 0) {
                    state = 1;
                } else if(strcmp(argv[loc], "-mass") == 0) {
                    state = 2;
                }  else if(strcmp(argv[loc], "-vel") == 0) {
                    state = 3;
                } else {
                    if(state == 1) {
                        int fixity;
                        if(Tcl_GetInt(interp, argv[loc], &fixity) != TCL_OK) {
                            opserr<<"WARNING: invalid constraint value "<<argv[loc];
                            opserr<<" -- PFEM2D discretize triangle\n";
                            return TCL_ERROR; 
                        }
                        fix[numfix++] = fixity;
                    } else if(state==2) {
                        double m;
                        if(Tcl_GetDouble(interp, argv[loc], &m) != TCL_OK) {
                            opserr<<"WARNING: invalid mass "<<argv[loc];
                            opserr<<" -- PFEM2D discretize triangle\n";
                            return TCL_ERROR; 
                        }
                        mass(nummass++) = m;
                    } else if(state==3) {
                        double v;
                        if(Tcl_GetDouble(interp, argv[loc], &v) != TCL_OK) {
                            opserr<<"WARNING: invalid vel "<<argv[loc];
                            opserr<<" -- PFEM2D discretize triangle\n";
                            return TCL_ERROR; 
                        }
                        vel(numvel++) = v;
                    }
                }
            
                loc++;
            }
            mass.resize(nummass);
            vel.resize(numvel);
            int endnode = startnode;

            int res = theMesher2D.discretize(startnode,x1,y1,x2,y2,x3,y3,n1,n2,ndf,fix,vel,mass,
                                             theDomain,endnode);
            if(res < 0) {
                opserr << "WARNING: failed to discretize\n";
                return TCL_ERROR;
            }

            Tcl_SetObjResult(interp, Tcl_NewIntObj(endnode));

        } else if(strcmp(argv[2], "circle")==0 || strcmp(argv[2], "Circle")==0) {

            if(argc < 11) {
                opserr<<"WARNING: wrong num of args -- ";
                opserr<<"PFEM2D discretize circle startnode xc yc r1 r2 nc nr ndf ";
                opserr<<"<-fix (constrValues) -mass (massValues)>\n";
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
            Vector mass(argc-loc), vel(argc-loc);
            ID fix;
            int state=-1, numfix=0, nummass=0, numvel=0;
            while(loc < argc) {
                if(strcmp(argv[loc], "-fix") == 0) {
                    state = 1;
                } else if(strcmp(argv[loc], "-mass") == 0) {
                    state = 2;
                } else if(strcmp(argv[loc], "-vel") == 0) {
                    state = 3;
                } else {
                    if(state == 1) {
                        int fixity;
                        if(Tcl_GetInt(interp, argv[loc], &fixity) != TCL_OK) {
                            opserr<<"WARNING: invalid constraint value "<<argv[loc];
                            opserr<<" -- PFEM2D discretize triangle\n";
                            return TCL_ERROR; 
                        }
                        fix[numfix++] = fixity;
                    } else if(state==2) {
                        double m;
                        if(Tcl_GetDouble(interp, argv[loc], &m) != TCL_OK) {
                            opserr<<"WARNING: invalid mass "<<argv[loc];
                            opserr<<" -- PFEM2D discretize triangle\n";
                            return TCL_ERROR; 
                        }
                        mass(nummass++) = m;
                    } else if(state==3) {
                        double v;
                        if(Tcl_GetDouble(interp, argv[loc], &v) != TCL_OK) {
                            opserr<<"WARNING: invalid vel "<<argv[loc];
                            opserr<<" -- PFEM2D discretize triangle\n";
                            return TCL_ERROR; 
                        }
                        vel(numvel++) = v;
                    }
                }
            
                loc++;
            }
            mass.resize(nummass);
            vel.resize(numvel);
            int endnode = startnode;

            int res = theMesher2D.discretize(startnode,xc,yc,r1,r2,nc,nr,ndf,fix,vel,mass,
                                             theDomain,endnode);
            if(res < 0) {
                opserr << "WARNING: failed to discretize\n";
                return TCL_ERROR;
            }

            Tcl_SetObjResult(interp, Tcl_NewIntObj(endnode));


        } else if(strcmp(argv[2], "frame")==0 || strcmp(argv[2], "Frame")==0) {

            if(argc < 9) {
                opserr<<"WARNING: wrong num of args -- ";
                opserr<<"PFEM2D discretize frame startnode type num nth nthfloor ndf ";
                opserr<<"<-fix (constrValues) -mass (massValues)>\n";
                return TCL_ERROR;
            }
            char type;
            int nth,nthfloor,ndf,startnode,num;
            if(Tcl_GetInt(interp, argv[3], &startnode) != TCL_OK) {
                opserr<<"WARNING: invalid startnod "<<argv[3]<<" -- PFEM2D discretize frame\n";
                return TCL_ERROR; 
            }
            type = argv[4][0];
            if(Tcl_GetInt(interp, argv[5], &num) != TCL_OK) {
                opserr<<"WARNING: invalid num "<<argv[5]<<" -- PFEM2D discretize frame\n";
                return TCL_ERROR; 
            }
            if(Tcl_GetInt(interp, argv[6], &nth) != TCL_OK) {
                opserr<<"WARNING: invalid nth "<<argv[6]<<" -- PFEM2D discretize frame\n";
                return TCL_ERROR; 
            }
            if(Tcl_GetInt(interp, argv[7], &nthfloor) != TCL_OK) {
                opserr<<"WARNING: invalid nthfloor "<<argv[7]<<" -- PFEM2D discretize frame\n";
                return TCL_ERROR; 
            }
            if(Tcl_GetInt(interp, argv[8], &ndf) != TCL_OK) {
                opserr<<"WARNING: invalid ndf "<<argv[8]<<" -- PFEM2D discretize frame\n";
                return TCL_ERROR; 
            }

            int loc = 9;
            Vector mass(argc-loc), vel(argc-loc);
            ID fix;
            int state=-1, numfix=0, nummass=0, numvel=0;
            while(loc < argc) {
                if(strcmp(argv[loc], "-fix") == 0) {
                    state = 1;
                } else if(strcmp(argv[loc], "-mass") == 0) {
                    state = 2;
                } else if(strcmp(argv[loc], "-vel") == 0) {
                    state = 3;
                } else {
                    if(state == 1) {
                        int fixity;
                        if(Tcl_GetInt(interp, argv[loc], &fixity) != TCL_OK) {
                            opserr<<"WARNING: invalid constraint value "<<argv[loc];
                            opserr<<" -- PFEM2D discretize frame\n";
                            return TCL_ERROR; 
                        }
                        fix[numfix++] = fixity;
                    } else if(state==2) {
                        double m;
                        if(Tcl_GetDouble(interp, argv[loc], &m) != TCL_OK) {
                            opserr<<"WARNING: invalid mass "<<argv[loc];
                            opserr<<" -- PFEM2D discretize frame\n";
                            return TCL_ERROR; 
                        }
                        mass(nummass++) = m;
                    } else if(state==3) {
                        double v;
                        if(Tcl_GetDouble(interp, argv[loc], &v) != TCL_OK) {
                            opserr<<"WARNING: invalid vel "<<argv[loc];
                            opserr<<" -- PFEM2D discretize frame\n";
                            return TCL_ERROR; 
                        }
                        vel(numvel++) = v;
                    }
                }
            
                loc++;
            }
            mass.resize(nummass);
            vel.resize(numvel);
            int endnode = startnode;
            int res = theMesher2D.discretize(startnode,type,num,nth,nthfloor,ndf,fix,vel,mass,
                                             theDomain,endnode);
            if(res < 0) {
                opserr << "WARNING: failed to discretize\n";
                return TCL_ERROR;
            }

            Tcl_SetObjResult(interp, Tcl_NewIntObj(endnode));

        }

    
            
    } else if(strcmp(argv[1], "doTriangulation") == 0) {
        if(argc < 4) {
            opserr << "WARNING: wrong num of args -- ";
            opserr << "PFEM2D doTriangulation alpha -nodes (start1 end1 ...)  ";
            opserr << "<-addnodes (start1 end1 ...) -PFEMElement2D starteletag rho mu b1 b2 (thk)>\n";
            return TCL_ERROR;
        }
        double alpha;
        if(Tcl_GetDouble(interp, argv[2], &alpha) != TCL_OK) {
            opserr<<"WARNING: invalid alpha "<<argv[2]<<" -- PFEM2D doTriangulation\n";
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
            } else if(strcmp(argv[loc], "-PFEMElement2D") == 0) {
                states.push_back(3);
                locs.push_back(loc);
            }
        }
        locs.push_back(argc);

        std::vector<int> nodes, addnodes;
        std::vector<double> pfemparams(4);
        double thk = 1.0;
        int starteletag;
        bool pfem = false;
        for(int i=0; i<(int)states.size(); i++) {
            int state = states[i];
            int l0 = locs[i]+1;
            int l1 = locs[i+1];
            if(state == 1) {          // nodes
                if(l1-l0<2) {
                    opserr<<"WARNING: wrong no of arguments -- PFEM2D doTriangulation\n";
                    return TCL_ERROR;
                }
                for(int loc=l0; loc<l1; loc++) {
                    int tag;
                    if(Tcl_GetInt(interp, argv[loc], &tag) != TCL_OK) {
                        opserr<<"WARNING: invalid node tag "<<argv[loc];
                        opserr<<" -- PFEM2D doTriangulation\n";
                        return TCL_ERROR; 
                    }
                    nodes.push_back(tag);
                }
            } else if(state == 2) {   // addnodes
                if(l1-l0<2) {
                    opserr<<"WARNING: wrong no of arguments -- PFEM2D doTriangulation\n";
                    return TCL_ERROR;
                }
                for(int loc=l0; loc<l1; loc++) {
                    int tag;
                    if(Tcl_GetInt(interp, argv[loc], &tag) != TCL_OK) {
                        opserr<<"WARNING: invalid node tag "<<argv[loc];
                        opserr<<" -- PFEM2D doTriangulation\n";
                        return TCL_ERROR; 
                    }
                    addnodes.push_back(tag);
                }
            } else if(state == 3) {   // PFEMElement2D
                if(l1-l0<5) {
                    opserr<<"WARNING: wrong no of arguments -- PFEM2D doTriangulation\n";
                    return TCL_ERROR;
                }
                pfem = true;
                if(Tcl_GetInt(interp, argv[l0], &starteletag) != TCL_OK) {
                    opserr<<"WARNING: invalid element tag "<<argv[l0];
                    opserr<<" -- PFEM2D doTriangulation\n";
                    return TCL_ERROR; 
                }
                for(int loc=l0+1; loc<l0+5; loc++) {
                    if(Tcl_GetDouble(interp, argv[loc], &pfemparams[loc-l0-1]) != TCL_OK) {
                        opserr<<"WARNING: invalid parameter "<<argv[loc];
                        opserr<<" -- PFEM2D doTriangulation\n";
                        return TCL_ERROR; 
                    }
                }
                if(l1-l0>5) {
                    if(Tcl_GetDouble(interp, argv[l0+5], &thk) != TCL_OK) {
                        opserr<<"WARNING: invalid thickness "<<argv[l0+5];
                        opserr<<" -- PFEM2D doTriangulation\n";
                        return TCL_ERROR; 
                    }
                }
            }
        }

        // triangulation
        int res;
        ID eles;
        if(pfem) {
            res = theMesher2D.doTriangulation(starteletag,nodes,
                                              alpha,addnodes,theDomain,
                                              pfemparams[0],pfemparams[1],
                                              pfemparams[2],pfemparams[3],thk);
        } else {
            res = theMesher2D.doTriangulation(nodes,alpha,
                                              addnodes,theDomain,eles);
            if(res < 0) {
                opserr<<"WARNING: failed to do triangulation -- ";
                opserr<<" -- PFEM2D doTriangulation\n";
                return TCL_ERROR; 
            }
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
            opserr << "PFEM2D save filename <-structure (start1 end1 ...)> ";
            return TCL_ERROR;            
        }

        if(argc>5 && strcmp(argv[3], "-structure")==0) {
            int numnodes = argc-4;
            ID snodes(numnodes);
            for(int i=0; i<numnodes; i++) {
                if(Tcl_GetInt(interp, argv[4+i], &snodes(i)) != TCL_OK) {
                    opserr<<"WARNING: invalid strcutural node "<<argv[4+i]<<" -- PFEM2D save\n";
                    return TCL_ERROR; 
                }
            }
            theMesher2D.save(argv[2], snodes, theDomain);
        
        } else {
            ID snodes;
            theMesher2D.save(argv[2], snodes, theDomain);
        }


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

    } else if(strcmp(argv[1], "setFrame") == 0) {
        if(argc < 5) {
            opserr << "WARNING: wrong num of args -- ";
            opserr << "PFEM2D setFrame xbase ybase -span () -height ()\n";
            return TCL_ERROR;
        }
        double x1, y1;
        if(Tcl_GetDouble(interp, argv[2], &x1) != TCL_OK) {
            opserr<<"WARNING: invalid node x1 "<<argv[2];
            opserr<<" -- PFEM2D setFrame\n";
            return TCL_ERROR; 
        }
        if(Tcl_GetDouble(interp, argv[3], &y1) != TCL_OK) {
            opserr<<"WARNING: invalid node y1 "<<argv[3];
            opserr<<" -- PFEM2D setFrame\n";
            return TCL_ERROR; 
        }
        int loc = 4;
        std::vector<double> Lspan, Height;
        int state=-1;
        while(loc < argc) {
            if(strcmp(argv[loc], "-span") == 0) {
                state = 1;
            } else if(strcmp(argv[loc], "-height") == 0) {
                state = 2;
            } else {
                if(state == 1) {
                    double L;
                    if(Tcl_GetDouble(interp, argv[loc], &L) != TCL_OK) {
                        opserr<<"WARNING: invalid L "<<argv[loc];
                        opserr<<" -- PFEM2D discretize setFrame\n";
                        return TCL_ERROR; 
                    }
                    Lspan.push_back(L);
                } else if(state==2) {
                    double H;
                    if(Tcl_GetDouble(interp, argv[loc], &H) != TCL_OK) {
                        opserr<<"WARNING: invalid H "<<argv[loc];
                        opserr<<" -- PFEM2D discretize setFrame\n";
                        return TCL_ERROR; 
                    }
                    Height.push_back(H);
                }
            }
            
            loc++;
        }

        theMesher2D.setFrame(x1,y1,Lspan,Height);

    } else if(strcmp(argv[1], "removeOutBoundNodes") == 0) {

        if(argc < 2) {
            opserr << "WARNING: wrong num of args -- ";
            opserr << "PFEM2D removeOutBoundNodes <-nodes startnode end ...>\n";
            return TCL_ERROR;            
        }

        if(argc>4 && strcmp(argv[3], "-nodes")==0) {
            int numnodes = argc-3;
            ID nodes(numnodes);
            for(int i=0; i<numnodes; i++) {
                if(Tcl_GetInt(interp, argv[3+i], &nodes(i)) != TCL_OK) {
                    opserr<<"WARNING: invalid node "<<argv[3+i]<<" -- PFEM2D removeOutBoundNodes\n";
                    return TCL_ERROR; 
                }
            }
            theMesher2D.removeOutBoundNodes(nodes, theDomain);
        
        }

    } else if(strcmp(argv[1], "calculateForces") == 0) {

        if(argc < 2) {
            opserr<<"WARNING: wrong num of args -- ";
            opserr<<"PFEM2D calculateForces -boundary $nd1 ... ";
            opserr<<"-basenode $nd -dragdir $dx $dy -liftdir $dx $dy"; 
            return TCL_ERROR;            
        }
        int loc = 2, state = -1;
        std::vector<int> boundary;
        int basenode = 0;
        Vector dragdir(2), liftdir(2);
        while(loc < argc) {
            if(strcmp(argv[loc],"-boundary") == 0) {
                state = 0;
                loc++;
            } else if(strcmp(argv[loc],"-basenode") == 0) {
                loc++;
                if(Tcl_GetInt(interp, argv[loc], &basenode) != TCL_OK) {
                    opserr<<"WARNING: invalid node "<<argv[loc]<<" -- PFEM2D calculateForces\n";
                    return TCL_ERROR; 
                }
                loc++;
            } else if(strcmp(argv[loc],"-dragdir") == 0) {
                loc++;
                for(int i=0; i<2; i++) {
                    double dir = 0;
                    if(Tcl_GetDouble(interp, argv[loc], &dir) != TCL_OK) {
                        opserr<<"WARNING: invalid dir "<<argv[loc]<<" -- PFEM2D calculateForces\n";
                        return TCL_ERROR; 
                    }
                    loc++;
                    dragdir(i) = dir;
                }
            } else if(strcmp(argv[loc],"-liftdir") == 0) {
                loc++;
                for(int i=0; i<2; i++) {
                    double dir = 0;
                    if(Tcl_GetDouble(interp, argv[loc], &dir) != TCL_OK) {
                        opserr<<"WARNING: invalid dir "<<argv[loc]<<" -- PFEM2D calculateForces\n";
                        return TCL_ERROR; 
                    }
                    loc++;
                    liftdir(i) = dir;
                }
            } else if(state == 0) {
                int tag = 0;
                if(Tcl_GetInt(interp, argv[loc], &tag) != TCL_OK) {
                    opserr<<"WARNING: invalid tag "<<argv[loc]<<" -- PFEM2D calculateForces\n";
                    return TCL_ERROR; 
                }
                loc++;
                boundary.push_back(tag);
            }
        }

        Vector forces = theMesher2D.calculateForces(boundary,basenode,dragdir,liftdir,theDomain);

        for(int i=0; i<forces.Size(); i++) {
            char buffer[60];
            sprintf(buffer, "%35.20f ", forces(i));
            Tcl_AppendResult(interp, buffer, NULL);
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
