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
                                                                        
// $Revision: 1.3 $
// $Date: 2003-02-25 23:32:40 $
// $Source: /usr/local/cvs/OpenSees/SRC/database/main.cpp,v $
                                                                        
                                                                        
// File: ~/model/main.C
//
// Written: fmk 12/95
// Revised:
//
// Purpose: This file is a driver to create a 2d plane-frame
// model and to perform a linear static analysis on that model.
// 
//

#include <stdlib.h>

#include <OPS_Globals.h>
#include <Timer.h>
#include <FileDatastore.h>
#include <Quick2dFrame.h>
#include <Node.h>
#include <NodeIter.h>
#include <Element.h>
#include <ElementIter.h>

#include <Domain.h>
#include <ID.h>
#include <Vector.h>
#include <Matrix.h>
#include <FEM_ObjectBroker.h>

#include <StandardStream.h>
StandardStream sserr;
OPS_Stream &opserr = sserr;

int main(int argc, char **argv)
{

    //
    //	now create a domain and a modelbuilder
    //  and build the model
    //
    FEM_ObjectBroker theBroker;
    Domain *theDomain = new Domain();
    ModelBuilder  *theModelBuilder = new Quick2dFrame(*theDomain);

    FileDatastore *theDatabase 
      = new FileDatastore("/usr/users/fmk/g3/SRC/database/test1",*theDomain,theBroker);
    FileDatastore &theDb = *theDatabase;

  /*
    opserr << "TESTING IDs: \n";
    ID id1(2);
    id1(0) = 1; id1(1) = 2;


    ID id2(2);
    id2(0) = 3; id2(1) = 4;

    ID id3(2);
    id3(0) = 5; id3(1) = 6;

    ID id4(2);
    id4(0) = 7; id4(1) = 8;

    theDb.sendID(1,1,id1);
    theDb.sendID(1,2,id2);
    theDb.sendID(2,1,id3);
    theDb.sendID(2,2,id4);
    theDb.sendID(1,1,id4);
    theDb.sendID(1,2,id3);
    theDb.sendID(2,1,id2);
    theDb.sendID(2,2,id1);

    opserr << "RESULTS\n";    
    ID id5(2);
    theDb.recvID(1,1,id5);
    opserr << id5;
    theDb.recvID(1,2,id5);
    opserr << id5;
    theDb.recvID(2,1,id5);
    opserr << id5;
    theDb.recvID(2,2,id5);
    opserr << id5;
    theDb.recvID(1,1,id5);
    opserr << id5;
    theDb.recvID(1,2,id5);
    opserr << id5;
    theDb.recvID(2,1,id5);
    opserr << id5;
    theDb.recvID(2,2,id5);
    opserr << id5;

    theDb.sendID(2,2,id1);
    theDb.sendID(2,1,id2);
    theDb.sendID(1,2,id3);
    theDb.sendID(1,1,id4);

    theDb.recvID(1,1,id5);
    opserr << id5;
    theDb.recvID(1,2,id5);
    opserr << id5;
    theDb.recvID(2,1,id5);
    opserr << id5;
    theDb.recvID(2,2,id5);
    opserr << id5;

    theDb.sendID(4,4,id5);
    theDb.recvID(4,4,id5);
    opserr << id5;

    theDb.recvID(5,5,id5);
    opserr << id5;
    */
  
    /**************************
    opserr << "TESTING VECTORS: \n";

    Vector vect1(2);
    vect1(0) = 1.11111; vect1(1) = 2.222e3;

    Vector vect2(2);
    vect2(0) = 3.33333; vect2(1) = 4.2;

    Vector vect3(2);
    vect3(0) = 5; vect3(1) = 6;

    Vector vect4(2);
    vect4(0) = 7; vect4(1) = 8.8e12;

    theDb.sendVector(1,1,vect1);
    theDb.sendVector(1,2,vect2);
    theDb.sendVector(2,1,vect3);
    theDb.sendVector(2,2,vect4);

    opserr << "RESULTS\n";    
    Vector vect5(2);
    theDb.recvVector(1,1,vect5);
    opserr << vect5;
    theDb.recvVector(1,2,vect5);
    opserr << vect5;
    theDb.recvVector(2,1,vect5);
    opserr << vect5;
    theDb.recvVector(2,2,vect5);
    opserr << vect5;

    theDb.sendVector(2,2,vect1);
    theDb.sendVector(2,1,vect2);
    theDb.sendVector(1,2,vect3);
    theDb.sendVector(1,1,vect4);

    theDb.recvVector(1,1,vect5);
    opserr << vect5;
    theDb.recvVector(1,2,vect5);
    opserr << vect5;
    theDb.recvVector(2,1,vect5);
    opserr << vect5;
    theDb.recvVector(2,2,vect5);
    opserr << vect5;


    theDb.sendVector(4,4,vect5);
    theDb.recvVector(4,4,vect5);
    opserr << vect5;

    theDb.recvVector(5,5,vect5);
    opserr << vect5;
    **************************/
    /*
    opserr << "TESTING MATRICES: \n";

    Matrix mat1(2,2);
    mat1(0,0) = 1.11111; mat1(0,1) = 2.222e3;
    mat1(1,0) = 3.3; mat1(1,1) = 4.44;


    Matrix mat2(mat1);
    mat2 *= 10;
    theDb.sendMatrix(2,1,mat2);

    theDb.sendMatrix(1,1,mat1);

    Matrix mat3(2,2);
    theDb.recvMatrix(1,1,mat3);
    opserr << mat3;
    

    theDb.recvMatrix(2,1,mat3);
    opserr << mat3;

    theDb.sendMatrix(2,1,mat1);
    theDb.recvMatrix(2,1,mat3);
    opserr << mat3;
    */

//  theModelBuilder->buildFE_Model();
//  theDb.commitState(0);
  theDb.restoreState(0);


    // theDb.restoreElements(0);

    /*
    theDb.restoreNode(1,0);
    theDb.restoreNode(2,0);
    theDb.restoreNode(3,0);
    theDb.restoreNode(4,0);
    theDb.restoreElement(0,0);
    theDb.restoreElement(1,0);
    theDb.restoreElement(2,0);
    */

  opserr << *theDomain;

    delete theDatabase;

    exit(0);
}	
	
