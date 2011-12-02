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
                                                                        
// $Revision: 1.2 $
// $Date: 2006-01-12 23:37:19 $
// $Source: /usr/local/cvs/OpenSees/SRC/graph/partitioner/Metis.h,v $
                                                                        
                                                                        
// Written: fmk 
//
// Description: This file contains the class definition for Metis.
// Metis is a type of GraphPartitioner which uses 'METIS - Unstructured
// Graph Partitioning And Sparse Matrix Ordering System', developed by
// G. Karypis and V. Kumar at the University of Minnesota. The metis
// files are found in metis-2.0 which were downloaded.
//     This class provides the C++ interface for metis which will allow
// it to fit seamlessly into our system.
//
// What: "@(#) Metis.h, revA"

#ifndef Metis_h
#define Metis_h

#include <GraphPartitioner.h>
#include <GraphNumberer.h>

#ifndef _bool_h
#include <bool.h>
#endif

class Metis : public GraphPartitioner, public GraphNumberer
{
  public:
    Metis(int numParts =1);
    Metis(int Ptype, 
	  int Mtype, 
	  int coarsenTo,
	  int Rtype, 
	  int IPtype,
	  int numParts =1);
    ~Metis();

    bool setOptions(int Ptype, 
		    int Mtype,
		    int coarsenTo,
		    int Rtype, 
		    int IPtype);

    bool setDefaultOptions(void);
    
    int partition(Graph &theGraph, int numPart);

    // the follwing methods are if the object is to be used as a numberer
    const ID &number(Graph &theGraph, int lastVertex = -1);
    const ID &number(Graph &theGraph, const ID &lastVertices);

    int sendSelf(int commitTag, 
		 Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, 
		 FEM_ObjectBroker &theBroker);

    
  protected:

  private:
    bool checkOptions(void);
    
    int myPtype ; 	// package type: 
                        //	pmetis = 1
                        //	kmetis = 2

    int myMtype;     	// type of matching scheme: 
			//	random = 1
			//  	heavy edge = 2
			//    	light edge = 3
			//	heavy clique = 4
			//	modified heavy edge = 5
			//    	sorted random = 11
			//	sorted heavy edge =21
			// 	sorted modified heavy edge = 51
   
   int myCoarsenTo; 	// the number of vertices the graph is coarsened down to
                        //   	if pmetis default is 100
			//	if kmetis default is 2000
	
    int myRtype;     	// type of refinement policy:  
			//	greedy = 1
			//	kernighan-lin = 2
			//    	combo greedy and K-L = 3
			// 	boundary greedy = 11
			//    	boundary K-L = 12
			//	combo of boundary greedy and boundary K-L = 13,
			//    	no-refinement = 20

    int myIPtype; 	// type of bisection algo:  
	                //    	graph growing partition = 1,
			//    	greedy graph growing partition = 2,
			//    	spectral bisection = 3,
			//    	graph growing followed by K-L = 4
	
    bool defaultOptions;			    
    
    int   numPartitions; // needed if to be used as a numberer
    ID    *theRefResult;
};

#endif

