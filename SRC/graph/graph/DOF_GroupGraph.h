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
                                                                        
// $Revision: 1.1.1.1 $
// $Date: 2000-09-15 08:23:21 $
// $Source: /usr/local/cvs/OpenSees/SRC/graph/graph/DOF_GroupGraph.h,v $
                                                                        
                                                                        
// File: ~/graph/graph/DOF_GroupGraph.h
// 
// Written: fmk 
// Created: Sun Sept 15 11:47:47: 1996
// Revision: A
//
// Description: This file contains the class definition for DOF_GroupGraph.
// DOF_GroupGraph is a graph of the DOF_Groups in the domain. It is used by 
// the DOF_Numberer to assign equation numbers to the DOFs.
//
// What: "@(#) DOF_GroupGraph.h, revA"

#ifndef DOF_GroupGraph_h
#define DOF_GroupGraph_h

#include <Graph.h>

class AnalysisModel;

class DOF_GroupGraph: public Graph
{
  public:
    DOF_GroupGraph(AnalysisModel &theModel);
    ~DOF_GroupGraph();

  protected:
    
  private:
    AnalysisModel &myModel;
};

#endif

