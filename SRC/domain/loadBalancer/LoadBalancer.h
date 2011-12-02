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
// $Date: 2000-09-15 08:23:19 $
// $Source: /usr/local/cvs/OpenSees/SRC/domain/loadBalancer/LoadBalancer.h,v $
                                                                        
                                                                        
// File: ~/domain/loadBalancer/LoadBalancer.h
// 
// Written: fmk 
// Created: Fri Aug 29 17:43:25 1997
// Revision: A
//
// Description: This file contains the class definition for LoadBalancer.
// A LoadBalancer is an object used to partition a PartitionedDomain.
//
// What: "@(#) LoadBalancer.h, revA"

#ifndef LoadBalancer_h
#define LoadBalancer_h

#include <DomainPartitioner.h>
class Vector;


class LoadBalancer
{
  public:
    
    LoadBalancer();
    virtual  ~LoadBalancer();    

    virtual void setLinks(DomainPartitioner &thePartitioner);
    virtual int balance(Graph &theWeightedGraph) =0;

  protected:    
    DomainPartitioner *getDomainPartitioner(void);
	
  private:
    DomainPartitioner *theDomainPartitioner;    
};

#endif


