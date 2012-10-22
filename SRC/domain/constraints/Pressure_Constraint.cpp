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
// $Date: 2012-08-21 12:59:06 $
// $Source: /usr/local/cvs/OpenSees/SRC/domain/constraints/Pressure_Constraint.h,v $

// Written: Minjie
// Created: 08/12
// Revision: A                                                             
//
// Purpose: This file contains the implementation of class Pressure_Constraint.
//
// The class Pressure_Constraint interface:
//

#include <Pressure_Constraint.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>

Pressure_Constraint::Pressure_Constraint(int classTag)
    :DomainComponent(0,classTag), nodeConstrained(0), fluid(true)
{
}

Pressure_Constraint::Pressure_Constraint(int classTag, int nodeId, bool Fluid)
    :DomainComponent(nodeId,classTag), nodeConstrained(nodeId), fluid(Fluid)
{
}

Pressure_Constraint::Pressure_Constraint(int nodeId, bool Fluid)
    :DomainComponent(nodeId,CNSTRNT_TAG_Pressure_Constraint), nodeConstrained(nodeId),
     fluid(Fluid)
{
}

Pressure_Constraint::~Pressure_Constraint()
{
}

int Pressure_Constraint::getNodeConstrained()const 
{
    return nodeConstrained;
}

bool Pressure_Constraint::isFluid()const
{
    return fluid;
}

int Pressure_Constraint::sendSelf(int commitTag, Channel& theChannel)
{
    static Vector data(3);
    data(0) = this->getTag();
    data(1) = nodeConstrained;
    data(2) = fluid;

    int result = theChannel.sendVector(this->getDbTag(), commitTag, data);

    if (result != 0) {
        opserr << "WARNING Pressure_Constraint::sendSelf - error sending Vector data\n";
        return result;
    }

    return 0;
}

int Pressure_Constraint::recvSelf(int commitTag, Channel &theChannel, 
                                  FEM_ObjectBroker &theBroker)
{
    static Vector data(3);

    int result = theChannel.recvVector(this->getDbTag(), commitTag, data);
    if (result < 0) {
	opserr << "WARNING Pressure_Constraint::recvSelf - error receiving Vector data\n";
	return result;
    }

    // if o.k. set the data
    this->setTag((int)data(0));
    nodeConstrained = (int)data(1);
    fluid = (data(2) != 0);

    return 0;
}

void Pressure_Constraint::Print(OPS_Stream& s, int flag)
{
    s << "Pressure_Constraint: " << this->getTag() << "    ";
    s << "Constrained Node: " << nodeConstrained << "    ";
    if(fluid) {
        s << "This is " << "fluid node\n";
    } 
}
