/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 2001, The Regents of the University of California    **
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
** Reliability module developed by:                                   **
**   Terje Haukaas (haukaas@ce.berkeley.edu)                          **
**   Armen Der Kiureghian (adk@ce.berkeley.edu)                       **
**                                                                    **
** ****************************************************************** */
                                                                        
// $Revision: 1.5 $
// $Date: 2008-05-27 20:04:30 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/designPoint/FindDesignPointAlgorithm.cpp,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu)
//

#include <FindDesignPointAlgorithm.h>
#include <ReliabilityDomain.h>

 
FindDesignPointAlgorithm::FindDesignPointAlgorithm()
{

}

FindDesignPointAlgorithm::~FindDesignPointAlgorithm()
{

}


void
FindDesignPointAlgorithm::set_u(Vector& v)
{
	opserr << "FindDesignPointAlgorithm::set_u() - " << endln
           << " this function is for newalgorithm" << endln;
}

void
FindDesignPointAlgorithm::set_x(Vector& x)
{
	opserr << "FindDesignPointAlgorithm::set_x() - " << endln
           << " this function is for newalgorithm" << endln;
}

double
FindDesignPointAlgorithm::get_beta()
{
	opserr << "FindDesignPointAlgorithm::get_beta() - " << endln
           << " this function is for newalgorithm" << endln;
	return 0;
}

Matrix
FindDesignPointAlgorithm::getJacobian_x_u()
{
	opserr << "FindDesignPointAlgorithm::getJacobian_x_u() - " << endln
           << " this function is for newalgorithm" << endln;
	Matrix* zero=0;
	return (*zero);
}

int
FindDesignPointAlgorithm::getNumberOfSensAna()
{
	opserr << "FindDesignPointAlgorithm::getNumberOfSensAna() - " << endln
           << " this function is for newalgorithm" << endln;
	return 0;
}

double
FindDesignPointAlgorithm::get_check1_init()
{
	opserr << "FindDesignPointAlgorithm::get_check1_init() - " << endln
           << " this function is for newalgorithm" << endln;
	return 0;
}

double
FindDesignPointAlgorithm::get_check2_init()
{
	opserr << "FindDesignPointAlgorithm::get_check2_init() - " << endln
           << " this function is for newalgorithm" << endln;
	return 0;
}

double
FindDesignPointAlgorithm::get_check1_conv()
{
	opserr << "FindDesignPointAlgorithm::get_check1_conv() - " << endln
           << " this function is for newalgorithm" << endln;
	return 0;
}

double
FindDesignPointAlgorithm::get_check2_conv()
{
	opserr << "FindDesignPointAlgorithm::get_check2_conv() - " << endln
           << " this function is for newalgorithm" << endln;
	return 0;
}
