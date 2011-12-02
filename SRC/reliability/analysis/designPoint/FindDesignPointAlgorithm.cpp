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
                                                                        
// $Revision: 1.3 $
// $Date: 2008-03-13 22:30:16 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/designPoint/FindDesignPointAlgorithm.cpp,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu)
//

#include <FindDesignPointAlgorithm.h>

 
FindDesignPointAlgorithm::FindDesignPointAlgorithm()
{
	outputFile =0;
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
GradGEvaluator*
FindDesignPointAlgorithm::getGradGEvaluator()
{
	opserr << "FindDesignPointAlgorithm::getGradGEvaluator() - " << endln
           << " this function is for newalgorithm" << endln;
	return 0;
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
