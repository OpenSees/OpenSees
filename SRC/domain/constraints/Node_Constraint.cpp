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

// Revision : 1.0.0
// Initial Date : 2022-07-16
// Current Rev. : 2022-07-16


#include <DomainComponent.h>
#include <Vector.h>
#include <Node_Constraint.h>

Node_Constraint::Node_Constraint(int nodeTag, Vector& spVal, bool isLoadConstant = false) :DomainComponent(0, CNSTRNT_TAG_NODE_SP_Constraint), NodeTag(nodeTag), value(0), loadConstant(isLoadConstant)
{

}