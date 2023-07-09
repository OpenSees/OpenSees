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
// $Date: 2016-1-27  $
                                                                        
// Written: Minjie Zhu
//
// Description: This class defines the BackgroundMesh 
//

#include "BackgroundFixData.h"
#include <Node.h>
#include <SP_Constraint.h>

BackgroundFixData::BackgroundFixData()
    :min(),max(),fix(),sps()
{
}

BackgroundFixData::~BackgroundFixData()
{

}

void
BackgroundFixData::addInfo(const Vector& min, const Vector& max, const ID& fix)
{
    bool fixed = false;
    for (int i=0; i<fix.Size(); i++) {
	if (fix(i) != 0) {
	    fixed = true;
	    break;
	}
    }
    if (!fixed) return;
    
    this->min.push_back(min);
    this->max.push_back(max);
    this->fix.push_back(fix);
}

int
BackgroundFixData::tryFix(int ndtag, Domain& domain)
{
    // get node
    Node* node = domain.getNode(ndtag);
    if (node == 0) {
	opserr<<"WARNING: node "<<ndtag<<" not exist\n";
	return -1;
    }
    int ndf = node->getNumberDOF();

    // crds
    const Vector& crds = node->getCrds();

    // check each info
    int index = -1;
    for (int i=0; i<(int)min.size(); i++) {
	if (crds.Size() != min[i].Size()) {
	    opserr << "WARNING: ndm for the nodes and fix range are not compatible\n";
	    return -1;
	}
	if (crds.Size() != max[i].Size()) {
	    opserr << "WARNING: ndm for the nodes and fix range are not compatible\n";
	    return -1;
	}

	// if in the fix area
	bool fixed = true;
	for (int j=0; j<crds.Size(); j++) {
	    if (crds(j)<min[i](j) || crds(j)>max[i](j)) {
		fixed = false;
		break;
	    }
	}
	if (fixed) {
	    index = i;
	    break;
	}
    }

    // if not fixed
    if (index == -1) return 0;

    // set nodal states
    Vector disp = node->getTrialDisp();
    Vector vel = node->getTrialVel();
    Vector accel = node->getTrialVel();
    for (int i=0; i<ndf; i++) {
	if (i>=fix[index].Size()) break;
	if (fix[index](i) != 0) {
	    disp(i) = 0.0;
	    vel(i) = 0.0;
	    accel(i) = 0.0;
	}
    }
    node->setTrialDisp(disp);
    node->setTrialVel(vel);
    node->setTrialAccel(accel);
    node->commitState();

    // add SPs
    for (int i=0; i<ndf; i++) {
	if (i>=fix[index].Size()) break;
	if (fix[index](i) != 0) {
	    SP_Constraint* sp = new SP_Constraint(ndtag, i, 0, true);
	    if (sp == 0) {
		opserr<<"WARNING: run out of memory\n";
		return -1;
	    }
	    if (domain.addSP_Constraint(sp) == false) {
		opserr<<"WARNING: failed to add sp to domain\n";
		delete sp;
		return -1;
	    }
	    sps.push_back(sp);
	}
    }
	
    return 1;
}

void
BackgroundFixData::clear(Domain& domain)
{
    for (int i=0; i<(int)sps.size(); i++) {
	if (sps[i] != 0) {
	    domain.removeSP_Constraint(sps[i]->getTag());
	    delete sps[i];
	}
    }
    sps.clear();
}
