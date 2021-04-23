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
// $Date: 2010-04-23 22:50:19 $
// $Source: /usr/local/cvs/OpenSees/SRC/domain/constraints/BeamSolidMPC.cpp,v $
                                                                        
                                                                        
// File: ~/model/constraints/BS3d.C
//
// Written: Ziping Zhu,zzpxbkarl@gmail.com
// Revised:
//
// Purpose: This file contains the class definition for BS3d which interacts element
// of beam and solid in 3D.
// BSMPC  $rNodeTag1  $rNodeTag2  $cNodeTag1  $cNodeTag2…$cNodeTag_n $v1 $v2 $v3
// $rNodeTag1 is a node of Beam-column element in multiscale interface.
// $rNodeTag2 is a node of Beam-column element which is not in multiscale interface.
// These two nodes can determine the local coordination x'.
// $cNodeTag1  $cNodeTag2…$cNodeTag_n are nodes of solid element in multiscale interface.
// $v1 $v2 $v3 are the verctor to determine local coordination, which is not parallel to x'.
// local coordination: x'= $rNodeTag2 - $rNodeTag1, y'=v * x', z'=x' * y'

#include <OPS_Globals.h>
#include <stdlib.h>
#include <Domain.h>
#include <Node.h>
#include <MP_Constraint.h>
#include <Matrix.h>
#include <ID.h>
#include <BSMPC.h>
#include <math.h>

BSMPC::BSMPC(Domain &theDomain, int ndm,int nR1, int nR2, 
	 ID &nC,const Vector &vectorforlocal) {//R1,R2,C1,C2 are uesd to determine local coordination.

    // check constrainedNodes ID does not contain the retained node
    if (nC.getLocation(nR1) >= 0) {
      opserr << "BS3d::BS3d - " << 
	"retained node" << nR1 << "is in constrained node list\n";
      return;
    }	

    //R1 is a retained node in beam element which is near solid element, R2 in beam element is far from the solid element.
    //to define local x'y'z'.
    Node *nodeR1 = theDomain.getNode(nR1);
    if (nodeR1 == 0) {
      opserr << "BS3d::BS3d - " << 
	"retained Node" <<  nR1 <<  "not in domain\n";
      return;
    }
    Node *nodeR2 = theDomain.getNode(nR2);
    //Node *nodeC1 = theDomain.getNode(nC1);
    //Node *nodeC2 = theDomain.getNode(nC2);
    const Vector &crdR1 = nodeR1->getCrds();
    const Vector &crdR2 = nodeR2->getCrds();
    //const Vector &crdC1 = nodeC1->getCrds();
    //const Vector &crdC2 = nodeC2->getCrds();
    //Vector localx = crdC1 - crdR1;
    //Vector localy = crdC2 -crdR1;
    //Vector localz = crdR2 -crdR1;
    Vector localx = crdR2 -crdR1;
    Vector localy(3);
    localy(0) = vectorforlocal(1)*localx(2) - vectorforlocal(2)*localx(1);
    localy(1) = vectorforlocal(2)*localx(0) - vectorforlocal(0)*localx(2);
    localy(2) = vectorforlocal(0)*localx(1) - vectorforlocal(1)*localx(0);
    
    Vector localz(3);
    localz(0) = localx(1)*localy(2) - localx(2)*localy(1);
    localz(1) = localx(2)*localy(0) - localx(0)*localy(2);
    localz(2) = localx(0)*localy(1) - localx(1)*localy(0);
    //cosaerfa1=cos<x',x>, cosbeta1 = cos<x',y>, cosgama1 = cos<x',z>;
    //cosaerfa2 = cos<y',x>, cosbeta2 = cos<y',y>, cosgama2 = cos<y',z>;
    //cosaerfa3 = cos<z',x>, cosbeta3 = cos<z',y>, cosgama3 = cos<z',z>;
    //x' = x*cosaerfa1 + y*cosbeta1 + z*cosgama1;
    //y' = x*cosaerfa2 + y*cosbeta2 + z*cosgama2;
    //z' = x*cosaerfa3 + y*cosbeta3 + z*cosgama3;
    cosaerfa1 = localx(0)/sqrt(localx(0)*localx(0)+localx(1)*localx(1)+localx(2)*localx(2));
    cosbeta1 = localx(1)/sqrt(localx(0)*localx(0)+localx(1)*localx(1)+localx(2)*localx(2));
    cosgama1 = localx(2)/sqrt(localx(0)*localx(0)+localx(1)*localx(1)+localx(2)*localx(2));
    cosaerfa2 = localy(0)/sqrt(localy(0)*localy(0)+localy(1)*localy(1)+localy(2)*localy(2));
    cosbeta2 = localy(1)/sqrt(localy(0)*localy(0)+localy(1)*localy(1)+localy(2)*localy(2));
    cosgama2 = localy(2)/sqrt(localy(0)*localy(0)+localy(1)*localy(1)+localy(2)*localy(2));
    cosaerfa3 = localz(0)/sqrt(localz(0)*localz(0)+localz(1)*localz(1)+localz(2)*localz(2));
    cosbeta3 =localz(1)/sqrt(localz(0)*localz(0)+localz(1)*localz(1)+localz(2)*localz(2));
    cosgama3 =localz(2)/sqrt(localz(0)*localz(0)+localz(1)*localz(1)+localz(2)*localz(2));

    //to get constrians of each constrained nodes in solid element.
for (int i=0; i<nC.Size(); i++) {
    ID idC(3),idR(6);
	idC(0) = 0; idC(1) = 1; idC(2) = 2;
	idR(0) = 0; idR(1) = 1; idR(2) = 2; idR(3) = 3; idR(4) = 4; idR(5) = 5;

     Matrix T1(3,3); T1.Zero();//transformation matrix for displacement of constrained nodes in local xyz.
    Matrix T2(6,6); T2.Zero();//transformation matrix for displacement of retained nodes in local xyz.
    T1(0,0)= cosaerfa1; T1(0,1)= cosbeta1; T1(0,2)= cosgama1;
    T1(1,0)= cosaerfa2; T1(1,1)= cosbeta2; T1(1,2)= cosgama2;
    T1(2,0)= cosaerfa3; T1(2,1)= cosbeta3; T1(2,2)= cosgama3;
    T2(0,0)=cosaerfa1; T2(0,1)= cosbeta1; T2(0,2)= cosgama1;
    T2(1,0)=cosaerfa2; T2(1,1)= cosbeta2; T2(1,2)= cosgama2;
    T2(2,0)= cosaerfa3; T2(2,1)= cosbeta3; T2(2,2)=cosgama3;
    T2(3,3)= cosaerfa1; T2(3,4)= cosbeta1; T2(3,5)=cosgama1;
    T2(4,3)= cosaerfa2; T2(4,4)= cosbeta2; T2(4,5)=cosgama2;
    T2(5,3)= cosaerfa3; T2(5,4)= cosbeta3; T2(5,5)=cosgama3;

	 int ndC = nC(i);
      Node *nodeC = theDomain.getNode(ndC);
      if (nodeC != 0) {

	 // get constrained nodes coordinates in local.
	const Vector &crdC = nodeC->getCrds();
    double clocalx = (crdC(0)-crdR1(0))*cosaerfa1+ (crdC(1)-crdR1(1))*cosbeta1+ (crdC(2)-crdR1(2))*cosgama1; 
    double clocaly = (crdC(0)-crdR1(0))*cosaerfa2+ (crdC(1)-crdR1(1))*cosbeta2+ (crdC(2)-crdR1(2))*cosgama2;
    double clocalz = (crdC(0)-crdR1(0))*cosaerfa3+ (crdC(1)-crdR1(1))*cosbeta3+ (crdC(2)-crdR1(2))*cosgama3;

    Matrix ci(3,6);ci.Zero();//uci'=ci*ur1'
    ci(0,0)= 1; ci(0,4)= clocalz; ci(0,5)=-clocaly;
    ci(1,1)= 1; ci(1,3)= -clocalz;
    ci(2,2)= 1; ci(2,3)= clocaly; 

 //Ci=T1'*ci*T2;
    Matrix Ci(3,6);Ci.Zero();//Uci=Ci*Ur1
    Ci(0,0)= cosaerfa1*cosaerfa1 + cosaerfa2*cosaerfa2 + cosaerfa3*cosaerfa3;
    Ci(0,1)= cosaerfa1*cosbeta1 + cosaerfa2*cosbeta2 + cosaerfa3*cosbeta3;
    Ci(0,2)= cosaerfa1*cosgama1 + cosaerfa2*cosgama2 + cosaerfa3*cosgama3;
    Ci(0,3)= (-clocalz*cosaerfa2 + clocaly*cosaerfa3)*cosaerfa1 + (clocalz*cosaerfa1)*cosaerfa2 +(-clocaly*cosaerfa1)*cosaerfa3;
    Ci(0,4)= (-clocalz*cosaerfa2 + clocaly*cosaerfa3)*cosbeta1 + (clocalz*cosaerfa1)*cosbeta2 + (-clocaly*cosaerfa1)*cosbeta3;
    Ci(0,5)= (-clocalz*cosaerfa2 + clocaly*cosaerfa3)*cosgama1 + (clocalz*cosaerfa1)*cosgama2 + (-clocaly*cosaerfa1)*cosgama3;
    Ci(1,0)= cosbeta1*cosaerfa1 + cosbeta2*cosaerfa2 + cosbeta3*cosaerfa3;
    Ci(1,1)= cosbeta1*cosbeta1 + cosbeta2*cosbeta2 + cosbeta3*cosbeta3;
    Ci(1,2)= cosbeta1*cosgama1 + cosbeta2*cosgama2 + cosbeta3*cosgama3;
    Ci(1,3)= (-clocalz*cosbeta2 + clocaly*cosbeta3)*cosaerfa1 + (clocalz*cosbeta1)*cosaerfa2 + (-clocaly*cosbeta1)*cosaerfa3;
    Ci(1,4)= (-clocalz*cosbeta2 + clocaly*cosbeta3)*cosbeta1 + (clocalz*cosbeta1)*cosbeta2 + (-clocaly*cosbeta1)*cosbeta3;
    Ci(1,5)= (-clocalz*cosbeta2 + clocaly*cosbeta3)*cosgama1 + (clocalz*cosbeta1)*cosgama2 + (-clocaly*cosbeta1)*cosgama3;
    Ci(2,0)= cosgama1*cosaerfa1 + cosgama2*cosaerfa2 + cosgama3*cosaerfa3;
    Ci(2,1)= cosgama1*cosbeta1 + cosgama2*cosbeta2 + cosgama3*cosbeta3;
    Ci(2,2)= cosgama1*cosgama1 + cosgama2*cosgama2 + cosgama3*cosgama3;
    Ci(2,3)= (-clocalz*cosgama2 + clocaly*cosgama3)*cosaerfa1 + (clocalz*cosgama1)*cosaerfa2 + (-clocaly*cosgama1)*cosaerfa3;
    Ci(2,4)= (-clocalz*cosgama2 + clocaly*cosgama3)*cosbeta1 + (clocalz*cosgama1)*cosbeta2 + (-clocaly*cosgama1)*cosbeta3;
    Ci(2,5)= (-clocalz*cosgama2 + clocaly*cosgama3)*cosgama1 + (clocalz*cosgama1)*cosgama2 + (-clocaly*cosgama1)*cosgama3;

//get constrains.
    MP_Constraint *newC = new MP_Constraint(nR1, ndC, Ci, idC, idR);
    if (newC == 0) {
	 opserr << "BS3d::BS3d - ignoring constrained Node " << ndC << 
	   ", out of memory\n";
	} else {
	    // add the constraint to the domain
	  if (theDomain.addMP_Constraint(newC) == false) {
	  opserr << "BS3d::BS3d - ignoring constrained Node " << ndC << 
	   ", failed to add\n";
	   delete newC;}    
	  }
   }else{
   	opserr << "BS3d::BS3d nodeC is not exist.";
   }
  
}
}

BSMPC::~BSMPC()
{
    // does nothing
}
//
