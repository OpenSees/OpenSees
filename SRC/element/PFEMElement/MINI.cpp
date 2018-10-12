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

// $Revision $
// $Date$
// $URL$

// Written: Minjie Zhu (zhum@oregonstate.edu)
//
// Description: MINI element which mass matrix can be lumped
//

#include "MINI.h"
#include <elementAPI.h>
#include <Node.h>
#include <Domain.h>
#include <ElementResponse.h>
#include <math.h>
#include <NodeIter.h>
#include <Pressure_Constraint.h>

Matrix MINI::mat;
Vector MINI::vec;

void* OPS_MINI()
{
    int ndm = OPS_GetNDM();
    int num = OPS_GetNumRemainingInputArgs();
    if (ndm == 2) {
	if(num < 10) {
	    opserr<<"WARNING: insufficient number of arguments -- MINI:";
	    opserr<<"tag, nd1, nd2, nd3, rho, mu, b1, b2, thk, kappa\n";
	    return 0;
	}
	
    } else if (ndm == 3) {
	if(num < 12) {
	    opserr<<"WARNING: insufficient number of arguments -- MINI:";
	    opserr<<"tag, nd1, nd2, nd3, nd4, rho, mu, b1, b2, b3, thk, kappa\n";
	    return 0;
	}
    }

    // tag, nd1, nd2, nd3
    if (ndm == 2) {
	num = 4;
    } else if (ndm == 3) {
	num = 5;
    }
    int idata[5];
    if(OPS_GetIntInput(&num,idata)<0) {
	opserr << "WARNING: failed to read integers -- MINI\n";
	return 0;
    }

    // (rho, mu, b1, b2,b3 thinkness, kappa)
    if (ndm == 2) {
	num = 6;
    } else if (ndm == 3) {
	num = 7;
    }
    double data[7];
    if(OPS_GetDoubleInput(&num,data) < 0) {
	opserr << "WARNING: failed to read doubles -- MINI\n";
	return 0;
    }

    if (ndm == 2) {

	return new MINI(idata[0],idata[1],idata[2],idata[3],
			data[0],data[1],data[2],data[3],
			data[4],data[5]);
    } else if (ndm == 3) {
	return new MINI(idata[0],idata[1],idata[2],idata[3],idata[4],
			data[0],data[1],data[2],data[3],
			data[4],data[5],data[6]);
    }

    return 0;
}


MINI::MINI()
    :Element(0, ELE_TAG_MINI),
     ntags(), nodes(), data(), dofs(1), Jco()
{
}

MINI::MINI(int tag, int nd1, int nd2, int nd3,
	   double r, double m,
	   double bx, double by, 
	   double t, double ka)
    :Element(tag, ELE_TAG_MINI),
     ntags(7), nodes(7), data(6), dofs(8), Jco(3,3)
{
    if (ka <= 0) {
	ntags.resize(6);
	nodes.resize(6);
	dofs.resize(7);
    }
    ntags(0) = nd1;
    ntags(1) = nd1;
    ntags(2) = nd2;
    ntags(3) = nd2;
    ntags(4) = nd3;
    ntags(5) = nd3;
    if (ka > 0) {
	ntags(6) = nd3;
    }
    data(0) = r;
    data(1) = m;
    data(2) = t;
    data(3) = ka;
    data(4) = bx;
    data(5) = by;

}

MINI::MINI(int tag, int nd1, int nd2, int nd3, int nd4,
	   double r, double m,
	   double bx, double by, double bz, 
	   double t, double ka)
    :Element(tag, ELE_TAG_MINI),
     ntags(9), nodes(9), data(7), dofs(10), Jco(4,4)
{
    if (ka <= 0) {
	ntags.resize(8);
	nodes.resize(8);
	dofs.resize(9);
    }

    ntags(0) = nd1;
    ntags(1) = nd1;
    ntags(2) = nd2;
    ntags(3) = nd2;
    ntags(4) = nd3;
    ntags(5) = nd3;
    ntags(6) = nd4;
    ntags(7) = nd4;
    if (ka > 0) {
	ntags(8) = nd4;
    }
    data(0) = r;
    data(1) = m;
    data(2) = t;
    data(3) = ka;
    data(4) = bx;
    data(5) = by;
    data(6) = bz;
}

// destructor
MINI::~MINI()
{
    Domain* domain = this->getDomain();
    
    int nump = ntags.Size() / 2;
    for(int i=0; i<nump; i++) {
	Pressure_Constraint* pc = domain->getPressure_Constraint(ntags(2*i));
        if(pc != 0) {
            pc->disconnect(this->getTag());
        }
    }
    if (data(3) > 0) {
	Node* bnode = nodes.back();
	if (domain != 0) {
	    if (bnode != 0) {
		domain->removeNode(bnode->getTag());
		delete bnode;
	    }
	}
    }
}

// methods dealing with nodes and number of external dof
int
MINI::getNumExternalNodes() const
{
    return ntags.Size();
}

const ID &
MINI::getExternalNodes()
{
    return ntags;
}

Node **
MINI::getNodePtrs()
{
    return &nodes[0];
}

int
MINI::getNumDOF()
{
    return dofs(dofs.Size()-1);
}

// public methods to set the state of the element
int
MINI::revertToLastCommit()
{
    return 0;
}

int
MINI::update()
{
    // get Jacobian matrix
    Matrix J(Jco.noRows(),Jco.noCols());
    for (int b=0; b<J.noCols(); ++b) {
    	const Vector& coord = nodes[2*b]->getCrds();
    	const Vector& disp = nodes[2*b]->getTrialDisp();

	for (int a=0; a<J.noRows(); ++a) {
	    if (a == 0) {
		J(a,b) = 1.0;
	    } else {
		J(a,b) = coord(a-1) + disp(a-1);
	    }
	}
    }

    // cofactor of Jacobian matrix
    Matrix sub(J.noRows()-1,J.noCols()-1);
    for (int a=0; a<J.noRows(); ++a) {
	for (int b=0; b<J.noCols(); ++b) {
	    
	    // fill submatrix
	    for (int m=0; m<sub.noRows(); ++m) {
		for (int n=0; n<sub.noCols(); ++n) {
		    int mm = m;
		    int nn = n;
		    if (m >= a) {
			mm += 1;
		    }
		    if (n >= b) {
			nn += 1;
		    }
		    sub(m,n) = J(mm,nn);
		}
	    }

	    // det of sub
	    double detsub = det(sub);

	    // sign of cofactor
	    if ((a+b) % 2 > 0) {
		detsub *= -1.0;
	    }

	    // transpose cofactor
	    Jco(b,a) = detsub;
	}
    }

    double detJ = 0.0;
    for (int i=0; i<Jco.noRows(); ++i) {
	detJ += Jco(i,0);
    }

    if (fabs(detJ) <= 0) {
	opserr << "WARNING: J <= 0\n";
	return -1;
    }

    return 0;
}

int
MINI::commitState()
{
    return Element::commitState();
}

// public methods to obtain stiffness, mass, damping and residual information
const Matrix &
MINI::getTangentStiff()
{
    int ndf = getNumDOF();
    mat.resize(ndf,ndf);
    mat.Zero();
    return mat;
}

const Matrix &
MINI::getInitialStiff()
{
    int ndf = getNumDOF();
    mat.resize(ndf,ndf);
    mat.Zero();
    return mat;
}

const Matrix &
MINI::getMass()
{
    int ndf = getNumDOF();
    mat.resize(ndf,ndf);
    mat.Zero();


    int ndm = OPS_GetNDM();

    double rho = data(0);
    double thk = data(2);
    double kappa = data(3);
    double J = 0.0;
    for (int i=0; i<Jco.noRows(); ++i) {
	J += Jco(i,0);
    }

    // if (J <= 0) return mat;

    // lumped mass
    double m = 0, mb = 0.0;
    double mp = 0;
    if (ndm == 2) {
	m = rho*J*thk/6.0;
	if (kappa > 0) {
	    mb = rho*J*thk*27.0/120.0;
	    mp = J*thk/6.0/kappa;
	}
    } else if (ndm == 3) {
	m = rho*J*thk/24.0;
	if (kappa > 0) {
	    mb = rho*J*thk*256.0/5040.0;
	    mp = J*thk/24.0/kappa;
	}
    }

    // number of corner points
    int nump = ntags.Size()/2;

    // corner
    for (int a=0; a<=nump; ++a) {

	if (kappa<=0 && a==nump) break;

	int dof = dofs(2*a);

	// velocity
	for (int i=0; i<ndm; ++i) {
	    if (a < nump) {
		mat(dof+i,dof+i) = m;
	    } else {
		mat(dof+i,dof+i) = mb;
	    }
	}

	// pressure
	if (a < nump && kappa > 0) {
	    int dofp = dofs(2*a+1);
	    mat(dofp,dofp) = mp;
	}
    }

    return mat;
}


const Matrix &
MINI::getDamp()
{
    int ndf = getNumDOF();
    int ndm = OPS_GetNDM();
    mat.resize(ndf,ndf);
    mat.Zero();

    double rho = data(0);
    double mu = data(1);
    double thk = data(2);
    double kappa = data(3);
    double J = 0.0;
    for (int i=0; i<Jco.noRows(); ++i) {
	J += Jco(i,0);
    }

    if (J <= 0) return mat;

    // number of corner points
    int nump = ntags.Size()/2;

    // sum(b_i^2), sum(c_i^2), sum(b_i*c_i)
    double b2 = 0.0, c2 =0.0, d2=0.0, bc2=0.0, bd2=0.0, cd2=0.0;
    if (mu > 0) {
	for (int i=0; i<Jco.noCols(); ++i) {
	    b2 = Jco(i,1)*Jco(i,1);
	    c2 = Jco(i,1)*Jco(i,2);
	    bc2 = Jco(i,1)*Jco(i,2);
	    if (Jco.noRows() > 3) {
		d2 = Jco(i,3)*Jco(i,3);
		bd2 = Jco(i,1)*Jco(i,3);
		cd2 = Jco(i,2)*Jco(i,3);
	    }
	}
    }

    // Gt*inv(Mb)*G
    Matrix Gb(ndm, nump), L(nump,nump);
    double alpha=0, beta=0;
    if (ndm == 2) {
	alpha = thk/6.0;
	beta = -thk*27.0/120.0;
    } else if (ndm == 3) {
	alpha = thk/24.0;
	beta = -thk*256.0/5040.0;
    }
    for (int i=0; i<ndm; ++i) {
	for (int b=0; b<nump; ++b) {
	    Gb(i, b) = -beta*Jco(b,i+1);
	}
    }
    if (kappa <= 0) {
	double invMb = 120*ops_Dt/(rho*J*thk*27.0);
	L.addMatrixTransposeProduct(0.0, Gb,Gb,invMb);
    }

    // K
    double alpha2D = mu*thk/(6*J);
    double alpha3D = mu*thk/(18*J);
    double beta2D = mu*thk/(1080*J);
    double beta3D = mu*thk/(272160*J);
    for (int a=0; a<=nump; ++a) {
	if (kappa <= 0 && a==nump) break;
	for (int b=0; b<=nump; ++b) {
	    if (kappa <= 0 && b==nump) break;
	    // dofs
	    int dofa = dofs(2*a);
	    int dofb = dofs(2*b);

	    // Kab
	    if (a < nump && b < nump && mu > 0) {

		// corner
		if (ndm == 2) {
		    mat(dofa,dofb) = alpha2D*(4*Jco(a,1)*Jco(b,1)
					      +3*Jco(a,2)*Jco(b,2));
		    mat(dofa,dofb+1) = alpha2D*(3*Jco(a,2)*Jco(b,1)
						-2*Jco(a,1)*Jco(b,2));
		    mat(dofa+1,dofb) = alpha2D*(3*Jco(a,1)*Jco(b,2)
						-2*Jco(a,2)*Jco(b,1));
		    mat(dofa+1,dofb+1) = alpha2D*(4*Jco(a,2)*Jco(b,2)
						  +3*Jco(a,1)*Jco(b,1));
		} else if (ndm == 3) {
		    mat(dofa,dofb) = alpha3D*(4*Jco(a,1)*Jco(b,1)
					       +3*Jco(a,2)*Jco(b,2)
					       +3*Jco(a,3)*Jco(b,3));
		    mat(dofa+1,dofb+1) = alpha3D*(4*Jco(a,2)*Jco(b,2)
						  +3*Jco(a,1)*Jco(b,1)
						  +3*Jco(a,3)*Jco(b,3));
		    mat(dofa+2,dofb+2) = alpha3D*(4*Jco(a,3)*Jco(b,3)
						  +3*Jco(a,2)*Jco(b,2)
						  +3*Jco(a,1)*Jco(b,1));
		    mat(dofa,dofb+1) = alpha3D*(3*Jco(a,2)*Jco(b,1)
						-2*Jco(a,1)*Jco(b,2));
		    mat(dofa+1,dofb) = alpha3D*(3*Jco(a,1)*Jco(b,2)
						-2*Jco(a,2)*Jco(b,1));
		    mat(dofa,dofb+2) = alpha3D*(3*Jco(a,3)*Jco(b,1)
						-2*Jco(a,1)*Jco(b,3));
		    mat(dofa+1,dofb+2) = alpha3D*(3*Jco(a,3)*Jco(b,2)
						  -2*Jco(a,2)*Jco(b,3));
		    mat(dofa+2,dofb) = alpha3D*(3*Jco(a,1)*Jco(b,3)
					      -2*Jco(a,3)*Jco(b,1));
		    mat(dofa+2,dofb+1) = alpha3D*(3*Jco(a,2)*Jco(b,3)
						  -2*Jco(a,3)*Jco(b,2));
		}


	    } else if (a == nump && b == nump && mu > 0) {
		
		// bubble
		if (ndm == 2) {
		    
		    mat(dofa,dofb) = beta2D*(4*b2+3*c2);
		    mat(dofa,dofb+1) = beta2D*bc2;
		    mat(dofa+1,dofb) = beta2D*bc2;
		    mat(dofa+1,dofb+1) = beta2D*(4*c2+3*b2);
		    
		} else if (ndm == 3) {
		    mat(dofa,dofb) = beta3D*(4*b2+3*c2+3*d2);
		    mat(dofa+1,dofb+1) = beta3D*(4*c2+3*b2+3*d2);
		    mat(dofa+2,dofb+2) = beta3D*(4*d2+3*c2+3*b2);
		    mat(dofa,dofb+1) = beta3D*bc2;
		    mat(dofa,dofb+2) = beta3D*bd2;
		    mat(dofa+1,dofb) = beta3D*bc2;
		    mat(dofa+1,dofb+2) = beta3D*cd2;
		    mat(dofa+2,dofb) = beta3D*bd2;
		    mat(dofa+2,dofb+1) = beta3D*cd2;
		}
	    }
	    
	    // Gab
	    if (b < nump) {
		int dofp = dofs(2*b+1);
		
		//-G and Gt
		for (int i=0; i<ndm; ++i) {

		    if (a < nump) {
			mat(dofa+i, dofp) = -alpha*Jco(a,i+1);
			mat(dofp, dofa+i) = alpha*Jco(a,i+1);
		    } else if (a == nump) {
			mat(dofa+i, dofp) = -beta*Jco(b,i+1);
			mat(dofp, dofa+i) = beta*Jco(b,i+1);
		    }
		}
	    }

	    // Gt*inv(Mb)*G
	    if (a < nump && b < nump && kappa <= 0) {
		int dofpa = dofs(2*a+1);
		int dofpb = dofs(2*b+1);
		mat(dofpa,dofpb) = L(a,b);
	    }
	    
	}
    }

    return mat;
}

// methods for applying loads
int
MINI::addInertiaLoadToUnbalance(const Vector &accel)
{
    return 0;
}

// methods for obtaining resisting force (force includes elemental loads)
const Vector &
MINI::getResistingForce()
{
    int ndf = getNumDOF();
    vec.resize(ndf);
    vec.Zero();
    return vec;
}

const Vector &
MINI::getResistingForceIncInertia()
{
    int ndf = getNumDOF();
    vec.resize(ndf);
    vec.Zero();

    int ndm = OPS_GetNDM();

    double rho = data(0);
    double thk = data(2);
    double kappa = data(3);
    double J = 0.0;
    for (int i=0; i<Jco.noRows(); ++i) {
	J += Jco(i,0);
    }

    if (J <= 0) return vec;

    // external force
    double f = 0;
    // double fb = 0;
    if (ndm == 2) {
	f = rho*J*thk/6.0;
	//fb = rho*J*thk*27.0/120.0;
    } else if (ndm == 3) {
	f = rho*J*thk/24.0;
	//fb = rho*J*thk*256.0/5040.0;
    }

    // number of corner points
    int nump = ntags.Size()/2;

    // -F
    for (int a=0; a<nump; ++a) {

	int dof = dofs(2*a);

	for (int i=0; i<ndm; ++i) {
	    // if (a < nump) {
	    vec(dof+i) = -f*data(4+i);
	    // } else if (a == nump) {
	    //vec(dof+i) = -fb*data(4+i);
	    //}
	}
    }

    // vdot, v
    Vector vdot(ndf), v(ndf);
    for (int a=0; a<=nump; ++a) {
	if (kappa <= 0 && a==nump) break;
	const Vector& vel = nodes[2*a]->getTrialVel();
	const Vector& accel = nodes[2*a]->getTrialAccel();

	int dof = dofs(2*a);

	for (int i=0; i<ndm; ++i) {
	    vdot(dof+i) = accel(i);
	    v(dof+i) = vel(i);
	}

	if (a < nump) {
	    int dofp = dofs(2*a+1);
	    const Vector& p = nodes[2*a+1]->getTrialVel();
	    const Vector& pdot = nodes[2*a+1]->getTrialAccel();
	    vdot(dofp) = pdot(0);
	    v(dofp) = p(0);
	}
    }

    // -r = M*vdot+K*v-F
    vec.addMatrixVector(1.0, this->getMass(), vdot, 1.0);
    vec.addMatrixVector(1.0, this->getDamp(), v, 1.0);

    return vec;
}

// MovableObject
const char *
MINI::getClassType() const
{
    return "MINI";
}

int
MINI::sendSelf(int commitTag, Channel &theChannel)
{
    return 0;
}

int
MINI::recvSelf(int commitTag, Channel &theChannel,
	       FEM_ObjectBroker &theBroker)
{
    return 0;
}

// DomainComponent
void
MINI::setDomain(Domain *theDomain)
{

    // set domain
    this->DomainComponent::setDomain(theDomain);
    if(theDomain == 0) return;

    // ndm
    int ndm = OPS_GetNDM();

    // ndf
    int ndf = 0;
    int eletag = this->getTag();
    Vector bcrds(ndm);

    // get corner nodes
    int nump = ntags.Size() / 2;
    for (int i=0; i<nump; ++i) {
	
	// get velocity node
    	nodes[2*i] = theDomain->getNode(ntags(2*i));
    	if (nodes[2*i] == 0) {
    	    opserr << "WARNING: node "<<ntags(2*i);
    	    opserr << "does not exist -- MINI\n";
    	    return;
    	}
    	if (nodes[2*i]->getNumberDOF() < ndm) {
    	    opserr << "WARNING: node "<<ntags(2*i);
    	    opserr << "ndf < ndm -- MINI\n";
    	    return;
    	}
    	if (nodes[2*i]->getCrds().Size() < ndm) {
    	    opserr << "WARNING: node "<<ntags(2*i);
    	    opserr << "crds.size() < ndm -- MINI\n";
    	    return;
    	}

	// center node
	const Vector& crds = nodes[2*i]->getCrds();
	const Vector& disp = nodes[2*i]->getTrialDisp();
	for (int j=0; j<ndm; j++) {
	    bcrds(j) += crds(j) + disp(j);
	}

	// ndf
	dofs(2*i) = ndf;
	ndf += nodes[2*i]->getNumberDOF();

	// get pc
	int pndf = 1;
        Pressure_Constraint* pc = theDomain->getPressure_Constraint(ntags(2*i));
        if(pc != 0) {
            pc->setDomain(theDomain);
        } else {
            pc = new Pressure_Constraint(ntags(2*i), pndf);
            if(pc == 0) {
                opserr<<"WARNING: no enough memory for PC -- ";
                opserr<<"MINI::setDomain "<<eletag<<"\n";
                return;
            }
            if(theDomain->addPressure_Constraint(pc) == false) {
                opserr<<"WARNING: failed to add PC to domain -- ";
                opserr<<"MINI::setDomain "<<eletag<<"\n";
                delete pc;
                return;
            }
        }

	// connect
        pc->connect(eletag);

        // get pressure node
        nodes[2*i+1] = pc->getPressureNode();
        if(nodes[2*i+1] == 0) {
            opserr<<"WARNING: pressure node does not exist ";
            opserr<<"in MINI - setDomain() "<<eletag<<"\n ";
            return;
        }
        ntags(2*i+1) = nodes[2*i+1]->getTag();
	dofs(2*i+1) = ndf;
        ndf += nodes[2*i+1]->getNumberDOF();
    }

    dofs(2*nump) = ndf;

    if (data(3) > 0) {
	// create bubble node
	NodeIter& theNodes = theDomain->getNodes();
	Node* theNode = theNodes();
	int bidx = ntags.Size() - 1;
	ntags(bidx) = 0;
	if (theNode != 0) {
	    ntags(bidx) = theNode->getTag();
	}
	ntags(bidx)--;
	
	bcrds /= nump;
	if (ndm == 2) {
	    nodes[bidx] = new Node(ntags(bidx),ndm,bcrds(0),bcrds(1));
	} else if (ndm == 3) {
	    nodes[bidx] = new Node(ntags(bidx),ndm,bcrds(0),bcrds(1),bcrds(2));
	}
	if (nodes[bidx] == 0) {
	    opserr<<"WARNING: run out of memory in creating node\n";
	    return;
	}
	if (theDomain->addNode(nodes[bidx]) == false) {
	    opserr<<"WARNING: failed to add node to domain\n";
	    delete nodes[bidx];
	    nodes[bidx] = 0;
	}
	
	ndf += nodes[bidx]->getNumberDOF();
	dofs(dofs.Size()-1) = ndf;
    }
}

// TaggedObject
void
MINI::Print(OPS_Stream &s, int flag)
{
    s << this->getClassType() <<"\n";
    s << "tag : " << this->getTag() << "\n";
    s << "nodes : ";
    for (int i=0; i<ntags.Size(); ++i) {
	s << ntags(i)<<" ";
    }
    s << "\n";
    s << "rho : " << data(0) << "\n";
    s << "mu : " << data(1) << "\n";
    s << "body force: ";
    for (int i=4; i<data.Size(); ++i) {
	s << data(i) << " ";
    }
    s << "\n";
    s << "thickness : " << data(2) << "\n";
    s << "kappa : " << data(3) << "\n";
}

int
MINI::displaySelf(Renderer &, int mode, float fact,
		  const char **displayModes,
		  int numModes)
{
    return 0;
}

double
MINI::det(const Matrix& M)
{
    if (M.noCols() != M.noRows()) {
	return 0.0;
    }

    // number of dimensions
    int num = M.noCols();

    // 2D
    if (num == 2) {
	return M(0,0)*M(1,1)-M(0,1)*M(1,0);
    }

    // 3D
    double res = 0.0;
    for (int j=0; j<num; ++j) {

	
	double res1 = 1.0;
	double res2 = 1.0;
	for (int i=0; i<num; ++i) {

	    // main diagonal
	    int m = i;
	    if (m >= num) m-=num;
	    int n = j+i;
	    if (n >= num) n-=num;
	    res1 *= M(m,n);

	    // off diagonal
	    m = i;
	    if (m >= num) m-=num;
	    n = j-i;
	    if (n < 0) n+=num;
	    res2 *= M(m,n);
	}

	res += res1;
	res -= res2;
	
    }

    return res;
    
}
