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
// $Date: 2012-08-31 11:36:36 $
// $Source: /usr/local/cvs/OpenSees/SRC/system_of_eqn/linearSOE/sparseGEN/PFEMLinSOE.h,v $

// Written: Minjie Zhu
// Created: August 2012
//
// Description: This file contains the class definition for PFEMLinSOE
// PFEMLinSOE is a subclass of SparseGenColLinSOE. It stores the matrix equation
// Ax=b using the sparse column-compacted storage scheme for storing the 
// matrix A. It solves the equations using the Fractional Step Method in PFEM.
//
// What: "@(#) PFEMLinSOE.h, revA"



#include <stdlib.h>
#include <PFEMLinSOE.h>
#include <PFEMSolver.h>
#include <Matrix.h>
#include <Graph.h>
#include <Vertex.h>
#include <VertexIter.h>
#include <math.h>
#include <Element.h>
#include <ElementIter.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <iostream>
using std::nothrow;
#include <Pressure_Constraint.h>
#include <Pressure_ConstraintIter.h>
#include <Domain.h>
#include <Node.h>
#include <NodeIter.h>
#include <DOF_Group.h>
#include <AnalysisModel.h>
#include <BackgroundMesh.h>
#include <elementAPI.h>
#ifdef _PARALLEL_INTERPRETERS
#include <mpi.h>
#endif

PFEMLinSOE::PFEMLinSOE(PFEMSolver &the_Solver)
    :LinearSOE(the_Solver, LinSOE_TAGS_PFEMLinSOE),
     M(0), Gft(0), Git(0), L(0), Qt(0),
     X(), B(), Mhat(), Mf(),
     dofType(), dofID(), assemblyFlag(0), stage(0)
{
    the_Solver.setLinearSOE(*this);
}


PFEMLinSOE::PFEMLinSOE()
    :LinearSOE(LinSOE_TAGS_PFEMLinSOE),
     M(0), Gft(0), Git(0), L(0), Qt(0),
     X(), B(), Mhat(), Mf(),
     dofType(), dofID(), assemblyFlag(0), stage(0)
{

}

PFEMLinSOE::PFEMLinSOE(int classTag)
    :LinearSOE(classTag),
     M(0), Gft(0), Git(0), L(0), Qt(0),
     X(), B(), Mhat(), Mf(),
     dofType(), dofID(), assemblyFlag(0), stage(0)
{

}


PFEMLinSOE::PFEMLinSOE(PFEMSolver &the_Solver, int classTag)
    :LinearSOE(the_Solver, classTag),
     M(0), Gft(0), Git(0), L(0), Qt(0),
     X(), B(), Mhat(), Mf(),
     dofType(), dofID(), assemblyFlag(0), stage(0)
{

}



    
PFEMLinSOE::~PFEMLinSOE()
{
    if(M != 0) cs_spfree(M);
    if(Gft != 0) cs_spfree(Gft);
    if(Git != 0) cs_spfree(Git);
    if(L != 0) cs_spfree(L);
    if(Qt != 0) cs_spfree(Qt);
}

int
PFEMLinSOE::solve(void)
{
    int res = LinearSOE::solve();
    if (res < 0) {
        return res;
    }

    assemblyFlag = 1;
    return 0;
}

int
PFEMLinSOE::getNumEqn(void) const
{
    return X.Size();
}

int 
PFEMLinSOE::setSize(Graph &theGraph)
{
    int result = 0;
    int size = theGraph.getNumVertex();
    if (size <= 0) {
	opserr << "WARNING: size<=0 -- PFEMLinSOE::setSize\n";
	return -1;
    }

    // resize vectors
    B.resize(size);
    X.resize(size);

    // zero the vectors
    B.Zero();
    X.Zero();

    // set Dof IDs
    int Ssize, Fsize, Isize, Psize, Pisize;
    result = this->setDofIDs(size, Ssize, Fsize, Isize, Psize, Pisize);

    // set matrix IDs
    result = this->setMatIDs(theGraph, Ssize, Fsize, Isize, Psize, Pisize);

    // reset flags
    assemblyFlag = 0;

    BackgroundMesh& bgmesh = OPS_getBgMesh();
    bool pressureonce = bgmesh.isPressureOnce();
    stage = 0;
    if (pressureonce) {
        stage = 1;
    }
    
    // invoke setSize() on the Solver    
    LinearSOESolver *the_Solver = this->getSolver();
    int solverOK = the_Solver->setSize();
    if (solverOK < 0) {
	opserr << "WARNING:PFEMLinSOE::setSize :";
	opserr << " solver failed setSize()\n";
	return solverOK;
    }
    return result;
}

int 
PFEMLinSOE::addA(const Matrix &m, const ID &id, double fact)
{
    // check for a quick return 
    if (fact == 0.0)  
	return 0;

    int idSize = id.Size();
    int size = X.Size();
    
    // check that m and id are of similar size
    if (idSize != m.noRows() && idSize != m.noCols()) {
        opserr << "PFEMLinSOE::addA() ";
        opserr << " - Matrix and ID not of similar sizes\n";
        return -1;
    }

    bool hasFluid = !skipFluid();

    int Ssize = M->n - Git->n;

    if (fact == 1.0) { // do not need to multiply 
        for (int i=0; i<idSize; i++) {
            int col = id(i);
            if(col>=size || col<0) continue;
            int coltype = dofType(col);         // column type
            int colid = dofID(col);             // column id

            if(coltype == 4) {                  // diganol of Mhat
                if (hasFluid) Mhat(colid) += m(i,i);
            } else if(coltype == 1) {           // diganol of Mf
                if (hasFluid) Mf(colid) += m(i,i);
            }

            if(coltype==4 || coltype<0) continue;

            for (int j=0; j<idSize; j++) {
                int row = id(j);
                if(row>=size || row<0) continue;
                cs* mat = 0;

                int rowtype = dofType(row);     // row type
                int rowid = dofID(row);         // row id
                int cid = colid;

                // get right matrix
                if(rowtype==0 && coltype==0) {                 // Ms
                    mat = M;
                } else if(rowtype==2 && coltype==2) {          // Ms
                    mat = M;
                    rowid += Ssize;
                    cid += Ssize;
                } else if(rowtype==0 && coltype==2) {          // Msi
                    mat = M;
                    cid += Ssize;
                } else if(rowtype==2 && coltype==0) {          // Mis
                    mat = M;
                    rowid += Ssize;
                } else if(rowtype==3 && coltype==1) {          // Gft
                    if (hasFluid) mat = Gft;
                } else if(rowtype==3 && coltype==2) {          // Git
                    mat = Git;
                } else if(rowtype==3 && coltype==3) {          // L
                    if (hasFluid) mat = L;
                } else if(rowtype==4 && coltype==3) {          // Qt
                    if (hasFluid) mat = Qt;
                }

                if(mat == 0) continue;

                // find place in mat
                for(int k=mat->p[cid]; k<mat->p[cid+1]; k++) {
                    if(mat->i[k] == rowid) {
                        mat->x[k] += m(j,i);
                        break;
                    }
                }
            }  // for j		
        }  // for i
    } else {
        for (int i=0; i<idSize; i++) {
            int col = id(i);
            if(col>=size || col<0) continue;
            int coltype = dofType(col);         // column type
            int colid = dofID(col);             // column id

            if(coltype == 4) {                  // diganol of Mhat
                if (hasFluid) Mhat(colid) += fact*m(i,i);
            } else if(coltype == 1) {           // diganol of Mf
                if (hasFluid) Mf(colid) += fact*m(i,i);
            }

            if(coltype==4 || coltype<0) continue;

            for (int j=0; j<idSize; j++) {
                int row = id(j);
                if(row>=size || row<0) continue;
                cs* mat = 0;

                int rowtype = dofType(row);     // row type
                int rowid = dofID(row);         // row id
                int cid = colid;

                // get right matrix
                if(rowtype==0 && coltype==0) {                 // Ms
                    mat = M;
                } else if(rowtype==2 && coltype==2) {          // Ms
                    mat = M;
                    rowid += Ssize;
                    cid += Ssize;
                } else if(rowtype==0 && coltype==2) {          // Msi
                    mat = M;
                    cid += Ssize;
                } else if(rowtype==2 && coltype==0) {          // Mis
                    mat = M;
                    rowid += Ssize;
                } else if(rowtype==3 && coltype==1) {          // Gft
                    if (hasFluid) mat = Gft;
                } else if(rowtype==3 && coltype==2) {          // Git
                    mat = Git;
                } else if(rowtype==3 && coltype==3) {          // L
                    if (hasFluid) mat = L;
                } else if(rowtype==4 && coltype==3) {          // Qt
                    if (hasFluid) mat = Qt;
                }

                if(mat == 0) continue;

                // find place in mat
                for(int k=mat->p[cid]; k<mat->p[cid+1]; k++) {
                    if(mat->i[k] == rowid) {
                        mat->x[k] += fact*m(j,i);
                        break;
                    }
                }
            }  // for j		
        }  // for i
    }
    return 0;
}

    
int 
PFEMLinSOE::addB(const Vector &v, const ID &id, double fact)
{
    // check for a quick return 
    if (fact == 0.0)  return 0;

    int idSize = id.Size();    
    int size = X.Size();

    // check that m and id are of similar size
    if (idSize != v.Size() ) {
	opserr << "PFEMLinSOE::addB() ";
	opserr << " - Vector and ID not of similar sizes\n";
	return -1;
    }    

    if (fact == 1.0) { // do not need to multiply if fact == 1.0
	for (int i=0; i<idSize; i++) {
	    int pos = id(i);
	    if (pos <size && pos >= 0)
		B(pos) += v(i);
	}
    } else if (fact == -1.0) { // do not need to multiply if fact == -1.0
	for (int i=0; i<idSize; i++) {
	    int pos = id(i);
	    if (pos <size && pos >= 0)
		B(pos) -= v(i);
	}
    } else {
	for (int i=0; i<idSize; i++) {
	    int pos = id(i);
	    if (pos <size && pos >= 0)
		B(pos) += v(i) * fact;
	}
    }	

    return 0;
}


int
PFEMLinSOE::setB(const Vector &v, double fact)
{
    // check for a quick return 
    if (fact == 0.0)  return 0;

    if (v.Size() != B.Size()) {
	opserr << "WARNING BandGenLinSOE::setB() -";
	opserr << " incompatible sizes " << B.Size() << " and " << v.Size() << endln;
	return -1;
    }
    
    if (fact == 1.0) { // do not need to multiply if fact == 1.0
        B.Zero();
        B += v;
    } else if (fact == -1.0) {
        B.Zero();
        B -= v;
    } else {
        B.Zero();
        B += v;
        B *= fact;
    }	
    return 0;
}

void 
PFEMLinSOE::zeroA(void)
{
    for (int i = 0; i < M->nzmax; i++)
        M->x[i] = 0.0;

    for (int i = 0; i < Git->nzmax; i++)
        Git->x[i] = 0.0;

    // zero fluid part
    bool hasFluid = !skipFluid();
    if (hasFluid) {
        for (int i = 0; i < Gft->nzmax; i++)
            Gft->x[i] = 0.0;
        for (int i = 0; i < L->nzmax; i++)
            L->x[i] = 0.0;
        for (int i = 0; i < Qt->nzmax; i++)
            Qt->x[i] = 0.0;

        Mhat.Zero();
        Mf.Zero();
    }
}
	
void 
PFEMLinSOE::zeroB(void)
{
    B.Zero();
}

void 
PFEMLinSOE::setX(int loc, double value)
{
    if (loc < X.Size() && loc >=0)
	X(loc) = value;
}

void 
PFEMLinSOE::setX(const Vector &x)
{
    X.Zero();
    X += x;

}

const Vector &
PFEMLinSOE::getX(void)
{
    return X;
}

const Vector &
PFEMLinSOE::getB(void)
{
    return B;
}

double 
PFEMLinSOE::normRHS(void)
{
    return B.Norm();
}    


int
PFEMLinSOE::setPFEMSolver(PFEMSolver &newSolver)
{
    newSolver.setLinearSOE(*this);
    
    if (X.Size() != 0) {
	int solverOK = newSolver.setSize();
	if (solverOK < 0) {
	    opserr << "WARNING:PFEMLinSOE::setSolver :";
	    opserr << "the new solver could not setSeize() - staying with old\n";
	    return -1;
	}
    }
    
    return this->LinearSOE::setSolver(newSolver);
}


int 
PFEMLinSOE::sendSelf(int cTag, Channel &theChannel)
{
    return 0;
}

int 
PFEMLinSOE::recvSelf(int cTag, Channel &theChannel, 
                     FEM_ObjectBroker &theBroker)  
{
    return 0;
}

int 
PFEMLinSOE::setDofIDs(int size,int& Ssize, int&Fsize, int& Isize,int& Psize,int& Pisize)
{
    if(theModel == 0) {
        opserr << "Analysis model has not been linked - PFEMLinSOE::setDofIDs()\n";
        return -1;
    }

    dofType.resize(size);
    dofID.resize(size);
    dofType.Zero();
    dofID.Zero();
    Domain* domain = theModel->getDomainPtr();
    if(domain == 0) {
        opserr<<"WARNING: no domain is set for the model";
        opserr<<" -- PFEMLinSOE::setDofIDs\n";
        return -1;
    }

    // loop through pcs
    Pressure_ConstraintIter& thePCs = domain->getPCs();
    Pressure_Constraint* thePC = 0;
    Ssize = Fsize = Isize = Psize = Pisize = 0;
    int Isosize = 0;
    while((thePC = thePCs()) != 0) {
        // int ptag = thePC->getPressureNode();
        int ntag = thePC->getTag();
        // Node* pnode = domain->getNode(ptag);
        Node* pnode = thePC->getPressureNode();
        // if(pnode == 0) {
        //     opserr<<"WARNING: pressure node "<<ptag<<" does not exists ";
        //     opserr<<" -- PFEMLinSOE::setDofIDs\n";
        //     return -1;
        // }
        Node* nnode = domain->getNode(ntag);
        if(nnode == 0) {
            opserr<<"WARNING: pressure node "<<ntag<<" does not exists ";
            opserr<<" -- PFEMLinSOE::setDofIDs\n";
            return -1;
        }
        // DOF_Group* pDOF = pnode->getDOF_GroupPtr();
        DOF_Group* pDOF = 0;
        if(pnode != 0) {
            int ptag = pnode->getTag();
            pDOF = pnode->getDOF_GroupPtr();
            if(pDOF == 0) {
                opserr<<"WARNING: no DOF_Group is found with node "<<ptag;
                opserr<<" -- PFEMLinSOE::setDofIDs\n";
                return -1;
            }
        }
        // if(pDOF == 0) {
        //     opserr<<"WARNING: no DOF_Group is found with node "<<ptag;
        //     opserr<<" -- PFEMLinSOE::setDofIDs\n";
        //     return -1;
        // }
        DOF_Group* nDOF = nnode->getDOF_GroupPtr();
        if(nDOF == 0) {
            opserr<<"WARNING: no DOF_Group is found with node "<<ntag;
            opserr<<" -- PFEMLinSOE::setDofIDs\n";
            return -1;
        }
        // const ID& pid = pDOF->getID();
        const ID& nid = nDOF->getID();

        // pressure nodes
        if(pnode != 0) {
            const ID& pid = pDOF->getID();

	    if(thePC->isFreeSurf() && thePC->isFluid()) {
		for(int i=0; i<pid.Size(); i++) {
                    if(pid(i) >= 0) {
                        dofType(pid(i)) = -1;    
                        dofID(pid(i)) = -1;
                    }
                }
	    } else if(thePC->isFluid() || thePC->isInterface()) {
                if(pid(0) >= 0) {
                    dofType(pid(0)) = 3;          // pressure
                    dofID(pid(0)) = Psize++;
                }
                for(int i=1; i<pid.Size(); i++) {
                    if(pid(i) >= 0) {
                        dofType(pid(i)) = 4;      // pressure gradient
                        dofID(pid(i)) = Pisize++;     
                    }
                }
            } else {
                for(int i=0; i<pid.Size(); i++) {
                    if(pid(i) >= 0) {
                        dofType(pid(i)) = -1;     // which are not connected to any PFEM elements
                        dofID(pid(i)) = -1;
                    }
                }

            }
        }

        // momentum nodes
        int ndm = nnode->getCrds().Size();
        if(thePC->isInterface()) {
            for(int i=0; i<ndm; i++) {
                if(nid(i) >= 0) {
                    dofType(nid(i)) = 2;      // interface momentum 
                    dofID(nid(i)) = Isize++;
                } else {
                }
            }
        } else if(thePC->isFluid()) {
            for(int i=0; i<ndm; i++) {
                if(nid(i) >= 0) {
                    dofType(nid(i)) = 1;      // fluid momentum 
                    dofID(nid(i)) = Fsize++;
                } else {
                }
            }
        } else if(thePC->isIsolated()) {
            for(int i=0; i<nid.Size(); i++) {
                if(nid(i) >= 0) {
                    dofType(nid(i)) = -1;     // which are not connected to any elements 
                    dofID(nid(i)) = -1;
                }
		Isosize++;
            }
        }
    }
    
    for(int col=0; col<size; col++) {         
        if(dofType(col) == 0) {               // structure momentum 
            dofID(col) = Ssize++;
        }
    }

#ifdef _PARALLEL_INTERPRETERS
    int myid;
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    if(myid == 0) {
	opserr<<"Ssize = "<<Ssize<<", ";
	opserr<<"Fsize = "<<Fsize<<", ";
	opserr<<"Isize = "<<Isize<<", ";
	opserr<<"Psize = "<<Psize<<", ";
	opserr<<"Isosize = "<<Isosize<<"\n";
    }

#else
    opserr<<"Ssize = "<<Ssize<<", ";
    opserr<<"Fsize = "<<Fsize<<", ";
    opserr<<"Isize = "<<Isize<<", ";
    opserr<<"Psize = "<<Psize<<", ";
    opserr<<"Isosize = "<<Isosize<<"\n";
#endif
	
    return 0;
}

int
PFEMLinSOE::setMatIDs(Graph& theGraph, int Ssize, int Fsize, int Isize, int Psize, int Pisize)
{
    cs* M1 = cs_spalloc(Ssize+Isize, Ssize+Isize, 1, 1, 1);
    cs* Gft1 = cs_spalloc(Psize, Fsize, 1, 1, 1);
    cs* Git1 = cs_spalloc(Psize, Isize, 1, 1, 1);
    cs* L1 = cs_spalloc(Psize, Psize, 1, 1, 1);
    cs* Qt1 = cs_spalloc(Pisize, Psize, 1, 1, 1);
    Mhat.resize(Pisize); Mhat.Zero();
    Mf.resize(Fsize); Mf.Zero();

    Vertex* theVertex = 0;
    int size = X.Size();
    for (int a=0; a<size; a++) {    // columns
        theVertex = theGraph.getVertexPtr(a);
        if (theVertex == 0) {
            opserr << "WARNING:PFEMLinSOE::setSize :";
            opserr << " vertex " << a << " not in graph!\n";
            break;
        }

        int col = theVertex->getTag();  // column 
        int coltype = dofType(col);     // column type
        int colid = dofID(col);         // column id
        if(coltype==4 || coltype<0) continue;      // don't need this column

        // diagnol terms
        if(coltype == 0) {                      // structure momentum
            cs_entry(M1, colid, colid, 0.0);   
        } else if(coltype == 2) {               // interface momentum
            cs_entry(M1, colid+Ssize, colid+Ssize, 0.0);   
        } else if(coltype == 3) {               // pressure
            cs_entry(L1, colid, colid, 0.0);    
        }

        // off diagnol terms
        const ID &theAdjacency = theVertex->getAdjacency();
        int idSize = theAdjacency.Size();
        for (int i=0; i<idSize; i++) {       // rows
            int row = theAdjacency(i);       // row 
            int rowtype = dofType(row);      // row type
            int rowid = dofID(row);          // row id
            
            if(rowtype==0 && coltype==0) {                // Ms
                cs_entry(M1, rowid, colid, 0.0);
            } else if(rowtype==2 && coltype==2) {         // Mi
                cs_entry(M1, rowid+Ssize, colid+Ssize, 0.0);
            } else if(rowtype==0 && coltype==2) {         // Msi
                cs_entry(M1, rowid, colid+Ssize, 0.0);
            } else if(rowtype==2 && coltype==0) {         // Mis
                cs_entry(M1, rowid+Ssize, colid, 0.0);
            } else if(rowtype==3 && coltype==1) {        // Gft
                cs_entry(Gft1, rowid, colid, 0.0);
            } else if(rowtype==3 && coltype==2) {        // Git
                cs_entry(Git1, rowid, colid, 0.0);
            } else if(rowtype==3 && coltype==3) {        // L
                cs_entry(L1, rowid, colid, 0.0);
            } else if(rowtype==4 && coltype==3) {        // Qt
                cs_entry(Qt1, rowid, colid, 0.0);
            }
        }
    }

    // convert to compressed format
    if(M != 0) cs_spfree(M);
    M = cs_compress(M1);
    cs_spfree(M1);

    if(Gft != 0) cs_spfree(Gft);
    Gft = cs_compress(Gft1);
    cs_spfree(Gft1);

    if(Git != 0) cs_spfree(Git);
    Git = cs_compress(Git1);
    cs_spfree(Git1);

    if(L != 0) cs_spfree(L);
    L = cs_compress(L1);
    cs_spfree(L1);

    if(Qt != 0) cs_spfree(Qt);
    Qt = cs_compress(Qt1);
    cs_spfree(Qt1);

    // reorder rows
    // cs* mats[5] = {M,Gft,Git,L,Qt};
    // for (int i=0; i<5; i++) {
    // 	cs* mat = mats[i];
    // 	for (int j=0; j<mat->n; j++) {
    // 	    ID col(0, mat->p[j+1]-mat->p[j]);
    // 	    for (int k=mat->p[j]; k<mat->p[j+1]; k++) {
    // 		col.insert(mat->i[k]);
    // 	    }
    // 	    int index = 0;
    // 	    for (int k=mat->p[j]; k<mat->p[j+1]; k++) {
    // 		mat->i[k] = col[index++];
    // 	    }
    // 	}
    // }

    return 0;
}

bool
PFEMLinSOE::isFluidID(const ID &id) const
{
    bool fluid = true;
    for (int i = 0; i < id.Size(); ++i) {
        if (dofType(id(i))==0 || dofType(id(i))==2) {
            fluid = false;
            break;
        }
    }

    return fluid;
}

bool
PFEMLinSOE::skipFluid() const
{
    BackgroundMesh& bgmesh = OPS_getBgMesh();
    return assemblyFlag==1 && bgmesh.isDispOn()==false && bgmesh.isFastAssembly();
}

void PFEMLinSOE::saveK(OPS_Stream& output) {
    if (M == 0) return;
    output << "sparse matrix <" << M->m << ", " << M->n << "> with "
           << M->nzmax << " entries\n";

    // save the matrix
    for (int j = 0; j < M->n; ++j) {
        for (int k = M->p[j]; k < M->p[j + 1]; ++k) {
            output << "    " << M->i[k] << "    " << j << "    ("
                   << M->x[k] << ")\n";
        }
    }
}