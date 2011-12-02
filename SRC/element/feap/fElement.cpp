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
// $Date: 2000-09-15 08:23:20 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/feap/fElement.cpp,v $
                                                                        
                                                                        
// File: ~/element/fortran/fElement.C
// 
// Written: fmk 
// Created: 07/98
// Revision: A
//
// Description: This file contains the implementation for the fElement class.
//
// What: "@(#) fElement.C, revA"

#include "fElement.h"
#include <Domain.h>
#include <Node.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <UniaxialMaterial.h>
#include <Renderer.h>

#include <math.h>
#include <stdlib.h>

// initialise all class wise pointers to 0 and numfElements to 0
Matrix **fElement::fElementM;
Vector **fElement::fElementV;
double *fElement::r;
double *fElement::s;
double *fElement::ul;
double *fElement::xl;
double *fElement::tl;
int    *fElement::ix;
int    fElement::numfElements(0);

#define MAX_NST 64

// constructor:
//  responsible for allocating the necessary space needed by each object
//  and storing the tags of the fElement end nodes.
fElement::fElement(int tag, 
		   int classTag,
		   int EleType,
		   int sizeD, int NEN,
		   int NDM, int NDF,
		   int numNh1, int numNh3)
:Element(tag,classTag), nh1(numNh1), nh3(numNh3), h(0), eleType(EleType),
  theNodes(0), u(0), nen(NEN), ndf(NDF), ndm(NDM), d(0), data(0), 
 connectedNodes(0),nrCount(0), theLoad(0)
  
{
    // allocate space for h array
    if (nh1 < 0) nh1 = 0;
    if (nh3 < 0) nh3 = 0;
    if (nh1 != 0 || nh3 != 0) {
	int sizeH = 2*nh1 + nh3;
	h = new double[sizeH];
	if (h == 0) {
	    cerr << "FATAL: fElement::fElement() - eleTag: " << tag;
	    cerr << " ran out of memory creating h of size " << 2*nh1+nh3 << endl;
	    exit(-1);
	}	    
	for (int i=0; i<sizeH; i++) h[i] = 0.0;
    }
    connectedNodes = new ID(NEN);
    d = new double[sizeD];
    for (int i=0; i<sizeD; i++) d[i] = 0.0;
    data = new Vector(d, sizeD);

    // allocate space for static varaibles on creation of first instance
    if (numfElements == 0) {
	fElementM = new Matrix *[MAX_NST+1];
	fElementV = new Vector *[MAX_NST+1];
	s = new double[(MAX_NST+1)*(MAX_NST+1)];
	r = new double[MAX_NST+1];
	ul = new double[(MAX_NST+1)*6];
	xl = new double[MAX_NST+1];
	tl = new double[MAX_NST+1];
	ix = new int[MAX_NST+1];	

	// check space was available -- otherwise exit
	if (fElementM == 0 || fElementV == 0 || ix == 0 ||
	    r == 0 || s == 0 || ul == 0 || xl == 0 || tl == 0) {

	    cerr << "FATAL: fElement::fElement() - eleTag: " << tag;
	    cerr << " ran out of memory initialising static stuff\n";
	    exit(-1);	    
	}
	
	for (int i=0; i<MAX_NST+1; i++) {
	    fElementM[i] = 0;
	    fElementV[i] = 0;
	}
        fElementM[0] = new Matrix(1,1); // dummy for error
	fElementV[0] = new Vector(1);
    }
    
    // increment number of elements
    numfElements++;
}


fElement::fElement(int tag, 
		   int classTag,
		   int EleType,
		   int sizeD, int NEN,
		   int NDM, int NDF, int iow)
:Element(tag,classTag), nh1(0), nh3(0), h(0), eleType(EleType),
  theNodes(0), u(0), nen(NEN), ndf(NDF), ndm(NDM), d(0), data(0), 
  connectedNodes(0), nrCount(0), theLoad(0)
{
    connectedNodes = new ID(NEN);
    d = new double[sizeD];
    data = new Vector(d, sizeD);
    if (d == 0 || data == 0) {
	cerr << "FATAL: fElement::fElement() - eleTag: " << tag;
	cerr << " ran out of memory creating d of size " << sizeD << endl;
	exit(-1);
    }	    
    for (int i=0; i<sizeD; i++) d[i] = 0.0;    

    // invoke the elmt() routine with isw == 1 to read in the element data
    // and create room for the h array stuff needed by the element
    this->invokefInit(1, iow); 
    // allocate space for h array
    if (nh1 < 0) nh1 = 0;
    if (nh3 < 0) nh3 = 0;
    if (nh1 != 0 || nh3 != 0) {
	int sizeH = 2*nh1+nh3;
	h = new double[sizeH];
	if (h == 0) {
	    cerr << "FATAL: fElement::fElement() - eleTag: " << this->getTag();
	    cerr << " ran out of memory creating h of size " << sizeH << endl;
	    exit(-1);
	}	    
	else
	    for (int i=0; i<sizeH; i++) h[i] = 0.0;
    }

    // allocate space for static varaibles on creation of first instance
    if (numfElements == 0) {
	fElementM = new Matrix *[MAX_NST+1];
	fElementV = new Vector *[MAX_NST+1];
	s = new double[(MAX_NST+1)*(MAX_NST+1)];
	r = new double[MAX_NST+1];
	ul = new double[(MAX_NST+1)*6];
	xl = new double[MAX_NST+1];
	tl = new double[MAX_NST+1];
	ix = new int[MAX_NST+1];	

	// check space was available -- otherwise exit
	if (fElementM == 0 || fElementV == 0 || ix == 0 ||
	    r == 0 || s == 0 || ul == 0 || xl == 0 || tl == 0) {

	    cerr << "FATAL: fElement::fElement() - eleTag: " << tag;
	    cerr << " ran out of memory initialising static stuff\n";
	    exit(-1);	    
	}
	
	for (int i=0; i<MAX_NST+1; i++) {
	    fElementM[i] = 0;
	    fElementV[i] = 0;
	}
        fElementM[0] = new Matrix(1,1); // dummy for error
	fElementV[0] = new Vector(1);
    }

    // increment number of elements
    numfElements++;
}

// constructor:
//   invoked by a FEM_ObjectBroker - blank object that recvSelf needs
//   to be invoked upon
fElement::fElement(int classTag)
:Element(0, classTag), nh1(0), nh3(0), h(0),
 theNodes(0), u(0), nen(0), ndf(0), ndm(0), d(0), data(0), connectedNodes(0),
 theLoad(0)
{
    // does nothing
}

//  destructor
//     delete must be invoked on any objects created by the object
//     and on the matertial object.
fElement::~fElement()
{
    // clear up any space allocated for individual object

    if (h != 0)
	delete [] h;
    if (u != 0)
	delete [] u;
    if (theNodes != 0)
	delete [] theNodes;

    if (data != 0)
	delete  data;
    if (connectedNodes != 0)
	delete  connectedNodes;
    if (d != 0)
	delete [] d;
    if (theLoad != 0)
      delete theLoad;

    // if last element - clear up space allocated

    numfElements --;    
    if (numfElements == 0) {
	for (int i=0; i<MAX_NST+1; i++) {
	    if (fElementM[i] != 0) delete fElementM[i];
	    if (fElementV[i] != 0) delete fElementV[i];
	}
	delete [] fElementM;
	delete [] fElementV;
	delete [] s;
	delete [] r;
	delete [] ul;
	delete [] xl;
	delete [] tl;
	delete [] ix;	
    }
}

int
fElement::getNumExternalNodes(void) const
{
    return connectedNodes->Size();
}

const ID &
fElement::getExternalNodes(void)
{
    return *connectedNodes;
}

int
fElement::getNumDOF(void)
{
    return ndf*nen;
}

void
fElement::setDomain(Domain *theDomain)
{
    // check Domain is not null - invoked when object removed from a domain
    if (theDomain == 0) {
	ndm = 0;
	nen = 0;
	ndf = 0;
	if (theNodes != 0) {
	    delete [] theNodes;
	    theNodes = 0;
	}
	return;
    }    

    // set up the pointers to the nodes
    const ID &theNodeTags = this->getExternalNodes();

    int numNodes = theNodeTags.Size();
    theNodes = new Node *[numNodes];
    for (int i=0; i<numNodes; i++) {
	Node *theNode = theDomain->getNode(theNodeTags(i));
	if (theNode == 0) {
	    cerr << "WARNING fElement::setDomain(Domain *theDomain) - node: ";
	    cerr << theNodeTags(i) << " does not exist in domain for ele " << *this;
	    ndm = 0; nen = 0; ndf = 0;
	    return;
	}
	// set up the pointer to the node
	theNodes[i] = theNode;

	
	// check the dimension and number of dof at the node 
	// is the same as all the other nodes of the element
	if (i == 0) {
	    const Vector &crds = theNode->getCrds();
	    ndm = crds.Size();	              // ndm = dimesion of mesh
	    ndf = theNode->getNumberDOF();    // ndf = number of dof at nodes
	} else {
	    const Vector &crds = theNode->getCrds();	
	    if (ndm != crds.Size()) {
		cerr << "WARNING fElement::setDomain(Domain *theDomain) - node: ";
		cerr << theNodeTags(i) << " not in correct dimension " << *this;
		ndm = 0; nen = 0; ndf = 0;
		return;		
	    }
	    if (ndf != theNode->getNumberDOF()) {
		cerr << "WARNING fElement::setDomain(Domain *theDomain) - node: ";
		cerr << theNodeTags(i) << " does not have correct #DOF " << *this;
		ndm = 0; nen = 0; ndf = 0;
		return;		
	    }
	}
    }

    
    // call the DomainComponent class method THIS IS VERY IMPORTANT
    this->DomainComponent::setDomain(theDomain);
    
    // set nen - the number of element nodes
    nen = numNodes;
    
    // allocate memory for u
    int nst = ndf*nen;    
    u = new double[nst]; 
    if (u == 0) {
	cerr << "WARNING fElement::setDomain(Domain *theDomain) -  ";
	cerr << " ran out of memory creating u of size: " << nen*ndf << *this;
	ndm = 0; nen = 0; ndf = 0;
	return;			
    }
    // zero u
    for (int ii=0; ii<nst; ii++) 
	u[ii] = 0.0;


    theLoad = new Vector(nst);
    if (theLoad == 0) {
	g3ErrorHandler->fatal("Truss::setDomain - truss %d %s %d\n",
			      this->getTag(), 
			      "out of memory creating vector of size",
			      nst);	
      return;
    }          

    // allocate the Matrix and Vector objects if none yet for this size nst
   if (fElementM[nst] == 0) {
       fElementM[nst] = new Matrix(s,nst,nst);
       fElementV[nst] = new Vector(r,nst);       

       if ((fElementM[nst] == 0) || (fElementV[nst] == 0)) {
	   cerr << "WARNING fElement::setDomain(Domain *theDomain) -  ";
	   cerr << " ran out of memory creating Matrix and Vector for " << *this;
	   ndm = 0; nen = 0; ndf = 0;
	   return;			
       }
   }
}


int
fElement::commitState()
{
    if (nh1 != 0) 
	for (int i=0; i<nh1; i++)
	    h[i] = h[i+nh1];

    nrCount = 0;
    return 0;
}

int
fElement::revertToLastCommit()
{
    if (nh1 != 0) 
	for (int i=0; i<nh1; i++)
	    h[i+nh1] = h[i];    

    nrCount = 0;
    return 0;
}

int
fElement::revertToStart()
{
    // does nothing
    return 0;
}


const Matrix &
fElement::getTangentStiff(void)
{

    // check for quick return
    if (nen == 0)
	return (*fElementM[0]);
    
    // get the current load factor
    Domain *theDomain=this->getDomain();
    double dm = theDomain->getCurrentTime();
    
    // set ctan, ior and iow
    double ctan[3];
    ctan[0] = 1.0; ctan[1] = 0.0; ctan[2] = 0.0;
    int ior = 0; int iow = 0;
    
    // call the ready routine to set ul, xl, tl and ix, NH1, NH2 and NH3
    int nstR = this->readyfRoutine(false);

    // zero the matrix
    fElementM[nstR]->Zero();    
    
    // invoke the fortran subroutine
    int isw = 3; 
    int nstI = this->invokefRoutine(ior, iow, ctan, isw);
    
    // check nst is as determined in readyfRoutine()
    if (nstI != nstR) {
	cerr << "FATAL fElement::getTangentStiff() problems with incompatable nst";
	cerr << " ready: " << nstR << " invoke: " << nstI << endl;
	exit(-1);
    }

    // return the matrix

    return *(fElementM[nstR]);

}


const Matrix &
fElement::getSecantStiff(void)
{
    // check for quick return
    if (nen == 0)
	return (*fElementM[0]);
    
    // get the current load factor
    Domain *theDomain=this->getDomain();
    double dm = theDomain->getCurrentTime();
    
    // set ctan, ior and iow
    double ctan[3];
    ctan[0] = 1.0; ctan[1] = 0.0; ctan[2] = 0.0;
    int ior = 0; int iow = 0;
    
    // call the ready routine to set ul, xl, tl and ix, NH1, NH2 and NH3
    int NH1, NH2, NH3;    
    int nstR = this->readyfRoutine(false);
    
    // zero the matrix
    fElementM[nstR]->Zero();    
    
    // invoke the fortran subroutine
    int isw = 3; 
    int nstI = this->invokefRoutine(ior, iow, ctan, isw);

    // check nst is as determined in readyfRoutine()
    if (nstI != nstR) {
	cerr << "FATAL fElement::getTangentStiff() problems with incompatable nst";
	cerr << " ready: " << nstR << " invoke: " << nstI << endl;
	exit(-1);
    }
    
    // return the matrix
    return *(fElementM[nstR]);
}
    
const Matrix &
fElement::getDamp(void)
{
    // check for quick return
    if (nen == 0)
	return (*fElementM[0]);
    
    // get the current load factor
    Domain *theDomain=this->getDomain();
    double dm = theDomain->getCurrentTime();
    
    // set ctan, ior and iow
    double ctan[3];
    ctan[0] = 0.0; ctan[1] = 1.0; ctan[2] = 0.0;
    int ior = 0; int iow = 0;
    
    // call the ready routine to set ul, xl, tl and ix, NH1, NH2 and NH3
    int NH1, NH2, NH3;    
    int nstR = this->readyfRoutine(true);
    
    // zero the matrix
    fElementM[nstR]->Zero();    
    
    // invoke the fortran subroutine
    int isw = 3; int nst = nen*ndf; int n = this->getTag();
    int nstI = this->invokefRoutine(ior, iow, ctan, isw);
    
    // check nst is as determined in readyfRoutine()
    if (nstI != nstR) {
	cerr << "FATAL fElement::getTangentStiff() problems with incompatable nst";
	cerr << " ready: " << nstR << " invoke: " << nstI << endl;
	exit(-1);
    }
    
    // return the matrix
    return *(fElementM[nstR]);
}


const Matrix &
fElement::getMass(void)
{
    // check for quick return
    if (nen == 0)
	return (*fElementM[0]);
    
    // get the current load factor
    Domain *theDomain=this->getDomain();
    double dm = theDomain->getCurrentTime();
    
    // set ctan, ior and iow
    double ctan[3];
    ctan[0] = 0.0; ctan[1] = 0.0; ctan[2] = 1.0;
    int ior = 0; int iow = 0;
    
    // call the ready routine to set ul, xl, tl and ix, NH1, NH2 and NH3
    int NH1, NH2, NH3;    
    int nstR = this->readyfRoutine(true);
    
    // zero the matrix and vector (consistant and lumped)
    fElementM[nstR]->Zero();    
    fElementV[nstR]->Zero();        
    
    // invoke the fortran subroutine
    int isw = 5; int nst = nen*ndf; int n = this->getTag();
    int nstI = this->invokefRoutine(ior, iow, ctan, isw);
    
    // check nst is as determined in readyfRoutine()
    if (nstI != nstR) {
	cerr << "FATAL fElement::getTangentStiff() problems with incompatable nst";
	cerr << " ready: " << nstR << " invoke: " << nstI << endl;
	exit(-1);
    }
    
    // return the matrix
    return *(fElementM[nstR]);
}



void 
fElement::zeroLoad(void)
{
  // does nothing now
  if (theLoad != 0)
    theLoad->Zero();
}

int
fElement::addLoad(const Vector &load)
{
  // does nothing now
  if (theLoad != 0)
    (*theLoad) += load;
	return 0;
}


const Vector &
fElement::getResistingForce()
{		
    // check for quick return
    if (nen == 0)
	return (*fElementV[0]);
    
    // get the current load factor
    Domain *theDomain=this->getDomain();
    double dm = theDomain->getCurrentTime();
    
    // set ctan, ior and iow
    double ctan[3];
    ctan[0] = 0.0; ctan[1] = 1.0; ctan[2] = 0.0;
    int ior = 0; int iow = 0;
    
    // call the ready routine to set ul, xl, tl and ix, NH1, NH2 and NH3
    int NH1, NH2, NH3;    
    int nstR = this->readyfRoutine(false);

    // zero the vector
    fElementV[nstR]->Zero();        
    
    // invoke the fortran subroutine
    int isw = 6; int nst = nen*ndf; int n = this->getTag();
    int nstI = this->invokefRoutine(ior, iow, ctan, isw);
    
    // check nst is as determined in readyfRoutine()
    if (nstI != nstR) {
	cerr << "FATAL fElement::getTangentStiff() problems with incompatable nst";
	cerr << " ready: " << nstR << " invoke: " << nstI << endl;
	exit(-1);
    }

    // negate the sign of the loads -- feap elements return -ku
    (*fElementV[nstR]) *= -1.0;
    
    // add the applied loads from other sources
    (*fElementV[nstR]) -= *theLoad;
    
    // return the matrix
    return *(fElementV[nstR]);    
}


const Vector &
fElement::getResistingForceIncInertia()
{	
    // check for quick return
    if (nen == 0)
	return (*fElementV[0]);
    
    // get the current load factor
    Domain *theDomain=this->getDomain();
    double dm = theDomain->getCurrentTime();
    
    // set ctan, ior and iow
    double ctan[3];
    ctan[0] = 0.0; ctan[1] = 0.0; ctan[2] = 0.0;
    int ior = 0; int iow = 0;
    
    // call the ready routine to set ul, xl, tl and ix, NH1, NH2 and NH3
    int NH1, NH2, NH3;    
    int nstR = this->readyfRoutine(true);

    // zero the vector
    fElementV[nstR]->Zero();        
    
    // invoke the fortran subroutine
    int isw = 6; int nst = nen*ndf; int n = this->getTag();
    int nstI = this->invokefRoutine(ior, iow, ctan, isw);
    
    // check nst is as determined in readyfRoutine()
    if (nstI != nstR) {
	cerr << "FATAL fElement::getTangentStiff() problems with incompatable nst";
	cerr << " ready: " << nstR << " invoke: " << nstI << endl;
	exit(-1);
    }
    
    // return the matrix
    return *(fElementV[nstR]);    
}


int
fElement::sendSelf(int commitTag, Channel &theChannel)
{
    return 0;
}

int
fElement::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    return 0;
}


int
fElement::displaySelf(Renderer &theViewer, int displayMode, float fact)
{
    return 0;
}


void
fElement::Print(ostream &s, int flag)
{
    int ior = 0; int iow = 1;    
    if (s == cerr || s == cout)
	ior = -1;
    else  {
	s << "fElement::Print() - can only print to cerr or cout at present\n";
	ior = -1;
    }

    // get the current load factor
    Domain *theDomain=this->getDomain();
    double dm = theDomain->getCurrentTime();
    
    // set ctan, ior and iow
    double ctan[3];
    ctan[0] = 0.0; ctan[1] = 0.0; ctan[2] = 0.0;

    
    // call the ready routine to set ul, xl, tl and ix, NH1, NH2 and NH3
    int NH1, NH2, NH3;    
    int nstR = this->readyfRoutine(false);

    // invoke the fortran subroutine
    int isw = 4; int nst = nen*ndf; int n = this->getTag();
    int nstI = this->invokefRoutine(ior, iow, ctan, isw);
}

#ifdef _WIN32

extern "C" int _stdcall GETCOMMON(int *mynh1, int *mynh3, int *sizeH, 
				  double *myh);


extern "C" int _stdcall FILLCOMMON(int *mynen, double *mydm, int *myn, 
				   int *myior, int *myiow, int *mynh1, 
				   int *mynh2, int *mynh3, int *sumnh, 
				   double *myh, double *myctan,
				   int *nrCount);

extern "C" int _stdcall ELMT01(double *d, double *ul, double *xl, int *ix, 
			       double *tl, double *s, double *r, int *ndf, 
			       int *ndm, int *nst, int *isw);

extern "C" int _stdcall ELMT02(double *d, double *ul, double *xl, int *ix, 
			       double *tl, double *s, double *r, int *ndf, 
			       int *ndm, int *nst, int *isw);
		       
extern "C" int _stdcall ELMT03(double *d, double *ul, double *xl, int *ix, 
			       double *tl, double *s, double *r, int *ndf, 
			       int *ndm, int *nst, int *isw);
		       
extern "C" int _stdcall ELMT04(double *d, double *ul, double *xl, int *ix, 
			       double *tl, double *s, double *r, int *ndf, 
			       int *ndm, int *nst, int *isw);
		       
extern "C" int _stdcall ELMT05(double *d, double *ul, double *xl, int *ix, 
			       double *tl, double *s, double *r, int *ndf, 
			       int *ndm, int *nst, int *isw);		       

#define getcommon_ 	GETCOMMON
#define fillcommon_	FILLCOMMON
#define elmt01_		ELMT01
#define elmt02_		ELMT02
#define elmt03_		ELMT03
#define elmt04_		ELMT03
#define elmt05_		ELMT05

#else
extern "C" int getcommon_(int *mynh1, int *mynh3, int *sizeH, double *myh);


extern "C" int fillcommon_(int *mynen, double *mydm, int *myn, int *myior, 
                           int *myiow, int *mynh1, int *mynh2, int *mynh3,
                           int *sumnh, double *myh, double *myctan,
			   int *nrCount);

extern "C" int elmt01_(double *d, double *ul, double *xl, int *ix, double *tl, 
                       double *s, double *r, int *ndf, int *ndm, int *nst, 
		       int *isw);

extern "C" int elmt02_(double *d, double *ul, double *xl, int *ix, double *tl, 
                       double *s, double *r, int *ndf, int *ndm, int *nst, 
		       int *isw);
		       
extern "C" int elmt03_(double *d, double *ul, double *xl, int *ix, double *tl, 
                       double *s, double *r, int *ndf, int *ndm, int *nst, 
		       int *isw);
		       
extern "C" int elmt04_(double *d, double *ul, double *xl, int *ix, double *tl, 
                       double *s, double *r, int *ndf, int *ndm, int *nst, 
		       int *isw);
		       
extern "C" int elmt05_(double *d, double *ul, double *xl, int *ix, double *tl, 
                       double *s, double *r, int *ndf, int *ndm, int *nst, 
		       int *isw); 
#endif	       

int
fElement::invokefRoutine(int ior, int iow, double *ctan, int isw)
{
    // fill the common blocks
    // determine position in h of nh1, nh2 and nh3 - remember Fortarn indexing
    int NH1, NH2, NH3;
    if (nh1 != 0) { 
	NH1 = 1; 
	NH2 = nh1 + NH1; 
	NH3 = nh1 + NH2; 
    } else {
	NH1 = 1;
	NH2 = 1;
	NH3 = 1;
    }
    
    int NDM = ndm;
    int NDF = ndf;
    
    int n = this->getTag();
    int sum = 2*nh1 + nh3;    
    int count = nrCount;


    double dm = 0.0; // load factor

    fillcommon_(&nen, &dm, &n, &ior, &iow, &NH1, &NH2, &NH3, &sum, 
		h, ctan, &count);

    // invoke the fortran subroutine

    int nst = nen*ndf;
    if (nst != 0) {
	if (eleType == 1)
	    elmt01_(d,ul,xl,ix,tl,s,r,&NDF,&NDM,&nst,&isw);
	else if (eleType == 2)
	    elmt02_(d,ul,xl,ix,tl,s,r,&NDF,&NDM,&nst,&isw);	    
	else if (eleType == 3)
	    elmt03_(d,ul,xl,ix,tl,s,r,&NDF,&NDM,&nst,&isw);	    
	else if (eleType == 4)
	    elmt04_(d,ul,xl,ix,tl,s,r,&NDF,&NDM,&nst,&isw);	    
	else if (eleType == 5) 
	    elmt05_(d,ul,xl,ix,tl,s,r,&ndf,&NDM,&nst,&isw);	    
	else {
	    cerr << "fElement::invokefRoutine() unknown element type ";
	    cerr << eleType << endl;
	}
	
	// now copy the stuff from common block to h array
	getcommon_(&NH1,&NH3,&sum,h);
    }

    return nst;
}


int
fElement::invokefInit(int isw, int iow)
{
    // fill the common blocks
    // determine position in h of nh1, nh2 and nh3 - remember Fortarn indexing
    int NH1 =0;
    int NH2 =0;
    int NH3 =0;

    int NDM = ndm;
    int NDF = ndf;   
    double ctan[3];    

    int n = this->getTag();
    int sum = 0;    
    int ior = 0;
    int count = nrCount;

    double dm = 0.0;
    
    fillcommon_(&nen, &dm, &n, &ior, &iow, &NH1, &NH2, &NH3, &sum, 
		h, ctan, &count);

    // invoke the fortran subroutine

    int nst = nen*ndf;
    if (nst != 0) {
	if (eleType == 1)
	    elmt01_(d,ul,xl,ix,tl,s,r,&NDF,&NDM,&nst,&isw);
	else if (eleType == 2)
	    elmt02_(d,ul,xl,ix,tl,s,r,&NDF,&NDM,&nst,&isw);	    
	else if (eleType == 3)
	    elmt03_(d,ul,xl,ix,tl,s,r,&NDF,&NDM,&nst,&isw);	    
	else if (eleType == 4)
	    elmt04_(d,ul,xl,ix,tl,s,r,&NDF,&NDM,&nst,&isw);	    
	else if (eleType == 5)
	    elmt05_(d,ul,xl,ix,tl,s,r,&NDF,&NDM,&nst,&isw);	    
	else {
	    cerr << "fElement::invokefRoutine() unknown element type ";
	    cerr << eleType << endl;
	}

	if (nst < 0) {
	    cerr << "FATAL: fElement::fElement() - eleTag: " << this->getTag();
	    cerr << " ran out of memory creating h of size " << nst << endl;
	    exit(-1);
	}
    }

    // now get the size of the state info needed by the element
    sum = 0;
    getcommon_(&NH1,&NH3,&sum,h);
    nh1 = NH1; nh3=NH3;
	return 0;
}

int
fElement::readyfRoutine(bool incInertia)
{
    // determine nst 
    int nst = ndf*nen;

    // loop over nodes - fill in xl, ul, ix as we go
    int posUl = 0;
    int posXl = 0;
    for (int j=0; j<nen; j++) {
	Node *theNode = theNodes[j];

        // add the node tag to ix
	ix[j] = theNode->getTag();
        
        // add displacement, velocity, accel and  increments to ul
	// Note: we get nodal vel and accel only if inertia is true, this
	// will save memory in static analysis -- look at Node implementation	
	const Vector &trialDisp = theNode->getTrialDisp();
	const Vector &commitDisp = theNode->getDisp();        
        const Vector &crd = theNode->getCrds();

	// add the coordinates to xl		
	int crdSize = crd.Size();

	for (int i=0; i<crdSize; i++) {
	    xl[posXl] = crd(i);	    
	    posXl++;
	}
	
	if (incInertia == true) { 
	    const Vector &trialVel = theNode->getTrialVel();
	    const Vector &trialAccel = theNode->getTrialAccel();
	    const Vector &commitVel = theNode->getVel();        
	    for (int i=0; i<trialDisp.Size(); i++) {
		double trialDispI = trialDisp(i);
		ul[posUl] = trialDispI;
		ul[posUl+nst] = trialDispI - commitDisp(i);
		ul[posUl+2*nst] = trialDispI - u[posUl];	    	    
		ul[posUl+3*nst] = trialVel(i);
		ul[posUl+4*nst] = trialAccel(i);	    
		ul[posUl+5*nst] = commitVel(i);	    	    		
		u[posUl] = trialDispI; // u(k-1) on next call
		posUl++;				    
	    }
	} else {
	    for (int i=0; i<trialDisp.Size(); i++) {
		double trialDispI = trialDisp(i);
		ul[posUl] = trialDispI;
		ul[posUl+nst] = trialDispI - commitDisp(i);
		ul[posUl+2*nst] = trialDispI - u[posUl];	    	    
		ul[posUl+3*nst] = 0.0;
		ul[posUl+4*nst] = 0.0;	    
		ul[posUl+5*nst] = 0.0;	    	    
		u[posUl] = trialDispI;	    
		posUl++;
	    }
	}
    }

    // check we have a matrix and vector created for an object of this size
    if (fElementM[nst] == 0) {
	fElementM[nst] = new Matrix(s,nst,nst);
    	fElementV[nst] = new Vector(r,nst);
    
	if (fElementM[nst] == 0 || fElementV[nst] == 0) {
	    cerr << "FATAL fElement::getTangentStiff() nst: " << nst;
	    cerr << "ran out of memory\n";
	    exit(-1);
	}  
    }
    return nst;
}



int
fElement::update()
{
    // determine nst 
    int nst = ndf*nen;

    // increment the newton-raphson count -- needed for Prof. Fillipou's element
    nrCount++;
    
    return 0;
}




