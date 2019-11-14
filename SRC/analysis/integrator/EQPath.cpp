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
                                                                        
// Written: M. Salehi opensees.net@gmail.com
// Created: 02/19
// Revision: A

//refs
//Structural Engineeringand Mechanics   Volume 48, Number 6, December25 2013, pages 849 - 878
//DOI: https://doi.org/10.12989/sem.2013.48.6.849	
//
//Comprehensive evaluation of structural geometrical nonlinear solution techniques Part I : Formulation and characteristics of the methods
//M.Rezaiee - Pajand, M.Ghalishooyan and M.Salehi - Ahmadabad
//FULLTEXT : https://www.researchgate.net/publication/264146397_Comprehensive_evaluation_of_structural_geometrical_nonlinear_solution_techniques_Part_I_Formulation_and_characteristics_of_the_methods


//Structural Engineeringand Mechanics   Volume 48, Number 6, December25 2013, pages 879 - 914
//DOI: https://doi.org/10.12989/sem.2013.48.6.879	
//
//Comprehensive evaluation of structural geometrical nonlinear solution techniques Part II : Comparing efficiencies of the methods
//M.Rezaiee - Pajand, M.Ghalishooyan and M.Salehi - Ahmadabad
//FULLTEXT : https://www.researchgate.net/publication/263361974_Comprehensive_evaluation_of_structural_geometrical_nonlinear_solution_techniques_Part_II_Comparing_efficiencies_of_the_methods

#include <EQPath.h>
#include <AnalysisModel.h>
#include <LinearSOE.h>
#include <Vector.h>
#include <Channel.h>
#include <math.h>
#include "ElementIter.h"
#include "Domain.h"
#include "Node.h"
#include "Element.h"

EQPath::EQPath(double arcLen,int method)
:StaticIntegrator(INTEGRATOR_TAGS_EQPath),
 arclen(arcLen), 
 uq(0),uq0(0), ur(0), du(0),du0(0),uqn(0), q(0),type(method),m(1),dl(0),changed(0),nitr(0)
 {
	
	
  
}

EQPath::~EQPath()
{
    // delete any vector object created
    if (uq != 0)
	delete uq;
    if (uq0 != 0)
	   delete uq0;
    if (uqn != 0)
	   delete uqn;
    if (ur != 0)
	delete ur;
    if (du != 0)
	delete du;
    if (du0 != 0)
	   delete du0;
    if (q != 0)
	delete q;
}

int
EQPath::newStep(void)
{
    // get pointers to AnalysisModel and LinearSOE
    AnalysisModel *theModel = this->getAnalysisModel();
    LinearSOE *theLinSOE = this->getLinearSOE();    
    if (theModel == 0 || theLinSOE == 0) {
	opserr << "WARNING EQPath::newStep() ";
	opserr << "No AnalysisModel or LinearSOE has been set\n";
	return -1;
    }

    // get the current load factor
    double currentLambda = theModel->getCurrentDomainTime();


    // determine dUhat
    this->formTangent();
    theLinSOE->setB(*q);
    if (theLinSOE->solve() < 0) {
      opserr << "EQPath::newStep(void) - failed in solver\n";
      return -1;
    }

    // for GDC method we need this
    if(uqn==0 && uq0!=0)
    {
	   uqn=new Vector(uq0->Size());
	   (*uqn)=(*uq0);
    }
    else if(uq0!=0)
    {
	   (*uqn)=(*uq0);
    }

    uq0=new Vector(du->Size());
    (*uq0) = theLinSOE->getX();
    int size = theModel->getNumEqn();
    // determine delta lambda(1) == dlambda
    
   double a=(*du)^(*uq0);
 
   if(a>=0)
	  sign=1;
   else
	  sign=-1;
    
    du->Zero();

    double dLambda = sign*arclen/uq0->Norm();

    (*du)=dLambda*(*uq0);

    du0=new Vector(du->Size());
    (*du0)=(*du);

    currentLambda += dLambda;
    dl+=dLambda;
    
    // update model with delta lambda and delta U
    theModel->incrDisp(*du);    
    theModel->applyLoadDomain(currentLambda);    
    if (theModel->updateDomain() < 0) {
      opserr << "EQPath::newStep - model failed to update for new dU\n";
      return -1;
    }
    
    //opserr << "new step : current dl " << dLambda << "\n" ;

    nitr=0;
    if (m != 1)
	   changed--;
    if (changed == 0)
	   m = 1;

    return 0;
}

int
EQPath::update(const Vector &dU)
{
    AnalysisModel *theModel = this->getAnalysisModel();
    LinearSOE *theLinSOE = this->getLinearSOE();    
    if (theModel == 0 || theLinSOE == 0) {
	opserr << "WARNING EQPath::update() ";
	opserr << "No AnalysisModel or LinearSOE has been set\n";
	return -1;
    }

    nitr++;

    (*ur) = dU;

    theLinSOE->setB(*q);
    theLinSOE->solve();
    (*uq) = theLinSOE->getX();    

    // determine delta lambda(i)
    double dLambda=0;
	   if (type==1) // minimum residual disp
	   {
		  double a = (*ur)^(*uq);
		  double b = (*uq)^(*uq);
		  if (b == 0) {
			 opserr << "EQPath::update() - zero denominator\n";
			 return -1;
		  }	  
		  dLambda = -a/b;
	   }
	   else if (type==2) // normal plain
	   {
		  double a = (*du0)^(*ur);
		  double b = (*du0)^(*uq);
		  if (b == 0) {
			 opserr << "EQPath::update() - zero denominator\n";
			 return -1;
		  }

		  dLambda = -a/b;
	   }
	   else if (type==3) // update normal plain
	   {
		  double a = (*du)^(*ur);
		  double b = (*du)^(*uq);
		  if (b == 0) {
			 opserr << "EQPath::update() - zero denominator\n";
			 return -1;
		  }

		  dLambda = -a/b;
		 
	   }
	   else if (type==4) // cylindrical arc-length
	   {
		  double A = (*uq)^(*uq);
		  double B = 2 * (((*du) + (*ur)) ^ (*uq));
		  double c1=((*du) + (*ur)) ^ (*du);
		  double c2=((*du) + (*ur)) ^(*ur);
		  double C = c1 + c2 - arclen * arclen;
		  double delta = B * B - 4 * A * C;

		  dLambda = 0;
		  if (delta < 0)
		  {
			 opserr << "EQPath::update() - negetive denominator\n";
			 return -1;
		  }
		  else if (delta == 0)
		  {
			 dLambda = -B / 2 / A;
		  }
		  else
		  {
			 double sl1 = (-B + pow(delta,0.5)) / 2 / A;
			 double sl2 = (-B - pow(delta,0.5)) / 2 / A;
			 double aa1=(*du) ^ (*ur);
			 double aa2=(*du) ^ (*du);
			 double aa3=(*du) ^ (*uq);
			 double costl1 = aa1 + aa2 + sl1 * aa3;
			 double costl2 = aa1 + aa2 + sl2 * aa3;
			 dLambda = sl1;
			 if (costl2 > costl1)
				dLambda = sl2;
		  }

	   }
	   else if (type==5) // modified normal flow
	   {
		  double a = (*ur)^(*uq);
		  double b = (*uq)^(*uq);
		  if (b == 0) {
			 opserr << "EQPath::update() - zero denominator\n";
			 return -1;
		  }

		  dLambda = -a/b;
	   }

	   else if (type==6) // GDC
	   {
		  double a,b;
		  if(uqn==0)
		  {
			 // first increment use MRD
			 a = (*ur)^(*uq);
			 b = (*uq)^(*uq);
		  }
		  else
		  {
			 a = (*ur)^(*uqn);
			 b = (*uq)^(*uqn);
		  }
		  
	   
		  if (b == 0) {
		  opserr << "EQPath::update() - zero denominator\n";
		  return -1;
		  }
		  dLambda = -a/b;
		  
	   }
	   else if (type==7) // Modified Update Normal Plane
	   {
		  double p1=(*uq)^(*uq);
		  double p2=(*du)^(*uq);
		  double p3=(*ur)^(*uq);
		  double p4=(*ur)^(*du);
		  double p5=(*ur)^(*ur);
		  double dl0=-p3/p1; //Chan constraint


		  double A = p1;//(*uq)^(*uq);
		  double B = p2+2*p3;//(((*du) + (*ur)+ (*ur)) ^ (*uq));
		  double C=p4+p5;//((*du) + (*ur)) ^ (*ur);
		  double delta = B * B - 4 * A * C;
		  //opserr << "negative denominator - alpha = " << 0 <<"\n";
		  dLambda = 0;
		  if (delta < 0)
		  {
			 Vector *v1=new Vector(ur->Size());
			 Vector *v2=new Vector(ur->Size());
			 (*v2)=(*ur);
			 v2->addVector(1,*uq,dl0);
			 (*v1)=(*du);
			 v1->addVector(1,*v2,1);
			 double l1=v1->Norm();
			 double l2=v2->Norm();
			 double alpha=(C-B*B/4/A)/l1/l2;
			 alpha+=0.1*(1-alpha);
			 //opserr << "negative denominator - alpha = " << alpha <<"\n";
			 delta = B * B - 4 * A * (C-alpha*l1*l2);
			 //dLambda=dl0;
		  }

		  if (delta == 0)
		  {
			 dLambda = -B / 2 / A;
		  }
		  else
		  {
			 double sl1 = (-B + pow(delta,0.5)) / 2 / A;
			 double sl2 = (-B - pow(delta,0.5)) / 2 / A;
			 double aa1=(*du) ^ (*ur);
			 double aa2=(*du) ^ (*du);
			 double aa3=(*du) ^ (*uq);
			 double costl1 = aa1 + aa2 + sl1 * aa3;
			 double costl2 = aa1 + aa2 + sl2 * aa3;
			 dLambda = sl1;
			 if (costl2 > costl1)
				dLambda = sl2;
		  }

	   }
	   else if (type==8) // normal plain
	   {
		  double p1=(*uq)^(*uq);
		  double p2=(*du)^(*uq);
		  double p3=(*ur)^(*uq);
		  if (p1 == 0) {
			 opserr << "EQPath::update() - zero denominator\n";
			 return -1;
		  }

		  dLambda = -(p2+p3)/p1;
	   }
	   else if (type==9) // local minimum residual displacement
	   {
		   AnalysisModel *am=this->getAnalysisModel();
		   Domain *theDomain=theModel->getDomainPtr();
		   Element *elePtr;
		   ElementIter &theElemIter = theDomain->getElements();    
		   while ((elePtr = theElemIter()) != 0)
		   {
			   Node	**nodes=elePtr->getNodePtrs();
			   
		   }
	   }
	   else {
		  opserr << "WARNING EQPath::update() ";
		  opserr << "Unknown update method has been set\n";
		  return -1;
	   }



    Vector *sd=new Vector(ur->Size());
    // determine delta U(i)
    (*sd) = (*ur);
    sd->addVector(1.0, *uq,dLambda);
    

    if (type==5)
    {
	   double fac = -((*sd)^(*uq))/((*uq0)^(*uq0));	  
	   sd->addVector(1.0,*uq,fac);
    }

    (*du) += (*sd);
    dl += dLambda;

    double currentLambda = theModel->getCurrentDomainTime();
    currentLambda+=dLambda;
    // update the model
    theModel->incrDisp(*sd);    
    theModel->applyLoadDomain(currentLambda);    

    if (theModel->updateDomain() < 0) {
      opserr << "EQPath::update - model failed to update for new dU\n";
      return -1;
    }
    
    // set the X soln in linearSOE to be deltaU for convergence Test
    theLinSOE->setX(*sd);

    
    return 0;
}

int 
EQPath::domainChanged(void)
{
    // we first create the Vectors needed
    AnalysisModel *theModel = this->getAnalysisModel();
    LinearSOE *theLinSOE = this->getLinearSOE();    
    if (theModel == 0 || theLinSOE == 0) {
	opserr << "WARNING EQPath::update() ";
	opserr << "No AnalysisModel or LinearSOE has been set\n";
	return -1;
    }    

    int size = theModel->getNumEqn(); // ask model in case N+1 space

    if (uq == 0 || uq->Size() != size) { // create new Vector
	if (uq != 0)
	    delete uq;   // delete the old
	uq = new Vector(size);
	if (uq == 0 || uq->Size() != size) { // check got it
	    opserr << "FATAL EQPath::domainChanged() - ran out of memory for";
	    opserr << " uq Vector of size " << size << endln;
	    exit(-1);
	}
    }

    if (du == 0 || du->Size() != size) { // create new Vector
	if (du != 0)
	    delete du;   // delete the old
	du = new Vector(size);
	if (du == 0 || du->Size() != size) { // check got it
	    opserr << "FATAL EQPath::domainChanged() - ran out of memory for";
	    opserr << " du Vector of size " << size << endln;
	    exit(-1);
	}
    }

    
    if (ur == 0 || ur->Size() != size) { // create new Vector
	if (ur != 0)
	    delete ur;   // delete the old
	ur = new Vector(size);
	if (ur == 0 || ur->Size() != size) { // check got it
	    opserr << "FATAL EQPath::domainChanged() - ran out of memory for";
	    opserr << " deltaU Vector of size " << size << endln;
	    exit(-1);
	}
    }

    if (q == 0 || q->Size() != size) { 
	if (q != 0)
	    delete q;  
	q = new Vector(size);
	if (q == 0 || q->Size() != size) { 
	    opserr << "FATAL EQPath::domainChanged() - ran out of memory for";
	    opserr << " q Vector of size " << size << endln;
	    exit(-1);
	}
    }

    // now we have to determine phat
    // do this by incrementing lambda by 1, applying load
    // and getting phat from unbalance.
    double currentLambda = theModel->getCurrentDomainTime();
    currentLambda += 1.0;
    theModel->applyLoadDomain(currentLambda);    
    this->formUnbalance(); // NOTE: this assumes unbalance at last was 0
    (*q) = theLinSOE->getB();
    currentLambda -= 1.0;
    theModel->setCurrentDomainTime(currentLambda);    


    // check there is a reference load
    int haveLoad = 0;
    for (int i=0; i<size; i++)
      if ( (*q)(i) != 0.0 ) {
	haveLoad = 1;
	i = size;
      }

    if (haveLoad == 0) {
      opserr << "WARNING ArcLength::domainChanged() - zero reference load";
      return -1;
    }

    
    return 0;
}

int
EQPath::sendSelf(int cTag,
		    Channel &theChannel)
{
  Vector data(3);
  data(0) = arclen;
  data(1) = dl;
  data(2) = sign;
  
  if (theChannel.sendVector(this->getDbTag(), cTag, data) < 0) {
      opserr << "EQPath::sendSelf() - failed to send the data\n";
      return -1;
  }
  return 0;
}


int
EQPath::recvSelf(int cTag,
		    Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  Vector data(3);
  if (theChannel.recvVector(this->getDbTag(), cTag, data) < 0) {
      opserr << "EQPath::sendSelf() - failed to send the data\n";
      return -1;
  }      

  // set the data

  arclen = data(0);
  dl = data(1);
  sign = data(2);
  
  return 0;
}

void
EQPath::Print(OPS_Stream &s, int flag)
{
    AnalysisModel *theModel = this->getAnalysisModel();
    if (theModel != 0) {
	double cLambda = theModel->getCurrentDomainTime();
	s << "\t EQPath - currentLambda: " << cLambda <<"\n";
	s << "\t EQPath - arcLength: " << arclen <<"\n";
	s << "\t EQPath - sign: " << sign <<"\n";
    } else 
	s << "\t EQPath - no associated AnalysisModel\n";
}

