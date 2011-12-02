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
// $Date: 2000-09-15 08:23:22 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/ParallelMaterial.cpp,v $
                                                                        
                                                                        
// File: ~/material/ParallelModel.C
//
// Written: fmk 
// Created: 07/98
// Revision: A
//
// Description: This file contains the class definition for 
// ParallelModel. ParallelModel is an aggregation
// of UniaxialMaterial objects all considered acting in parallel.
//
// What: "@(#) ParallelModel.C, revA"

#include <ParallelMaterial.h>
#include <ID.h>
#include <Vector.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <stdlib.h>


ParallelMaterial::ParallelMaterial(
				 int tag, 
				 int num, 
				 UniaxialMaterial ** theMaterialModels,
                                 double min, double max)
:UniaxialMaterial(tag,MAT_TAG_ParallelMaterial),
 trialStrain(0.0), numMaterials(num), otherDbTag(0), theModels(0), 
 epsmin(min), epsmax(max), Cfailed(0), Tfailed(0)
{
    // create an array (theModels) to store copies of the MaterialModels
    theModels = new UniaxialMaterial *[num];

    if (theModels == 0) {
	cerr << "FATAL ParallelMaterial::ParallelMaterial() ";
	cerr << " ran out of memory for array of size: " << num << "\n";
	exit(-1);
    }

    // into the newly created array store a ponter to a copy
    // of the UniaxialMaterial stored in theMaterialModels
    for (int i=0; i<num; i++) {
	theModels[i] = theMaterialModels[i]->getCopy();
    }
}



// this constructor is used for a ParallelMaterailModel object that
// needs to be constructed in a remote actor process. recvSelf() needs
// to be called on this object
ParallelMaterial::ParallelMaterial()
:UniaxialMaterial(0,MAT_TAG_ParallelMaterial),
 trialStrain(0.0), numMaterials(0), otherDbTag(0), theModels(0), 
 epsmin(NEG_INF_STRAIN), epsmax(POS_INF_STRAIN), Cfailed(0), Tfailed(0)
{

}


ParallelMaterial::~ParallelMaterial()
{
    // invoke the destructor on each MaterialModel object
    for (int i=0; i<numMaterials; i++)
	delete theModels[i];

    // now we can delete the array
    if (theModels != 0) // just in case blank constructor called and no recvSelf()
	delete [] theModels;
}



int 
ParallelMaterial::setTrialStrain(double strain, double strainRate)
{
    // set the trialStrain and the trialStrain in each of the
    // local MaterialModel objects 
	Tfailed = Cfailed;
    trialStrain = strain;
	trialStrainRate = strainRate;

    if (trialStrain < epsmin || trialStrain > epsmax)
        Tfailed = 1;

    if (!Tfailed) {
        for (int i=0; i<numMaterials; i++)
	    theModels[i]->setTrialStrain(strain, strainRate);
    }

    return 0;
}


double 
ParallelMaterial::getStrain(void)
{
    return trialStrain;
}

double 
ParallelMaterial::getStrainRate(void)
{
    return trialStrainRate;
}

double 
ParallelMaterial::getStress(void)
{
    // get the stress = sum of stress in all local MaterialModel objects
    double stress = 0.0;
    if (!Tfailed) {
        for (int i=0; i<numMaterials; i++)
	    stress +=theModels[i]->getStress();
    }

    return stress;
}



double 
ParallelMaterial::getTangent(void)
{
    // get the tangent = sum of tangents in all local MaterialModel objects    
    double E = 0.0;
    if (!Tfailed) {
        for (int i=0; i<numMaterials; i++)
	    E +=theModels[i]->getTangent();    
    }

    return E;
}

double 
ParallelMaterial::getDampTangent(void)
{
    // get the damp tangent = sum of damp tangents in all local MaterialModel objects    
    double eta = 0.0;
    if (!Tfailed) {
        for (int i=0; i<numMaterials; i++)
	    eta +=theModels[i]->getDampTangent();    
    }

    return eta;
}

double ParallelMaterial::getSecant ()
{
	// get the secant = sum of secants in all local MaterialModel objects    
    double S = 0.0;
    if (!Tfailed) {
        for (int i=0; i<numMaterials; i++)
	    S +=theModels[i]->getSecant();    
    }

    return S;
}

int 
ParallelMaterial::commitState(void)
{
    Cfailed = Tfailed;

    // invoke commitState() on each of local MaterialModel objects
    for (int i=0; i<numMaterials; i++)
	if (theModels[i]->commitState() != 0) {
	    cerr << "WARNING ParallelMaterial::commitState() ";
	    cerr << "MaterialModel failed to commitState():" ;
	    theModels[i]->Print(cerr);
	}
    
    return 0;    
}

int 
ParallelMaterial::revertToLastCommit(void)
{
    Tfailed = Cfailed;

    // invoke commitState() on each of local MaterialModel objects
    for (int i=0; i<numMaterials; i++)
	if (theModels[i]->revertToLastCommit() != 0) {
	    cerr << "WARNING ParallelMaterial::revertToLastCommit() ";
	    cerr << "MaterialModel failed to revertToLastCommit():" ;
	    theModels[i]->Print(cerr);
	}
    
    return 0;    
}


int 
ParallelMaterial::revertToStart(void)
{
    Cfailed = 0;

    // invoke commitState() on each of local MaterialModel objects
    for (int i=0; i<numMaterials; i++)
	if (theModels[i]->revertToStart() != 0) {
	    cerr << "WARNING ParallelMaterial::revertToStart() ";
	    cerr << "MaterialModel failed to revertToStart():" ;
	    theModels[i]->Print(cerr);
	}
    
    return 0;    
}



UniaxialMaterial *
ParallelMaterial::getCopy(void)
{
    ParallelMaterial *theCopy = new 
       ParallelMaterial(this->getTag(),numMaterials,theModels,epsmin,epsmax);

    theCopy->trialStrain = trialStrain;
	theCopy->trialStrainRate = trialStrainRate;
    theCopy->Cfailed = Cfailed;
    theCopy->Tfailed = Tfailed;

    return theCopy;
}


int 
ParallelMaterial::sendSelf(int cTag, Channel &theChannel)
{
    Vector data(6);

    // check to see if we have a database tag for the other ID,
    // the object requires two database tags as sends two ID's.
    if (otherDbTag == 0) 
      otherDbTag = theChannel.getDbTag();

    data(0) = this->getTag();
    data(1) = numMaterials;
    data(2) = otherDbTag;
    data(3) = epsmin;
    data(4) = epsmax;
    data(5) = Cfailed;

    theChannel.sendVector(this->getDbTag(), cTag, data);
    
    // now create an ID containing the class tags and dbTags of all
    // the MaterialModel objects in this ParallelMaterial
    // then send each of the MaterialModel objects
    ID classTags(numMaterials*2);
    for (int i=0; i<numMaterials; i++) {
	classTags(i) = theModels[i]->getClassTag();
	int matDbTag = theModels[i]->getDbTag();
	if (matDbTag == 0) {
	  matDbTag  = theChannel.getDbTag();
	  if (matDbTag != 0)
	    theModels[i]->setDbTag(matDbTag);
	}
	classTags(i+numMaterials) = matDbTag;
    }
    theChannel.sendID(otherDbTag, cTag, classTags);
    for (int j=0; j<numMaterials; j++)
	theModels[j]->sendSelf(cTag, theChannel);
    
    return 0;
}

int 
ParallelMaterial::recvSelf(int cTag, Channel &theChannel, 
				FEM_ObjectBroker &theBroker)
{
    Vector data(6);
    theChannel.recvVector(this->getDbTag(), cTag, data);
    this->setTag(int(data(0)));
    numMaterials = int(data(1));
    otherDbTag = int(data(2));
    epsmin = data(3);
    epsmax = data(4);
    Cfailed = int(data(5));

    // create and receive an ID for the classTags and dbTags of the local 
    // MaterialModel objects
    ID classTags(numMaterials*2);
    theChannel.recvID(otherDbTag, cTag, classTags);
    theModels = new UniaxialMaterial *[numMaterials];

    // create space for the objects
    if (theModels == 0) {
	cerr << "FATAL ParallelMaterial::recvSelf() ";
	cerr << " ran out of memory for array of size: " << numMaterials << "\n";
	exit(-1);
    }    

    // now for each of the MaterialModel objects, create a new object
    // and invoke recvSelf() on it
    for (int i=0; i<numMaterials; i++) {
	UniaxialMaterial *theMaterialModel = 
	    theBroker.getNewUniaxialMaterial(classTags(i));
	if (theMaterialModel != 0) {
	    theModels[i] = theMaterialModel;
	    theMaterialModel->setDbTag(classTags(i+numMaterials));
	}
	else {
	    cerr << "FATAL ParallelMaterial::recvSelf() ";
	    cerr << " could not get a UniaxialMaterial \n";
	    exit(-1);
	}    	    
	theMaterialModel->recvSelf(cTag, theChannel, theBroker);
    }
    return 0;
}

void 
ParallelMaterial::Print(ostream &s, int flag)
{
    s << "Parallel tag: " << this->getTag() << endl;
    if (epsmin != NEG_INF_STRAIN)
      s << "  epsmin: " << epsmin << endl;
    if (epsmax != POS_INF_STRAIN)
      s << "  epsmax: " << epsmax << endl;
    for (int i=0; i<numMaterials; i++) {
      s << ' ';
      theModels[i]->Print(s, flag);
    }
    
}
