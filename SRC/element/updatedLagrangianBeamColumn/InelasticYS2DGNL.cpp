// InelasticYS2DGNL.cpp
//////////////////////////////////////////////////////////////////////

#include <Domain.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
//#include <Renderer.h>
#include <math.h>
#include <stdlib.h>

#ifdef _NOGRAPHICS

#else
#ifdef _GLX    // Boris Jeremic added 23Oct2002
#include <OpenGLRenderer.h>
#endif         // Boris Jeremic added 23Oct2002
#endif

#include "InelasticYS2DGNL.h"
#include <Renderer.h>
//#include <WindowManager.h>
#include <PlainMap.h>
#include <Information.h>
#include <ElementResponse.h>
#include <string.h>
//#define debug  1
//#define fdebug 1
//#define pdebug 1
#define updateDebug 0
#define plastkDebug 0

#define ERROR 1e-8
#define TAG_InelasticYS2DGNL -1

// Vector InelasticYS2DGNL::trialForce(6);
Vector InelasticYS2DGNL::elasticForce(6);
Vector InelasticYS2DGNL::F1(6);
Vector InelasticYS2DGNL::F2(6);
Vector InelasticYS2DGNL::Fs(6);
double InelasticYS2DGNL::storage(0);

const int InelasticYS2DGNL::INSIDE(-1);
const int InelasticYS2DGNL::WITHIN(0);
const int InelasticYS2DGNL::OUTSIDE(1);

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

// Element specific to PM interaction only
InelasticYS2DGNL::InelasticYS2DGNL(int tag, // double a, double e, double i,
                                   int Nd1, int Nd2,
				   YieldSurface_BC *ysEnd1,
				   YieldSurface_BC *ysEnd2,
				   int rf_algo, bool islinear, double rho)
  :UpdatedLagrangianBeam2D(tag, TAG_InelasticYS2DGNL, Nd1, Nd2, islinear),
   end1G(6,1), end2G(6,1), Stiff(6,6),
   forceRecoveryAlgo(rf_algo), forceRecoveryAlgo_orig(rf_algo),
   end1Damage(false), end2Damage(false), split_step(false)
{
  
  {
    debug  = 0;
    fdebug = 0;
    pdebug = 0;
    ydebug = 0;
    statusDebug = 0;
  }
  
  
  
  if(ysEnd1==0) {
      opserr << "WARNING - InelasticYS2DGNL(): ys1 = 0" << endln;
  } else {
    ys1 =  ysEnd1->getCopy();
    ys1->setTransformation(2, 0, -1,  1);
    ys1->setEleInfo(getTag(), 1);
    // Hogging moment at end1 => Negative moment
    // Positive disp at node 1 (compression) => Positive axial force
  }
  
  if(ysEnd2==0) {
      opserr << "WARNING - InelasticYS2DGNL(): ys2 = 0" << endln;
  } else {
    
    ys2 =  ysEnd2->getCopy();
    ys2->setTransformation(5, 3, 1, -1);
    ys2->setEleInfo(getTag(), 2);
		// Sagging moment at end2 positive => Positive moment
    // Positive disp at node 2 (tension) => Negative axial force
  }
  
  pView = 0;
  end1Plastify = false;
  end2Plastify = false;
  end1Plastify_hist = false;
  end2Plastify_hist = false;
  
  init = false;
}


InelasticYS2DGNL::~InelasticYS2DGNL()
{
  //~ if(theMap != NULL)
  	//~ delete theMap;
}

// do everything here
int InelasticYS2DGNL::update()
{
    if (L == 0)
	    return 0;

	ys1->update();
	ys2->update();		
			
//	split_step = false;

//////////////////////////////////////////////////////////////////////
// Step 1: Trial Elastic Forces
//////////////////////////////////////////////////////////////////////

	 // Get the local elastic stiffness matrix, store in Stiff
     getLocalStiff(Stiff);

     // Add internal geometric stiffness matrix
     addInternalGeomStiff(Stiff);

     // Get incremental local displacements  trial-conv.
     getIncrNaturalDisp(disp);
	 //getIncrLocalDisp(disp);   IMPORTANT - Do not change!!

     // Compute local incremental force
     force = Stiff*disp;

Vector trial_force(6);

     // Compute total trial local force - store in trialForce
     trial_force = eleForce_hist + force;	 
	 computeTrueEleForce(trial_force);
	 checkSpecialCases();
	 return 0;
}


int InelasticYS2DGNL::computeTrueEleForce(Vector &trialForce)
{
	 // ys1->displayForcePoint(trialForce);
	 // opserr << "\a";

	 // opserr << "Trial Force = " << trialForce; opserr << "\a";

//////////////////////////////////////////////////////////////////////
// Step 2: Check and plastify ends if required
//////////////////////////////////////////////////////////////////////
	// if(forceRecoveryAlgo == -1)
	{
		int plastify = this->plasticPredictor(trialForce);
		if(!plastify) // did not plastify
			return 0;
	}
	/*else
	{
		int plastify = this->elasticCorrector(trialForce); // needs debugging
		if(!plastify) // did not plastify
			return 0;
	}*/

//////////////////////////////////////////////////////////////////////
// Step 3: Drift control
//////////////////////////////////////////////////////////////////////

    if(end1Plastify)
    {

    	int d1 = ys1->getTrialForceLocation(eleForce);
    	if(d1 == OUTSIDE)
		{
    		ys1->setToSurface(eleForce, ys1->RadialReturn);
		//cout <<" ys1 outside\n";
		}
    	else
		{
    	    ys1->setToSurface(eleForce, ys1->ConstantYReturn);
		//cout << "ys1 inside \n";
		}
    }

    if(end2Plastify)
    {
    	int d2 = ys2->getTrialForceLocation(eleForce);
    	if(d2 == OUTSIDE)
		{
    		ys2->setToSurface(eleForce, ys2->RadialReturn);
		//cout << "ys2 outside\n";
		}
    	else
		{
    		ys2->setToSurface(eleForce, ys2->ConstantYReturn);
		//cout << "ys2 inside\n";
		}
    }

//////////////////////////////////////////////////////////////////////
// Step 4: Force Balance  first Axial, then Shear - still have to do dF for elastic case
//////////////////////////////////////////////////////////////////////
    forceBalance(eleForce, 1);

    /*if(
    	(eleForce(0) > 0 && eleForce_hist(0) < 0) ||
    	(eleForce(0) < 0 && eleForce_hist(0) > 0)
      )
    {
    	opserr << "NOTE: Axial force sign changed from converged state\n";
     // nothing can be done about that
//    	eleForce(0) =  eleForce_hist(0);
//    	eleForce(3) =  eleForce_hist(3);
    }*/
	return 0;
}

void InelasticYS2DGNL::checkSpecialCases(void)
{
	if(fabs(eleForce(0)) < ERROR && fabs(eleForce(3)) < ERROR)
	{
		eleForce(0) = 0.0;
		eleForce(3) = 0.0;
		return;
	}
    // check if the 2 force-points are on the same side in YS    
    if(sign(eleForce(0)) == sign(eleForce(3))) // special case
    {
		opserr << "oops 1: element " << getTag() << " okay \n";
		opserr << eleForce;
		// ys is already evolved
		
		// Stiff = Kt now
		// getLocalStiff(Stiff);
		// addInternalGeomStiff(Stiff);
		getIncrNaturalDisp(disp);
		force = Stiff*disp;
		eleForce = eleForce_hist + force;
		// elasticCorrector(trialForce, ys1->ConstantYReturn); - untested
		// trial force vector is stored in eleForce
		bool end1Drifts, end2Drifts;
		this->checkEndStatus(end1Drifts, end2Drifts, eleForce);
		// Have to use constantY only
        if(end1Plastify)
			ys1->setToSurface(eleForce, ys1->ConstantYReturn);
        if(end2Plastify)
			ys2->setToSurface(eleForce, ys2->ConstantYReturn);
		
        // the two P's should be same now, so Pavg should be = P's
        forceBalance(eleForce, 1); // shear
        if(sign(eleForce(0)) == sign(eleForce(3)))
        {
			opserr << "oops 2: element " << getTag() << " not okay \n";
			opserr << eleForce;
		}
	}

	if(updateDebug)
	{
    	if(pView)
		{
    		pView->clearImage();
    		pView->startImage();
	    	ys1->displaySelf(*pView, 1, 1);
	    	ys1->displayForcePoint(eleForce, 2);
	    	
	    	ys2->displaySelf(*pView, 1, 1);

	    	ys2->displayForcePoint(eleForce, 2);
	    	
	    	pView->doneImage();
	    	opserr << "Trial Force points plotted \n";
	    	opserr << "\a";
		}
	}
			
}
// update

// no longer valid
int InelasticYS2DGNL::elasticCorrector(Vector &trialForce, int algo)
{
bool	end1Drifts, end2Drifts;

	this->checkEndStatus(end1Drifts, end2Drifts, trialForce);

	if(!end1Plastify && !end2Plastify) // still elastic
	{
		// opserr << "InelasticYS2DGNL::elasticCorrector - Trial Force Elastic" << endln;
		eleForce = trialForce;
		return 0;
	}
	
	// if any end shoot through, force point is retracted
	// else if they drift, lamda is computed and the
	// surface modified, force point is retracted back to
	// the surface using the forceRecoveryAlgo

	if(end1Plastify)
	{
		this->plastifyOneEnd(1, ys1, trialForce, disp, Stiff, eleForce, algo);
	}

    if(end2Plastify)
	{
		this->plastifyOneEnd(2, ys2, trialForce, disp, Stiff, eleForce, algo);
	}

	return 1;
}


int InelasticYS2DGNL::plasticPredictor(Vector &trialForce)
{

bool	end1Drifts, end2Drifts;
Vector totalForce(6);

	this->checkEndStatus(end1Drifts, end2Drifts, trialForce);
		
	if(end1Plastify && !end2Plastify)
	{
		// this updates both Kt and F
		this->plastifyOneEnd(1, ys1, trialForce, disp, Stiff, eleForce, -1);
	}
	else if(end2Plastify && !end1Plastify)
	{
		// this updates both Kt and F
		this->plastifyOneEnd(2, ys2, trialForce, disp, Stiff, eleForce, -1);
	}
	else if(end1Plastify && end2Plastify)
	{
		// opserr << "drift status both plastify = " << end1Drifts << " " << end2Drifts << endln; opserr << "\a";
	
   	   if( (end1Drifts && !end2Drifts))
     	{     		
     		this->splitStep(2, ys2, ys1, trialForce, Stiff, eleForce);
    	}
   	    else if ((end2Drifts && !end1Drifts))
   	    {   	    	
         	 this->splitStep(1, ys1, ys2, trialForce, Stiff, eleForce);
	    }
    	 else
			this->plastifyBothEnds(trialForce, disp, Stiff, eleForce);

//   	   if( (end1Drifts && !end2Drifts) || (end2Drifts && !end1Drifts))
//     	{     		
//         	this->elasticCorrector(trialForce, -10);
//    	}
//	    else
//			this->plastifyBothEnds(trialForce, disp, Stiff, eleForce);
			
					
	 	/* if(!end1Drifts && !end2Drifts)
	 	{
	 	 	// order does not matter,
	 	 	// only the force point is retracted
	 	 	this->plastifyOneEnd(1, ys1, trialForce, disp, Stiff, eleForce, -1);
	 	 	this->plastifyOneEnd(2, ys2, trialForce, disp, Stiff, eleForce, -1);
	 	}
	 	else if(end1Drifts && !end2Drifts)
	 	{
	 	 	// end that drifts modifies Kt and F
	 	 	// then call the end that shoots through -
	 	 	// retracts the force point - order matters
	 	 	this->plastifyOneEnd(1, ys1, trialForce, disp, Stiff, eleForce, -1);
	 	 	this->plastifyOneEnd(2, ys2, trialForce, disp, Stiff, eleForce, -1);
	 	}
	 	else if(end2Drifts && !end1Drifts)

	 	{
	 	 	// end that drifts modifies Kt and F
	 	 	// then call the end that shoots through -
	 	 	// retracts the force point - order matters
	 	 	this->plastifyOneEnd(2, ys2, trialForce, disp, Stiff, eleForce, -1);
	 	 	this->plastifyOneEnd(1, ys1, trialForce, disp, Stiff, eleForce, -1);
	 	}
	 	else if(end1Drifts && end2Drifts)
	 	{
	 	 	this->plastifyBothEnds(trialForce, disp, Stiff, eleForce);
	 	}
	 	else
	 	{
	 	 	opserr << "InelasticYS2DGNL::predictor() - This condition should not happen\n";
	 	 	opserr << "\a";
	 	} */
	 	
	}
	else
	{
		if(!end1Plastify && !end2Plastify)
		{
		 	// both ends elastic
		 	eleForce = trialForce;
		 	return 0;
		}
		else
		{
			opserr << "InelasticYS2DGNL::predictor() - didn't think of this condition\n";

			opserr << "\a";
		}
	}
	
	// opserr << "RETURNING FROM PLASTIC PREDICTOR \n"; opserr << "\a";
	
	return 1;
}


// if an end plastifies, it can either shoot through or drift

void InelasticYS2DGNL::checkEndStatus(bool &end1drifts, bool &end2drifts, Vector &trialForce)
{
int driftOld;

	end1Plastify = false;
	end2Plastify = false;

	// check the status for end1
	int d1 = ys1->getTrialForceLocation(trialForce);
    if(d1 != INSIDE)
    {
    	end1Plastify = true;
    	
    	driftOld = ys1->getCommitForceLocation();
    	if(driftOld == INSIDE)
		{
		 	end1drifts = false;
		 	if(statusDebug)
		 		opserr << "checkEndStatus(..) ["<<getTag()<<"] - End 1 shoots through\n";
		}
		else if(driftOld == WITHIN)
		{
			// in this case if the current point is within or
    		// outside, eiher way, force point is taken as
    		// drifted
		 	end1drifts = true;
		 	if(statusDebug)
		 		opserr << "checkEndStatus(..) ["<<getTag()<<"] - End 1 drifts\n";
		 	
		}
		else
		{
		 	opserr << "WARNING - checkEndStatus end1 force_hist outside ["<<getTag()<<"]\n";
		 	//cin.get();
			//ys1->setToSurface(eleForce, ys1->ConstantYReturn);
		}
    }
    else
    {
    	if(statusDebug)
    	{
        	driftOld = ys1->getCommitForceLocation();
        	if(driftOld != INSIDE)
        	{
        		double drift = ys1->getTrialDrift(trialForce);
        		opserr << "checkEndStatus(..) ["<<getTag()<<"] - End 1 unloading "<<drift<<" \n";
			// opserr << "\a";
        	}
        	else
        		opserr << "checkEndStatus(..) ["<<getTag()<<"] - End 1 remains elastic\n";
    	}
    }

    // check the status for end2
	int d2 = ys2->getTrialForceLocation(trialForce);
    if(d2 != INSIDE)
    {
    	end2Plastify = true;
    	
    	driftOld = ys2->getCommitForceLocation();
    	if(driftOld == INSIDE)
		{
		 	end2drifts = false;

		 	if(statusDebug)
		 		opserr << "checkEndStatus(..) ["<<getTag()<<"] - End 2 shoots through\n";
		}
		else if(driftOld == WITHIN)
		{
		 	end2drifts = true;
		 	if(statusDebug)
		 		opserr << "checkEndStatus(..) ["<<getTag()<<"] - End 2 drifts\n";
		}
		else
		{


		 	opserr << "WARNING - checkEndStatus end2 force_hist outside ["<<getTag()<<"]\n";

		 	// opserr << "\a";
			// this condition could happen due to force balance at
			// privous step


		}
    }
    else
    {
    	if(statusDebug)
    	{
    		driftOld = ys2->getCommitForceLocation();
    		if(driftOld != INSIDE)
    		{
    			double drift = ys2->getTrialDrift(trialForce);
    			opserr << "checkEndStatus(..) ["<<getTag()<<"] - End 2 unloading "<<drift<<"\a \n";
			//cin.get();
        	}
        	else
        		opserr << "checkEndStatus(..) ["<<getTag()<<"] - End 2 remains elastic\n";

    	}
    }

//    if(statusDebug)
//    	opserr << "\a";
}


void InelasticYS2DGNL::forceBalance(Vector &eleforce, int algo)
{
    int sign1 = 1, sign2 = 1;
	if(eleforce(0) <0) sign1 = -1;
	if(eleforce(3) <0) sign2 = -1;
	double Pavg = (fabs(eleforce(0)) + fabs(eleforce(3)))/2;

	double Pmin = min_(fabs(eleforce(0)), fabs(eleforce(3)));
	double Pmax = max_(fabs(eleforce(0)), fabs(eleforce(3)));

    switch (algo)
	{
		case 1: // use Pavg
		{
			eleforce(0) = Pavg*sign1;
			eleforce(3) = Pavg*sign2;
			break;
		}
		case 2: // use Pmin
		{
			eleforce(0) = Pmin*sign1;
			eleforce(3) = Pmin*sign2;
			break;
		}

		case 3: // use Pmax
		{
			eleforce(0) = Pmax*sign1;
			eleforce(3) = Pmax*sign2;
			break;
		}
		default:
			opserr << "InelasticYS2DGNL::forceBalance - unkown algo\n";
			break;
	}


	if(end1Plastify)
	{
		ys1->setToSurface(eleforce, ys1->ConstantYReturn);

	}

	if(end2Plastify)
	{
		ys2->setToSurface(eleforce, ys2->ConstantYReturn);

	}


	// shear force balance
    eleforce(1) = (eleforce(2) + eleforce(5))/L;
	eleforce(4) = -eleforce(1);


}

const Matrix &InelasticYS2DGNL::getTangentStiff(void)
{
	// opserr << " getTangentStiff Called \n";

	if(!init)
	{
		this->update();
		init = true;
	}

//	if(forceRecoveryAlgo != -1) // do not use Kt
//		return this->Element2D02::getTangentStiff();
	
	transformToGlobal(Stiff);
	
//	opserr << "Kt = " << Stiff;
	
	return Stiff;
}



const Vector &InelasticYS2DGNL::getResistingForce(void)
{
	// opserr << "getResistingForce Called \n";

	if(!init)
	{
		//~ eleForce.Zero();
		this->update();
		init = true;
	}
    // check for quick return
    if (L == 0)
	    return ZeroVector;
	m_Iter++;

	// untested code !!
	bool needBal = false;
    if(ys1->hModel->freezeEvolution)
    {
		needBal = true;
		for(int i=0; i<3; i++)
			eleForce(i) = eleForce_hist(i);
	}
    if(ys2->hModel->freezeEvolution)
    {
		needBal = true;
		for(int i=3; i<6; i++)
			eleForce(i) = eleForce_hist(i);
	}
	if(needBal)
		this->forceBalance(eleForce, 1);
	
	// determine the ele end forces in global coords - want -F into rForce
    force(0) =  cs*eleForce(0) - sn*eleForce(1);
    force(1) =  sn*eleForce(0) + cs*eleForce(1);
    force(2) =  eleForce(2);
    force(3) =  cs*eleForce(3) - sn*eleForce(4);
    force(4) =  sn*eleForce(3) + cs*eleForce(4);
    force(5) =  eleForce(5);

    if(pdebug)
    {
    	opserr << "Returning Force \n";
    	opserr << force;
    }

    storage = 0;
    if(getTag()==1 || getTag() ==3)
    {
    	// opserr << "EleForce - global [" << getTag() << "]\n";

    	// for(int i = 0; i<6; i++)

    	// opserr << force[2] << " ";
    	storage += force[2];
    	//cerr << "\n";
    }

   // opserr << m_Iter << " RF = " << force; //cin.get();
//    if(getTag()==3)
//    	opserr << "\n" << storage;
	return force;
}




void InelasticYS2DGNL::plastifyOneEnd(int end, YieldSurface_BC *ys,  Vector &trial_force,
								Vector &incrDisp, Matrix &K, Vector &total_force, int algo)
{

// only one end is plastifying
// need to calculate lambda = G'Ke d(Del)_in / G'(Ke+Kp)G
//                          = G' dF_in / G'(Ke+Kp)G
// F_trial ( = F_hist + dF_in ) is outside the surface
// set F_trial to surface (modify dF_in)
// no need to re-balance the axial force or shear
// for 1 end plastifying G' is such that the numerator is uncoupled from
// the other end

	// opserr << "InelasticYS2DGNL::plastifyOneEnd " << endln;
	if(plastkDebug)
		opserr << "----------------------------------------------------------------------"
		     << endln;

Vector trialForce(6);
	// copy trial_force so that it does not get modified
	trialForce =  trial_force;

Vector surfaceForce(6);
Matrix G(6,1);
bool use_Kr = true;

// case: if it shoots through
//         use dF return to surface

	// driftnew is either outside or within (that's why plastify end was called)

	// either case, if previous point was elastic, implies a partly elastic load step
	// and remaining inelastic
	int driftOld = ys->getCommitForceLocation();

	if(driftOld == INSIDE)
	{
		use_Kr = false;		
		surfaceForce = trialForce;




		ys->setToSurface(surfaceForce, ys->RadialReturn);  //dFReturn, ConstantYReturn, RadialReturn
		// this part can get crazy!!
		// opserr << "    -------------- HERE? ------------ \n";
		// opserr << "surface drift = " << ys->getTrialDrift(surfaceForce) << endln; opserr << "\a";

		ys->getTrialGradient(G, surfaceForce);

		if(plastkDebug)
		{
			opserr << "Element (" << getTag() << ") plastifyEnd shoots through end: " << end << "\n";
		}
	}
	// Now we know that force point has drifted from the surface
	else if(driftOld != WITHIN)
	{
	 	opserr << "WARNING: InelasticYS2DGNL::plastifyOneEnd = " << end << " - driftOld outside [" << getTag()<<"]\n";
	 	opserr << "\a";
	}
	else
	{
		ys->getCommitGradient(G);
		surfaceForce =  eleForce_hist;

		if(plastkDebug)
			opserr << "Element (" << getTag() << ") plastifyEnd (" << end << ") drifts from surface\n";
	}

	/*if(driftOld==INSIDE)
	{
		surfaceForce = trial_force;
		ys->setToSurface(surfaceForce, ys->RadialReturn);
		total_force = surfaceForce;
		//cout << "returning surface force\n"; opserr << "\a";
		return;
	}*/

	
Vector dF_in(6);
	dF_in = trialForce - surfaceForce;
	
Matrix Ktp(6,6); //Ke(6,6); // want to leave 'K' unmodified for now

	int drift_test = ys->getTrialForceLocation(surfaceForce);
//	opserr << " force_hist = " << eleForce_hist;
//	opserr << " drift_test = " << drift_test << endln;

	Ktp = K; //includes Ke + Kg
	ys->addPlasticStiffness(Ktp);


	Matrix KI = G^(Ktp*G);

	double inv = 1/KI(0,0);

	Vector Lm = G^(dF_in);

	Lm = Lm*inv;
	double lamda = Lm(0);
	if(fabs(lamda) < ERROR) lamda = 0.0; // to get rid of -1e-15 etc

	if(lamda < 0)
	{
		//cout << "lamda = " << lamda << endln;
		//cin.get();
		use_Kr = false;
		lamda = 0.0;
	}

	Vector delP(6);
	delP(0) = G(0,0);
	delP(1) = G(1,0);

	delP(2) = G(2,0);
	delP(3) = G(3,0);
	delP(4) = G(4,0);
	delP(5) = G(5,0);
	
	delP = delP*lamda;	
int grow;
	if(algo != 20)	
		grow = ys->modifySurface(lamda, surfaceForce, G);
	else
	{
		grow = ys->modifySurface(lamda, surfaceForce, G, 1);
		use_Kr = false;
	}

	if(grow < 0)
		forceRecoveryAlgo = ys->ConstantYReturn;
	else
		forceRecoveryAlgo = forceRecoveryAlgo_orig;

Vector dF_t(6);
	dF_t = dF_in - K*delP;

	if(split_step)
		total_force = surfaceForce + dF_t; //otherwise convergence problem
	else
		total_force = surfaceForce + dF_in; // faster convergence
	

	// ys->displayForcePoint(total_force, 3);

	if(plastkDebug)
    {
		opserr << " InelasticYS2DGNL::plastifyOneEnd - tag = " << getTag() << "\n";

    	opserr << "lamda = " << lamda << "\n";
	opserr << " G = " << G << endln;
	opserr << "delP = " << delP <<endln;
		// opserr << "\a";
    }
	
	// now we update the stiffness
	if(algo == -10) 
		use_Kr = false;

	if(split_step)
		use_Kr = false;	
		
	// opserr << "Plastify One End, not using Kr\n";
		
	if(use_Kr)
	{
	Matrix Kr(6,6);
		Kr = K*G*(G^K)*inv;
		Stiff = Stiff - Kr;
		
		// F1 = (Stiff)*incrDisp + surfaceForce;
		// ys->displayForcePoint(F1, 3); 
		// opserr << "\a";
	}
	// else
	//	 opserr << "NOT USING Kr " << endln;

}


void InelasticYS2DGNL::splitStep(int end_shoot, YieldSurface_BC *ys_shoots, YieldSurface_BC *ys_drifts,
	               Vector &trial_force, Matrix &K, Vector &total_force)
{
// first we need to find a factor by which to split the step

//	opserr << " InelasticYS2DGNL::splitStep(" << end_shoot << ")\n"; // opserr << "\a";
	
	split_step = true;
	
Vector F1(6);
	F1 =  trial_force;
	ys_shoots->setToSurface(F1, ys1->dFReturn);
	
	// opserr << "trial force location " << ys_shoots->getTrialForceLocation(F1) << endln; opserr << "\a";
	
int Pi = 0, Mi = 2;
int p_shoot = 0, m_shoot = 2;
int p_drift = 3, m_drift = 5;

	if(end_shoot==2)
	{
	 	Pi = 3;  p_shoot = 3; p_drift = 0;
	 	Mi = 5;  m_shoot = 5; m_drift = 2;
	}
		
   	double num = pow((F1(Pi) - eleForce_hist(Pi)),2)
   				 + pow((F1(Mi) - eleForce_hist(Mi)),2);
					
   	num = sqrt(num);
		
   	double denom = pow((trial_force(Pi) - eleForce_hist(Pi)),2)
   				 + pow((trial_force(Mi) - eleForce_hist(Mi)),2);
	
	denom = sqrt(denom);

double t= num/denom;

// now drift the other end by the same amount
	Vector trialForce2(6), step_force(6);
	
	
	
	//cout << "Splitting step by - " << t << endl ;
	//cout << "Original trial force = " << trial_force;
	trialForce2 = eleForce_hist + t*(trial_force - eleForce_hist);
	//cout << "Drift one end by trial force = " <<  trialForce2;
	
Vector f_surface = eleForce_hist;
int count =0;
	// while(1)
	{
    	this->driftOneEnd(ys_drifts, trialForce2, f_surface, K, step_force);
    	

    	this->forceBalance(step_force, 1);    	
    	// opserr << "After force balance, force = " << step_force;    	
    	// opserr << "Drift both ends, trialForce = " <<  trial_force
    	// << "surface force = " << step_force;
    	
    	/* Vector df =   trial_force - step_force;
    	// opserr << "df = " << df;
    	
    	double ratio =  fabs(df(m_shoot)/df(m_drift));
    	// opserr << "ratio = " << ratio << endln; // opserr << "\a";
    	
    	if(ratio > 0.75)
    		break;
    	count++;
    	
    	if(count > 100)
    	{
    		opserr << "WARNING: InelasticYS2DGNL::splitStep - artificial hardening not converging\n";
    		break;
    	}
    	
    	f_surface = step_force;   */
	}
	
//	opserr << "Count = " << count << endln; // opserr << "\a";
// now drift both ends by the remainder amount	
	
	trialForce2 = step_force + (1-t)*(trial_force - eleForce_hist);
	
	if(ys1->getTrialForceLocation(trialForce2) <0)
		opserr << "oops - 1\n";
	if(ys2->getTrialForceLocation(trialForce2) <0)
		opserr << "oops - 2\n";	
	
	// check the status first!
//	this->driftBothEnds(trialForce2, step_force, K, total_force);
	this->driftBothEnds(trialForce2, step_force, K, eleForce);	

}// split step	


void InelasticYS2DGNL::driftOneEnd(YieldSurface_BC *ys, Vector &trialForce, Vector &surfaceForce,
					 Matrix &K, Vector &total_force)
{
	
Matrix G(6,1);	
	//!!
	// ys->setToSurface(surfaceForce, ys->ConstantYReturn);
	ys->getTrialGradient(G, surfaceForce);	
	
Vector dF_in(6);
	dF_in = trialForce - surfaceForce;
	
Matrix Ktp(6,6); 
	
	Ktp = K; //includes Ke + Kg
	ys->addPlasticStiffness(Ktp);

	Matrix KI = G^(Ktp*G);

	double inv = 1/KI(0,0);

	Vector Lm = G^(dF_in);

	Lm = Lm*inv;
	double lamda = Lm(0);
	if(fabs(lamda) < ERROR) lamda = 0.0; // to get rid of -1e-15 etc

	if(lamda < 0)
	{
		//cout << "lamda = " << lamda << endln;
		//cin.get();

		lamda = 0.0;
	}

	Vector delP(6);
	delP(0) = G(0,0);
	delP(1) = G(1,0);
	delP(2) = G(2,0);
	delP(3) = G(3,0);
	delP(4) = G(4,0);
	delP(5) = G(5,0);
	
	delP = delP*lamda;	
	
	int grow = ys->modifySurface(lamda, surfaceForce, G);

	if(grow < 0)
		forceRecoveryAlgo = ys->ConstantYReturn;
	else
		forceRecoveryAlgo = forceRecoveryAlgo_orig;

Vector dF_t(6);
	dF_t = dF_in - K*delP;

//	total_force = surfaceForce + dF_t;
	total_force = surfaceForce + dF_in;
	
	bool use_Kr = false;
	if(use_Kr)
	{
		
	Matrix Kr(6,6);
		Kr = K*G*(G^K)*inv;
	Matrix Ks = Stiff;
		Stiff = Ks - Kr;		
	}
}


void InelasticYS2DGNL::driftBothEnds(Vector &trialForce, Vector &surfaceForce,
						  Matrix &K, Vector &total_force)
{
Matrix G1(6,1), G2(6,1), G(6, 2);

	ys1->getTrialGradient(G1, surfaceForce);
	ys2->getTrialGradient(G2, surfaceForce);
	

	for(int i=0; i<6; i++)
    {
     	G(i, 0) = G1(i,0);
     	G(i, 1) = G2(i,0);
    }


	
Vector dF_in(6);
	dF_in = trialForce - surfaceForce;

Matrix Ktp(6,6); 
	
	Ktp = K; //includes Ke + Kg

	
	ys1->addPlasticStiffness(Ktp);
	ys2->addPlasticStiffness(Ktp);

	Matrix KI = G^(Ktp*G);
	
	Vector Lm(2);
	Lm(0) = G1(0,0)*dF_in(0) + G1(2,0)*dF_in(2);
	Lm(1) = G2(3,0)*dF_in(3) + G2(5,0)*dF_in(5);	
	
	//cout << "dF = "  << dF_in;
	//cout << "G  = "  << G;
	//cout << "Lm = " << Lm;
	//cout << "KI = " << KI;
	
	double mx = max_(Lm(0), Lm(1));
	double mn = min_(Lm(0), Lm(1));
	double avg = (mx + mn)/2;
	
//	 if(    fabs(mn/mx)*100 < 25)
	{
//		opserr << "BALANCING Lm"<<endl;
//		Lm(0) = mx;
//		Lm(1) = mx;
		
//	    if(Lm(0)==mn)
//	    	Lm(0) = 0.0;
//	    if(Lm(1)==mn)
//	    	Lm(1) = 0.0;		
	}
	
	Lm = Lm/KI;
	
	//cout << " Lamda = " << Lm; //cin.get();
	
	double lamda1 = Lm(0);
	double lamda2 = Lm(1);
	
	if(fabs(lamda1) < ERROR) lamda1 = 0.0; // to get rid of -1e-15 etc

	if(fabs(lamda2) < ERROR) lamda2 = 0.0; // to get rid of -1e-15 etc	

	if(lamda1 < 0 || lamda2 < 0)
	{
	 	//cout << "lamda1 = " << lamda1 << " lamda2 = " << lamda2 << endln;
	 	//cin.get();	
	 	

	 	if(lamda1 < 0.0) lamda1 = 0.0; //fabs(lamda1);
	 	if(lamda2 < 0.0) lamda2 = 0.0; //fabs(lamda2);
	}
	
	int grow1 = ys1->modifySurface(lamda1, surfaceForce, G1, 1);
	int grow2 = ys2->modifySurface(lamda2, surfaceForce, G2, 1);
	
Vector delP(6);
	delP(0) = G(0,0)*lamda1;
	delP(1) = G(1,0)*lamda1;
	delP(2) = G(2,0)*lamda1;
	delP(3) = G(3,1)*lamda2;
	delP(4) = G(4,1)*lamda2;
	delP(5) = G(5,1)*lamda2;

Vector dF_t(6);
	dF_t = dF_in - K*delP;

//	total_force = surfaceForce + dF_t;	
	total_force = surfaceForce + dF_in;
	
	// now we update the stiffness

	bool use_Kr = false;
	
	if(use_Kr)
	{
		Matrix Kr(6,6);
		// opserr << "  ----- USING KR -----" << endln;
		Matrix inv(2,2);

		inv(0,0) =    KI(1,1);
		inv(0,1) = -1*KI(0,1);
		inv(1,0) = -1*KI(1,0);
		inv(1,1) =    KI(0,0);

		double det = KI(0,0)*KI(1,1) - KI(1,0)*KI(0,1);

		if(fabs(det) < 1e-8) det = 1e-8;

		for(int i=0; i< 2; i++)

		{
			for(int j=0; j<2; j++)
			{
				inv(i,j) = inv(i,j)/det;
			}
		}

		Matrix Gn1 = G*inv;
		Matrix K1 = K*Gn1;
		Kr = K1*(G^K);
		

		Matrix Ks = Stiff;
		

		Stiff = Ks - Kr;
	}


}



// called only when force point at both the ends has drifted from the surface
void InelasticYS2DGNL::plastifyBothEnds(Vector &trial_force, Vector &incrDisp,
									 Matrix &K, Vector &total_force)
{

	// opserr << "InelasticYS2DGNL::plastifyBothEnds " << endln;
	if(plastkDebug)
		opserr << "----------------------------------------------------------------------"
		     << endln;

Vector trialForce(6);
	// copy trial_force so that it does not get modified
	trialForce =  trial_force;

Vector surfaceForce(6);
Matrix G1(6,1), G2(6,1), G(6, 2);
bool use_Kr = true;
bool end1drifts = true;
bool end2drifts = true;
	
	if(split_step)
		use_Kr = false;
int end = 1;
	
	// check for ys1 if the ysfp has shot through or drifted	
	// case: if it shoots through
	//         use dF return to surface
	// driftnew is either outside or within (that's why plastify end was called)
	// either case, if previous point was elastic, implies a partly elastic load step

	// and remaining inelastic
int driftOld = ys1->getCommitForceLocation();
	
	if(driftOld == INSIDE)
	{
		use_Kr = false;
		end1drifts = false;

		for(int i=0; i< 3; i++)
			surfaceForce(i) = trialForce(i);
		
			ys1->setToSurface(surfaceForce, ys1->RadialReturn);  //dFReturn, ConstantYReturn, RadialReturn
			ys1->getTrialGradient(G1, surfaceForce);

		if(plastkDebug)
		{
			opserr << "Element (" << getTag() << ") plastifyBothEnds shoots through end: " << end << "\n";
		}
	}
	// Now we know that force point has drifted from the surface
	else if(driftOld != WITHIN)
	{
	 	opserr << "WARNING: InelasticYS2DGNL::plastifyBothEnds = " << end << " - driftOld outside [" << getTag()<<"]\n";
	 	opserr << "\a";
	}
	else
	{
		ys1->getCommitGradient(G1);
		
		for(int i=0; i< 3; i++)		
			surfaceForce(i) =  eleForce_hist(i);

		if(plastkDebug)
			opserr << "Element (" << getTag() << ") plastifyBothEnds (" << end << ") drifts from surface\n";
	}
// end1 stuff over	
	
	
end = 2;
	
	// check for ys2 if the ysfp has shot through or drifted	
	// case: if it shoots through
	//         use dF return to surface
	// driftnew is either outside or within (that's why plastify end was called)
	// either case, if previous point was elastic, implies a partly elastic load step
	// and remaining inelastic
	driftOld = ys2->getCommitForceLocation();


	if(driftOld == INSIDE)
	{
		use_Kr = false;
		end2drifts = false;
		
		for(int i=3; i<6; i++)
			surfaceForce(i) = trialForce(i);

			ys2->setToSurface(surfaceForce, ys2->RadialReturn);  //dFReturn, ConstantYReturn, RadialReturn
		// this part can get crazy!!  - FIXED
		// actually, if dF return is used, trial force in different iterations can change
		// radial-return might be better to use		
	     // opserr << "surface drift = " << ys2->getTrialDrift(surfaceForce) << endln; opserr << "\a";

			ys2->getTrialGradient(G2, surfaceForce);

		if(plastkDebug)

		{
			opserr << "Element (" << getTag() << ") plastifyBothEnds shoots through end: " << end << "\n";
		}
	}
	// Now we know that force point has drifted from the surface
	else if(driftOld != WITHIN)
	{
	 	opserr << "WARNING: InelasticYS2DGNL::plastifyBothEnds = " << end << " - driftOld outside [" << getTag()<<"]\n";
	 	opserr << "\a";
	}
	else
	{
		ys2->getCommitGradient(G2);
		
		for(int i=3; i< 6; i++)		
			surfaceForce(i) =  eleForce_hist(i);

		if(plastkDebug)
			opserr << "Element (" << getTag() << ") plastifyBothEnds (" << end << ") drifts from surface\n";
	}
// end2 stuff over	
	
	
	/*if(!end1drifts && !end2drifts)
	{
		surfaceForce = trial_force;
		ys1->setToSurface(surfaceForce, ys1->ConstantYReturn);
		ys2->setToSurface(surfaceForce, ys2->ConstantYReturn);
		total_force = surfaceForce;
		//cout << "returning surface force\n"; opserr << "\a";
		return;
	}*/
	
	// now check that axial force is still the same, else take avg.



bool force_bal = false;
	
	if( fabs(surfaceForce(0)) != fabs(surfaceForce(3)))
	{
		//cout << "FORCE IMBALANCE \n";
		force_bal = true;

		this->forceBalance(surfaceForce, 1);
		
		ys1->setToSurface(surfaceForce, ys1->ConstantYReturn);
		ys2->setToSurface(surfaceForce, ys2->ConstantYReturn);
		
		ys1->getTrialGradient(G1, surfaceForce);		
		ys2->getTrialGradient(G2, surfaceForce);
	}
	
	for(int i=0; i<6; i++)
    {
     	G(i, 0) = G1(i,0);
     	G(i, 1) = G2(i,0);
    }
	
    //cout << "Matrix G = " << G;
	// opserr << "G1 = " << G1;
	// opserr << "G2 = " << G2;
	
// now to make the shear term consistent
// do axial force balancing later	
	
	// shear force balance - not needed now
    //~ surfaceForce(1) = (surfaceForce(2) + surfaceForce(5))/L;
	//~ surfaceForce(4) = -surfaceForce(1);
	
	
	
Vector dF_in(6);
	dF_in = trialForce - surfaceForce;
	
	/*if(force_bal)
	{
		opserr << "df_in before f_bal = " << dF_in;
		double avg = (dF_in(2) + dF_in(5)) / 2;
		opserr << "avg = " << avg << endln;
		
		dF_in(2) = avg;

		dF_in(5) = avg;


		
		// this->forceBalance(dF_in,1);
	}*/

	// opserr << "trialForce   = " << trialForce;
	// opserr << "surfaceForce = " << surfaceForce;
	// opserr << "dF_in        = " << dF_in << endln;

Matrix Ktp(6,6); //Ke(6,6); // want to leave 'K' unmodified for now

//	int drift_test = ys->getTrialForceLocation(surfaceForce);
//	opserr << " force_hist = " << eleForce_hist;
//	opserr << " drift_test = " << drift_test << endln;

	Ktp = K; //includes Ke + Kg

	
	if(end1drifts)
		ys1->addPlasticStiffness(Ktp);
	if(end2drifts)
		ys2->addPlasticStiffness(Ktp);

	Matrix KI = G^(Ktp*G);
	//Vector Lm = G^(dF_in);
	
	Vector Lm(2);
	Lm(0) = G1(0,0)*dF_in(0) + G1(2,0)*dF_in(2);
	Lm(1) = G2(3,0)*dF_in(3) + G2(5,0)*dF_in(5);	
		
	/*double mx = max_(Lm(0), Lm(1));
	double mn = min_(Lm(0), Lm(1));
	double avg = (mx + mn)/2;
	
	if(    (mn/mx)*100 < 50)
	{
		Lm(0) = mn;
		Lm(1) = mn;
	}*/
	
	//cout << "dF = " << dF_in;
	//cout << "Lm = " << Lm;
	//cout << "KI = " << KI;

	
	Lm = Lm/KI;
	
	//cout << " Lamda = " << Lm;  // opserr << "\a";
	
	double lamda1 = Lm(0);
	double lamda2 = Lm(1);
	
	
	// opserr << "KI = " << KI;
	// opserr << "lamda = " << lamda1 << " " << lamda2 << endln;


	
	if(fabs(lamda1) < ERROR) lamda1 = 0.0; // to get rid of -1e-15 etc
	if(fabs(lamda2) < ERROR) lamda2 = 0.0; // to get rid of -1e-15 etc	

	// happens when one end is shooting through and other drifting
	if(lamda1 < 0 || lamda2 < 0)
	{
		//~ opserr << "lamda = " << lamda1 << " " << lamda2 << endln;
		// opserr << "\a";
		use_Kr = false;
		
		//~ if(lamda1 < 0)
			//~ lamda1 = fabs(lamda1);
		//~ if(lamda2 < 0)
			//~ lamda2 = fabs(lamda1);
		
		if(lamda1 < 0)
			lamda1 = 0.0;
		if(lamda2 < 0)
			lamda2 = 0.0;

	}
			
	// Given the force point modify the surface
	int grow1 = ys1->modifySurface(lamda1, surfaceForce, G1);
	int grow2 = ys2->modifySurface(lamda2, surfaceForce, G2);
	
	// int factor = 1;
	//~ if(grow1 <= 0)


			//~ factor -= 0.5;
	//~ if(grow2 <= 0)
			//~ factor -= 0.5;

		
	// separate this also
	if(grow1 < 0 || grow2 < 0)
	{
		forceRecoveryAlgo = ys1->ConstantYReturn;
	}
	else // need to check every-time
		forceRecoveryAlgo = forceRecoveryAlgo_orig;
	
	Vector delP(6);
	delP(0) = G(0,0)*lamda1;
	delP(1) = G(1,0)*lamda1;
	delP(2) = G(2,0)*lamda1;
	delP(3) = G(3,1)*lamda2;
	delP(4) = G(4,1)*lamda2;
	delP(5) = G(5,1)*lamda2;


Vector dF_t(6);
	dF_t = dF_in - K*delP;

	total_force = surfaceForce + dF_t;
	
	//ys1->displayForcePoint(total_force, 3);
	//ys2->displayForcePoint(total_force, 3);
	
	if(plastkDebug)
    {
		opserr << " InelasticYS2DGNL::plastifyOneEnd - tag = " << getTag() << "\n";
    	opserr << "lamda = " << lamda1 << " " << lamda2 << "\n";
		opserr << " G = " << G << endln;
		opserr << "delP = " << delP <<endln;
		// opserr << "\a";
    }

	// now we update the stiffness

	Matrix Kr(6,6);

	//use_Kr = false;
	
	if(use_Kr)
	{
		Matrix inv(2,2);

		inv(0,0) =    KI(1,1);
		inv(0,1) = -1*KI(0,1);
		inv(1,0) = -1*KI(1,0);
		inv(1,1) =    KI(0,0);

		double det = KI(0,0)*KI(1,1) - KI(1,0)*KI(0,1);

		if(fabs(det) < 1e-8) det = 1e-8;


		for(int i=0; i< 2; i++)
		{
			for(int j=0; j<2; j++)
			{
				inv(i,j) = inv(i,j)/det;
			}
		}

		Matrix Gn1 = G*inv;
		Matrix K1 = K*Gn1;
		Kr = K1*(G^K);
		//Kr = K*G*(G^K)*inv;
		// opserr << "factor = " << factor << endln; opserr << "\a";
		
		Stiff = Stiff - Kr;
		
		// opserr << Stiff;
		
		// F1 = (Stiff)*incrDisp + surfaceForce;
		// ys->displayForcePoint(F1, 3); 
		// opserr << "\a";
	}
//	else
//		opserr << "NOT USING Kr " << endln;
	
	
}




int InelasticYS2DGNL::commitState()
{

	 if(pdebug)
		opserr << " ############# commit ############ ["<< getTag() << "]\n";
	//cin.get();
	/*bool end1drifts, end2drifts;
	this->checkEndStatus(end1drifts, end2drifts, eleForce);

	if(plastkDebug)
	{
		opserr << "End1Plastify = " << end1Plastify << ", End2Plastify = " << end2Plastify << "\n";
		opserr << "End1Drifts   = " << end1drifts   << ", End2Drifts   = " << end2drifts   << "\n";
		opserr << " ------------------------------------\n\n";
	}*/
	//cin.get();
//	if(getTag()==3)
//	{
//		opserr << eleForce;
//		opserr << "\a";
//	}

    split_step = false;
	this->UpdatedLagrangianBeam2D::commitState();

	if(end1Plastify)
		end1Damage = true;

	if(end2Plastify)
		end2Damage = true;

	//!! eleForce for shoot through goes outside
	ys1->commitState(eleForce);
	ys2->commitState(eleForce);

	end1Plastify_hist =  end1Plastify;
	end2Plastify_hist =  end2Plastify;

	if(pView)
	{
    	pView->clearImage();
    	pView->startImage();
    	ys1->displaySelf(*pView, 1, 1);
    	ys2->displaySelf(*pView, 1, 1);
    	pView->doneImage();
	}

	// opserr << "--- commit ---\n"; opserr << "\a";

	//if(getTag()==3) opserr << storage << "000000";

	return 0;
}



//////////////////////////////////////////////////////////////////////
// Pure virtual abstract methods
//////////////////////////////////////////////////////////////////////

void InelasticYS2DGNL::getLocalMass(Matrix &M)
{
	if(massDof < 0)
    {
		opserr << "Element2dGNL::getMass - Distributed mass not implemented\n";
	    M.Zero();
    }
    else if(massDof == 0)//this cond. is taken care of already
    {
        M.Zero();
    }
    else
    {
        M.Zero();
	    M(0,0) = M(1,1) = M(2,2) = M(3,3) = M(4,4) = M(5,5) = massDof;
    }

}

/*
void InelasticYS2DGNL::getLocalStiff(Matrix &K)
{
 double	EIbyL = E*Iz/L;

    K(0, 1) = K(0, 2) = K(0, 4) = K(0, 5)=0;
    K(1, 0) = K(1, 3) =0;
    K(2, 0) = K(2, 3) =0;
    K(3, 1) = K(3, 2) = K(3, 4) = K(3, 5)=0;
    K(4, 0) = K(4, 3) =0;
    K(5, 0) = K(5, 3) =0;

  	K(0,0) = K(3,3) = (A/Iz)*(EIbyL);
  	K(0,3) = K(3,0) = (-A/Iz)*(EIbyL);
  	K(1,1) = K(4,4) = (12/(L*L))*(EIbyL);
  	K(1,4) = K(4,1) = (-12/(L*L))*(EIbyL);
  	K(1,2) = K(2,1) = K(1,5) = K(5,1) = (6/L)*(EIbyL);
  	K(2,4) = K(4,2) = K(4,5) = K(5,4) = (-6/L)*(EIbyL);
  	K(2,2) = K(5,5) = 4*(EIbyL);
  	K(2,5) = K(5,2) = 2*(EIbyL);


}//getLocalStiff
*/

//////////////////////////////////////////////////////////////////////
// Print/Render  Send/Recv
//////////////////////////////////////////////////////////////////////
void InelasticYS2DGNL::createView(char *title, double scale, int x, int y, int cx, int cy, char displaytype)
{
	displayType = displaytype;


#ifdef _NOGRAPHICS

#else
#ifdef _GLX // Boris Jeremic added 23Oct2002
theMap = new PlainMap();
pView =  new OpenGLRenderer(title, x, y, cx, cy, *theMap);

 if(pView){
   pView->setVRP(0.0, 0.0, 0.0);
   pView->setVPN(0.0, 0.0, 1.0);
   pView->setVUP(0.0, 1.0, 0.0);
   pView->setFillMode("wire");             // wire mode
   pView->setPlaneDist(1.0, -1.0);
   pView->setPRP(0.0, 0.0, 10.0);
   pView->setPortWindow(-1, 1, -1, 1);  // use the whole window
   
   pView->setViewWindow(-scale, scale, -scale, scale);
   
   pView->clearImage();
   pView->startImage();
   

   ys1->setView(pView);
   ys2->setView(pView);
   
   ys1->displaySelf(*pView, 10, 1);
   ys2->displaySelf(*pView, 10, 1);
pView->doneImage();
 
 }
 else
   opserr << "WARNING: InelasticYS2DGNL::createView - Renderer not available\n";
#endif   // Boris Jeremic added 23Oct2002
#endif

}



/*
void InelasticYS2DGNL::createView(char *title, WindowManager *theWM, char displaytype)
{
	displayType = displaytype;
	theMap = new PlainMap();
	pView =	theWM->getRenderer(title, *theMap);
	//pView = new OpenGLRenderer(title, *theMap);
	if(!pView)
		opserr << "WARNING: InelasticYS2DGNL::createView - Renderer not available\n";
	//theWM->setRenderer(*pView);	
		
	if(pView)
	{
    	pView->setViewWindow(-3.5, 3.5, -3.5, 3.5);

        pView->clearImage();
    	pView->startImage();

    	ys1->setView(pView);
		ys2->setView(pView);

    	ys1->displaySelf(*pView, 1, 1);
    	ys2->displaySelf(*pView, 1, 1);
    	pView->doneImage();
	}
}
*/

void InelasticYS2DGNL::Print(OPS_Stream &s, int flag)
{
    s << "\nElement No: " << this->getTag();
    s << " type: InelasticYS2DGNL  iNode: " << connectedExternalNodes(0);
    s << " jNode: " << connectedExternalNodes(1);
    //s << "\nElement Forces ";
}

int InelasticYS2DGNL::sendSelf(int commitTag, Channel &theChannel)
{
	return -1;
}

int InelasticYS2DGNL::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
	return -1;
}


Response* InelasticYS2DGNL::setResponse(const char **argv, int argc, Information &eleInformation)
{
Response *suResponse=0;

	suResponse = this->UpdatedLagrangianBeam2D::setResponse(argv, argc, eleInformation);

	if(suResponse != 0)
		return suResponse;

    if (strcmp(argv[0],"ysVisual") == 0)
	{
		suResponse =  new ElementResponse(this, DISPLAY_YS);
    }

	return suResponse;
}


int InelasticYS2DGNL::getResponse(int responseID, Information &eleInformation)
{
int res = this->UpdatedLagrangianBeam2D::getResponse(responseID, eleInformation);

	if(res != -1)
		return res;

	if(responseID == DISPLAY_YS)
		res = responseID;

	return res;
}



int InelasticYS2DGNL::displaySelf(Renderer &theViewer, int displayMode, float fact)
{
    // first determine the two end points of the element based on
    // the display factor (a measure of the distorted image)
    // store this information in 2 3d vectors v1 and v2
    //cerr << "Inside display self mode: " << displayMode << " fact " << fact <<"\n";

	if(displayMode == DISPLAY_YS)
	{
		ys1->setView(&theViewer);
		ys2->setView(&theViewer);

		ys1->displaySelf(theViewer, 1, 1);
		ys2->displaySelf(theViewer, 1, 1);
		return 0;
	}


  this->UpdatedLagrangianBeam2D::displaySelf(theViewer, displayMode, fact);

  const Vector &end1Crd = end1Ptr->getCrds();
  const Vector &end2Crd = end2Ptr->getCrds();
  const Vector &end1Disp = end1Ptr->getTrialDisp();
  const Vector &end2Disp = end2Ptr->getTrialDisp();

	Vector v1(3);
    Vector v2(3);
	Vector vc(3);

	Vector rgb(3);
	rgb(0) = 0;
	rgb(1) = 0.9;
	rgb(2) = 0;

    for (int i=0; i<2; i++) { //!! i < 3
	v1(i) = end1Crd(i)+end1Disp(i)*fact;
	v2(i) = end2Crd(i)+end2Disp(i)*fact;
    }
	double e = 0.05;

	if (displayMode == 1) //theViewer.drawLine(v1, v2, 1, 1);
	{
		if(end1Damage && !end1Plastify)
		{
			vc(2) = v1(2);
			vc(0) = v1(0) + e*(v2(0) - v1(0));


			vc(1) = v1(1) + e*(v2(1) - v1(1));

			theViewer.drawPoint(vc, rgb, 3);
		}

		if(end2Damage && !end2Plastify)
		{
			vc(2) = v2(2);
			vc(0) = v2(0) + e*(v1(0) - v2(0));
			vc(1) = v2(1) + e*(v1(1) - v2(1));

			theViewer.drawPoint(vc, rgb, 3);
		}


		if(!end1Plastify && !end2Plastify) return 0;

		rgb(0) = 1;
		rgb(1) = 0;
		rgb(2) = 0;
		if(end1Plastify)
		{
			vc(2) = v1(2);

			vc(0) = v1(0) + e*(v2(0) - v1(0));
			vc(1) = v1(1) + e*(v2(1) - v1(1));

			theViewer.drawPoint(vc, rgb, 3);
		}
		if(end2Plastify)
		{
			vc(2) = v2(2);
			vc(0) = v2(0) + e*(v1(0) - v2(0));
			vc(1) = v2(1) + e*(v1(1) - v2(1));

			theViewer.drawPoint(vc, rgb, 3);

		}
	}
	return 0;

}




/* if(split_step && !end1Plastify)
	{
	 	ys1->setToSurface(trialForce, ys1->ConstantYReturn);
	 	trialForce(2) = trialForce(2)*1.00001;
	    this->forceBalance(trialForce, 1);
	}
	
	if(split_step && !end2Plastify)
	{
	 	ys2->setToSurface(trialForce, ys1->ConstantYReturn);
	 	trialForce(5) = trialForce(5)*1.00001;
	    this->forceBalance(trialForce, 1);
	}
    this->checkEndStatus(end1Drifts, end2Drifts, trialForce);
*/






/*
// called only when force point at both the ends has drifted from the surface
void InelasticYS2DGNL::plastifyBothEnds(Vector &trial_force, Vector &incrDisp,
									 Matrix &K, Vector &total_force)
{
Vector trialForce(6);

	// copy trial_force so that it does not get modified
	trialForce =  trial_force;
    Matrix G(6, 2);
    Matrix Kt(6,6), Ktp(6,6), Ke(6,6);

    Ke 	=  K;
    Ktp	= Ke;
    Kt 	= Ktp;

	ys1->addPlasticStiffness(Ktp);
	ys2->addPlasticStiffness(Ktp);
	
	Matrix ysG1(6,1), ysG2(6,1);
	Vector elasticForce(6);

	elasticForce =  eleForce_hist;
	
//	int driftOld = ys1->getTrialForceLocation(elasticForce);
//    if(driftOld == OUTSIDE)
//	{
//	 	opserr << "WARNING: InelasticYS2DGNL::plastifyBothEnds - end1 driftOld outside\n";

//	 	opserr << "\a";
//	}
	// In a trial state, elasticForce may not exactly be on surface, esp. in case
	// of hardening - it may actually be inside, from second trial step onwards
	// since phi = phi_hist + dR
	// set elasticForce to surface anyway, there might a minute difference in G,
	// but computation of phi will be accurate, and the force-balance step will
	// take care of any difference in F_total ( = F_on_surface + Kt*ddel)
	
	//!! ys1->setToSurface(elasticForce, ys1->ConstantYReturn);
	// this causes convergence problems, especially if ys is shrinking
	
//	driftOld = ys2->getTrialForceLocation(elasticForce);
//    if(driftOld == OUTSIDE)
//	{
//	 	opserr << "WARNING: InelasticYS2DGNL::plastifyBothEnds - end2 driftOld outside\n";
//	 	opserr << "\a";
//	}
	//!! ys2->setToSurface(elasticForce, ys2->ConstantYReturn);
	
    ys1->getCommitGradient(ysG1);
    ys2->getCommitGradient(ysG2);

    for(int i=0; i<6; i++)
    {
     	G(i, 0) = ysG1(i,0);
     	G(i, 1) = ysG2(i,0);
    }

	Matrix KI = G^(Ktp*G);
	Matrix inv(2,2);

	inv(0,0) =    KI(1,1);
	inv(0,1) = -1*KI(0,1);
	inv(1,0) = -1*KI(1,0);
	inv(1,1) =    KI(0,0);

	double det = KI(0,0)*KI(1,1) - KI(1,0)*KI(0,1);

	if(fabs(det) < 1e-8) det = 1e-8;

	for(int i=0; i< 2; i++)
	{
		for(int j=0; j<2; j++)
		{
			inv(i,j) = inv(i,j)/det;
		}
	}

	// Calculate modified stiffness
	Matrix G1 = G*inv;
	Matrix K1 = Ke*G1;
	Matrix Kr = K1*(G^Ke);
	Kt = Ke - Kr;

	Vector din = incrDisp;
	total_force = Kt*din;


	// add to the elastic force point on the surface
	total_force +=  elasticForce;

	Matrix dinc(6,1);
	for(int i=0; i<6; i++)
		dinc(i,0) = din(i);

	Matrix Lm = inv*(G^Ke)*dinc;

	double lamda1 = Lm(0,0), lamda2 = Lm(1,0);
    if(fabs(lamda1) < ERROR) lamda1 = 0; // to get rid of -1e-15 etc
    if(fabs(lamda2) < ERROR) lamda2 = 0; // to get rid of -1e-15 etc

	// Given the force points drifted, modify the surface
	//~ int grow1 = ys1->modifySurface(lamda1, elasticForce);
	//~ int grow2 = ys2->modifySurface(lamda2, elasticForce);

	//~ if(grow1 < 0 || grow2 < 0)
		//~ forceRecoveryAlgo = ys1->ConstantYReturn;
    //~ else
		//~ forceRecoveryAlgo = forceRecoveryAlgo_orig;


}
*/

