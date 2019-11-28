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

// Written: Saumyasuchi Das
//          May 2013, University of Canterbury
// Updated: Chris McGann
//          June 2015, Washington State University

#include <StressDensityModel2D.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Information.h>
#include <Parameter.h>

// full constructor
StressDensityModel2D::StressDensityModel2D(int tag, double constDensity, double initialVoidRatio, double constA,
                                               double exponentN, double poissonRatio, double constAlpha1, 
                                               double constBeta1, double constAlpha2, double constBeta2, 
                                               double constAlpha3, double constBeta3, double constDegradation,
                                               double constMumin, double constMucyclic,	double constDilatancyStrain,
                                               double constMumax, double constPatm, double constsslvoidatP1,
                                               double constsslvoidatP2, double constsslvoidatP3, double constsslvoidatP4,
                                               double constsslvoidatP5, double constsslvoidatP6, double constsslvoidatP7,
                                               double constsslvoidatP8, double constsslvoidatP9, double constsslvoidatP10,
                                               double consthslvoid, double constP1, double constP2, double constP3, 
                                               double constP4, double constP5, double constP6, double constP7, 
                                               double constP8, double constP9, double constP10, double constRxx,
                                               double constRyy, double constRxy)
  : StressDensityModel(tag, ND_TAG_StressDensityModel2D, constDensity, initialVoidRatio, constA, exponentN,	poissonRatio,
                          constAlpha1, constBeta1, constAlpha2, constBeta2, constAlpha3, constBeta3, constDegradation,
                          constMumin, constMucyclic, constDilatancyStrain, constMumax, constPatm, constsslvoidatP1,
                          constsslvoidatP2, constsslvoidatP3, constsslvoidatP4, constsslvoidatP5, constsslvoidatP6,
                          constsslvoidatP7, constsslvoidatP8, constsslvoidatP9, constsslvoidatP10, consthslvoid,
                          constP1, constP2, constP3, constP4, constP5, constP6, constP7, constP8, constP9, constP10,
                          constRxx, constRyy, constRxy),
    stressCurrent(3),
    stressNext(3),
    strainCurrent(3),
    strainNext(3),
    initialTangent(3,3),
    currentTangent(3,3)
{
    // initialise variables
    this->initialise();
    
    // get the initial material tangent
	this->CalInitialTangent();
    // set current tangent as initial to start
	currentTangent = initialTangent;
}

// null constructor
StressDensityModel2D::StressDensityModel2D()
  : StressDensityModel(),
    stressCurrent(3),
    stressNext(3),
    strainCurrent(3),
    strainNext(3),
    initialTangent(3,3),
    currentTangent(3,3)
{
    this->initialise();
}

// destructor
StressDensityModel2D::~StressDensityModel2D()
{
}

NDMaterial *
StressDensityModel2D::getCopy(void)
{
	/*StressDensityModel2D *theCopy = 0;
	if(theCopy){
		*theCopy=*this;
		return theCopy;
	} else {
		opserr<<"StressDensityModel2D::getCopy failed to get copy: " << endln;
	  	return 0;
	}*/
    StressDensityModel2D *theCopy;
    theCopy = new StressDensityModel2D();
    *theCopy = *this;
    return theCopy;

}

const char *
StressDensityModel2D::getType(void) const 
{
	return "PlaneStrain";
}

int 
StressDensityModel2D::getOrder( ) const 
{
	return 3;
}

int 
StressDensityModel2D::setTrialStrain(const Vector &strain_from_element) 
{
    strainNext = strain_from_element;
	this->GetCurrentStress(); // calculates the trial stress and current tangent
	return 0;
}

int 
StressDensityModel2D::setTrialStrain(const Vector &v, const Vector &r)
{
	return this->setTrialStrain(v);
}

const Matrix &
StressDensityModel2D:: getTangent(void) 
{
    return currentTangent;
}

const Matrix &
StressDensityModel2D::getInitialTangent(void) 
{
    this->CalInitialTangent();
    return initialTangent;
}

const Vector &
StressDensityModel2D::getStress(void)
{
    return stressNext;
}

const Vector &
StressDensityModel2D::getStrain(void)
{
    return strainCurrent;
}

int 
StressDensityModel2D::commitState(void) 
{
    //Commit stress and strain
	stressCurrent=stressNext;
	strainCurrent=strainNext;

	//Commit hardening parameters
	for(int i=0;i<(Nsurface*7+5);i++) {
        hard_para_real[i] = _hard_para_real[i]; 
    }
	for(int i=0;i<2;i++) {
        hard_para_int[i]= _hard_para_int[i];
    }

	return 0;
}

int 
StressDensityModel2D::sendSelf(int commitTag, Channel &theChannel)
{
    int res = 0;

    // place data in a vector
    static Vector vData(318);

	vData(0)  = this->getTag();
	vData(1)  = theStage;
	vData(2)  = hard_para_int[0];
	vData(3)  = hard_para_int[1];
    vData(4)  = theDensity;
    vData(5)  = modelParameter[0];
    vData(6)  = modelParameter[1];
    vData(7)  = modelParameter[2];
    vData(8)  = modelParameter[3];
    vData(9)  = modelParameter[4];
    vData(10) = modelParameter[5];
    vData(11) = modelParameter[6];
    vData(12) = modelParameter[7];
    vData(13) = modelParameter[8];
    vData(14) = modelParameter[9];
    vData(15) = modelParameter[10];
    vData(16) = modelParameter[11];
    vData(17) = modelParameter[12];
    vData(18) = modelParameter[13];
    vData(19) = modelParameter[14];
    vData(20) = modelParameter[15];

    vData(21) = sslVoidratio[0]; vData(31) = hslVoidratio[0]; vData(41) = refPressure[0];
    vData(22) = sslVoidratio[1]; vData(32) = hslVoidratio[1]; vData(42) = refPressure[1];
    vData(23) = sslVoidratio[2]; vData(33) = hslVoidratio[2]; vData(43) = refPressure[2];
    vData(24) = sslVoidratio[3]; vData(34) = hslVoidratio[3]; vData(44) = refPressure[3];
    vData(25) = sslVoidratio[4]; vData(35) = hslVoidratio[4]; vData(45) = refPressure[4];
    vData(26) = sslVoidratio[5]; vData(36) = hslVoidratio[5]; vData(46) = refPressure[5];
    vData(27) = sslVoidratio[6]; vData(37) = hslVoidratio[6]; vData(47) = refPressure[6];
    vData(28) = sslVoidratio[7]; vData(38) = hslVoidratio[7]; vData(48) = refPressure[7];
    vData(29) = sslVoidratio[8]; vData(39) = hslVoidratio[8]; vData(49) = refPressure[8];
    vData(30) = sslVoidratio[9]; vData(40) = hslVoidratio[9]; vData(50) = refPressure[9];

    vData(51) = refOrigin[0]; vData(54) = stressCurrent(0); vData(57) = strainCurrent(0);
    vData(52) = refOrigin[1]; vData(55) = stressCurrent(1); vData(58) = strainCurrent(1);
    vData(53) = refOrigin[2]; vData(56) = stressCurrent(2); vData(59) = strainCurrent(2);

    vData(60) = currentTangent(0,0); vData(63) = currentTangent(0,1); vData(66) = currentTangent(0,2);
    vData(61) = currentTangent(1,0); vData(64) = currentTangent(1,1); vData(67) = currentTangent(1,2);
    vData(62) = currentTangent(2,0); vData(65) = currentTangent(2,1); vData(68) = currentTangent(2,2);

    for (int i=0;i<Nsurface*7+5;i++) {
        vData(69+i) = hard_para_real[i];
    }

	res = theChannel.sendVector(this->getDbTag(), commitTag, vData);
	if (res < 0) {
      opserr << "StressDensityModel::sendSelf() - failed to send vData\n";
	  return -1;
	}

	return 0;
}

int 
StressDensityModel2D::recvSelf(int commitTag, Channel &theChannel,FEM_ObjectBroker &theBroker)
{
    int res = 0;

    // place data in a vector
    static Vector vData(318);

	res = theChannel.recvVector(this->getDbTag(), commitTag, vData);
	if (res < 0) {
		opserr << "StressDensityModel::recvSelf() - failed to recv vData\n";
		return -1;
    }
	
    this->setTag((int)vData(0));
	theStage           = (int)vData(1);
	hard_para_int[0]   = (int)vData(2);
	hard_para_int[1]   = (int)vData(3);
    theDensity         = vData(4);
    modelParameter[0]  = vData(5);
    modelParameter[1]  = vData(6);
    modelParameter[2]  = vData(7);
    modelParameter[3]  = vData(8);
    modelParameter[4]  = vData(9);
    modelParameter[5]  = vData(10);
    modelParameter[6]  = vData(11);
    modelParameter[7]  = vData(12);
    modelParameter[8]  = vData(13);
    modelParameter[9]  = vData(14);
    modelParameter[10] = vData(15);
    modelParameter[11] = vData(16);
    modelParameter[12] = vData(17);
    modelParameter[13] = vData(18);
    modelParameter[14] = vData(19);
    modelParameter[15] = vData(20);

    sslVoidratio[0] = vData(21); hslVoidratio[0] = vData(31); refPressure[0] = vData(41);
    sslVoidratio[1] = vData(22); hslVoidratio[1] = vData(32); refPressure[1] = vData(42);
    sslVoidratio[2] = vData(23); hslVoidratio[2] = vData(33); refPressure[2] = vData(43);
    sslVoidratio[3] = vData(24); hslVoidratio[3] = vData(34); refPressure[3] = vData(44);
    sslVoidratio[4] = vData(25); hslVoidratio[4] = vData(35); refPressure[4] = vData(45);
    sslVoidratio[5] = vData(26); hslVoidratio[5] = vData(36); refPressure[5] = vData(46);
    sslVoidratio[6] = vData(27); hslVoidratio[6] = vData(37); refPressure[6] = vData(47);
    sslVoidratio[7] = vData(28); hslVoidratio[7] = vData(38); refPressure[7] = vData(48);
    sslVoidratio[8] = vData(29); hslVoidratio[8] = vData(39); refPressure[8] = vData(49);
    sslVoidratio[9] = vData(30); hslVoidratio[9] = vData(40); refPressure[9] = vData(50);

    refOrigin[0] = vData(51); stressCurrent(0) = vData(54); strainCurrent(0) = vData(57);
    refOrigin[1] = vData(52); stressCurrent(1) = vData(55); strainCurrent(1) = vData(58);
    refOrigin[2] = vData(53); stressCurrent(2) = vData(56); strainCurrent(2) = vData(59);

    currentTangent(0,0) = vData(60); currentTangent(0,1) = vData(63); currentTangent(0,2) = vData(66);
    currentTangent(1,0) = vData(61); currentTangent(1,1) = vData(64); currentTangent(1,2) = vData(67);
    currentTangent(2,0) = vData(62); currentTangent(2,1) = vData(65); currentTangent(2,2) = vData(68);

    for (int i=0;i<Nsurface*7+5;i++) {
        hard_para_real[i] = vData(69+i);
    }
	
	return 0;
}

// ******************************* PRIVATE METHODS **********************************

void
StressDensityModel2D::initialise() 
{
    // initialise hardening parameters to zeros
    for (int i=0;i<7*Nsurface+5;i++) {
        hard_para_real[i]=0.0;
    }
	for (int i=0;i<2;i++) {
        hard_para_int[i]=0;
    }

    // initialise Vector and Matrix variables
    stressCurrent.Zero();
    stressNext.Zero();
    strainCurrent.Zero();
    strainNext.Zero();
    initialTangent.Zero();
    currentTangent.Zero();

    // initialise FORTRAN in/out variables
    for (int i=0;i<16;i++) {
        _model_parameter[i] = 0.0;
    }
    for (int i=0;i<10;i++) {
        _ssl_void_ratio[i] = 0.0;
        _hsl_void_ratio[i] = 0.0;
        _ssl_pressure[i] = 0.0;
        _hsl_pressure[i] = 0.0;
    }
    for (int i=0;i<3;i++) {
        _stress_current[i] = 0.0;
        _strain_current[i] = 0.0;
        _strain_next[i] = 0.0;
    }
    for (int i=0;i<7*Nsurface+5;i++) {
        _hard_para_real[i]=0.0;
    }
	for (int i=0;i<2;i++) {
        _hard_para_int[i]=0;
    }
}

/*-----------------------------------------------------------------------------------
This method calculates current stress from a given current strain but does not update
the model. The intial state remains the same during iteration. The intial state is
stored in class members. 

The class members are nor directly passed to the FORTRAN subroutine. They are first
copied in temporary variables (of static construct and hence are the same for all 
instances of this class) and then passed to FORTRAN subroutine.

If the state is committed, the class members are updated from these temporary
variables; else not.

The temporary variables have an underscore as the prefix.
-----------------------------------------------------------------------------------*/
void
StressDensityModel2D::GetCurrentStress(void){

	// Elastic 
	// ---------

	if (theStage!=1){ 
		stressNext=stressCurrent+currentTangent*(strainNext-strainCurrent);
		return;
	}

	// Plastic
	// ---------
	
	// Copy the temporary variables

	// Change the sign of the input data from OPENSEES
	// Compressive (normal) stress and strain is +ve

	for(int i=0;i<2;i++)_stress_current[i] = -stressCurrent(i); 
	for(int i=2;i<3;i++)_stress_current[i] =  stressCurrent(i); 

	for(int i=0;i<2;i++)_strain_current[i] = -strainCurrent(i);  
	for(int i=2;i<3;i++)_strain_current[i] =  strainCurrent(i)/2.;	//convert to true strain

	for(int i=0;i<2;i++)_strain_next[i] = -strainNext(i);
	for(int i=2;i<3;i++) {
        _strain_next[i] =  strainNext(i)/2.;
    }

	for(int i=0;i<16;i++) {
        _model_parameter[i] = modelParameter[i];
    }

	for(int i=0;i<10;i++) {
        _ssl_void_ratio[i] = sslVoidratio[i];
    }

	for(int i=0;i<10;i++) {
        _ssl_pressure[i] = refPressure[i];
    }

	for(int i=0;i<10;i++) {
        _hsl_void_ratio[i] = hslVoidratio[i];
    }

	for(int i=0;i<10;i++) {
        _hsl_pressure[i] = refPressure[i];
    }

	for(int i=0;i<Nsurface*7+5;i++) {
        _hard_para_real[i] = hard_para_real[i];
    }

	for(int i=0;i<2;i++) {
        _hard_para_int[i] = hard_para_int[i];
    }

	// in Fortran the double-scripted arrays will be transposed
	for(int i=0;i<3;i++) {
		for(int j=0;j<3;j++) _tangent[j][i]=currentTangent(i,j);
    }
	
	//Fortran subroutine for stress integration...
    sdm2d_(_stress_current, 
	       _strain_current, 
	       _strain_next, 
	       _model_parameter, 
	       _ssl_void_ratio, 
	       _ssl_pressure, 
	       _hsl_void_ratio, 
	       _hsl_pressure, 
	       _hard_para_real,
	       _hard_para_int,
	       _tangent);

	//Update the stress variable only as the state is not committed
	for(int i=0;i<2;i++) stressNext(i) = -_stress_current[i];
	for(int i=2;i<3;i++) stressNext(i) =  _stress_current[i];
		
	// in Fortran the double-scripted arrays will be transposed
	// in this section tangent is calculated for the Full-Newton-Raphson scheme
	for(int i=0;i<3;i++) {
		for(int j=0;j<3;j++) currentTangent(i,j) = _tangent[j][i];
    }

}

void 
StressDensityModel2D::CalInitialTangent(void)
{
    double nu, G, A, m, eo, patm;

    A    = modelParameter[1];
    m    = modelParameter[2];
    eo   = modelParameter[0];
    patm = modelParameter[15];

    // assume p = patm for initial shear modulus
	G = A*patm*(2.17 - eo)*(2.17 - eo)/(1 + eo)*(pow(1,m));
	nu = modelParameter[3]; 

	initialTangent(0,0) = 2*G/(1-2*nu)*(1-nu);
	initialTangent(0,1) = 2*G/(1-2*nu)*nu;
	initialTangent(0,2) = 0;


	initialTangent(1,0) = 2*G/(1-2*nu)*nu;
	initialTangent(1,1) = 2*G/(1-2*nu)*(1-nu);
	initialTangent(1,2) = 0;
	

	initialTangent(2,0) = 0;
	initialTangent(2,1) = 0;
	initialTangent(2,2) = G;
	
}
