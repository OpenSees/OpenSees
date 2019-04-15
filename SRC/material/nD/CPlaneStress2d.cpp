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
** ****************************************************************** */
                                                                        
// $Revision: 1.7 $
// $Date: 2008/10/20 22:23:03 $
// Written: Tesser, Talledo

//  Send strains in following format :
// 
//     strain_vec = {   eps_00
//                      eps_11
//                    2 eps_01   }   <--- note the 2

#include <CPlaneStress2d.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>

//static vectors and matrices
Vector CPlaneStress2d :: strain_vec(3) ;
Vector CPlaneStress2d :: stress_vec(3) ;
Matrix CPlaneStress2d :: tangent_matrix(3,3) ;
Matrix CPlaneStress2d :: InitialStiffness2d(3,3);

void CPlaneStress2d :: zero ( ) 
{
	rnp   = r0p;
	rnn   = r0n;
	rnp1p = rnp;
	rnp1n = rnn;
	pStrain2d_n.Zero();
	pStrain2d_n1.Zero();
}

//null constructor
CPlaneStress2d ::  CPlaneStress2d( ) : 
Concrete( ),
stress2d(3),
strain2d(3),
pStrain2d_n(3),
pStrain2d_n1(3),
Stiffness2d(3,3)
{
	E    = 0.0;
	ni   = 0.0;
	f01d = 0.0;
	f02d = 0.0;
	f0p  = 0.0;
	beta = 0.0;
	An   = 0.0;
	Bn   = 0.0;
	Ap   = 0.0;

	lambda  = 0.0;
	shear_p = 0.0;
	bulk_p  = 0.0;
	k_p     = 0.0;
	r0n     = 0.0;
	r0p     = 0.0;

	dn  = 0.0;
	dp  = 0.0;
	//dn2  = 0.0;
	//dp2  = 0.0;
	rnn  = r0n;
	rnp1n= rnn;
	rnp  = r0p;
	rnp1p= rnp;
	//rnn2 = r0n;
	//rnp2 = r0p;

	this->zero( );
	plastic_damage2D( );
}


//full constructor
CPlaneStress2d :: 
CPlaneStress2d(int tag,double E_i,double ni_i,double f01d_i,double f02d_i,
			   double f0p_i,double beta_i,double An_i,double Bn_i,double Ap_i,
			   double gammaC_i, double dchem_i, double dchemp_i, int def_i, double dpMax_i, double dnMax_i, bool srfCompr_i,bool isDchemVar_i,double eps_u_i) :
Concrete(tag, ND_TAG_CPlaneStress2d,E_i,ni_i,f01d_i,f02d_i,f0p_i,beta_i,An_i,Bn_i,Ap_i,gammaC_i,dchem_i,dchemp_i,def_i,dpMax_i,dnMax_i,srfCompr_i,isDchemVar_i,eps_u_i),
stress2d(3),
strain2d(3),
pStrain2d_n(3),
pStrain2d_n1(3),
Stiffness2d(3,3)
{
	E    = E_i;
	ni   = ni_i;
	f01d = f01d_i;
	f02d = f02d_i;
	f0p  = f0p_i;
	beta = beta_i;
	An   = An_i;
	Bn   = Bn_i;
	Ap   = Ap_i;
	
	lambda  = E*ni/((1+ni)*(1-2*ni));
	shear_p = E/(2*(1+ni));
	bulk_p  = lambda+2/3.0*shear_p;
	k_p     = sqrt(2.0)*(f02d-f01d)/(2*f02d-f01d);
	// 10/01/2013 Diego Talledo: Added different equivalent tension definitions
	if (usedEquivalentTensionDefinition == ORIGINAL)
	{
		r0n		= sqrt(sqrt(3.0)*(k_p-sqrt(2.0))*f01d/3);
		r0p     = sqrt(f0p/sqrt(E));
	}
	else if (usedEquivalentTensionDefinition == HOMOGENEOUS)
	{
		r0n     = sqrt(sqrt(3.0)*(k_p-sqrt(2.0))*f01d/3);
		r0p     = sqrt(f0p);
	}
	else if (usedEquivalentTensionDefinition == COMPDYN)
	{
		r0n     = sqrt(3.0)*(k_p-sqrt(2.0))*f01d/3;
		r0p     = f0p;
	}

	dn    = 0.0;
	dp    = 0.0;
	rnn   = r0n;
	rnp1n = rnn;
	rnp   = r0p;
	rnp1p = rnp;

	this->zero( );
	plastic_damage2D( );
	this->doInitialStiffness();
}


//elastic constructor
CPlaneStress2d :: 
CPlaneStress2d(int tag,double E_i,double ni_i ) :
Concrete(tag, ND_TAG_CPlaneStress2d,E_i,ni_i),
stress2d(3),
strain2d(3),
pStrain2d_n(3),
pStrain2d_n1(3),
Stiffness2d(3,3)
{
	E    = E_i;
	ni   = ni_i;
	f01d = -1.0e16*E;
	f02d = -1.0e16*E;
	f0p  = 1.0e16*E;
	beta = 0.0;
	An   = 1.0;
	Bn   = 1.0;
	Ap   = 1.0;

	lambda  = E*ni/((1+ni)*(1-2*ni));
	shear_p = E/(2*(1+ni));
	bulk_p  = lambda+2*shear_p/3;
	k_p     = 0.0;
	// 10/01/2013 Diego Talledo: Added different equivalent tension definitions
	if (usedEquivalentTensionDefinition == ORIGINAL)
	{
		r0n		= sqrt(sqrt(3.0)*(k_p-sqrt(2.0))*f01d/3);
		r0p     = sqrt(f0p/sqrt(E));
	}
	else if (usedEquivalentTensionDefinition == HOMOGENEOUS)
	{
		r0n     = sqrt(sqrt(3.0)*(k_p-sqrt(2.0))*f01d/3);
		r0p     = sqrt(f0p);
	}
	else if (usedEquivalentTensionDefinition == COMPDYN)
	{
		r0n     = sqrt(3.0)*(k_p-sqrt(2.0))*f01d/3;
		r0p     = f0p;
	}

	dn  = 0.0;
	dp  = 0.0;
	rnn = r0n;
	rnp1n=rnn;
	rnp = r0p;
	rnp1p=rnp;

	//dn2  = 0.0;
	//dp2  = 0.0;
	//rnn2 = r0n;
	//rnp2 = r0p;
}



//destructor
CPlaneStress2d :: ~CPlaneStress2d( ) 
{ 

} 


//make a clone of this material
NDMaterial* CPlaneStress2d :: getCopy( ) 
{ 
	CPlaneStress2d  *clone;
	clone = new CPlaneStress2d() ;   //new instance of this class
	*clone = *this ;          //asignment to make copy
	return clone ;
}



//send back type of material
const char* CPlaneStress2d :: getType( ) const 
{
  return "PlaneStress2d" ;
}


//send back order of strain in vector form
int CPlaneStress2d :: getOrder( ) const 
{ 
  return 3 ; 
} 


int CPlaneStress2d :: setTrialStrain( const Vector &strain_from_element) 
{
	strain2d.Zero();
	strain2d = strain_from_element;
	this -> plastic_damage2D();
	return 0 ;
}


//unused trial strain functions
int CPlaneStress2d :: setTrialStrain( const Vector &v, const Vector &r )
{ 
   return this->setTrialStrain( v ) ;
} 

int CPlaneStress2d :: setTrialStrainIncr( const Vector &v ) 
{
  static Vector newStrain(3);
  newStrain(0) = strain(0,0) + v(0);
  newStrain(1) = strain(1,1) + v(1);
  newStrain(2) = strain(0,1) + v(2);
  
  return this->setTrialStrain(newStrain);  
}

int CPlaneStress2d :: setTrialStrainIncr( const Vector &v, const Vector &r ) 
{
    return this->setTrialStrainIncr(v);
}


//send back the strain
const Vector& CPlaneStress2d :: getStrain( ) 
{
  strain_vec(0) = strain2d(0) ;
  strain_vec(1) = strain2d(1) ;
  strain_vec(2) = strain2d(2) ;

  return strain_vec ;
} 


//send back the stress 
const Vector& CPlaneStress2d :: getStress( ) 
{
  stress_vec(0) = stress2d(0) ;
  stress_vec(1) = stress2d(1) ;
  stress_vec(2) = stress2d(2) ;
  return stress_vec ;
}

//send back the tangent 
const Matrix& CPlaneStress2d :: getTangent( ) 
{
	tangent_matrix = Stiffness2d;
	return tangent_matrix ;
} 


//send back the tangent 
const Matrix& CPlaneStress2d :: getInitialTangent( ) 
{
	this->doInitialStiffness();

	tangent_matrix = InitialStiffness2d;

	return tangent_matrix ;
} 

int 
CPlaneStress2d :: commitState( ) 
{
	rnp = rnp1p;
	rnn = rnp1n;
	pStrain2d_n = pStrain2d_n1;
	
	return 0;
}

int 
CPlaneStress2d :: revertToLastCommit( )
{
	return 0;
}


int 
CPlaneStress2d :: revertToStart( ) 
{
	this->zero( ) ;
	return 0;
}

int
CPlaneStress2d :: sendSelf (int commitTag, Channel &theChannel)
{
  // we place all the data needed to define material and it's state
  // int a vector object
  static Vector data(15);
  int cnt = 0;
  data(cnt++) = this->getTag();
  data(cnt++) = E    ;
  data(cnt++) = ni   ;
  data(cnt++) = f01d ;
  data(cnt++) = f02d ;
  data(cnt++) = f0p  ;
  data(cnt++) = beta ;
  data(cnt++) = An   ;
  data(cnt++) = Bn   ;
  data(cnt++) = Ap   ;
  data(cnt++) = rnn  ;
  data(cnt++) = rnp  ;
  for (int i=0; i<3; i++)
      data(cnt++) = pStrain2d_n(i);

  // send the vector object to the channel
  if (theChannel.sendVector(this->getDbTag(), commitTag, data) < 0) {
    opserr << "CPlaneStress2d::sendSelf - failed to send vector to channel\n";
    return -1;
  }
  return 0;
}

int
CPlaneStress2d :: recvSelf (int commitTag, Channel &theChannel,FEM_ObjectBroker &theBroker)
{

  // recv the vector object from the channel which defines material param and state
  static Vector data(15);
  if (theChannel.recvVector(this->getDbTag(), commitTag, data) < 0) {
    opserr << "CPlaneStress2d::recvSelf - failed to recv vector from channel\n";
    return -1;
  }

  // set the material parameters and state variables
  int cnt = 0;
  this->setTag(data(cnt++));
  E    = data(cnt++);
  ni   = data(cnt++);
  f01d = data(cnt++);
  f02d = data(cnt++);
  f0p  = data(cnt++);
  beta = data(cnt++);
  An   = data(cnt++);
  Bn   = data(cnt++);
  Ap   = data(cnt++);
  rnn  = data(cnt++);
  rnp  = data(cnt++);
  for (int i=0; i<3; i++)
      pStrain2d_n(i) = data(cnt++);
  rnp1n=rnn;
  rnp1p=rnp;

  return 0;
}

void CPlaneStress2d :: plastic_damage2D( )
{
  const double tolerance = (1.0e-14)*f0p ;
  int i,j,k;
  Vector senp1t(3);
  Matrix Stress2d(2,2);
  double V[2][2];
  double dg[2],dgn[2],dgp[2];
  double sigoct,tauoct,taun,taup;
  double alfa,norms,pint,lambdap;
  double st[3],dea[3],ls[3],ses[3],dep[3];
  double supp1[2][2],supp2[2][2],trans[3];
  double g, rhoQ, rhoP, alfasp, alfasn, rhoL, thetaL;

  for (i=0;i<3;i++)
	  dep[i]=0.0;

  this->doInitialStiffness();
  senp1t = strain2d - pStrain2d_n;
  senp1t = InitialStiffness2d*senp1t;
    
  if (beta>0.0)
  {
	  solve_eig(senp1t(0),senp1t(2),senp1t(2),senp1t(1),V,dg);
	  for (i=0;i<2;i++)
		  dgn[i]=(dg[i]-fabs(dg[i]))/2.0;
	  sigoct=(dgn[0]+dgn[1])/3;
	  tauoct=sqrt((dgn[0]-dgn[1])*(dgn[0]-dgn[1])+dgn[0]*dgn[0]+dgn[1]*dgn[1])/3;
	// 10/01/2013 Diego Talledo: Added different equivalent tension definitions
	if ((usedEquivalentTensionDefinition == ORIGINAL) || (usedEquivalentTensionDefinition == HOMOGENEOUS))
	{
		taun=sqrt(3.0)*(k_p*sigoct+tauoct);
		if (taun >= 0)
			taun = sqrt(taun);
		else
			taun = 0;
	}
	else if (usedEquivalentTensionDefinition == COMPDYN)
	{
		taun=sqrt(3.0)*(k_p*sigoct+tauoct);
	}
	  
	  for (i=0;i<2;i++)
		  dgp[i]=(dg[i]+fabs(dg[i]))/2.0;
	  taup = ((dgp[0]*dgp[0]+dgp[1]*dgp[1]))/E-2*dgp[0]*dgp[1]*ni/E;
	// 10/01/2013 Diego Talledo: Added different equivalent tension definitions
	if (usedEquivalentTensionDefinition == ORIGINAL)
	{
		taup=sqrt(sqrt(taup));
	}
	else if  (usedEquivalentTensionDefinition == HOMOGENEOUS)
	{
		taup = sqrt(sqrt(taup*E));
	}
	else if (usedEquivalentTensionDefinition == COMPDYN)
	{
		taup=sqrt(taup*E);
	}

	g=((taup/rnp)*(taup/rnp))+((taun/rnn)*(taun/rnn))-1;
	if (g>tolerance)
	{
		rhoQ = sqrt(taup*taup+taun*taun);
		rhoP = rnp*rnn*sqrt((taun*taun+taup*taup)/((taun*rnp)*(taun*rnp)+(taup*rnn)*(taup*rnn)));
		if (rnn >= rnp)
		{
			if (rhoP<rnp)
				rhoP =rnp;
			if (rhoP>rnn)
				rhoP = rnn;
		}
		else if (rnn < rnp)
		{
			if (rhoP>rnp)
				rhoP =rnp;
			if (rhoP<rnn)
				rhoP = rnn;
		}
		alfa=rhoQ/rhoP;
		// valida solo per definizioni 2 e 4
		if ((usedEquivalentTensionDefinition == ORIGINAL) || (usedEquivalentTensionDefinition == HOMOGENEOUS))
			alfa*=alfa;

		// compute dea as usual
		for (i=0;i<3;i++)
			st[i]  = senp1t[i]*(1-1/alfa);
		dea[0] = (1/E)*st[0]-(ni/E)*st[1];
		dea[1] = (1/E)*st[1]-(ni/E)*st[0];
		dea[2] = (1/shear_p)*st[2];
		norms=senp1t[0]*senp1t[0]+senp1t[1]*senp1t[1]+2*senp1t[2]*senp1t[2];
		norms=sqrt(norms);
		for (i=0;i<3;i++)
			ls[i]=senp1t[i]/norms;
		pint = 0.0;
		for (i=0;i<3;i++)
			pint+=ls[i]*dea[i];
		if (pint>0)
		{
			lambdap=1-beta*E*pint/norms; 
			for (i=0;i<3;i++)
				ses[i]=senp1t[i]*lambdap;
			// 10/01/2013 Diego Talledo: Added different equivalent tension definitions
			if ((usedEquivalentTensionDefinition == ORIGINAL) || (usedEquivalentTensionDefinition == HOMOGENEOUS))
			{
				taun*=sqrt(lambdap);
				taup*=sqrt(lambdap);
			}
			else if (usedEquivalentTensionDefinition == COMPDYN)
			{
				taun*=lambdap;
				taup*=lambdap;
			}
			g=(taup/rnp)*(taup/rnp)+(taun/rnn)*(taun/rnn)-1;
			if (g>tolerance)
			{
				for (i=0;i<3;i++)
					senp1t[i]=ses[i];
				trans[0] = (1/E)*ls[0]-(ni/E)*ls[1];
				trans[1] = (1/E)*ls[1]-(ni/E)*ls[0];
				trans[2] = (1/shear_p)*ls[2];
				for (i=0;i<3;i++)
				{
					dep[i]=beta*E*pint*trans[i];  
					pStrain2d_n1(i)=pStrain2d_n(i)+dep[i]; 
				};
			}	
			else
				pStrain2d_n1 = pStrain2d_n; 
		} 
		else 
			pStrain2d_n1 = pStrain2d_n;
	  } 
	  else 
		  pStrain2d_n1 = pStrain2d_n;
  };
  solve_eig(senp1t(0),senp1t(2),senp1t(2),senp1t(1),V,dg);
  for (i=0;i<2;i++)
	  dgn[i]=(dg[i]-fabs(dg[i]))/2;
  sigoct=(dgn[0]+dgn[1])/3; 
  tauoct=sqrt((dgn[0]-dgn[1])*(dgn[0]-dgn[1])+dgn[0]*dgn[0]+dgn[1]*dgn[1])/3; 
	// 10/01/2013 Diego Talledo: Added different equivalent tension definitions
	if ((usedEquivalentTensionDefinition == ORIGINAL) || (usedEquivalentTensionDefinition == HOMOGENEOUS))
	{
		taun=sqrt(3.0)*(k_p*sigoct+tauoct);
		if (taun >= 0)
			taun = sqrt(taun);
		else
			taun = 0;
	}
	else if (usedEquivalentTensionDefinition == COMPDYN)
	{
		taun=sqrt(3.0)*(k_p*sigoct+tauoct);
	}
	for (i=0;i<2;i++)
		dgp[i]=(dg[i]+fabs(dg[i]))/2;
	taup = ((dgp[0]*dgp[0]+dgp[1]*dgp[1]))/E-2*dgp[0]*dgp[1]*ni/E;
	// 10/01/2013 Diego Talledo: Added different equivalent tension definitions
	if (usedEquivalentTensionDefinition == ORIGINAL)
	{
		taup=sqrt(sqrt(taup));
	}
	else if  (usedEquivalentTensionDefinition == HOMOGENEOUS)
	{
		taup = sqrt(sqrt(taup*E));
	}
	else if (usedEquivalentTensionDefinition == COMPDYN)
	{
		taup=sqrt(taup*E);
	}
	g=(taup/rnp)*(taup/rnp)+(taun/rnn)*(taun/rnn)-1;
	if (g>tolerance)
	{
		rhoQ = sqrt(taup*taup+taun*taun); 
		rhoP = rnp*rnn*sqrt((taun*taun+taup*taup)/((taun*rnp)*(taun*rnp)+(taup*rnn)*(taup*rnn))); 
		if (rnn >= rnp) 
		{ 
			if (rhoP<rnp)
				rhoP =rnp;
			if (rhoP>rnn)
				rhoP = rnn;
		}
		else if (rnn < rnp)
		{
			if (rhoP>rnp)
				rhoP =rnp;
			if (rhoP<rnn)
				rhoP = rnn;
		}
		alfa=rhoQ/rhoP; 
		thetaL = atan((rnp*rnp)/(rnn*rnn));
		rhoL=sqrt((rnp*rnp*rnn*rnn)/(rnn*rnn*sin(thetaL)*sin(thetaL)+rnp*rnp*cos(thetaL)*cos(thetaL))); 
		if (((rhoP>rhoL) && (rhoP<=rnn)) || ((rhoP>=rnn) && (rhoP<rhoL))) {
			alfasp=1+(alfa-1)*(rnn-rhoP)/(rnn-rhoL);
			rnp1p=rnp*alfasp;
			rnp1n=sqrt((rnp1p*rnp1p*taun*taun)/(rnp1p*rnp1p-taup*taup)); 
		} else 
		{
			alfasn=1+(alfa-1)*(rhoP-rnp)/(rhoL-rnp);
			rnp1n=rnn*alfasn; 
			rnp1p=sqrt((rnp1n*rnp1n*taup*taup)/(rnp1n*rnp1n-taun*taun));
		}
		//if (rnp1n<tolerance)
		//	dn = 0.0;
		//else
		//	dn=1-r0n/rnp1n*(1-An)-An*exp(Bn*(1-rnp1n/r0n));
		//// 11/01/2013 Diego Talledo: Added Negative Damage Limit
		///*if (dn>=1.0)
		//	dn=1.0;*/
		//if (dn>=dnMax)
		//	dn = dnMax;
		//if (dn<0.0)
		//	dn=0.0;
		//if (rnp1p<tolerance)
		//	dp = 0.0;
		//else
		//	dp=1-((r0p*r0p)/(rnp1p*rnp1p))*exp(Ap*(1-(rnp1p*rnp1p)/(r0p*r0p)));
		////Original leo 26/11/2012
		///*if (dp >=1.0)
		//	dp=1.0;*/
		//// 11/01/2013 Diego Talledo: Added Negative Damage Limit
		//if (dp>=dpMax)
		//	dp = dpMax;
		//if (dp <0)
		//	dp=0.0;
	} else
	{
		rnp1n=rnn;
		rnp1p=rnp;
		//if (rnp1n<tolerance)
		//	dn = 0.0;
		//else
		//	dn=1-r0n/rnp1n*(1-An)-An*exp(Bn*(1-rnp1n/r0n)); 
		//// 11/01/2013 Diego Talledo: Added Negative Damage Limit
		///*if (dn>=1.0)
		//	dn=1.0;*/
		//if (dn>=dnMax)
		//	dn = dnMax;
		//if (dn<0.0)
		//	dn=0.0;
		//if (rnp1p<tolerance)
		//	dp = 0.0;
		//else
		//	dp=1-((r0p*r0p)/(rnp1p*rnp1p))*exp(Ap*(1-(rnp1p*rnp1p)/(r0p*r0p))); 
		////Original leo 26/11/2012
		///*if (dp >=1.0)
		//	dp=1.0;*/
		//// 11/01/2013 Diego Talledo: Added Negative Damage Limit
		//if (dp>=dpMax)
		//	dp = dpMax;
		//if (dp <0)
		//	dp=0.0;
	};

	if (rnp1n<tolerance)
		dn = 0.0;
	else
	{
		// 10/01/2013 Diego Talledo: Added different equivalent tension definitions
		if ((usedEquivalentTensionDefinition == ORIGINAL) || (usedEquivalentTensionDefinition == HOMOGENEOUS))
		{
			dn=1-r0n/rnp1n*(1-An)-An*exp(Bn*(1-rnp1n/r0n)); 
		}
		else if (usedEquivalentTensionDefinition == COMPDYN)
		{
			dn=1-(sqrt(r0n))/(sqrt(rnp1n))*(1-An)-An*exp(Bn*(1-(sqrt(rnp1n))/(sqrt(r0n)))); 
		}
	}
	// 11/01/2013 Diego Talledo: Added Negative Damage Limit
	/*if (dn>=1.0)
		dn=1.0;*/
	if (dn>=dnMax)
		dn = dnMax;
	if (dn<0.0)
		dn=0.0;
	if (rnp1p<tolerance)
		dp = 0.0;
	else
	{
		dp=1-((r0p*r0p)/(rnp1p*rnp1p))*exp(Ap*(1-(rnp1p*rnp1p)/(r0p*r0p))); 
		// 10/01/2013 Diego Talledo: Added different equivalent tension definitions
		if ((usedEquivalentTensionDefinition == ORIGINAL) || (usedEquivalentTensionDefinition == HOMOGENEOUS))
		{
			dp=1-((r0p*r0p)/(rnp1p*rnp1p))*exp(Ap*(1-(rnp1p*rnp1p)/(r0p*r0p)));  
		}
		else if (usedEquivalentTensionDefinition == COMPDYN)
		{
			dp=1-((r0p)/(rnp1p))*exp(Ap*(1-(rnp1p)/(r0p))); 
		}
	}
	//Original leo 26/11/2012
	/*if (dp >=1.0)
		dp=1.0;*/
	// 11/01/2013 Diego Talledo: Added Negative Damage Limit
	if (dp>=dpMax)
		dp = dpMax;
	if (dp <0)
		dp=0.0;

	// 13/01/2013 Diego Talledo: Added Shear Retention Factor
	if (gammaC <= 0.0) {
		SRF12 = 0.0;
	}
	else {
		SRF12 = 1-abs(strain2d(2))/(2*gammaC);  //gammaC = epsilonREF_1,2
		if (SRF12 <= 0.0) {
			SRF12 = 0.0;
		}
	}

  for (i=0;i<2;i++)
	  for (j=0;j<2;j++)
	  {
		  supp1[i][j]=0.0;
		  supp2[i][j]=0.0;
		  supp1[i][j]=supp1[i][j]+V[i][j]*dgn[j];
		  supp2[i][j]=supp2[i][j]+V[i][j]*dgp[j];
	  };
  for (i=0;i<2;i++)
	  for (j=0;j<2;j++)
	  {
		  Dn[i][j]=0.0;
		  Dp[i][j]=0.0;
		  for (k=0;k<2;k++)
		  {
			  Dn[i][j]=Dn[i][j]+supp1[i][k]*V[j][k];
			  Dp[i][j]=Dp[i][j]+supp2[i][k]*V[j][k];
		  };
	  };

	// 11/03/2013 Diego Talledo: Added Environmental Chemical Damage
	// Compute dnstar and dpstar
	dnstar = dn + dchemn - (dchemn * dn);
	//dpstar = dp + dchem - (dchem * dp);
	dpstar = dp + dchemp - (dchemp * dp);

	// 14/07/2014 Diego Talledo: Aggiunto dchemn variabile
	// dchem variabile: retta tra dchemn e 1. Con 1 raggiunto per epsilon pari a eps_u
	if (isDchemVariableSelected) {
		double dchem_var = (1-dchemn)/tau_n_eps_u*rnp1n+dchemn;
		dnstar = dn + dchem_var - (dchem_var * dn);
	}
	// Fine prova dchem variabile

  //compute Cauchy stress
  Stress2d.Zero();
  for (i=0;i<2;i++)
	  for (j=0;j<2;j++)
	  {

		  //Stress2d(i,j)=(1-dn)*Dn[i][j]+(1-dp)*Dp[i][j];
		  // 11/03/2013 Diego Talledo: Added Environmental Chemical Damage
		  Stress2d(i,j)=(1-dnstar)*Dn[i][j]+(1-dpstar)*Dp[i][j];
		  // 13/01/2013 Diego Talledo: Added Shear Retention Factor
		  if (((i==0) && (j==1)) || ((i==1) && (j==0)))
			  if (!srfCompr) // 11/03/2013 Diego Talledo: Apply SRF also to compression.
				  Stress2d(i,j)=(1-dnstar)*Dn[i][j]+(1-(1-SRF12)*dpstar)*Dp[i][j];
			  else
				  Stress2d(i,j)=(1-(1-SRF12)*dnstar)*Dn[i][j]+(1-(1-SRF12)*dpstar)*Dp[i][j];
	  }
  stress2d(0) = Stress2d(0,0);
  stress2d(1) = Stress2d(1,1);
  stress2d(2) = Stress2d(0,1);
  //compute stiffness matrix
  if (dn > 0.0 && dn < 0.99999)
	  Stiffness2d = (1-dn) * InitialStiffness2d;
  else
	  Stiffness2d = InitialStiffness2d;
  return ;
} 

void CPlaneStress2d :: doInitialStiffness( )
{
	InitialStiffness2d.Zero();
	InitialStiffness2d(0,0) = E/(1-(ni*ni));
	InitialStiffness2d(1,1) = InitialStiffness2d(0,0);
	InitialStiffness2d(2,2) = shear_p;
	InitialStiffness2d(0,1) = E*ni/(1-(ni*ni));
	InitialStiffness2d(1,0) = InitialStiffness2d(0,1);
	return ;
}

void CPlaneStress2d :: solve_eig(double A, double B, double C, double D, double V[2][2], double d[2] )
{
	double tolerance = 0.1e-20;
	double lambda1  = 0.0;
	double lambda2 = 0.0;
	double v1x = 0.0;
	double v1y = 0.0;
	double v2x = 0.0;
	double v2y = 0.0;

	if(B*C <= tolerance  )
	{
		lambda1 = A; v1x = 1; v1y = 0;
		lambda2 = D; v2x = 0; v2y = 1;
		d[0] = lambda1;
		d[1] = lambda2;
		V[0][0] = v1x;
		V[1][0] = v1y;
		V[0][1] = v2x;
		V[1][1] = v2y;
		return;
	}
	double tr = A + D;
	double det = A * D - B * C;
	double S = sqrt( (tr/2)*(tr/2) - det );
	lambda1 = tr/2 + S;
	lambda2 = tr/2 - S;
	double S2 = ((A-D)/2)*((A-D)/2) + B * C;
	double SS = 0.0;
	if (S2 > 0.0)
		SS = sqrt(S2);
	if( A - D < 0 ) {
		v1x = C;
		v1y = - (A-D)/2 + SS;
		v2x = + (A-D)/2 - SS;
		v2y = B;
	} else {
		v2x = C;
		v2y = - (A-D)/2 - SS;
		v1x = + (A-D)/2 + SS;
		v1y = B;
	}
	double n1 = sqrt(v1x*v1x+v1y*v1y);  
	v1x /= n1;
	v1y /= n1;
	double n2 = sqrt(v2x*v2x+v2y*v2y);
	v2x /= n2;
	v2y /= n2;

	d[0] = lambda1;
	d[1] = lambda2;

	V[0][0] = v1x;
	V[1][0] = v1y;
	V[0][1] = v2x;
	V[1][1] = v2y;
}

// Diego Gennaio 2015 : Per implementazione Danno Globale
const Vector & CPlaneStress2d::getElasticFreeEnergy()
{
	static Vector tmp(2); // Psi0- Psi0+

	// Variabili temporanee - Sigma Elastica + Decomposizione Spettrale
	Vector sigmaElastic(3);
	double V[2][2];
	double dg[2],dgn[2],dgp[2];
	double supp1[2][2],supp2[2][2],trans[3];

	this->doInitialStiffness();
	sigmaElastic = strain2d - pStrain2d_n;
	sigmaElastic = InitialStiffness2d*sigmaElastic;
    
    solve_eig(sigmaElastic(0),sigmaElastic(2),sigmaElastic(2),sigmaElastic(1),V,dg);
	for (int i=0; i<2; i++)
		dgn[i]=(dg[i]-fabs(dg[i]))/2.0;
	// Calcolare DGP XXX
	for (int i=0;i<2;i++)
		dgp[i]=(dg[i]+fabs(dg[i]))/2.0;

	// Bring it to the reference system
	for (int i=0; i<2; i++)
	  for (int j=0; j<2; j++)
	  {
		  supp1[i][j]=0.0;
		  supp2[i][j]=0.0;
		  supp1[i][j]=supp1[i][j]+V[i][j]*dgn[j];
		  supp2[i][j]=supp2[i][j]+V[i][j]*dgp[j];
	  };
  for (int i=0; i<2; i++)
	  for (int j=0; j<2; j++)
	  {
		  Dn[i][j]=0.0;
		  Dp[i][j]=0.0;
		  for (int k=0; k<2; k++)
		  {
			  Dn[i][j]=Dn[i][j]+supp1[i][k]*V[j][k];
			  Dp[i][j]=Dp[i][j]+supp2[i][k]*V[j][k];
		  };
	  };

	// Psi0-
	tmp(0) = 0.5*(Dn[0][0]*(strain2d(0)-pStrain2d_n(0)) + Dn[1][1]*(strain2d(1)-pStrain2d_n(1)) + Dn[0][1]*(strain2d(2)-pStrain2d_n(2)));
	// Psi9+
	tmp(1) = 0.5*(Dp[0][0]*(strain2d(0)-pStrain2d_n(0)) + Dp[1][1]*(strain2d(1)-pStrain2d_n(1)) + Dp[0][1]*(strain2d(2)-pStrain2d_n(2)));

	return tmp;
}

const Vector & CPlaneStress2d::getDamagedFreeEnergy()
{
	static Vector tmp(2);
	const Vector &ElFreeEnergy = this->getElasticFreeEnergy();

	tmp(0) = ElFreeEnergy(0) * (1-this->dnstar);
	tmp(1) = ElFreeEnergy(1) * (1-this->dpstar);

	return tmp;
}



