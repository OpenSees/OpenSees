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
                                                                        
// Written: Diego Turello(*), Alborz Ghofrani and Pedro Arduino
//			Sep 2017, University of Washington
//          (*) Universidad Nacional de Córdoba, FCEFyN. Depto Estructuras.
//              Universidad Tecnológica Nacional, GIMNI.
//              CONICET
// 
// Description: This file contains the implementation for the Borja material class.
// MultiaxialCyclicPlasticity for clays 
// Borja R.I, Amies, A.P.Multiaxial Cyclic Plasticity Model for Clays,
// ASCE J.Geotech.Eng.Vol 120, No 6, 1051 - 1070
//            

#include <J2CyclicBoundingSurface.h>


#include <Information.h>
#include <Parameter.h>

#include <string.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <elementAPI.h>

//parameters

char  unsigned      J2CyclicBoundingSurface::m_ElastFlag = 1;


void* OPS_J2CyclicBoundingSurface()
{
    int numdata = OPS_GetNumRemainingInputArgs();

    if (numdata < 9) {
	opserr << "WARNING: Insufficient arguements\n";
	opserr << "Want: nDMaterial J2CyclicBoundingSurface tag? G? K? su? rho? h? m? k_in? beta?\n";
	return 0;
    }

    int tag; 
	
    numdata = 1;
    if (OPS_GetIntInput(&numdata,&tag) < 0) {
	opserr << "WARNING invalid J2CyclicBoundingSurface tag\n";
	return 0;
    }

    double data[8] = {0,0,0,0,0,0,0,0};
    numdata = OPS_GetNumRemainingInputArgs();
    if (numdata != 8) {
		opserr << "WARNING error in  J2CyclicBoundingSurface number of arg incorrect\n";
		return 0;
    }
    if (OPS_GetDoubleInput(&numdata,data)) {
	opserr << "WARNING invalid J2CyclicBoundingSurface double inputs\n";
	return 0;
    }

    NDMaterial* mat = new J2CyclicBoundingSurface(tag,data[0],data[1],data[2],data[3],data[4],data[5],data[6],data[7]);
    if (mat == 0) {
	opserr << "WARNING: failed to create J2CyclicBoundingSurface material\n";
	return 0;
    }

    return mat;
}


//null constructor
J2CyclicBoundingSurface::J2CyclicBoundingSurface() :
	NDMaterial()
    {

	}


//full constructor
J2CyclicBoundingSurface::J2CyclicBoundingSurface(int  tag,
		double G,
		double K,
		double su,
		double rho,
		double h,
		double m,
		double k_in,
        double beta)
:
NDMaterial(tag, ND_TAG_J2CyclicBoundingSurface),
m_sigma0_n(6), m_sigma0_np1(6),  m_stress_n(6), m_stress_np1(6), m_strain_np1(6), m_strain_n(6), m_Ktan_np1(6,6), m_Kelas(6,6)
{
  m_young = 9.*K*G/(3.*K+G);
  m_poiss = (3.*K-2.*G)/2./ (3. * K + G);
  //opserr << m_poiss << endln;
  if (m_poiss > 0.5) {
	  opserr << "\n ERROR! J2CyclicBoundingSurface Poiss can not be grater than 0.50!" << endln;
	  return;
  }
  m_su    = su;
  m_bulk  = K ; // E / 3 / (1-2*nu)
  m_shear = G;  // E / 2 / (1 + nu);
  m_R = sqrt(2.)*su; // sqrt(2.)*su; with su from simple shear test //sqrt(8. / 3.)*su from undrained unconfined triaxial test
  m_density     = rho;
  m_kappa_inf   = k_in;
  m_h_par = h;
  m_m_par = m;
  m_beta  = beta;
  m_kappa_n = m_kappa_inf;
  m_xi_n = 2.*m_shear;

  m_isElast2Plast = false;

  doInitialTangent();

 /* plastic_integrator();*/
}


//destructor
J2CyclicBoundingSurface :: ~J2CyclicBoundingSurface( ) 
{  } 


//--------------------Plasticity-------------------------------------

void J2CyclicBoundingSurface::zero()
{
}

void J2CyclicBoundingSurface::integrate()
{
	// Force elastic response
	if (m_ElastFlag == 0) {
		elastic_integrator();
		//opserr << "elastic" << endln;
	}
	// ElastoPlastic response
	else {
		J2CyclicBoundingSurface::plastic_integrator();
		//opserr << "plastic" << endln;
	}
}

void J2CyclicBoundingSurface::elastic_integrator()
{
    Vector delta_strain = m_strain_np1 - m_strain_n; //delta strain for the step
	m_stress_np1 = m_stress_n + m_Kelas * delta_strain;
	
}


//plastic integration routine
void J2CyclicBoundingSurface :: plastic_integrator( )
{
  const double tol_rel = (1.0e-6) ;
  Matrix ACO(6, 6);
  ACO(0, 0) = 1.0; ACO(1, 1) = 1.0; ACO(2, 2) = 1.0; ACO(3, 3) = 0.5; ACO(4, 4) = 0.5; ACO(5, 5) = 0.5;
  Vector AVI(6);
  AVI(0) = 1.0; AVI(1) = 1.0; AVI(2) = 1.0;

  Vector delta_strain  = m_strain_np1 - m_strain_n; //delta strain for the step
   
  double vol_delta_strain  = trace(delta_strain) ;        // volumetric strain
  Vector dev_delta_strain  = dev_strain_op(delta_strain); // deviatoric strain

  Vector dev_stress_n(6)  ; //deviatoric stress step n
  Vector dev_stress_np1(6); //deviatoric stress step n+1
  Vector dev_sigma0_np1(6); //deviatoric stress 0 step n+1

  dev_stress_n   = dev_stress_op(m_stress_n);
  // initialize m_sigma0_np1
  m_sigma0_np1   = m_sigma0_n;
  dev_sigma0_np1 = dev_stress_op(m_sigma0_np1);

  double H_n ;
  double H_np1;
   
  double cond = -1;
  double norm_dev_stress_n   = sqrt(inner_product_stress(dev_stress_n, dev_stress_n));
  double norm_dev_sigma0_np1 = sqrt(inner_product_strain(dev_sigma0_np1, dev_sigma0_np1));
  Vector auxVec1;
  //****************************************
  // Condition over kappa
  if (m_isElast2Plast) {
	  //m_kappa_n = m_kappa_inf;
	  m_kappa_n = abs(m_R / norm_dev_stress_n - 1.0);
	  H_n       = H(m_kappa_n);
	  m_xi_n    = 2.0*m_shear * H_n / (H_n + 3.0*m_shear);
	  m_kappa_np1 = m_kappa_n;
	  m_xi_np1    = m_xi_n;
	  m_isElast2Plast = false;
  }
  else {
	  // Esta condicion no estoy seguro este bien
	  //if (m_kappa_n == 0.0) {
		 // m_kappa_n = m_kappa_inf;
	  //}
	  //if (m_xi_n == 0.0) {
		 // m_xi_n = 2. * m_shear;
	  //}
	  m_kappa_np1 = m_kappa_n;
	  m_xi_np1    = m_xi_n;
	  //m_xi_np1 = 2. * m_shear;
	  //****************************************
	  // check loading/unloading
	  if (norm_dev_stress_n > 0 || norm_dev_sigma0_np1 > 0)
	  {
		  auxVec1 = (-(1 + m_kappa_n)*dev_stress_n - m_kappa_n*(1 + m_kappa_n)*(dev_stress_n - dev_sigma0_np1)) /
			  (
				  inner_product_stress(dev_stress_n, (dev_stress_n - dev_sigma0_np1)) + m_kappa_n*inner_product_stress((dev_stress_n - dev_sigma0_np1), (dev_stress_n - dev_sigma0_np1))
				  );
		  cond = inner_product_contra_covariant(auxVec1, dev_delta_strain);
		  if (cond > 0.)
		  {
			  m_sigma0_np1 = m_stress_n;
			  m_kappa_n    = m_kappa_inf;
			  m_kappa_np1  = m_kappa_n;
			  H_n          = H(m_kappa_n);
			  m_xi_np1     = 2. * m_shear;
			  dev_sigma0_np1 = dev_stress_op(m_sigma0_np1);
		  }
	  }
  }
  // 
  // Newton for the material, solve for n_xi_np1 and n_kappa_np1
  // 
  //  Initial Residual 
  H_n   = H(m_kappa_n);
  H_np1 = H(m_kappa_np1);

  Vector auxVec2 = dev_stress_n + m_xi_np1 * ACO * dev_delta_strain + m_kappa_np1 * (dev_stress_n + m_xi_np1 * ACO * dev_delta_strain - dev_sigma0_np1);
  double auxVec2_norm = sqrt(inner_product_stress(auxVec2, auxVec2));

  Vector res(2); double res_norm;
  //res(0) = m_xi_np1 + 3. * m_shear * m_xi_np1 * ((1. - m_beta) / H_n + m_beta / H_np1) - 2. * m_shear;
  //res(1) = auxVec2_norm - m_R;

  res(0) = m_xi_np1/2./ m_shear + 3./2. * m_xi_np1 * ((1. - m_beta) / H_n + m_beta / H_np1) - 1.0;
  res(1) = auxVec2_norm/ m_R - 1.0;

  res_norm = vector_norm(res);
  
  // Initialize variables for the Newton
  int iteration_counter = 0;
  const int max_iterations = 100;
  Matrix Ktan(2,2);

  Vector incVar(2);
  double tol_material = tol_rel*res_norm;
  double constNew = 1.0;

  // Main while in the Newton loop
  while (res_norm>(tol_material+1.e-10) &&  iteration_counter<max_iterations && auxVec2_norm>1e-10)
  {
	  iteration_counter ++ ;
	  //Ktan(0, 0) = 1. + 3. * m_shear * ((1. - m_beta) / H_n + m_beta / H_np1);
	  //Ktan(0, 1) = 3. * m_shear*m_xi_np1*(-m_beta / pow(H_np1, 2) * m_m_par*m_h_par*pow(m_kappa_np1, (m_m_par - 1.)));
	  //Ktan(1, 0) = inner_product(auxVec2, (dev_delta_strain + m_kappa_np1*dev_delta_strain)) / auxVec2_norm; //vector_norm(auxVec2)
	  //Ktan(1, 1) = inner_product(auxVec2, (dev_stress_n   + m_xi_np1   *dev_delta_strain-dev_sigma0_np1)) / auxVec2_norm; //vector_norm(auxVec2)

	  Ktan(0, 0) = 1./2./ m_shear + 3./2. * ((1. - m_beta) / H_n + m_beta / H_np1);
	  Ktan(0, 1) = 3./2. * m_xi_np1*(-m_beta / pow(H_np1, 2) * m_m_par*m_h_par*pow(m_kappa_np1, (m_m_par - 1.)));
	  Ktan(1, 0) = inner_product_contra_covariant(auxVec2, (dev_delta_strain + m_kappa_np1*dev_delta_strain)) / auxVec2_norm / m_R; //vector_norm(auxVec2)
	  Ktan(1, 1) = inner_product_stress(auxVec2, (dev_stress_n + m_xi_np1 * ACO * dev_delta_strain - dev_sigma0_np1)) / auxVec2_norm / m_R; //vector_norm(auxVec2)
	  //opserr << "res = " << res << endln;
	  //opserr << "Ktan = " << Ktan << endln;
	  //
	  // Solve the system
	  Ktan.Solve(res, incVar);
	  if ((m_xi_np1<constNew*incVar(0)) || (m_kappa_np1<constNew*incVar(1)) ) {
		  constNew = constNew / 2.0;
	  }
	  // Actualize variable
	  incVar = -constNew * incVar;
	  //opserr << "incVar = " << incVar << endln;
	  m_xi_np1    = m_xi_np1    + incVar(0);
	  m_kappa_np1 = m_kappa_np1 + incVar(1);

	  // m_xi_np1 should be greater than 0, puts a value 10% of the maximun
	  // if (m_xi_np1<=0) { m_xi_np1 = 0.1 * 2. * m_shear; }	 
	  // m_xi_np1 should be less than 2*shear, puts the maximun value
	  // if (m_xi_np1>2*m_shear) { m_xi_np1 = 2. * m_shear; }
	  // kappa should be possitive
	  // if (m_kappa_np1<0){ m_kappa_np1 = 1.0e-5;}

	  // Actualize hardening
	  H_np1 = H(m_kappa_np1);

	  // Calcule sigma_np1
	  m_stress_np1 = m_stress_n + m_bulk * trace(delta_strain) * AVI + m_xi_np1 * ACO * dev_delta_strain;
	  
	  //m_stress_np1(0) = m_stress_n(0) + m_bulk*trace(delta_strain) + m_xi_np1*dev_delta_strain(0);
	  //m_stress_np1(1) = m_stress_n(1) + m_bulk*trace(delta_strain) + m_xi_np1*dev_delta_strain(1);
	  //m_stress_np1(2) = m_stress_n(2) + m_bulk*trace(delta_strain) + m_xi_np1*dev_delta_strain(2);
	  //m_stress_np1(3) = m_stress_n(3) + m_xi_np1*dev_delta_strain(3);
	  //m_stress_np1(4) = m_stress_n(4) + m_xi_np1*dev_delta_strain(4);
	  //m_stress_np1(5) = m_stress_n(5) + m_xi_np1*dev_delta_strain(5);

	  dev_stress_np1 = dev_stress_op(m_stress_np1);
	  //  Residual 
	  auxVec2 = dev_stress_n + m_xi_np1 * ACO * dev_delta_strain + m_kappa_np1 * (dev_stress_n + m_xi_np1 * ACO * dev_delta_strain - dev_sigma0_np1);
	  auxVec2_norm = sqrt(inner_product_stress(auxVec2, auxVec2));

	  //res(0) = m_xi_np1 + 3. * m_shear*m_xi_np1*((1. - m_beta) / H_n + m_beta / H_np1) - 2. * m_shear;
	  //res(1) = auxVec2_norm - m_R;

	  res(0) = m_xi_np1 / 2. / m_shear + 3. / 2. * m_xi_np1 * ((1. - m_beta) / H_n + m_beta / H_np1) - 1.0;
	  res(1) = auxVec2_norm / m_R - 1.0;

	  res_norm = vector_norm(res);
	   
  }
  if (iteration_counter == max_iterations) 
  {
	  opserr << "Material model does not converge after  " << iteration_counter << " iterations/n " << endln;
  }

  return ;
} 

// DT Functions
// Trace Operator
double J2CyclicBoundingSurface :: trace(Vector V)
{
	return V(0) + V(1) + V(2);
}
// Deviatoric strain operator
Vector J2CyclicBoundingSurface::dev_strain_op(Vector V)
{
	double vol_strain = trace(V);
	for (int i = 0; i < 3; i++)
		V(i) = V(i) - 1. / 3. * vol_strain;
	return V;
}
// Deviatoric stress operator
Vector J2CyclicBoundingSurface::dev_stress_op(Vector V)
{
	double vol_stress = 1./3.*trace(V);
	for (int i = 0; i < 3; i++)
		V(i) = V(i) - vol_stress;
	return V;
}
// Inner product stress
double J2CyclicBoundingSurface::inner_product_stress(Vector x, Vector y)
{
	double inner = 0.0;
	//int numberOfComponents = sizeof(x);
	if (x.Size() == 3) {
		for (int i = 0; i < x.Size(); i++)
			inner = inner + x(i)*y(i);
	}
	else {
		for (int i = 0; i < x.Size(); i++)
			if (i<3){inner = inner + x(i)*y(i);}
			else { inner = inner + 2*x(i)*y(i); }
	}
	return inner;
}
// Inner product strain
double J2CyclicBoundingSurface::inner_product_strain(Vector x, Vector y)
{
	double inner = 0.0;
	//int numberOfComponents = sizeof(x);
	if (x.Size() == 3) {
		for (int i = 0; i < x.Size(); i++)
			inner = inner + x(i)*y(i);
	}
	else {
		for (int i = 0; i < x.Size(); i++)
			if (i<3) { inner = inner + x(i)*y(i); }
			else { inner = inner + 0.5 * x(i)*y(i); }
	}
	return inner;
}
// Inner product contra covariant
double J2CyclicBoundingSurface::inner_product_contra_covariant(Vector x, Vector y)
{
	double inner = 0.0;
	//int numberOfComponents = sizeof(x);
	for (int i = 0; i < x.Size(); i++)
        inner = inner + x(i)*y(i);

	return inner;
}

// Norm of a vector
double J2CyclicBoundingSurface::vector_norm(Vector x)
{
	double vector_norm = sqrt(inner_product_contra_covariant(x, x));
	return vector_norm;
}

void J2CyclicBoundingSurface::index_map(int matrix_index, int & i, int & j)
{
}





//hardening function
double J2CyclicBoundingSurface :: H( double kappa ) 
{
//  H(kappa) = h*kappa^m
	return    m_h_par*pow(kappa, m_m_par);
}


void J2CyclicBoundingSurface::Print(OPS_Stream & s, int flag)
{
}

NDMaterial*
J2CyclicBoundingSurface::getCopy (void)
{
  opserr << "J2CyclicBoundingSurface::getCopy -- subclass responsibility\n"; 
  exit(-1);
  return 0;
}

const char*
J2CyclicBoundingSurface::getType (void) const
{
    //opserr << "J2CyclicBoundingSurface::getType -- subclass responsibility\n";
    //exit(-1);
    //return 0;
	return "ThreeDimensional";
}

int
J2CyclicBoundingSurface::getOrder (void) const
{
    opserr << "J2CyclicBoundingSurface::getOrder -- subclass responsibility\n";
    exit(-1);
    return 0;
}


NDMaterial * J2CyclicBoundingSurface::getCopy(const char * type)
{
	
	if (strcmp(type, "ThreeDimensional") == 0 || strcmp(type, "3D") == 0) {
		J2CyclicBoundingSurface *clone;
		clone = new J2CyclicBoundingSurface(this->getTag(), m_shear, m_bulk, m_su, m_density, m_h_par, m_m_par, m_kappa_inf, m_beta);
		return clone;
	}
	else {
		opserr << "J2CyclicBoundingSurface::getCopy failed to get copy: " << type << endln;
		return 0;
	}
}

int
J2CyclicBoundingSurface::commitState( ) 
{
	m_sigma0_n = m_sigma0_np1;
	m_stress_n = m_stress_np1;
	m_kappa_n  = m_kappa_np1;
	m_xi_n     = m_xi_np1;
	// Actualize m_strain_n
	m_strain_n = m_strain_np1; //Guarda m_strain_n para calcular el delta_strain del proximo paso

	return 0;
}

int 
J2CyclicBoundingSurface::revertToLastCommit( ) 
{
  return 0;
}


int 
J2CyclicBoundingSurface::revertToStart( ) {

  // added: C.McGann, U.Washington for InitialStateAnalysis
  if (ops_InitialStateAnalysis) {
	// do nothing, keep state variables from last step
  } else {
	// normal call for revertToStart (not initialStateAnalysis)
    this->zero( ) ;
  }

  return 0;
}

int
J2CyclicBoundingSurface::setParameter(const char **argv, int argc,
				      Parameter &param)
{
	int theMaterialTag;
	theMaterialTag = atoi(argv[1]);

	if (theMaterialTag == this->getTag()) {

		if (strcmp(argv[0], "updateMaterialStage") == 0) {     // enforce elastic/elastoplastic response
			return param.addObject(1, this);
		}
		else if (strcmp(argv[0], "materialState") == 0) {     // enforce elastic/elastoplastic response
			return param.addObject(2, this);
		}
	}
	return -1;
}

int 
J2CyclicBoundingSurface::updateParameter(int responseID, Information &info)
{ 
	// called updateMaterialStage in tcl file
	if (responseID == 1) {
		m_ElastFlag = info.theInt;
		m_isElast2Plast = true;
		return 0;
	}
	// called materialState in tcl file
	else if (responseID == 2) {
		m_ElastFlag = (int)info.theDouble;
		m_isElast2Plast = true;
		return 0;
	}
	return -1;
}

int
J2CyclicBoundingSurface::activateParameter(int paramID)
{
  parameterID = paramID;

  return 0;
}

int
J2CyclicBoundingSurface::sendSelf(int commitTag, Channel &theChannel)
{
  // we place all the data needed to define material and it's state
  // int a vector object
  //static Vector data(10+9);
  //int cnt = 0;
  //data(cnt++) = this->getTag();
  //data(cnt++) = bulk;
  //data(cnt++) = shear;
  //data(cnt++) = sigma_0;
  //data(cnt++) = sigma_infty;
  //data(cnt++) = delta;
  //data(cnt++) = Hard;
  //data(cnt++) = eta;
  //data(cnt++) = rho;

  //data(cnt++) = xi_n;

  //for (int i=0; i<3; i++) 
  //  for (int j=0; j<3; j++) 
  //    data(cnt++) = epsilon_p_n(i,j);


  //// send the vector object to the channel
  //if (theChannel.sendVector(this->getDbTag(), commitTag, data) < 0) {
  //  opserr << "J2CyclicBoundingSurface::sendSelf - failed to send vector to channel\n";
  //  return -1;
  //}

  return 0;
}

int
J2CyclicBoundingSurface::recvSelf (int commitTag, Channel &theChannel, 
			 FEM_ObjectBroker &theBroker)
{
	//// recv the vector object from the channel which defines material param and state
	//static Vector data(10+9);
	//if (theChannel.recvVector(this->getDbTag(), commitTag, data) < 0) {
	// opserr << "J2CyclicBoundingSurface::recvSelf - failed to recv vector from channel\n";
	//return -1;
	// }

	//// set the material parameters and state variables
	//int cnt = 0;
	//this->setTag(data(cnt++));
	//bulk = data(cnt++);
	// shear = data(cnt++);
	// sigma_0 = data(cnt++);
	//sigma_infty = data(cnt++);
	//delta = data(cnt++);
	// Hard = data(cnt++);
	//eta = data(cnt++);
  //rho = data(cnt++);
  //
  //xi_n = data(cnt++);
  //
  //for (int i=0; i<3; i++)
  //  for (int j=0; j<3; j++) 
  //    epsilon_p_n(i,j) = data(cnt++);
  //
  //epsilon_p_nplus1 = epsilon_p_n;
  //xi_nplus1        = xi_n;
  //
  return 0;
}

// get the strain and integrate plasticity equations
int
J2CyclicBoundingSurface::setTrialStrain(const Vector &strain_from_element)
{
	m_strain_np1 =  strain_from_element; // -1.0 is for geotechnical sign convention
	//opserr << "strain = " << strain_from_element;
	//this->plastic_integrator();
	this->integrate();
	
	return 0;
}

// unused trial strain functions
int
J2CyclicBoundingSurface::setTrialStrain(const Vector &v, const Vector &r)
{
	return this->setTrialStrain(v);
}

// send back the strain
const Vector&
J2CyclicBoundingSurface::getStrain()
{
	return m_strain_np1; // -1.0 is for geotechnical sign convention
}

// send back the stress 
const Vector&
J2CyclicBoundingSurface::getStress()
{
	// this->integrate();
	return m_stress_np1; // -1.0 is for geotechnical sign convention
}

// set up for initial elastic
void
J2CyclicBoundingSurface::doInitialTangent()
{
	double lambda = 2.0*m_poiss*m_shear / (1.0 - 2.0*m_poiss);

	Matrix SSOIT(6, 6), FOAT1(6, 6);
	for (unsigned int i = 0; i < 3; i++)
		for (unsigned int j = 0; j < 3; j++)
		{
			SSOIT(i, j) = 1;
		}

	FOAT1(0, 0) = 1.0; FOAT1(1, 1) = 1.0; FOAT1(2, 2) = 1.0;
	FOAT1(3, 3) = 0.5; FOAT1(4, 4) = 0.5; FOAT1(5, 5) = 0.5;

	m_Kelas = m_bulk*SSOIT + 2 * m_shear*FOAT1;
	//opserr << "m_Kelas :" << m_Kelas << endln;

	return;

}

// send back the tangent 
const Matrix&
J2CyclicBoundingSurface::getTangent()
{
	// Force elastic response
	if (m_ElastFlag == 0) {
		return	m_Kelas;

	}
	// ElastoPlastic response
	else {
		// only the symetric part
		double lambda = 2.0*m_poiss*m_shear / (1.0 - 2.0*m_poiss);

		Matrix SSOIT(6, 6), FOAT1(6, 6);
		for (unsigned int i = 0; i < 3; i++)
			for (unsigned int j = 0; j < 3; j++)
			{
				SSOIT(i, j) = 1;
			}

		FOAT1(0, 0) = 1.0; FOAT1(1, 1) = 1.0; FOAT1(2, 2) = 1.0;
		FOAT1(3, 3) = 0.5; FOAT1(4, 4) = 0.5; FOAT1(5, 5) = 0.5;

	m_Ktan_n   = m_bulk*SSOIT + m_xi_n*FOAT1;

	//return m_Ktan_n;
	return m_Kelas;
	}
		
}

// send back the tangent 
const Matrix&
J2CyclicBoundingSurface::getInitialTangent()
{
	return m_Kelas;
}
