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

#ifndef J2CyclicBoundingSurface_h
#define J2CyclicBoundingSurface_h

#include <stdio.h> 
#include <stdlib.h> 
#include <math.h> 

#include <Vector.h>
#include <Matrix.h>
#include <NDMaterial.h>


class J2CyclicBoundingSurface : public NDMaterial {

	//-------------------Declarations-------------------------------

public:

	//null constructor
	J2CyclicBoundingSurface();

	//full constructor
	J2CyclicBoundingSurface(int tag,
		double G,
		double K,
		double su,
		double rho,
		double h,
		double m,
		double h0,
		double chi,
		double beta);

	J2CyclicBoundingSurface(int tag, int classTag,
		double G,
		double K,
		double su,
		double rho,
		double h,
		double m,
		double h0,
		double chi,
		double beta);

	//destructor
	virtual ~J2CyclicBoundingSurface();

	virtual const char *getClassType(void) const { return "J2CyclicBoundingSurface"; };

	virtual NDMaterial* getCopy(const char *type);

	//swap history variables
	virtual int commitState();

	//revert to last saved state
	virtual int revertToLastCommit();

	//revert to start
	virtual int revertToStart();

	//sending and receiving
	virtual int sendSelf(int commitTag, Channel &theChannel);
	virtual int recvSelf(int commitTag, Channel &theChannel,
		FEM_ObjectBroker &theBroker);

	//print out material data
	void Print(OPS_Stream &s, int flag = 0);

	virtual NDMaterial *getCopy(void);
	virtual const char *getType(void) const;
	virtual int getOrder(void) const;

	double getRho(void) { return m_density; }



	virtual int setParameter(const char **argv, int argc, Parameter &param);
	virtual int updateParameter(int responseID, Information &info);
	virtual int activateParameter(int paramID);

	virtual const Matrix& getDampTangent();

	virtual int setTrialStrain(const Vector &strain_from_element);
	virtual int setTrialStrain(const Vector &v, const Vector &r);

	// send back the strain
	virtual const Vector&	 getStrain();

	// send back the stress 
	virtual const Vector&	  getStress();

	// send back the tangent 
	virtual const Matrix&	 getTangent();

	// send back the tangent 
	virtual const Matrix&  getInitialTangent();


protected:
	static char unsigned m_ElastFlag;	// 1: enforce elastic response

//material parameters
	double m_su;          // undrained shear strength
	double m_R;           // radius of the bounding surface
	double m_bulk;        // bulk modulus
	double m_shear;       // shear modulus
	double m_density;     // material density
	double m_h_par;       // exponential hardening parameter
	double m_m_par;       // exponential hardening parameter
	double m_h0_par;
	double m_beta;        // Integration scheme parameter beta = 0, explicit. beta = 1, implicit. beta = 0.5 mid point rule


	//internal variables
	Vector m_sigma0_n;           // sigma0 time n
	Vector m_sigma0_np1;         // sigma0 time n+1
	double m_kappa_n;            // kappa  time n
	double m_kappa_np1;          // kappa  time n+1
	double m_psi_n;              // hardening variable time n
	double m_psi_np1;            // hardening variable time n+1
	double m_kappa_inf;          // kappa inf  

	//material response 
	Vector m_stress_n;           //stress vector time n
	Vector m_stress_np1;         //stress vector time n+1
	Vector m_stress_vis_n;
	Vector m_stress_vis_n1;
	Vector m_stress_t_n1;
	Matrix m_Cep;
	Matrix m_Ce;

	// viscous behavior
	double m_chi;
	Matrix m_D;

	//material input
	Vector m_strain_n;           //strain vector time n
	Vector m_strain_np1;         //strain vector time n+1
	Vector m_strainRate_n;
	Vector m_strainRate_n1;

	bool m_isElast2Plast;

	//zero internal variables
	void zero();

	//plasticity integration routine
	void	integrate();
	void	elastic_integrator();
	void	plastic_integrator();
	void    viscoElastic_integrator();

	void    calcInitialTangent();

	//hardening function
	double H(double kappa);

	double trace(Vector V);

	Vector getDevPart(Vector V);

	double inner_product(Vector x, Vector y, int type);

	double vector_norm(Vector x, int type);

	Vector convert_to_stressLike(Vector v);

	double small;
	bool debugFlag;

}; //end of J2CyclicBoundingSurface declarations
#endif
