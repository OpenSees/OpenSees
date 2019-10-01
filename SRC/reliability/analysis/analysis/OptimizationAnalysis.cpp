/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 2001, The Regents of the University of California    **
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
** Reliability module developed by:                                   **
**   Terje Haukaas (haukaas@ce.berkeley.edu)                          **
**   Armen Der Kiureghian (adk@ce.berkeley.edu)                       **
**                                                                    **
** ****************************************************************** */
                                                                        
// $Revision: 1.2 $
// $Date: 2006-12-06 22:32:23 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/analysis/OptimizationAnalysis.cpp,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu)
//

#include <OptimizationAnalysis.h>
#include <ReliabilityDomain.h>
#include <ReliabilityAnalysis.h>
#include <Vector.h>
#include <Matrix.h>
#include <MatrixOperations.h>
#include <NormalRV.h>
#include <RandomVariable.h>
#include <math.h>

#include <fstream>
#include <iomanip>
#include <iostream>
using std::ifstream;
using std::ios;
using std::setw;
using std::setprecision;
using std::setiosflags;


OptimizationAnalysis::OptimizationAnalysis(ReliabilityDomain *passedReliabilityDomain,
						   TCL_Char *passedFileName)
:ReliabilityAnalysis()
{
	theReliabilityDomain = passedReliabilityDomain;
	strcpy(fileName,passedFileName);
}


OptimizationAnalysis::~OptimizationAnalysis()
{
  
}



int 
OptimizationAnalysis::analyze(void)
{

	// Alert the user that the FORM analysis has started
	opserr << "OptimizationAnalysis is running ... " << endln;


	// Assumptions
	// - Have start values for the design parameters x


	// Algorithm parameters
	double alphapar = 0.5;
	double betapar = 0.8;
	double gammapar = 2.0;
	double deltapar = 1.0;

	// Bound on probability of failure
	double p0 = 0.00134989803163;
	double beta0 = 3.0;

	// Related probability quantities...
	double pmax = p0;
	double pmin = 0.00001*p0;
	double betavarmax = -5.56;
	double betavarmin = beta0;

	double covtarget = 0.11;

	// Number of limit-state functions
	int K = 4; 

	// Declare 'tvector'
	Vector tvector(K,1);
	for (int k=0; k<K; k++) {
		for (int j=0; j<K; j++) {
			tvector(k,j) = 1.0;
		}
	}

	// Even more parameters...
	double targetcost = 0.0;
	int maxouteriter = 25;  
	int maxinneriter = 150;
	int N = 50000;


	// Number of random variables
	int size_u = 8; 

	// Discretization of u-space (one matrix per limit-state function!)
	Matrix Omega1(1,size_u);
	Matrix Omega2(1,size_u);
	Matrix Omega3(1,size_u);
	Matrix Omega4(1,size_u);

	double betavar = beta0;

	// Evaluate cost (objective function)
	cost = f0(x,betavar); 

	// Evaluate deterministic constraints
	convio = fjs(x,betavar,betavarmax,betavarmin);

	termflag = 1;
	for (int i=1; i<1; i++) {
   
		int ok = mooa(x_out, betavar, maxconstr, maxoverballs,u,
					 x,
					 alphapar,
					 betapar,
					 gammapar,
					 deltapar,
					 tvector,
					 Omega,
					 maxinneriter,
					 size_u,
					 assump,
					 betavar,
					 betavarmax,
					 betavarmin);



		// compute prob of failure
		[pf,beta,cov_pf,numsim,cputime,counter] = impsamp(x,u,i*N,0,covtarget - (0.004*i),assump);


		// Heuristic update of parameters t  
		for (int kk=0; kk<K; kk++) {
			tvector(kk) = tvector(kk)*betavar/beta;
		}

	}


	// Open output file
	ofstream outputFile( fileName, ios::out );

	outputFile << "#######################################################################" << endln;
	outputFile << "#  FORM ANALYSIS RESULTS, LIMIT-STATE FUNCTION NUMBER "
		<<setiosflags(ios::left)<<setprecision(1)<<setw(4)<<lsf <<"            #" << endln;
	outputFile << "#                                                                     #" << endln;
	outputFile << "#  No convergence!                                                    #" << endln;
	outputFile << "#                                                                     #" << endln;
	outputFile << "#######################################################################" << endln << endln << endln;

	// Clean up
	outputFile.close();

	// Print summary of results to screen (more here!!!)
	opserr << "OptimizationAnalysis completed." << endln;

	return 0;
}








double 
OptimizationAnalysis::f0(Vector x, double betavar)
{
	double Cs = 50.0;
	double Cc = 1.0;
	double alpha = 0.1; 
	double Y3 = 18.30; 

	double X19 = 2*Y3*(1/x(7-1) + 1/x(8-1) + 1/x(9-1))/6;

	double initialcost = 3*Cs*Y3*x(1-1)/4 + Cs*X19*x(6-1)*(x(3-1) + x(5-1) - alpha + 0.5*x(4-1)) + Cc*Y3*(x(2-1)*x(3-1) + x(4-1)*x(5-1));

	double costoffailure = 500*initialcost;

	double cost = initialcost + costoffailure*normcdf(-betavar);

	return cost;
}


Vector 
OptimizationAnalysis::fjs(Vector x, double betavar, double betavarmax, double betavarmin)
{

	double alpha = 0.1;

	double X20b = (x(3-1) + x(5-1) - alpha)/2.0;
	double X20c = 24.0*0.0254;

	Vector fjvector(25);

	fjvector(1-1) = x(7-1) - X20b; 
	fjvector(2-1) = x(7-1) - X20c; 
	fjvector(3-1) = x(8-1) - X20b; 
	fjvector(4-1) = x(8-1) - X20c;
	fjvector(5-1) = x(9-1) - X20b; 
	fjvector(6-1) = x(9-1) - X20c;

	fjvector(7-1) = 0.5*x(4-1) - x(3-1);  
	fjvector(8-1) = x(2-1) - 4*x(4-1);
	fjvector(9-1) = x(4-1) - x(2-1);
	fjvector(10-1) = x(2) - 1.22;
	fjvector(11-1) = 0.15 - x(3-1);
	fjvector(12-1) = 0.15 - x(4-1);
	fjvector(13-1) = x(5-1)/x(4-1) - 4;
	fjvector(14-1) = 1 - x(6-1)/0.0001;
	fjvector(15-1) = -x(7-1);
	fjvector(16-1) = -x(8-1);
	fjvector(17-1) = -x(9-1);
	fjvector(18-1) = x(3-1) + x(5-1) - 1.2;

	fjvector(19-1) = -g(x,zeros(8,1),9,assump);

	// Assuming tension force in the steel is balanced by a concrete compression zone IN THE FLANGE
//	fjvector(20-1) = -g(x,zeros(8,1),5,assump);
//	fjvector(21-1) = -g(x,zeros(8,1),7,assump);
//	fjvector(22-1) = g(x,zeros(8,1),10,assump);          

	// Assuming tension force in steel is balanced by a concrete compression zone that propagates into the THE WEB
	fjvector(20-1) = -g(x,zeros(8,1),6,assump);
	fjvector(21-1) = -g(x,zeros(8,1),8,assump);
	fjvector(22-1) = g(x,zeros(8,1),11,assump);


	fjvector(23-1) = -x(5-1);
	fjvector(24-1) = betavar - betavarmax;
	fjvector(25-1) = -betavar + betavarmin;


}


double 
OptimizationAnalysis::g(Vector x, Vector u, int k) 
{

	// Probability transformations
	Vector v(8);
	v(1-1) = u(1-1)*0.15*413.4*pow(10.0,6.0) + 413.4*pow(10.0,6.0);       
	v(2-1) = u(2-1)*0.15*27.56*pow(10.0,6.0) + 27.56*pow(10.0,6.0);   
	v(3-1) = u(3-1)*0.2*13.57*pow(10.0,3.0) + 13.57*pow(10.0,3.0);      
	v(4-1) = u(4-1)*0.243*929*pow(10.0,3.0) + 929*pow(10.0,3.0);      
	v(5-1) = u(5-1)*0.243*138.31*pow(10.0,3.0) + 138.31*pow(10.0,3.0); 
	v(6-1) = u(6-1)*0.243*183.39*pow(10.0,3.0) + 183.39*pow(10.0,3.0); 
	v(7-1) = u(7-1)*0.243*228.51*pow(10.0,3.0) + 228.51*pow(10.0,3.0); 
	v(8-1) = u(8-1)*0.1*22.74*pow(10.0,3.0) + 22.74*pow(10.0,3.0);     

	double Y3 = 18.30;
	double alpha = 0.1;

	double X15 = (2*x(4-1)/0.0254*(x(3-1)/0.0254 + x(5-1)/0.0254 - alpha/0.0254)*sqrt(v(2-1)/(6.89*pow(10.0,3.0))))*4.45; 
	double X16 = x(6-1)*v(1-1)*(x(3-1) + x(5-1) - alpha)/x(7-1);
	double X17 = x(6-1)*v(1-1)*(x(3-1) + x(5-1) - alpha)/x(8-1);
	double X18 = x(6-1)*v(1-1)*(x(3-1) + x(5-1) - alpha)/x(9-1);


	double gval;

	switch (k) {
	case 1:
		double X14 = 0.85*v(2-1)*(x(2-1) - x(4-1))*x(3-1)/v(1-1);
		double X10 = (x(1-1) - X14)*v(1-1)/(0.85*v(2-1)*x(4-1));
		double X11 = (x(1-1) - X14)*v(1-1)*((x(3-1) + x(5-1) - alpha) - X10/2) + X14*v(1-1)*((x(3-1) + x(5-1) - alpha) - x(3-1)/2.0);
		gval = 1 - v(4-1)/X11 - (v(3-1)*Y3^2/8)/X11 - ((x(2-1)*x(3-1) + x(4-1)*x(5-1))*v(8-1)*Y3^2/8)/X11;
		break;
	case 2:  
		double cap = X15 + X16;
		gval = 1 - v(5-1)/cap - (v(3-1)*Y3/6)/cap - ((x(2-1)*x(3-1) + x(4-1)*x(5-1))*v(8-1)*Y3/6)/cap;
		break;
	case 3:
		double cap = X15 + X17;
		gval = 1 - v(6-1)/cap - (v(3-1)*Y3/3)/cap - ((x(2-1)*x(3-1) + x(4-1)*x(5-1))*v(8-1)*Y3/3)/cap;
		break;
	case 4:
		double cap = X15 + X18;
		gval = 1 - v(7-1)/cap - (v(3-1)*Y3/2)/cap - ((x(2-1)*x(3-1) + x(4-1)*x(5-1))*v(8-1)*Y3/2)/cap;
		break;
	case 5:
		double X12 = x(1-1)/(x(2-1)*(x(3-1) + x(5-1) - alpha));
		double X13 = (0.85*0.85*v(2-1)/v(1-1))*(87/(87 + v(1-1)/(6.89*pow(10.0,6.0))));
		gval = -X12 + 0.75*X13;
		break;
	case 6:
		double X12 = x(1-1)/(x(4-1)*(x(3-1) + x(5-1) - alpha));
		double X14 = 0.85*v(2-1)*(x(2-1) - x(4-1))*x(3-1)/v(1-1); 
		double X13 = ((0.85*0.85*v(2-1)/v(1-1))*87/(87 + v(1-1)/(6.89*pow(10.0,6.0))) + X14/(x(4-1)*(x(3-1) + x(5-1) - alpha)))*x(4-1)/x(2-1);
		gval = -X12 + 0.75*X13;
		break;
	case 7:
		double X21 = 200/(1000*v(1-1)/(6.89*pow(10.0,6.0)));  
		double X12 = x(1-1)/(x(2-1)*(x(3-1) + x(5-1) - alpha));
		gval = -X21 + X12;
		break;
	case 8:
		double X21 = 200/(1000*v(1-1)/(6.89*pow(10.0,6.0)));  
		double X12 = x(1-1)/(x(4-1)*(x(3-1) + x(5-1) - alpha));
		gval = -X21 + X12;
		break;
	case 9:
		gval = -x(6-1)*(v(1-1)/pow(10.0,3.0)/6.89)/x(9-1)/2/x(4-1)*(1/sqrt(v(2-1)/(6.89*pow(10.0,3.0)))) + 4; 
		break;
	case 10:
		gval = 1 - (0.85*v(2-1)*x(2-1)*x(3-1))/(v(1-1)*x(1-1));
		break;
	case 11:
		gval = -1 + (0.85*v(2-1)*x(2-1)*x(3-1))/(v(1-1)*x(1-1));
		break;
	default:
		break;
	}

	return gval;

}



int 
OptimizationAnalysis::mooa(Vector x_out, double betavar, double maxconstr, double maxoverballs, Vector u,
						   Vector x,
						   double alpha,
						   double beta,
						   double gammapar,
						   double delta,
						   Vector tvector,
						   Matrix Omega,
						   int maxiter,
						   int size_u,
						   int assump,
						   double betavar,
						   double betavarmax,
						   double betavarmin) 
{





Omega_basic = Omega;

epsilon = 0.1/(N*N);

K = length(tvector); // Number of random variables...

innersol = [];
psi_0_bar = [];
psi_omega = [];

Omeganames = fieldnames(Omega);

	Vector uinitial(nrv);


	for (int N=1; N<maxiter+1; N++) {
   
	   // Compute approximate solution to inner problem(s)
	   eps_inner = epsilon/K;
   
	   // For each limit-state function...
//	   for (int kk=1; kk<=K; kk++) {
		  innersol_temp = phinner(x,
								  uinitial,
								  alpha,
								  beta,
								  gammapar,
								  delta,
								  tvector(kk),
								  eps_inner,
								  kk,
								  assump,
								  betavar);
//	   }
   
	   uinitial = innersol_temp(1:m,:);  
	   innersol = [innersol innersol_temp];
   
	   // Construct new set omega_N
	   fjvalues = fjs(x,assump,betavar,betavarmax,betavarmin);
	   psi_0_bar = [psi_0_bar max(fjvalues)];
   
	   Omega = Omega_basic;

	   eps_vector = eps_constdrop(N);

	   psi_omega = [psi_omega  max([0; max(innersol(m + 1,(N - 1)*K + 1:N*K)); psi_0_bar(N)])];

   
	   for n = 1:N

		  psival = psi_omega(n);

		  for kk = 1:K

			 ctr = innersol(size_u + 1,(n - 1)*K + kk);

			 if abs(ctr - psival) < 10^(-15) & psival > eps_vector(n) 

				kfield = char(Omeganames(kk));

				oldvalues = getfield(Omega,kfield);

				Omega = setfield(Omega,kfield,[oldvalues innersol(1:m,(n - 1)*K + kk)]);

			 end
		  end
	   end
   
	   [x,betavar,constr,theta,maxoverballs] = 
	   phouter(x,
			   alpha,
			   beta,
			   gammapar,
			   Omega,
			   K,
			   epsilon,
			   assump,
			   betavar,
			   betavarmax,
			   betavarmin,
			   tvector);
   
   
   
	   newres = [x;betavar;constr;maxoverballs;theta];
      
	   H = [H newres];
   
	   N = N + 1;
	   epsilon = eps_outer(N);
   
	}

last_u = uinitial;

}


