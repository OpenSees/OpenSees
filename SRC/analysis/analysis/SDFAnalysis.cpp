/* ****************************************************************** **
**    OpenSees System for Earthquake Engineering Simulation           **
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

#include <PathSeries.h>

#include <math.h>
#include <elementAPI.h>

int OPS_sdfResponse()
{
  int numData = 1;
  double uresidual = 0.0;
  double umaxprev = 0.0;
  
  int numOptionalArgs = 0;
  int numArgs = OPS_GetNumRemainingInputArgs();
  while (OPS_GetNumRemainingInputArgs() > 0) {
    std::string type = OPS_GetString();
    if (type == "-uresid" || type == "-uresidual") {
      numOptionalArgs++;
      if (OPS_GetNumRemainingInputArgs() > 0) {
	numData = 1;
	if (OPS_GetDoubleInput(&numData, &uresidual) < 0) {
	  opserr << "ERROR sdfResponse - failed to read uresidual" << endln;
	  return 0;
	}
	numOptionalArgs++;
      }
    }
    if (type == "-umaxprev") {
      numOptionalArgs++;
      if (OPS_GetNumRemainingInputArgs() > 0) {
	numData = 1;
	if (OPS_GetDoubleInput(&numData, &umaxprev) < 0) {
	  opserr << "ERROR sdfResponse - failed to read umaxprev" << endln;
	  return 0;
	}
	numOptionalArgs++;
      }
    }
  }

  if (numArgs > 0) {
    OPS_ResetCurrentInputArg(-numArgs);
  }
  numArgs = numArgs - numOptionalArgs;
    
    if (numArgs < 7 || numArgs > 8) {
	opserr << "Incorrect number of arguments to sdfResponse --";
	opserr << "m, zeta, k, Fy, alpha, dtF, filename, dt, <-uresidual, uresid, -umaxprev, umaxp>" << endln;
	opserr << "m, zeta, k, Fy, alpha, tsTag, dt, <-uresidual, uresid, -umaxprev, umaxp>" << endln;	
	return -1;
    }

    double data[5] = {0,0,0,0,0};
    numData = 5;
    if (OPS_GetDoubleInput(&numData, &data[0]) < 0) {
	opserr << "WARNING sdfResponse -- could not read input \n";
	return -1;
   } 

    bool timeSeries = false;
    TimeSeries *accelSeries = 0;
    // Time series tag specified
    if (numArgs == 7) {
      int tsTag = 0;
      numData = 1;
      if (OPS_GetIntInput(&numData, &tsTag) < 0) {
	opserr << "WARNING sdfResponse -- could not read timeSeries tag" << endln;
	return -1;
      }

      accelSeries = OPS_getTimeSeries(tsTag);
      if (accelSeries == 0) {
	opserr << "WARNING invalid accel series: " << tsTag;
	opserr << " sdfResponse\n";
	return -1;
      }
      
      timeSeries = true;
    }

    // dtF and filename, create internal time series
    if (numArgs == 8) {
      double dtF = 0.0;

      numData = 1;
      if (OPS_GetDoubleInput(&numData, &dtF) < 0) {
	opserr << "WARNING sdfResponse -- could not read input \n";
	return -1;
      }

      const char *filename = OPS_GetString();

      accelSeries = new PathSeries(0, filename, dtF);      
    }

    if (accelSeries == 0) {
      opserr << "ERROR sdfResponse -- failed to allocate time series" << endln;
      if (!timeSeries && accelSeries != 0)
	delete accelSeries;
      return -1;
    }
    
    double dt = 0.0;
    numData = 1;
    if (OPS_GetDoubleInput(&numData, &dt) < 0) {
	opserr << "WARNING sdfResponse -- could not read dt input \n";
	return -1;
    }

    
    double m = data[0];
    double zeta = data[1];
    double k = data[2];
    double Fy = data[3];
    double alpha = data[4];
    //double dtF = data[5];
    //double dt = data[6];
    //double uresidual = data[7];
    //double umaxprev = data[8];
    
    double gamma = 0.5;
    double beta = 0.25;
    double tol = 1.0e-8;
    int maxIter = 10;

    double c = zeta*2*sqrt(k*m);
    double Hkin = alpha/(1.0-alpha)*k;

    double p0 = 0.0;
    double u0 = uresidual;
    double v0 = 0.0;
    double fs0 = 0.0;
    double a0 = (p0-c*v0-fs0)/m;

    double a1 = m/(beta*dt*dt) + (gamma/(beta*dt))*c;
    double a2 = m/(beta*dt) + (gamma/beta-1.0)*c;
    double a3 = (0.5/beta-1.0)*m + dt*(0.5*gamma/beta-1.0)*c;

    double au = 1.0/(beta*dt*dt);
    double av = 1.0/(beta*dt);
    double aa = 0.5/beta-1.0;

    double vu = gamma/(beta*dt);
    double vv = 1.0-gamma/beta;
    double va = dt*(1-0.5*gamma/beta);
    
    double kT0 = k;

    double umax = fabs(umaxprev);
    double amax = 0.0; double tamax = 0.0;
    double up = uresidual; double up0 = up;
    int i = 0;
    double ft, u=0, du, v, a, fs, zs, ftrial, kT, kTeff, dg, phat, R, R0, accel;
    double time = accelSeries->getStartTime();
    double Tend = accelSeries->getDuration();
    while (time < Tend) {

      ft = accelSeries->getFactor(time);
      
	i++;
    
	u = u0;
      
	fs = fs0;
	kT = kT0;
	up = up0;
      
	phat = ft + a1*u0 + a2*v0 + a3*a0;
      
	R = phat - fs - a1*u;
	R0 = R;
	if (R0 == 0.0) {
	    R0 = 1.0;
	}
    
	int iter = 0;

	while (iter < maxIter && fabs(R/R0) > tol) {
	    iter++;

	    kTeff = kT + a1;

	    du = R/kTeff;

	    u = u + du;

	    fs = k*(u-up0);
	    zs = fs-Hkin*up0;
	    ftrial = fabs(zs)-Fy;
	    if (ftrial > 0) {
		dg = ftrial/(k+Hkin);
		if (fs < 0) {
		    fs = fs + dg*k;
		    up = up0 - dg;
		} else {
		    fs = fs - dg*k;
		    up = up0 + dg;
		}
		kT = k*Hkin/(k+Hkin);
	    } else {
		kT = k;
	    }
      
	    R = phat - fs - a1*u;
	}

	v = vu*(u-u0) + vv*v0 + va*a0;
	a = au*(u-u0) - av*v0 - aa*a0;

	u0 = u;
	v0 = v;
	a0 = a;
	fs0 = fs;
	kT0 = kT;
	up0 = up;

	if (fabs(u) > umax) {
	    umax = fabs(u);
	}
	if (fabs(a) > amax) {
	    amax = fabs(a);
	    tamax = i*dt;
	}

	time += dt;
    }
  
    double output[] = {umax, u, up, amax, tamax};
    numData = 5;
    if (OPS_SetDoubleOutput(&numData,output, false) < 0) {
	opserr << "WARNING: failed to set output -- sdfResponse\n";
	return -1;
    }

    if (!timeSeries && accelSeries != 0)
      delete accelSeries;
    
    return 0;
}
