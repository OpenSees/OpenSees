/* *********************************************************************
**    Module:	PyLiq1.h 
**
**    Purpose:	Provide a p-y material that gets pore pressure from a
**				specified element that contains a PorousFluidSolid.
**              
**
**    Developed by Ross W. Boulanger
**
** Copyright @ 2002 The Regents of the University of California (The Regents). All Rights Reserved.
**
** The Regents grants permission, without fee and without a written license agreement, for (a) use, 
** reproduction, modification, and distribution of this software and its documentation by educational, 
** research, and non-profit entities for noncommercial purposes only; and (b) use, reproduction and 
** modification of this software by other entities for internal purposes only. The above copyright 
** notice, this paragraph and the following three paragraphs must appear in all copies and modifications 
** of the software and/or documentation.
**
** Permission to incorporate this software into products for commercial distribution may be obtained 
** by contacting the University of California 
** Office of Technology Licensing 
** 2150 Shattuck Avenue #510, 
** Berkeley, CA 94720-1620, 
** (510) 643-7201.
**
** This software program and documentation are copyrighted by The Regents of the University of California. 
** The Regents does not warrant that the operation of the program will be uninterrupted or error-free. The 
** end-user understands that the program was developed for research purposes and is advised not to rely 
** exclusively on the program for any reason.
**
** IN NO EVENT SHALL REGENTS BE LIABLE TO ANY PARTY FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR 
** CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS, ARISING OUT OF THE USE OF THIS SOFTWARE AND ITS 
** DOCUMENTATION, EVEN IF REGENTS HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. REGENTS GRANTS 
** NO EXPRESS OR IMPLIED LICENSE IN ANY PATENT RIGHTS OF REGENTS BUT HAS IMPLEMENTED AN INDIVIDUAL 
** CONTRIBUTOR LICENSE AGREEMENT FOR THE OPENSEES PROJECT AT THE UNIVERSITY OF CALIFORNIA, BERKELEY 
** TO BENEFIT THE END USER.
**
** REGENTS SPECIFICALLY DISCLAIMS ANY WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES 
** OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. THE SOFTWARE AND ACCOMPANYING DOCUMENTATION,
** IF ANY, PROVIDED HEREUNDER IS PROVIDED "AS IS". REGENTS HAS NO OBLIGATION TO PROVIDE MAINTENANCE, 
** SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS.
**
** ****************************************************************** */

// $Revision: 1.0
// $Date: 2002/5/15
// $Source: /OpenSees/SRC/material/uniaxial/PyLiq1.h

#ifndef PYLIQ1_H
#define PYLIQ1_H

// Written: RWB
// Created: May 2002
//
// Description: This file contains the class definition for PyLiq1.
// 

#include <UniaxialMaterial.h>
#include <Domain.h>
//#include <InitialStateAnalysisWrapper.h>
#include "PySimple1.h"
#include <iostream>
#include <string>

class TimeSeries;

class PyLiq1 : public PySimple1
{
  public:
  PyLiq1(int tag, int classtag, int soilType, double pult, double y50, double drag,
	 double dashpot, double pRes, int solidElem1, int solidElem2, Domain *theDomain);
  PyLiq1(int tag, int classtag, int soilType, double pult, double y50, double drag,
	 double dashpot, double pRes, Domain *theDomain, TimeSeries *theSeries);
  PyLiq1();
  ~PyLiq1();

  int setTrialStrain(double y, double yRate); 
  double getStrain(void);          
  double getStress(void);
  double getTangent(void);
  double getStrainRate(void);
  double getDampTangent(void);
  double getInitialTangent(void);
  
  int commitState(void);
  int revertToLastCommit(void);    
  int revertToStart(void);        
  
  UniaxialMaterial *getCopy(void);
  
  int sendSelf(int commitTag, Channel &theChannel);  
  int recvSelf(int commitTag, Channel &theChannel, 
	       FEM_ObjectBroker &theBroker);
  
  //  Command for initiating vertConsolStress from TclUpdateMaterialStageCommand
  int updateParameter(int snum, Information &eleInformation);
  int setParameter(const char **argv, int argc, Parameter &param);
  
  
  void Print(OPS_Stream &s, int flag =0);
  
  protected:
    
  private:

  // Residual p (other parameters in PySimple1 base class)
  double pRes;
  
  // Committed and trial values for p, y, and ru
  double Ty;
  double Cy;
  double Tp;
  double Cp;
  double Tangent;
  double maxTangent;
  double Tru;
  double Cru;
  double Hru;
  
  // Solid element from which pore pressures are obtained, domain pointer
  // and stage information to get the initial vertical effective stress.
  int solidElem1, solidElem2, theSeriesTag;
  double meanConsolStress, meanStress;
  static int loadStage;
  int lastLoadStage;
  std::string elemFlag; 
  Domain *theDomain;
  TimeSeries *theSeries;
  
  // Initial tangent
  double initialTangent;
  
  // Function for obtaining effective stresses from adjoining solid soil elements
  double getEffectiveStress(void);
  double getEffectiveStress(TimeSeries *theSeries);
  static Vector stressV3;
};

#endif // PYLIQ1_H



