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
                                                                        
// $Revision: 1.0 $
// $Date: 2012-05-24 22:03:16 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/ConcreteS.cpp,v $

//
// Yuli Huang (yulihuang@gmail.com) & Xinzheng Lu (luxz@tsinghua.edu.cn)
//
// Simple Concrete Model, plane stress
// 


#include <ConcreteS.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <MaterialResponse.h>
#include <elementAPI.h>

void * OPS_ADD_RUNTIME_VPV(OPS_ConcreteS)
{
    int argc = OPS_GetNumRemainingInputArgs() + 2;
    if (argc < 8) {
	opserr << "WARNING insufficient arguments\n";
	opserr << "Want: nDMaterial ConcreteS tag? E? nu? fc? ft? Es?" << endln;
	return 0;
    }

    int tag;
    int numdata = 1;
    if (OPS_GetIntInput(&numdata, &tag) < 0) {
	opserr << "WARNING invalid nDMaterial ConcreteS tag" << endln;
	return 0;
    }

    // double E, nu, fc, ft, Es;
    double data[5];
    numdata = 5;
    if (OPS_GetDoubleInput(&numdata, data) < 0) {
	opserr << "WARNING invalid double inputs" << endln;
	opserr << "ConcreteS: " << tag << endln;
	return 0;
    }

    return new ConcreteS( tag, data[0], data[1], data[2], data[3], data[4]);
}

//null constructor
ConcreteS::ConcreteS( ) : 
NDMaterial(0, ND_TAG_ConcreteS ), 
strain0(3), strain(3), stress0(3), stress(3), stressd(3),
tangent(3,3),eTangent(3,3),
cStrain0(0),cStrain(0)
{ }


//full constructor
ConcreteS::ConcreteS(int tag, double rE, double rnu, double rfc, double rft, double rEs) :
NDMaterial( tag, ND_TAG_ConcreteS ),
strain0(3), strain(3), stress0(3), stress(3), stressd(3),
tangent(3,3),eTangent(3,3),
cStrain0(0),cStrain(0),
E(rE), nu(rnu), fc(fabs(rfc)), ft(rft), Es(fabs(rEs))
{
  setInitials();
}


//destructor
ConcreteS::~ConcreteS( ) 
{ 

} 

void ConcreteS::setInitials()
{
  double fac;
  eTangent.Zero();
  eTangent(0,0) = 1.0;
  eTangent(0,1) = nu;
  eTangent(1,0) = nu;
  eTangent(1,1) = 1.0;
  eTangent(2,2) = 0.5 * (1.0 - nu);
  fac = E / (1.0 - nu * nu);
  eTangent *= fac;
  tangent = eTangent;
  Ep = 10 * E;
  EmEp1 = 1.0 / (E + Ep);
  beta = - (Ep + 11 * Es) / ft;
//  Ep = E * Es / (E + Es)
//  EmEp1 = 1.0 / (E - Ep) = (E + Es) / E^2;
//  EmEp1 = (E + Es) / (E * E);
//  beta = -(E + 2 * Es) / ft;
//  Es = -E/2;
//  EmEp1 = 0.5 / E;
  ftmin = 1.0e-3 * ft;
}

//make a clone of this material
NDMaterial*
ConcreteS::getCopy( ) 
{
  ConcreteS *clone ;   //new instance of this class

  clone = new ConcreteS( this->getTag(), E, nu, fc, ft, Es);

  return clone ;
}


//make a clone of this material
NDMaterial* 
ConcreteS::getCopy( const char *type ) 
{
  return this->getCopy( ) ;
}


//send back order of strain in vector form
int 
ConcreteS::getOrder( ) const
{
  return 3 ;
}


const char*
ConcreteS::getType( ) const 
{
  return "ConcreteS" ; 
}



//swap history variables
int 
ConcreteS::commitState( ) 
{
  stress0 = stress;
  strain0 = strain;
  cStrain0 = cStrain;
  return 0;
}


//revert to last saved state
int 
ConcreteS::revertToLastCommit( )
{
  return 0;
}


//revert to start
int
ConcreteS::revertToStart( )
{
  strain0.Zero();
  strain.Zero();
  stress0.Zero();
  stress.Zero();
  cStrain0 = 0.0;
  cStrain = 0.0;
  return 0;
}

//receive the strain
int 
ConcreteS::setTrialStrain( const Vector &strainFromElement )
{
  static Matrix CfCf(3,3);
  static Vector flow(3), Cf(3);
  double vStress, vStress1, yieldFunc;
  double fCf, sigm, sigd, theta;
  double ps1, ps2, psmax, tStrain, eps;
  strain(0) = strainFromElement(0) ;
  strain(1) = strainFromElement(1) ;
  strain(2) = strainFromElement(2) ;

  stress = stress0 + eTangent * (strain - strain0);
  
  tangent = eTangent;
  
  vStress = sqrt(  stress(0) * stress(0)
                 - stress(0) * stress(1)
                 + stress(1) * stress(1)
                 + stress(2) * stress(2) * 3.0);
  
  yieldFunc = vStress - fc;
  sigm = 0.5 * (stress(0) + stress(1));
  sigd = 0.5 * (stress(0) - stress(1));
  ps2 = sqrt(sigd * sigd + stress(2) * stress(2));
  ps1 = sigm + ps2;
  ps2 = sigm - ps2;
  if (yieldFunc > 0.0 && ps1 <= ft && ps2 <= ft)
  {
    vStress1 = 1.0 / vStress;
    flow(0) = (stress(0) - 0.5 * stress(1)) * vStress1;
    flow(1) = (stress(1) - 0.5 * stress(0)) * vStress1;
    flow(2) = 3.0 * stress(2) * vStress1;
    Cf = eTangent * flow;
    fCf = flow^Cf;
    eps = yieldFunc / fCf;
    stress -= eTangent * (eps * flow);
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++)
        CfCf(i,j) = Cf(i) * Cf(j);
    tangent -= CfCf / (fCf + 1.e-3 * E);
  }
  sigm = 0.5 * (stress(0) + stress(1));
  sigd = 0.5 * (stress(0) - stress(1));
  theta = atan2(stress(2), sigd);
  ps2 = sqrt(sigd * sigd + stress(2) * stress(2));
  ps1 = sigm + ps2;
  ps2 = sigm - ps2;
  
  if (ps2 > 0.0 && ps1 < -fc) ps1 = -fc;
  if (ps1 > 0.0 && ps2 < -fc) ps2 = -fc;
  psmax = ft + Ep * cStrain0;
//  if (psmax < ftmin) psmax = ftmin;
  yieldFunc = ps1 - psmax;

  if (yieldFunc > 0.0)
  {
    cStrain = cStrain0 + yieldFunc * EmEp1;
    ps1 = ft + Ep * cStrain;
  }

  yieldFunc = ps2 - psmax;

  if (yieldFunc > 0.0)
  {
    tStrain = cStrain0 + yieldFunc * EmEp1;
    ps2 = ft + Ep * tStrain;
    if (cStrain < tStrain) cStrain = tStrain;
  }
  
  sigm = 0.5 * (ps1 + ps2);
  sigd = 0.5 * fabs(ps1 - ps2);
  
  stress(1) = sigd * cos(theta);
  stress(0) = sigm + stress(1);
  stress(1) = sigm - stress(1);
  stress(2) = sigd * sin(theta);

  stressd = stress;   
  
  double damage = exp(beta * cStrain);
  if (ps1 > 0.0) ps1 *= damage;
  if (ps2 > 0.0) ps2 *= damage;

  sigm = 0.5 * (ps1 + ps2);
  sigd = 0.5 * fabs(ps1 - ps2);
  
  stressd(1) = sigd * cos(theta);
  stressd(0) = sigm + stressd(1);
  stressd(1) = sigm - stressd(1);
  stressd(2) = sigd * sin(theta);

//  if (ps1 > 0.0 || ps2 > 0.0)
//  {
//    double damage = exp(beta * cStrain);
//    stressd *= damage;
////    tangent *= damage;
//  }

  return 0;
}


//send back the strain
const Vector& 
ConcreteS::getStrain( )
{
  return strain ;
}


//send back the stress 
const Vector&  
ConcreteS::getStress( )
{
  return stressd ;
}


//send back the tangent 
const Matrix&  
ConcreteS::getTangent( )
{
  return tangent ;
}

const Matrix&  
ConcreteS::getInitialTangent
( )
{
  return eTangent ;
}


//print out data
void  
ConcreteS::Print( OPS_Stream &s, int flag )
{
  s << "ConcreteS Material tag: " << this->getTag() << endln ; 
  s << "  E:  " << E  << " ";
  s << "  nu: " << nu << " ";
  s << "  fc: " << fc << " ";
  s << "  ft: " << ft << " ";
  s << "  Es: " << Es << " ";

}


int 
ConcreteS::sendSelf(int commitTag, Channel &theChannel) 
{
  int res = 0, cnt = 0;

  static Vector data(13);

  data(cnt++) = this->getTag();
  data(cnt++) = E;
  data(cnt++) = nu;
  data(cnt++) = fc;
  data(cnt++) = ft;
  data(cnt++) = Es;
  data(cnt++) = cStrain0;

  int i;
  for (i = 0; i < 3; i++) 
    data(cnt++) = strain0(i);

  for (i = 0; i < 3; i++) 
    data(cnt++) = stress0(i);

   res = theChannel.sendVector(this->getDbTag(), commitTag, data);
   if (res < 0) 
      opserr << "ConcreteS::sendSelf() - failed to send data" << endln;

   return res;
}

int 
ConcreteS::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  int res = 0, cnt = 0;

  static Vector data(13);

  res = theChannel.recvVector(this->getDbTag(), commitTag, data);
  if (res < 0) {
   opserr << "ConcreteS::recvSelf -- could not recv Vector" << endln;
   return res;
  }

  this->setTag(int(data(cnt++)));
  E = data(cnt++);
  nu = data(cnt++);
  fc = data(cnt++);
  ft = data(cnt++);
  Es = data(cnt++);
  cStrain0 = data(cnt++);

  setInitials();

  int i;
  for (i = 0; i < 3; i++)
    strain0(i) = data(cnt++);

  for (i = 0; i < 3; i++)
    stress0(i) = data(cnt++);

  return res;
}
 
