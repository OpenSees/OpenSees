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
// $Date: 2012-06-08 22:03:16 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/PlaneStressUserMaterial.cpp,v $

//
// Yuli Huang (yulihuang@gmail.com) & Xinzheng Lu (luxz@tsinghua.edu.cn)
//
// Plane Stress User Defined Material
//

#include <PlaneStressUserMaterial.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <MaterialResponse.h>
#include <elementAPI.h>

#ifdef _HAVE_PSUMAT

#ifdef _WIN32

#else
#define PSUMAT psumat_
#endif

extern "C" 
{
  void PSUMAT(int *nstatev, int *nprops, double *props,
              double *stress, double *strain0, double *strain1, double *dstrain,
              double *statev, double *tangent);
}
#else
void PSUMAT(int *nstatev, int *nprops, double *props,
	    double *stress, double *strain0, double *strain1, double *dstrain,
	    double *statev, double *tangent)
{
  opserr << "PSUMAT - NOT DEFINED IN THIS VERSION, SOURCE CODE RESTRICTED\n";
}
#endif

void* OPS_PlaneStressUserMaterial()
{
    int argc = OPS_GetNumRemainingInputArgs() + 2;
    if (argc < 6) {
	opserr << "WARNING: Insufficient arguements\n";
	opserr << "Want: nDMaterial PlaneStressUserMaterial tag? nstatevs? nprops? prop1? ... propn?" << endln;
	return 0;
    }

    // int tag, nstatevs, nprops;
    int idata[3];
    int numdata = 3;
    if (OPS_GetIntInput(&numdata,idata) < 0) {
	opserr << "WARNING invalid nDMaterial PlaneStressUserMaterial int inputs" << endln;
	return 0;
    }
    int tag = idata[0];
    int nstatevs = idata[1];
    int nprops = idata[2];
    if (nstatevs < 1) nstatevs = 1;
    if (nprops < 1) nprops = 1;

    if (OPS_GetNumRemainingInputArgs() < nprops) {
	opserr << "WARNING insufficient arguments\n";
	return 0;
    }
    double *props;
    props = new double[nprops];
    if (OPS_GetDoubleInput(&nprops , props) < 0) {
	opserr << "WARNING invalid prop" << endln;
	opserr << "PlaneStressUserMaterial: " << tag << endln;
	return 0;
    }
    void* temp = new PlaneStressUserMaterial( tag, nstatevs, nprops, props);
    if (props != 0) delete props;

    return temp;
}


//null constructor
PlaneStressUserMaterial::PlaneStressUserMaterial( ) : 
NDMaterial(0, ND_TAG_PlaneStressUserMaterial ), 
strain0(3), strain(3), stress0(3), stress(3),
tangent(3,3), eTangent(3,3), 
vprops(0), statev0(0), statev(0),
props(0), statevdata(0),
nstatevs(0), nprops(0)
{ 

}


//full constructor
PlaneStressUserMaterial::PlaneStressUserMaterial(int tag, int istatevs, int iprops, double *rprops) :
NDMaterial( tag, ND_TAG_PlaneStressUserMaterial ),
strain0(3), strain(3), stress0(3), stress(3),
tangent(3,3), eTangent(3,3),
statev0(0), statev(0),
statevdata(0),
nstatevs(istatevs), nprops(iprops)
{
  props = new double[nprops];
  for (int i = 0; i < nprops; i++)
  {
    props[i] = rprops[i];
  }
  vprops = new Vector(props, nprops);

  //  vprops->setData(props, nprops);
  for (int i=0; i<9; i++)
    tangentdata[i]=0;


  statevdata = new double[nstatevs];
  statev0 = new Vector(istatevs);
  statev  = new Vector(istatevs);
  setInitials();
}


//destructor
PlaneStressUserMaterial::~PlaneStressUserMaterial( ) 
{ 
  if (props != 0)
    delete [] props;
  if (vprops != 0)
    delete vprops;

  if (statevdata != 0)
    delete [] statevdata;

  if (statev0 != 0)
    delete statev0;
  if (statev != 0)
    delete statev;
} 

void PlaneStressUserMaterial::setInitials()
{
  int i = 0;
  for (i = 0; i < 3; i++)
  {
    stressdata[i]  = 0.0;
    strain0data[i] = 0.0;
    straindata[i]  = 0.0;
    dstraindata[i] = 0.0;
  }
  for (i = 0; i < nstatevs; i++) statevdata[i] = 0.0;

  PSUMAT(&nstatevs, &nprops, props,
         stressdata, strain0data, straindata, dstraindata,
         statevdata, tangentdata);
  
  for (i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) {
      eTangent(i,j) = tangentdata[j + 3 * i];
    }
  
  tangent = eTangent;
}

//make a clone of this material
NDMaterial*
PlaneStressUserMaterial::getCopy( ) 
{
  PlaneStressUserMaterial *clone ;   //new instance of this class
  clone = new PlaneStressUserMaterial(this->getTag(), nstatevs, nprops, props);
  return clone ;
}


//make a clone of this material
NDMaterial* 
PlaneStressUserMaterial::getCopy( const char *type ) 
{
  if ((strcmp(type, "PlaneStress") == 0) ||
      (strcmp(type, "PlaneStress2D") == 0))
    return this->getCopy( ) ;
  else
    return 0;
}


//send back order of strain in vector form
int 
PlaneStressUserMaterial::getOrder( ) const
{
  return 3 ;
}


const char*
PlaneStressUserMaterial::getType( ) const 
{
  return "PlaneStress" ; 
}



//swap history variables
int 
PlaneStressUserMaterial::commitState( ) 
{
  stress0 = stress;
  strain0 = strain;
  (*statev0) = (*statev);
  return 0;
}


//revert to last saved state
int 
PlaneStressUserMaterial::revertToLastCommit( )
{
  return 0;
}


//revert to start
int
PlaneStressUserMaterial::revertToStart( )
{
  strain0.Zero();
  strain.Zero();
  stress0.Zero();
  stress.Zero();
  statev0->Zero();
  statev->Zero();
  return 0;
}

//receive the strain
int 
PlaneStressUserMaterial::setTrialStrain( const Vector &strainFromElement )
{
  strain(0) = strainFromElement(0) ;
  strain(1) = strainFromElement(1) ;
  strain(2) = strainFromElement(2) ;

  tangent = eTangent;

  int i = 0;
  for (i = 0; i < 3; i++)
  {
    stressdata[i]  = stress0(i);
    strain0data[i] = strain0(i);
    straindata[i]  = strain(i);
    dstraindata[i] = strain(i) - strain0(i);
  }
  for (i = 0; i < nstatevs; i++) statevdata[i] = (*statev0)(i);

  PSUMAT(&nstatevs, &nprops, props,
         stressdata, strain0data, straindata, dstraindata,
         statevdata, tangentdata);

  stress.setData(stressdata, 3);
  statev->setData(statevdata, nstatevs);
  for (i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      tangent(i,j) = tangentdata[j + 3 * i];

  return 0;
}


//send back the strain
const Vector& 
PlaneStressUserMaterial::getStrain( )
{
  return strain ;
}


//send back the stress 
const Vector&  
PlaneStressUserMaterial::getStress( )
{
  return stress ;
}


//send back the tangent 
const Matrix&  
PlaneStressUserMaterial::getTangent( )
{
  return tangent ;
}

const Matrix&  
PlaneStressUserMaterial::getInitialTangent
( )
{
  return eTangent ;
}


//print out data
void  
PlaneStressUserMaterial::Print( OPS_Stream &s, int flag )
{
    if (flag == OPS_PRINT_PRINTMODEL_MATERIAL) {
        s << "PlaneStressUserMaterial tag: " << this->getTag() << endln;
        for (int i = 0; i < nprops; i++)
            s << "prop" << i << ":  " << props[i] << " ";
    }
    
    if (flag == OPS_PRINT_PRINTMODEL_JSON) {
        s << "\t\t\t{";
        s << "\"name\": \"" << this->getTag() << "\", ";
        s << "\"type\": \"PlaneStressUserMaterial\", ";
        s << "\"properties\": [";
        for (int i = 0; i < nprops-1; i++)
            s << props[i] << ", ";
        s << props[nprops-1] << "]}";
    }
}


int 
PlaneStressUserMaterial::sendSelf(int commitTag, Channel &theChannel) 
{
  int res = 0;

  int dataTag = this->getDbTag();

  static ID idData(3);

  idData(0) = this->getTag();
  idData(1) = nstatevs;
  idData(2) = nprops;

  res = theChannel.sendID(dataTag, commitTag, idData);
  if (res < 0) 
    opserr << "PlaneStressUserMaterial::sendSelf() - failed to send ID" << endln;

  res = theChannel.sendVector(dataTag, commitTag, strain0);
  if (res < 0) 
    opserr << "PlaneStressUserMaterial::sendSelf() - failed to send data" << endln;

  res = theChannel.sendVector(dataTag, commitTag, stress0);
  if (res < 0) 
    opserr << "PlaneStressUserMaterial::sendSelf() - failed to send data" << endln;

  res = theChannel.sendVector(dataTag, commitTag, *statev0);
  if (res < 0) 
    opserr << "PlaneStressUserMaterial::sendSelf() - failed to send data" << endln;

  res = theChannel.sendVector(dataTag, commitTag, *vprops);
  if (res < 0) 
    opserr << "PlaneStressUserMaterial::sendSelf() - failed to send data" << endln;

  return res;
}

int 
PlaneStressUserMaterial::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  int res = 0;

  int dataTag = this->getDbTag();

  static ID idData(3);

  res = theChannel.recvID(dataTag, commitTag, idData);
  if (res < 0) {
   opserr << "PlaneStressUserMaterial::recvSelf -- could not recv ID" << endln;
   return res;
  }

  this->setTag(idData(0));
  
  if (nstatevs != idData(1))
  {
    nstatevs = idData(1);
    if (statev0 != 0) delete statev0;
    statev0 = new Vector(nstatevs);
    if (statev != 0) delete statev;
    statev = new Vector(nstatevs);
    if (statevdata != 0) delete statevdata;
    statevdata = new double[nstatevs];
  }
  
  if (nprops != idData(2))
  {
    nprops = idData(2);
    if (vprops != 0) delete vprops;
    vprops = new Vector(nprops);
    if (props != 0) delete props;
    props = new double[nprops];
  }

  res = theChannel.recvVector(dataTag, commitTag, strain0);
  if (res < 0) {
   opserr << "PlaneStressUserMaterial::recvSelf -- could not recv data" << endln;
   return res;
  }

  res = theChannel.recvVector(dataTag, commitTag, stress0);
  if (res < 0) {
   opserr << "PlaneStressUserMaterial::recvSelf -- could not recv data" << endln;
   return res;
  }

  res = theChannel.recvVector(dataTag, commitTag, *statev0);
  if (res < 0) {
   opserr << "PlaneStressUserMaterial::recvSelf -- could not recv data" << endln;
   return res;
  }

  res = theChannel.recvVector(dataTag, commitTag, *vprops);
  if (res < 0) {
   opserr << "PlaneStressUserMaterial::recvSelf -- could not recv data" << endln;
   return res;
  }

  setInitials();

  return res;
}
 
