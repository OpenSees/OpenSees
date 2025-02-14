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
// $Date: 2021-07-02 14:29:01 $
// $Source: /usr/local/cvs/OpenSees/SRC/damping/URDDamping.h,v $

// Revised by: Y Tian
// Created: 02/2020
// Revision: A
//
// Description: This file contains the definition for the URDDamping class.
// URDDamping provides the abstraction of an elemental damping imposition
// providing user-define damping over a frequency range
//
// Reference:
//

// What: "@(#) URDDamping.cpp, revA"

#include <math.h>
#include <Vector.h>
#include <Matrix.h>
#include <Domain.h>
#include <Channel.h>
#include <elementAPI.h>
#include <URDDamping.h>
#include <FEM_ObjectBroker.h>
#include <ID.h>
//#include <algorithm>
#include <iostream>
#include <fstream>
#include <iomanip>
using namespace std;

//extern StaticAnalysis *theStaticAnalysis;

// constructor:
URDDamping::URDDamping(int tag, int nfreq, Matrix *etaf, double tol, double t1, double t2, TimeSeries *f, int ptag, int iter):
Damping(tag, DMP_TAG_URDDamping),
nComp(0), nFilter(0),
numfreq(nfreq), dptol(tol), ta(t1), td(t2), fac(f), prttag(ptag), maxiter(iter),
alpha(0), omegac(0), omegaetaf(0), qL(0), qLC(0), qd(0), qdC(0), q0(0), q0C(0),
Freqlog(0), Fredif(0), Freqk(0), Freqb(0)
{
  etaFreq = new Matrix(*etaf);
  Initialize();
}


// constructor:
URDDamping::URDDamping(int tag, int nfreq, Matrix* etaf, double tol, double t1, double t2, TimeSeries *f, int nF, Vector *a, Vector *w, Vector *ef, int ptag, int iter):
Damping(tag, DMP_TAG_URDDamping),
nComp(0), nFilter(0),
numfreq(nfreq), dptol(tol), ta(t1), td(t2), fac(f), prttag(ptag), maxiter(iter),
alpha(0), omegac(0), omegaetaf(0), qL(0), qLC(0), qd(0), qdC(0), q0(0), q0C(0),
Freqlog(0), Fredif(0), Freqk(0), Freqb(0)
{
  etaFreq = new Matrix(*etaf);
  if (nF > 0 && a->Size() == nF && w->Size() == nF && ef->Size() == nF)
  {
    nFilter = nF;
    alpha = new Vector(*a);
    omegac = new Vector(*w);
    omegaetaf = new Vector(*ef);
  }
  else
  {
    Initialize();
  }
}


// constructor:
// invoked by a FEM_ObjectBroker, recvSelf() needs to be invoked on this object.
URDDamping::URDDamping():
Damping(0, DMP_TAG_URDDamping),
nComp(0), nFilter(0),
numfreq(0), etaFreq(0), dptol(0.0), ta(0.0), td(0.0), fac(0), prttag(0), maxiter(0),
alpha(0), omegac(0), omegaetaf(0), qL(0), qLC(0), qd(0), qdC(0), q0(0), q0C(0),
Freqlog(0), Fredif(0), Freqk(0), Freqb(0)
{

}


// destructor:
URDDamping::~URDDamping() 
{
  if (fac) delete fac;
  if (alpha) delete alpha;
  if (omegac) delete omegac;
  if (omegaetaf) delete omegaetaf;
  if (qL) delete qL;
  if (qLC) delete qLC;
  if (qd) delete qd;
  if (qdC) delete qdC;
  if (q0) delete q0;
  if (q0C) delete q0C;
  if (Freqlog) delete Freqlog;
  if (Fredif) delete Fredif;
  if (Freqk) delete Freqk;
  if (Freqb) delete Freqb;
  if (etaFreq) delete etaFreq;
}

int
URDDamping::Initialize(void)
{
   // Update order
  double tmpf, tempeta;
  for (int i = 0; i < numfreq-1; i++)
  {
      for (int j = i+1; j < numfreq; j++)
      {
          if ((*etaFreq)(i, 0) > (*etaFreq)(j, 0)) {
              tmpf = (*etaFreq)(j, 0);
              tempeta = (*etaFreq)(j, 1);
              (*etaFreq)(j, 0) = (*etaFreq)(i, 0);
              (*etaFreq)(j, 1) = (*etaFreq)(i, 1);
              (*etaFreq)(i, 0) = tmpf;
              (*etaFreq)(i, 1) = tempeta;
          }
          else if((*etaFreq)(i, 0) == (*etaFreq)(j, 0)) {
              opserr << "WARNING invalid frequency series - freqi cannot equal freqj\n";
              return  1;
          }
      }
  }

  //Calculate interpolation parameters
  Freqlog = new Vector(numfreq);
  Fredif = new Vector(numfreq - 1);
  Freqk = new Vector(numfreq - 1);
  Freqb = new Vector(numfreq - 1);
  for (int i = 0; i < numfreq; i++)
  {
     (*Freqlog)(i) = log10((*etaFreq)(i, 0));
  }
  for (int i = 0; i < numfreq - 1; i++)
  {
      (*Fredif)(i) = (*Freqlog)(i + 1) - (*Freqlog)(0);
      (*Freqk)(i) = ((*etaFreq)(i + 1, 1) - (*etaFreq)(i, 1)) / ((*Freqlog)(i + 1) - (*Freqlog)(i));
      (*Freqb)(i) = (*etaFreq)(i, 1) - (*Freqk)(i) * (*Freqlog)(i);
  }
  
  int countfreq1, countfreq2;
  int omigacount1, omigacount2;
  int tagFindFilter = 0;
  int iter = 0;
  double delta1 = dptol; //0.05; //relative damping ratio error threshold
  double delta2 = 0.0002; //absolute damping ratio error threshold
  double f1log = (*Freqlog)(0);
  double f2log = (*Freqlog)(numfreq-1);
  nFilter = 1;
  //Calculate alpha
  while (iter < maxiter) 
  {
      ++iter;
      if ((tagFindFilter == 0)) {
          ++nFilter;
      }
      else if ((tagFindFilter == 1) || (tagFindFilter == 2)) {
          nFilter = nFilter + numfreq - 1;
      }
      
      omegac = new Vector(nFilter);
      alpha = new Vector(nFilter);
      omegaetaf = new Vector(nFilter);
      
      //Calculate omegac
      if (tagFindFilter == 0) {
          //Method #0
          double dfreq = (f2log - f1log) / (nFilter - 1);
          Vector *omegactag = new Vector(nFilter);
          omigacount2 = 0;
          for (int i = 0; i < numfreq - 1; i++)
          {
              omigacount1 = floor((*Fredif)(i) / dfreq) + 1;
              if (omigacount1 > omigacount2) {
                  for (int j = omigacount2; j < omigacount1; j++)
                  {
                      (*omegactag)(j) = i;
                  }
                  omigacount2 = omigacount1;
              }
          }
          (*omegactag)(0) = 0;
          (*omegactag)(nFilter - 1) = numfreq - 2;
          tmpf = f1log - dfreq;
          for (int i = 0; i < nFilter; i++)
          {
              tmpf += dfreq;
              (*omegaetaf)(i) = (*Freqk)((*omegactag)(i)) * tmpf + (*Freqb)((*omegactag)(i));
          }

          for (int i = 0; i < nFilter; ++i)
          {
              (*omegac)(i) = 6.28318530718 * pow(10.0, f1log + i * dfreq);
          }
          delete omegactag;
      }
      else if ((tagFindFilter == 1) || (tagFindFilter == 2)) {
          //Method #1
          for (int i = 0; i < numfreq - 1; i++) {
              for (int j = 0; j < iter; j++) {
                  (*omegac)(i * iter + j) = (*Freqlog)(i) + ((*Freqlog)(i + 1) - (*Freqlog)(i)) / iter * j;
                  (*omegaetaf)(i * iter + j) = (*Freqk)(i) * (*omegac)(i * iter + j) + (*Freqb)(i);
                  (*omegac)(i * iter + j) = 6.28318530718 * pow(10.0, (*omegac)(i * iter + j));
              }
          }
          (*omegac)(nFilter - 1) = (*Freqlog)(numfreq - 1);
          (*omegac)(nFilter - 1) = 6.28318530718 * pow(10.0, (*omegac)(nFilter - 1));
          (*omegaetaf)(nFilter - 1) = (*etaFreq)(numfreq-1, 1);
      }
      

      Vector* y = new Vector(nFilter);
      Matrix* X = new Matrix(nFilter, nFilter);
      //Calculate y and X
      if ((tagFindFilter == 0) || (tagFindFilter == 1))
      {
          //Prepare interpolated damping ratio for integral
          double df = 0.005;
          int nf = ceil((f2log - f1log) / df) + 1;
          df = (f2log - f1log) / (nf - 1);
          Vector* tmptag = new Vector(nf);
          Vector* tmpeta = new Vector(nf);
          countfreq2 = 0;

          for (int i = 0; i < numfreq - 1; i++)
          {
              countfreq1 = floor((*Fredif)(i) / df) + 1;
              if (countfreq1 > countfreq2) {
                  for (int j = countfreq2; j < countfreq1; j++)
                  {
                      (*tmptag)(j) = i;
                  }
                  countfreq2 = countfreq1;
              }
          }
          (*tmptag)(0) = 0;
          (*tmptag)(nf - 1) = numfreq - 2;

          tmpf = f1log - df;
          for (int i = 0; i < nf; i++)
          {
              tmpf += df;
              (*tmpeta)(i) = (*Freqk)((*tmptag)(i)) * tmpf + (*Freqb)((*tmptag)(i));
          }
          delete tmptag;

          for (int i = 0; i < nf; ++i)
          {
              double omega = 6.28318530718 * pow(10.0, f1log + i * df);
              for (int j = 0; j < nFilter; ++j)
              {
                  double wjn = omega / (*omegac)[j];
                  double phij = 2.0 * wjn / (1.0 + wjn * wjn);
                  (*y)(j) += (phij * (*tmpeta)(i) * (*omegaetaf)(j));
                  for (int k = 0; k < nFilter; ++k)
                  {
                      double wkn = omega / (*omegac)[k];
                      double phik = 2.0 * wkn / (1.0 + wkn * wkn);
                      (*X)(j, k) += (phij * phik * (*omegaetaf)(j) * (*omegaetaf)(k));
                  }
              }
          }
          *alpha = (*y) / (*X);

          delete y;
          delete X;
          bool converged = true;

          for (int i = 0; i < nf; ++i)
          {
              double omega = 6.28318530718 * pow(10.0, f1log + i * df);
              double err = 0.0;
              for (int j = 0; j < nFilter; ++j)
              {
                  double wjn = omega / (*omegac)[j];
                  double phij = 2.0 * wjn / (1.0 + wjn * wjn);
                  err += (*alpha)(j) * phij * (*omegaetaf)(j);
              }
              err -= (*tmpeta)(i);
              if ((*tmpeta)(i) > 0.005) {
                  if (err / (*tmpeta)(i) > delta1 || err / (*tmpeta)(i) < -delta1)
                  {
                      converged = false;
                      break;
                  }
              }
              else {
                  if (err > delta2 || err < -delta2)
                  {
                      converged = false;
                      break;
                  }
              }
          }
          delete tmpeta;

          //Check out-zone
          for (int i = 1; i < 8; ++i)
          {
              double omega = 6.28318530718 * pow(10.0, f1log - i * log10(2.0));
              double err = 0.0;
              for (int j = 0; j < nFilter; ++j)
              {
                  double wjn = omega / (*omegac)[j];
                  double phij = 2.0 * wjn / (1.0 + wjn * wjn);
                  err += (*alpha)(j) * phij * (*omegaetaf)(j);
              }
              if (err < 0.0) {
                  converged = false;
                  break;
              }
              else if (err > (*omegaetaf)(0)) {
                  converged = false;
                  break;
              }
          }
          for (int i = 1; i < 8; ++i)
          {
              double omega = 6.28318530718 * pow(10.0, f2log + i * log10(2.0));
              double err = 0.0;
              for (int j = 0; j < nFilter; ++j)
              {
                  double wjn = omega / (*omegac)[j];
                  double phij = 2.0 * wjn / (1.0 + wjn * wjn);
                  err += (*alpha)(j) * phij * (*omegaetaf)(j);
              }
              if (err < 0.0) {
                  converged = false;
                  break;
              }
              else if (err > (*omegaetaf)(nFilter - 1)) {
                  converged = false;
                  break;
              }
          }

          if (converged) {
              if (prttag == 1) 
              {
                  opserr << "damping URD: nFilter " << nFilter << "\n";
                  opserr << "damping URD: alpha: " << "\n";
                  for (int j = 0; j < nFilter; ++j) {
                      cout << setiosflags(ios::fixed) << setprecision(16) << (*alpha)(j) << " ";
                  }
                  opserr << "\n";
                  opserr << "damping URD: nf " << nf << "\n";
                  opserr << "damping URD: omega_cn: " << "\n";
                  for (int j = 0; j < nFilter; ++j) {
                      cout << setiosflags(ios::fixed) << setprecision(16) << (*omegac)(j) << " ";
                  }
                  opserr << "\n";
                  opserr << "damping URD: ita_omega_cn: " << "\n";
                  for (int j = 0; j < nFilter; ++j) {
                      cout << setiosflags(ios::fixed) << setprecision(16) << (*omegaetaf)(j) << " ";
                  }
                  opserr << "\n";
              }
              
              //system("pause");
              break;
          }
          else if (iter == maxiter) {
              if (tagFindFilter == 0) {
                  tagFindFilter = 2; //Give up method 1
                  nFilter = 1;
                  iter = 0;
                  opserr << "Change nFilter calculation method to Method #"<< tagFindFilter <<"\n";
              }
              else if (tagFindFilter == 1) {
                  opserr << "Error: Damping model not converged\n";
                  for (int i = 0; i < numfreq; i++) {
                      (*etaFreq)(i, 1) /= (2.0);
                  }
                  opserr << (*etaFreq) << "\n";
                  return  1;
                  break;
                  //tagFindFilter = 2;
                  //nFilter = 1;
                  //iter = 0;
                  //opserr << "Change nFilter calculation method to Method #"<< tagFindFilter <<"\n";
              }
          }
      }
      else if (tagFindFilter == 2)
      {
          for (int i = 0; i < nFilter; ++i)
          {
              (*y)(i) = (*omegaetaf)(i);
              for (int j = 0; j < nFilter; ++j)
              {
                  double wjn = (*omegac)[i] / (*omegac)[j];
                  double phij = 2.0 * wjn / (1.0 + wjn * wjn);
                  (*X)(i, j) = (*omegaetaf)(j) * phij;
              }
          }
          *alpha = (*y) / (*X);
          delete y;
          delete X;    
          bool converged = true;
          //check out-zone
          for (int i = 1; i < 8; ++i)
          {
              double omega = 6.28318530718 * pow(10.0, f1log - i * log10(2.0));
              double err = 0.0;
              for (int j = 0; j < nFilter; ++j)
              {
                  double wjn = omega / (*omegac)[j];
                  double phij = 2.0 * wjn / (1.0 + wjn * wjn);
                  err += (*alpha)(j) * phij * (*omegaetaf)(j);
              }
              if (err < 0.0) {
                  converged = false;
                  break;
              }
              else if (err > (*omegaetaf)(0)) {
                  converged = false;
                  break;
              }
          }
          for (int i = 1; i < 8; ++i)
          {
              double omega = 6.28318530718 * pow(10.0, f2log + i * log10(2.0));
              double err = 0.0;
              for (int j = 0; j < nFilter; ++j)
              {
                  double wjn = omega / (*omegac)[j];
                  double phij = 2.0 * wjn / (1.0 + wjn * wjn);
                  err += (*alpha)(j) * phij * (*omegaetaf)(j);
              }
              if (err < 0.0) {
                  converged = false;
                  break;
              }
              else if (err > (*omegaetaf)(nFilter - 1)) {
                  converged = false;
                  break;
              }
          }
          if (converged) {
              if (prttag == 1) 
              {
                 opserr << "damping URD: nFilter " << nFilter << "\n";
                 opserr << "damping URD: alpha: " << "\n";
                 for (int j = 0; j < nFilter; ++j) {
                     cout << setiosflags(ios::fixed) << setprecision(16) << (*alpha)(j) << " ";
                 }
                 opserr << "\n";
                 opserr << "damping URD: omega_cn: " << "\n";
                 for (int j = 0; j < nFilter; ++j) {
                     cout << setiosflags(ios::fixed) << setprecision(16) << (*omegac)(j) << " ";
                 }
                 opserr << "\n";
                 opserr << "damping URD: ita_omega_cn: " << "\n";
                 for (int j = 0; j < nFilter; ++j) {
                     cout << setiosflags(ios::fixed) << setprecision(16) << (*omegaetaf)(j) << " ";
                 } 
                 opserr << "\n";
              }
             
              //system("pause");
              break;
          }
          else if (iter == maxiter) 
          {              
              tagFindFilter = 1;
              nFilter = 1;
              iter = 0;
              opserr << "Change nFilter calculation method to Method #" << tagFindFilter << "\n";
              
          }
      }
      
  }
  
  return 0;
}

int
URDDamping::commitState(void)
{
  *qdC = *qd;
  *q0C = *q0;
  *qLC = *qL;
  return 0;
}


int
URDDamping::revertToLastCommit(void)
{
  *qd = *qdC;
  *q0 = *q0C;
  *qL = *qLC;
  return 0;
}


int
URDDamping::revertToStart(void)
{
  (*qd).Zero();
  (*qdC).Zero();
  (*q0).Zero();
  (*q0C).Zero();
  (*qL).Zero();
  (*qLC).Zero();
  return 0;
}


int 
URDDamping::setDomain(Domain *domain, int nC)
{
  theDomain = domain;
  nComp = nC;
  
  qd = new Vector(nComp);
  qdC = new Vector(nComp);
  q0 = new Vector(nComp);
  q0C = new Vector(nComp);
  qL = new Matrix(nComp, nFilter);
  qLC = new Matrix(nComp, nFilter);
  
  return 0;
}


int
URDDamping::update(Vector q)
{
  double t = theDomain->getCurrentTime();
  double dT = theDomain->getDT();
  StaticAnalysis **theStaticAnalysis = OPS_GetStaticAnalysis();
  if (*theStaticAnalysis)
  {
    *q0 = q;
    (*qd).Zero();
    for (int i = 0; i < nFilter; ++i)
      for (int j = 0; j < nComp; ++j)
        (*qL)(j,i) = q(j);
  }
  else if (dT > 0.0)
  {
    *q0 = q;
    (*qd).Zero();
    if (t < td)
    {
      if (t > ta)
      {
        for (int i = 0; i < nFilter; ++i)
        {
          double dTomegac = dT * (*omegac)(i);
          double cd = 4.0 * (*alpha)(i) * (*omegaetaf)(i) / (2.0 + dTomegac);
          double c0 = dTomegac / (2.0 + dTomegac);
          double cL = (2.0 - dTomegac) / (2.0 + dTomegac);
          for (int j = 0; j < nComp; ++j)
          {
            (*qd)(j) += cd * ((*q0C)(j) + (*q0)(j) - 2.0 * (*qLC)(j,i));
            (*qL)(j,i) = c0 * ((*q0C)(j) + (*q0)(j)) + cL * (*qLC)(j,i);
          }
        }
        *qd -= *qdC;
      }
      else
      {
        for (int i = 0; i < nFilter; ++i)
          for (int j = 0; j < nComp; ++j)
            (*qL)(j,i) = q(j);
      }
      if (fac) *qd *= fac->getFactor(t);
    }
  }
  return 0;
}

const Vector &
URDDamping::getDampingForce(void)
{
  return (*qd);
}

double URDDamping::getStiffnessMultiplier(void)
{
  double t = theDomain->getCurrentTime();
  double dT = theDomain->getDT();
  double km = 0.0;
  StaticAnalysis **theStaticAnalysis = OPS_GetStaticAnalysis();
  if (!*theStaticAnalysis && dT > 0.0 && t > ta && t < td)
  {
    for (int i = 0; i < nFilter; ++i)
    {
        km += 4.0 * (*alpha)(i) * (*omegaetaf)(i) / (2.0 + (*omegac)(i) * dT); 
    }
    if (fac) km *= fac->getFactor(t);
  }
  return 1.0 + km;
}

Damping *URDDamping::getCopy(void)
{
  // create a new instance of URDDamping 

  URDDamping *theCopy;

  theCopy = new URDDamping(this->getTag(), numfreq, etaFreq, dptol, ta, td, fac, nFilter, alpha, omegac, omegaetaf, prttag, maxiter);

  return theCopy;
}


int 
URDDamping::sendSelf(int cTag, Channel &theChannel)
{
  int dbTag = this->getDbTag();
  static ID idData(2);
  static Vector data(7);
  static Matrix *data2;

  if (fac)
  {
    idData(0) = fac->getClassTag();
    int seriesDbTag = fac->getDbTag();
    if (seriesDbTag == 0)
    {
      seriesDbTag = theChannel.getDbTag();
      fac->setDbTag(seriesDbTag);
    }
    idData(1) = seriesDbTag;
  }
  else
    idData(0) = -1;

  

  data(0) = this->getTag();
  data(1) = numfreq;
  data(2) = dptol;
  data(3) = ta;
  data(4) = td;
  data(5) = prttag;
  data(6) = maxiter;
  (*data2) = (*etaFreq);

  int res = theChannel.sendID(dbTag, cTag, idData);
  res += theChannel.sendVector(dbTag, cTag, data);
  res += theChannel.sendMatrix(dbTag, cTag, *data2);
  if (res < 0)
  {
    opserr << " URDDamping::sendSelf() - data could not be sent\n" ;
    return -1;
  }

  if (fac)
  {
    res = fac->sendSelf(cTag, theChannel);
    if (res < 0)
    {
      opserr << " URDDamping::sendSelf() - failed to send factor series\n";
      return res;
    }
  }

  return 0;
}


int 
URDDamping::recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  int dbTag = this->getDbTag();
  static ID idData(2);
  static Vector data(7);
  static Matrix *data2;

  int res = theChannel.recvID(dbTag, cTag, idData);
  res += theChannel.recvVector(dbTag, cTag, data);
  res += theChannel.recvMatrix(dbTag, cTag, *data2);
  if (res < 0) {opserr << " URDDamping::recvSelf() - data could not be received\n" ;
    return -1;
  }

  int seriesClassTag = idData(0);
  if (seriesClassTag != -1)
  {
    int seriesDbTag = idData(1);
    if (fac == 0 || fac->getClassTag() != seriesClassTag)
    {
      if (fac != 0)
        delete fac;
      fac = theBroker.getNewTimeSeries(seriesClassTag);
      if (fac == 0)
      {
        opserr << "GroundMotion::recvSelf - could not create a Series object\n";
        return -2;
      }
    }
    fac->setDbTag(seriesDbTag);
    res = fac->recvSelf(cTag, theChannel, theBroker);
    if (res < 0)
    {
      opserr << "URDDamping::recvSelf() - factor series could not be received\n";
      return res;
    }
  }
    
  this->setTag((int)data(0));
  numfreq = data(1);
  dptol = data(2);
  ta = data(3);
  td = data(4);
  prttag = data(5);
  maxiter = data(6);
  (*etaFreq) = (*data2);

  Initialize();
  return 0;
}

void
URDDamping::Print(OPS_Stream &s, int flag)
{
  if (flag == OPS_PRINT_CURRENTSTATE)
  {
    s << "\nDamping: " << this->getTag() << " Type: URDDamping";
    s << "\tnumber of frequencies: " << numfreq << endln;
    s << ", \"frequency\": [";
    for (int i = 0; i < numfreq; ++i)
    {
        s << (*etaFreq)(i, 0);
    }
    s << "]";
    s << ", \"loss factor\": [";
    for (int i = 0; i < numfreq; ++i)
    {
        s << (*etaFreq)(i, 1);
    }
    s << "]";
    s << "\ttol: " << dptol << endln;
    s << "\tactivation time: " << ta << endln;
    s << "\tdeactivation time: " << td << endln;
  }

  if (flag == OPS_PRINT_PRINTMODEL_JSON)
  {
    s << "\t\t\t{\"name\": \"" << this->getTag() << "\", \"type\": \"URDDamping\"";
    s << ", \"number of frequencies\": [" << numfreq << "]";
    s << ", \"frequency\": [";
    for (int i = 0; i < numfreq; ++i)
    {
        s << (*etaFreq)(i,0) ;
    } 
    s << "]";
    s << ", \"loss factor\": [";
    for (int i = 0; i < numfreq; ++i)
    {
        s << (*etaFreq)(i, 1);
    }
    s << "]";
    s << ", \"tol\": [" << dptol << "]";
    s << ", \"activation time\": [" << ta << "]";
    s << ", \"deactivation time\": [" << td << "]";
    s << "}";
  }
}
