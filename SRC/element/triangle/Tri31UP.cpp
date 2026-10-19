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

// ============================================================================
// 2026 By Pedro Arduino @ University of Washington, WA.
// parduino@uw.edu
// ============================================================================
#include "Tri31UP.h"
#include <elementAPI.h>
#include <Domain.h>
#include <Node.h>
#include <NDMaterial.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <ElementResponse.h>
#include <ElementalLoad.h>
#include <Information.h>
#include <Parameter.h>
#include <Renderer.h>
#include <cmath>
#include <cstring>
#include <cfloat>

void *OPS_Tri31UP()
{
    static bool announced=false;
    if (!announced) {
        announced=true;
        opserr << "Tri31UP element - Written: P.Arduino, U.Washington\n";
    }
    if (OPS_GetNDM()!=2 || OPS_GetNDF()!=3 || OPS_GetNumRemainingInputArgs()<12) {
        opserr << "Tri31UP requires ndm=2, ndf=3: element Tri31UP tag n1 n2 n3 matTag "
               << "thickness fBulk fDensity kx ky voidRatio alpha <bx by> <-lumped>\n";
        return 0;
    }
    int ids[5], n=5;
    double d[9]={0};
    if (OPS_GetIntInput(&n,ids)<0) return 0;
    n=7;
    if (OPS_GetDoubleInput(&n,d)<0) return 0;
    bool lumped=false;
    int nb=0;
    while (OPS_GetNumRemainingInputArgs()>0) {
        const int before=OPS_GetNumRemainingInputArgs();
        double value;
        n=1;
        if (OPS_GetDoubleInput(&n,&value)==0) {
            if (nb==2) { opserr << "Tri31UP: too many body components\n"; return 0; }
            d[7+nb++]=value;
        } else {
            if (OPS_GetNumRemainingInputArgs()<before) OPS_ResetCurrentInputArg(-1);
            const char *option=OPS_GetString();
            if (std::strcmp(option,"-lumped")!=0 || lumped) {
                opserr << "Tri31UP: unknown or repeated option " << option << endln;
                return 0;
            }
            lumped=true;
        }
    }
    if (nb!=0 && nb!=2) { opserr << "Tri31UP: supply both body components\n"; return 0; }
    for (int i=0;i<9;++i) if (!std::isfinite(d[i])) {
        opserr << "Tri31UP: physical inputs must be finite\n"; return 0;
    }
    if (d[0]<=0 || d[1]<=0 || d[2]<0 || d[3]<0 || d[4]<0 || d[5]<0 || d[6]<0) {
        opserr << "Tri31UP: thickness and bulk must be positive; other physical inputs nonnegative\n";
        return 0;
    }
    for (int i=0;i<3;++i) for (int j=0;j<i;++j) if (ids[i+1]==ids[j+1]) {
        opserr << "Tri31UP: repeated node\n"; return 0;
    }
    NDMaterial *mat=OPS_getNDMaterial(ids[4]);
    if (!mat) { opserr << "Tri31UP: material not found\n"; return 0; }
    NDMaterial *copy=mat->getCopy("PlaneStrain");
    if (!copy) { opserr << "Tri31UP: material must support PlaneStrain\n"; return 0; }
    const bool compatible=copy->getOrder()==3;
    delete copy;
    if (!compatible) { opserr << "Tri31UP: expected three plane-strain components\n"; return 0; }
    Tri31UP *element=new Tri31UP(ids[0],ids[1],ids[2],ids[3],*mat,
        d[0],d[1],d[2],d[3],d[4],d[5],d[6],d[7],d[8],lumped);
    // Domain::addElement does not propagate setDomain/update failure in release builds.
    element->setDomain(OPS_GetDomain());
    if (!element->getDomain()) { delete element; return 0; }
    return element;
}

Tri31UP::Tri31UP()
    : Element(0,ELE_TAG_Tri31UP), nodes(NumNodes), material(0),
      thickness(1), fBulk(1), fDensity(0), porosity(0), alpha(0), volume(0),
      lumped(false), loaded(false), valid(false), grad(NumNodes,Dim), B(Stress,SolidDofs),
      stiffness(Dofs,Dofs), initial(Dofs,Dofs), mass(Dofs,Dofs), damp(Dofs,Dofs),
      residual(Dofs), unbalance(Dofs)
{
    for (int i=0;i<NumNodes;++i) node[i]=0;
    for (int a=0;a<Dim;++a) perm[a]=body[a]=appliedBody[a]=0;
}
Tri31UP::Tri31UP(int tag, int n1, int n2, int n3, NDMaterial &mat,
    double thick, double bulk, double fluidDensity, double kx, double ky,
    double voidRatio, double stabilization, double bx, double by, bool lump)
    : Tri31UP()
{
    setTag(tag);
    nodes(0)=n1; nodes(1)=n2; nodes(2)=n3;
    material=mat.getCopy("PlaneStrain");
    thickness=thick; fBulk=bulk; fDensity=fluidDensity;
    porosity=voidRatio/(1+voidRatio); alpha=stabilization; lumped=lump;
    perm[0]=kx; perm[1]=ky; body[0]=bx; body[1]=by;
}
Tri31UP::~Tri31UP() { delete material; }

void Tri31UP::setDomain(Domain *domain)
{
    valid=false; volume=0;
    DomainComponent::setDomain(0);
    for (int i=0;i<NumNodes;++i) node[i]=0;
    if (!domain || !material || material->getOrder()!=Stress || thickness<=0) return;
    for (int i=0;i<NumNodes;++i) {
        node[i]=domain->getNode(nodes(i));
        if (!node[i] || node[i]->getNumberDOF()!=Ndf || node[i]->getCrds().Size()!=Dim) {
            opserr << "Tri31UP " << getTag() << ": requires three 2D nodes with three DOFs\n";
            return;
        }
    }
    // J=[x2-x1,x3-x1]; use translated coordinates and accept either orientation.
    double J[Dim][Dim], scale=0;
    for (int a=0;a<Dim;++a) for (int j=0;j<Dim;++j) {
        J[a][j]=node[j+1]->getCrds()(a)-node[0]->getCrds()(a);
        if (!std::isfinite(J[a][j])) return;
        if (std::fabs(J[a][j])>scale) scale=std::fabs(J[a][j]);
    }
    const double det=J[0][0]*J[1][1]-J[0][1]*J[1][0];
    if (!std::isfinite(det) || scale==0 || std::fabs(det)<=64*DBL_EPSILON*scale*scale) {
        opserr << "Tri31UP " << getTag() << ": degenerate triangle\n"; return;
    }
    volume=std::fabs(det)*thickness/2;
    if (!std::isfinite(volume)) return;
    grad(1,0)=J[1][1]/det; grad(1,1)=-J[0][1]/det;
    grad(2,0)=-J[1][0]/det; grad(2,1)=J[0][0]/det;
    for (int a=0;a<Dim;++a) grad(0,a)=-grad(1,a)-grad(2,a);
    B.Zero();
    for (int i=0;i<NumNodes;++i) {
        B(0,Dim*i)=grad(i,0); B(1,Dim*i+1)=grad(i,1);
        B(2,Dim*i)=grad(i,1); B(2,Dim*i+1)=grad(i,0);
    }
    valid=true;
    DomainComponent::setDomain(domain);
}
int Tri31UP::update()
{
    if (!valid) return -1;
    Vector u(SolidDofs), strain(Stress);
    for (int i=0;i<NumNodes;++i) for (int a=0;a<Dim;++a) u(Dim*i+a)=node[i]->getTrialDisp()(a);
    strain.addMatrixVector(0,B,u,1);
    return material->setTrialStrain(strain);
}
int Tri31UP::commitState() {
    int result=Element::commitState();
    return material ? result+material->commitState() : -1;
}
int Tri31UP::revertToLastCommit() { return material ? material->revertToLastCommit() : -1; }
int Tri31UP::revertToStart() { return material ? material->revertToStart() : -1; }

const Matrix &Tri31UP::formStiffness(const Matrix &tangent)
{
    Matrix solid(SolidDofs,SolidDofs);
    solid.addMatrixTripleProduct(0,B,tangent,volume);
    stiffness.Zero();
    for (int i=0;i<NumNodes;++i) for (int j=0;j<NumNodes;++j)
        for (int a=0;a<Dim;++a) for (int b=0;b<Dim;++b)
            stiffness(Ndf*i+a,Ndf*j+b)=solid(Dim*i+a,Dim*j+b);
    return stiffness;
}
const Matrix &Tri31UP::getTangentStiff() {
    if (!valid) { stiffness.Zero(); return stiffness; }
    return formStiffness(material->getTangent());
}
const Matrix &Tri31UP::getInitialStiff() {
    initial.Zero();
    if (valid) initial=formStiffness(material->getInitialTangent());
    return initial;
}
const Matrix &Tri31UP::getMass()
{
    mass.Zero();
    if (!valid) return mass;
    const double rho=material->getRho(); // mixture density, as in SSPquadUP
    for (int i=0;i<NumNodes;++i) for (int j=0;j<NumNodes;++j) {
        const double nij=volume*(i==j ? 2.0 : 1.0)/12;
        const double solidMass=rho*(lumped ? (i==j ? volume/3 : 0) : nij);
        for (int a=0;a<Dim;++a) mass(Ndf*i+a,Ndf*j+a)=solidMass;
        double dot=0;
        for (int a=0;a<Dim;++a) dot+=grad(i,a)*grad(j,a);
        // Never row-sum lump the pressure gradient stabilization: its row sums vanish.
        mass(Ndf*i+Dim,Ndf*j+Dim)=-porosity/fBulk*nij-alpha*volume*dot;
    }
    return mass;
}
const Matrix &Tri31UP::getDamp()
{
    damp.Zero();
    if (!valid) return damp;
    if (alphaM!=0) {
        const Matrix &M=getMass();
        for (int i=0;i<Dofs;++i) if (i%Ndf!=Dim)
            for (int j=0;j<Dofs;++j) if (j%Ndf!=Dim) damp(i,j)+=alphaM*M(i,j);
    }
    if (betaK!=0) damp.addMatrix(1,getTangentStiff(),betaK);
    if (betaK0!=0) damp.addMatrix(1,getInitialStiff(),betaK0);
    if (betaKc!=0 && Kc) damp.addMatrix(1,*Kc,betaKc);
    for (int i=0;i<NumNodes;++i) for (int j=0;j<NumNodes;++j) {
        double h=0;
        for (int a=0;a<Dim;++a) {
            damp(Ndf*i+a,Ndf*j+Dim)=-volume*grad(i,a)/3;
            damp(Ndf*j+Dim,Ndf*i+a)=-volume*grad(i,a)/3;
            h+=grad(i,a)*perm[a]*grad(j,a);
        }
        damp(Ndf*i+Dim,Ndf*j+Dim)=-volume*h;
    }
    return damp;
}
void Tri31UP::zeroLoad() {
    unbalance.Zero(); loaded=false;
    for (int a=0;a<Dim;++a) appliedBody[a]=0;
}
int Tri31UP::addLoad(ElementalLoad *load, double factor)
{
    int type;
    const Vector &data=load->getData(type,factor);
    if (type!=LOAD_TAG_SelfWeight) {
        opserr << "Tri31UP: unsupported elemental load; use SelfWeight or consistent nodal loads\n";
        return -1;
    }
    loaded=true;
    for (int a=0;a<Dim;++a) appliedBody[a]+=body[a]*factor*data(a);
    return 0;
}
int Tri31UP::addInertiaLoadToUnbalance(const Vector &accel)
{
    if (!valid) return -1;
    Vector r(Dofs);
    for (int i=0;i<NumNodes;++i) {
        const Vector &rv=node[i]->getRV(accel);
        if (rv.Size()!=Ndf) return -1;
        for (int a=0;a<Dim;++a) r(Ndf*i+a)=rv(a);
        r(Ndf*i+Dim)=0;
    }
    unbalance.addMatrixVector(1,getMass(),r,-1);
    return 0;
}
const Vector &Tri31UP::getResistingForce()
{
    residual.Zero();
    if (!valid) return residual;
    Vector force(SolidDofs);
    force.addMatrixTransposeVector(0,B,material->getStress(),volume);
    const double *g=loaded ? appliedBody : body;
    for (int i=0;i<NumNodes;++i) for (int a=0;a<Dim;++a) {
        residual(Ndf*i+a)=force(Dim*i+a)-material->getRho()*volume*g[a]/3;
        residual(Ndf*i+Dim)+=volume*fDensity*grad(i,a)*perm[a]*g[a];
    }
    residual.addVector(1,unbalance,-1);
    return residual;
}
const Vector &Tri31UP::getResistingForceIncInertia()
{
    getResistingForce();
    if (!valid) return residual;
    Vector velocity(Dofs), acceleration(Dofs);
    for (int i=0;i<NumNodes;++i) for (int a=0;a<Ndf;++a) {
        velocity(Ndf*i+a)=node[i]->getTrialVel()(a);
        acceleration(Ndf*i+a)=node[i]->getTrialAccel()(a);
    }
    residual.addMatrixVector(1,getMass(),acceleration,1);
    residual.addMatrixVector(1,getDamp(),velocity,1);
    return residual;
}

int Tri31UP::sendSelf(int commitTag, Channel &channel)
{
    if (!material) return -1;
    if (!material->getDbTag()) material->setDbTag(channel.getDbTag());
    ID id(7);
    id(0)=getTag();
    for (int i=0;i<NumNodes;++i) id(i+1)=nodes(i);
    id(4)=material->getClassTag(); id(5)=material->getDbTag(); id(6)=lumped ? 1 : 0;
    Vector data(26);
    data(0)=thickness; data(1)=fBulk; data(2)=fDensity; data(3)=porosity; data(4)=alpha;
    for (int a=0;a<Dim;++a) {
        data(5+a)=perm[a]; data(7+a)=body[a]; data(9+a)=appliedBody[a];
    }
    data(11)=loaded ? 1 : 0;
    data(12)=alphaM; data(13)=betaK; data(14)=betaK0; data(15)=betaKc;
    for (int i=0;i<Dofs;++i) data(16+i)=unbalance(i);
    data(25)=Kc ? 1 : 0;
    int result=channel.sendID(getDbTag(),commitTag,id);
    if (result<0) return result;
    result=channel.sendVector(getDbTag(),commitTag,data);
    if (result<0) return result;
    if (Kc) {
        result=channel.sendMatrix(getDbTag(),commitTag,*Kc);
        if (result<0) return result;
    }
    return material->sendSelf(commitTag,channel);
}
int Tri31UP::recvSelf(int commitTag, Channel &channel, FEM_ObjectBroker &broker)
{
    Domain *previousDomain=getDomain();
    valid=false;
    ID id(7); Vector data(26);
    int result=channel.recvID(getDbTag(),commitTag,id);
    if (result<0) return result;
    result=channel.recvVector(getDbTag(),commitTag,data);
    if (result<0) return result;
    setTag(id(0));
    for (int i=0;i<NumNodes;++i) { nodes(i)=id(i+1); node[i]=0; }
    lumped=id(6)!=0;
    thickness=data(0); fBulk=data(1); fDensity=data(2); porosity=data(3); alpha=data(4);
    for (int a=0;a<Dim;++a) {
        perm[a]=data(5+a); body[a]=data(7+a); appliedBody[a]=data(9+a);
    }
    loaded=data(11)!=0;
    setRayleighDampingFactors(data(12),data(13),data(14),data(15));
    for (int i=0;i<Dofs;++i) unbalance(i)=data(16+i);
    if (data(25)!=0) {
        if (!Kc) Kc=new Matrix(Dofs,Dofs);
        result=channel.recvMatrix(getDbTag(),commitTag,*Kc);
        if (result<0) return result;
    } else if (Kc) { delete Kc; Kc=0; }
    if (!material || material->getClassTag()!=id(4)) {
        delete material;
        material=broker.getNewNDMaterial(id(4));
    }
    if (!material) return -1;
    material->setDbTag(id(5));
    volume=0;
    DomainComponent::setDomain(0);
    result=material->recvSelf(commitTag,channel,broker);
    if (result>=0 && previousDomain) setDomain(previousDomain);
    return result;
}
void Tri31UP::Print(OPS_Stream &s, int flag) {
    s << "Tri31UP " << getTag() << " nodes: " << nodes << "thickness: " << thickness
      << " volume: " << volume << " alpha: " << alpha << " (experimental P1/P1 plane strain)" << endln;
}
int Tri31UP::displaySelf(Renderer &viewer, int mode, float fact, const char **modes, int numModes)
{
    if (!valid) return -1;
    Matrix coords(NumNodes,3);
    Vector values(NumNodes), point(3);
    double value=0;
    if (mode>0 && mode<=Stress) value=material->getStress()(mode-1);
    for (int i=0;i<NumNodes;++i) {
        node[i]->getDisplayCrds(point,fact,mode);
        for (int a=0;a<3;++a) coords(i,a)=point(a);
        values(i)=value;
    }
    return viewer.drawPolygon(coords,values,getTag());
}
Response *Tri31UP::setResponse(const char **argv, int argc, OPS_Stream &s)
{
    if (argc<1 || !material) return 0;
    if (!std::strcmp(argv[0],"force") || !std::strcmp(argv[0],"forces"))
        return new ElementResponse(this,1,Vector(Dofs));
    if (!std::strcmp(argv[0],"stress") || !std::strcmp(argv[0],"stress2D3"))
        return new ElementResponse(this,2,Vector(Stress));
    if (!std::strcmp(argv[0],"strain") || !std::strcmp(argv[0],"strain2D3"))
        return new ElementResponse(this,3,Vector(Stress));
    if (!std::strcmp(argv[0],"mass")) return new ElementResponse(this,4,Matrix(Dofs,Dofs));
    if (!std::strcmp(argv[0],"damping")) return new ElementResponse(this,5,Matrix(Dofs,Dofs));
    if (!std::strcmp(argv[0],"stiffness")) return new ElementResponse(this,6,Matrix(Dofs,Dofs));
    if (!std::strcmp(argv[0],"initialStiffness")) return new ElementResponse(this,7,Matrix(Dofs,Dofs));
    if (!std::strcmp(argv[0],"dynamicForce")) return new ElementResponse(this,8,Vector(Dofs));
    if (!std::strcmp(argv[0],"material")) {
        int skip=(argc>2 && !std::strcmp(argv[1],"1")) ? 2 : 1;
        return material->setResponse(argv+skip,argc-skip,s);
    }
    return material->setResponse(argv,argc,s);
}
int Tri31UP::getResponse(int id, Information &info)
{
    if (!material) return -1;
    switch (id) {
    case 1: return info.setVector(getResistingForce());
    case 2: return info.setVector(material->getStress());
    case 3: return info.setVector(material->getStrain());
    case 4: return info.setMatrix(getMass());
    case 5: return info.setMatrix(getDamp());
    case 6: return info.setMatrix(getTangentStiff());
    case 7: return info.setMatrix(getInitialStiff());
    case 8: return info.setVector(getResistingForceIncInertia());
    default: return -1;
    }
}
int Tri31UP::setParameter(const char **argv, int argc, Parameter &param)
{
    if (argc<1 || !material) return -1;
    const char *names[]={"xPerm","yPerm","alpha"};
    for (int i=0;i<3;++i) if (!std::strcmp(argv[0],names[i])) return param.addObject(i+1,this);
    return material->setParameter(argv,argc,param);
}
int Tri31UP::updateParameter(int id, Information &info)
{
    if (id>=1 && id<=3) {
        if (!std::isfinite(info.theDouble) || info.theDouble<0) return -1;
        if (id==3) alpha=info.theDouble;
        else perm[id-1]=info.theDouble;
        return 0;
    }
    return -1;
}
