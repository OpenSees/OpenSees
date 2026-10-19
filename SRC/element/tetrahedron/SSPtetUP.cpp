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
// 2026 By Pedro Arduino @ University of Washington, WA and
// Marlio Molina Manrique @ Politecnico di Torino, Turin, Italy
// parduino@uw.edu | marlio.molina@polito.it
// ============================================================================
// Please read detailed description in SSPtetUP.h.
// ============================================================================

#include "SSPtetUP.h"
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
#include <cmath>
#include <cstring>
#include <cfloat>

void *OPS_SSPtetUP()
{
    static int num_SSPtetUP = 0;
    if (num_SSPtetUP == 0) {
        num_SSPtetUP++;
        opserr << "SSPtetUP element - Written: P.Arduino, M.Molina_Manrique, U.Washington\n";
    }

    const int count = OPS_GetNumRemainingInputArgs();
    if (OPS_GetNDM()!=3 || OPS_GetNDF()!=4 || count<13) {
        opserr << "SSPtetUP requires ndm=3, ndf=4: element SSPtetUP tag n1 n2 n3 n4 matTag "
               << "fBulk fDensity kx ky kz voidRatio alpha <bx by bz> <-lumped>\n";
        return 0;
    }
    int ids[6], n=6;
    double d[10]={0};
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
            if (nb==3) { opserr << "SSPtetUP: too many body components\n"; return 0; }
            d[7+nb++]=value;
        } else {
            if (OPS_GetNumRemainingInputArgs()<before) OPS_ResetCurrentInputArg(-1);
            const char *option=OPS_GetString();
            if (std::strcmp(option,"-lumped")!=0 || lumped) {
                opserr << "SSPtetUP: unknown or repeated option " << option << endln;
                return 0;
            }
            lumped=true;
        }
    }
    if (nb!=0 && nb!=3) { opserr << "SSPtetUP: supply all three body components\n"; return 0; }
    for (int i=0;i<10;++i) if (!std::isfinite(d[i])) return 0;
    if (d[0]<=0 || d[1]<0 || d[2]<0 || d[3]<0 || d[4]<0 || d[5]<0 || d[6]<0) {
        opserr << "SSPtetUP: bulk must be positive; density, permeability, void ratio and alpha nonnegative\n";
        return 0;
    }
    for (int i=0;i<4;++i) for (int j=0;j<i;++j) if (ids[i+1]==ids[j+1]) {
        opserr << "SSPtetUP: repeated node\n"; return 0;
    }
    NDMaterial *mat=OPS_getNDMaterial(ids[5]);
    if (!mat) { opserr << "SSPtetUP: material not found\n"; return 0; }
    // Check before construction so unsupported materials never terminate the process.
    NDMaterial *copy=mat->getCopy("ThreeDimensional");
    if (!copy) return 0;
    const bool compatible=copy->getOrder()==6;
    delete copy;
    if (!compatible) return 0;
    SSPtetUP *element=new SSPtetUP(ids[0],ids[1],ids[2],ids[3],ids[4],*mat,
        d[0],d[1],d[2],d[3],d[4],d[5],d[6],d[7],d[8],d[9],lumped);
    // Domain::addElement does not propagate setDomain/update failure in release
    // builds. Validate here so a bad geometry cannot silently enter the model.
    element->setDomain(OPS_GetDomain());
    if (!element->getDomain()) { delete element; return 0; }
    return element;
}

SSPtetUP::SSPtetUP()
    : Element(0,ELE_TAG_SSPtetUP), nodes(4), material(0),
      fBulk(1), fDensity(0), porosity(0), alpha(0), volume(0),
      lumped(false), loaded(false), valid(false), grad(4,3), B(6,12),
      stiffness(16,16), initial(16,16), mass(16,16), damp(16,16),
      residual(16), unbalance(16)
{
    for (int i=0;i<4;++i) node[i]=0;
    for (int a=0;a<3;++a) perm[a]=body[a]=appliedBody[a]=0;
}

SSPtetUP::SSPtetUP(int tag, int n1, int n2, int n3, int n4, NDMaterial &mat,
    double bulk, double fluidDensity, double kx, double ky, double kz,
    double voidRatio, double stabilization, double bx, double by, double bz, bool lump)
    : SSPtetUP()
{
    setTag(tag);
    nodes(0)=n1; nodes(1)=n2; nodes(2)=n3; nodes(3)=n4;
    material=mat.getCopy("ThreeDimensional");
    fBulk=bulk; fDensity=fluidDensity; porosity=voidRatio/(1+voidRatio);
    alpha=stabilization; lumped=lump;
    perm[0]=kx; perm[1]=ky; perm[2]=kz;
    body[0]=bx; body[1]=by; body[2]=bz;
}

SSPtetUP::~SSPtetUP() { delete material; }

void SSPtetUP::setDomain(Domain *domain)
{
    valid=false;
    volume=0;
    DomainComponent::setDomain(0);
    for (int i=0;i<4;++i) node[i]=0;
    if (!domain || !material) return;
    for (int i=0;i<4;++i) {
        node[i]=domain->getNode(nodes(i));
        if (!node[i] || node[i]->getNumberDOF()!=4 || node[i]->getCrds().Size()!=3) {
            opserr << "SSPtetUP " << getTag() << ": requires four 3D nodes with four DOFs\n";
            return;
        }
    }
    // J=[x2-x1,x3-x1,x4-x1]. Translate first for robustness to large offsets.
    double J[3][3], cofactor[3][3], scale=0;
    for (int a=0;a<3;++a) for (int j=0;j<3;++j) {
        J[a][j]=node[j+1]->getCrds()(a)-node[0]->getCrds()(a);
        if (std::fabs(J[a][j])>scale) scale=std::fabs(J[a][j]);
    }
    for (int a=0;a<3;++a) for (int j=0;j<3;++j)
        cofactor[a][j]=J[(a+1)%3][(j+1)%3]*J[(a+2)%3][(j+2)%3]
                      -J[(a+1)%3][(j+2)%3]*J[(a+2)%3][(j+1)%3];
    double det=0;
    for (int j=0;j<3;++j) det+=J[0][j]*cofactor[0][j];
    if (!std::isfinite(det) || scale==0 || std::fabs(det)<=64*DBL_EPSILON*scale*scale*scale) {
        opserr << "SSPtetUP " << getTag() << ": degenerate tetrahedron\n";
        return;
    }
    // Both orientations are accepted; integration uses positive physical volume.
    volume=std::fabs(det)/6;
    grad.Zero(); B.Zero();
    for (int a=0;a<3;++a) for (int i=1;i<4;++i) {
        grad(i,a)=cofactor[a][i-1]/det;
        grad(0,a)-=grad(i,a);
    }
    for (int i=0;i<4;++i) {
        B(0,3*i)=grad(i,0); B(1,3*i+1)=grad(i,1); B(2,3*i+2)=grad(i,2);
        B(3,3*i)=grad(i,1); B(3,3*i+1)=grad(i,0);
        B(4,3*i+1)=grad(i,2); B(4,3*i+2)=grad(i,1);
        B(5,3*i)=grad(i,2); B(5,3*i+2)=grad(i,0);
    }
    valid=true;
    DomainComponent::setDomain(domain);
}

int SSPtetUP::update()
{
    if (!valid) return -1;
    Vector u(12), strain(6);
    for (int i=0;i<4;++i) for (int a=0;a<3;++a) u(3*i+a)=node[i]->getTrialDisp()(a);
    strain.addMatrixVector(0,B,u,1);
    return material->setTrialStrain(strain);
}
int SSPtetUP::commitState() {
    int result=Element::commitState();
    return material ? result+material->commitState() : -1;
}
int SSPtetUP::revertToLastCommit() { return material ? material->revertToLastCommit() : -1; }
int SSPtetUP::revertToStart() { return material ? material->revertToStart() : -1; }

const Matrix &SSPtetUP::formStiffness(const Matrix &tangent)
{
    Matrix solid(12,12);
    solid.addMatrixTripleProduct(0,B,tangent,volume);
    stiffness.Zero();
    for (int i=0;i<4;++i) for (int j=0;j<4;++j)
        for (int a=0;a<3;++a) for (int b=0;b<3;++b)
            stiffness(4*i+a,4*j+b)=solid(3*i+a,3*j+b);
    return stiffness;
}
const Matrix &SSPtetUP::getTangentStiff() {
    if (!valid) { stiffness.Zero(); return stiffness; }
    return formStiffness(material->getTangent());
}
const Matrix &SSPtetUP::getInitialStiff() {
    initial.Zero();
    if (valid) initial=formStiffness(material->getInitialTangent());
    return initial;
}

const Matrix &SSPtetUP::getMass()
{
    mass.Zero();
    if (!valid) return mass;
    const double rho=material->getRho(); // supplied material density is mixture density
    for (int i=0;i<4;++i) for (int j=0;j<4;++j) {
        const double nij=volume*(i==j ? 2.0 : 1.0)/20;
        const double solidMass=rho*(lumped ? (i==j ? volume/4 : 0) : nij);
        for (int a=0;a<3;++a) mass(4*i+a,4*j+a)=solidMass;
        double dot=0;
        for (int a=0;a<3;++a) dot+=grad(i,a)*grad(j,a);
        // Lump only the displacement mass: row-sum lumping destroys stabilization.
        mass(4*i+3,4*j+3)=-porosity/fBulk*nij-alpha*volume*dot;
    }
    return mass;
}

const Matrix &SSPtetUP::getDamp()
{
    damp.Zero();
    if (!valid) return damp;
    // Rayleigh damping applies only to the displacement block.
    if (alphaM!=0) {
        const Matrix &M=getMass();
        for (int i=0;i<16;++i) if (i%4!=3)
            for (int j=0;j<16;++j) if (j%4!=3) damp(i,j)+=alphaM*M(i,j);
    }
    if (betaK!=0) damp.addMatrix(1,getTangentStiff(),betaK);
    if (betaK0!=0) damp.addMatrix(1,getInitialStiff(),betaK0);
    if (betaKc!=0 && Kc) damp.addMatrix(1,*Kc,betaKc);
    for (int i=0;i<4;++i) for (int j=0;j<4;++j) {
        double h=0;
        for (int a=0;a<3;++a) {
            damp(4*i+a,4*j+3)=-volume*grad(i,a)/4;
            damp(4*j+3,4*i+a)=-volume*grad(i,a)/4;
            h+=grad(i,a)*perm[a]*grad(j,a);
        }
        damp(4*i+3,4*j+3)=-volume*h;
    }
    return damp;
}

void SSPtetUP::zeroLoad() {
    unbalance.Zero(); loaded=false;
    for (int a=0;a<3;++a) appliedBody[a]=0;
}
int SSPtetUP::addLoad(ElementalLoad *load, double factor)
{
    int type;
    const Vector &data=load->getData(type,factor);
    if (type!=LOAD_TAG_BrickSelfWeight && type!=LOAD_TAG_SelfWeight) {
        opserr << "SSPtetUP: unsupported elemental load\n"; return -1;
    }
    loaded=true;
    for (int a=0;a<3;++a)
        appliedBody[a]+=body[a]*factor*(type==LOAD_TAG_SelfWeight ? data(a) : 1);
    return 0;
}
int SSPtetUP::addInertiaLoadToUnbalance(const Vector &accel)
{
    if (!valid) return -1;
    Vector r(16);
    for (int i=0;i<4;++i) {
        const Vector &rv=node[i]->getRV(accel);
        if (rv.Size()!=4) return -1;
        for (int a=0;a<3;++a) r(4*i+a)=rv(a);
        r(4*i+3)=0; // ground excitation never acts on pressure storage
    }
    unbalance.addMatrixVector(1,getMass(),r,-1);
    return 0;
}
const Vector &SSPtetUP::getResistingForce()
{
    residual.Zero();
    if (!valid) return residual;
    Vector force(12);
    force.addMatrixTransposeVector(0,B,material->getStress(),volume);
    const double *g=loaded ? appliedBody : body;
    for (int i=0;i<4;++i) {
        for (int a=0;a<3;++a) {
            residual(4*i+a)=force(3*i+a)-material->getRho()*volume*g[a]/4;
            residual(4*i+3)+=volume*fDensity*grad(i,a)*perm[a]*g[a];
        }
    }
    residual.addVector(1,unbalance,-1);
    return residual;
}
const Vector &SSPtetUP::getResistingForceIncInertia()
{
    getResistingForce();
    if (!valid) return residual;
    Vector velocity(16), acceleration(16);
    for (int i=0;i<4;++i) for (int a=0;a<4;++a) {
        velocity(4*i+a)=node[i]->getTrialVel()(a);
        acceleration(4*i+a)=node[i]->getTrialAccel()(a);
    }
    residual.addMatrixVector(1,getMass(),acceleration,1);
    residual.addMatrixVector(1,getDamp(),velocity,1);
    return residual;
}

int SSPtetUP::sendSelf(int commitTag, Channel &channel)
{
    if (!material) return -1;
    if (!material->getDbTag()) material->setDbTag(channel.getDbTag());
    ID id(8);
    id(0)=getTag();
    for (int i=0;i<4;++i) id(i+1)=nodes(i);
    id(5)=material->getClassTag(); id(6)=material->getDbTag(); id(7)=lumped ? 1 : 0;
    Vector data(35);
    data(0)=fBulk; data(1)=fDensity; data(2)=porosity; data(3)=alpha;
    for (int a=0;a<3;++a) {
        data(4+a)=perm[a]; data(7+a)=body[a]; data(10+a)=appliedBody[a];
    }
    data(13)=loaded ? 1 : 0;
    data(14)=alphaM; data(15)=betaK; data(16)=betaK0; data(17)=betaKc;
    for (int i=0;i<16;++i) data(18+i)=unbalance(i);
    data(34)=Kc ? 1 : 0;
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
int SSPtetUP::recvSelf(int commitTag, Channel &channel, FEM_ObjectBroker &broker)
{
    Domain *previousDomain=getDomain();
    valid=false;
    ID id(8); Vector data(35);
    int result=channel.recvID(getDbTag(),commitTag,id);
    if (result<0) return result;
    result=channel.recvVector(getDbTag(),commitTag,data);
    if (result<0) return result;
    setTag(id(0));
    for (int i=0;i<4;++i) { nodes(i)=id(i+1); node[i]=0; }
    lumped=id(7)!=0;
    fBulk=data(0); fDensity=data(1); porosity=data(2); alpha=data(3);
    for (int a=0;a<3;++a) {
        perm[a]=data(4+a); body[a]=data(7+a); appliedBody[a]=data(10+a);
    }
    loaded=data(13)!=0;
    setRayleighDampingFactors(data(14),data(15),data(16),data(17));
    for (int i=0;i<16;++i) unbalance(i)=data(18+i);
    if (data(34)!=0) {
        if (!Kc) Kc=new Matrix(16,16);
        result=channel.recvMatrix(getDbTag(),commitTag,*Kc);
        if (result<0) return result;
    } else if (Kc) { delete Kc; Kc=0; }
    if (!material || material->getClassTag()!=id(5)) {
        delete material;
        material=broker.getNewNDMaterial(id(5));
    }
    if (!material) return -1;
    material->setDbTag(id(6));
    valid=false; volume=0;
    DomainComponent::setDomain(0);
    result=material->recvSelf(commitTag,channel,broker);
    // Domain::recvSelf reuses existing elements when geometry has not changed;
    // in that path it does not call addElement/setDomain again.
    if (result>=0 && previousDomain) setDomain(previousDomain);
    return result;
}

void SSPtetUP::Print(OPS_Stream &s, int flag) {
    s << "SSPtetUP " << getTag() << " nodes: " << nodes << "volume: " << volume
      << " alpha: " << alpha << " (experimental P1/P1)" << endln;
}
Response *SSPtetUP::setResponse(const char **argv, int argc, OPS_Stream &s)
{
    if (argc<1 || !material) return 0;
    if (!std::strcmp(argv[0],"force") || !std::strcmp(argv[0],"forces"))
        return new ElementResponse(this,1,Vector(16));
    if (!std::strcmp(argv[0],"stress") || !std::strcmp(argv[0],"stress3D6"))
        return new ElementResponse(this,2,Vector(6));
    if (!std::strcmp(argv[0],"strain") || !std::strcmp(argv[0],"strain3D6"))
        return new ElementResponse(this,3,Vector(6));
    if (!std::strcmp(argv[0],"mass")) return new ElementResponse(this,4,Matrix(16,16));
    if (!std::strcmp(argv[0],"damping")) return new ElementResponse(this,5,Matrix(16,16));
    if (!std::strcmp(argv[0],"stiffness")) return new ElementResponse(this,6,Matrix(16,16));
    if (!std::strcmp(argv[0],"material")) {
        // Accept both "material stress" and "material 1 stress" (one point).
        int skip=(argc>2 && !std::strcmp(argv[1],"1")) ? 2 : 1;
        return material->setResponse(argv+skip,argc-skip,s);
    }
    return material->setResponse(argv,argc,s);
}
int SSPtetUP::getResponse(int id, Information &info)
{
    if (!material) return -1;
    switch (id) {
    case 1: return info.setVector(getResistingForce());
    case 2: return info.setVector(material->getStress());
    case 3: return info.setVector(material->getStrain());
    case 4: return info.setMatrix(getMass());
    case 5: return info.setMatrix(getDamp());
    case 6: return info.setMatrix(getTangentStiff());
    default: return -1;
    }
}
int SSPtetUP::setParameter(const char **argv, int argc, Parameter &param)
{
    if (argc<1 || !material) return -1;
    const char *names[]={"xPerm","yPerm","zPerm","alpha"};
    for (int i=0;i<4;++i) if (!std::strcmp(argv[0],names[i])) return param.addObject(i+1,this);
    return material->setParameter(argv,argc,param);
}
int SSPtetUP::updateParameter(int id, Information &info)
{
    if (id>=1 && id<=4) {
        if (!std::isfinite(info.theDouble) || info.theDouble<0) return -1;
        if (id==4) alpha=info.theDouble;
        else perm[id-1]=info.theDouble;
        return 0;
    }
    return -1;
}
