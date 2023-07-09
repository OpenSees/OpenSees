#include <Python.h>
#include <elementAPI.h>
#include <ElasticMaterial.cpp>

#include <UniaxialMaterial.h>
#include <string.h>
#include <string>
#include <ID.h>

#include <StandardStream.h>
StandardStream sserr;
OPS_Stream *opserrPtr = &sserr;

#include <SimulationInformation.h>
SimulationInformation simulationInfo;
SimulationInformation *theSimulationInfoPtr = 0;

#include <FE_Datastore.h>
FE_Datastore *theDatabase  =0;

#include <Domain.h>
static Domain theDomain;

#include "PythonModelBuilder.h"
#include "PythonAnalysisBuilder.h"

static UniaxialMaterial *theTestingUniaxialMaterial =0;
// static double ops_strain = 0;

extern UniaxialMaterial *OPS_ParseUniaxialMaterialCommand(const char *matType);
extern "C" void OPS_ResetCommandLine(int nArgs, int cArg, PyObject *args);

static PyObject *ops_UniaxialMaterialCommand(PyObject *self, PyObject *args)
{
  OPS_ResetCommandLine(PyTuple_Size(args), 0, args);
  int numberArgs = OPS_GetNumRemainingInputArgs();

  if (numberArgs < 2) {
    PyErr_SetString(PyExc_RuntimeError, "ERROR too few arguments: uniaxialMaterial type? tag? args.");
    return NULL;
  }

  const char *matType = OPS_GetString();
  UniaxialMaterial *theMaterial = OPS_ParseUniaxialMaterialCommand(matType);

  if (theMaterial == 0) {
    PyErr_SetString(PyExc_RuntimeError, "ERROR could not create uniaxialMaterial.");
    return NULL;
  }
  
  // Now add the material to the modelBuilder
  if (OPS_addUniaxialMaterial(theMaterial) == false) {
    PyErr_SetString(PyExc_RuntimeError, "ERROR could not add uniaaialMaterial.");
    delete theMaterial; // invoke the material objects destructor, otherwise mem leak
    return NULL;
  }

  PyObject *ret = Py_BuildValue("d", 0.0);
  return ret;  
}

static PyObject *ops_testUniaxialMaterial(PyObject *self, PyObject *args)
{
  OPS_ResetCommandLine(PyTuple_Size(args), 0, args);
  int numberArgs = OPS_GetNumRemainingInputArgs();

  if (numberArgs !=  1) {
    PyErr_SetString(PyExc_RuntimeError, "testUniaxialMaterial - You must provide a material tag.");
    return NULL;
  }

  int tag;
  int numData = 1;
  OPS_GetIntInput(&numData, &tag);

  // if (theTestingUniaxialMaterial != 0) 
  //   delete theTestingUniaxialMaterial;

  theTestingUniaxialMaterial=OPS_getUniaxialMaterial(tag);

  if (theTestingUniaxialMaterial == 0) {
      PyErr_SetString(PyExc_ValueError,"testUniaxialMaterial - Material Not Found.");
      return NULL;
  }
  PyObject *ret = Py_BuildValue("d", 0.0);
  return ret;  
}


static PyObject *ops_setStrain(PyObject *self, PyObject *args)
{
  int argc = PyTuple_Size(args);
  if (argc == 0) {
    return NULL;
    PyErr_SetString(PyExc_RuntimeError, "You must provide a strain value.");
  }

  PyObject *o = PyTuple_GetItem(args,0);
  if (!PyFloat_Check(o)) {
    PyErr_SetString(PyExc_ValueError,"setStrain expected a floating point value");
    return NULL;
  }

  double opsStrain = PyFloat_AS_DOUBLE(o);
  if (theTestingUniaxialMaterial !=0 ) {
    theTestingUniaxialMaterial->setTrialStrain(opsStrain);
    theTestingUniaxialMaterial->commitState();
  } else {
    PyErr_SetString(PyExc_ValueError,"setStrain WARNING no active UniaxialMaterial - use testUniaxialMaterial command");    
    return NULL;
  }
  PyObject *ret = Py_BuildValue("d", opsStrain);
  return ret;
}


static PyObject *ops_getStrain(PyObject *self, PyObject *args)
{
  if (theTestingUniaxialMaterial !=0) {
    double strain = theTestingUniaxialMaterial->getStrain();
    PyObject *ret = Py_BuildValue("d", strain);
    return ret;
  } else {
    PyErr_SetString(PyExc_ValueError,"getStrain no active uniaxialMaerial - use testUniaxialMaterial function");
    return NULL;
  }
}


static PyObject *ops_getStress(PyObject *self, PyObject *args)
{
  if (theTestingUniaxialMaterial !=0) {
    double stress = theTestingUniaxialMaterial->getStress();
    PyObject *ret = Py_BuildValue("d", stress);
    return ret;
  } else {
    PyErr_SetString(PyExc_ValueError,"getStress no active uniaxialMaerial - use testUniaxialMaterial function");
    return NULL;
  }
}


static PyObject *ops_getTangent(PyObject *self, PyObject *args)
{
  if (theTestingUniaxialMaterial !=0) {
    double tangent = theTestingUniaxialMaterial->getTangent();
    PyObject *ret = Py_BuildValue("d", tangent);
    return ret;
  } else {
    PyErr_SetString(PyExc_ValueError,"getTangent no active uniaxialMaerial - use testUniaxialMaterial function");
    return NULL;
  }
}

// static char module_docstring[] =
//   "This module add OpenSees Commands to Python";
// static char opensees_docstring[] =
//   "Blah";
static char uniaxialMaterial_docstring[] =
    "* uniaxialMaterial(matType,matTag,...)\n\n\
Create a uniaxial material with given type, tag and arguments\n\n\
matType -- a string of the material name\n\
matTag -- an integer identifying the material\n\
The type of material created and the additional arguments required \n\
depends on the matType? provided in the command:\n\n\
* uniaxialMaterial('Elastic',matTag,E[,eta[,Eneg]])\n\n\
E -- tangent\n\
eta -- damping tangent (optional, default=0.0)\n\
Eneg -- tangent in compression (optional, default=E)\n\n\
* uniaxialMaterial('ElasticPP',matTag,E,epsyP[,epsyN[,eps0]])\n\n\
E -- tangent\n\
epsyP -- strain or deformation at which material reaches plastic state in tension\n\
epsyN -- strain or deformation at which material reaches plastic state \n\
         in compression.(optional, default is tension value)\n\
eps0 -- initial strain (optional, default: zero)\n\n\
* uniaxialMaterial('Parallel',matTag,tag1,tag2,...)\n\n\
$tag1 $tag2 ... -- identification tags of materials making up the material model\n\n\
* uniaxialMaterial('Steel01',matTag,Fy,E0,b[,a1,a2,a3,a4])\n\n\
Fy -- yield strength\n\
E0 -- initial elastic tangent\n\
b -- strain-hardening ratio (ratio between post-yield tangent and initial elastic tangent)\n\
a1 -- isotropic hardening parameter, increase of compression yield envelope as proportion of \n\
      yield strength after a plastic strain of a2*(Fy/E0). (optional)\n\
a2 -- isotropic hardening parameter (see explanation under a1). (optional).\n\
a3 -- isotropic hardening parameter, increase of tension yield envelope as \n\
      proportion of yield strength after a plastic strain of a4*(Fy/E0). (optional)\n\
a4 -- isotropic hardening parameter (see explanation under a3). (optional)\n\n\
* uniaxialMaterial('Concrete01',matTag,fpc,epsc0,fpcu,epsU)\n\n\
fpc -- concrete compressive strength at 28 days (compression is negative)\n\
epsc0 -- concrete strain at maximum strength\n\
fpcu -- concrete crushing strength\n\
epsU -- concrete strain at crushing strength\n\n\
* uniaxialMaterial('Series',matTag,tag1,tag2,...)\n\n\
matTag -- integer tag identifying material\n\
tag1,tag2,... -- identification tags of materials making up the material model\n\n\
*uniaxialMaterial('Hysteretic',matTag,s1p,e1p,s2p,e2p[,s3p,e3p],s1n,e1n,s2n,e2n[,s3n,e3n],\n\
                  pinchX,pinchY,damage1,damage2[,beta])\n\n\
s1p,e1p -- stress and strain (or force & deformation) at first point of the envelope \n\
           in the positive direction\n\
s2p,e2p -- stress and strain (or force & deformation) at second point of the envelope \n\
           in the positive direction\n\
s3p,e3p -- stress and strain (or force & deformation) at third point of the envelope \n\
           in the positive direction (optional)\n\
s1n,e1n -- stress and strain (or force & deformation) at first point of the envelope \n\
           in the negative direction\n\
s2n,e2n -- stress and strain (or force & deformation) at second point of the envelope \n\
           in the negative direction\n\
s3n,e3n -- stress and strain (or force & deformation) at third point of the envelope \n\
           in the negative direction (optional)\n\
pinchx -- pinching factor for strain (or deformation) during reloading\n\
pinchy -- pinching factor for stress (or force) during reloading\n\
damage1 -- damage due to ductility: D1(mu-1)\n\
damage2 -- damage due to energy: D2(Eii/Eult)\n\
beta -- power used to determine the degraded unloading stiffness based on ductility, \n\
        mu-beta (optional, default=0.0)\n\n\
* uniaxialMaterial('Cable',tag,presetress,E,effUnitWeight,Lelement)\n\n";
static char testUMaterial_docstring[] =
    "* testUniaxialMaterial(matTag)\n\n\
Test a uniaxial material which has been defined with matTag\n";
static char setStrain_docstring[] =
    "* setStrain(double)\n\n\
Test the uniaxial material which has been set with testUniaxialMaterial()\n";
static char getStrain_docstring[] =
    "* getStrain() -> double\n\n\
Test the uniaxial material which has been set with testUniaxialMaterial()\n";
static char getStress_docstring[] =
    "* getStress() -> double\n\n\
Test the uniaxial material which has been set with testUniaxialMaterial()\n";
static char getTangent_docstring[] =
    "* getTangent() -> double\n\n\
Test the uniaxial material which has been set with testUniaxialMaterial()\n";
static char model_docstring[] =
    "* model('basic','-ndm',ndm[,'-ndf',ndf])\n\n\
Set the basic model builder with spatial dimension of model and number of DOF\n\n\
ndm -- spatial dimension of problem (1,2 or 3)\n\
ndf -- number of DOF at node, default value dpends on value of ndm\n\
       ndm=1 -> ndf=1\n\
       ndm=2 -> ndf=3\n\
       ndm=3 -> ndf=6\n";
static char node_docstring[] =
    "* node(nodeTag,(ndm coords)[,'-mass',(ndf mass)][,'-disp',(ndf disp)][,'-vel',(ndf vel)]\
[,'-dispLoc',(ndf dispLoc)]) -> nodeTag\n\n\
Create a Node object DOF\n\n\
nodeTag -- integer tag identifying node\n\
coords -- nodal coordinates (ndm arguments)\n\
mass -- nodal mass (ndf arguments)\n\
disp -- nodal displacement (ndf arguments)\n\
vel -- nodal velocity (ndf arguments)\n\
dispLoc -- nodal display location (ndm arguments)\n";
static char fix_docstring[] =
    "* fix(nodeTag,(ndf constrValues))\n\n\
Construct single-point homogeneous boundary constraints\n\n\
nodeTag -- integer tag identifying the node to be constrained\n\
constrValues -- ndf constraint values (0 or 1) corresponding to the ndf dof\n\
                0 unconstrained (or free)\n\
                1 constrained (or fix)\n";
static char element_docstring[] =
    "* element(eleType,eleTag,...) -> eleTag\n\n\
Construct an Element object\n\n\
eleType -- a string of the element name\n\
eleTag -- integer tag identifying the element\n\n\
The type of element created and the additional arguments required depends on \n\
the eleType? provided in the command:\n\n\
* element('truss',eleTag,iNode,jNode,A,matTag[,'-rho',rho][,-cMass,cFlag][,-doRayleigh $rFlag])\n\n\
iNode,jNode -- end nodes\n\
A -- cross-sectional area of element\n\
matTag -- tag associated with previously-defined UniaxialMaterial\n\
rho -- mass per unit length, optional, default = 0.0\n\
cFlag -- consistent mass flag, optional, default = 0\n\
         cFlag = 0 lumped mass matrix (default)\n\
         cFlag = 1 consistent mass matrix\n\
rFlag -- Rayleigh damping flag, optional, default = 0\n\
         rFlag = 0 NO RAYLEIGH DAMPING (default)\n\
         rFlag = 1 include Rayleigh damping\n\n\
* element('zeroLengthSection',eleTag,iNode,jNode,secTag[,'-orient',x1,x2,x3,yp1,yp2,yp3][,'-doRayleigh',rFlag])\n\n\
iNode,jNode -- end nodes\n\
secTag -- tag associated with previously-defined Section object\n\
x1,x2,x3 -- vector components in global coordinates defining local x-axis (optional)\n\
yp1,yp2,yp3 -- vector components in global coordinates defining vector yp which lies \n\
               in the local x-y plane for the element. (optional)\n\
rFlag -- optional, default = 1\n\
         rFlag = 0 no Rayleigh damping\n\
         rFlag = 1 include Rayleigh damping (default)\n\n\
* element('forceBeamColumn',eleTag,iNode,jNode,transfTag,beamintegrationTag[,'-iter',maxIter,tol][,'-mass',massDens])\n\n\
iNode, jNode -- end nodes\n\
transfTag -- identifier for previously-defined coordinate-transformation (CrdTransf) object\n\
beamintegrationTag -- identifier for previously-defined BeamIntegrationRule object\n\
massDens -- element mass density (per unit length), from which a lumped-mass matrix is \n\
            formed (optional, default=0.0)\n\
maxIters -- maximum number of iterations to undertake to satisfy element compatibility\n\
            (optional, default=10)\n\
tol -- tolerance for satisfaction of element compatibility (optional, default=10e-12)\n\n\
* element('elasticBeamColumn',eleTag,iNode,jNode,A,E,Iz,transfTag[,'-alpha',alpha][,'-depth',depth][,'-mass',massDens][,'-cMass'])\n\
For two-dimensional problem\n\n\
* element('elasticBeamColumn',eleTag,iNode,jNode,A,E,G,J,Iy,Iz,transfTag[,'-mass',massDens][,'-cMass'])\n\
* element('elasticBeamColumn',eleTag,iNode,jNode,secTag,transfTag[,'-mass',massDens][,'-cMass'])\n\
For three-dimensional problem\n\n\
iNode,jNode -- end nodes\n\
A -- cross-sectional area of element\n\
E -- Young's Modulus\n\
G -- Shear Modulus\n\
J -- torsional moment of inertia of cross section\n\
Iz -- second moment of area about the local z-axis\n\
Iy -- second moment of area about the local y-axis\n\
transfTag -- identifier for previously-defined coordinate-transformation (CrdTransf) object\n\
massDens -- element mass per unit length (optional, default = 0.0)\n\
-cMass -- to form consistent mass matrix (optional, default = lumped mass matrix)\n\
alpha -- temperature expansion coefficient\n\
depth -- depth of section\n\
secTag -- section for 3D problem\n";

static char timeSeries_docstring[] =
    "* timeSeries(tsType,tsTag,...) -> tsTag\n\n\
Construct an TimeSeries object\n\n\
tsType -- a string of the TimeSeries name\n\
tsTag -- integer tag identifying the timeseries\n\n\
The type of time series created and the additional arguments required depends on the \n\
seriesType? provided in the command:\n\n\
* timeSeries('Constant',tsTag[,'-factor',factor])\n\n\
* timeSeries('Linear',tsTag[,'-factor',factor])\n\n\
* timeSeries('Trig',tsTag,tStart,tEnd,period[,'-factor',factor]['-shift',shift])\n\n\
* timeSeries('Triangle',tsTag,tStart,tEnd,period[,'-factor',factor]['-shift',shift])\n\n\
* timeSeries('Rectangular',tsTag,tStart,tEnd[,'-factor',factor])\n\n\
* timeSeries('Pulse',tsTag,tStart,tEnd,period['-width',width]['-shift',shift][,'-factor',factor])\n\n\
tStart -- starting time of non-zero load factor\n\
tEnd -- ending time of non-zero load factor\n\
period -- characteristic period of sine wave\n\
shift -- phase shift in radians (optional, default=0.0)\n\
factor -- the load factor amplification factor (optional, default=1.0)\n\
pulseWidth -- pulse width as a fraction of the period (optional, default = 0.5)\n\n\
* timeSeries('Path',tsTag[,'-dt',dt][,'-values',...][,'-factor',factor])\n\
* timeSeries('Path',tsTag[,'-dt',dt][,'-filePath',filePath][,'-factor',factor])\n\
* timeSeries('Path',tsTag[,'-dt',dt]['-time',...][,'-values',...][,'-factor',factor])\n\
* timeSeries('Path',tsTag[,'-dt',dt][,'-fileTime',fileTime][,'-filePath',filePath][,'-factor',factor])\n\n\
tag -- unique tag among TimeSeries objects.\n\
filePath -- file containing the load factors values\n\
fileTime -- file containing the time values for corresponding load factors\n\
dT -- time interval between specified points.\n\
'-time',... -- load factor values\n\
'-values',... -- time values\n\
factor -- optional, a factor to multiply load factors by (default = 1.0)\n";
static char pattern_docstring[] =
    "* pattern(patternType,patternTag,...) -> patternTag\n\n\
Construct an LoadPattern object and set the current pattern.\n\n\
patternType -- a string of the LoadPattern name\n\
patternTag -- integer tag identifying the pattern\n\n\
The type of pattern created and the additional arguments required depends on the \n\
patternType? provided in the command:\n\n\
* pattern('Plain',patternTag,tsTag[,'-fact',factor]\n\n\
tsTag -- the tag of the time series to be used in the load pattern\n\
factor -- constant factor (optional, default=1.0)\n\n\
* pattern('MultipleSupport',patternTag[,'-fact',factor]\n\n\
* pattern('UniformExcitation',patternTag,dir[,'-accel',tsTag][,'-vel',tsTag][,'-disp',tsTag][,'-fact',factor])\n\n"
    "tsTag -- tag of TimeSeries object created using timeSeries command.\n"
    "cFactor -- constant factor (optional, default=1.0)";

static char load_docstring[] =
    "* load(nodeTag,(ndf LoadValues)[,'-const'][,'-pattern',patternTag])-> patternTag\n\n\
Construct a Nodalload object and add it to the current pattern or specified pattern.\n\n\
nodeTag -- integer tag identifying the node to which load is applied\n\
LoadValues -- ndf reference load values\n\
'-const' -- make it a constant load\n\
patternTag -- specify a pattern to add, otherwise the load is added to the current pattern\n";
static char sp_docstring[] =
    "* sp(nodeTag,dofTag,dofValue[,'-const'][,'-pattern',patternTag]) -> spTag\n\n\
Construct a sing-point constraint object and add it to the current pattern or specified pattern.\n\n\
nodeTag -- integer tag identifying the node to which load is applied\n\
dofTag -- the dof at the node to which constraint is applied (1 trhough ndf)\n\
dofValue -- reference constraint value\n\
'-const' -- make it a constant constraint\n\
patternTag -- specify a pattern to add, otherwise the load is added to the current pattern\n";
static char section_docstring[] =
    "* section(secType,secTag,...) -> secTag\n\n\
Construct a SectionForceDeformation object, hereto referred to as Section, \n\
which represents force-deformation (or resultant stress-strain) relationships at \n\
beam-column and plate sample points.\n\n\
secType -- a string of the name of Section\n\
secTag -- integer tag identifying the section\n\n\
The type of section created and the additional arguments required depends on the \n\
secType? provided in the command:\n\n\
* section('Elastic',secTag,E,A,Iz[,G,alphaY])\n\
* section('Elastic',secTag,E,A,Iz,Iy,G,J[,alphaY,alphaZ])\n\n\
E -- Young's Modulus\n\
A -- cross-sectional area of section\n\
Iz -- second moment of area about the local z-axis\n\
Iy -- second moment of area about the local y-axis (required for 3D analysis)\n\
G -- Shear Modulus (optional for 2D analysis, required for 3D analysis)\n\
J -- torsional moment of inertia of section (required for 3D analysis)\n\
alphaY -- shear shape factor along the local y-axis (optional)\n\
alphaZ -- shear shape factor along the local z-axis (optional)\n\n\
* section('Fiber',secTag)\n\n\
* section('NDFiber',secTag)\n\n";
static char nDMaterial_docstring[] =
    "* nDMaterial(matType,matTag,...) -> matTag\n\n\
Construct an NDMaterial object which represents the stress-strain relationship \n\
at the gauss-point of a continuum element.\n\n\
matType -- a string of the name of nDMaterial\n\
matTag -- integer tag identifying the nDMaterial\n\n\
The type of ndMaterial created and the additional arguments required depends on the \n\
matType? provided in the command:\n\n\
* nDMaterial('ElasticIsotropic',matTag,E,v[,rho])\n\n\
E -- elastic Modulus\n\
v -- Poisson's ratio\n\
rho -- mass density, optional default = 0.0.\n\n\
* nDMaterial('ElasticOrthotropic',matTag,Ex,Ey,Ez,vxy,vyz,vzx,Gxy,Gyz,Gzx[,rho])\n\n\
Ex,Ey,Ez -- elastic modulii in three mutually perpendicular directions (x, y, and z)\n\
vxy,vyz,vzx -- Poisson's ratios\n\
Gxy,Gyz,Gzx -- shear modulii\n\
rho -- mass density, optional default = 0.0.\n\n";
static char wipeAnalysis_docstring[] =
    "* wipeAnalysis()\n\n\
Wipe all the analysis\n";
static char wipeReliability_docstring[] =
    "* wipeReliability()\n\n\
Wipe all the reliability analysis and random variables\n";
static char wipe_docstring[] =
    "* wipe()\n\n\
Wipe all the analysis and model\n";
static char Print_docstring[] =
    "* Print([filename])\n\n\
Print all objects in the domain\n\n\
filename -- name of file to which data will be sent. overwrtes existing file. default is \n\
to print to stderr\n\n\
* Print([filename][,'-node'[,'-flag',flag][,nd1,...]]])\n\n\
flag -- integer flag to be sent to the print() method, depending node and element type\n\
nd1,... -- integer tags of nodes to be printed. default is to print all\n\n\
* Print([filename][,'-ele'[,'-flag',flag][,ele1,...]]])\n\n\
ele1,... -- integer tags of elements to be printed. default is to print all\n";
static char system_docstring[] =
    "* system(systemType,...)\n\n\
Construct the LinearSOE and LinearSolver objects to store and solve the system of equations\n\n\
systemType -- a string of the solver name\n\n\
The type of LinearSOE created and the additional arguments required depends on the \n\
systemType? provided in the command:\n\n\
* system('BandGeneral')\n\n\
* system('BandSPD')\n\n\
* system('SuperLU'[,'-piv'][,'-np',np][,'-npRow',npRow][,'-npCol',npCol])\n\n";
static char numberer_docstring[] =
    "* numberer(numbererType,...)\n\n\
Construct the DOF_Numberer objects to determine the mapping between equations numbers and dofs\n\n\
numbererType -- a string of the Numberer name\n\n\
The type of DOF_Numberer created and the additional arguments required depends on the \n\
numbererType? provided in the command:\n\n\
* numberer('Plain')\n\n\
* numberer('RCM')\n\n\
* numberer('AMD')\n\n";
static char constraints_docstring[] =
    "* constraints(constraintType,...)\n\n\
Construct the ConstraintHandler objects to determine how the constraint equations are enforced\n\n\
constraintType -- a string of the Handler name\n\n\
The following contain information about numbererType? and the args required for each of the \n\
available constraint handler types:\n\n\
* constraints('Plain')\n\n\
* constraints('Lagrange'[,alphaS,alphaM])\n\n\
alphaS -- factor on single points. optional, default = 1.0\n\
alphaM -- factor on multi-points, optional, default = 1.0;\n\n\
* constraints('Penalty',alphaS,alphaM)\n\n\
alphaS -- factor on single points\n\
alphaM -- factor on multi-points\n\n\
* constraints('Transformation')\n\n";
static char algorithm_docstring[] =
    "* algorithm(algoType,...)\n\n\
Construct the SolutionAlgorithm objects to determine the sequence of steps taken to \n\
solve the non-linear equations\n\n\
algoType -- a string of the Algorithm name\n\n\
The type of solution algorithm created and the additional arguments required depends on \n\
the algorithmType? provided in the command:\n\n\
* algorithm('Linear'[,'-initial'][,'-secent'][,'-factorOnce'])\n\n\
-secant -- optional flag to indicate to use secant stiffness\n\
-initial -- optional flag to indicate to use initial stiffness\n\
-factorOnce -- optional flag to indicate to only set up and factor matrix once\n\n\
* algorithm('Newton'[,'-initial'][,'-secent'][,'-initialThenCurrent'])\n\n\
-secant -- optional flag to indicate to use secant stiffness\n\
-initial -- optional flag to indicate to use initial stiffness\n\
-initialThenCurrent -- optional flag to indicate to use initial stiffness on first step,\n\
                       then use current stiffness for subsequent steps \n\n";
static char integrator_docstring[] =
    "* integrator(integType,...)\n\n\
Construct the Integrator objects to determine the meaning of the terms in the system \n\
of equation object Ax=B. \n\n\
integType -- a string of the Integrator name\n\n\
The type of integrator used in the analysis is dependent on whether it is a \n\
static analysis or transient analysis.\n\n\
Static integrators:\n\n\
* integrator('LoadControl',lambda[,numIter,miniLambda,maxLambda])\n\n\
lambda -- the load factor increment\n\
numIter -- the number of iterations the user would like to occur in the solution \n\
           algorithm. Optional, default = 1.0.\n\
minLambda -- the min stepsize the user will allow. optional, default = lambda\n\
maxLambda -- the max stepsize the user will allow. optional, default = lambda\n\n\
* integrator('DisplacementControl',node,dof,incr[,numIter,dUmin,dUmax])\n\n\
node -- node whose response controls solution\n\
dof -- degree of freedom at the node, valid options: 1 through ndf at node.\n\
incr -- first displacement increment dUdof\n\
numIter -- the number of iterations the user would like to occur in the solution algorithm. Optional, default = 1.0.\n\
dUmin -- the min stepsize the user will allow. optional, default dU0\n\
dUmax -- the max stepsize the user will allow. optional, default dU0\n\n\
Transient Integrators:\n\n\
* integrator('Newmark',gamma,beta)\n\n";
static char analysis_docstring[] =
    "* analysis(analysisType)\n\n\
Construct the Analysis object, which defines what type of analysis is to be performed\n\
analysisType -- a string identifying type of analysis object. Currently 4 valid options:\n\n\
* analysis('Static')\n\
* analysis('Transient')\n\
* analysis('VariableTransient')\n\
* analysis('PFEM',dtmin,dtmax[,ratio])\n\n\
dtmin -- minimum time increment\n\
dtmin -- maximum time increment\n\
ratio -- the ratio to adjust time step when last time step is failed\n";
static char test_docstring[] =
    "* test(testType,...)\n\n\
Construct a ConvergenceTest object. Certain SolutionAlgorithm objects require a \n\
ConvergenceTest object to determine if convergence has been achieved at the end of\n\
an iteration step. The convergence test is applied to the matrix equation, AX=B stored\n\
in the LinearSOE. The type of convergence test created and the additional arguments \n\
required depends on the testType? provided in the command.\n\n\
The following contain information about testType? and the args required for each of\n\
the available test types:\n\n\
* test('NormUnbalance',tol,iter[,flag[,nType[,maxincr[,maxTol]]]])\n\
* test('NormDispIncr',tol,iter[,flag[,nType[,maxTol]]])\n\n\
This command is used to construct a convergence test which uses the norm \n\
of the right hand side of the matrix equation to determine if convergence has \n\
been reached. What the right-hand-side of the matrix equation is depends on \n\
integraor and constraint handler chosen. Usually, though not always, it is equal \n\
to the unbalanced forces in the system.\n\n\
tol -- the tolerance criteria used to check for convergence\n\
iter --	the max number of iterations to check before returning failure condition\n\
flag --	optional print flag, default is 0. valid options:\n\
         0 print nothing\n\
         1 print information on norms each time test() is invoked\n\
         2 print information on norms and number of iterations at end of successful test\n\
         4 at each step it will print the norms and also the dU and R(U) vectors.\n\
         5 if it fails to converge at end of numIter it will print an error message BUT RETURN A SUCEESSFULL test\n\
nType -- optional type of norm, default is 2. (0 = max-norm, 1 = 1-norm, 2 = 2-norm, ...)\n\
maxincr -- max number of iteration that the norm increases\n\
maxTol -- max tolerance that can reach";
static char analyze_docstring[] =
    "* analyze(numStep[,dt[,dtmin,dtmax,Jd]]) -> 0 (successful) or < 0 (failed)\n\n\
This command is used to perform the analysis\n\n\
numStep -- number of analysis steps to perform\n\
dt -- time step increment. Required if transient or variable transient analysis\n\
dtmin,dtmax -- minimum and maximum time steps. \n\
               Required if a variable time step transient analysis was specified.\n\
               The variable transient analysis will change current time step if last \n\
               analysis step took more or less iterations than this to converge.\n\
Jd -- Required if a variable time step transient analysis was specified.\n";
static char nodeDisp_docstring[] =
    "* nodeDisp(nodeTag[,dof]) -> disp list\n\n\
Returns the current displacement at a specified node.\n\
If optional dof is not provided, an array containing all displacement components \n\
for every dof at the node is returned.\n\n\
nodeTag -- integer tag identifying node\n\
dof -- specific dof at the node (1 through ndf)\n\n\
EXAMPLE:\n\
disp = nodeDisp(nodeTag, 1)\n";
static char fiber_docstring[] =
    "* fiber(yLoc,zLoc,A,matTag) -> secTag\n\n\
Construct a single fiber and add it to the enclosing FiberSection or NDFiberSection\n\n\
yLoc -- y coordinate of the fiber in the section (local coordinate system)\n\
zLoc -- z coordinate of the fiber in the section (local coordinate system)\n\
A -- area of the fiber.\n\
matTag -- material tag associated with this fiber (UniaxialMaterial tag for a FiberSection and NDMaterial tag for use in an NDFiberSection).\n";
static char patch_docstring[] =
    "* patch(type,matTag...) -> secTag\n\n\
Generate a number of fibers over a cross-sectional area.\n\n\
type -- a string of patch type\n\
matTag -- tag of previously defined material (UniaxialMaterial tag for a \n\
          FiberSection or NDMaterial tag for use in an NDFiberSection)\n\n\
* patch('quad',matTag,numSubdivIJ,numSubdivJK,yI,zI,yJ,zJ,yK,zK,yL,zL)\n\n\
numSubdivIJ -- number of subdivisions (fibers) in the IJ direction.\n\
numSubdivJK -- number of subdivisions (fibers) in the JK direction.\n\
yI,zI -- y & z-coordinates of vertex I (local coordinate system)\n\
yJ,zJ -- y & z-coordinates of vertex J (local coordinate system)\n\
yK,zK -- y & z-coordinates of vertex K (local coordinate system)\n\
yL,zL -- y & z-coordinates of vertex L (local coordinate system)\n\n\
* patch('rect',matTag,numSubdivY,numSubdivZ,yI,zI,yJ,zJ)\n\n\
numSubdivY -- number of subdivisions (fibers) in the local y direction.\n\
numSubdivZ -- number of subdivisions (fibers) in the local z direction.\n\
yI,zI -- y & z-coordinates of vertex I (local coordinate system)\n\
yJ,zJ -- y & z-coordinates of vertex J (local coordinate system)\n\n\
* patch('circ',matTag,numSubdivCirc,numSubdivRad,yCenter,zCenter,intRad,extRad,startAng,endAng)\n\n\
numSubdivCirc -- number of subdivisions (fibers) in the circumferential direction\n\
numSubdivRad -- number of subdivisions (fibers) in the radial direction.\n\
yCenter,zCenter -- y & z-coordinates of the center of the circle\n\
intRad -- internal radius\n\
extRad -- external radius\n\
startAng -- starting angle\n\
endAng -- ending angle\n";
static char layer_docstring[] =
    "* layer(type,matTag,...) -> secTag\n\n\
Generate a number of fibers along a line or a circular arc.\n\n\
type -- a string of patch type\n\
matTag -- tag of previously defined material (UniaxialMaterial tag for a \n\
          FiberSection or NDMaterial tag for use in an NDFiberSection)\n\n\
* layer('straight',matTag,numFiber,areaFiber,yStart,zStart,yEnd,zEnd)\n\n\
numFibers -- number of fibers along line\n\
areaFiber -- area of each fiber\n\
yStart,zEnd -- y and z-coordinates of first fiber in line (local coordinate system)\n\
yEnd,zEnd -- y and z-coordinates of last fiber in line (local coordinate system)\n\n\
* layer('circ',matTag,numFiber,areaFiber,yCenter,zCenter,radius[,startAng,endAng])\n\n\
numFiber -- number of fibers along arc\n\
areaFiber -- area of each fiber\n\
yCenter,zCenter -- y and z-coordinates of center of circular arc\n\
radius -- radius of circlular arc\n\
startAng -- starting angle (optional, default = 0.0)\n\
endAng -- ending angle (optional, default = 360.0 - 360/numFiber)\n";
static char getLoadFactor_docstring[] =
    "* getLoadFactor(patternTag) -> double\n\n\
Get load factor of the pattern\n\n\
patternTag -- integer tag of a pattern\n";
static char geomTransf_docstring[] =
    "* geomTransf(transfType,transfTag,...)\n\n\
Construct a coordinate-transformation (CrdTransf) object, which transforms beam \n\
element stiffness and resisting force from the basic system to the global-coordinate system. \n\n\
transfType -- a string of transformation name\n\
transfTag -- an integer of tag\n\n\
* geomTransf('Linear',transfTag[,'-jntOffset',dXi,dYi,dXj,dYj])\n\
* geomTransf('PDelta',transfTag[,'-jntOffset',dXi,dYi,dXj,dYj])\n\
* geomTransf('Corotational',transfTag[,'-jntOffset',dXi,dYi,dXj,dYj])\n\
For 2D problem.\n\n\
* geomTransf('Linear',transfTag,vecxzX,vecxzY,vecxzZ[,'-jntOffset',dXi,dYi,dZi,dXj,dYj,dZj])\n\
* geomTransf('PDelta',transfTag,vecxzX,vecxzY,vecxzZ[,'-jntOffset',dXi,dYi,dZi,dXj,dYj,dZj])\n\
* geomTransf('Corotational',transfTag,vecxzX,vecxzY,vecxzZ[,'-jntOffset',dXi,dYi,dZi,dXj,dYj,dZj])\n\
For 3D problem.\n\n\
vecxzX,vecxzY,vecxzZ -- X, Y, and Z components of vecxz, the vector used to define the local x-z \n\
                        plane of the local-coordinate system. The local y-axis is defined by \n\
                        taking the cross product of the vecxz vector and the x-axis. These components \n\
                        are specified in the global-coordinate system X,Y,Z and define a vector that is \n\
                        in a plane parallel to the x-z plane of the local-coordinate system. These items \n\
                        need to be specified for the three-dimensional problem.\n\
dXi,dYi,dZi -- joint offset values, offsets specified with respect to the global coordinate system for \n\
               element-end node i (the number of arguments depends on the dimensions of the current model).\n\
               The offset vector is oriented from node i to node j as shown in a figure below. (optional)\n\
dXj,dYj,dZj -- joint offset values, offsets specified with respect to the global coordinate system for \n\
               element-end node j (the number of arguments depends on the dimensions of the current model).\n\
               The offset vector is oriented from node j to node i as shown in a figure below. (optional)\n";

static char beamIntegration_docstring[] =
    "* beamIntegrationRule(type,integrationTag,...) -> integrationTag\n\n\
Construct a BeamIntegration object for ForceBasedBeamColumn\n\n\
type -- a string of BeamIntegration type\n\
integrationTag -- an integer of BeamIntegration tag\n\n\
* beamIntegrationRule('Lobatto',tag,secTag,N)\n\
* beamIntegrationRule('Legendre',tag,secTag,N)\n\
* beamIntegrationRule('NewtonCotes',tag,secTag,N)\n\
* beamIntegrationRule('Radau',tag,secTag,N)\n\
* beamIntegrationRule('Trapezoidal',tag,secTag,N)\n\
* beamIntegrationRule('CompositeSimpson',tag,secTag,N)\n\n\
secTag -- section tag\n\
N -- number of integration points\n\n\
* beamIntegrationRule('UserDefined',tag,N,secTag1,...,pt1,...,wt1,...)\n\
* beamIntegrationRule('FixedLocation',tag,N,secTag1,...,pt1,...)\n\
* beamIntegrationRule('LowOrder',tag,N,secTag1,...,pt1,...,wt1,...)\n\
* beamIntegrationRule('MidDistance',tag,N,secTag1,...,pt1,...)\n\n\
N -- number of integration points\n\
secTag1,... -- N section tags\n\
pt1,... -- N integration locations\n\
wt1,... -- N weights (for 'LowOrder', Nc<N weights)\n\n\
* beamIntegrationRule('UserHinge',tag,secTagE,npL,secTagL1,...,ptL1,...,wtL1,...,npR,secTagR1,...,ptR1,...,wtR1,...)\n\n\
secTagE -- section tag for element interior\n\
npL-- number of integration points at left end\n\
secTagL1,... -- npL section tags\n\
ptL1,... -- npL integration locations\n\
wtL1,... -- npL weights\n\
npR-- number of integration points at right end\n\
secTagR1,... -- npR section tags\n\
ptR1,... -- npR integration locations\n\
wtR1,... -- npR weights\n\n\
* beamIntegrationRule('HingeMidpoint',tag,secTagI,lpI,secTagJ,lpJ,secTagE)\n\
* beamIntegrationRule('HingeRadau',tag,secTagI,lpI,secTagJ,lpJ,secTagE)\n\
* beamIntegrationRule('HingeRadauTwo',tag,secTagI,lpI,secTagJ,lpJ,secTagE)\n\
* beamIntegrationRule('HingeEndpoint',tag,secTagI,lpI,secTagJ,lpJ,secTagE)\n\n\
secTagI -- section tag for end I\n\
lpI -- plastic hinge length at end I\n\
secTagJ -- section tag for end J\n\
lpJ -- plastic hinge length at end J\n\
secTagE -- section tag for element interior\n";

static char mass_docstring[] =
    "* mass(nodeTag,(ndf massValues))\n\n"
    "Set the mass at a node\n\n"
    "nodeTag -- integer tag identifying node whose mass is set\n"
    "massValues -- ndf nodal mass values corresponding to each DOF\n";

static char loadConst_docstring[] =
    "* loadConst(['-time',pseudoTime])\n\n"
    "Set the loads constant in the domain and to also set the time in the domain.\n\n"
    "pseudoTime -- Time domain is to be set to (optional)\n";

static char groundMotion_docstring[] =
    "* groundMotion(gmTag,gmType,...) -> gmTag\n\n"
    "Construct a plain GroundMotion object\n\n"
    "gmTag -- unique tag among ground motions in load pattern\n"
    "gmType -- a string of ground motion type\n\n"
    "* groundMotion(gmTag,'Plain'[,'-accel',tsTag][,'-vel',tsTag][,'-disp',tsTag][,'-fact',factor])\n\n"
    "tsTag -- tag of TimeSeries object created using timeSeries command."
    "cFactor -- constant factor (optional, default=1.0)\n\n"
    "* groundMotion(gmTag,'Interpolated',gmTag1,gmTag2,...[,'-fact',fact1,fact2,...])\n\n"
    "gmTag1,gmTag2,... -- the tags of existing ground motions in pattern to be used for interpolation.\n"
    "fact1,fact2,... -- the interpolation factors.\n";

static char rayleigh_docstring[] =
    "* rayleigh(alphaM,betaK,betaK0,betaKc)\n\n"
    "Assign damping to all previously-defined elements and nodes\n\n"
    "alphaM -- factor applied to elements or nodes mass matrix\n"
    "betaK -- factor applied to elements current stiffness matrix.\n"
    "betaKinit -- factor applied to elements initial stiffness matrix.\n"
    "betaKcomm -- factor applied to elements committed stiffness matrix.\n";

static char eigen_docstring[] =
    "* eigen([solver,]numEigenvalues) -> list eigen values\n\n"
    "The eigenvectors are stored at the nodes and can be printed out using\n"
    "the nodeEigenvector command, or the Print command. The default eigensolver\n"
    "is able to solve only for N-1 eigenvalues, where N is the number of\n"
    "inertial DOFs. When running into this limitation the -fullGenLapack \n"
    "solver can be used instead of the default Arpack solver.\n\n"
    "numEigenvalues -- number of eigenvalues required\n"
    "solver -- optional string detailing type of solver: -genBandArpack,\n"
    "          -symmBandLapack, -fullGenLapack (default: -genBandArpack)";

static char eleLoad_docstring[] =
    "* eleLoad(['-ele',eleTag1,eleTag2,...][,'-type',...][,'-range',eleTag1,eleTag2]) -> patternTag\n\n"
    "Construct an ElementalLoad object and add it to the enclosing LoadPattern\n\n"
    "eleTag1,eleTag2,... -- tag of previously defined element\n\n"
    "* eleLoad(,'-type','-beamUniform',Wy[,Wx])\n"
    "* eleLoad(,'-type','-beamPoint',Py,xL[,Px])\n\n"
    "Two-dimentional problem\n\n"
    "* eleLoad(,'-type','-beamUniform',Wy,Wz[,Wx])\n"
    "* eleLoad(,'-type','-beamPoint',Py,Pz,xL[,Px])\n\n"
    "Three-dimentional problem\n\n"
    "Wx -- mag of uniformily distributed ref load acting in direction along member length\n"
    "Wy -- mag of uniformily distributed ref load acting in local y direction of element\n"
    "Wz -- mag of uniformily distributed ref load acting in local z direction of element\n"
    "Py -- mag of ref point load acting in direction along member length\n"
    "Py -- mag of ref point load acting in local y direction of element\n"
    "Pz -- mag of ref point load acting in local z direction of element\n"
    "xL -- location of point load relative to node I, prescribed as fraction of element length";

static PyMethodDef methodsOpenSees[] = {
  // model methods
  {"uniaxialMaterial", ops_UniaxialMaterialCommand, METH_VARARGS, uniaxialMaterial_docstring},
  {"testUniaxialMaterial", ops_testUniaxialMaterial, METH_VARARGS, testUMaterial_docstring},
  {"setStrain", ops_setStrain, METH_VARARGS, setStrain_docstring},
  {"getStrain", ops_getStrain, METH_VARARGS, getStrain_docstring},
  {"getStress", ops_getStress, METH_VARARGS, getStress_docstring},
  {"getTangent", ops_getTangent, METH_VARARGS, getTangent_docstring},
  {"model", ops_Model, METH_VARARGS, model_docstring},
  {"node", ops_addNode, METH_VARARGS, node_docstring},
  {"fix", ops_addHomogeneousBC, METH_VARARGS, fix_docstring},
  {"element", ops_addElement, METH_VARARGS, element_docstring},
  {"timeSeries", ops_addTimeSeries, METH_VARARGS, timeSeries_docstring},
  {"pattern", ops_addPattern, METH_VARARGS, pattern_docstring},
  {"load", ops_addNodalLoad, METH_VARARGS, load_docstring},
  {"sp", ops_addSP, METH_VARARGS, sp_docstring},
  {"section", ops_addSection, METH_VARARGS, section_docstring},
  {"nDMaterial", ops_addNDMaterial, METH_VARARGS, nDMaterial_docstring},
  {"fiber", ops_addFiber, METH_VARARGS, fiber_docstring},
  {"patch", ops_addPatch, METH_VARARGS, patch_docstring},
  {"layer", ops_addLayer, METH_VARARGS, layer_docstring},
  {"geomTransf", ops_addGeomTransf, METH_VARARGS, geomTransf_docstring},
  {"beamIntegrationRule", ops_addBeamIntegrationRule, METH_VARARGS, beamIntegration_docstring},
  {"mass", ops_addNodalMass, METH_VARARGS, mass_docstring},
  {"groundMotion", ops_addGroundMotion, METH_VARARGS, groundMotion_docstring},
  {"eleLoad", ops_addElementalLoad, METH_VARARGS, eleLoad_docstring},
  // analysis methods
  {"wipeAnalysis", ops_wipeAnalysis, METH_VARARGS, wipeAnalysis_docstring},
  {"wipeReliability", ops_wipeReliability, METH_VARARGS, wipeReliability_docstring},
  {"wipe", ops_wipeModel, METH_VARARGS, wipe_docstring},
  {"Print", ops_printModel, METH_VARARGS, Print_docstring},
  {"system", ops_specifySOE, METH_VARARGS, system_docstring},
  {"numberer", ops_specifyNumberer, METH_VARARGS, numberer_docstring},
  {"constraints", ops_specifyConstraintHandler, METH_VARARGS, constraints_docstring},
  {"algorithm", ops_specifyAlgorithm, METH_VARARGS, algorithm_docstring},
  {"integrator", ops_specifyIntegrator, METH_VARARGS, integrator_docstring},
  {"analysis", ops_specifyAnalysis, METH_VARARGS, analysis_docstring},
  {"test", ops_specifyCTest, METH_VARARGS, test_docstring},
  {"analyze", ops_analyzeModel, METH_VARARGS, analyze_docstring},
  {"nodeDisp", ops_nodeDisp, METH_VARARGS, nodeDisp_docstring},
  {"getLoadFactor", ops_getLoadFactor, METH_VARARGS, getLoadFactor_docstring},
  {"loadConst", ops_setLoadConst, METH_VARARGS, loadConst_docstring},
  {"rayleigh", ops_rayleighDamping, METH_VARARGS, rayleigh_docstring},
  {"eigen", ops_eigenAnalysis, METH_VARARGS, eigen_docstring},
  {NULL, NULL, 0, NULL}
};

#if PY_MAJOR_VERSION==2
void addOpenSeesCommands(void) {
  PyImport_AddModule("opensees");
  Py_InitModule("opensees", methodsOpenSees);
}

PyMODINIT_FUNC
initopensees(void)
{
    (void) Py_InitModule("opensees", methodsOpenSees);
}

#elif PY_MAJOR_VERSION==3
static struct PyModuleDef opsmodule = {
   PyModuleDef_HEAD_INIT,
   "opensees",   /* name of module */
   NULL, /* module documentation, may be NULL */
   -1,       /* size of per-interpreter state of the module,
                or -1 if the module keeps state in global variables. */
   methodsOpenSees
};


PyMODINIT_FUNC
PyInit_opensees(void)
{
    return PyModule_Create(&opsmodule);
}

void addOpenSeesCommands(void) {
  PyModule_Create(&opsmodule);
}

#endif


/* ********************************* *
 * elementAPI FUNCTIONS for Python   *
 *********************************** */

static int numberArgs = 0;
static int currentArg = 0;
PyObject *currentArgs = 0;

extern "C"
void OPS_ResetCommandLine(int nArgs, int cArg, PyObject *args) {
  numberArgs = nArgs;
  currentArg = cArg;
  currentArgs = args;
}


extern "C"   
int OPS_GetNumRemainingInputArgs()
{
  return numberArgs - currentArg;
}

extern "C"   
int OPS_GetIntInput(int *numData, int*data)
{
  if ((numberArgs - currentArg) >= *numData) {
    for (int i=0; i<*numData; i++) {
      PyObject *o = PyTuple_GetItem(currentArgs,currentArg);
      
#if PY_MAJOR_VERSION==2
      if (!PyInt_Check(o)) {
	  //PyErr_SetString(PyExc_ValueError,"invalid integer");
	  return -1;
      }
#elif PY_MAJOR_VERSION==3
      if (!PyLong_Check(o)) {
	  //PyErr_SetString(PyExc_ValueError,"invalid integer");
	  return -1;
      }
#endif   

      
#if PY_MAJOR_VERSION==2
      long value = PyInt_AS_LONG(o);
#elif PY_MAJOR_VERSION==3
      long value = PyLong_AsLong(o);
#endif
      
      data[i] = value;
      currentArg++;
    }
  
    return 0;
  }
  return -1;
}


extern "C" 
int OPS_GetDoubleInput(int *numData, double *data)
{
  int numD = *numData;
  if ((numberArgs - currentArg) >= numD) {
    for (int i=0; i<numD; i++) {
      PyObject *o = PyTuple_GetItem(currentArgs,currentArg);
      if (!PyFloat_Check(o)) {
	  
#if PY_MAJOR_VERSION==2
	if (!PyInt_Check(o)) {
	    //PyErr_SetString(PyExc_ValueError,"invalid double");
	    return -1;
	}
#elif PY_MAJOR_VERSION==3
	if (!PyLong_Check(o)) {
	    //PyErr_SetString(PyExc_ValueError,"invalid double");
	    return -1;
	}
#endif
	
#if PY_MAJOR_VERSION==2
	data[i] = PyInt_AS_LONG(o);
#elif PY_MAJOR_VERSION==3
	data[i] = PyLong_AsLong(o);
#endif
	currentArg++;
	continue;
      }
      data[i] = PyFloat_AS_DOUBLE(o);
      currentArg++;
    }
    return 0;  
  }

  return -1;
}

extern "C" 
const char *OPS_GetString(void)
{
  PyObject *o = PyTuple_GetItem(currentArgs,currentArg);
  currentArg++;
#if PY_MAJOR_VERSION==2
  return PyString_AS_STRING(o);
#elif PY_MAJOR_VERSION==3
  return PyUnicode_AS_DATA(o);
#endif
}

extern "C" int OPS_GetNDM()
{
    PythonModelBuilder& pBuilder = OPS_GetPythonModelBuilder();
    if(pBuilder.isSet()) return pBuilder.getNDM();
    return 0;
}

extern "C" int OPS_GetNDF()
{
    PythonModelBuilder& pBuilder = OPS_GetPythonModelBuilder();
    if(pBuilder.isSet()) return pBuilder.getNDF();
    return 0;
}

extern UniaxialMaterial *
OPS_GetUniaxialMaterial(int matTag) {
    return OPS_getUniaxialMaterial(matTag);
}

extern Domain*
OPS_GetDomain()
{
    return &theDomain;
}
