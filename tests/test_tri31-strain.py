try:
   import opensees as ops
except ModuleNotFoundError:
   import openseespy.opensees as ops
from math import isclose

# A single Tri31 in plane stress. The element is a constant strain triangle, so
# a uniform stress state gives strains that are exact, not approximate.

E = 1000.0
nu = 0.25
G = E/(2*(1+nu))

th = 1.5
a = 2.0
b = 3.0

sigma = 5.0
tau = 4.0

def uniaxial_tension():
   # Supports remove the rigid body modes only, and the applied node force is
   # statically equivalent to sigma11 = sigma with sigma22 = sigma12 = 0.
   ops.wipe()
   ops.model('basic','-ndm',2,'-ndf',2)

   ops.node(1,0.0,0.0); ops.fix(1,1,1)
   ops.node(2,a,0.0)
   ops.node(3,0.0,b); ops.fix(3,1,0)

   ops.nDMaterial('ElasticIsotropic',1,E,nu)
   ops.element('tri31',1,1,2,3,th,'PlaneStress',1)

   ops.timeSeries('Constant',1)
   ops.pattern('Plain',1,1)
   ops.load(2,sigma*b*th/2,0.0)

   ops.analysis('Static','-noWarnings')
   ops.analyze(1)

def pure_shear():
   # Node force statically equivalent to sigma12 = tau with sigma11 = sigma22 = 0.
   ops.wipe()
   ops.model('basic','-ndm',2,'-ndf',2)

   ops.node(1,0.0,0.0); ops.fix(1,1,1)
   ops.node(2,a,0.0); ops.fix(2,0,1)
   ops.node(3,0.0,b)

   ops.nDMaterial('ElasticIsotropic',1,E,nu)
   ops.element('tri31',1,1,2,3,th,'PlaneStress',1)

   ops.timeSeries('Constant',1)
   ops.pattern('Plain',1,1)
   ops.load(3,tau*a*th/2,0.0)

   ops.analysis('Static','-noWarnings')
   ops.analyze(1)

def test_axial_strain():
   uniaxial_tension()
   strains = ops.eleResponse(1,'strains')

   assert isclose(strains[0],sigma/E)
   assert isclose(strains[1],-nu*sigma/E)
   assert isclose(strains[2],0.0,abs_tol=1e-12)

def test_shear_strain():
   pure_shear()
   strains = ops.eleResponse(1,'strains')

   assert isclose(strains[0],0.0,abs_tol=1e-12)
   assert isclose(strains[1],0.0,abs_tol=1e-12)
   assert isclose(strains[2],tau/G)

def test_strain_alias():
   uniaxial_tension()
   strains = ops.eleResponse(1,'strains')

   assert len(strains) == 3
   assert ops.eleResponse(1,'strain') == strains
