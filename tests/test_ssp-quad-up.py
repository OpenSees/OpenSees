try:
   import opensees as ops
except ModuleNotFoundError:
   import openseespy.opensees as ops
from math import isclose, pi, sin, exp
import pytest

# SSPquadUP with a massless solid (rho = 0), as in a quasi-static
# consolidation analysis. The pore fluid of bulk modulus Kf fills the
# porosity n = e/(1 + e), so the storage coefficient is S = n/Kf; Kf is low
# enough here for S to change the answer.

E = 1.0e4
nu = 0.3
e = 2.0/3.0
n = e/(1 + e)
Kf = 2.0e4
perm = 1.0e-6
q = 10.0

M = E*(1 - nu)/((1 + nu)*(1 - 2*nu))  # constrained modulus
u0 = q/(1 + n*M/Kf)                    # undrained pore pressure

def transient(constraints='Plain'):
   ops.constraints(constraints)
   ops.numberer('RCM')
   ops.system('BandGeneral')
   ops.test('NormDispIncr',1.0e-12,10,0)
   ops.algorithm('Linear')
   ops.integrator('Newmark',0.5,0.25)
   ops.analysis('Transient')

def test_undrained_pressure():
   # One element, laterally confined and sealed, loaded on top. With Newmark
   # gamma = 1/2 the storage balance S*p + div(u) = 0 holds exactly at every
   # step, so the pressure is the undrained q/(1 + n*M/Kf) for any time step.
   # The pressure DOFs are tied, as in an undrained element test: the element's
   # one-point permeability and stabilization leave the pressure hourglass
   # mode free when no DOF is drained.
   ops.wipe()
   ops.model('basic','-ndm',2,'-ndf',3)

   ops.node(1,0.0,0.0); ops.fix(1,1,1,0)
   ops.node(2,1.0,0.0); ops.fix(2,1,1,0)
   ops.node(3,1.0,1.0); ops.fix(3,1,0,0)
   ops.node(4,0.0,1.0); ops.fix(4,1,0,0)
   for node in (2,3,4):
      ops.equalDOF(1,node,3)

   ops.nDMaterial('ElasticIsotropic',1,E,nu,0.0)
   ops.element('SSPquadUP',1,1,2,3,4,1,1.0,Kf,1.0,perm,perm,e,1.0e-6)

   ops.timeSeries('Constant',1)
   ops.pattern('Plain',1,1)
   ops.load(3,0.0,-q/2,0.0)
   ops.load(4,0.0,-q/2,0.0)

   transient('Transformation')
   for step in range(3):
      assert ops.analyze(1,0.1) == 0
      # pore pressure is the velocity of the pressure DOF
      for node in (1,2,3,4):
         assert isclose(ops.nodeVel(node,3),u0,rel_tol=1e-9)
      assert isclose(ops.nodeDisp(3,2),-n/Kf*u0,rel_tol=1e-9)

def test_consolidation():
   # A column drained at the top and loaded at t = 0 follows Terzaghi's
   # solution with the diffusivity c = perm/(n/Kf + 1/M) of a compressible
   # pore fluid, starting from the undrained pressure u0.
   H, w, ny = 1.0, 0.1, 20
   c = perm/(n/Kf + 1/M)

   ops.wipe()
   ops.model('basic','-ndm',2,'-ndf',3)
   for j in range(ny + 1):
      for k,x in enumerate((0.0,w)):
         ops.node(2*j + k + 1,x,j*H/ny)
         ops.fix(2*j + k + 1,1,1 if j == 0 else 0,1 if j == ny else 0)

   ops.nDMaterial('ElasticIsotropic',1,E,nu,0.0)
   for j in range(ny):
      ops.element('SSPquadUP',j + 1,2*j + 1,2*j + 2,2*j + 4,2*j + 3,1,1.0,Kf,1.0,perm,perm,e,0.0)

   ops.timeSeries('Constant',1)
   ops.pattern('Plain',1,1)
   ops.load(2*ny + 1,0.0,-q*w/2,0.0)
   ops.load(2*ny + 2,0.0,-q*w/2,0.0)

   transient()
   T = 0.1
   assert ops.analyze(50,T/50*H*H/c) == 0

   modes = [pi*(2*m + 1)/2 for m in range(200)]
   def pressure(z):
      return u0*sum(2/Mm*sin(Mm*z/H)*exp(-Mm*Mm*T) for Mm in modes)
   U = 1 - sum(2/(Mm*Mm)*exp(-Mm*Mm*T) for Mm in modes)
   settlement = (q - u0)*H/M + u0*H/M*U

   for j in range(ny):
      assert abs(ops.nodeVel(2*j + 1,3) - pressure(H - j*H/ny)) < 1.0e-2*u0
   assert abs(-ops.nodeDisp(2*ny + 1,2) - settlement) < 5.0e-3*q*H/M

# (Pup, Plow, Pleft, Pright), the nodes of the loaded side, its outward normal
SIDES = {
   'lower': ((0.0,10.0,0.0,0.0),(1,2),(0.0,-1.0)),
   'upper': ((10.0,0.0,0.0,0.0),(3,4),(0.0,1.0)),
   'left': ((0.0,0.0,10.0,0.0),(4,1),(-1.0,0.0)),
   'right': ((0.0,0.0,0.0,10.0),(2,3),(1.0,0.0)),
}

@pytest.mark.parametrize('side',SIDES)
def test_side_pressure(side):
   # A pressure P on one side of a unit square acts along the outward normal
   # and loads each node of that side with P/2, on its displacement DOFs only.
   # Every DOF but one, off the loaded side, is fixed, so the reactions
   # return that load with the opposite sign.
   P,loaded,normal = SIDES[side]

   ops.wipe()
   ops.model('basic','-ndm',2,'-ndf',3)
   for node,(x,y) in enumerate(((0,0),(1,0),(1,1),(0,1)),1):
      ops.node(node,float(x),float(y))
      ops.fix(node,1,1,1)
   free = next(node for node in (1,2,3,4) if node not in loaded)
   ops.remove('sp',free,1)

   ops.nDMaterial('ElasticIsotropic',1,E,nu,0.0)
   ops.element('SSPquadUP',1,1,2,3,4,1,1.0,Kf,1.0,perm,perm,e,0.0,0.0,0.0,*P)

   ops.timeSeries('Constant',1)
   ops.pattern('Plain',1,1)

   ops.constraints('Plain')
   ops.numberer('Plain')
   ops.system('FullGeneral')
   ops.test('NormDispIncr',1.0e-12,10,0)
   ops.algorithm('Linear')
   ops.integrator('LoadControl',1.0)
   ops.analysis('Static')
   assert ops.analyze(1) == 0
   ops.reactions()

   for node in (1,2,3,4):
      share = 5.0 if node in loaded else 0.0
      expected = (-share*normal[0],-share*normal[1],0.0)
      for dof in (1,2,3):
         assert isclose(ops.nodeReaction(node,dof),expected[dof - 1],abs_tol=1e-12)
