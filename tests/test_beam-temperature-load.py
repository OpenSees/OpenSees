try:
   import opensees as ops
except ModuleNotFoundError:
   import openseespy.opensees as ops
from math import isclose
import pytest

# Beam2dTempLoad on the 2d elastic beam-columns.  The temperature is linear
# through the depth d and along the length, so the free thermal strain is an
# axial part alpha*Tc, Tc the mid-depth temperature, and a linear curvature
# k0 = -alpha*(Ttop - Tbot)/d.  Reference values are closed forms of the beam
# equations under that imposed strain, except where two models must agree.

E = 2.0e8
A = 0.02
I = 3.0e-4
L = 6.0

alpha = 1.2e-5
d = 0.4

G = 8.0e7
Av = 8.0e-4
phi = 12*E*I/(G*Av*L**2)

# (Ttop, Tbot) at node I, then at node J
FIELDS = [(25.0,5.0,25.0,5.0),
          (10.0,10.0,40.0,0.0),
          (30.0,-10.0,0.0,20.0)]
FIELD_IDS = ['uniform','rising','reversing']

# the largest gradient above sets the scale of the absolute tolerances
DISP = alpha*40.0/d*L**2
FORCE = E*I*alpha*40.0/d

def near(a,b,scale):
   return isclose(a,b,rel_tol=1e-9,abs_tol=1e-9*scale)

def curvature(Ttop,Tbot):
   # a hotter top bends the member concave down
   return -alpha*(Ttop - Tbot)/d

def end_curvatures(T):
   return curvature(T[0],T[1]),curvature(T[2],T[3])

def mean_temperature(T):
   # mean of the mid-depth temperature along the length
   return sum(T)/4

def analyze():
   ops.analysis('Static','-noWarnings')
   ops.analyze(1)
   ops.reactions()

def elastic_beam(T,fixI,fixJ,release=0,shear=False):
   ops.wipe()
   ops.model('basic','-ndm',2,'-ndf',3)

   ops.node(1,0.0,0.0); ops.fix(1,*fixI)
   ops.node(2,L,0.0)
   if any(fixJ):
      ops.fix(2,*fixJ)

   ops.geomTransf('Linear',1)
   props = [A,E,I,G,Av] if shear else [A,E,I]
   ops.element('elasticBeamColumn',1,1,2,*props,1,
               '-alpha',alpha,'-depth',d,'-release',release)

   ops.timeSeries('Constant',1)
   ops.pattern('Plain',1,1)
   ops.eleLoad('-ele',1,'-type','-beamTemp',*T)

   analyze()

def mod_elastic_beam(T,fixJ,K=(4.0,4.0,2.0),depth=True):
   ops.wipe()
   ops.model('basic','-ndm',2,'-ndf',3)

   ops.node(1,0.0,0.0); ops.fix(1,1,1,1)
   ops.node(2,L,0.0)
   if any(fixJ):
      ops.fix(2,*fixJ)

   ops.geomTransf('Linear',1)
   options = ['-alpha',alpha] + (['-d',d] if depth else [])
   ops.element('ModElasticBeam2d',1,1,2,A,E,I,*K,1,*options)

   ops.timeSeries('Constant',1)
   ops.pattern('Plain',1,1)
   ops.eleLoad('-ele',1,'-type','-beamTemp',*T)

   analyze()

def reactions():
   return [ops.nodeReaction(1,2),ops.nodeReaction(1,3),
           ops.nodeReaction(2,2),ops.nodeReaction(2,3)]

def clamped_reactions(T,p):
   # Clamped against bending.  For k0 = c0 + c1*x/L the end rotation and end
   # deflection conditions of a member with shear flexibility phi give
   #   MI =  EI*(c0 + c1*phi/(2*(1 + phi)))
   #   MJ = -EI*(c0 + c1*(2 + phi)/(2*(1 + phi)))
   # and phi = 0 leaves the Euler-Bernoulli moments EI*k0(0) and -EI*k0(L).
   kI,kJ = end_curvatures(T)
   c0 = kI
   c1 = kJ - kI
   MI = E*I*(c0 + c1*p/(2*(1 + p)))
   MJ = -E*I*(c0 + c1*(2 + p)/(2*(1 + p)))
   V = (MI + MJ)/L
   return [V,MI,-V,MJ]

@pytest.mark.parametrize('shear',[False,True])
@pytest.mark.parametrize('T',FIELDS,ids=FIELD_IDS)
def test_cantilever(T,shear):
   # Statically determinate, so the tip follows the free thermal strain for
   # any shear flexibility: the tip rotation is the integral of k0 and the
   # tip deflection its first moment about the tip.
   elastic_beam(T,(1,1,1),(0,0,0),shear=shear)
   kI,kJ = end_curvatures(T)

   assert near(ops.nodeDisp(2,1),alpha*L*mean_temperature(T),DISP)
   assert near(ops.nodeDisp(2,2),L**2*(2*kI + kJ)/6,DISP)
   assert near(ops.nodeDisp(2,3),L*(kI + kJ)/2,DISP)

@pytest.mark.parametrize('shear',[False,True])
@pytest.mark.parametrize('T',FIELDS,ids=FIELD_IDS)
def test_simply_supported(T,shear):
   # Determinate as well: w'' = k0 with w(0) = w(L) = 0 fixes the end
   # rotations, and no support reacts to a load that applies no force.
   elastic_beam(T,(1,1,0),(0,1,0),shear=shear)
   kI,kJ = end_curvatures(T)

   assert near(ops.nodeDisp(1,3),-L*(2*kI + kJ)/6,DISP)
   assert near(ops.nodeDisp(2,3),L*(kI + 2*kJ)/6,DISP)
   assert near(ops.nodeReaction(1,2),0.0,FORCE)
   assert near(ops.nodeReaction(2,2),0.0,FORCE)

@pytest.mark.parametrize('shear',[False,True])
@pytest.mark.parametrize('T',FIELDS,ids=FIELD_IDS)
def test_clamped(T,shear):
   elastic_beam(T,(1,1,1),(0,1,1),shear=shear)
   expected = clamped_reactions(T,phi if shear else 0.0)

   for got,ref in zip(reactions(),expected):
      assert near(got,ref,FORCE)
   assert near(ops.nodeDisp(2,1),alpha*L*mean_temperature(T),DISP)

@pytest.mark.parametrize('release',[1,2,3])
@pytest.mark.parametrize('shear',[False,True])
@pytest.mark.parametrize('T',FIELDS,ids=FIELD_IDS)
def test_release(T,shear,release):
   # A released end is a hinge.  Between clamped nodes the released element
   # must carry what the unreleased one carries with that end pinned instead,
   # which leaves the condensation to the global solution.
   pinI = release in (1,3)
   pinJ = release in (2,3)
   elastic_beam(T,(1,1,0 if pinI else 1),(0,1,0 if pinJ else 1),shear=shear)
   expected = reactions()

   elastic_beam(T,(1,1,1),(0,1,1),release=release,shear=shear)
   for got,ref in zip(reactions(),expected):
      assert near(got,ref,FORCE)

@pytest.mark.parametrize('T',FIELDS,ids=FIELD_IDS)
def test_fully_clamped(T):
   # Two elements clamped at both ends, the temperature continuous at the
   # middle node.  The member stays straight, M = -EI*k0 throughout, and the
   # middle node shifts by alpha*l*(TcI - TcJ)/4 under the compression
   # EA*alpha*mean(Tc), l = L/2 being the element length.
   ops.wipe()
   ops.model('basic','-ndm',2,'-ndf',3)

   ops.node(1,0.0,0.0); ops.fix(1,1,1,1)
   ops.node(2,L/2,0.0)
   ops.node(3,L,0.0); ops.fix(3,1,1,1)

   ops.geomTransf('Linear',1)
   ops.element('elasticBeamColumn',1,1,2,A,E,I,1,'-alpha',alpha,'-depth',d)
   ops.element('elasticBeamColumn',2,2,3,A,E,I,1,'-alpha',alpha,'-depth',d)

   Tm = [(T[0] + T[2])/2,(T[1] + T[3])/2]
   ops.timeSeries('Constant',1)
   ops.pattern('Plain',1,1)
   ops.eleLoad('-ele',1,'-type','-beamTemp',T[0],T[1],*Tm)
   ops.eleLoad('-ele',2,'-type','-beamTemp',*Tm,T[2],T[3])

   analyze()
   kI,kJ = end_curvatures(T)
   TcI = (T[0] + T[1])/2
   TcJ = (T[2] + T[3])/2
   MI = E*I*kI
   MJ = -E*I*kJ

   assert near(ops.nodeDisp(2,1),alpha*(L/2)*(TcI - TcJ)/4,DISP)
   assert near(ops.nodeDisp(2,2),0.0,DISP)
   assert near(ops.nodeDisp(2,3),0.0,DISP)
   assert near(ops.nodeReaction(1,1),E*A*alpha*mean_temperature(T),FORCE)
   assert near(ops.nodeReaction(1,2),(MI + MJ)/L,FORCE)
   assert near(ops.nodeReaction(1,3),MI,FORCE)
   assert near(ops.nodeReaction(3,3),MJ,FORCE)

@pytest.mark.parametrize('T',FIELDS,ids=FIELD_IDS)
def test_mod_cantilever(T):
   # Standard stiffness coefficients reduce ModElasticBeam2d to the elastic
   # beam, so the tip follows the free thermal strain.
   mod_elastic_beam(T,(0,0,0))
   kI,kJ = end_curvatures(T)

   assert near(ops.nodeDisp(2,1),alpha*L*mean_temperature(T),DISP)
   assert near(ops.nodeDisp(2,2),L**2*(2*kI + kJ)/6,DISP)
   assert near(ops.nodeDisp(2,3),L*(kI + kJ)/2,DISP)

@pytest.mark.parametrize('T',FIELDS,ids=FIELD_IDS)
def test_mod_clamped(T):
   # The clamped element carries the standard fixed-end forces whatever its
   # stiffness modifiers, as for its other element loads.
   mod_elastic_beam(T,(0,1,1),K=(4.4,4.4,2.2))

   for got,ref in zip(reactions(),clamped_reactions(T,0.0)):
      assert near(got,ref,FORCE)

@pytest.mark.parametrize('element',['elasticBeamColumn','ModElasticBeam2d'])
def test_without_depth(element):
   # With no depth given there is no gradient to act on, and a uniform
   # temperature change only lengthens the member.
   T = [15.0]*4
   if element == 'elasticBeamColumn':
      ops.wipe()
      ops.model('basic','-ndm',2,'-ndf',3)
      ops.node(1,0.0,0.0); ops.fix(1,1,1,1)
      ops.node(2,L,0.0)
      ops.geomTransf('Linear',1)
      ops.element('elasticBeamColumn',1,1,2,A,E,I,1,'-alpha',alpha)
      ops.timeSeries('Constant',1)
      ops.pattern('Plain',1,1)
      ops.eleLoad('-ele',1,'-type','-beamTemp',*T)
      analyze()
   else:
      mod_elastic_beam(T,(0,0,0),depth=False)

   assert near(ops.nodeDisp(2,1),alpha*L*15.0,DISP)
   assert near(ops.nodeDisp(2,2),0.0,DISP)
   assert near(ops.nodeDisp(2,3),0.0,DISP)
