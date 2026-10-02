try:
   import opensees as ops
except ModuleNotFoundError:
   import openseespy.opensees as ops
from math import isclose
import pytest

# A unit cube of mass density rho under a body force b. Displacements are
# fixed laterally everywhere and vertically at the base, and pressures are
# fixed everywhere, so a static analysis returns the element's body force as
# reactions summing to -rho*V*b. The three components differ, so the order
# in which they are read is checked too.

rho = 2.0
b = (1.5, -2.0, -9.81)

# arguments after the eight nodes, up to the optional body force
ARGS = {
   'brickUP': (1, 2.0e4, 1.0, 1.0e-6, 1.0e-6, 1.0e-6),
   'bbarBrickUP': (1, 2.0e4, 1.0, 1.0e-6, 1.0e-6, 1.0e-6),
   'SSPbrickUP': (1, 2.0e4, 1.0, 1.0e-6, 1.0e-6, 1.0e-6, 0.5, 0.0),
}

FORMS = {
   'numbers': list(b),
   'numbers then -lumped': list(b) + ['-lumped'],
   '-lumped then numbers': ['-lumped'] + list(b),
   'numeric strings': [str(v) for v in b],
}

def total_reaction(element, options):
   ops.wipe()
   ops.model('basic','-ndm',3,'-ndf',4)

   corners = [(0,0,0),(1,0,0),(1,1,0),(0,1,0),(0,0,1),(1,0,1),(1,1,1),(0,1,1)]
   for tag,(x,y,z) in enumerate(corners,1):
      ops.node(tag,float(x),float(y),float(z))
      ops.fix(tag,1,1,1 if z == 0 else 0,1)

   ops.nDMaterial('ElasticIsotropic',1,1.0e4,0.3,rho)
   ops.element(element,1,*range(1,9),*ARGS[element],*options)

   ops.timeSeries('Constant',1)
   ops.pattern('Plain',1,1)

   ops.analysis('Static','-noWarnings')
   ops.analyze(1)
   ops.reactions()

   return [sum(ops.nodeReaction(tag,dof) for tag in range(1,9)) for dof in (1,2,3)]

@pytest.mark.parametrize('form',FORMS)
@pytest.mark.parametrize('element',ARGS)
def test_body_force(element,form):
   for got,bi in zip(total_reaction(element,FORMS[form]),b):
      assert isclose(got,-rho*bi,rel_tol=1e-12)
