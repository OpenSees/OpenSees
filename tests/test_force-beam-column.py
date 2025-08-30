try:
   import opensees as ops
except ModuleNotFoundError:
   import openseespy.opensees as ops
from math import isclose

L = 48
E = 29000
A = 20
I = 800
P = 10

def end_loaded_cantilever():
   ops.wipe()
   ops.model('basic','-ndm',2,'-ndf',3)

   ops.node(1,0,0); ops.fix(1,1,1,1)
   ops.node(2,L,0)

   ops.section('Elastic',1,E,A,I)
   ops.beamIntegration('Lobatto',1,1,3)

   ops.geomTransf('Linear',1)

   ops.element('forceBeamColumn',1,1,2,1,1)

   ops.timeSeries('Constant',1)
   ops.pattern('Plain',1,1)
   ops.load(2,0,P,0)

   ops.analysis('Static','-noWarnings')
   ops.analyze(1)
   ops.reactions()

def test_deflection():
   end_loaded_cantilever()
   assert isclose(ops.nodeDisp(2,2),P*L**3/(3*E*I))

def test_rotation():
   end_loaded_cantilever()
   assert isclose(ops.nodeDisp(2,3),P*L**2/(2*E*I))

def test_force():
   end_loaded_cantilever()
   assert isclose(ops.nodeReaction(1,2),-P)

def test_moment():
   end_loaded_cantilever()
   assert isclose(ops.nodeReaction(1,3),-P*L)
