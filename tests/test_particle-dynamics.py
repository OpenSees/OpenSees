try:
   import opensees as ops
except ModuleNotFoundError:
   import openseespy.opensees as ops
from math import isclose 

m = 1.0
F = 1.0
t = 1.0

a = F/m

def define_model():
   ops.wipe()
   ops.model('basic','-ndm',1,'-ndf',1)

   ops.node(1,0); ops.mass(1,m)

   ops.timeSeries('Linear',1)
   ops.pattern('Plain',1,1)
   ops.load(1,F)

   ops.integrator('Newmark',0.5,1.0/6)
   ops.analysis('Transient','-noWarnings')
   Nsteps = 10
   dt = t/Nsteps
   ops.analyze(Nsteps,dt)

def test_acceleration():
   define_model()

   assert isclose(ops.nodeAccel(1,1),a)

def test_velocity():
   define_model()

   assert isclose(ops.nodeVel(1,1),0.5*a*t)

def test_displacement():
   define_model()

   assert isclose(ops.nodeDisp(1,1),a*t*t/6)

if __name__ == '__main__':
   test_displacement()
