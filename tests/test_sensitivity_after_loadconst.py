try:
   import opensees as ops
except ModuleNotFoundError:
   import openseespy.opensees as ops
import math
import pytest

# DDM sensitivity of a transient analysis that follows a static gravity stage.
#
# The usual workflow applies gravity in a static analysis, holds it with
# loadConst, and then runs a transient analysis. loadConst freezes the load
# factor of the gravity pattern at its value at that time. The load
# sensitivity of a frozen pattern is therefore
#
#     P * dlambda_frozen/dh + lambda_frozen * dP/dh,
#
# where dlambda_frozen/dh is the factor sensitivity of the time series at the
# time the pattern was frozen. It is zero for any parameter that is not a
# parameter of the time series (a material parameter, a nodal load), and
# nonzero for a factor of the time series itself.
#
# The test compares DDM with central finite differences for a fiber-section
# column (Concrete02 and Steel01) under ground acceleration. Gravity is
#   'static'    applied in a static LoadControl stage, held with loadConst
#   'transient' ramped in a first transient stage, held with loadConst
#   'constant'  a constant load inside the transient analysis (reference)
# for three kinds of parameter: material parameters, the magnitude of the
# gravity load, and the factor of the time series of the gravity pattern.
# The gravity pattern uses a series factor of 2 and half the load, so that
# the frozen factor (2) differs from its sensitivity to the series factor (1).
#
# The static stage is run with material parameters only: during a static
# LoadControl stage the sensitivity to load parameters (the load magnitude
# and the series factor) is not propagated consistently with the load factor
# of the step, which is a separate matter from loadConst.

NOMINAL = {'fc': -30.0, 'epsc0': -0.002, 'fcu': -6.0, 'epscu': -0.012,
           'rat': 0.1, 'ft': 3.0, 'Ets': 1500.0}
GRAVITY = -200.0e3     # times the series factor
SERIES_FACTOR = 2.0
DT = 0.01
NUM_STEPS = 150
RAMP_STEPS = 50

# (kind, name, nominal value)
PARAMETERS = [('material', 'fc', NOMINAL['fc']),
              ('material', 'Ets', NOMINAL['Ets']),
              ('load', 'gravity', GRAVITY),
              ('series', 'factor', SERIES_FACTOR)]
MODES = ['static', 'transient', 'constant']
CASES = [(mode, *p) for mode in MODES for p in PARAMETERS
         if not (mode == 'static' and p[0] != 'material')]


def run(mode, kind=None, name=None, perturbation=0.0, sensitivity=False):
   values = dict(NOMINAL)
   gravity = GRAVITY
   factor = SERIES_FACTOR
   if kind == 'material':
      values[name] += perturbation
   elif kind == 'load':
      gravity += perturbation
   elif kind == 'series':
      factor += perturbation

   ops.wipe()
   ops.model('basic', '-ndm', 2, '-ndf', 3)
   ops.node(1, 0.0, 0.0); ops.fix(1, 1, 1, 1)
   ops.node(2, 0.0, 1500.0)
   ops.uniaxialMaterial('Concrete02', 1, *values.values())
   ops.uniaxialMaterial('Steel01', 2, 420.0, 200000.0, 0.01)
   ops.section('Fiber', 1)
   ops.patch('rect', 1, 20, 1, -150.0, -150.0, 150.0, 150.0)
   ops.layer('straight', 2, 2, 300.0, 110.0, -100.0, 110.0, 100.0)
   ops.layer('straight', 2, 2, 300.0, -110.0, -100.0, -110.0, 100.0)
   ops.geomTransf('Linear', 1)
   ops.beamIntegration('Lobatto', 1, 1, 5)
   ops.element('dispBeamColumn', 1, 1, 2, 1, 1)
   ops.mass(2, 40.0, 40.0, 1.0e-9)

   ops.timeSeries('Constant' if mode == 'constant' else 'Linear', 1, '-factor', factor)
   ops.pattern('Plain', 1, 1)
   ops.load(2, 0.0, gravity, 0.0)

   if sensitivity:
      if kind == 'material':
         ops.parameter(1, 'element', 1, 'material', 1, name)
      elif kind == 'load':
         ops.parameter(1, 'loadPattern', 1, 'loadAtNode', 2, '2')
      elif kind == 'series':
         ops.parameter(1, 'loadPattern', 1, 'timeSeries', 'factor')

   def common():
      ops.constraints('Plain'); ops.numberer('Plain'); ops.system('BandGeneral')
      ops.test('NormDispIncr', 1.0e-12, 100); ops.algorithm('Newton')

   ops.rayleigh(2.1, 0.0, 0.0, 0.0)
   if mode == 'static':
      common()
      ops.integrator('LoadControl', 0.1)
      ops.analysis('Static')
      if sensitivity:
         ops.sensitivityAlgorithm('-computeAtEachStep')
      for _ in range(10):
         assert ops.analyze(1) == 0
      ops.loadConst('-time', 0.0)
      ops.wipeAnalysis()
   elif mode == 'transient':
      # gravity ramped over pseudo time 0 to 0.5 in a transient analysis
      common()
      ops.integrator('Newmark', 0.5, 0.25)
      ops.analysis('Transient')
      if sensitivity:
         ops.sensitivityAlgorithm('-computeAtEachStep')
      for _ in range(RAMP_STEPS):
         assert ops.analyze(1, DT) == 0
      ops.loadConst('-time', 0.0)
      ops.wipeAnalysis()

   accel = [9810.0 * 0.5 * math.sin(2.0 * math.pi * 1.5 * k * DT) * min(1.0, k * DT)
            for k in range(NUM_STEPS + 10)]
   ops.timeSeries('Path', 2, '-dt', DT, '-values', *accel)
   ops.pattern('UniformExcitation', 2, 1, '-accel', 2)
   common()
   ops.integrator('Newmark', 0.5, 0.25)
   ops.analysis('Transient')
   if sensitivity:
      ops.sensitivityAlgorithm('-computeAtEachStep')

   U, S = [], []
   for step in range(NUM_STEPS):
      assert ops.analyze(1, DT) == 0, f'analysis failed at step {step}'
      U.append(ops.nodeDisp(2, 1))
      S.append(ops.sensNodeDisp(2, 1, 1) if sensitivity else 0.0)
   return U, S


def relative_errors(mode, kind, name, nominal, h_rel=1.0e-6):
   _, ddm = run(mode, kind, name, sensitivity=True)
   h = h_rel * abs(nominal)
   Up, _ = run(mode, kind, name, h)
   Um, _ = run(mode, kind, name, -h)
   fd = [(a - b) / (2.0 * h) for a, b in zip(Up, Um)]
   scale = max(abs(x) for x in fd)
   # largest error relative to the peak, and largest error relative to the
   # local value wherever the local value exceeds 1% of the peak
   peak = max(abs(a - b) for a, b in zip(ddm, fd)) / scale
   local = max((abs(a - b) / abs(b) for a, b in zip(ddm, fd) if abs(b) > 0.01 * scale), default=0.0)
   return peak, local


@pytest.mark.parametrize('mode,kind,name,nominal', CASES, ids=[f'{c[0]}-{c[2]}' for c in CASES])
def test_transient_ddm_matches_finite_difference(mode, kind, name, nominal):
   peak, local = relative_errors(mode, kind, name, nominal)
   assert peak < 1.0e-6 and local < 1.0e-4, \
      f'{name}: relative difference {peak:.2e} of the peak, {local:.2e} locally'
