try:
   import opensees as ops
except ModuleNotFoundError:
   import openseespy.opensees as ops
import pytest

# DDM response sensitivity of Concrete02IS with respect to its eight material
# parameters, verified against central finite differences in two models: a
# material-level model with trusses and a cantilever column with a fiber
# section.
#
# Concrete02IS is Concrete02 with a user-specified initial stiffness E0 and a
# Popovics curve before the peak in compression,
#
#    stress = fc*eta*r/(r - 1 + eta^r),  eta = strain/epsc0,
#    r = E0/(E0 - fc/epsc0),
#
# with r = 400 when fc/epsc0 is not smaller than E0. The nominal E0 below
# differs from 2*fc/epsc0, so the exponent r is not 2.
#
# Truss model. Two parallel trusses of unit length and unit area connect a
# fixed node to a free node: one is Concrete02IS, the other is elastic and
# stiffer than the concrete, so that load control with plain Newton converges
# through the whole cyclic history. The history visits every branch of the
# model: the Popovics curve, the linear descending branch, the residual
# plateau, unloading and reloading in compression, and cracking, softening
# and reloading in tension.
#
# Plain Newton is used on purpose: the sensitivity solve reuses the tangent
# left in the system of equations by the solution algorithm, so the algorithm
# must form the current tangent at every iteration.

# Concrete02IS is piecewise smooth. At a corner of the envelopes or at a
# change of branch the response has only one-sided derivatives, and DDM
# returns the derivative of the branch that is active. The load histories
# below are chosen so that no converged state lies on such a point, where a
# centered finite difference would straddle two branches.

# Internal sign convention of Concrete02IS: compression values are negative
NOMINAL = {
   'E0': 25000.0,
   'fc': -30.0,
   'epsc0': -0.002,
   'fcu': -6.0,
   'epscu': -0.012,
   'rat': 0.1,
   'ft': 3.0,
   'Ets': 1500.0,
}
NAMES = list(NOMINAL)

E_ELASTIC = 60000.0

# Load reversal points and the number of equal load steps between them
LOAD_PEAKS = [-179.0, 18.5, -381.0, 61.0, -553.0, 49.0, -847.0, 25.0]
STEPS_PER_LEG = 40

# Relative perturbation for the central finite differences
FD_STEP = 1.0e-6
# Allowed difference between DDM and finite difference at each step:
# RELATIVE_TOLERANCE times the local finite difference, plus a floor of
# ABSOLUTE_FLOOR times the largest finite difference in the history. The
# floor only matters where the sensitivity passes through zero
RELATIVE_TOLERANCE = 1.0e-4
ABSOLUTE_FLOOR = 1.0e-6


def assert_sensitivities_match(parameter, ddm, fd):
   scale = max(abs(value) for value in fd)
   assert scale > 0.0

   for step, (d, f) in enumerate(zip(ddm, fd)):
      allowed = RELATIVE_TOLERANCE*abs(f) + ABSOLUTE_FLOOR*scale
      assert abs(d - f) <= allowed, \
         f'{parameter}: step {step}, DDM {d}, finite difference {f}'


def define_material(values):
   ops.uniaxialMaterial('Concrete02IS', 1, *[values[name] for name in NAMES])


def build_model(values, peaks=LOAD_PEAKS):
   ops.wipe()
   ops.model('basic', '-ndm', 1, '-ndf', 1)

   ops.node(1, 0.0); ops.fix(1, 1)
   ops.node(2, 1.0)

   define_material(values)
   ops.uniaxialMaterial('Elastic', 2, E_ELASTIC)

   ops.element('truss', 1, 1, 2, 1.0, 1)
   ops.element('truss', 2, 1, 2, 1.0, 2)

   # Piecewise linear load history through the reversal points. The last
   # value is repeated so that round-off in the accumulated time cannot
   # take the final step past the end of the series
   times = list(range(len(peaks) + 2))
   loads = [0.0] + list(peaks) + [peaks[-1]]
   ops.timeSeries('Path', 1, '-time', *times, '-values', *loads)
   ops.pattern('Plain', 1, 1)
   ops.load(2, 1.0)


def define_analysis(compute_sensitivity):
   ops.constraints('Plain')
   ops.numberer('Plain')
   ops.system('BandGeneral')
   ops.test('NormUnbalance', 1.0e-11, 100)
   ops.algorithm('Newton')
   ops.integrator('LoadControl', 1.0/STEPS_PER_LEG)
   ops.analysis('Static')

   if compute_sensitivity:
      ops.sensitivityAlgorithm('-computeAtEachStep')


def advance(num_steps, num_parameters=0):
   # Runs num_steps load steps. Returns the displacement, the stress in the
   # Concrete02IS truss and the DDM displacement sensitivities (one value
   # per registered parameter) after each step
   displacements = []
   stresses = []
   sensitivities = []
   for step in range(num_steps):
      assert ops.analyze(1) == 0, f'analysis failed at step {step}'
      displacements.append(ops.nodeDisp(2, 1))
      stresses.append(ops.basicForce(1)[0])
      sensitivities.append([ops.sensNodeDisp(2, 1, tag)
                            for tag in range(1, num_parameters + 1)])

   return displacements, stresses, sensitivities


def run_history(values, parameter=None, updated_value=None, num_steps=None,
                peaks=LOAD_PEAKS):
   # Returns the displacement history and, when a parameter name is given,
   # the DDM displacement sensitivity history for that parameter. When
   # updated_value is given, the parameter is set to it with updateParameter
   # before the analysis and no sensitivity is computed
   build_model(values, peaks)

   if parameter is not None:
      ops.parameter(1, 'element', 1, parameter)
   if updated_value is not None:
      ops.updateParameter(1, updated_value)

   compute_sensitivity = parameter is not None and updated_value is None
   define_analysis(compute_sensitivity)

   if num_steps is None:
      num_steps = len(peaks)*STEPS_PER_LEG

   displacements, _, sensitivities = advance(
      num_steps, 1 if compute_sensitivity else 0)

   return displacements, [values[0] for values in sensitivities if values]


def finite_difference(values, parameter, peaks=LOAD_PEAKS):
   step = FD_STEP*abs(values[parameter])

   plus = dict(values); plus[parameter] += step
   minus = dict(values); minus[parameter] -= step

   u_plus, _ = run_history(plus, peaks=peaks)
   u_minus, _ = run_history(minus, peaks=peaks)

   return [(up - um)/(2.0*step) for up, um in zip(u_plus, u_minus)]


def test_truss_history_reaches_all_branches():
   # Guards the tests below: the nominal history must reach the residual
   # plateau in compression, soften in tension and crack completely,
   # otherwise some branches of the model would go unchecked
   build_model(NOMINAL)
   define_analysis(False)
   u, stress, _ = advance(len(LOAD_PEAKS)*STEPS_PER_LEG)

   fc = NOMINAL['fc']
   fcu = NOMINAL['fcu']
   steps = range(1, len(u))

   # Past the peak stress and onto the residual plateau
   assert min(u) < NOMINAL['epscu']
   assert min(stress) < 0.99*fc
   assert any(abs(s - fcu) <= 1.0e-12*abs(fcu) for s in stress)
   # Unloading in compression: stress and strain both increase
   assert any(u[k] > u[k-1] and stress[k-1] < stress[k] < 0.0 for k in steps)
   # Tension softening: stress decreases while the strain increases
   assert any(u[k] > u[k-1] and 0.0 < stress[k] < stress[k-1] for k in steps)
   # Completely cracked after having carried tension. The residual stress
   # after cracking is exactly zero, as in Concrete02
   assert any(stress[k] == 0.0 and max(stress[:k]) > 0.0 for k in steps)


@pytest.mark.parametrize('parameter', NAMES)
def test_truss_ddm_matches_finite_difference(parameter):
   _, ddm = run_history(NOMINAL, parameter)
   fd = finite_difference(NOMINAL, parameter)

   assert_sensitivities_match(parameter, ddm, fd)


# When fc/epsc0 is not smaller than E0 the Popovics exponent is fixed at
# r = 400 and no longer depends on fc, epsc0 and E0
CAPPED_EXPONENT = dict(NOMINAL, E0=14000.0)


@pytest.mark.parametrize('parameter', ['E0', 'fc', 'epsc0'])
def test_capped_popovics_exponent_ddm_matches_finite_difference(parameter):
   assert CAPPED_EXPONENT['fc']/CAPPED_EXPONENT['epsc0'] >= \
      CAPPED_EXPONENT['E0']

   _, ddm = run_history(CAPPED_EXPONENT, parameter)
   fd = finite_difference(CAPPED_EXPONENT, parameter)

   assert_sensitivities_match(parameter, ddm, fd)


def test_all_parameters_together_match_single_parameter_runs():
   # With all eight parameters defined in the same analysis, each one must
   # give the sensitivity it gives when it is the only parameter. This
   # checks that the sensitivity history of one parameter is not mixed
   # with that of another
   build_model(NOMINAL)
   for tag, name in enumerate(NAMES, 1):
      ops.parameter(tag, 'element', 1, name)
   define_analysis(True)
   _, _, together = advance(len(LOAD_PEAKS)*STEPS_PER_LEG, len(NAMES))

   for index, name in enumerate(NAMES):
      _, alone = run_history(NOMINAL, name)
      scale = max(abs(value) for value in alone)
      assert scale > 0.0
      for step, (values, expected) in enumerate(zip(together, alone)):
         assert abs(values[index] - expected) <= 1.0e-12*scale, \
            f'{name}: step {step}'


# Load history with three legs of constant load: on the compression envelope
# before the peak, after partial unloading, and on the softening branch in
# tension. The strain does not change in those steps, so the sensitivity
# must stay constant as well
HOLD_PEAKS = [-134.0, -134.0, -60.0, -60.0, -170.0, 33.0, 33.0]


@pytest.mark.parametrize('parameter', NAMES)
def test_constant_load_steps_match_finite_difference(parameter):
   _, ddm = run_history(NOMINAL, parameter, peaks=HOLD_PEAKS)
   fd = finite_difference(NOMINAL, parameter, peaks=HOLD_PEAKS)

   assert_sensitivities_match(parameter, ddm, fd)


ALIASES = [
   ('fpc', 'fc'),
   ('epsco', 'epsc0'),
   ('fpcu', 'fcu'),
   ('epsU', 'epscu'),
   ('epsu', 'epscu'),
   ('lambda', 'rat'),
]


@pytest.mark.parametrize('alias, parameter', ALIASES)
def test_parameter_alias(alias, parameter):
   # Alternative parameter names address the same material property
   _, expected = run_history(NOMINAL, parameter, num_steps=2*STEPS_PER_LEG)
   _, ddm = run_history(NOMINAL, alias, num_steps=2*STEPS_PER_LEG)

   assert max(abs(value) for value in expected) > 0.0
   assert ddm == expected


def test_reset_clears_sensitivity_history():
   # After reset (revertToStart) the same load history must give the same
   # sensitivities again, which requires the sensitivity history kept by
   # the material to be cleared together with its state
   num_steps = 3*STEPS_PER_LEG

   build_model(NOMINAL)
   ops.parameter(1, 'element', 1, 'E0')
   define_analysis(True)
   u_first, _, ddm_first = advance(num_steps, 1)

   ops.reset()
   assert ops.nodeDisp(2, 1) == 0.0

   u_second, _, ddm_second = advance(num_steps, 1)

   assert max(abs(values[0]) for values in ddm_first) > 0.0
   assert u_second == u_first
   assert ddm_second == ddm_first


# updateParameter. A parameter changed with updateParameter before the
# analysis must give the same response as a material created with the new
# value. Two load legs are enough to reach unloading from the compression
# envelope and cracking in tension.
UPDATE_FACTOR = 1.1
UPDATE_STEPS = 2*STEPS_PER_LEG


@pytest.mark.parametrize('parameter', NAMES)
def test_update_parameter_matches_new_material(parameter):
   updated = dict(NOMINAL); updated[parameter] *= UPDATE_FACTOR
   expected, _ = run_history(updated, num_steps=UPDATE_STEPS)
   nominal, _ = run_history(NOMINAL, num_steps=UPDATE_STEPS)

   u, _ = run_history(NOMINAL, parameter, updated_value=updated[parameter],
                      num_steps=UPDATE_STEPS)

   # The update must change the response, and must reproduce the response
   # of the material created with the new value
   assert max(abs(a - b) for a, b in zip(expected, nominal)) > 0.0
   for step, (a, b) in enumerate(zip(u, expected)):
      assert a == pytest.approx(b, rel=1.0e-12, abs=1.0e-15), \
         f'{parameter}: step {step}'


@pytest.mark.parametrize('parameter', ['fc', 'epsc0', 'fcu', 'epscu'])
def test_update_parameter_keeps_compression_negative(parameter):
   # Compression values are stored as negative numbers. A positive value
   # passed to updateParameter is given the same sign as in the constructor
   updated = dict(NOMINAL); updated[parameter] *= UPDATE_FACTOR
   expected, _ = run_history(updated, num_steps=UPDATE_STEPS)

   u, _ = run_history(NOMINAL, parameter,
                      updated_value=abs(updated[parameter]),
                      num_steps=UPDATE_STEPS)

   for step, (a, b) in enumerate(zip(u, expected)):
      assert a == pytest.approx(b, rel=1.0e-12, abs=1.0e-15), \
         f'{parameter}: step {step}'


# Column model. A cantilever column modeled with one displacement-based
# beam-column element and a fiber section: Concrete02IS fibers over the depth
# and two layers of elastic steel. A constant axial compression is combined
# with a cyclic lateral load at the top. At the base section the extreme
# concrete fibers pass the peak stress in compression and crack, soften and
# reload in tension, with different fibers on different branches at the same
# time. The parameter is attached to every Concrete02IS fiber of every
# section.
#
# Units: N, mm, MPa
COLUMN_LENGTH = 1500.0
SECTION_DEPTH = 300.0
SECTION_WIDTH = 300.0
NUM_CONCRETE_FIBERS = 20
STEEL_LAYER_DEPTH = 110.0
STEEL_BAR_AREA = 600.0
E_STEEL = 200000.0
AXIAL_LOAD = -300.0e3

# Lateral load reversal points and the number of load steps between them
LATERAL_PEAKS = [100.0e3, -140.0e3, 180.0e3, -220.0e3, 150.0e3]
LATERAL_STEPS_PER_LEG = 40


def build_column(values):
   ops.wipe()
   ops.model('basic', '-ndm', 2, '-ndf', 3)

   ops.node(1, 0.0, 0.0); ops.fix(1, 1, 1, 1)
   ops.node(2, 0.0, COLUMN_LENGTH)

   define_material(values)
   ops.uniaxialMaterial('Elastic', 2, E_STEEL)

   h = 0.5*SECTION_DEPTH
   b = 0.5*SECTION_WIDTH
   ys = STEEL_LAYER_DEPTH
   ops.section('Fiber', 1)
   ops.patch('rect', 1, NUM_CONCRETE_FIBERS, 1, -h, -b, h, b)
   ops.layer('straight', 2, 2, STEEL_BAR_AREA, ys, -100.0, ys, 100.0)
   ops.layer('straight', 2, 2, STEEL_BAR_AREA, -ys, -100.0, -ys, 100.0)

   ops.geomTransf('Linear', 1)
   ops.beamIntegration('Lobatto', 1, 1, 5)
   ops.element('dispBeamColumn', 1, 1, 2, 1, 1)

   # Constant axial compression
   ops.timeSeries('Constant', 1)
   ops.pattern('Plain', 1, 1)
   ops.load(2, 0.0, AXIAL_LOAD, 0.0)

   # Cyclic lateral load; the last value is repeated for the same reason
   # as in the truss model
   times = list(range(len(LATERAL_PEAKS) + 2))
   loads = [0.0] + LATERAL_PEAKS + [LATERAL_PEAKS[-1]]
   ops.timeSeries('Path', 2, '-time', *times, '-values', *loads)
   ops.pattern('Plain', 2, 2)
   ops.load(2, 1.0, 0.0, 0.0)


def run_column_history(values, parameter=None):
   # Returns the lateral displacement history at the top of the column, the
   # DDM sensitivity history when a parameter name is given, and the lowest
   # and highest extreme-fiber strains reached at the base section
   build_column(values)

   if parameter is not None:
      ops.parameter(1, 'element', 1, 'material', 1, parameter)

   ops.constraints('Plain')
   ops.numberer('Plain')
   ops.system('BandGeneral')
   # The displacement increment test makes the last tangent formed by
   # Newton the tangent of the converged state, which the sensitivity
   # solve then reuses
   ops.test('NormDispIncr', 1.0e-11, 50)
   ops.algorithm('Newton')
   ops.integrator('LoadControl', 1.0/LATERAL_STEPS_PER_LEG)
   ops.analysis('Static')

   if parameter is not None:
      ops.sensitivityAlgorithm('-computeAtEachStep')

   displacements = []
   sensitivities = []
   min_strain = 0.0
   max_strain = 0.0
   for step in range(len(LATERAL_PEAKS)*LATERAL_STEPS_PER_LEG):
      assert ops.analyze(1) == 0, f'analysis failed at step {step}'
      displacements.append(ops.nodeDisp(2, 1))
      if parameter is not None:
         sensitivities.append(ops.sensNodeDisp(2, 1, 1))

      axial, curvature = ops.eleResponse(1, 'section', 1, 'deformation')
      edge = 0.5*SECTION_DEPTH*abs(curvature)
      min_strain = min(min_strain, axial - edge)
      max_strain = max(max_strain, axial + edge)

   return displacements, sensitivities, (min_strain, max_strain)


@pytest.mark.parametrize('parameter', NAMES)
def test_fiber_section_ddm_matches_finite_difference(parameter):
   u, ddm, (min_strain, max_strain) = run_column_history(NOMINAL, parameter)

   # The extreme fibers must pass the peak stress in compression and the
   # end of tension softening, otherwise some branches would go unchecked
   assert min_strain < NOMINAL['epsc0']
   assert max_strain > NOMINAL['ft']*(1.0/NOMINAL['Ets'] +
                                      1.0/NOMINAL['E0'])

   step = FD_STEP*abs(NOMINAL[parameter])
   plus = dict(NOMINAL); plus[parameter] += step
   minus = dict(NOMINAL); minus[parameter] -= step
   u_plus, _, _ = run_column_history(plus)
   u_minus, _, _ = run_column_history(minus)
   fd = [(up - um)/(2.0*step) for up, um in zip(u_plus, u_minus)]

   assert_sensitivities_match(parameter, ddm, fd)
