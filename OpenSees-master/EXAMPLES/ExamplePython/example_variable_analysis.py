import sys
sys.path.insert(0, '../../SRC/interpreter/')
# sys.path.insert(0, '../build/lib/')
import opensees as ops


ops.model('basic', '-ndm', 1, '-ndf', 1)
ops.uniaxialMaterial('Elastic', 1, 3000.0)

ops.node(1, 0.0)
ops.node(2, 72.0)

ops.fix(1, 1)

ops.element('Truss', 1, 1, 2, 10.0, 1)
ops.timeSeries('Linear', 1)
ops.pattern('Plain', 1, 1)
ops.load(2, 100.0)

ops.constraints('Transformation')
ops.numberer('ParallelPlain')
ops.test('NormDispIncr', 1e-6, 6, 2)
ops.system('ProfileSPD')
ops.integrator('Newmark', 0.5, 0.25)
# ops.analysis('Transient')
ops.algorithm('Linear')
ops.analysis('VariableTransient')


ops.analyze(5, 0.0001, 0.00001, 0.001, 10)
time = ops.getTime()
print(f'time: ', ops.getTime())
approx_vtime = 0.0001 + 0.001  # One step at target, then one step at maximum
assert 0.99 < time / approx_vtime < 1.01,  (time,  approx_vtime)
ops.setTime(0.0)
# Can still run a non-variable analysis - since analyze function has multiple dispatch.
ops.analyze(5, 0.0001)
time = ops.getTime()
print(f'time: ', ops.getTime())
approx_vtime = 0.0001 * 5  # variable transient is not active so time should be dt * 5
# If variable transient is not active then time would be 0.0005
assert 0.99 < time / approx_vtime < 1.01,  (time,  approx_vtime)
