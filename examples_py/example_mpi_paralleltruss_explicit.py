# import custom_openseespy.opensees as ops
import opensees as ops

pid = ops.getPID()
print('pid: ', pid)
np = ops.getNP()
print('np: ', np)
ops.start()
a = open('nps.txt', 'w')
a.write(str(np))
a.close()
if np < 2:
    exit()


ops.model('basic', '-ndm', 2, '-ndf', 2)
ops.uniaxialMaterial('Elastic', 1, 3000.0)

if pid == 0:
    ops.node(1, 0.0, 0.0)
    ops.node(4, 72.0, 96.0)

    ops.fix(1, 1, 1)
    ops.mass(4, 100.0, 100.0)

    ops.element('Truss', 1, 1, 4, 10.0, 1)
    ops.timeSeries('Linear', 1)
    ops.pattern('Plain', 1, 1)
    ops.load(4, 100.0, -50.0)

else:
    ops.node(2, 144.0, 0.0)
    ops.node(3, 168.0, 0.0)
    ops.node(4, 72.0, 96.0)

    ops.fix(2, 1, 1)
    ops.fix(3, 1, 1)
    ops.mass(4, 100.0, 100.0)

    ops.element('Truss', 2, 2, 4, 5.0, 1)
    ops.element('Truss', 3, 3, 4, 5.0, 1)

ops.constraints('Transformation')
ops.numberer('ParallelPlain')
ops.test('NormDispIncr', 1e-6, 6, 2)
ops.algorithm('Linear')
etype = 'central_difference'
etype = 'explicit_difference'  # Comment out this line to run with central difference
if etype == 'central_difference':
    ops.system('Mumps')
    ops.integrator('CentralDifference')
else:
    # ops.system('Mumps')
    ops.system('MPIDiagonal')  # Can use Mumps here but not sure if it scales as well
    ops.integrator('ExplicitDifference')
ops.analysis('Transient')
for i in range(30):
    print(f'######################################## run {i} ##')
    ops.analyze(1, 0.000001)
print('PPP')
ops.analyze(20, 0.00001)

print(pid, ' Node 4: ', [ops.nodeCoord(4), ops.nodeDisp(4)])
print(pid, " COMPLETED")

