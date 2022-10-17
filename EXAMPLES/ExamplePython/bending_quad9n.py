import openseespy.opensees as ops
# import opensees as ops

ops.wipe()

ops.model('basic', '-ndm', 2, '-ndf', 2)

L = 5.
H = 1.

thk = 0.01

P = 100.
E = 200.e6
nu = 0.3

ops.nDMaterial('ElasticIsotropic', 1, E, nu)

ops.node(1, 0., 0.)
ops.node(2, 5., 0.)
ops.node(3, 5., 1.)
ops.node(4, 0., 1.)
ops.node(5, 2.5, 0.)
ops.node(6, 5., .5)
ops.node(7, 2.5, 1.)
ops.node(8, 0., .5)
ops.node(9, 2.5, .5)  # comment for quad8n element

ops.element('quad9n', 1, 1, 2, 3, 4, 5, 6, 7, 8, 9, thk, 'PlaneStress', 1)
# ops.element('quad8n', 1, 1, 2, 3, 4, 5, 6, 7, 8, thk, 'PlaneStress', 1)

ops.fix(1, 1, 1)
ops.fix(4, 1, 0)
ops.fix(8, 1, 0)

ops.timeSeries('Linear', 1)
ops.pattern('Plain', 1, 1)
ops.load(2, P, 0.)
ops.load(3, -P, 0.)

ops.analysis('Static')

ops.analyze(1)

ops.printModel()

# verification:
# tip vertical displacement (node 2 and 3) = 0.0075
# bottom Gauss Point stress_xx = 46475.8
# bottom extreme stress_xx (extrapolated) = 60000.0

exit()
