import opensees as ops

# units: kip, in

# remove existing model
ops.wipe()

# set modelbuilder
ops.model('basic', '-ndm', 2, '-ndf', 2)

# create nodes
ops.node(1, 0.0, 0.0, "-disp", 0.0, 0.0, "-vel", 0.0,0.0, "-mass", 0.0,0.0)
ops.node(2, 144.0,  0.0)
ops.node(3, 168.0,  0.0)
ops.node(4,  72.0, 96.0)

# set boundary condition
ops.fix(1, 1, 1)
ops.fix(2, 1, 1)
ops.fix(3, 1, 1)

# define materials
ops.uniaxialMaterial("Elastic", 1, 3000)

# define elements
ops.element("Truss",1,1,4,10.0,1)
ops.element("Truss",2,2,4,5.0,1)
ops.element("Truss",3,3,4,5.0,1)

# create TimeSeries
ops.timeSeries("Linear", 1)

# create a plain load pattern
ops.pattern("Plain", 1, 1, "-fact", 1.0)
ops.load(4, 100, -50)

# print model
#ops.Print()

# create SOE
ops.system("BandSPD")

# create DOF number
ops.numberer("RCM")

# create constraint handler
ops.constraints("Plain")

# create algorithm
ops.algorithm("Linear")

# create integrator
ops.integrator("LoadControl", 1.0)

# create analysis object
ops.analysis("Static")

# perform the analysis
ops.analyze(1)

# print results
print "node 4 displacement: ", ops.nodeDisp(4)
ops.Print('node',4)
ops.Print('ele')



