import opensees as ops
area = 548.4 * 1e-6
gamma = 46.12 / area
E = 1.31 * 1e11
l = 20
h = 2
x = [-10,0,10]
y = [0,0,0]
z = [h,0,h]
alfa = 6.5e-6
cambiodetemp = 0.
w3 = -gamma * area
rho = 0
L0 = 10.2878
sag = .6
def Opensees_model():	
    ops.wipe()
    ops.model('basic', '-ndm', 3, '-ndf', 3)
    nnode = 0
    for inode in range(3):
        nnode += 1
        ops.node(nnode, x[inode], y[inode], z[inode])
        ops.fix(nnode, 1, 1, 1)

    errorTol = 1e-6
    NSubSteps = 30
    num_ele = 0
    for nele in range(1, nnode, 1):
        num_ele += 1
        if(num_ele ==1):
            ops.element('CatenaryCable', num_ele, nele, nele + 1, w3, E, area, sag, alfa, cambiodetemp, rho,errorTol, NSubSteps, 0, "-maxsag")
        else:
            ops.element('CatenaryCable', num_ele, nele, nele + 1, w3, E, area,  L0, alfa, cambiodetemp, rho,errorTol, NSubSteps, 0)
        
    ops.timeSeries('Constant', 1)
    ops.pattern('Plain', 1, 1)

    NSteps = 20
    ops.recorder('Node', '-file', "disp.txt", '-time', '-nodeRange', 1, nnode, '-dof', 1, 2, 3, 'disp')
    ops.recorder('Element', '-file', "forces.txt", '-time', '-eleRange', 1, num_ele, 'force')

    ops.system('UmfPack')  # ('FullGeneral')
    ops.constraints('Plain')
    ops.numberer('Plain')
    ops.test('NormDispIncr', 1.0e-6, 200, 0)
    ops.integrator('LoadControl', 1.0 / NSteps)
    ops.algorithm('RaphsonNewton')
    ops.analysis('Static')
    ops.analyze(NSteps)
#
    Node_list = ops.getNodeTags()
    Node_Rection = []
    Node_Coord = []
    for inode in Node_list:
        Node_Coord.append(ops.nodeCoord(inode))
    
    ops.reactions()
    Node_Rection.append(ops.nodeReaction(Node_list[0]))
    Node_Rection.append(ops.nodeReaction(Node_list[-1]))
    return Node_list, Node_Coord, Node_Rection
    
if __name__ == "__main__":
    Node_list, Node_Coord, Node_Rection = Opensees_model()  
    print(Node_list)
    print(Node_Coord)
    print(Node_Rection)    
