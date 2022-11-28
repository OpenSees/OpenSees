import openseespy.opensees as ops

def DiscretizeMember(ndI,ndJ,numEle,eleType,integrTag,transfTag,nodeTag,eleTag):

    nodeList = []
    eleList = []

    # A lazy workaround in case someone uses elasticBeamColumn, which
    # has different input order than the material nonlinear elements
    if eleType == 'elasticBeamColumn':
        eleType = 'elasticForceBeamColumn'
    
    if numEle <= 1:
        ops.element(eleType,eleTag,ndI,ndJ,transfTag,integrTag)
        eleList.append(eleTag)
        return eleList,nodeList

    Xi = ops.nodeCoord(ndI,'X')
    Yi = ops.nodeCoord(ndI,'Y')
    Xj = ops.nodeCoord(ndJ,'X')
    Yj = ops.nodeCoord(ndJ,'Y')

    dX = (Xj-Xi)/numEle
    dY = (Yj-Yi)/numEle

    threeD = True
    if len(ops.nodeCoord(ndI)) < 3:
        threeD = False
    else:
        Zi = ops.nodeCoord(ndI,'Z')
        Zj = ops.nodeCoord(ndJ,'Z')
        dZ = (Zj-Zi)/numEle

    nodes = [None]*(numEle+1)
    nodes[0] = ndI
    nodes[numEle] = ndJ

    for i in range(1,numEle):
        if threeD:
            ops.node(nodeTag,Xi+i*dX,Yi+i*dY,Zi+i*dZ)
        else:
            ops.node(nodeTag,Xi+i*dX,Yi+i*dY)
        nodes[i] = nodeTag
        nodeList.append(nodeTag)
        nodeTag = nodeTag+1

    #print(eleType,eleTag,ndI,nodes[1],transfTag,integrTag)
    ops.element(eleType,eleTag,ndI,nodes[1],transfTag,integrTag)
    eleList.append(eleTag)
    eleTag = eleTag + 1

    for i in range(1,numEle-1):
        ops.element(eleType,eleTag,nodes[i],nodes[i+1],transfTag,integrTag)
        eleList.append(eleTag)
        eleTag = eleTag + 1

    ops.element(eleType,eleTag,nodes[numEle-1],ndJ,transfTag,integrTag)
    eleList.append(eleTag)
    eleTag = eleTag + 1

    return eleList,nodeList
