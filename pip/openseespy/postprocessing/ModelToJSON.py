import json

def WriteModelToJSON(filename,indent=-1,outformat='.6g'):
    
    modelData = {'ndm': ops.getNDM()[0], 'bounds': ops.nodeBounds()}
    NDF = ops.getNDF()[0]
    
    nodes = ops.getNodeTags()
    Nnodes = len(nodes)
    fixedNodes = ops.getFixedNodes()
    
    nodeData = {}
    for nd in nodes:
        nodeData[nd] = {'crds': [float(f"{item:{outformat}}") for item in ops.nodeCoord(nd)]}
        if nd in fixedNodes:
            fixedDOFs = ops.getFixedDOFs(nd)
            fixities = [0]*NDF
            for dof in fixedDOFs:
                fixities[dof-1] = -1;
            nodeData[nd]['fix'] = fixities
                
    modelData['nodes'] = nodeData
        
    eles = ops.getEleTags()
    Neles = len(eles)

    eleTypes = {}
    for ele in eles:
        eleTypes[ops.getEleClassTags(ele)[0]] = ops.eleType(ele)
    modelData['elementTypes'] = eleTypes
    
    modelData['sectionTypes'] = {1: 'rect'}
    modelData['sections'] = {}
    
    # A unique burner tag
    #
    paramTag = 1
    paramTags = ops.getParamTags()
    if len(paramTags) > 0:
        paramTag = max(paramTags) + 1
    
    for ele in eles:
        if ops.eleType(ele) in ['ElasticBeam3d']:
            ops.parameter(paramTag,'element',ele,'A')
            A = ops.getParamValue(paramTag)
            ops.remove('parameter',paramTag)
            ops.parameter(paramTag,'element',ele,'Iz')
            Iz = ops.getParamValue(paramTag)
            ops.remove('parameter',paramTag)
            ops.parameter(paramTag,'element',ele,'Iy')
            Iy = ops.getParamValue(paramTag)
            ops.remove('parameter',paramTag)            
            h = (12*Iz/A)**0.5
            b = (12*Iy/A)**0.5
            h = float(f"{h:{outformat}}")
            b = float(f"{b:{outformat}}")
            modelData['sections'][ele] = {'type': 1, 'data': [b,h]}
            
    elementData = {}
    for ele in eles:
        elementData[ele] = {'nodes': ops.eleNodes(ele), 'type': ops.getEleClassTags(ele)[0]}
        
        sectionData = {}
        eleType = ops.eleType(ele)
        if eleType in ['ElasticBeam2d','ElasticBeam3d']:
            sectionData['yaxis'] = [float(f"{item:{outformat}}") for item in ops.eleResponse(ele,'yaxis')]
            sectionData['zaxis'] = [float(f"{item:{outformat}}") for item in ops.eleResponse(ele,'zaxis')]
            sectionData['tag'] = ele
            elementData[ele]['section'] = sectionData
        if eleType in ['ZeroLength','TwoNodeLink']:
            sectionData['yaxis'] = [float(f"{item:{outformat}}") for item in ops.eleResponse(ele,'yaxis')]
            sectionData['zaxis'] = [float(f"{item:{outformat}}") for item in ops.eleResponse(ele,'zaxis')]
            sectionData['xaxis'] = [float(f"{item:{outformat}}") for item in ops.eleResponse(ele,'xaxis')]
            sectionData['mats'] = [int(item) for item in ops.eleResponse(ele,'materials')]
            sectionData['dirs'] = [int(item) for item in ops.eleResponse(ele,'directions')]
            elementData[ele]['section'] = sectionData
    modelData['elements'] = elementData

    patternData = {}
    patterns = ops.getPatterns()
    Npatterns = len(patterns)
    for pattern in patterns:

        patternData[pattern] = {}

        patternData[pattern]['nodalLoad'] = {}
        nodalLoads = ops.getNodeLoadTags(pattern)
        nodalLoadData = ops.getNodeLoadData(pattern)
        loc = 0
        for nd in nodalLoads:
            patternData[pattern]['nodalLoad'][nd] = [float(f"{item:{outformat}}") for item in nodalLoadData[loc:loc+NDF]]
            loc += NDF

        patternData[pattern]['eleLoad'] = {}
        eleLoads = ops.getEleLoadTags(pattern)
        eleLoadClassTags = ops.getEleLoadClassTags(pattern)
        eleLoadData = ops.getEleLoadData(pattern)
        loc = 0
        for iele,ele in enumerate(eleLoads):
            if eleLoadClassTags[iele] == 5: # Beam3DUniformLoad
                Ndata = 3
                wy,wz,wx = eleLoadData[loc:loc+Ndata]
                patternData[pattern]['eleLoad'][ele] = [float(f"{item:{outformat}}") for item in [wx,wy,wz]]
            loc += Ndata
    modelData['patterns'] = patternData
    
    outjson = open(filename,'w')
    if indent < 0:
        outjson.write(json.dumps(modelData,separators=(',', ':')) + '\n')
    else:
        outjson.write(json.dumps(modelData,indent=indent,separators=(',', ':')) + '\n')
    outjson.close()
