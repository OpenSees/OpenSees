#!/usr/bin/env python
import h5py
import numpy as np
import os, sys

class Xdmf:
  def __init__(self, filename):
    self.filename = filename
    self.cellMap  = {1 : {1 : 'Polyvertex', 2 : 'Polyline'}, 2 : {3 : 'Triangle', 4 : 'Quadrilateral'}, 3 : {4 : 'Tetrahedron', 6: 'Wedge', 8 : 'Hexahedron'}}
    self.typeMap  = {'scalar' : 'Scalar', 'vector' : 'Vector', 'tensor' : 'Tensor6', 'matrix' : 'Matrix'}
    self.typeExt  = {2 : {'vector' : ['x', 'y'], 'tensor' : ['xx', 'yy', 'xy']}, 3 : {'vector' : ['x', 'y', 'z'], 'tensor' : ['xx', 'yy', 'zz', 'xy', 'yz', 'xz']}}
    return

  def writeHeader(self, fp, hdfFilename):
    fp.write('''\
<?xml version="1.0" ?>
<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" [
<!ENTITY HeavyData "%s">
]>
''' % os.path.basename(hdfFilename))
    fp.write('\n<Xdmf>\n  <Domain Name="domain">\n')
    return

  def writeCells(self, fp, topologyPath, numCells, numCorners):
    fp.write('''\
    <DataItem Name="cells"
	      ItemType="Uniform"
	      Format="HDF"
	      NumberType="Float" Precision="8"
	      Dimensions="%d %d">
      &HeavyData;:/%s/cells
    </DataItem>
''' % (numCells, numCorners, topologyPath))
    return

  def writeVertices(self, fp, geometryPath, numVertices, spaceDim):
    fp.write('''\
    <DataItem Name="vertices"
	      Format="HDF"
	      Dimensions="%d %d">
      &HeavyData;:/%s/vertices
    </DataItem>
    <!-- ============================================================ -->
''' % (numVertices, spaceDim, geometryPath))
    return

  def writeLocations(self, fp, numParticles, spaceDim):
    fp.write('''\
    <DataItem Name="xcoord"
	      Format="HDF"
	      Dimensions="%d">
      &HeavyData;:/particles/xcoord
    </DataItem>
''' % (numParticles))
    if spaceDim == 1: return
    fp.write('''\
    <DataItem Name="ycoord"
	      Format="HDF"
	      Dimensions="%d">
      &HeavyData;:/particles/ycoord
    </DataItem>
''' % (numParticles))
    if spaceDim == 2: return
    fp.write('''\
    <DataItem Name="zcoord"
	      Format="HDF"
	      Dimensions="%d">
      &HeavyData;:/particles/zcoord
    </DataItem>
''' % (numParticles))
    return

  def writeTimeGridHeader(self, fp, time):
    fp.write('''\
    <Grid Name="TimeSeries" GridType="Collection" CollectionType="Temporal">
      <Time TimeType="List">
        <DataItem Format="XML" NumberType="Float" Dimensions="%d">
          ''' % (len(time)))
    fp.write(' '.join(map(str, map(float, time))))
    fp.write('''
        </DataItem>
      </Time>
''')
    return

  def writeSpaceGridHeader(self, fp, numCells, numCorners, cellDim, spaceDim):
    fp.write('''\
      <Grid Name="domain" GridType="Uniform">
	<Topology
	   TopologyType="%s"
	   NumberOfElements="%d">
	  <DataItem Reference="XML">
	    /Xdmf/Domain/DataItem[@Name="cells"]
	  </DataItem>
	</Topology>
	<Geometry GeometryType="%s">
	  <DataItem Reference="XML">
	    /Xdmf/Domain/DataItem[@Name="vertices"]
	  </DataItem>
	</Geometry>
''' % (self.cellMap[cellDim][numCorners], numCells, "XYZ" if spaceDim > 2 else "XY"))
    return

  def writeFieldSingle(self, fp, numSteps, timestep, spaceDim, name, f, domain):
    if len(f[1].shape) > 2:
      dof = f[1].shape[1]
      bs  = f[1].shape[2]
    elif len(f[1].shape) > 1:
      if numSteps > 1:
        dof = f[1].shape[1]
        bs  = 1
      else:
        dof = f[1].shape[0]
        bs  = f[1].shape[1]
    else:
      dof = f[1].shape[0]
      bs  = 1
    fp.write('''\
	<Attribute
	   Name="%s"
	   Type="%s"
	   Center="%s">
          <DataItem ItemType="HyperSlab"
		    Dimensions="1 %d %d"
		    Type="HyperSlab">
            <DataItem
	       Dimensions="3 3"
	       Format="XML">
              %d 0 0
              1 1 1
              1 %d %d
	    </DataItem>
	    <DataItem
	       DataType="Float" Precision="8"
	       Dimensions="%d %d %d"
	       Format="HDF">
	      &HeavyData;:%s
	    </DataItem>
	  </DataItem>
	</Attribute>
''' % (f[0], self.typeMap[f[1].attrs['vector_field_type']], domain, dof, bs, timestep, dof, bs, numSteps, dof, bs, name))
    return

  def writeFieldComponents(self, fp, numSteps, timestep, spaceDim, name, f, domain):
    vtype = f[1].attrs['vector_field_type']
    if len(f[1].shape) > 2:
      dof    = f[1].shape[1]
      bs     = f[1].shape[2]
      cdims  = '1 %d 1' % dof
      dims   = '%d %d %d' % (numSteps, dof, bs)
      stride = '1 1 1'
      size   = '1 %d 1' % dof
    else:
      dof    = f[1].shape[0]
      bs     = f[1].shape[1]
      cdims  = '%d 1' % dof
      dims   = '%d %d' % (dof, bs)
      stride = '1 1'
      size   = '%d 1' % dof
    for c in range(bs):
      ext = self.typeExt[spaceDim][vtype][c]
      if len(f[1].shape) > 2: start  = '%d 0 %d' % (timestep, c)
      else:                   start  = '0 %d' % c
      fp.write('''\
	<Attribute
	   Name="%s"
	   Type="Scalar"
	   Center="%s">
          <DataItem ItemType="HyperSlab"
		    Dimensions="%s"
		    Type="HyperSlab">
            <DataItem
	       Dimensions="3 %d"
	       Format="XML">
              %s
              %s
              %s
	    </DataItem>
	    <DataItem
	       DataType="Float" Precision="8"
	       Dimensions="%s"
	       Format="HDF">
	      &HeavyData;:%s
	    </DataItem>
	  </DataItem>
	</Attribute>
''' % (f[0]+'_'+ext, domain, cdims, len(f[1].shape), start, stride, size, dims, name))
    return

  def writeField(self, fp, numSteps, timestep, cellDim, spaceDim, name, f, domain):
    ctypes = ['tensor', 'matrix']
    if spaceDim == 2 or cellDim != spaceDim: ctypes.append('vector')
    if f[1].attrs['vector_field_type'] in ctypes:
      self.writeFieldComponents(fp, numSteps, timestep, spaceDim, name, f, domain)
    else:
      self.writeFieldSingle(fp, numSteps, timestep, spaceDim, name, f, domain)
    return

  def writeSpaceGridFooter(self, fp):
    fp.write('      </Grid>\n')
    return

  def writeParticleGridHeader(self, fp, numParticles, spaceDim):
    fp.write('''\
      <Grid Name="particle_domain" GridType="Uniform">
	<Topology TopologyType="Polyvertex" NodesPerElement="%d" />
	<Geometry GeometryType="%s">
	  <DataItem Reference="XML">/Xdmf/Domain/DataItem[@Name="xcoord"]</DataItem>
	  <DataItem Reference="XML">/Xdmf/Domain/DataItem[@Name="ycoord"]</DataItem>
	  %s
	</Geometry>
''' % (numParticles, "X_Y_Z" if spaceDim > 2 else "X_Y", "<DataItem Reference=\"XML\">/Xdmf/Domain/DataItem[@Name=\"zcoord\"]</DataItem>" if spaceDim > 2 else ""))
    return

  def writeParticleField(self, fp, numParticles, spaceDim):
    fp.write('''\
    <Attribute Name="particles/icoord">
      <DataItem Name="icoord"
                Format="HDF"
                Dimensions="%d">
                &HeavyData;:/particles/icoord
      </DataItem>
    </Attribute>''' % (numParticles))
    return

  def writeTimeGridFooter(self, fp):
    fp.write('    </Grid>\n')
    return

  def writeFooter(self, fp):
    fp.write('  </Domain>\n</Xdmf>\n')
    return

  def write(self, hdfFilename, topologyPath, numCells, numCorners, cellDim, geometryPath, numVertices, spaceDim, time, vfields, cfields, numParticles):
    useTime = not (len(time) < 2 and time[0] == -1)
    with file(self.filename, 'w') as fp:
      self.writeHeader(fp, hdfFilename)
      # Field information
      self.writeCells(fp, topologyPath, numCells, numCorners)
      self.writeVertices(fp, geometryPath, numVertices, spaceDim)
      if useTime: self.writeTimeGridHeader(fp, time)
      for t in range(len(time)):
        self.writeSpaceGridHeader(fp, numCells, numCorners, cellDim, spaceDim)
        for vf in vfields: self.writeField(fp, len(time), t, cellDim, spaceDim, '/vertex_fields/'+vf[0], vf, 'Node')
        for cf in cfields: self.writeField(fp, len(time), t, cellDim, spaceDim, '/cell_fields/'+cf[0], cf, 'Cell')
        self.writeSpaceGridFooter(fp)
      if useTime: self.writeTimeGridFooter(fp)
      if numParticles:
        self.writeLocations(fp, numParticles, spaceDim)
        if useTime: self.writeTimeGridHeader(fp, time)
        for t in range(len(time)):
          self.writeParticleGridHeader(fp, numParticles, spaceDim)
          self.writeSpaceGridFooter(fp)
        if useTime: self.writeTimeGridFooter(fp)
      self.writeFooter(fp)
    return

def generateXdmf(hdfFilename, xdmfFilename = None):
  if xdmfFilename is None:
    xdmfFilename = os.path.splitext(hdfFilename)[0] + '.xmf'
  # Read mesh
  h5          = h5py.File(hdfFilename, 'r')
  if 'viz' in h5 and 'geometry' in h5['viz']:
    geomPath  = 'viz/geometry'
    geom      = h5['viz']['geometry']
  else:
    geomPath  = 'geometry'
    geom      = h5['geometry']
  if 'viz' in h5 and 'topology' in h5['viz']:
    topoPath  = 'viz/topology'
    topo      = h5['viz']['topology']
  else:
    topoPath  = 'topology'
    topo      = h5['topology']
  vertices    = geom['vertices']
  numVertices = vertices.shape[0]
  spaceDim    = vertices.shape[1]
  cells       = topo['cells']
  numCells    = cells.shape[0]
  numCorners  = cells.shape[1]
  cellDim     = topo['cells'].attrs['cell_dim']
  if 'time' in h5:
    time      = np.array(h5['time']).flatten()
  else:
    time      = [-1]
  vfields     = []
  cfields     = []
  if 'vertex_fields' in h5: vfields = h5['vertex_fields'].items()
  if 'cell_fields' in h5: cfields = h5['cell_fields'].items()
  numParticles = 0
  if 'particles' in h5: numParticles = h5['particles']['xcoord'].shape[0]

  # Write Xdmf
  Xdmf(xdmfFilename).write(hdfFilename, topoPath, numCells, numCorners, cellDim, geomPath, numVertices, spaceDim, time, vfields, cfields, numParticles)
  h5.close()
  return

if __name__ == '__main__':
  for f in sys.argv[1:]:
    generateXdmf(f)
