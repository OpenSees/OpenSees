"""PetscBinaryIO
===============

Provides
  1. PETSc-named objects Vec, Mat, and IS that inherit numpy.ndarray
  2. A class to read and write these objects from PETSc binary files.

The standard usage of this module should look like:

  >>> import PetscBinaryIO
  >>> io = PetscBinaryIO.PetscBinaryIO()
  >>> objects = io.readBinaryFile('file.dat')

or

  >>> import PetscBinaryIO
  >>> import numpy
  >>> vec = numpy.array([1., 2., 3.]).view(PetscBinaryIO.Vec)
  >>> io = PetscBinaryIO.PetscBinaryIO()
  >>> io.writeBinaryFile('file.dat', [vec,])

to read in objects one at a time use such as

  >>> import PetscBinaryIO
  >>> io = PetscBinaryIO.PetscBinaryIO()
  >>> fh = open('file.dat')
  >>> objecttype = io.readObjectType(fh)
  >>> if objecttype == 'Vec':
  >>>   v = io.readVec(fh)

   Note that one must read in the object type first and then call readVec(), readMat() etc.


See also PetscBinaryIO.__doc__ and methods therein.
"""

import numpy as np
import functools

try:
    basestring                  # Python-2 has basestring as a common parent of unicode and str
except NameError:
    basestring = str            # Python-3 is unicode through and through

def update_wrapper_with_doc(wrapper, wrapped):
    """Similar to functools.update_wrapper, but also gets the wrapper's __doc__ string"""
    wdoc = wrapper.__doc__

    functools.update_wrapper(wrapper, wrapped)
    if wdoc is not None:
        if wrapper.__doc__ is None:
            wrapper.__doc__ = wdoc
        else:
            wrapper.__doc__ = wrapper.__doc__ + wdoc
    return wrapper

def wraps_with_doc(wrapped):
    """Similar to functools.wraps, but also gets the wrapper's __doc__ string"""
    return functools.partial(update_wrapper_with_doc, wrapped=wrapped)

def decorate_with_conf(f):
    """Decorates methods to take kwargs for precisions."""
    @wraps_with_doc(f)
    def decorated_f(self, *args, **kwargs):
        """
        Additional kwargs:
          precision: 'single', 'double', '__float128' for scalars
          indices: '32bit', '64bit' integer size
          complexscalars: True/False

          Note these are set in order of preference:
            1. kwargs if given here
            2. PetscBinaryIO class __init__ arguments
            3. PETSC_DIR/PETSC_ARCH defaults
        """

        changed = False
        old_precision = self.precision
        old_indices = self.indices
        old_complexscalars = self.complexscalars

        try:
            self.precision = kwargs.pop('precision')
        except KeyError:
            pass
        else:
            changed = True

        try:
            self.indices = kwargs.pop('indices')
        except KeyError:
            pass
        else:
            changed = True

        try:
            self.complexscalars = kwargs.pop('complexscalars')
        except KeyError:
            pass
        else:
            changed = True

        if changed:
            self._update_dtypes()

        result = f(self, *args, **kwargs)

        if changed:
            self.precision = old_precision
            self.indices = old_indices
            self.complexscalars = old_complexscalars
            self._update_dtypes()

        return result
    return decorated_f


class DoneWithFile(Exception): pass


class Vec(np.ndarray):
    """Vec represented as 1D numpy array

    The best way to instantiate this class for use with writeBinaryFile()
    is through the numpy view method:

    vec = numpy.array([1,2,3]).view(Vec)
    """
    _classid = 1211214


class MatDense(np.matrix):
    """Mat represented as 2D numpy array

    The best way to instantiate this class for use with writeBinaryFile()
    is through the numpy view method:

    mat = numpy.array([[1,0],[0,1]]).view(Mat)
    """
    _classid = 1211216


class MatSparse(tuple):
    """Mat represented as CSR tuple ((M, N), (rowindices, col, val))

    This should be instantiated from a tuple:

    mat = MatSparse( ((M,N), (rowindices,col,val)) )
    """
    _classid = 1211216
    def __repr__(self):
        return 'MatSparse: %s'%super(MatSparse, self).__repr__()


class IS(np.ndarray):
    """IS represented as 1D numpy array

    The best way to instantiate this class for use with writeBinaryFile()
    is through the numpy "view" method:

    an_is = numpy.array([3,4,5]).view(IS)
    """
    _classid = 1211218


class PetscBinaryIO(object):
    """Reader/Writer class for PETSc binary files.

    Note that by default, precisions for both scalars and indices, as well as
    complex scalars, are picked up from the PETSC_DIR/PETSC_ARCH configuration
    as set by environmental variables.

    Alternatively, defaults can be overridden at class instantiation, or for
    a given method call.
    """

    _classid = {1211216:'Mat',
                1211214:'Vec',
                1211218:'IS',
                1211219:'Bag',
                1211213:'Real'}

    def __init__(self, precision=None, indices=None, complexscalars=None):
        if (precision is None) or (indices is None) or (complexscalars is None):
            import petsc_conf
            defaultprecision, defaultindices, defaultcomplexscalars = petsc_conf.get_conf()
            if precision is None:
                if defaultprecision is None:
                    precision = 'double'
                else:
                    precision = defaultprecision

            if indices is None:
                if defaultindices is None:
                    indices = '32bit'
                else:
                    indices = defaultindices

            if complexscalars is None:
                if defaultcomplexscalars is None:
                    complexscalars = False
                else:
                    complexscalars = defaultcomplexscalars

        self.precision = precision
        if self.precision == '__float128' :
            raise RuntimeError('__float128 (quadruple) precision is not properly supported. One may use double precision by using -binary_write_double in PETSc and precision=\'double\' here')
        self.indices = indices
        self.complexscalars = complexscalars
        self._update_dtypes()

    def _update_dtypes(self):
        if self.indices == '64bit':
            self._inttype = np.dtype('>i8')
        else:
            self._inttype = np.dtype('>i4')

        if self.precision == '__float128':
            nbyte = 16
        elif self.precision == 'single':
            nbyte = 4
        else:
            nbyte = 8

        if self.complexscalars:
            name = 'c'
            nbyte = nbyte * 2 # complex scalar takes twice as many bytes
        else:
            name = 'f'

        self._scalartype = '>{0}{1}'.format(name, nbyte)

    @decorate_with_conf
    def readReal(self, fh):
        """Reads a single real from a binary file handle, must be called after readObjectType()."""

        try:
            vals = np.fromfile(fh, dtype=self._scalartype, count=1)
        except MemoryError:
            raise IOError('Inconsistent or invalid real data in file')
        if (len(vals) is 0):
            raise IOError('Inconsistent or invalid real data in file')
        return vals

    @decorate_with_conf
    def readVec(self, fh):
        """Reads a PETSc Vec from a binary file handle, must be called after readObjectType()."""

        nz = np.fromfile(fh, dtype=self._inttype, count=1)[0]
        try:
            vals = np.fromfile(fh, dtype=self._scalartype, count=nz)
        except MemoryError:
            raise IOError('Inconsistent or invalid Vec data in file')
        if (len(vals) is 0):
            raise IOError('Inconsistent or invalid Vec data in file')
        return vals.view(Vec)

    @decorate_with_conf
    def writeVec(self, fh, vec):
        """Writes a PETSc Vec to a binary file handle."""

        metadata = np.array([Vec._classid, len(vec)], dtype=self._inttype)
        metadata.tofile(fh)
        vec.astype(self._scalartype).tofile(fh)
        return

    @decorate_with_conf
    def readMatSparse(self, fh):
        """Reads a PETSc Mat, returning a sparse representation of the data. Must be called after readObjectType()

        (M,N), (I,J,V) = readMatSparse(fid)

        Input:
          fid : file handle to open binary file.
        Output:
          M,N : matrix size
          I,J : arrays of row and column for each nonzero
          V: nonzero value
        """

        try:
            M,N,nz = np.fromfile(fh, dtype=self._inttype, count=3)
            I = np.empty(M+1, dtype=self._inttype)
            I[0] = 0
            rownz = np.fromfile(fh, dtype=self._inttype, count=M)
            np.cumsum(rownz, out=I[1:])
            assert I[-1] == nz

            J = np.fromfile(fh, dtype=self._inttype,    count=nz)
            assert len(J) == nz
            V = np.fromfile(fh, dtype=self._scalartype, count=nz)
            assert len(V) == nz
        except (AssertionError, MemoryError, IndexError):
            raise IOError('Inconsistent or invalid Mat data in file')

        return MatSparse(((M, N), (I, J, V)))

    @decorate_with_conf
    def writeMatSparse(self, fh, mat):
        """Writes a Mat into a PETSc binary file handle"""

        ((M,N), (I,J,V)) = mat
        metadata = np.array([MatSparse._classid,M,N,I[-1]], dtype=self._inttype)
        rownz = I[1:] - I[:-1]

        assert len(J.shape) == len(V.shape) == len(I.shape) == 1
        assert len(J) == len(V) == I[-1] == rownz.sum()
        assert (rownz > -1).all()

        metadata.tofile(fh)
        rownz.astype(self._inttype).tofile(fh)
        J.astype(self._inttype).tofile(fh)
        V.astype(self._scalartype).tofile(fh)
        return

    @decorate_with_conf
    def readMatDense(self, fh):
        """Reads a PETSc Mat, returning a dense represention of the data, must be called after readObjectType()"""

        try:
            M,N,nz = np.fromfile(fh, dtype=self._inttype, count=3)
            I = np.empty(M+1, dtype=self._inttype)
            I[0] = 0
            rownz = np.fromfile(fh, dtype=self._inttype, count=M)
            np.cumsum(rownz, out=I[1:])
            assert I[-1] == nz

            J = np.fromfile(fh, dtype=self._inttype, count=nz)
            assert len(J) == nz
            V = np.fromfile(fh, dtype=self._scalartype, count=nz)
            assert len(V) == nz

        except (AssertionError, MemoryError, IndexError):
            raise IOError('Inconsistent or invalid Mat data in file')

        mat = np.zeros((M,N), dtype=self._scalartype)
        for row in range(M):
            rstart, rend = I[row:row+2]
            mat[row, J[rstart:rend]] = V[rstart:rend]
        return mat.view(MatDense)

    @decorate_with_conf
    def readMatSciPy(self, fh):
        from scipy.sparse import csr_matrix
        (M, N), (I, J, V) = self.readMatSparse(fh)
        return csr_matrix((V, J, I), shape=(M, N))

    @decorate_with_conf
    def writeMatSciPy(self, fh, mat):
        from scipy.sparse import csr_matrix
        if hasattr(mat, 'tocsr'):
            mat = mat.tocsr()
        assert isinstance(mat, csr_matrix)
        V = mat.data
        M,N = mat.shape
        J = mat.indices
        I = mat.indptr
        return self.writeMatSparse(fh, (mat.shape, (mat.indptr,mat.indices,mat.data)))

    @decorate_with_conf
    def readMat(self, fh, mattype='sparse'):
        """Reads a PETSc Mat from binary file handle, must be called after readObjectType()

        optional mattype: 'sparse" or 'dense'

        See also: readMatSparse, readMatDense
        """

        if mattype == 'sparse':
            return self.readMatSparse(fh)
        elif mattype == 'dense':
            return self.readMatDense(fh)
        elif mattype == 'scipy.sparse':
            return self.readMatSciPy(fh)
        else:
            raise RuntimeError('Invalid matrix type requested: choose sparse/dense/scipy.sparse')

    @decorate_with_conf
    def readIS(self, fh):
        """Reads a PETSc Index Set from binary file handle, must be called after readObjectType()"""

        try:
            nz = np.fromfile(fh, dtype=self._inttype, count=1)[0]
            v = np.fromfile(fh, dtype=self._inttype, count=nz)
            assert len(v) == nz
        except (MemoryError,IndexError):
            raise IOError('Inconsistent or invalid IS data in file')
        return v.view(IS)

    @decorate_with_conf
    def writeIS(self, fh, anis):
        """Writes a PETSc IS to binary file handle."""

        metadata = np.array([IS._classid, len(anis)], dtype=self._inttype)
        metadata.tofile(fh)
        anis.astype(self._inttype).tofile(fh)
        return

    @decorate_with_conf
    def readObjectType(self, fid):
        """Returns the next object type as a string in the file"""
        try:
              header = np.fromfile(fid, dtype=self._inttype, count=1)[0]
        except (MemoryError, IndexError):
              raise DoneWithFile
        try:
              objecttype = self._classid[header]
        except KeyError:
              raise IOError('Invalid PetscObject CLASSID or object not implemented for python')
        return objecttype

    @decorate_with_conf
    def readBinaryFile(self, fid, mattype='sparse'):
        """Reads a PETSc binary file, returning a tuple of the contained objects.

        objects = self.readBinaryFile(fid, **kwargs)

        Input:
          fid : either file name or handle to an open binary file.

        Output:
          objects : tuple of objects representing the data in numpy arrays.

        Optional:
          mattype :
            'sparse': Return matrices as raw CSR: (M, N), (row, col, val).
            'dense': Return matrices as MxN numpy arrays.
            'scipy.sparse': Return matrices as scipy.sparse objects.
        """

        close = False

        if isinstance(fid, basestring):
            fid = open(fid, 'rb')
            close = True

        objects = []
        try:
            while True:
                objecttype = self.readObjectType(fid)

                if objecttype == 'Vec':
                    objects.append(self.readVec(fid))
                elif objecttype == 'IS':
                    objects.append(self.readIS(fid))
                elif objecttype == 'Mat':
                    objects.append(self.readMat(fid,mattype))
                elif objecttype == 'Real':
                    objects.append(self.readReal(fid))
                elif objecttype == 'Bag':
                    raise NotImplementedError('Bag Reader not yet implemented')
        except DoneWithFile:
            pass
        finally:
            if close:
                fid.close()

        return tuple(objects)

    @decorate_with_conf
    def writeBinaryFile(self, fid, objects):
        """Writes a PETSc binary file containing the objects given.

        readBinaryFile(fid, objects)

        Input:
          fid : either file handle to an open binary file, or filename.
          objects : list of objects representing the data in numpy arrays,
                    which must be of type Vec, IS, MatSparse, or MatSciPy.
        """
        close = False
        if isinstance(fid, basestring):
            fid = open(fid, 'wb')
            close = True

        for petscobj in objects:
            if (isinstance(petscobj, Vec)):
                self.writeVec(fid, petscobj)
            elif (isinstance(petscobj, IS)):
                self.writeIS(fid, petscobj)
            elif (isinstance(petscobj, MatSparse)):
                self.writeMatSparse(fid, petscobj)
            elif (isinstance(petscobj, MatDense)):
                if close:
                    fid.close()
                raise NotImplementedError('Writing a dense matrix is not yet supported')
            else:
                try:
                    self.writeMatSciPy(fid, petscobj)
                except AssertionError:
                    if close:
                        fid.close()
                    raise TypeError('Object %s is not a valid PETSc object'%(petscobj.__repr__()))
        if close:
            fid.close()
        return
