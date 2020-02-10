import unittest
from petsc4py import PETSc
import os
from PetscBinaryIO import *

class TestPetscBinaryIO(unittest.TestCase):
    def setUp(self):
        try:
            os.remove('test.dat')
        except OSError:
            pass
        try:
            os.remove('test.dat.info')
        except OSError:
            pass

    def test_VecRead(self):
        """Test reading a Vec"""
        array = np.array([1.1, 2.2, 3.3])
        vec = PETSc.Vec().createSeq(3)
        vec[...] = array
        viewer = PETSc.Viewer().createBinary('test.dat', PETSc.Viewer.Mode.W)
        vec.view(viewer)
        viewer.destroy()
        vec.destroy()

        result, = PetscBinaryIO().readBinaryFile('test.dat')
        self.assertTrue(np.allclose(array, result))

    def test_VecWrite(self):
        """Test writing a Vec"""
        array = np.array([1.1, 2.2, 3.3])
        PetscBinaryIO().writeBinaryFile('test.dat', [array.view(Vec),])

        vec = PETSc.Vec().createSeq(3)
        vec.set(0.)
        viewer = PETSc.Viewer().createBinary('test.dat', PETSc.Viewer.Mode.R)
        vec.load(viewer)
        viewer.destroy()

        self.assertTrue(np.allclose(array, vec[...]))
        vec.destroy()

    def test_ISRead(self):
        """Test reading an IS"""
        indices = np.array([3,4,5])
        anis = PETSc.IS().createGeneral(list(indices))
        viewer = PETSc.Viewer().createBinary('test.dat', PETSc.Viewer.Mode.W)
        anis.view(viewer)
        viewer.destroy()
        anis.destroy()

        result, = PetscBinaryIO().readBinaryFile('test.dat')
        self.assertTrue((indices == result).all())

    def test_MatRead(self):
        """Test reading a Mat"""
        mat = PETSc.Mat().createAIJ(2)
        mat[0,0] = 1.1
        mat[0,1] = 2.1
        mat[1,1] = 3.1
        mat.assemble()

        vals = np.array([1.1,2.1,3.1])
        counts = np.array([0,2,3])
        cols = np.array([0,1,1])

        viewer = PETSc.Viewer().createBinary('test.dat', PETSc.Viewer.Mode.W)
        mat.view(viewer)
        viewer.destroy()
        mat.destroy()

        result, = PetscBinaryIO().readBinaryFile('test.dat')
        self.assertTrue(np.allclose(vals, result[1][2]))
        self.assertTrue((counts == result[1][0]).all())
        self.assertTrue((cols == result[1][1]).all())
        self.assertTrue((2,2) == result[0])

    def test_MatWrite(self):
        """Test writing a Mat"""
        vals = np.array([1.1,2.1,3.1])
        counts = np.array([0,2,3])
        cols = np.array([0,1,1])
        mat = MatSparse(((2,2),(counts,cols,vals)))

        dense = np.array([1.1,2.1,0.0,3.1])

        PetscBinaryIO().writeBinaryFile('test.dat', [mat,])

        mat = PETSc.Mat().createAIJ(2)
        viewer = PETSc.Viewer().createBinary('test.dat', PETSc.Viewer.Mode.R)
        mat.load(viewer)
        viewer.destroy()

        self.assertTrue(np.allclose(dense, mat[:,:].ravel()))
        mat.destroy()


if __name__ == '__main__':
    unittest.main()
    try:
        os.remove('test.dat')
    except OSError:
        pass
    try:
        os.remove('test.dat.info')
    except OSError:
        pass
