/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
** University of California, Berkeley, is strictly prohibited.  See   **
** file 'COPYRIGHT'  in main directory for information on usage and   **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
** Developed by:                                                      **
**   Frank McKenna (fmckenna@ce.berkeley.edu)                         **
**   Gregory L. Fenves (fenves@ce.berkeley.edu)                       **
**   Filip C. Filippou (filippou@ce.berkeley.edu)                     **
**                                                                    **
** ****************************************************************** */

// Description: This file contains the implementation for PetscSOE

#include <limits>


#include <petsc.h>
#include <PetscSOE.h>
#include <PetscSolver.h>
#include <petscmat.h>
#include <petscvec.h>
#include <petscsys.h>

#include <Matrix.h>
#include <Graph.h>
#include <Vertex.h>
#include <VertexIter.h>
#include <f2c.h>
#include <math.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>


#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string>
#include <sstream>
#include <petscpc.h>


using namespace std;




PetscSOE::PetscSOE(PetscSolver &theSOESolver, MatType matType, int bs)
  : LinearSOE(theSOESolver, LinSOE_TAGS_PetscSOE),
    isFactored(0), size(0), processID(0), numProcesses(0), B(0), X(0),
    indices(0), vectX(0), vectB(0), A(0), M(0),  x(0), b(0), blockSize(bs),
    numChannels(0), theChannels(0),  localCol(0), mType(matType)
{
  theSOESolver.setLinearSOE(*this);
  init_done = false;
  startRow = -1;
  endRow = -1;
}

int
PetscSOE::getNumEqn(void) const
{
  return size;
}

PetscSOE::~PetscSOE()
{
  cout << "PetscSOE::~PetscSOE()\n";
  if (theChannels != 0)
  {
    delete [] theChannels;
  }

  if (localCol != 0)
    for (int i = 0; i < numChannels; i++)
      if (localCol[i] != 0)
      {
        delete localCol[i];
      }

  delete [] localCol;

  if (vectX != 0)
  {
    delete vectX;
  }

  if (vectB != 0)
  {
    delete vectB;
  }

  if (B != 0)
  {
    delete [] B;
  }

  if (X != 0)
  {
    delete [] X;
  }

  // invoke the petsc destructors
  if (A != 0)
  {
    MatDestroy(&A);
  }

  if (M != 0)
  {
    MatDestroy(&M);
  } ;

  if (b != 0)
  {
    VecDestroy(&b);
  }

  if (x != 0)
  {
    VecDestroy(&x);
  }

  if (processID >= 0)
  {
    cout << "Processor " << processID << " printing PETSc log\n";
    PetscViewer    viewer;
    
    PetscViewerASCIIGetStdout(PETSC_COMM_WORLD, &viewer);
    PetscLogView(viewer);
    PetscFinalize();
  }


}


int
PetscSOE::setSize(Graph &theGraph)
{
  int result = 0;
  int ierr = 0;




  if (!init_done)
  {

    MPI_Comm_size(MPI_COMM_WORLD, &numProcesses);
    MPI_Comm_rank(MPI_COMM_WORLD, &processID);

    //Check if petsc_options.txt is present
    if (FILE *file = fopen(PETSCSOE_PETSCOPTIONS_FILENAME, "r"))
    {
      fclose(file);
      cout << "        + PetscSOE::setSize() [" << processID <<  "] : Found petsc_options.txt file. \n";
      PetscInitialize(0, PETSC_NULL, PETSCSOE_PETSCOPTIONS_FILENAME, "PetscSOE::setSize - Initialize with file");
      PetscPopSignalHandler();
    }
    else
    {
      PetscInitialize(0, PETSC_NULL, (char *)0, "PetscSOE::setSize - Initialize");
      PetscPopSignalHandler();
    }

    PetscLogAllBegin();


    init_done = true;
  }


  //
  // first determine system size
  //

  PETSCSOE_DEBUGOUT << "        + PetscSOE::setSize() [" << processID <<  "] : Getting local max dof tag\n";

  // Size of the system is the maximum DOF tag + 1

  // Each processor visits its DOF graph and determines the maximum DOF tag
  // for its domain.
  Vertex* theVertex;
  VertexIter& theVertices = theGraph.getVertices();

  #ifdef max // This is needed because, to my frustration, someone, somewhere overloaded max()
  #undef max // Angrily, J. Abell
  #endif

  int mindoftag = std::numeric_limits<int>::max();
  int maxdoftag = 0;

  while ((theVertex = theVertices()) != 0)
  {
    int tag = theVertex->getTag();

    if (tag > maxdoftag)
    {
      maxdoftag = tag;
    }
    if (tag < mindoftag)
    {
      mindoftag = tag;
    }
  }
  size = maxdoftag;


  // Then, all maximumg tags (size) are sent to P0 which gets the global maximum.

  PETSCSOE_DEBUGOUT << "        + PetscSOE::setSize() [" << processID <<  "] : Determining global maximum dof tag\n";
  isFactored = 0;
  static ID data(1);
  if (processID != 0)
  {
    Channel *theChannel = theChannels[0];

    data(0) = size;
    theChannel->sendID(0, 0, data);
    theChannel->recvID(0, 0, data);

    size = data(0);
  }
  else
  {

    for (int j = 0; j < numChannels; j++)
    {
      Channel *theChannel = theChannels[j];
      theChannel->recvID(0, 0, data);

      if (data(0) > size)
      {
        size = data(0);
      }
    }

    data(0) = size;

    for (int j = 0; j < numChannels; j++)
    {
      Channel *theChannel = theChannels[j];
      theChannel->sendID(0, 0, data);
    }
  }

  size = size + 1; // vertices numbered 0 through n-1, and n is the size of the system. therefore, add 1

  PETSCSOE_DEBUGOUT << "        + PetscSOE::setSize() [" << processID <<  "] : size = "  << size << "\n";


  PETSCSOE_DEBUGOUT << "        + PetscSOE::setSize() [" << processID <<  "] : Allocating memory\n";

  // invoke the petsc destructors
  if (A != 0)
  {
    MatDestroy(&A);
  }

  if (b != 0)
  {
    VecDestroy(&b);
  }

  if (x != 0)
  {
    VecDestroy(&x);
  }

  //
  // now we create the opensees vector objects
  //

  // delete the old vectors
  if (B != 0)
  {
    delete [] B;
  }

  if (X != 0)
  {
    delete [] X;
  }

  // create the new
  B = new double[size];
  X = new double[size];

  if (B == 0 || X == 0)
  {
    cerr << "WARNING PetscSOE::PetscSOE :";
    cerr << " ran out of memory for vectors (size) (";
    cerr << size << ") \n";
    size = 0;
    result = -1;
  }

  // zero the vectors
  for (int j = 0; j < size; j++)
  {
    B[j] = 0;
    X[j] = 0;
  }

  if (vectX != 0)
  {
    delete vectX;
  }

  if (vectB != 0)
  {
    delete vectB;
  }

  vectX = new Vector(X, size);

  if (vectX == NULL)
  {
    std::cerr << "PetscSOE::PetscSOE : can not allocate memory for vectX \n";
  }

  vectB = new Vector(B, size);


  if (vectB == NULL)
  {
    std::cerr << "PetscSOE::PetscSOE : can not allocate memory for vectB \n";
  }


  // Determine start and end rows for each processor
  // ================================================

  PETSCSOE_DEBUGOUT << "        + PetscSOE::setSize() [" << processID <<  "] : Distributing SOE rows amongst processors\n";

  // Create array to store number of dofs in each processor (allgather)
  int nlocaldofs[numProcesses];

  //Fill with -1 : indicating unused component
  for (int i = 0; i < numProcesses; i++)
  {
    nlocaldofs[i] = -1;
  }

  //Make component corresponding to local component = the number of DOFS in this proc
  nlocaldofs[processID] = theGraph.getNumVertex();

  int my_nlocaldofs = nlocaldofs[processID];
  // cout << "Processor " << processID_world << " owns " << nlocaldofs[processID_world] << " DOFS.\n";

  //Now gather across all processes so we can determine a reasonable partition of
  MPI_Allgather(&my_nlocaldofs, 1, MPI_INT, &nlocaldofs, 1, MPI_INT,
                MPI_COMM_WORLD);

  // Count the excess, ie. dofs which are shared across processes. These need to be
  // evenly distributed.
  int excess = 0;
  for (int i = 0; i < numProcesses; i++)
  {
    excess += nlocaldofs[i];
  }
  excess -= size;

  // Remove excess dofs uniformly from all processes
  for (int i = 0; i < excess; i++)
  {
    nlocaldofs[i % (numProcesses - 1) + 1]--;
  }


  int startRow_vec[numProcesses];
  int endRow_vec[numProcesses];

  // Assign start and end rows
  startRow_vec[0] = 0; // P0 gets nothing
  endRow_vec[0] = 0;

  int dof = 0;
  for (int p = 0; p < numProcesses; p++)
  {
    startRow_vec[p] = dof;
    endRow_vec[p] = dof + nlocaldofs[p] - 1;
    dof += nlocaldofs[p];
  }
  startRow = startRow_vec[processID];
  endRow = endRow_vec[processID];

  //
  // now create petsc matrix and vector objects
  //

  //
  // Call Petsc VecCreate & MatCreate; NOTE: using previously allocated storage for vectors
  //
  //


  if (processID >= 0)
  {

    PETSCSOE_DEBUGOUT << "        + PetscSOE::setSize() - [" << processID <<  "] : Initializing PETSc\n";
    // PetscOptionsGetInt(PETSC_NULL, "-n", &size, PETSC_NULL);

    ierr = MatCreate(PETSC_COMM_WORLD, &A);
    CHKERRQ(ierr);


    ierr = MatSetType(A, MATMPIAIJ);
    CHKERRQ(ierr);


    int ndofs = nlocaldofs[processID];

    // Performance opportunity.
    //Can use 2nd and 3rd parameters to customize which rows belong to what processor.
    // ierr = MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, size, size);
    ierr = MatSetSizes(A, ndofs, ndofs, size, size);
    // CHKERRQ(ierr);


    // Disable auto-ownership
    // MatGetOwnershipRange(A, &startRow, &endRow);

    PETSCSOE_DEBUGOUT << "        + PetscSOE::setSize() - Process " << processID << " owns rows from " << startRow << " to " << endRow << endl;

    //Define custom ownership
    PetscInt* d_nnz[numProcesses];
    PetscInt* o_nnz[numProcesses];
    PetscInt d_nnz_global[ndofs];
    PetscInt o_nnz_global[ndofs];

    for (int proc = 0; proc < numProcesses; proc++)
    {
      int ndofs_proc = nlocaldofs[proc];
      d_nnz[proc] = new PetscInt[ndofs_proc];
      o_nnz[proc] = new PetscInt[ndofs_proc];
    }

   
    // Initialize diagonal to at least one nonzero element (the diagonal) for current proc
    // Initialize off-diagonal to zero
    for (int proc = 0; proc < numProcesses; proc++)
    {
      PetscInt dval = 0;
      if (proc == processID)
      {
        dval = 1;
      }

      for (int row = 0; row < nlocaldofs[proc]; row++)
      {
        d_nnz[proc][row] = dval;
        o_nnz[proc][row] = 0;
      }
    }

    //Iterate the vertices again to determine the number of nonzeros at each dof....
    Vertex* theVertex;
    VertexIter& theVertices = theGraph.getVertices();

    while ((theVertex = theVertices()) != 0)
    {
      int dof = theVertex->getTag();
      const ID& adj = theVertex->getAdjacency();
      //Iterate over adjacent dofs and determine whether they are in the "diagonal" block of the PETSc matrix
      //  or not....
      // Diagonal here means that the adjacent dof also belongs to this processor
      for (int i = 0; i < adj.Size(); i++)
      {
        int col = adj(i);

        for (int proc = 0; proc < numProcesses; proc++)
        {
          // If the dof number belongs to this processor
          if ( startRow_vec[proc] <= dof && dof <= endRow_vec[proc]   )
          {
            int row = dof - startRow_vec[proc];
            if ( startRow_vec[proc] <= col && col <= endRow_vec[proc]   )
            {
              d_nnz[proc][row ] += 1;
            }
            else
            {
              o_nnz[proc][row ] += 1;
            }
          }

          // Check the symmetric condition of this (dofs are symmetricly connected!)
          // If the dof number belongs to this processor
          if ( startRow_vec[proc] <= col && col <= endRow_vec[proc]   )
          {
            int row = col - startRow_vec[proc];
            if ( startRow_vec[proc] <= dof && dof <= endRow_vec[proc]   )
            {
              d_nnz[proc][row ] += 1;
            }
            else
            {
              o_nnz[proc][row ] += 1;
            }
          }
        }
      }
    }

    // Now reduce the arrays to one global one in the local processor
    for (int proc = 0; proc < numProcesses; proc++)
    {
      MPI_Reduce(
        d_nnz[proc],
        d_nnz_global,
        nlocaldofs[proc],
        MPI_INT,
        MPI_SUM,
        proc , // Because this is the process number in the petsc_comm communicator
        PETSC_COMM_WORLD);
      MPI_Reduce(
        o_nnz[proc],
        o_nnz_global,
        nlocaldofs[proc],
        MPI_INT,
        MPI_SUM,
        proc , // Because this is the process number in the petsc_comm communicator
        PETSC_COMM_WORLD);
    }




    for (int row = 0; row < nlocaldofs[processID]; row++)
    {
      if (d_nnz_global[row] < PETSCSOE_MIN_DNNZ)
      {
        d_nnz_global[row] =  PETSCSOE_MIN_DNNZ;
      }
      if (o_nnz_global[row] < PETSCSOE_MIN_ONNZ)
      {
        o_nnz_global[row] =  PETSCSOE_MIN_ONNZ;
      }
      if (d_nnz_global[row] > ndofs)
      {
        d_nnz_global[row] =  ndofs;
      }
      if (o_nnz_global[row] > ndofs)
      {
        o_nnz_global[row] =  ndofs;
      }
    }


    // ierr = MatCreateAIJ(PETSC_COMM_WORLD, ndofs, ndofs, size, size, 0, d_nnz, 0, o_nnz, &A);
    // CHKERRQ(ierr);

    // PetscErrorCode  MatMPIAIJSetPreallocation(Mat B,PetscInt d_nz,const PetscInt d_nnz[],PetscInt o_nz,const PetscInt o_nnz[])
    // B   - the matrix
    // d_nz    - number of nonzeros per row in DIAGONAL portion of local submatrix (same value is used for all local rows)
    // d_nnz   - array containing the number of nonzeros in the various rows of the DIAGONAL portion of the local submatrix (possibly different for each row) or NULL (PETSC_NULL_INTEGER in Fortran), if d_nz is used to specify the nonzero structure. The size of this array is equal to the number of local rows, i.e 'm'. For matrices that will be factored, you must leave room for (and set) the diagonal entry even if it is zero.
    // o_nz    - number of nonzeros per row in the OFF-DIAGONAL portion of local submatrix (same value is used for all local rows).
    // o_nnz   - array containing the number of nonzeros in the various rows of the OFF-DIAGONAL portion of the local submatrix (possibly different for each row) or NULL (PETSC_NULL_INTEGER in Fortran), if o_nz is used to specify the nonzero structure. The size of this array is equal to the number of local rows, i.e 'm'.

    MatMPIAIJSetPreallocation(A, 0, d_nnz_global, 0, o_nnz_global);

    // ierr = MatMPIAIJSetPreallocation(A, 650, PETSC_NULL, 650, PETSC_NULL);
    // ierr = MatMPIAIJSetPreallocation(A, 650, PETSC_NULL, 650, PETSC_NULL);

    // CHKERRQ(ierr);
    // ierr = MatSeqAIJSetPreallocation(A, 650 , PETSC_NULL);
    // CHKERRQ(ierr);


    MatSetOption(A, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);



    ierr = VecCreateMPIWithArray(PETSC_COMM_WORLD, blockSize, ndofs , size, &X[startRow], &x);
    CHKERRQ(ierr);
    ierr = VecCreateMPIWithArray(PETSC_COMM_WORLD, blockSize, ndofs , size, &B[startRow], &b);
    CHKERRQ(ierr);

    for (int proc = 0; proc < numProcesses; proc++)
    {
      delete [] d_nnz[proc];
      delete [] o_nnz[proc];
    }

    // cout << "Preallocation done! " << endl;
  }
  else
  {
    // cout << "Process " << processID << " owns rows from " << startRow << " to " << endRow << endl;
  }

  // invoke setSize() on the Solver
  LinearSOESolver *tSolver = this->getSolver();
  int solverOK = tSolver->setSize();

  if (solverOK < 0)
  {
    cerr << "WARNING:PetscSOE::setSize :";
    cerr << " solver failed setSize()\n";
    return solverOK;
  }

  // MatGetOwnershipRange(A, &startRow, &endRow);


  PETSCSOE_DEBUGOUT << "Process " << processID << " owns rows from " << startRow << " to " << endRow << endl;


  return result;
}



int
PetscSOE::setSize(int MaxDOFtag)
{
  int result = -1;
  cerr << "PetscSOE::setSize(int MaxDOFtag) -- DEPRECATED! Should not get called\n";

  return result;

}


int
PetscSOE::addA(const Matrix &m, const ID &id, double fact)
{

  if (processID >= 0)
  {
    isFactored = 0;

    // check for a quick return
    if (fact == 0.0)
    {
      return 0;
    }


    // check that m and id are of similar size
    int idSize = id.Size();

    if (idSize != m.noRows() && idSize != m.noCols())
    {
      cerr << "PetscSOE::addA() - Matrix and ID not of similar sizes\n";
      return -1;
    }

    int n = id.Size();


    int indx_ij[n ];

    double values[n * n];


    // const int * iddata = id.getData();
    // const double * matdata = m.getData();
    int pos = 0;
    int nsmall = 0;

    for (int i = 0; i < n; i++)
    {
      if (id(i) >= 0)
      {
        indx_ij[nsmall] = id(i);
        nsmall++;
      }
    }

    int ii = 0;
    for (int i = 0; i < n; i++)
    {
      int jj = 0;
      if (id(i) >= 0)
      {
        for (int j = 0; j < n; j++) // Row-major, like C
        {
          if (id(j) >= 0)
          {
            values[ii * nsmall + jj] = m(i, j) * fact;
            pos++;
            jj++;
          }
        }
        ii++;
      }
    }

    int ierr = MatSetValues(A, nsmall, indx_ij, nsmall, indx_ij, values, ADD_VALUES);

    CHKERRQ(ierr);
  }

  return 0;
}


int
PetscSOE::addB(const Vector &v, const ID &id, double fact)
{
  if (processID >= 0)
  {
    // check for a quick return
    if (fact == 0.0)
    {
      return 0;
    }


    // check that m and id are of similar size
    int idSize = id.Size();

    if (idSize != v.Size() )
    {
      cerr << "addB::addB()	- Vector and ID not of similar sizes\n";
      return -1;
    }

    if (fact == 1.0)   // do not need to multiply if fact == 1.0
    {
      for (int i = 0; i < idSize; i++)
      {
        int pos = id(i);

        if (pos < size && pos >= 0)
        {
          B[pos] += v(i);
          if (B[pos] != B[pos])
          {
            cerr << "PetscSOE::addB() - rank = " << processID << "- NaN found! At pos=" << pos << endl;
            return -1;
          }
        }
      }
    }
    else if (fact == -1.0)
    {
      for (int i = 0; i < idSize; i++)
      {
        int pos = id(i);

        if (pos < size && pos >= 0)
        {
          B[pos] -= v(i);
        }
      }
    }
    else
    {
      for (int i = 0; i < idSize; i++)
      {
        int pos = id(i);

        if (pos < size && pos >= 0)
        {
          B[pos] += v(i) * fact;
        }
      }
    }
  }

  return 0;
}


int
PetscSOE::setB(const Vector &v, double fact)
{
  if (processID >= 0)
  {
    // check for a quick return
    if (fact == 0.0)
    {
      return 0;
    }


    if (v.Size() != size)
    {
      cerr << "WARNING BandGenLinSOE::setB() -";
      cerr << " incomptable sizes " << size << " and " << v.Size() << endln;
      return -1;
    }

    if (fact == 1.0)   // do not need to multiply if fact == 1.0
    {
      for (int i = 0; i < size; i++)
      {
        B[i] = v(i);
      }
    }
    else if (fact == -1.0)
    {
      for (int i = 0; i < size; i++)
      {
        B[i] = -v(i);
      }
    }
    else
    {
      for (int i = 0; i < size; i++)
      {
        B[i] = v(i) * fact;
      }
    }
  }
  return 0;
}


void
PetscSOE::zeroA(void)
{
  isFactored = 0;
  if ( processID >= 0 )
  {
    MatZeroEntries(A);
  }
}

void
PetscSOE::zeroB(void)
{
  // if ( processID_world > 0 )
  // {

  double *Bptr = B;

  for (int i = 0; i < size; i++)
  {
    *Bptr++ = 0;
  }

  // VecZeroEntries(b);
  // }
}


const Vector &
PetscSOE::getX(void)
{
  if (vectX == 0)
  {
    cerr << "FATAL PetscSOE::getX - vectX == 0!";
    exit(-1);
  }

  return *vectX;
}


const Vector &
PetscSOE::getB(void)
{
  if (vectB == 0)
  {
    cerr << "FATAL PetscSOE::getB - vectB == 0!";
    exit(-1);
  }

  return *vectB;
}


double
PetscSOE::normRHS(void)
{
  this->getB();
  double norm = 0.0;
  double *Bptr = B;

  for (int i = 0; i < size; i++)
  {
    double Yi = *Bptr++;
    norm += Yi * Yi;
  }

  return sqrt(norm);
}


void
PetscSOE::setX(int loc, double value)
{
  if (loc < size && loc >= 0)
  {
    X[loc] = value;
  }
}

void
PetscSOE::setX(const Vector &xData)
{
  if (xData.Size() == size && vectX != 0)
  {
    *vectX = xData;
  }
}

int
PetscSOE::setSolver(PetscSolver &newSolver)
{
  newSolver.setLinearSOE(*this);

  if (size != 0)
  {
    int solverOK = newSolver.setSize();

    if (solverOK < 0)
    {
      cerr << "WARNING:PetscSOE::setSolver :";
      cerr << "the new solver could not setSeize() - staying with old\n";
      return solverOK;
    }
  }

  return this->LinearSOE::setSolver(newSolver);
}


int
PetscSOE::sendSelf(int cTag, Channel &theChannel)
{

  // check if already using this object
  bool found = false;

  for (int i = 0; i < numChannels; i++)
    if (theChannels[i] == &theChannel)
    {
      // sendID = i + 1;
      found = true;
    }

  // if new object, enlarge Channel pointers to hold new channel * & allocate new ID
  if (found == false)
  {
    int nextNumChannels = numChannels + 1;
    Channel **nextChannels = new Channel *[nextNumChannels];

    if (nextNumChannels == 0)
    {
      cerr << "DistributedBandGenLinSOE::sendSelf() - failed to allocate channel array of size: " <<
           nextNumChannels << endln;
      return -1;
    }

    for (int i = 0; i < numChannels; i++)
    {
      nextChannels[i] = theChannels[i];
    }

    nextChannels[numChannels] = &theChannel;
    numChannels = nextNumChannels;

    if (theChannels != 0)
    {
      delete [] theChannels;
    }

    theChannels = nextChannels;

    if (localCol != 0)
    {
      delete [] localCol;
    }

    localCol = new ID *[numChannels];

    if (localCol == 0)
    {
      cerr << "DistributedBandGenLinSOE::sendSelf() - failed to allocate id array of size: " <<
           nextNumChannels << endln;
      return -1;
    }

    for (int i = 0; i < numChannels; i++)
    {
      localCol[i] = 0;
    }

  }

  return 0;
}


int
PetscSOE::recvSelf(int cTag, Channel &theChannel,
                   FEM_ObjectBroker &theBroker)
{
  numChannels = 1;
  theChannels = new Channel *[1];
  theChannels[0] = &theChannel;

  localCol = new ID *[numChannels];

  for (int i = 0; i < numChannels; i++)
  {
    localCol[i] = 0;
  }

  return 0;
}

int
PetscSOE::setChannels(int nChannels, Channel **theC)
{
  numChannels = nChannels;

  if (theChannels != 0)
  {
    delete [] theChannels;
  }

  theChannels = new Channel *[numChannels];

  for (int i = 0; i < numChannels; i++)
  {
    theChannels[i] = theC[i];
  }


  localCol = new ID *[nChannels];

  for (int i = 0; i < numChannels; i++)
  {
    localCol[i] = 0;
  }

  return 0;
}




int
PetscSOE::addM(const Matrix &m, const ID &id, double fact)
{

  isFactored = 0;

  // check for a quick return
  if (fact == 0.0)
  {
    return 0;
  }


  // check that m and id are of similar size
  int idSize = id.Size();

  if (idSize != m.noRows() && idSize != m.noCols())
  {
    std::cerr << "PetscSOE::addM() - Matrix and ID not of similar sizes\n";
    return -1;
  }

  int n = id.Size();
  int row;
  int col;
  double value;

  for (int i = 0; i < n; i++)
  {
    row = id(i);

    if (row >= 0)
    {
      for (int j = 0; j < n; j++)
      {
        col = id(j);

        if (col >= 0)
        {
          value = m(i, j) * fact;
          int ierr = MatSetValues(M, 1, &row, 1, &col, &value, ADD_VALUES);

          if (ierr)
          {
            cerr << processID << " " << row << " " << col << endln;
          }

          CHKERRQ(ierr);
        }
      }
    }
  }

  return 0;
}


int
PetscSOE::addK(const Matrix &m, const ID &id, double fact)
{

  isFactored = 0;

  // check for a quick return
  if (fact == 0.0)
  {
    return 0;
  }


  // check that m and id are of similar size
  int idSize = id.Size();

  if (idSize != m.noRows() && idSize != m.noCols())
  {
    std::cerr << "PetscSOE::addM() - Matrix and ID not of similar sizes\n";
    return -1;
  }

  int n = id.Size();
  int row;
  int col;
  double value;

  for (int i = 0; i < n; i++)
  {
    row = id(i);

    if (row >= 0)
    {
      for (int j = 0; j < n; j++)
      {
        col = id(j);

        if (col >= 0)
        {
          value = m(i, j) * fact;
          int ierr = MatSetValues(K, 1, &row, 1, &col, &value, ADD_VALUES);

          if (ierr)
          {
            cerr << processID << " " << row << " " << col << endln;
          }

          CHKERRQ(ierr);
        }
      }
    }
  }

  return 0;
}


void
PetscSOE::zeroM(void)
{
  isFactored = 0;
  MatZeroEntries(M);
}



void
PetscSOE::zeroK(void)
{
  isFactored = 0;
  MatZeroEntries(K);
}




int
PetscSOE::setEigenSize(Graph &theGraph)
{

  return -1;
}
