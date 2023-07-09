/* 
 * @author: George Petropoulos <gnp>
 *
 * @Description: This file contains the class definition for MPIDiagonalSOE 
 *               MPIDiagonalSOE is a subclass of LinearSOE. It stores a diagonal system
 *               MPIDiagonalSOE is a subclass of LinearSOE. It stores a diagonal system
 *               of equation, i.e. just the diagonal
 *               
 * @Date: 12/05
 *
 * Copyright: ALL RIGHTS RESERVED BY AUTHOR
 *
 */


#ifndef MPIDiagonalSOE_h
#define MPIDiagonalSOE_h

#include <mpi.h>

#include <map>
#include <LinearSOE.h>
#include <Vector.h>
#include <ID.h>
#include <stdlib.h>
#include <iostream>

class MPIDiagonalSolver;

class MPIDiagonalSOE : public LinearSOE
{
  public:
    MPIDiagonalSOE(MPIDiagonalSolver &theSolver);
    ~MPIDiagonalSOE();

    int getNumEqn(void) const;
    int setSize(Graph &theGraph);
    //void testSize(Graph &theGraph);
    int addA(const Matrix &, const ID &, double fact = 1.0);
    int addB(const Vector &, const ID &, double fact = 1.0);    
    int setB(const Vector &, double fact = 1.0);        
    
    void zeroA(void);
    void zeroB(void);

    const Vector& getpartofA(Vector& At, const ID& ids);

    void setX(int loc, double value);
    void setX(const Vector &x);
    
    const Vector &getX(void);
    const Vector &getB(void);
    double normRHS(void);

    int setDiagonalSolver(MPIDiagonalSolver &newSolver);    
    
    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);    
    int setChannels(int nChannels, Channel **theC);

    int setAnalysisModel(AnalysisModel &theModel);

    
    void intersections(ID& arrayA, ID& arrayB, int sizeA, int sizeB, int* numShared, int* sharedDOFs);
    void quickSort(ID& array, int array_size);
    void q_sort(ID& array, int left, int right);

    void dontUpdateA();

    friend class MPIDiagonalSolver;
    
  protected:
    
  private:
    int actualNeighbors;
    //
    bool updateA;
    int maxNeighbors;
    int* myNeighborsSizes; 
    int* myNeighbors;
    int* posLocKey;
    int size;
    double *A, *B, *X , *sharedA, *sharedB, *dataShared, *maxSharedA, *maxSharedB;
    Vector* vectX;
    Vector* vectB;
    Vector* vectA;
    bool isAfactored;

    int *myDOFsArray, *myDOFsSharedArray, *maxDOFsSharedArray;

    std::map<int,int*> myActualNeighborsSharedDOFs;
    std::map<int,double*> myActualNeighborsSharedBs;
    std::map<int,double*> myActualNeighborsBsToSend;
    int processID;
    int numProcesses;
    int numChannels;
    Channel **theChannels;
    ID **localCol;

    ID myDOFs;
    ID myDOFsShared;
    int numShared,maxShared;
    AnalysisModel *theModel;
};


#endif



