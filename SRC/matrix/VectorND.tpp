
  int
  addVector(const T thisFact, const Vector &other, const T otherFact) {
    if (otherFact == 0.0 && thisFact == 1.0)
      return 0; 

    else if (thisFact == 1.0) {
      // want: this += other * otherFact
      double *dataPtr = values;
      double *otherDataPtr = other.theData;
      if (otherFact == 1.0) { 
        // no point doing a multiplication if otherFact == 1.0
        for (int i=0; i<N; i++) 
          *dataPtr++ += *otherDataPtr++;
      } else if (otherFact == -1.0) { 
        // no point doing a multiplication if otherFact == 1.0
        for (int i=0; i<N; i++) 
          *dataPtr++ -= *otherDataPtr++;
      } else 
        for (int i=0; i<N; i++) 
          *dataPtr++ += *otherDataPtr++ * otherFact;

    } else if (thisFact == 0.0) {
        // want: this = other * otherFact
        double *dataPtr = values;
        double *otherDataPtr = other.theData;
        if (otherFact == 1.0) { // no point doing a multiplication if otherFact == 1.0
          for (int i=0; i<N; i++) 
            *dataPtr++ = *otherDataPtr++;
        } else if (otherFact == -1.0) { // no point doing a multiplication if otherFact == 1.0
          for (int i=0; i<N; i++) 
            *dataPtr++ = -(*otherDataPtr++);
        } else 
          for (int i=0; i<N; i++) 
            *dataPtr++ = *otherDataPtr++ * otherFact;
    } else {
      // want: this = this * thisFact + other * otherFact
      double *dataPtr = values;
      double *otherDataPtr = other.theData;
      if (otherFact == 1.0) { // no point doing a multiplication if otherFact == 1.0
        for (int i=0; i<N; i++) {
          double value = *dataPtr * thisFact + *otherDataPtr++;
          *dataPtr++ = value;
        }
      } else if (otherFact == -1.0) { // no point doing a multiplication if otherFact == 1.0
        for (int i=0; i<N; i++) {
          double value = *dataPtr * thisFact - *otherDataPtr++;
          *dataPtr++ = value;
        }
      } else 
        for (int i=0; i<N; i++) {
          double value = *dataPtr * thisFact + *otherDataPtr++ * otherFact;
          *dataPtr++ = value;
      }
    }

    // successfull
    return 0;
  }

  int
  addVector(const T thisFact, const VectorND<N> &other, const T otherFact) {
    if (otherFact == 0.0 && thisFact == 1.0)
      return 0; 

    else if (thisFact == 1.0) {
      // want: this += other * otherFact
      double *dataPtr = values;
      const double * otherDataPtr = other.values;
      if (otherFact == 1.0) { // no point doing a multiplication if otherFact == 1.0
        for (int i=0; i<N; i++) 
          *dataPtr++ += *otherDataPtr++;
      } else if (otherFact == -1.0) { // no point doing a multiplication if otherFact == 1.0
        for (int i=0; i<N; i++) 
          *dataPtr++ -= *otherDataPtr++;
      } else 
        for (int i=0; i<N; i++) 
          *dataPtr++ += *otherDataPtr++ * otherFact;

    } else if (thisFact == 0.0) {
        // want: this = other * otherFact
        double *dataPtr = values;
        const double *otherDataPtr = other.values;
        if (otherFact == 1.0) {
          // no point doing a multiplication if otherFact == 1.0
          for (int i=0; i<N; i++) 
            *dataPtr++ = *otherDataPtr++;
        } else if (otherFact == -1.0) {
          // no point doing a multiplication if otherFact == 1.0
          for (int i=0; i<N; i++) 
            *dataPtr++ = -(*otherDataPtr++);
        } else 
          for (int i=0; i<N; i++) 
            *dataPtr++ = *otherDataPtr++ * otherFact;
    } else {
      // want: this = this * thisFact + other * otherFact
      double *dataPtr = values;
      const double *otherDataPtr = other.values;
      if (otherFact == 1.0) { // no point doing a multiplication if otherFact == 1.0
        for (int i=0; i<N; i++) {
          double value = *dataPtr * thisFact + *otherDataPtr++;
          *dataPtr++ = value;
        }
      } else if (otherFact == -1.0) { // no point doing a multiplication if otherFact == 1.0
        for (int i=0; i<N; i++) {
          double value = *dataPtr * thisFact - *otherDataPtr++;
          *dataPtr++ = value;
        }
      } else 
        for (int i=0; i<N; i++) {
          double value = *dataPtr * thisFact + *otherDataPtr++ * otherFact;
          *dataPtr++ = value;
      }
    }

    // successfull
    return 0;
  }


  template <int NC>
  inline int
  addMatrixVector(double thisFact, const MatrixND<N, NC, double> &m, const Vector &v, double otherFact)
  {
    // check the sizes are compatable
    assert(NC == v.sz);

    // see if quick return
    if (thisFact == 1.0 && otherFact == 0.0)
      return 0;

    else {
      int incr = 1,
             i = N,
             n = NC;
       DGEMV("N", &i, &n,
             &otherFact,
             &m.values[0][0], &i,
             v.theData, &incr,
             &thisFact,
             values,   &incr);
      // successfull
      return 0;
    } 
  }

  template <int NR>
  inline int
  addMatrixTransposeVector(double thisFact, const MatrixND<NR, N, double> &m, const Vector &v, double otherFact)
  {
    // check the sizes are compatable
    assert(NR == v.sz);

    // see if quick return
    if (thisFact == 1.0 && otherFact == 0.0)
      return 0;

    else {
      int incr = 1,
             i = NR,
             n = N;
      DGEMV("T", &i, &n,
            &otherFact,
            &m.values[0][0], &i,
            v.theData, &incr,
            &thisFact,
            values,   &incr);
      // successfull
      return 0;
    } 
  }



  inline int
  addMatrixVector(const double thisFact, const Matrix &m, const Vector &v, const double otherFact)
  {
    // check the sizes are compatable
    assert(N == m.noRows());
    assert(m.noCols() == v.sz);

    // see if quick return
    if (thisFact == 1.0 && otherFact == 0.0)
      return 0;

#ifdef VECTOR_BLAS
    else if (v.sz > 10) {
      int incr = 1,
             i = m.numRows,
             n = m.numCols;
      return
        DGEMV("N", &i, &n,
              &otherFact,
              m.data, &i,
              v.theData, &incr,
              &thisFact,
              values,   &incr);
    }
#endif

    else if (thisFact == 1.0) {

      // want: this += m * v * otherFact
      if (otherFact == 1.0) { // no point doing multiplication if otherFact = 1.0
        int otherSize = v.sz;
        double *matrixDataPtr = m.data;
        const double *otherDataPtr = v.theData;
        for (int i=0; i<otherSize; i++) {
          double otherData = *otherDataPtr++;
          for (int j=0; j < N; j++)
            values[j] += *matrixDataPtr++ * otherData;
        }
      } 
      else if (otherFact == -1.0) { // no point doing multiplication if otherFact = -1.0
        int otherSize = v.sz;
        double *matrixDataPtr = m.data;
        const double *otherDataPtr = v.theData;
        for (int i=0; i<otherSize; i++) {
          double otherData = *otherDataPtr++;
          for (int j=0; j < N; j++)
            values[j] -= *matrixDataPtr++ * otherData;
        }
      } 
      else { // have to do the multiplication
        int otherSize = v.sz;
        double *matrixDataPtr = m.data;
        double *otherDataPtr = v.theData;
        for (int i=0; i<otherSize; i++) {
          double otherData = *otherDataPtr++ * otherFact;
          for (int j=0; j < N; j++)
            values[j] += *matrixDataPtr++ * otherData;
        }
      }
    }

    else if (thisFact == 0.0) {
      
      // want: this = m * v * otherFact
      for (int i=0; i < N; i++)
        values[i] = 0.0;

      if (otherFact == 1.0) { 
        // avoid multiplication when otherFact = 1.0
        int otherSize = v.sz;
        const double *matrixDataPtr = m.data;
        const double *otherDataPtr = v.theData;
        for (int i=0; i<otherSize; i++) {
          const double otherData = *otherDataPtr++;
          for (int j=0; j < N; j++)
            values[j] += *matrixDataPtr++ * otherData;
        }
      } 

      else if (otherFact == -1.0) { 
        // avoid multiplication when otherFact = -1.0
        int otherSize = v.sz;
        double *matrixDataPtr = m.data;
        const double *otherDataPtr = v.theData;
        for (int i=0; i<otherSize; i++) {
          const double otherData = *otherDataPtr++;
          for (int j=0; j < N; j++)
            values[j] -= *matrixDataPtr++ * otherData;
        }

      } else {
        int otherSize = v.sz;
        double *matrixDataPtr = m.data;
        double *otherDataPtr = v.theData;
        for (int i=0; i<otherSize; i++) {
          double otherData = *otherDataPtr++ * otherFact;
          for (int j=0; j < N; j++)
            values[j] += *matrixDataPtr++ * otherData;
        }
      }
    }

    else {

      // want: this = this * thisFact + m * v * otherFact
      for (int i=0; i<N; i++)
        values[i] *= thisFact;

      if (otherFact == 1.0) { // no point doing multiplication if otherFact = 1.0
        int otherSize = v.sz;
        double *matrixDataPtr = m.data;
        double *otherDataPtr = v.theData;
        for (int i=0; i<otherSize; i++) {
          double otherData = *otherDataPtr++;
          for (int j=0; j < N; j++)
            values[j] += *matrixDataPtr++ * otherData;
        }
      } else if (otherFact == -1.0) { // no point doing multiplication if otherFact = 1.0
        int otherSize = v.sz;
        double *matrixDataPtr = m.data;
        double *otherDataPtr = v.theData;
        for (int i=0; i<otherSize; i++) {
          double otherData = *otherDataPtr++;
          for (int j=0; j < N; j++)
            values[j] -= *matrixDataPtr++ * otherData;
        }
      } else {
        int otherSize = v.sz;
        double *matrixDataPtr = m.data;
        double *otherDataPtr = v.theData;
        for (int i=0; i<otherSize; i++) {
          double otherData = *otherDataPtr++ * otherFact;
          for (int j=0; j < N; j++)
            values[j] += *matrixDataPtr++ * otherData;
        }
      }
    }

    // successfull
    return 0;
  }

