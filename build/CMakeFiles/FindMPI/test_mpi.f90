      program hello
      use mpi_f08
      implicit none
      integer(kind=MPI_INTEGER_KIND) ierror
      call MPI_INIT(ierror)
      call MPI_FINALIZE(ierror)
      end program
