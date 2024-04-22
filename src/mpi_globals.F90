!!--------------------------------------------------------------------------!
!! libNEGF: a general library for Non-Equilibrium Green's functions.        !
!! Copyright (C) 2012                                                       !
!!                                                                          !
!! This file is part of libNEGF: a library for                              !
!! Non Equilibrium Green's Function calculation                             !
!!                                                                          !
!! Developers: Alessandro Pecchia, Gabriele Penazzi                         !
!! Former Conctributors: Luca Latessa, Aldo Di Carlo                        !
!!                                                                          !
!! libNEGF is free software: you can redistribute it and/or modify          !
!! it under the terms of the GNU Lesse General Public License as published  !
!! by the Free Software Foundation, either version 3 of the License, or     !
!! (at your option) any later version.                                      !
!!                                                                          !
!!  You should have received a copy of the GNU Lesser General Public        !
!!  License along with libNEGF.  If not, see                                !
!!  <http://www.gnu.org/licenses/>.                                         !
!!--------------------------------------------------------------------------!


module mpi_globals

#:if defined("MPI")
  use mpi
  use omp_lib, only : omp_get_max_threads, omp_get_num_threads
  use libmpifx_module!, only : mpifx_comm
  use ln_precision
  private
#:endif

  integer, public ::  numprocs = 1
  integer, public ::  id = 0
  logical, public ::  id0 = .true.

#:if defined("MPI")
  public :: negf_mpi_init
  public :: negf_cart_init
  public :: check_cart_comm
  public :: test_mpifx_reduce
  public :: test_bare_reduce

  contains

    subroutine negf_mpi_init(energyComm, ioProc)
      type(mpifx_comm) :: energyComm
      logical, optional :: ioProc

      id =energyComm%rank
      numprocs = energyComm%size

      if (present(ioProc)) then
        id0 = ioProc
      else
        id0 = (id == 0)
      end if

    end subroutine negf_mpi_init

    ! Initialize a 2D cartesian grid
    !
    ! Order: dim 1: k; dim 2: E
    ! It must be periodic in k for all-to-all communications
    ! CAVEAT:
    ! All processes MUST have the same number of points in K and E
    ! For E it is used to compute where E +/- wq are located
    ! For K it is used to compute where another q is placed
    !
    subroutine negf_cart_init(inComm, nk, cartComm, energyComm, kComm, bareCartComm, barekComm)
      !> Input communicator
      type(mpifx_comm), intent(in) :: inComm
      !> Number of processors for k
      integer, intent(in) :: nk
      !> Output 2D cartesian communicator
      type(mpifx_comm) :: cartComm
      !> Output communicator for the energy sub-grid
      type(mpifx_comm), intent(out) :: energyComm
      !> Output communicator for the k sub-grid
      type(mpifx_comm), intent(out) :: kComm

      !> Output communicators of type int for TiberCAD
      integer, intent(out), optional :: bareCartComm, barekComm

      integer :: outComm
      integer :: ndims = 2
      integer :: dims(2)
      logical :: periods(2) = .false.
      logical :: remain_dims(2)
      integer :: nE
      logical :: reorder = .true.
      integer :: mpierr

      if (mod(inComm%size,nk) /=0 ) then
        stop "Error in cart_init: cannot build a 2D cartesian grid with incompatible sizes"
      end if

      !call check_omp_mpi(inComm, mpierr)

      nE = inComm%size/nk
      dims(1)=nk; dims(2)=nE
      periods(1) = .true.

      call MPI_CART_CREATE(inComm%id, ndims, dims, periods, reorder, outComm, mpierr)
      call cartComm%init(outComm, mpierr)
      ! Global master id=0 node as writing node
      id0 = (cartComm%rank == 0)
      if (present(bareCartComm)) bareCartComm = outComm

      ! Extract sub-communicators
      remain_dims(:) = [.false., .true.]
      call MPI_CART_SUB(cartComm%id, remain_dims, outComm, mpierr)
      call energyComm%init(outComm, mpierr)
      id = energyComm%rank
      numprocs = energyComm%size

      remain_dims(:) = [.true., .false.]
      call MPI_CART_SUB(cartComm%id, remain_dims, outComm, mpierr)
      call kComm%init(outComm, mpierr)
      if (present(barekComm)) barekComm = outComm

      ! print*, "DEBUG: inside negf_cart_init: energyComm: ", energyComm
      ! print*, "DEBUG: inside negf_cart_init: kComm: ", kComm
      ! print*, "DEBUG: inside negf_cart_init: cartComm: ", cartComm

    end subroutine negf_cart_init

    subroutine check_cart_comm(cartComm, mpierror)
      !> Input 2d cartesian communicator
      type(mpifx_comm), intent(in) :: cartComm
      !> output error
      integer, intent(out) :: mpierror

      integer :: coords(2)

      mpierror = 0

      call MPI_Cart_coords(cartComm%id, 0, 2, coords, mpierror)

    end subroutine check_cart_comm

   ! subroutine check_omp_mpi(comm, mpierror)
   !   type(mpifx_comm), intent(in) :: comm
   !   integer, intent(out) :: mpierror

   !   integer :: newcomm, info, nprocs_shared, phys_cores, num_threads

   !   call MPI_COMM_SPLIT_TYPE(comm%id, MPI_COMM_TYPE_SHARED, 1, info, newcomm, mpierror)
   !   if (mpierror /= 0) then
   !      stop "ERROR in MPI_COMM_SPLIT_TYPE"
   !   end if
   !   call MPI_COMM_SIZE(newcomm, nprocs_shared, mpierror)

   !   phys_cores = get_num_cores()
   !   num_threads = omp_get_max_threads()

   !   if (comm%rank==0) then
   !     print*, "Number of physical cores:", phys_cores
   !     print*, "Number of processors on same shared memory:",nprocs_shared
   !     print*, "Maximum number of OMP threads:",num_threads
   !   end if
   !   if (num_threads*nprocs_shared > phys_cores) then
   !     call omp_set_max_threads(phys_cores/nprocs_shared)
   !   end if

   ! end subroutine check_omp_mpi

    subroutine test_mpifx_reduce(inComm, name)
      type(mpifx_comm), intent(in) :: inComm
      character(len=*), intent(in) :: name

      real(dp) :: test_array(5)
      
      test_array = 1.0_dp

      print*, ""
      print*, "TEST REDUCE (MPIFX), communicator:", name
      print*, name, "%id", inComm%id
      print*, name, "%size", inComm%size
      print*, name, "%rank", inComm%rank
      print*, name, "%leadrank", inComm%leadrank
      print*, name, "%lead", inComm%lead  
      print*, "Array before calling reduce:"
      print*, test_array
      print*, "Calling reduce..."
      call mpifx_reduceip(inComm, test_array, MPI_SUM)
      print*, "Array after calling reduce:"
      print*, test_array
      print*, "END TEST REDUCE"
      print*, ""

    end subroutine test_mpifx_reduce


    subroutine test_bare_reduce(inComm, name)
      integer, intent(in) :: inComm
      character(len=*), intent(in) :: name

      integer :: comm_size, rank, root0, error0, error
      real(dp) :: test_array(5)
      test_array = 1.0_dp

      root0 = 0

      call mpi_comm_size(inComm, comm_size, error0)
      call handle_errorflag(error0, "mpi_comm_size() in mpifx_comm_init()", error)
      if (error0 /= 0) then
        return
      end if
      call mpi_comm_rank(inComm, rank, error0)
      call handle_errorflag(error0, "mpi_comm_rank() in mpifx_comm_init()", error)
      if (error0 /= 0) then
        return
      end if

      print*, ""
      print*, "TEST REDUCE (PURE MPI), communicator: ", name
      print*, name, "%id", inComm
      print*, name, "%size", comm_size
      print*, name, "%rank", rank
      print*, "Array before calling reduce:"
      print*, test_array
      print*, "Calling reduce..."
      call mpi_reduce(MPI_IN_PLACE, test_array, size(test_array), MPI_DOUBLE_PRECISION, MPI_SUM, &
      &root0, inComm, error0)
      print*, "Array after calling reduce:"
      print*, test_array
      print*, "END TEST REDUCE"
      print*, ""
    
    end subroutine test_bare_reduce

    
    subroutine handle_errorflag(error0, msg, error)

      !> Error flag as returned by some routine.
      integer, intent(in) :: error0
  
      !>  Msg to print out, if program is stopped.
      character(*), intent(in) :: msg
  
      !> Optional error flag.
      !!
      !! If present, error0 is passed to it, otherwise if error0 was not zero, the
      !! error message in msg is printed and the program is stopped.
      !!
      integer, intent(out), optional :: error
  
      integer :: aborterror
  
      if (present(error)) then
        error = error0
      elseif (error0 /= 0) then
        write(*, "(A)") "Operation failed!"
        write(*, "(A)") msg
        write(*, "(A,I0)") "Error: ", error0
        call mpi_abort(MPI_COMM_WORLD, MPIFX_UNHANDLED_ERROR, aborterror)
        if (aborterror /= 0) then
          write(*, "(A)") "Stopping code with 'mpi_abort' did not succeed, trying 'stop' instead"
          stop 1
        end if
      end if
  
    end subroutine handle_errorflag

#:endif

end module mpi_globals
