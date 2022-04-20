module scba
  use ln_precision, only : dp
  use interactions
  use mat_def, only : z_CSR, z_DNS, create, destroy
  use sparsekit_drv, only : trace, clone

  type TScbaDriver
    !> Keep track of SCBA iteration
    integer :: scba_iter = 0

    !> SCBA Tolerance (Exact meaning may depend on model)
    real(dp) :: scba_tol = 1.0d-7

    !> SCBA error set for output
    real(dp) :: scba_err = 0.0_dp

    !> Holds the value of the previous Gn
    type(z_CSR), allocatable :: Mat_old

    !> Holds the value of the previous current
    real(dp), dimension(:), allocatable :: J_old

    !> internal status
    logical :: converged = .false.

    !> internal status
    logical :: do_write = .false.

    contains
    procedure :: init => scba_init
    procedure :: destroy => scba_destroy
    procedure :: set_scba_iter
    procedure :: is_converged
    procedure :: check_Mat_convergence
    procedure :: check_J_convergence
  end type TScbaDriver

  contains

  !> Initialize the SCBA Driver
  subroutine scba_init(this, tol, dowrite)
    class(TScbaDriver) :: this
    real(dp), intent(in) :: tol
    logical, intent(in) :: dowrite
    this%scba_tol = tol
    this%do_write = dowrite
    allocate(this%Mat_old)
  end subroutine scba_init

  !> Release internal space
  subroutine scba_destroy(this)
    class(TScbaDriver) :: this
    integer :: ii
    if (allocated(this%Mat_old%nzval)) call destroy(this%Mat_old)
    if (allocated(this%Mat_old)) deallocate(this%Mat_old)
    if (allocated(this%J_old)) deallocate(this%J_old)
  end subroutine scba_destroy


  !> Sets the scba_iteration in all interactions
  subroutine set_scba_iter(this, iter, interactArray)
    class(TScbaDriver) :: this
    integer, intent(in) :: iter
    type(TInteractionArray) :: interactArray(:)

    integer :: ii
    this%scba_iter = iter
    do ii = 1, size(interactArray)
      call interactArray(ii)%inter%set_scba_iter(iter)
    end do
  end subroutine set_scba_iter


  !> Check convergence on Gn matrix blocks
  subroutine check_Mat_convergence(this, Mat)
    class(TScbaDriver) :: this
    type(z_CSR), intent(in) :: Mat

    real(dp) :: maxerr
    integer :: ii

    if (.not.allocated(this%Mat_old)) then
      stop 'ERROR: TScbaDriver must be initialized first'
    end if
    if (.not.allocated(this%Mat_old%nzval)) then
      call clone(Mat, this%Mat_old)
      return
    end if

    maxerr = maxval(abs(Mat%nzval-this%Mat_old%nzval))
    if (maxerr < this%scba_tol) then
      this%converged = .true.
      if (this%do_write) then
         write(*,*) "SCBA loop converged in",this%scba_iter,&
                  & " iterations with error", maxerr
      end if
    else
      if (this%do_write) then
         write(*,*) "SCBA iteration",this%scba_iter," error", maxerr
      end if
    end if
    call destroy(this%Mat_old)
    call clone(Mat, this%Mat_old)
  end subroutine check_Mat_convergence

  !> Check convergence on the 1st layer current
  subroutine check_J_convergence(this, J)
    class(TScbaDriver) :: this
    real(dp), dimension(:), intent(in) :: J

    real(dp) :: error

    if (.not.allocated(this%J_old)) then
      allocate(this%J_old, source=J)
      return
    end if

    error = maxval(abs(this%J_old - J))
    if (error < this%scba_tol) then
      this%converged = .true.
      if (this%do_write) then
         write(*,*) "SCBA loop converged in",this%scba_iter,&
                  & " iterations with error", maxerr
      end if
    else
      if (this%do_write) then
         write(*,*) "SCBA iteration",this%scba_iter," error", maxerr
      end if
    end if
    this%J_old = J
  end subroutine check_J_convergence

  function is_converged(this) result (conv)
    class(TScbaDriver) :: this
    logical :: conv

    conv = this%converged
  end function is_converged

end module scba

