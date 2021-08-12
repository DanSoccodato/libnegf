module iterative_gpu
  use ln_precision
  use ln_constants, only : pi
  use ln_allocation
  use mat_def
  use ln_structure, only : TStruct_Info
  use lib_param, only : MAXNCONT, Tnegf, intarray, cusolverDnHandle, cublasHandle
  use mpi_globals, only : id, numprocs, id0
  use sparsekit_drv, only: csr2blk_sod
  use clock
  use cudautils
  use, intrinsic :: ieee_arithmetic

  implicit none
  private 

  public :: calculate_gsmr_blocks
  public :: calculate_gsml_blocks
  public :: calculate_Gr_tridiag_blocks
  public :: calculate_Gn_tridiag_blocks
  public :: calculate_single_transmission_2_contacts
  public :: calculate_single_transmission_N_contacts
  public :: build_ESH_onGPU 

  interface calculate_gsmr_blocks
     module procedure  calculate_gsmr_blocks_sp
     module procedure  calculate_gsmr_blocks_dp
  end interface calculate_gsmr_blocks

  interface calculate_gsml_blocks
     module procedure  calculate_gsml_blocks_sp
     module procedure  calculate_gsml_blocks_dp
  end interface calculate_gsml_blocks

  interface calculate_Gr_tridiag_blocks
     module procedure calculate_Gr_tridiag_blocks_sp
     module procedure calculate_Gr_tridiag_blocks_dp
  end interface calculate_Gr_tridiag_blocks

  interface calculate_Gn_tridiag_blocks
     module procedure calculate_Gn_tridiag_blocks_sp
     module procedure calculate_Gn_tridiag_blocks_dp
  end interface calculate_Gn_tridiag_blocks

  interface calculate_single_transmission_2_contacts
     module procedure calculate_single_transmission_2_contacts_sp
     module procedure calculate_single_transmission_2_contacts_dp
  end interface calculate_single_transmission_2_contacts

  interface calculate_single_transmission_N_contacts
     module procedure calculate_single_transmission_N_contacts_sp
     module procedure calculate_single_transmission_N_contacts_dp
  end interface calculate_single_transmission_N_contacts

  interface build_ESH_onGPU 
!     module procedure build_ESH_onGPU_sp 
     module procedure build_ESH_onGPU_dp
  end interface build_ESH_onGPU 
  
  interface get_tun_mask 
     module procedure get_tun_mask_sp 
     module procedure get_tun_mask_dp
  end interface get_tun_mask 

contains

  subroutine calculate_gsmr_blocks_sp(negf,ESH,sbl,ebl,gsmr,keepall)

    !In/Out
    type(c_DNS), dimension(:), intent(inout) :: gsmr
    type(Tnegf), intent(in) :: negf
    type(c_DNS), dimension(:,:), intent(in) :: ESH
    integer, intent(in) :: sbl,ebl                 ! start block, end block
    logical, intent(in), optional :: keepall

    !Work
    type(CublasHandle) :: hh
    type(CusolverDnHandle) :: hhsol
    complex(sp), parameter :: one = (1.0_sp, 0.0_sp)
    complex(sp), parameter :: mone = (-1.0_sp, 0.0_sp)
    complex(sp), parameter :: zero = (0.0_sp, 0.0_sp)
    type(c_DNS) :: work1, work2
    integer :: nrow, ncol
    integer :: i, nbl, istat
    real(sp) :: sum1
    logical :: keep

    if (sbl.lt.ebl) return

    nbl = size(ESH,1)
    if (nbl.eq.1) return

    keep = .true.
    if (present(keepall)) then
       keep = keepall
    end if

    nrow=ESH(sbl,sbl)%nrow

    hh = negf%hcublas
    hhsol = negf%hcusolver

    call createAll(gsmr(sbl),nrow,nrow)
    call copyToGPU(ESH(sbl,sbl))
    call inverse_gpu(hh, hhsol, ESH(sbl,sbl), gsmr(sbl), istat)
    call deleteGPU(ESH(sbl,sbl))

    do i=sbl-1,ebl,-1

       call createAll(work1, ESH(i,i+1)%nrow, gsmr(i+1)%ncol)    
       call copyToGPU(ESH(i,i+1))
       call matmul_gpu(hh, one, ESH(i,i+1), gsmr(i+1), zero, work1)
       call deleteGPU(ESH(i,i+1))

       if (.not.keep) then
          call destroyAll(gsmr(i+1))
       end if

       call createAll(work2, work1%nrow, ESH(i+1,i)%ncol)
       call copyToGPU(ESH(i+1,i))
       call matmul_gpu(hh, one, work1, ESH(i+1,i), zero, work2)
       call deleteGPU(ESH(i+1,i))
       call destroyAll(work1)

       call createAll(work1,  ESH(i,i)%nrow, ESH(i,i)%ncol)
       call copyToGPU(ESH(i,i))
       call matsum_gpu(hh, one, ESH(i,i), mone, work2, work1)
       call deleteGPU(ESH(i,i))
       call destroyAll(work2)

       call createAll(gsmr(i), work1%nrow, work1%ncol)
       call inverse_gpu(hh, hhsol, work1, gsmr(i), istat)    
       call destroyAll(work1)

    end do

  end subroutine calculate_gsmr_blocks_sp

  subroutine calculate_gsmr_blocks_dp(negf,ESH,sbl,ebl,gsmr,keepall)

    !In/Out
    type(z_DNS), dimension(:), intent(inout) :: gsmr
    type(Tnegf), intent(in) :: negf
    type(z_DNS), dimension(:,:), intent(in) :: ESH
    integer, intent(in) :: sbl,ebl                 ! start block, end block
    logical, intent(in), optional :: keepall

    !Work
    type(CublasHandle) :: hh
    type(CusolverDnHandle) :: hhsol
    complex(dp), parameter :: one = (1.0_dp, 0.0_dp)
    complex(dp), parameter :: mone = (-1.0_dp, 0.0_dp)
    complex(dp), parameter :: zero = (0.0_dp, 0.0_dp)
    type(z_DNS) :: work1, work2
    integer :: nrow, ncol
    integer :: i, nbl, istat
    real(dp) :: sum1
    logical :: keep

    if (sbl.lt.ebl) return

    nbl = size(ESH,1)
    if (nbl.eq.1) return

    keep = .true.
    if (present(keepall)) then
       keep = keepall
    end if


    nrow=ESH(sbl,sbl)%nrow

    hh = negf%hcublas
    hhsol = negf%hcusolver

    call createAll(gsmr(sbl),nrow,nrow)
    call copyToGPU(ESH(sbl,sbl))
    call inverse_gpu(hh, hhsol, ESH(sbl,sbl), gsmr(sbl), istat)
    call deleteGPU(ESH(sbl,sbl))

    do i=sbl-1,ebl,-1

       call createAll(work1, ESH(i,i+1)%nrow, gsmr(i+1)%ncol)    
       call copyToGPU(ESH(i,i+1))
       call matmul_gpu(hh, one, ESH(i,i+1), gsmr(i+1), zero, work1)
       call deleteGPU(ESH(i,i+1))

       if (.not.keep) then
          call destroyAll(gsmr(i+1))
       end if

       call createAll(work2, work1%nrow, ESH(i+1,i)%ncol)
       call copyToGPU(ESH(i+1,i))
       call matmul_gpu(hh, one, work1, ESH(i+1,i), zero, work2)
       call deleteGPU(ESH(i+1,i))
       call destroyAll(work1)

       call createAll(work1,  ESH(i,i)%nrow, ESH(i,i)%ncol)
       call copyToGPU(ESH(i,i))
       call matsum_gpu(hh, one, ESH(i,i), mone, work2, work1)
       call deleteGPU(ESH(i,i))
       call destroyAll(work2)

       call createAll(gsmr(i), work1%nrow, work1%ncol)
       call inverse_gpu(hh, hhsol, work1, gsmr(i), istat)    
       call destroyAll(work1)

    end do

  end subroutine calculate_gsmr_blocks_dp

  subroutine calculate_gsml_blocks_sp(negf,ESH,sbl,ebl,gsml)

    !In/Out
    type(c_DNS), dimension(:), intent(inout) :: gsml
    type(Tnegf), intent(in) :: negf
    type(c_DNS), dimension(:,:), intent(in) :: ESH
    integer, intent(in) :: sbl,ebl                 ! start block, end block

    !Work
    type(CublasHandle) :: hh
    type(CusolverDnHandle) :: hhsol
    complex(sp), parameter :: one = (1.0_sp, 0.0_sp)
    complex(sp), parameter :: mone = (-1.0_sp, 0.0_sp)
    complex(sp), parameter :: zero = (0.0_sp, 0.0_sp)
    type(c_DNS) :: work1, work2
    integer :: nrow, ncol
    integer :: i, nbl, istat
    real(sp) :: sum1

    if (sbl.gt.ebl) return

    nbl = size(ESH,1)
    if (nbl.eq.1) return

    nrow=ESH(sbl,sbl)%nrow

    hh = negf%hcublas
    hhsol = negf%hcusolver

    call createAll(gsml(sbl),nrow,nrow)
    call copyToGPU(ESH(sbl,sbl))
    call inverse_gpu(hh, hhsol, ESH(sbl,sbl), gsml(sbl), istat)
    call deleteGPU(ESH(sbl,sbl))

    do i=sbl+1,ebl

       call createAll(work1, ESH(i,i-1)%nrow, gsml(i-1)%ncol)    
       call copyToGPU(ESH(i,i-1))
       call matmul_gpu(hh, one, ESH(i,i-1), gsml(i-1), zero, work1)
       call deleteGPU(ESH(i,i-1))

       call createAll(work2, work1%nrow, ESH(i-1,i)%ncol)
       call copyToGPU(ESH(i-1,i))
       call matmul_gpu(hh, one, work1, ESH(i-1,i), zero, work2)
       call deleteGPU(ESH(i-1,i))
       call destroyAll(work1)

       call createAll(work1,  ESH(i,i)%nrow, ESH(i,i)%ncol)
       call copyToGPU(ESH(i,i))
       call matsum_gpu(hh, one, ESH(i,i), mone, work2, work1)
       call deleteGPU(ESH(i,i))
       call destroyAll(work2)

       call createAll(gsml(i), work1%nrow, work1%ncol)
       call inverse_gpu(hh, hhsol, work1, gsml(i), istat)    
       call destroyAll(work1)

    end do

  end subroutine calculate_gsml_blocks_sp

  subroutine calculate_gsml_blocks_dp(negf,ESH,sbl,ebl,gsml)
    !In/Out
    type(z_DNS), dimension(:), intent(inout) :: gsml
    type(Tnegf), intent(in) :: negf
    type(z_DNS), dimension(:,:), intent(in) :: ESH
    integer, intent(in) :: sbl,ebl                 ! start block, end block

    !Work
    type(CublasHandle) :: hh
    type(CusolverDnHandle) :: hhsol
    complex(dp), parameter :: one = (1.0_dp, 0.0_dp)
    complex(dp), parameter :: mone = (-1.0_dp, 0.0_dp)
    complex(dp), parameter :: zero = (0.0_dp, 0.0_dp)
    type(z_DNS) :: work1, work2
    integer :: nrow, ncol
    integer :: i, nbl, istat
    real(dp) :: sum1

    if (sbl.gt.ebl) return

    nbl = size(ESH,1)
    if (nbl.eq.1) return

    nrow=ESH(sbl,sbl)%nrow

    hh = negf%hcublas
    hhsol = negf%hcusolver

    call createAll(gsml(sbl),nrow,nrow)
    call copyToGPU(ESH(sbl,sbl))
    call inverse_gpu(hh, hhsol, ESH(sbl,sbl), gsml(sbl), istat)
    call deleteGPU(ESH(sbl,sbl))

    do i=sbl+1,ebl

       call createAll(work1, ESH(i,i-1)%nrow, gsml(i-1)%ncol)    
       call copyToGPU(ESH(i,i-1))
       call matmul_gpu(hh, one, ESH(i,i-1), gsml(i-1), zero, work1)
       call deleteGPU(ESH(i,i-1))

       call createAll(work2, work1%nrow, ESH(i-1,i)%ncol)
       call copyToGPU(ESH(i-1,i))
       call matmul_gpu(hh, one, work1, ESH(i-1,i), zero, work2)
       call deleteGPU(ESH(i-1,i))
       call destroyAll(work1)

       call createAll(work1,  ESH(i,i)%nrow, ESH(i,i)%ncol)
       call copyToGPU(ESH(i,i))
       call matsum_gpu(hh, one, ESH(i,i), mone, work2, work1)
       call deleteGPU(ESH(i,i))
       call destroyAll(work2)

       call createAll(gsml(i), work1%nrow, work1%ncol)
       call inverse_gpu(hh, hhsol, work1, gsml(i), istat)    
       call destroyAll(work1)

    end do

  end subroutine calculate_gsml_blocks_dp

  subroutine calculate_Gr_tridiag_blocks_sp(negf,ESH,gsml,gsmr,Gr,sbl,ebl)
    !In/Out
    type(c_DNS), dimension(:,:), intent(inout) :: Gr
    type(c_DNS), dimension(:), intent(in) :: gsmr, gsml
    type(Tnegf), intent(in) :: negf
    type(c_DNS), dimension(:,:), intent(in) :: ESH
    integer, intent(in) :: sbl
    integer, intent(in), optional :: ebl

    !Work
    type(CublasHandle) :: hh
    type(CusolverDnHandle) ::hhsol
    integer :: i,nrow,nbl
    type(c_DNS), target :: work1, work2, work3
    complex(sp), parameter :: one = (1.0_sp, 0.0_sp)
    complex(sp), parameter :: mone = (-1.0_sp, 0.0_sp)
    complex(sp), parameter :: zero = (0.0_sp, 0.0_sp)
    integer :: istat

    nbl = size(ESH,1)
    hh = negf%hcublas
    hhsol = negf%hcusolver

    if (sbl.gt.nbl) return
    if (sbl.lt.1) return

    if (.not.present(ebl)) then
       if (nbl.eq.1) then
          call createAll(Gr(sbl,sbl), ESH(sbl,sbl)%nrow, ESH(sbl,sbl)%ncol)

          call copyToGPU(ESH(sbl,sbl))
          call inverse_gpu(hh, hhsol, ESH(sbl,sbl), Gr(sbl,sbl), istat)
          call deleteGPU(ESH(sbl,sbl))
       else
          call createAll(work1, ESH(sbl,sbl)%nrow, ESH(sbl,sbl)%ncol)

          call copyToGPU(ESH(sbl,sbl))
          call copy_mat_gpu(hh, ESH(sbl,sbl), work1)
          call deleteGPU(ESH(sbl,sbl))

          if (sbl+1.le.nbl) then
             call createAll(work2, ESH(sbl,sbl+1)%nrow, gsmr(sbl+1)%ncol)
             call copyToGPU(ESH(sbl,sbl+1))
             call matmul_gpu(hh, one, ESH(sbl,sbl+1), gsmr(sbl+1), zero, work2)
             call deleteGPU(ESH(sbl,sbl+1))

             call createAll(work3, work2%nrow, ESH(sbl+1,sbl)%ncol)
             call copyToGPU(ESH(sbl+1,sbl))
             call matmul_gpu(hh, one, work2, ESH(sbl+1,sbl), zero, work3)
             call deleteGPU(ESH(sbl+1,sbl))

             call copyToGPU(ESH(sbl,sbl))
             call matsum_gpu(hh, one, ESH(sbl,sbl), mone, work3, work1)
             call deleteGPU(ESH(sbl,sbl))

             call destroyAll(work2)
             call destroyAll(work3)
          end if
          if (sbl-1.ge.1) then
             call createAll(work2, ESH(sbl,sbl-1)%nrow, gsml(sbl-1)%ncol)
             call copyToGPU(ESH(sbl,sbl-1))
             call matmul_gpu(hh, one, ESH(sbl,sbl-1), gsml(sbl-1), zero, work2)
             call deleteGPU(ESH(sbl,sbl-1))

             call createAll(work3, work2%nrow, ESH(sbl-1,sbl)%ncol)
             call copyToGPU(ESH(sbl-1,sbl))
             call matmul_gpu(hh, one, work2, ESH(sbl-1,sbl), zero, work3)
             call deleteGPU(ESH(sbl-1,sbl))

             call copyToGPU(ESH(sbl,sbl))
             call matsum_gpu(hh, one, ESH(sbl,sbl), mone, work3, work1)
             call deleteGPU(ESH(sbl,sbl))

             call destroyAll(work2)
             call destroyAll(work3)
          end if

          call createAll(Gr(sbl,sbl), work1%nrow, work1%ncol)
          call inverse_gpu(hh, hhsol, work1, Gr(sbl,sbl), istat)
          call destroyAll(work1)
       endif
       return
    endif
    !***
    !Diagonal, Subdiagonal and Superdiagonal blocks
    !***
    if ((ebl.ge.sbl).and.(ebl.gt.1).and.(sbl.gt.1)) THEN
       do i=sbl,ebl,1
          call createAll(work1, gsmr(i)%nrow, ESH(i,i-1)%ncol)
          call createAll(Gr(i,i-1), work1%nrow, Gr(i-1,i-1)%ncol)
          
          call copyToGPU(ESH(i,i-1))
          call matmul_gpu(hh, one, gsmr(i), ESH(i,i-1), zero, work1)
          call deleteGPU(ESH(i,i-1))

          call matmul_gpu(hh, mone, work1, Gr(i-1,i-1), zero, Gr(i,i-1))

          call destroyAll(work1)
          call createAll(work2, ESH(i-1,i)%nrow, gsmr(i)%ncol)
          call createAll(Gr(i-1,i), Gr(i-1,i-1)%nrow, work2%ncol)

          call copyToGPU(ESH(i-1,i))
          call matmul_gpu(hh, one, ESH(i-1,i), gsmr(i), zero, work2)
          call deleteGPU(ESH(i-1,i))

          call matmul_gpu(hh, mone, Gr(i-1,i-1), work2, zero, Gr(i-1,i))

          call createAll(work1, Gr(i,i-1)%nrow, work2%ncol)
          call matmul_gpu(hh, mone, Gr(i,i-1), work2, zero, work1)

          call destroyAll(work2)
          call createAll(Gr(i,i), gsmr(i)%nrow, gsmr(i)%ncol)
          call matsum_gpu(hh, one, gsmr(i), one, work1, Gr(i,i))
          call destroyAll(work1)
       end do
    else
       do i=sbl,ebl,-1
          call createAll(work1, gsml(i)%nrow, ESH(i,i+1)%ncol)
          call createAll(Gr(i,i+1),work1%nrow,Gr(i+1,i+1)%ncol)

          call copyToGPU(ESH(i,i+1))
          call matmul_gpu(hh, one, gsml(i), ESH(i,i+1), zero, work1)
          call deleteGPU(ESH(i,i+1))

          call matmul_gpu(hh, mone, work1, Gr(i+1,i+1), zero, Gr(i,i+1))

          call destroyAll(work1)
          call createAll(work2, ESH(i+1,i)%nrow, gsml(i)%ncol)
          call createAll(Gr(i+1,i), Gr(i+1,i+1)%nrow, work2%ncol)

          call copyToGPU(ESH(i+1,i))
          call matmul_gpu(hh, one, ESH(i+1,i), gsml(i), zero, work2)
          call deleteGPU(ESH(i+1,i))

          call matmul_gpu(hh, mone, Gr(i+1,i+1), work2, zero, Gr(i+1,i))

          call createAll(work1, Gr(i,i+1)%nrow, work2%ncol)
          call matmul_gpu(hh, mone, Gr(i,i+1), work2, zero, work1)

          call destroyAll(work2)
          call createAll(Gr(i,i), gsml(i)%nrow, gsml(i)%ncol)
          call matsum_gpu(hh, one, gsml(i), one, work1, Gr(i,i))
          call destroyAll(work1)
       end do
    endif
           
  end subroutine calculate_Gr_tridiag_blocks_sp

  subroutine calculate_Gr_tridiag_blocks_dp(negf,ESH,gsml,gsmr,Gr,sbl,ebl)
    !In/Out
    type(z_DNS), dimension(:,:), intent(inout) :: Gr
    type(z_DNS), dimension(:), intent(in) :: gsmr, gsml
    type(Tnegf), intent(in) :: negf
    type(z_DNS), dimension(:,:), intent(in) :: ESH
    integer, intent(in) :: sbl
    integer, intent(in), optional :: ebl

    !Work
    type(CublasHandle) :: hh
    type(CusolverDnHandle) :: hhsol
    integer :: i, nrow, nbl, istat
    type(z_DNS), target :: work1, work2, work3
    complex(dp), parameter :: one = (1.0_dp, 0.0_dp)
    complex(dp), parameter :: mone = (-1.0_dp, 0.0_dp)
    complex(dp), parameter :: zero = (0.0_dp, 0.0_dp)
    real(dp) :: summ

    nbl = size(ESH,1)
    hh = negf%hcublas
    hhsol = negf%hcusolver

    if (sbl.gt.nbl) return
    if (sbl.lt.1) return

    if (.not.present(ebl)) then
       if (nbl.eq.1) then
          call createAll(Gr(sbl,sbl), ESH(sbl,sbl)%nrow, ESH(sbl,sbl)%ncol)

          call copyToGPU(ESH(sbl,sbl))
          call inverse_gpu(hh, hhsol, ESH(sbl,sbl), Gr(sbl,sbl), istat)
          call deleteGPU(ESH(sbl,sbl))
       else
          call createAll(work1, ESH(sbl,sbl)%nrow, ESH(sbl,sbl)%ncol)

          call copyToGPU(ESH(sbl,sbl))
          call copy_mat_gpu(hh, ESH(sbl,sbl), work1)
          call deleteGPU(ESH(sbl,sbl))

          if (sbl+1.le.nbl) then
             call createAll(work2, ESH(sbl,sbl+1)%nrow, gsmr(sbl+1)%ncol)
             call copyToGPU(ESH(sbl,sbl+1))
             call matmul_gpu(hh, one, ESH(sbl,sbl+1), gsmr(sbl+1), zero, work2)
             call deleteGPU(ESH(sbl,sbl+1))

             call createAll(work3, work2%nrow, ESH(sbl+1,sbl)%ncol)
             call copyToGPU(ESH(sbl+1,sbl))
             call matmul_gpu(hh, one, work2, ESH(sbl+1,sbl), zero, work3)
             call deleteGPU(ESH(sbl+1,sbl))

             call copyToGPU(ESH(sbl,sbl))
             call matsum_gpu(hh, one, ESH(sbl,sbl), mone, work3, work1)
             call deleteGPU(ESH(sbl,sbl))

             call destroyAll(work2)
             call destroyAll(work3)
          end if
          if (sbl-1.ge.1) then
             call createAll(work2, ESH(sbl,sbl-1)%nrow, gsml(sbl-1)%ncol)
             call copyToGPU(ESH(sbl,sbl-1))
             call matmul_gpu(hh, one, ESH(sbl,sbl-1), gsml(sbl-1), zero, work2)
             call deleteGPU(ESH(sbl,sbl-1))

             call createAll(work3, work2%nrow, ESH(sbl-1,sbl)%ncol)
             call copyToGPU(ESH(sbl-1,sbl))
             call matmul_gpu(hh, one, work2, ESH(sbl-1,sbl), zero, work3)
             call deleteGPU(ESH(sbl-1,sbl))

             call copyToGPU(ESH(sbl,sbl))
             call matsum_gpu(hh, one, ESH(sbl,sbl), mone, work3, work1)
             call deleteGPU(ESH(sbl,sbl))

             call destroyAll(work2)
             call destroyAll(work3)
          end if

          call createAll(Gr(sbl,sbl), work1%nrow, work1%ncol)
          call inverse_gpu(hh, hhsol, work1, Gr(sbl,sbl), istat)
          call destroyAll(work1)
       endif
       return
    endif
    !***
    !Diagonal, Subdiagonal and Superdiagonal blocks
    !***
    if ((ebl.ge.sbl).and.(ebl.gt.1).and.(sbl.gt.1)) THEN
       do i=sbl,ebl,1
          call createAll(work1, gsmr(i)%nrow, ESH(i,i-1)%ncol)
          call createAll(Gr(i,i-1), work1%nrow, Gr(i-1,i-1)%ncol)
          
          call copyToGPU(ESH(i,i-1))
          call matmul_gpu(hh, one, gsmr(i), ESH(i,i-1), zero, work1)
          call deleteGPU(ESH(i,i-1))

          call matmul_gpu(hh, mone, work1, Gr(i-1,i-1), zero, Gr(i,i-1))

          call destroyAll(work1)
          call createAll(work2, ESH(i-1,i)%nrow, gsmr(i)%ncol)
          call createAll(Gr(i-1,i), Gr(i-1,i-1)%nrow, work2%ncol)

          call copyToGPU(ESH(i-1,i))
          call matmul_gpu(hh, one, ESH(i-1,i), gsmr(i), zero, work2)
          call deleteGPU(ESH(i-1,i))

          call matmul_gpu(hh, mone, Gr(i-1,i-1), work2, zero, Gr(i-1,i))

          call createAll(work1, Gr(i,i-1)%nrow, work2%ncol)
          call matmul_gpu(hh, mone, Gr(i,i-1), work2, zero, work1)

          call destroyAll(work2)
          call createAll(Gr(i,i), gsmr(i)%nrow, gsmr(i)%ncol)
          call matsum_gpu(hh, one, gsmr(i), one, work1, Gr(i,i))
          call destroyAll(work1)
       end do
    else
       do i=sbl,ebl,-1
          call createAll(work1, gsml(i)%nrow, ESH(i,i+1)%ncol)
          call createAll(Gr(i,i+1),work1%nrow,Gr(i+1,i+1)%ncol)

          call copyToGPU(ESH(i,i+1))
          call matmul_gpu(hh, one, gsml(i), ESH(i,i+1), zero, work1)
          call deleteGPU(ESH(i,i+1))

          call matmul_gpu(hh, mone, work1, Gr(i+1,i+1), zero, Gr(i,i+1))

          call destroyAll(work1)
          call createAll(work2, ESH(i+1,i)%nrow, gsml(i)%ncol)
          call createAll(Gr(i+1,i), Gr(i+1,i+1)%nrow, work2%ncol)

          call copyToGPU(ESH(i+1,i))
          call matmul_gpu(hh, one, ESH(i+1,i), gsml(i), zero, work2)
          call deleteGPU(ESH(i+1,i))

          call matmul_gpu(hh, mone, Gr(i+1,i+1), work2, zero, Gr(i+1,i))

          call createAll(work1, Gr(i,i+1)%nrow, work2%ncol)
          call matmul_gpu(hh, mone, Gr(i,i+1), work2, zero, work1)

          call destroyAll(work2)
          call createAll(Gr(i,i), gsml(i)%nrow, gsml(i)%ncol)
          call matsum_gpu(hh, one, gsml(i), one, work1, Gr(i,i))
          call destroyAll(work1)
       end do
    endif
           
  end subroutine calculate_Gr_tridiag_blocks_dp

  subroutine calculate_Gn_tridiag_blocks_sp(negf, ESH, SelfEneR, frm, ref, struct, gsml, gsmr, Gr, Gn) 
    !In/Out
    type(Tnegf), intent(in) :: negf
    type(c_DNS), dimension(:), intent(in) :: gsmr, gsml
    type(c_DNS), dimension(:,:), intent(in) :: ESH, Gr
    type(c_DNS), dimension(:,:), intent(inout) :: Gn
    type(c_DNS), dimension(:), intent(in) :: SelfEneR
    type(Tstruct_info), intent(in) :: struct
    real(sp), dimension(:), intent(in) :: frm
    integer, intent(in) :: ref

    !Work
    type(CublasHandle) :: hh
    Type(c_DNS) :: Gam
    type(c_DNS) :: work1, work2, work3
    integer :: i,j
    integer :: cb, nbl, ncont, istat
    complex(sp) :: frmdiff
    complex(sp), parameter :: one = (1.0_sp, 0.0_sp)
    complex(sp), parameter :: mone = (-1.0_sp, 0.0_sp)
    complex(sp), parameter :: zero = (0.0_sp, 0.0_sp)
    real(sp) :: somma

    ncont = struct%num_conts
    nbl = struct%num_PLs
    hh = negf%hcublas

    do j=1,ncont

       if (j.NE.ref .AND. ABS(frm(j)-frm(ref)).GT.EPS) THEN

          cb=struct%cblk(j) ! block corresponding to contact j    

          !Creating Gam = i (SE -  SE^+)
          call createAll(Gam, SelfEneR(j)%nrow, SelfEneR(j)%ncol)
          call spectral_gpu(hh, SelfEneR(j), Gam)

          frmdiff = cmplx(frm(j)-frm(ref),0.0_dp,dp)

          ! Computation of Gn(cb,cb) = Gr(cb,cb) Gam(cb) Gr(cb,cb)^+
          call createAll(work1, Gam%nrow, Gr(cb,cb)%nrow) !it's Gr^+ --> I take nrow
          call matmul_gpu(hh, frmdiff, Gam, Gr(cb,cb), zero, work1, 'dag_2nd')
          call matmul_gpu(hh, one, Gr(cb,cb), work1, zero, Gn(cb,cb))
          call destroyAll(work1)


          !***************************************************************
          !*** Gr(cb-1, cb) already allocated: must insert if in loop ****
          !***************************************************************

          do i=cb-1, 1, -1
             !Gr(i,cb) = - gL(i) ESH(i,i+1) Gr(i+1,cb)
             if (i .eq. cb-1) then
                call createAll(work1, ESH(i,i+1)%nrow, Gr(i+1,cb)%ncol)
             else 
                call createAll(work1, ESH(i,i+1)%nrow, Gr(i+1,cb)%ncol)   
                call createAll(Gr(i,cb), gsml(i)%nrow, work1%ncol)
             end if
             call copyToGPU(ESH(i,i+1))
             call matmul_gpu(hh, one, ESH(i,i+1), Gr(i+1,cb),zero, work1)
             call deleteGPU(ESH(i,i+1))
             call matmul_gpu(hh, mone, gsml(i), work1, zero, Gr(i,cb))
             call destroyAll(work1)

             !Gn(i, i) = Gr(i, cb) Gam(cb) Gr(i, cb)^+
             call createAll(work2, Gam%nrow, Gr(i,cb)%nrow)
             call matmul_gpu(hh, one, Gam, Gr(i,cb), zero, work2, 'dag_2nd')
             call matmul_gpu(hh, frmdiff, Gr(i,cb), work2, zero, Gn(i,i))

             !Gn(i+1,i)  = Gr(i+1, cb) Gam(cb) Gr(i, cb)^+
             call matmul_gpu(hh, frmdiff, Gr(i+1,cb), work2, zero, Gn(i+1,i))
             call destroyAll(work2)

             !Gn(i,i+1)  = Gr(i, cb) Gam(cb) Gr(i+1, cb)^+
             call createAll(work3, Gam%nrow, Gr(i+1,cb)%nrow)
             call matmul_gpu(hh, frmdiff, Gam, Gr(i+1,cb), zero, work3, 'dag_2nd')
             call matmul_gpu(hh, one, Gr(i,cb), work3, zero, Gn(i,i+1))
             call destroyAll(work3)
          end do

          !***************************************************************
          !*** Gr(cb+1, cb) already allocated: must insert if in loop ****
          !***************************************************************

          do i=cb+1, nbl
             !Gr(i,cb) = - gR(i) ESH(i,i-1) Gr(i-1,cb)
             if (i .eq. cb+1) then
                call createAll(work1, ESH(i,i-1)%nrow, Gr(i-1,cb)%ncol)

             else 
                call createAll(work1, ESH(i,i-1)%nrow, Gr(i-1,cb)%ncol)
                call createAll(Gr(i,cb), gsmr(i)%nrow, work1%ncol)
             end if
             call copyToGPU(ESH(i,i-1))
             call matmul_gpu(hh, one, ESH(i,i-1), Gr(i-1,cb),zero, work1)
             call deleteGPU(ESH(i,i-1))
             call matmul_gpu(hh, mone, gsmr(i), work1, zero, Gr(i,cb))
             call destroyAll(work1)


             !Gn(i, i) = Gr(i, cb) Gam(cb) Gr(i, cb)^+
             call createAll(work2, Gam%nrow, Gr(i,cb)%nrow)
             call matmul_gpu(hh, one, Gam, Gr(i,cb), zero, work2, 'dag_2nd')
             call matmul_gpu(hh, frmdiff, Gr(i,cb), work2, zero, Gn(i,i))

             !Gn(i-1,i)  = Gr(i-1, cb) Gam(cb) Gr(i, cb)^+
             call matmul_gpu(hh, frmdiff, Gr(i-1,cb), work2, zero, Gn(i-1,i))
             call destroyAll(work2)

             !Gn(i,i-1)  = Gr(i, cb) Gam(cb) Gr(i-1, cb)^+
             call createAll(work3, Gam%nrow, Gr(i-1,cb)%nrow)
             call matmul_gpu(hh, one, Gam, Gr(i-1,cb), zero, work3, 'dag_2nd')
             call matmul_gpu(hh, frmdiff, Gr(i,cb), work3, zero, Gn(i,i-1))
             call destroyAll(work3)
          end do

          do i=cb-2, 1, -1
             call destroyAll(Gr(i,cb))
          end do

          do i=cb+2, nbl
             call destroyAll(Gr(i,cb))
          end do
          call destroyAll(Gam)
       endif
    end do

  end subroutine calculate_Gn_tridiag_blocks_sp

  subroutine calculate_Gn_tridiag_blocks_dp(negf, ESH, SelfEneR, frm, ref, struct, gsml, gsmr, Gr, Gn)
    !In/Out
    type(Tnegf) :: negf
    type(z_DNS), dimension(:,:), intent(in) :: ESH, Gr
    type(z_DNS), dimension(:,:), intent(inout) :: Gn
    type(z_DNS), dimension(:), intent(in) :: gsmr, gsml
    type(z_DNS), dimension(:), intent(in) :: SelfEneR
    type(Tstruct_info), intent(in) :: struct
    real(dp), dimension(:), intent(in) :: frm
    integer, intent(in) :: ref

    !Work
    type(CublasHandle) :: hh
    Type(z_DNS) :: Gam
    type(z_DNS) :: work1, work2, work3
    integer :: i, j, istat
    integer :: cb, nbl, ncont
    complex(dp) :: frmdiff
    complex(dp), parameter :: one = (1.0_dp, 0.0_dp)
    complex(dp), parameter :: mone = (-1.0_dp, 0.0_dp)
    complex(dp), parameter :: zero = (0.0_dp, 0.0_dp)
    real(dp) :: somma

    ncont = struct%num_conts
    nbl = struct%num_PLs
    hh = negf%hcublas

    do j=1,ncont

       if (j.NE.ref .AND. ABS(frm(j)-frm(ref)).GT.EPS) THEN

          cb=struct%cblk(j) ! block corresponding to contact j    

          !Creating Gam = i (SE -  SE^+)
          call createAll(Gam, SelfEneR(j)%nrow, SelfEneR(j)%ncol)
          call spectral_gpu(hh, SelfEneR(j), Gam)

          frmdiff = cmplx(frm(j)-frm(ref),0.0_dp,dp)

          ! Computation of Gn(cb,cb) = Gr(cb,cb) Gam(cb) Gr(cb,cb)^+
          call createAll(work1, Gam%nrow, Gr(cb,cb)%nrow) !it's Gr^+ --> I take nrow
          call matmul_gpu(hh, frmdiff, Gam, Gr(cb,cb), zero, work1, 'dag_2nd')
          call matmul_gpu(hh, one, Gr(cb,cb), work1, zero, Gn(cb,cb))
          call destroyAll(work1)


          !***************************************************************
          !*** Gr(cb-1, cb) already allocated: must insert if in loop ****
          !***************************************************************

          do i=cb-1, 1, -1
             !Gr(i,cb) = - gL(i) ESH(i,i+1) Gr(i+1,cb)
             if (i .eq. cb-1) then
                call createAll(work1, ESH(i,i+1)%nrow, Gr(i+1,cb)%ncol)
             else 
                call createAll(work1, ESH(i,i+1)%nrow, Gr(i+1,cb)%ncol)   
                call createAll(Gr(i,cb), gsml(i)%nrow, work1%ncol)
             end if
             call copyToGPU(ESH(i,i+1))
             call matmul_gpu(hh, one, ESH(i,i+1), Gr(i+1,cb),zero, work1)
             call deleteGPU(ESH(i,i+1))
             call matmul_gpu(hh, mone, gsml(i), work1, zero, Gr(i,cb))
             call destroyAll(work1)

             !Gn(i, i) = Gr(i, cb) Gam(cb) Gr(i, cb)^+
             call createAll(work2, Gam%nrow, Gr(i,cb)%nrow)
             call matmul_gpu(hh, one, Gam, Gr(i,cb), zero, work2, 'dag_2nd')
             call matmul_gpu(hh, frmdiff, Gr(i,cb), work2, zero, Gn(i,i))

             !Gn(i+1,i)  = Gr(i+1, cb) Gam(cb) Gr(i, cb)^+
             call matmul_gpu(hh, frmdiff, Gr(i+1,cb), work2, zero, Gn(i+1,i))
             call destroyAll(work2)

             !Gn(i,i+1)  = Gr(i, cb) Gam(cb) Gr(i+1, cb)^+
             call createAll(work3, Gam%nrow, Gr(i+1,cb)%nrow)
             call matmul_gpu(hh, frmdiff, Gam, Gr(i+1,cb), zero, work3, 'dag_2nd')
             call matmul_gpu(hh, one, Gr(i,cb), work3, zero, Gn(i,i+1))
             call destroyAll(work3)
          end do

          !***************************************************************
          !*** Gr(cb+1, cb) already allocated: must insert if in loop ****
          !***************************************************************

          do i=cb+1, nbl
             !Gr(i,cb) = - gR(i) ESH(i,i-1) Gr(i-1,cb)
             if (i .eq. cb+1) then
                call createAll(work1, ESH(i,i-1)%nrow, Gr(i-1,cb)%ncol)

             else 
                call createAll(work1, ESH(i,i-1)%nrow, Gr(i-1,cb)%ncol)
                call createAll(Gr(i,cb), gsmr(i)%nrow, work1%ncol)
             end if
             call copyToGPU(ESH(i,i-1))
             call matmul_gpu(hh, one, ESH(i,i-1), Gr(i-1,cb),zero, work1)
             call deleteGPU(ESH(i,i-1))
             call matmul_gpu(hh, mone, gsmr(i), work1, zero, Gr(i,cb))
             call destroyAll(work1)


             !Gn(i, i) = Gr(i, cb) Gam(cb) Gr(i, cb)^+
             call createAll(work2, Gam%nrow, Gr(i,cb)%nrow)
             call matmul_gpu(hh, one, Gam, Gr(i,cb), zero, work2, 'dag_2nd')
             call matmul_gpu(hh, frmdiff, Gr(i,cb), work2, zero, Gn(i,i))

             !Gn(i-1,i)  = Gr(i-1, cb) Gam(cb) Gr(i, cb)^+
             call matmul_gpu(hh, frmdiff, Gr(i-1,cb), work2, zero, Gn(i-1,i))
             call destroyAll(work2)

             !Gn(i,i-1)  = Gr(i, cb) Gam(cb) Gr(i-1, cb)^+
             call createAll(work3, Gam%nrow, Gr(i-1,cb)%nrow)
             call matmul_gpu(hh, one, Gam, Gr(i-1,cb), zero, work3, 'dag_2nd')
             call matmul_gpu(hh, frmdiff, Gr(i,cb), work3, zero, Gn(i,i-1))
             call destroyAll(work3)
          end do

          do i=cb-2, 1, -1
             call destroyAll(Gr(i,cb))
          end do

          do i=cb+2, nbl
             call destroyAll(Gr(i,cb))
          end do
          call destroyAll(Gam)
       endif
    end do

  end subroutine calculate_Gn_tridiag_blocks_dp

  subroutine calculate_single_transmission_2_contacts_sp(negf,ni,nf,ESH,SelfEneR,cblk,tun_proj,Gr,TUN)
    type(Tnegf), intent(in) :: negf 
    integer, intent(in) :: ni,nf
    type(c_DNS), intent(in) :: SelfEneR(MAXNCONT)
    type(c_DNS), intent(in) :: ESH(:,:)
    integer, intent(in) :: cblk(:)
    type(intarray), intent(in) :: tun_proj
    type(c_DNS), dimension(:,:), intent(in) :: Gr
    real(sp), intent(out) :: TUN

    !Work variables
    type(CublasHandle) :: hh
    Integer :: ct1, bl1
    logical, dimension(:), allocatable :: tun_mask
    Type(c_DNS) :: work1, work2, GAM1_dns, GA, TRS, AA
    complex(sp), parameter :: j = (0.0_sp,1.0_sp)  ! CMPX unity
    complex(sp), parameter :: mj = (0.0_sp,-1.0_sp) 
    complex(sp), parameter :: one = (1.0_sp, 0.0_sp)
    complex(sp), parameter :: mone = (-1.0_sp, 0.0_sp)
    complex(sp), parameter :: zero = (0.0_sp, 0.0_sp)

    if (size(cblk).gt.2) then
       write(*,*) "ERROR: calculate_single_transmission_2_contacts is valid only for 2 contacts"
       TUN = 0.0_sp
       return
    endif

    hh = negf%hcublas

    !Arrange contacts in way that order between first and second is always the
    !same (always ct1 < ct2)

    if (cblk(ni).lt.cblk(nf)) then
       ct1=ni;
    else
       ct1=nf;
    endif

    bl1=cblk(ct1);

    ! Computes the Gamma matrices
    call createAll(GAM1_dns, SelfEneR(ct1)%nrow, SelfEneR(ct1)%ncol)
    call spectral_gpu(hh, SelfEneR(ct1), GAM1_dns)

    ! Work to compute transmission matrix (Gamma G Gamma G)
    call createAll(work1, GAM1_dns%nrow, Gr(bl1,bl1)%ncol)
    call createAll(work2, work1%nrow, GAM1_dns%ncol)
    call matmul_gpu(hh, one, GAM1_dns, Gr(bl1,bl1), zero, work1)
    call matmul_gpu(hh, one, work1, GAM1_dns, zero, work2)

    call destroyAll(work1)

    call createAll(work1, work2%nrow, Gr(bl1,bl1)%nrow)
    call matmul_gpu(hh, one, work2, Gr(bl1,bl1), zero, work1,'dag_2nd')
    call destroyAll(work2)

    call createAll(AA, Gr(bl1,bl1)%ncol, Gr(bl1,bl1)%nrow)
    call matsum_gpu(hh, j, Gr(bl1,bl1), mj, Gr(bl1,bl1), AA, 'dag_2nd')

    call createAll(work2, GAM1_dns%nrow, AA%ncol)
    call matmul_gpu(hh, one, GAM1_dns, AA, zero, work2)
    call destroyAll(GAM1_dns)
    call destroyAll(AA)

    call createAll(TRS, work1%nrow, work1%ncol)
    call matsum_gpu(hh, one, work2, mone, work1, TRS)
    call get_tun_mask(ESH, bl1, tun_proj, tun_mask)
    call trace_gpu(hh, TRS, TUN, tun_mask)
    call log_deallocate(tun_mask)

    call destroyAll(TRS)
    call destroyAll(work1)
    call destroyAll(work2)
  
  end subroutine calculate_single_transmission_2_contacts_sp

  subroutine calculate_single_transmission_2_contacts_dp(negf,ni,nf,ESH,SelfEneR,cblk,tun_proj,Gr,TUN)
    type(Tnegf), intent(in) :: negf 
    integer, intent(in) :: ni,nf
    type(z_DNS), intent(in) :: SelfEneR(MAXNCONT)
    type(z_DNS), intent(in) :: ESH(:,:)
    integer, intent(in) :: cblk(:)
    type(intarray), intent(in) :: tun_proj
    type(z_DNS), dimension(:,:), intent(in) :: Gr
    real(dp), intent(out) :: TUN

    !Work variables
    type(CublasHandle) :: hh
    Integer :: ct1, bl1
    logical, dimension(:), allocatable :: tun_mask
    Type(z_DNS) :: work1, work2, GAM1_dns, GA, TRS, AA
    complex(dp), parameter :: j = (0.0_dp,1.0_dp)  ! CMPX unity
    complex(dp), parameter :: mj = (0.0_dp,-1.0_dp) 
    complex(dp), parameter :: one = (1.0_dp, 0.0_dp)
    complex(dp), parameter :: mone = (-1.0_dp, 0.0_dp)
    complex(dp), parameter :: zero = (0.0_dp, 0.0_dp)
    real(dp) :: summ

    if (size(cblk).gt.2) then
       write(*,*) "ERROR: calculate_single_transmission_2_contacts is valid only for 2 contacts"
       TUN = 0.0_dp
       return
    endif

    hh = negf%hcublas

    !Arrange contacts in way that order between first and second is always the
    !same (always ct1 < ct2)

    if (cblk(ni).lt.cblk(nf)) then
       ct1=ni;
    else
       ct1=nf;
    endif

    bl1=cblk(ct1);

    ! Computes the Gamma matrices
    call createAll(GAM1_dns, SelfEneR(ct1)%nrow, SelfEneR(ct1)%ncol)
    call spectral_gpu(hh, SelfEneR(ct1), GAM1_dns)

    ! Work to compute transmission matrix (Gamma G Gamma G)
    call createAll(work1, GAM1_dns%nrow, Gr(bl1,bl1)%ncol)
    call createAll(work2, work1%nrow, GAM1_dns%ncol)
    call matmul_gpu(hh, one, GAM1_dns, Gr(bl1,bl1), zero, work1)
    call matmul_gpu(hh, one, work1, GAM1_dns, zero, work2)

    call destroyAll(work1)

    call createAll(work1, work2%nrow, Gr(bl1,bl1)%nrow)
    call matmul_gpu(hh, one, work2, Gr(bl1,bl1), zero, work1,'dag_2nd')
    call destroyAll(work2)

    call createAll(AA, Gr(bl1,bl1)%ncol, Gr(bl1,bl1)%nrow)
    call matsum_gpu(hh, j, Gr(bl1,bl1), mj, Gr(bl1,bl1), AA, 'dag_2nd')

    call createAll(work2, GAM1_dns%nrow, AA%ncol)
    call matmul_gpu(hh, one, GAM1_dns, AA, zero, work2)
    call destroyAll(GAM1_dns)
    call destroyAll(AA)

    call createAll(TRS, work1%nrow, work1%ncol)
    call matsum_gpu(hh, one, work2, mone, work1, TRS)
    call get_tun_mask(ESH, bl1, tun_proj, tun_mask)
    call trace_gpu(hh, TRS, TUN, tun_mask)
    call log_deallocate(tun_mask)

    call destroyAll(TRS)
    call destroyAll(work1)
    call destroyAll(work2)
  
  end subroutine calculate_single_transmission_2_contacts_dp

  subroutine calculate_single_transmission_N_contacts_sp(negf,ni,nf,ESH,SelfEneR,cblk,tun_proj,gsmr,Gr,TUN)
    type(Tnegf), intent(in) :: negf
    integer, intent(in) :: ni,nf
    type(c_DNS), intent(in) :: SelfEneR(MAXNCONT)
    type(c_DNS), intent(in) :: ESH(:,:)
    type(c_DNS), dimension(:),intent(in) :: gsmr
    type(c_DNS), dimension(:,:),intent(in) :: Gr
    integer, intent(in) :: cblk(:)
    type(intarray), intent(in) :: tun_proj
    real(sp), intent(out) :: TUN

    !Work variables
    type(CublasHandle) :: hh
    Integer :: ct1, ct2, bl1, bl2, i, nbl
    logical, dimension(:), allocatable :: tun_mask
    Type(c_DNS) :: work1, work2, GAM1_dns, GAM2_dns, TRS
    Real(sp) :: max
    complex(sp), parameter :: one = (1.0_sp, 0.0_sp)
    complex(sp), parameter :: mone = (-1.0_sp, 0.0_sp)
    complex(sp), parameter :: zero = (0.0_sp, 0.0_sp)

    !Arrange contacts in way that order between first and second is always the
    !same (always ct1 < ct2)

    if (cblk(ni).lt.cblk(nf)) then
       ct1=ni;ct2=nf;
    else
       ct1=nf;ct2=ni;
    endif

    bl1=cblk(ct1); bl2=cblk(ct2);
    nbl = size(cblk)
    hh = negf%hcublas
    ! in this way nt1 < nt2 by construction
    if ( nbl.gt.1 .and. (bl2-bl1).gt.1) then

       ! Compute column-blocks of Gr(i,bl1) up to i=bl2
       ! Gr(i,bl1) = -gr(i) T(i,i-1) Gr(i-1,bl1)
       do i = bl1+1, bl2
          !Checks whether previous block is non null.
          !If so next block is also null => TUN = 0
          call copyFromGPU(Gr(i-1,bl1))
          max=maxval(abs(Gr(i-1,bl1)%val))
          
          if (max.lt.EPS) then
             TUN = EPS*EPS !for log plots
             !Destroy also the block adjecent to diagonal since
             !this is not deallocated anymore in calling subroutine
             if (i.gt.(bl1+1)) call destroyAll(Gr(i-1,bl1))
             return
          endif

          !Checks whether block has been created, if not do it
          if (.not.allocated(Gr(i,bl1)%val)) then  

             call createAll(work1, gsmr(i)%nrow, ESH(i,i-1)%ncol)
             call createAll(Gr(i,bl1), work1%nrow, Gr(i-1,bl1)%ncol)

             call copyToGPU(ESH(i,i-1))
             call matmul_gpu(hh, mone, gsmr(i), ESH(i,i-1), zero, work1)
             call deleteGPU(ESH(i,i-1))

             call matmul_gpu(hh, one, work1, Gr(i-1,bl1) ,zero, Gr(i,bl1))
             call destroyAll(work1)

          endif

          ! avoid destroying blocks closer to diagonal
          if (i.gt.(bl1+2)) call destroyAll(Gr(i-1,bl1))
       end do

    endif
    ! Computes the Gamma matrices
    call createAll(GAM1_dns, SelfEneR(ct1)%nrow, SelfEneR(ct1)%ncol)
    call createAll(GAM2_dns, SelfEneR(ct2)%nrow, SelfEneR(ct2)%ncol)
    call spectral_gpu(hh, SelfEneR(ct1),GAM1_dns)
    call spectral_gpu(hh ,SelfEneR(ct2),GAM2_dns)

    ! Work to compute transmission matrix (Gamma2 Gr Gamma1 Ga)
    call createAll(work1, Gr(bl2,bl1)%nrow, GAM1_dns%ncol)
    call matmul_gpu(hh, one, Gr(bl2,bl1), GAM1_dns, zero, work1)

    call createAll(work2, GAM2_dns%nrow, work1%ncol)
    call matmul_gpu(hh, one, GAM2_dns, work1, zero, work2)

    call destroyAll(work1)
    call destroyAll(GAM2_dns)
    call destroyAll(GAM1_dns)

    call createAll(TRS, work2%nrow, Gr(bl2,bl1)%nrow)

    call matmul_gpu(hh, one, work2, Gr(bl2,bl1), zero, TRS, 'dag_2nd')
    call destroyAll(work2)
    if (bl2.gt.bl1+1) call destroyAll(Gr(bl2,bl1))

    call get_tun_mask(ESH, bl2, tun_proj, tun_mask)
    call trace_gpu(hh, TRS, TUN, tun_mask) 
    call log_deallocate(tun_mask)

    call destroyAll(TRS)

  end subroutine calculate_single_transmission_N_contacts_sp

  subroutine calculate_single_transmission_N_contacts_dp(negf,ni,nf,ESH,SelfEneR,cblk,tun_proj,gsmr,Gr,TUN)
    type(Tnegf), intent(in) :: negf
    integer, intent(in) :: ni,nf
    type(z_DNS), intent(in) :: SelfEneR(MAXNCONT)
    type(z_DNS), intent(in) :: ESH(:,:)
    type(z_DNS), dimension(:),intent(in) :: gsmr
    type(z_DNS), dimension(:,:),intent(in) :: Gr
    integer, intent(in) :: cblk(:)
    type(intarray), intent(in) :: tun_proj
    real(dp), intent(out) :: TUN

    !Work variables
    type(CublasHandle) :: hh
    Integer :: ct1, ct2, bl1, bl2, i, nbl
    logical, dimension(:), allocatable :: tun_mask
    Type(z_DNS) :: work1, work2, GAM1_dns, GAM2_dns, TRS
    Real(dp) :: max
    complex(dp), parameter :: one = (1.0_dp, 0.0_dp)
    complex(dp), parameter :: mone = (-1.0_dp, 0.0_dp)
    complex(dp), parameter :: zero = (0.0_dp, 0.0_dp)

    !Arrange contacts in way that order between first and second is always the
    !same (always ct1 < ct2)

    if (cblk(ni).lt.cblk(nf)) then
       ct1=ni;ct2=nf;
    else
       ct1=nf;ct2=ni;
    endif

    bl1=cblk(ct1); bl2=cblk(ct2);
    nbl = size(cblk)
    hh = negf%hcublas
    ! in this way nt1 < nt2 by construction
    if ( nbl.gt.1 .and. (bl2-bl1).gt.1) then

       ! Compute column-blocks of Gr(i,bl1) up to i=bl2
       ! Gr(i,bl1) = -gr(i) T(i,i-1) Gr(i-1,bl1)
       do i = bl1+1, bl2
          !Checks whether previous block is non null.
          !If so next block is also null => TUN = 0
          call copyFromGPU(Gr(i-1,bl1))
          max=maxval(abs(Gr(i-1,bl1)%val))
          
          if (max.lt.EPS) then
             TUN = EPS*EPS !for log plots
             !Destroy also the block adjecent to diagonal since
             !this is not deallocated anymore in calling subroutine
             if (i.gt.(bl1+1)) call destroyAll(Gr(i-1,bl1))
             return
          endif

          !Checks whether block has been created, if not do it
          if (.not.allocated(Gr(i,bl1)%val)) then  

             call createAll(work1, gsmr(i)%nrow, ESH(i,i-1)%ncol)
             call createAll(Gr(i,bl1), work1%nrow, Gr(i-1,bl1)%ncol)

             call copyToGPU(ESH(i,i-1))
             call matmul_gpu(hh, mone, gsmr(i), ESH(i,i-1), zero, work1)
             call deleteGPU(ESH(i,i-1))

             call matmul_gpu(hh, one, work1, Gr(i-1,bl1) ,zero, Gr(i,bl1))
             call destroyAll(work1)

          endif

          ! avoid destroying blocks closer to diagonal
          if (i.gt.(bl1+2)) call destroyAll(Gr(i-1,bl1))
       end do

    endif
    ! Computes the Gamma matrices
    call createAll(GAM1_dns, SelfEneR(ct1)%nrow, SelfEneR(ct1)%ncol)
    call createAll(GAM2_dns, SelfEneR(ct2)%nrow, SelfEneR(ct2)%ncol)
    call spectral_gpu(hh, SelfEneR(ct1),GAM1_dns)
    call spectral_gpu(hh ,SelfEneR(ct2),GAM2_dns)

    ! Work to compute transmission matrix (Gamma2 Gr Gamma1 Ga)
    call createAll(work1, Gr(bl2,bl1)%nrow, GAM1_dns%ncol)
    call matmul_gpu(hh, one, Gr(bl2,bl1), GAM1_dns, zero, work1)

    call createAll(work2, GAM2_dns%nrow, work1%ncol)
    call matmul_gpu(hh, one, GAM2_dns, work1, zero, work2)

    call destroyAll(work1)
    call destroyAll(GAM2_dns)
    call destroyAll(GAM1_dns)

    call createAll(TRS, work2%nrow, Gr(bl2,bl1)%nrow)

    call matmul_gpu(hh, one, work2, Gr(bl2,bl1), zero, TRS, 'dag_2nd')
    call destroyAll(work2)
    if (bl2.gt.bl1+1) call destroyAll(Gr(bl2,bl1))

    call get_tun_mask(ESH, bl2, tun_proj, tun_mask)
    call trace_gpu(hh, TRS, TUN, tun_mask) 
    call log_deallocate(tun_mask)

    call destroyAll(TRS)

  end subroutine calculate_single_transmission_N_contacts_dp

  subroutine build_ESH_onGPU_dp(negf, E, S, H, ESH)
    type(z_DNS), intent(inout), dimension(:,:) :: ESH
    type(z_DNS), intent(inout), dimension(:,:) :: S, H
    type(Tnegf), intent(in) :: negf
    complex(dp), intent(in) :: E
         
    type(z_CSR) :: H_csr, S_csr
    integer :: i, nbl
    type(cublasHandle) :: hh
    complex(dp), parameter :: one = (1.0_dp,0.0_dp )
    complex(dp), parameter :: mone = (-1.0_dp,0.0_dp )

    nbl = negf%str%num_PLs
    H_csr = negf%H
    S_csr = negf%S

    hh = negf%hcublas

    call csr2blk_sod(H_csr, H, negf%str%mat_PL_start)
    call csr2blk_sod(S_csr, S, negf%str%mat_PL_start)
    call copy_trid_toGPU(S)
    call copy_trid_toGPU(H)

    call createAll(ESH(1,1), S(1,1)%nrow, S(1,1)%ncol)
    call matsum_gpu(hh, E, S(1,1), mone, H(1,1), ESH(1,1))

    do i=2,nbl
       call createAll(ESH(i,i), S(i,i)%nrow, S(i,i)%ncol)
       call matsum_gpu(hh,E,S(i,i), mone, H(i,i), ESH(i,i))

       call createAll(ESH(i-1,i), S(i-1,i)%nrow, S(i-1,i)%ncol)
       call matsum_gpu(hh, E,S(i-1,i), mone, H(i-1,i), ESH(i-1,i))

       call createAll(ESH(i,i-1), S(i,i-1)%nrow, S(i,i-1)%ncol)
       call matsum_gpu(hh, E,S(i,i-1), mone, H(i,i-1), ESH(i,i-1))
    end do

  end subroutine build_ESH_onGPU_dp        

!  subroutine build_ESH_onGPU_sp(negf, E, S, H, ESH)
!    type(c_DNS), intent(inout), dimension(:,:) :: ESH
!    type(c_DNS), intent(inout), dimension(:,:) :: S, H
!    type(Tnegf), intent(in) :: negf
!    complex(sp), intent(in) :: E
!         
!    type(c_CSR) :: H_csr, S_csr
!    integer :: i, nbl
!    type(cublasHandle) :: hh
!    complex(sp), parameter :: one = (1.0_sp,0.0_sp )
!    complex(sp), parameter :: mone = (-1.0_sp,0.0_sp )
!
!    nbl = negf%str%num_PLs
!    H_csr = negf%H
!    S_csr = negf%S
!
!    hh = negf%hcublas
!
!    call csr2blk_sod(H_csr, H, negf%str%mat_PL_start)
!    call csr2blk_sod(S_csr, S, negf%str%mat_PL_start)
!    call copy_trid_toGPU(S)
!    call copy_trid_toGPU(H)
!
!    call createAll(ESH(1,1), S(1,1)%nrow, S(1,1)%ncol)
!    call matsum_gpu(hh, E, S(1,1), mone, H(1,1), ESH(1,1))
!
!    do i=2,nbl
!       call createAll(ESH(i,i), S(i,i)%nrow, S(i,i)%ncol)
!       call matsum_gpu(hh,E,S(i,i), mone, H(i,i), ESH(i,i))
!
!       call createAll(ESH(i-1,i), S(i-1,i)%nrow, S(i-1,i)%ncol)
!       call matsum_gpu(hh, E,S(i-1,i), mone, H(i-1,i), ESH(i-1,i))
!
!       call createAll(ESH(i,i-1), S(i,i-1)%nrow, S(i,i-1)%ncol)
!       call matsum_gpu(hh, E,S(i,i-1), mone, H(i,i-1), ESH(i,i-1))
!    end do
!
!  end subroutine build_ESH_onGPU_sp        

  subroutine get_tun_mask_sp(ESH,nbl,tun_proj,tun_mask)
    Type(c_DNS), intent(in) :: ESH(:,:)
    integer, intent(in) :: nbl
    type(intarray), intent(in) :: tun_proj
    logical, intent(out), allocatable :: tun_mask(:)    

    integer :: ii, istart, iend, ind

    call log_allocate(tun_mask, ESH(nbl,nbl)%nrow)

    if (allocated(tun_proj%indexes)) then
       tun_mask = .false.

       ! set the start/end indices of nbl
       ! NB: istart has offset -1 to avoid +/-1 operations
       istart = 0
       do ii = 1, nbl-1
          istart = istart + ESH(ii,ii)%nrow
       end do
       iend = istart + ESH(nbl,nbl)%nrow + 1  

       ! select the indices in tun_proj 
       do ii = 1, size(tun_proj%indexes)
          ind = tun_proj%indexes(ii)
          if (ind > istart .and. ind < iend) then
             tun_mask(ind - istart) = .true. 
          end if
       end do
    else
       tun_mask = .true.
    end if

  end subroutine get_tun_mask_sp

  subroutine get_tun_mask_dp(ESH,nbl,tun_proj,tun_mask)
    Type(z_DNS), intent(in) :: ESH(:,:)
    integer, intent(in) :: nbl
    type(intarray), intent(in) :: tun_proj
    logical, intent(out), allocatable :: tun_mask(:)    

    integer :: ii, istart, iend, ind

    call log_allocate(tun_mask, ESH(nbl,nbl)%nrow)

    if (allocated(tun_proj%indexes)) then
       tun_mask = .false.

       ! set the start/end indices of nbl
       ! NB: istart has offset -1 to avoid +/-1 operations
       istart = 0
       do ii = 1, nbl-1
          istart = istart + ESH(ii,ii)%nrow
       end do
       iend = istart + ESH(nbl,nbl)%nrow + 1  

       ! select the indices in tun_proj
       do ii = 1, size(tun_proj%indexes)
          ind = tun_proj%indexes(ii)
          if (ind > istart .and. ind < iend) then
             tun_mask(ind - istart) = .true. 
          end if
       end do
    else
       tun_mask = .true.
    end if

  end subroutine get_tun_mask_dp

end module iterative_gpu
