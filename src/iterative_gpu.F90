module iterative_gpu
  use ln_precision
  use ln_constants, only : pi
  use ln_allocation
  use mat_def
  use ln_structure, only : TStruct_Info
  use lib_param, only : MAXNCONT, Tnegf, intarray
  use mpi_globals, only : id, numprocs, id0
  use sparsekit_drv, only: csr2blk_sod
  use clock
  use cudautils
  use cublas_v2
  use cusolverDn
  use, intrinsic :: ieee_arithmetic

  implicit none
  private 

  public :: calculate_gsmr_blocks
  public :: calculate_Gr_tridiag_blocks
  public :: calculate_Gn_tridiag_blocks
  public :: calculate_single_transmission_2_contacts
  public :: calculate_single_transmission_N_contacts
  public :: check_convergence_trid
  public :: check_convergence_vec

  interface calculate_gsmr_blocks
     module procedure  calculate_gsmr_blocks_sp
     module procedure  calculate_gsmr_blocks_dp
  end interface calculate_gsmr_blocks

  interface calculate_Gr_tridiag_blocks
     module procedure calculate_Gr_tridiag_blocks_sp
     module procedure calculate_Gr_tridiag_blocks_dp
  end interface calculate_Gr_tridiag_blocks

!  interface calculate_Gn_tridiag_blocks
!     module procedure calculate_Gn_tridiag_blocks_sp
!     module procedure calculate_Gn_tridiag_blocks_dp
!  end interface calculate_Gn_tridiag_blocks

  interface calculate_single_transmission_2_contacts
     module procedure calculate_single_transmission_2_contacts_sp
     module procedure calculate_single_transmission_2_contacts_dp
  end interface calculate_single_transmission_2_contacts

  interface calculate_single_transmission_N_contacts
     module procedure calculate_single_transmission_N_contacts_sp
     module procedure calculate_single_transmission_N_contacts_dp
  end interface calculate_single_transmission_N_contacts

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
    call inverse_gpu(hh, hhsol, ESH(sbl,sbl)%val, gsmr(sbl)%val, istat)
    call deleteGPU(ESH(sbl,sbl))

    do i=sbl-1,ebl,-1

       call createAll(work1, ESH(i,i+1)%nrow, gsmr(i+1)%ncol)    
       call copyToGPU(ESH(i,i+1))
       call matmul_gpu(hh, one, ESH(i,i+1)%val, gsmr(i+1)%val, zero, work1%val)
       call deleteGPU(ESH(i,i+1))

       if (.not.keep) then
          call destroyAll(gsmr(i+1))
       end if

       call createAll(work2, work1%nrow, ESH(i+1,i)%ncol)
       call copyToGPU(ESH(i+1,i))
       call matmul_gpu(hh, one, work1%val, ESH(i+1,i)%val, zero, work2%val)
       call deleteGPU(ESH(i+1,i))

       call destroyAll(work1)

       call createAll(work1,  ESH(i,i)%nrow, ESH(i,i)%ncol)
       call copyToGPU(ESH(i,i))
       call add_cublas(hh, work1%val, ESH(i,i)%val, mone, work2%val)
       call deleteGPU(ESH(i,i))

       call destroyAll(work2)

       call createAll(gsmr(i), work1%nrow, work1%ncol)
       call inverse_gpu(hh, hhsol, work1%val, gsmr(i)%val, istat)    
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
    call inverse_gpu(hh, hhsol, ESH(sbl,sbl)%val, gsmr(sbl)%val, istat)
    call deleteGPU(ESH(sbl,sbl))

    do i=sbl-1,ebl,-1

       call createAll(work1, ESH(i,i+1)%nrow, gsmr(i+1)%ncol)    
       call copyToGPU(ESH(i,i+1))
       call matmul_gpu(hh, one, ESH(i,i+1)%val, gsmr(i+1)%val, zero, work1%val)
       call deleteGPU(ESH(i,i+1))

       if (.not.keep) then
          call destroyAll(gsmr(i+1))
       end if

       call createAll(work2, work1%nrow, ESH(i+1,i)%ncol)
       call copyToGPU(ESH(i+1,i))
       call matmul_gpu(hh, one, work1%val, ESH(i+1,i)%val, zero, work2%val)
       call deleteGPU(ESH(i+1,i))

       call destroyAll(work1)

       call createAll(work1,  ESH(i,i)%nrow, ESH(i,i)%ncol)
       call copyToGPU(ESH(i,i))
       call add_cublas(hh, work1%val, ESH(i,i)%val, mone, work2%val)
       call deleteGPU(ESH(i,i))

       call destroyAll(work2)

       call createAll(gsmr(i), work1%nrow, work1%ncol)
       call inverse_gpu(hh, hhsol, work1%val, gsmr(i)%val, istat)    
       call destroyAll(work1)

    end do

  end subroutine calculate_gsmr_blocks_dp

  subroutine calculate_Gr_tridiag_blocks_sp(negf,ESH,gsmr,Gr,sbl,ebl)
    !In/Out
    type(c_DNS), dimension(:,:), intent(inout) :: Gr
    type(c_DNS), dimension(:), intent(in) :: gsmr
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
          call inverse_gpu(hh, hhsol, ESH(sbl,sbl)%val, Gr(sbl,sbl)%val, istat)
          call deleteGPU(ESH(sbl,sbl))
       else
          call createAll(work1, ESH(sbl,sbl)%nrow, ESH(sbl,sbl)%ncol)
          call copyToGPU(ESH(sbl,sbl))
          call copy_mat_gpu(hh, ESH(sbl,sbl)%val, work1%val)
          call deleteGPU(ESH(sbl,sbl))

          if (sbl+1.le.nbl) then
             call createAll(work2, ESH(sbl,sbl+1)%nrow, gsmr(sbl+1)%ncol)
             call copyToGPU(ESH(sbl,sbl+1))
             call matmul_gpu(hh, one, ESH(sbl,sbl+1)%val, gsmr(sbl+1)%val, zero, work2%val)
             call deleteGPU(ESH(sbl,sbl+1))

             call createAll(work3, work2%nrow, ESH(sbl+1,sbl)%ncol)
             call copyToGPU(ESH(sbl+1,sbl))
             call matmul_gpu(hh, one, work2%val, ESH(sbl+1,sbl)%val, zero, work3%val)
             call deleteGPU(ESH(sbl+1,sbl))

             call copyToGPU(ESH(sbl,sbl))
             call add_cublas(hh, work1%val, ESH(sbl,sbl)%val, mone, work3%val)
             call deleteGPU(ESH(sbl,sbl))

             call destroyAll(work2)
             call destroyAll(work3)
          end if
          if (sbl-1.ge.1) then
             stop "Error: Gr_tridiag requires gsml"      
          end if

          call createAll(Gr(sbl,sbl), work1%nrow, work1%ncol)
          call inverse_gpu(hh, hhsol, work1%val, Gr(sbl,sbl)%val, istat)
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
          call matmul_gpu(hh, one, gsmr(i)%val, ESH(i,i-1)%val, zero, work1%val)
          call deleteGPU(ESH(i,i-1))
          call matmul_gpu(hh, mone, work1%val, Gr(i-1,i-1)%val, zero, Gr(i,i-1)%val)

          call destroyAll(work1)

          call createAll(work2, ESH(i-1,i)%nrow, gsmr(i)%ncol)
          call createAll(Gr(i-1,i), Gr(i-1,i-1)%nrow, work2%ncol)
          call copyToGPU(ESH(i-1,i))
          call matmul_gpu(hh, one, ESH(i-1,i)%val, gsmr(i)%val, zero, work2%val)
          call deleteGPU(ESH(i-1,i))
          call matmul_gpu(hh, mone, Gr(i-1,i-1)%val, work2%val, zero, Gr(i-1,i)%val)

          call createAll(work1, Gr(i,i-1)%nrow, work2%ncol)
          call matmul_gpu(hh, mone, Gr(i,i-1)%val, work2%val, zero, work1%val)

          call destroyAll(work2)
          call createAll(Gr(i,i), gsmr(i)%nrow, gsmr(i)%ncol)
          call add_cublas(hh, Gr(i,i)%val, gsmr(i)%val, one, work1%val)
          call destroyAll(work1)
       end do
    else
       stop "Error: Gr_tridiag requires gsml"        
    endif
           
  end subroutine calculate_Gr_tridiag_blocks_sp

  subroutine calculate_Gr_tridiag_blocks_dp(negf,ESH,gsmr,Gr,sbl,ebl)
    !In/Out
    type(z_DNS), dimension(:,:), intent(inout) :: Gr
    type(z_DNS), dimension(:), intent(in) :: gsmr
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
          call inverse_gpu(hh, hhsol, ESH(sbl,sbl)%val, Gr(sbl,sbl)%val, istat)
          call deleteGPU(ESH(sbl,sbl))
       else
          call createAll(work1, ESH(sbl,sbl)%nrow, ESH(sbl,sbl)%ncol)
          call copyToGPU(ESH(sbl,sbl))
          call copy_mat_gpu(hh, ESH(sbl,sbl)%val, work1%val)
          call deleteGPU(ESH(sbl,sbl))

          if (sbl+1.le.nbl) then
             call createAll(work2, ESH(sbl,sbl+1)%nrow, gsmr(sbl+1)%ncol)
             call copyToGPU(ESH(sbl,sbl+1))
             call matmul_gpu(hh, one, ESH(sbl,sbl+1)%val, gsmr(sbl+1)%val, zero, work2%val)
             call deleteGPU(ESH(sbl,sbl+1))

             call createAll(work3, work2%nrow, ESH(sbl+1,sbl)%ncol)
             call copyToGPU(ESH(sbl+1,sbl))
             call matmul_gpu(hh, one, work2%val, ESH(sbl+1,sbl)%val, zero, work3%val)
             call deleteGPU(ESH(sbl+1,sbl))

             call copyToGPU(ESH(sbl,sbl))
             call add_cublas(hh, work1%val, ESH(sbl,sbl)%val, mone, work3%val)
             call deleteGPU(ESH(sbl,sbl))

             call destroyAll(work2)
             call destroyAll(work3)
          end if
          if (sbl-1.ge.1) then
             stop "Error: Gr_tridiag requires gsml"      
          end if

          call createAll(Gr(sbl,sbl), work1%nrow, work1%ncol)
          call inverse_gpu(hh, hhsol, work1%val, Gr(sbl,sbl)%val, istat)
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
          call matmul_gpu(hh, one, gsmr(i)%val, ESH(i,i-1)%val, zero, work1%val)
          call deleteGPU(ESH(i,i-1))
          call matmul_gpu(hh, mone, work1%val, Gr(i-1,i-1)%val, zero, Gr(i,i-1)%val)

          call destroyAll(work1)

          call createAll(work2, ESH(i-1,i)%nrow, gsmr(i)%ncol)
          call createAll(Gr(i-1,i), Gr(i-1,i-1)%nrow, work2%ncol)
          call copyToGPU(ESH(i-1,i))
          call matmul_gpu(hh, one, ESH(i-1,i)%val, gsmr(i)%val, zero, work2%val)
          call deleteGPU(ESH(i-1,i))
          call matmul_gpu(hh, mone, Gr(i-1,i-1)%val, work2%val, zero, Gr(i-1,i)%val)

          call createAll(work1, Gr(i,i-1)%nrow, work2%ncol)
          call matmul_gpu(hh, mone, Gr(i,i-1)%val, work2%val, zero, work1%val)

          call destroyAll(work2)
          call createAll(Gr(i,i), gsmr(i)%nrow, gsmr(i)%ncol)
          call add_cublas(hh, Gr(i,i)%val, gsmr(i)%val, one, work1%val)
          call destroyAll(work1)
       end do
    else
       stop "Error: Gr_tridiag requires gsml"      
    endif
           
  end subroutine calculate_Gr_tridiag_blocks_dp

  subroutine check_convergence_trid(negf,T,nome,gpu)
    type(Tnegf), intent(in) :: negf      
    type(z_DNS), dimension(:,:), intent(in) :: T
    character(3), intent(in) :: nome
    logical, intent(in) :: gpu

    integer :: nbl, i
    type(CublasHandle) :: hh   
    real(dp) :: summ 
    character(1) :: ci, cim1

    nbl = size(T,1)
    hh = negf%hcublas
    if (gpu .eq. .true.) then
       !write(*,*) '~-~-~-~-',nome,' check convergence GPU: ~-~-~-~-'     
       call checksum(hh, T(1,1)%val, nome//'(1,1)')

       do i= 2,nbl-1
          write(ci,'(i1)') i
          write(cim1,'(i1)') i-1
          call checksum(hh, T(i,i)%val, nome//'('//ci//','//ci//')')
          call checksum(hh, T(i-1,i)%val, nome//'('//cim1//','//ci//')')
          call checksum(hh, T(i,i-1)%val, nome//'('//ci//','//cim1//')')
       end do
       !write(*,*) '~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-'     
    else
       !write(*,*) '~-~-~-~-',nome,' check convergence CPU: ~-~-~-~-'    
       summ = sum(ABS(T(1,1)%val)) 
       if (ieee_is_nan(summ)) write(*,*) 'CPU:   ',nome,'(',1,1,')=', summ

       do i= 2,nbl-1
          summ = sum(ABS(T(i,i)%val)) 
          if (ieee_is_nan(summ)) write(*,*) 'CPU:   ',nome,'(',i,i,')=', summ
          summ = sum(ABS(T(i-1,i)%val)) 
          if (ieee_is_nan(summ)) write(*,*) 'CPU:   ',nome,'(',i-1,i,')=', summ
          summ = sum(ABS(T(i,i-1)%val)) 
          if (ieee_is_nan(summ)) write(*,*) 'CPU:   ',nome,'(',i,i-1,')=', summ
       end do
       !write(*,*) '~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-'     
    endif     
  end subroutine check_convergence_trid

  subroutine check_convergence_vec(negf,T,nome,gpu)
    type(Tnegf), intent(in) :: negf      
    type(z_DNS), dimension(:), intent(in) :: T
    character(4), intent(in) :: nome
    logical, intent(in) :: gpu

    integer :: nbl, i
    type(CublasHandle) :: hh   
    real(dp) :: summ 
    character(1) :: ci

    nbl = size(T)
    hh = negf%hcublas
    if (gpu .eq. .true.) then
       !write(*,*) '~-~-~-~-',nome,' check convergence GPU: ~-~-~-~-'     

       do i= 2,nbl
          write(ci,'(i1)') i
          call checksum(hh, T(i)%val, nome//'('//ci//')')
       end do
       !write(*,*) '~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-'     
    else 
       !write(*,*) '~-~-~-~-',nome,' check convergence CPU: ~-~-~-~-'     

       do i= 2,nbl
          summ = sum(ABS(T(i)%val))
          if (ieee_is_nan(summ)) write(*,*) 'CPU:   ',nome,'(',i,')=', summ
       end do
       !write(*,*) '~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-'     
    endif
            
  end subroutine check_convergence_vec


  subroutine calculate_Gn_tridiag_blocks(negf,ESH,SelfEneR,frm,ref,struct,gsmr,Gr,Gn)
    type(TNegf), intent(in) :: negf    
    type(z_DNS), dimension(:,:), intent(in) :: ESH, Gr
    type(z_DNS), dimension(:), intent(in) :: SelfEneR, gsmr
    real(dp), dimension(:), intent(in) :: frm
    integer, intent(in) :: ref
    type(Tstruct_info), intent(in) :: struct
    type(z_DNS), dimension(:,:), intent(inout) :: Gn

    !Work
    type(CublasHandle) :: hh
    complex(dp), parameter :: one = (1.0_dp, 0.0_dp)
    complex(dp), parameter :: zero = (0.0_dp, 0.0_dp)
    complex(dp), parameter :: minusone = (-1.0_dp, 0.0_dp)
    type(z_DNS), dimension(:,:), allocatable :: Sigma_n
    type(z_DNS) :: work1, Gam
    complex(dp) :: frmdiff
    integer :: i, j
    integer :: nbl, ncont, cb

    ncont = struct%num_conts
    nbl = struct%num_PLs

    hh = negf%hcublas

    !build Sigma_n from SelfEneR
    call allocate_blk_dns(Sigma_n, nbl)
    call init_tridiag_blk(Sigma_n, ESH)

    if (allocated(negf%inter)) then
      call negf%inter%get_sigma_n(Sigma_n, negf%ie)
    end if

   call copy_trid_toGPU(Sigma_n)

    do j=1,ncont
      frmdiff = cmplx(frm(j) - frm(ref),0.0_dp,dp)
      if (j.NE.ref .AND. ABS(frmdiff).GT.EPS) THEN
        cb=struct%cblk(j) ! block corresponding to contact j
        call createAll(Gam,SelfEneR(j)%nrow,SelfEneR(j)%ncol)
        call spectral_gpu(SelfEner(j)%val, Gam%val)
        call add_cublas(hh, Sigma_n(cb,cb)%val, Sigma_n(cb,cb)%val, frmdiff, Gam%val)
        call destroyAll(Gam)
      endif
    end do

    call calculate_sigma_n()

    call createAll(work1,Sigma_n(1,1)%nrow, Gr(1,1)%nrow)
    call matmul_gpu(hh, one, Sigma_n(1,1)%val, Gr(1,1)%val, zero, work1%val, 'dag_2nd')
    call matmul_gpu(hh, one, Gr(1,1)%val, work1%val, one, Gn(1,1)%val)
    call destroyAll(work1)
    
    if (nbl .eq. 1) return

    !Explicit formulae:
    !Gn(i+1,i) = gsmr(i+1)*[Sigma(i+1,i)Ga(i,i) + Sigma(i+1,i+1)Ga(i+1,i) - Tr(i+1,i)Gn(i,i)]
    !Gn(i,i+1) = [Gr(i,i)Sigma(i,i+1) + Gr(i,i+1)Sigma(i+1,i+1) - Gn(i,i)Ta(i,i+1)] * gsma(i+1)
    !Use Hermitian property of Gn:
    !Gn(i,i+1) = Gn(i+1,i)^dag
    !Gn(i+1,i+1) = gsmr(i+1) * [Sigma(i+1,i)Ga(i,i+1) + Sigma(i+1,i+1)Ga(i+1,i+1) - Tr(i+1,i)Gn(i,i+1)]
    !Implementation exploits cumulative sum of prealloc_mult, C = C + A*B

    do i = 1, nbl-1

        call createAll(work1,Sigma_n(i+1,i)%nrow,Gr(i,i)%nrow)
        call matmul_gpu(hh, one, Sigma_n(i+1,i)%val, Gr(i,i)%val, zero, work1%val, 'dag_2nd')
        call matmul_gpu(hh, one, Sigma_n(i+1,i+1)%val, Gr(i,i+1)%val, one, work1%val, 'dag_2nd')

        call copyToGPU(ESH(i+1,i))
        call matmul_gpu(hh, minusOne, ESH(i+1,i)%val, Gn(i,i)%val, one, work1%val)
        call deleteGPU(ESH(i+1,i))

        call matmul_gpu(hh, one, gsmr(i+1)%val, work1%val, zero, Gn(i+1,i)%val)
        call destroyAll(work1)

        call dagger_gpu(Gn(i+1,i)%val,Gn(i,i+1)%val)

        call createAll(work1, Sigma_n(i+1,i)%nrow, Gr(i+1,i)%nrow)
        call matmul_gpu(hh, one,Sigma_n(i+1,i)%val, Gr(i+1,i)%val, zero, work1%val, 'dag_2nd')

        call matmul_gpu(hh, one,Sigma_n(i+1,i+1)%val, Gr(i+1,i+1)%val, one, work1%val, 'dag_2nd')

        call matmul_gpu(hh, minusOne, ESH(i+1,i)%val, Gn(i,i+1)%val, one, work1%val)


        call matmul_gpu(hh, one,gsmr(i+1)%val, work1%val, zero, Gn(i+1,i+1)%val)
        call destroyAll(work1)

    end do

    call destroy_tridiag_blk(Sigma_n)
    call deallocate_blk_dns(Sigma_n)
 
    contains
    ! Recursive calculation of Sigma_n:      
    ! gns(i+1) = gsmr(i+1) Sigma(i+1,i+1) gsmr(i+1)^dag
    ! Sigma(i,i) = Sigma(i,i) + Tr(i,i+1) gns(i+1) Ta(i+1,i) 
    !                         - Tr(i,i+1) gsmr(i+1) Sigma(i+1,i) 
    !                         - Sigma(i,i+1) gsmr^dag(i+1) Ta(i+1,i)]
    ! 
    subroutine calculate_sigma_n()      
      !Work
      type(z_DNS) :: work, gns 

      ! if nbl = 1 => Sigma_n(1,1) is ready
      if (nbl.eq.1) return
      !g^n(nbl) = gsmr(nbl) Sigma(nbl,nbl) gsma(nbl)
      call createAll(work1, gsmr(nbl)%nrow, Sigma_n(nbl,nbl)%ncol)
      call createAll(gns, work1%nrow, gsmr(nbl)%nrow)
      call matmul_gpu(hh, one, gsmr(nbl)%val, Sigma_n(nbl,nbl)%val, zero, work1%val)
      call matmul_gpu(hh, one, work1%val, gsmr(nbl)%val, zero, gns%val, 'dag_2nd')
      call destroyAll(work1)

      do i = nbl-1, 1, -1
        !work1 = Tr(i,i+1) gns(i+1) Ta(i+1,i)
        ! Tr(i,i+1) = ESH(i,i+1);  Ta(i+1,i) = ESH(i,i+1)^dag
        call createAll(work, ESH(i,i+1)%nrow, gns%ncol)

        call matmul_gpu(hh, one, ESH(i,i+1)%val, gns%val, zero, work%val)
        call matmul_gpu(hh, one, work%val, ESH(i,i+1)%val, one, Sigma_n(i,i)%val, 'dag_2nd')

        call destroyAll(work)
        call destroyAll(gns)

        !work2 = Sigma(i,i+1) gsmr^dag(i+1) Ta(i+1,i)
        call createAll(work, Sigma_n(i,i+1)%nrow, gsmr(i+1)%ncol)

        call matmul_gpu(hh, one, Sigma_n(i,i+1)%val, gsmr(i+1)%val, zero, work%val, 'dag_2nd')
        call matmul_gpu(hh, minusOne, work%val, ESH(i,i+1)%val, one, Sigma_n(i,i)%val, 'dag_2nd')

        call destroyAll(work)

        !work3 = ESH(i,i+1) gsmr(i+1) Sigma(i+1,i)
        call createAll(work, ESH(i,i+1)%nrow, gsmr(i+1)%ncol)

        call matmul_gpu(hh, one, ESH(i,i+1)%val, gsmr(i+1)%val, zero, work%val)
        call matmul_gpu(hh, minusOne, work%val, Sigma_n(i+1,i)%val, one, Sigma_n(i,i)%val)

        call destroyAll(work)

        if (i > 1) then
          !gns(i) = gsmr(i) * Sigma_n(i,i) * gsmr^dag(i)
          call createAll(work, gsmr(i)%nrow, Sigma_n(i,i)%ncol)
          call createAll(gns, work%nrow, gsmr(i)%nrow)

          call matmul_gpu(hh, one, gsmr(i)%val, Sigma_n(i,i)%val, zero, work%val)
          call matmul_gpu(hh, one, work%val, gsmr(i)%val, zero, gns%val, 'dag_2nd')

          call destroyAll(work)
        end if

      end do

    end subroutine calculate_sigma_n

  end subroutine calculate_Gn_tridiag_blocks

!  subroutine calculate_Gn_tridiag_blocks_sp(negf, ESH, SelfEneR, frm, ref, struct, gsml, gsmr, Gr, Gn) 
!    !In/Out
!    type(Tnegf), intent(in) :: negf
!    type(c_DNS), dimension(:), intent(in) :: gsmr, gsml
!    type(c_DNS), dimension(:,:), intent(in) :: ESH, Gr
!    type(c_DNS), dimension(:,:), intent(inout) :: Gn
!    type(c_DNS), dimension(:), intent(in) :: SelfEneR
!    type(Tstruct_info), intent(in) :: struct
!    real(sp), dimension(:), intent(in) :: frm
!    integer, intent(in) :: ref
!
!    !Work
!    type(CublasHandle) :: hh
!    Type(c_DNS) :: Gam, GA
!    type(c_DNS) :: work1, work2, work3
!    integer :: i,j
!    integer :: cb, nbl, ncont, istat
!    complex(sp) :: frmdiff
!    complex(sp), parameter :: one = (1.0_sp, 0.0_sp)
!    complex(sp), parameter :: mone = (-1.0_sp, 0.0_sp)
!    complex(sp), parameter :: zero = (0.0_sp, 0.0_sp)
!    real(sp) :: somma
!
!    ncont = struct%num_conts
!    nbl = struct%num_PLs
!    hh = negf%hcublas
!
!    do j=1,ncont
!
!       if (j.NE.ref .AND. ABS(frm(j)-frm(ref)).GT.EPS) THEN
!
!          cb=struct%cblk(j) ! block corresponding to contact j    
!
!          !Creating Gam = i (SE -  SE^+)
!          call createAll(Gam, SelfEneR(j)%nrow, SelfEneR(j)%ncol)
!          call spectral_gpu(SelfEneR(j)%val,Gam%val)
!
!          frmdiff = cmplx(frm(j)-frm(ref),0.0_dp,dp)
!
!          ! Computation of Gn(cb,cb) = Gr(cb,cb) Gam(cb) Gr(cb,cb)^+
!          call createAll(work1, Gam%nrow, Gr(cb,cb)%nrow) !it's Gr^+ --> I take nrow
!          call matmul_gpu(hh, frmdiff, Gam%val, Gr(cb,cb)%val, zero, work1%val, 'dag_2nd')
!          call matmul_gpu(hh, one, Gr(cb,cb)%val, work1%val, zero, Gn(cb,cb)%val)
!          call destroyAll(work1)
!          call destroyAll(GA)
!
!
!          !***************************************************************
!          !*** Gr(cb-1, cb) already allocated: must insert if in loop ****
!          !***************************************************************
!
!          do i=cb-1, 1, -1
!             !Gr(i,cb) = - gL(i) ESH(i,i+1) Gr(i+1,cb)
!             if (i .eq. cb-1) then
!                call createAll(work1, ESH(i,i+1)%nrow, Gr(i+1,cb)%ncol)
!             else 
!                call createAll(work1, ESH(i,i+1)%nrow, Gr(i+1,cb)%ncol)   
!                call createAll(Gr(i,cb), gsml(i)%nrow, work1%ncol)
!             end if
!             call copyToGPU(ESH(i,i+1))
!             call matmul_gpu(hh, one, ESH(i,i+1)%val, Gr(i+1,cb)%val,zero, work1%val)
!             call deleteGPU(ESH(i,i+1))
!             call matmul_gpu(hh, mone, gsml(i)%val, work1%val, zero, Gr(i,cb)%val)
!             call destroyAll(work1)
!
!             !Gn(i, i) = Gr(i, cb) Gam(cb) Gr(i, cb)^+
!             call createAll(work2, Gam%nrow, Gr(i,cb)%nrow)
!             call matmul_gpu(hh, one, Gam%val, Gr(i,cb)%val, zero, work2%val, 'dag_2nd')
!             call matmul_gpu(hh, frmdiff, Gr(i,cb)%val, work2%val, zero, Gn(i,i)%val)
!
!             !Gn(i+1,i)  = Gr(i+1, cb) Gam(cb) Gr(i, cb)^+
!             call matmul_gpu(hh, frmdiff, Gr(i+1,cb)%val, work2%val, zero, Gn(i+1,i)%val)
!             call destroyAll(work2)
!
!             !Gn(i,i+1)  = Gr(i, cb) Gam(cb) Gr(i+1, cb)^+
!             call createAll(work3, Gam%nrow, Gr(i+1,cb)%nrow)
!             call matmul_gpu(hh, frmdiff, Gam%val, Gr(i+1,cb)%val, zero, work3%val, 'dag_2nd')
!             call matmul_gpu(hh, one, Gr(i,cb)%val, work3%val, zero, Gn(i,i+1)%val)
!             call destroyAll(work3)
!          end do
!
!          !***************************************************************
!          !*** Gr(cb+1, cb) already allocated: must insert if in loop ****
!          !***************************************************************
!
!          do i=cb+1, nbl
!             !Gr(i,cb) = - gR(i) ESH(i,i-1) Gr(i-1,cb)
!             if (i .eq. cb+1) then
!                call createAll(work1, ESH(i,i-1)%nrow, Gr(i-1,cb)%ncol)
!
!             else 
!                call createAll(work1, ESH(i,i-1)%nrow, Gr(i-1,cb)%ncol)
!                call createAll(Gr(i,cb), gsmr(i)%nrow, work1%ncol)
!             end if
!             call copyToGPU(ESH(i,i-1))
!             call matmul_gpu(hh, one, ESH(i,i-1)%val, Gr(i-1,cb)%val,zero, work1%val)
!             call deleteGPU(ESH(i,i-1))
!             call matmul_gpu(hh, mone, gsmr(i)%val, work1%val, zero, Gr(i,cb)%val)
!             call destroyAll(work1)
!
!
!             !Gn(i, i) = Gr(i, cb) Gam(cb) Gr(i, cb)^+
!             call createAll(work2, Gam%nrow, Gr(i,cb)%nrow)
!             call matmul_gpu(hh, one, Gam%val, Gr(i,cb)%val, zero, work2%val, 'dag_2nd')
!             call matmul_gpu(hh, frmdiff, Gr(i,cb)%val, work2%val, zero, Gn(i,i)%val)
!
!             !Gn(i-1,i)  = Gr(i-1, cb) Gam(cb) Gr(i, cb)^+
!             call matmul_gpu(hh, frmdiff, Gr(i-1,cb)%val, work2%val, zero, Gn(i-1,i)%val)
!             call destroyAll(work2)
!
!             !Gn(i,i-1)  = Gr(i, cb) Gam(cb) Gr(i-1, cb)^+
!             call createAll(work3, Gam%nrow, Gr(i-1,cb)%nrow)
!             call matmul_gpu(hh, one, Gam%val, Gr(i-1,cb)%val, zero, work3%val, 'dag_2nd')
!             call matmul_gpu(hh, frmdiff, Gr(i,cb)%val, work3%val, zero, Gn(i,i-1)%val)
!             call destroyAll(work3)
!          end do
!
!          do i=cb-2, 1, -1
!             call destroyAll(Gr(i,cb))
!          end do
!
!          do i=cb+2, nbl
!             call destroyAll(Gr(i,cb))
!          end do
!          call destroyAll(Gam)
!       endif
!    end do
!
!  end subroutine calculate_Gn_tridiag_blocks_sp
!
!  subroutine calculate_Gn_tridiag_blocks_dp(negf, ESH, SelfEneR, frm, ref, struct, gsml, gsmr, Gr, Gn)
!    !In/Out
!    type(Tnegf) :: negf
!    type(z_DNS), dimension(:,:), intent(in) :: ESH, Gr
!    type(z_DNS), dimension(:,:), intent(inout) :: Gn
!    type(z_DNS), dimension(:), intent(in) :: gsmr, gsml
!    type(z_DNS), dimension(:), intent(in) :: SelfEneR
!    type(Tstruct_info), intent(in) :: struct
!    real(dp), dimension(:), intent(in) :: frm
!    integer, intent(in) :: ref
!
!    !Work
!    type(CublasHandle) :: hh
!    Type(z_DNS) :: Gam, GA
!    type(z_DNS) :: work1, work2, work3
!    integer :: i, j, istat
!    integer :: cb, nbl, ncont
!    complex(dp) :: frmdiff
!    complex(dp), parameter :: one = (1.0_dp, 0.0_dp)
!    complex(dp), parameter :: mone = (-1.0_dp, 0.0_dp)
!    complex(dp), parameter :: zero = (0.0_dp, 0.0_dp)
!    real(dp) :: somma
!
!    ncont = struct%num_conts
!    nbl = struct%num_PLs
!    hh = negf%hcublas
!
!    do j=1,ncont
!
!       if (j.NE.ref .AND. ABS(frm(j)-frm(ref)).GT.EPS) THEN
!
!          cb=struct%cblk(j) ! block corresponding to contact j    
!
!          !Creating Gam = i (SE -  SE^+)
!          call createAll(Gam, SelfEneR(j)%nrow, SelfEneR(j)%ncol)
!          call spectral_gpu(SelfEneR(j)%val,Gam%val)
!
!          frmdiff = cmplx(frm(j)-frm(ref),0.0_dp,dp)
!
!          ! Computation of Gn(cb,cb) = Gr(cb,cb) Gam(cb) Gr(cb,cb)^+
!          call createAll(work1, Gam%nrow, Gr(cb,cb)%nrow) !it's Gr^+ --> I take nrow
!          call matmul_gpu(hh, frmdiff, Gam%val, Gr(cb,cb)%val, zero, work1%val, 'dag_2nd')
!          call matmul_gpu(hh, one, Gr(cb,cb)%val, work1%val, zero, Gn(cb,cb)%val)
!          call destroyAll(work1)
!          call destroyAll(GA)
!
!
!          !***************************************************************
!          !*** Gr(cb-1, cb) already allocated: must insert if in loop ****
!          !***************************************************************
!
!          do i=cb-1, 1, -1
!             !Gr(i,cb) = - gL(i) ESH(i,i+1) Gr(i+1,cb)
!             if (i .eq. cb-1) then
!                call createAll(work1, ESH(i,i+1)%nrow, Gr(i+1,cb)%ncol)
!             else 
!                call createAll(work1, ESH(i,i+1)%nrow, Gr(i+1,cb)%ncol)   
!                call createAll(Gr(i,cb), gsml(i)%nrow, work1%ncol)
!             end if
!             call copyToGPU(ESH(i,i+1))
!             call matmul_gpu(hh, one, ESH(i,i+1)%val, Gr(i+1,cb)%val,zero, work1%val)
!             call deleteGPU(ESH(i,i+1))
!             call matmul_gpu(hh, mone, gsml(i)%val, work1%val, zero, Gr(i,cb)%val)
!             call destroyAll(work1)
!
!             !Gn(i, i) = Gr(i, cb) Gam(cb) Gr(i, cb)^+
!             call createAll(work2, Gam%nrow, Gr(i,cb)%nrow)
!             call matmul_gpu(hh, one, Gam%val, Gr(i,cb)%val, zero, work2%val, 'dag_2nd')
!             call matmul_gpu(hh, frmdiff, Gr(i,cb)%val, work2%val, zero, Gn(i,i)%val)
!
!             !Gn(i+1,i)  = Gr(i+1, cb) Gam(cb) Gr(i, cb)^+
!             call matmul_gpu(hh, frmdiff, Gr(i+1,cb)%val, work2%val, zero, Gn(i+1,i)%val)
!             call destroyAll(work2)
!
!             !Gn(i,i+1)  = Gr(i, cb) Gam(cb) Gr(i+1, cb)^+
!             call createAll(work3, Gam%nrow, Gr(i+1,cb)%nrow)
!             call matmul_gpu(hh, frmdiff, Gam%val, Gr(i+1,cb)%val, zero, work3%val, 'dag_2nd')
!             call matmul_gpu(hh, one, Gr(i,cb)%val, work3%val, zero, Gn(i,i+1)%val)
!             call destroyAll(work3)
!          end do
!
!          !***************************************************************
!          !*** Gr(cb+1, cb) already allocated: must insert if in loop ****
!          !***************************************************************
!
!          do i=cb+1, nbl
!             !Gr(i,cb) = - gR(i) ESH(i,i-1) Gr(i-1,cb)
!             if (i .eq. cb+1) then
!                call createAll(work1, ESH(i,i-1)%nrow, Gr(i-1,cb)%ncol)
!
!             else 
!                call createAll(work1, ESH(i,i-1)%nrow, Gr(i-1,cb)%ncol)
!                call createAll(Gr(i,cb), gsmr(i)%nrow, work1%ncol)
!             end if
!             call copyToGPU(ESH(i,i-1))
!             call matmul_gpu(hh, one, ESH(i,i-1)%val, Gr(i-1,cb)%val,zero, work1%val)
!             call deleteGPU(ESH(i,i-1))
!             call matmul_gpu(hh, mone, gsmr(i)%val, work1%val, zero, Gr(i,cb)%val)
!             call destroyAll(work1)
!
!
!             !Gn(i, i) = Gr(i, cb) Gam(cb) Gr(i, cb)^+
!             call createAll(work2, Gam%nrow, Gr(i,cb)%nrow)
!             call matmul_gpu(hh, one, Gam%val, Gr(i,cb)%val, zero, work2%val, 'dag_2nd')
!             call matmul_gpu(hh, frmdiff, Gr(i,cb)%val, work2%val, zero, Gn(i,i)%val)
!
!             !Gn(i-1,i)  = Gr(i-1, cb) Gam(cb) Gr(i, cb)^+
!             call matmul_gpu(hh, frmdiff, Gr(i-1,cb)%val, work2%val, zero, Gn(i-1,i)%val)
!             call destroyAll(work2)
!
!             !Gn(i,i-1)  = Gr(i, cb) Gam(cb) Gr(i-1, cb)^+
!             call createAll(work3, Gam%nrow, Gr(i-1,cb)%nrow)
!             call matmul_gpu(hh, one, Gam%val, Gr(i-1,cb)%val, zero, work3%val, 'dag_2nd')
!             call matmul_gpu(hh, frmdiff, Gr(i,cb)%val, work3%val, zero, Gn(i,i-1)%val)
!             call destroyAll(work3)
!          end do
!
!          do i=cb-2, 1, -1
!             call destroyAll(Gr(i,cb))
!          end do
!
!          do i=cb+2, nbl
!             call destroyAll(Gr(i,cb))
!          end do
!          call destroyAll(Gam)
!       endif
!    end do
!
!  end subroutine calculate_Gn_tridiag_blocks_dp

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
    call createAll(GA, Gr(bl1,bl1)%ncol, Gr(bl1,bl1)%nrow)
    call dagger_gpu(Gr(bl1,bl1)%val,GA%val)

    ! Computes the Gamma matrices
    call createAll(GAM1_dns, SelfEneR(ct1)%nrow, SelfEneR(ct1)%ncol)
    call spectral_gpu(SelfEneR(ct1)%val,GAM1_dns%val)

    ! Work to compute transmission matrix (Gamma G Gamma G)
    call createAll(work1, GAM1_dns%nrow, Gr(bl1,bl1)%ncol)
    call createAll(work2, work1%nrow, GAM1_dns%ncol)
    call matmul_gpu(hh, one, GAM1_dns%val, Gr(bl1,bl1)%val, zero, work1%val)
    call matmul_gpu(hh, one, work1%val, GAM1_dns%val, zero, work2%val)

    call destroyAll(work1)

    call createAll(work1, work2%nrow, GA%ncol)
    call matmul_gpu(hh, one, work2%val, GA%val, zero, work1%val)
    call destroyAll(work2)

    call createAll(AA, GA%nrow, GA%ncol)
    call add_gpu(j,Gr(bl1,bl1)%val,mj,GA%val,AA%val)
    call destroyAll(GA)

    call createAll(work2, GAM1_dns%nrow, AA%ncol)
    call matmul_gpu(hh, one, GAM1_dns%val, AA%val, zero, work2%val)
    call destroyAll(GAM1_dns)
    call destroyAll(AA)

    call createAll(TRS, work1%nrow, work1%ncol)
    call add_gpu(one,work2%val,mone,work1%val,TRS%val)
    call get_tun_mask(ESH, bl1, tun_proj, tun_mask)

    TUN = abs(real(trace_gpu(TRS%val, tun_mask)))
    call deleteGPU(tun_mask)
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
    complex(dp), parameter ::    j = (0.0_dp,1.0_dp)  ! CMPX unity
    complex(dp), parameter ::    mj = (0.0_dp,-1.0_dp)  
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
    call createAll(GA, Gr(bl1,bl1)%ncol, Gr(bl1,bl1)%nrow)
    call dagger_gpu(Gr(bl1,bl1)%val,GA%val)

    ! Computes the Gamma matrices
    call createAll(GAM1_dns, SelfEneR(ct1)%nrow, SelfEneR(ct1)%ncol)
    call spectral_gpu(SelfEneR(ct1)%val,GAM1_dns%val)

    ! Work to compute transmission matrix (Gamma G Gamma G)
    call createAll(work1, GAM1_dns%nrow, Gr(bl1,bl1)%ncol)
    call createAll(work2, work1%nrow, GAM1_dns%ncol)
    call matmul_gpu(hh, one, GAM1_dns%val, Gr(bl1,bl1)%val, zero, work1%val)
    call matmul_gpu(hh, one, work1%val, GAM1_dns%val, zero, work2%val)

    call destroyAll(work1)

    call createAll(work1, work2%nrow, GA%ncol)
    call matmul_gpu(hh, one, work2%val, GA%val, zero, work1%val)
    call destroyAll(work2)

    call createAll(AA, GA%nrow, GA%ncol)
    call add_gpu(j,Gr(bl1,bl1)%val,mj,GA%val,AA%val)
    call destroyAll(GA)

    call createAll(work2, GAM1_dns%nrow, AA%ncol)
    call matmul_gpu(hh, one, GAM1_dns%val, AA%val, zero, work2%val)
    call destroyAll(GAM1_dns)
    call destroyAll(AA)

    call createAll(TRS, work1%nrow, work1%ncol)
    call add_gpu(one,work2%val,mone,work1%val,TRS%val)
    call get_tun_mask(ESH, bl1, tun_proj, tun_mask)

    TUN = abs(real(trace_gpu(TRS%val, tun_mask)))
    call deleteGPU(tun_mask)
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
    Type(c_DNS) :: work1, work2, GAM1_dns, GAM2_dns, GA, TRS
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
             call matmul_gpu(hh, mone, gsmr(i)%val, ESH(i,i-1)%val, zero, work1%val)
             call deleteGPU(ESH(i,i-1))

             call matmul_gpu(hh, one, work1%val, Gr(i-1,bl1)%val ,zero, Gr(i,bl1)%val)
             call destroyAll(work1)

          endif

          ! avoid destroying blocks closer to diagonal
          if (i.gt.(bl1+2)) call destroyAll(Gr(i-1,bl1))
       end do

    endif
    ! Computes the Gamma matrices
    call createAll(GAM1_dns, SelfEneR(ct1)%nrow, SelfEneR(ct1)%ncol)
    call createAll(GAM2_dns, SelfEneR(ct2)%nrow, SelfEneR(ct2)%ncol)
    call spectral_gpu(SelfEneR(ct1)%val,GAM1_dns%val)
    call spectral_gpu(SelfEneR(ct2)%val,GAM2_dns%val)

    ! Work to compute transmission matrix (Gamma2 Gr Gamma1 Ga)
    call createAll(work1, Gr(bl2,bl1)%nrow, GAM1_dns%ncol)
    call matmul_gpu(hh, one, Gr(bl2,bl1)%val, GAM1_dns%val, zero, work1%val)

    call createAll(work2, GAM2_dns%nrow, work1%ncol)
    call matmul_gpu(hh, one, GAM2_dns%val, work1%val, zero, work2%val)

    call destroyAll(work1)
    call destroyAll(GAM2_dns)
    call destroyAll(GAM1_dns)

    call createAll(TRS, work2%nrow, Gr(bl2,bl1)%nrow)

    call createAll(GA, Gr(bl2,bl1)%ncol, Gr(bl2,bl1)%nrow)
    call dagger_gpu(Gr(bl2,bl1)%val,GA%val)
    call matmul_gpu(hh, one, work2%val, GA%val, zero, TRS%val)
    call destroyAll(work2)
    call destroyAll(GA)
    if (bl2.gt.bl1+1) call destroyAll(Gr(bl2,bl1))

    call get_tun_mask(ESH, bl2, tun_proj, tun_mask)
    TUN = abs(real(trace_gpu(TRS%val, tun_mask)))
    
    call deleteGPU(tun_mask) 
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
    Type(z_DNS) :: work1, work2, GAM1_dns, GAM2_dns, GA, TRS
    real(dp) :: max, summ
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
             call matmul_gpu(hh, mone, gsmr(i)%val, ESH(i,i-1)%val, zero, work1%val)
             call deleteGPU(ESH(i,i-1))

             call matmul_gpu(hh, one, work1%val, Gr(i-1,bl1)%val ,zero, Gr(i,bl1)%val)
             call destroyAll(work1)

          endif

          ! avoid destroying blocks closer to diagonal
          if (i.gt.(bl1+2)) call destroyAll(Gr(i-1,bl1))
       end do

    endif
    ! Computes the Gamma matrices
    call createAll(GAM1_dns, SelfEneR(ct1)%nrow, SelfEneR(ct1)%ncol)
    call createAll(GAM2_dns, SelfEneR(ct2)%nrow, SelfEneR(ct2)%ncol)
    call spectral_gpu(SelfEneR(ct1)%val,GAM1_dns%val)
    call spectral_gpu(SelfEneR(ct2)%val,GAM2_dns%val)

    ! Work to compute transmission matrix (Gamma2 Gr Gamma1 Ga)
    call createAll(work1, Gr(bl2,bl1)%nrow, GAM1_dns%ncol)
    call matmul_gpu(hh, one, Gr(bl2,bl1)%val, GAM1_dns%val, zero, work1%val)

    call createAll(work2, GAM2_dns%nrow, work1%ncol)
    call matmul_gpu(hh, one, GAM2_dns%val, work1%val, zero, work2%val)

    call destroyAll(work1)
    call destroyAll(GAM2_dns)
    call destroyAll(GAM1_dns)

    call createAll(TRS, work2%nrow, Gr(bl2,bl1)%nrow)

    call createAll(GA, Gr(bl2,bl1)%ncol, Gr(bl2,bl1)%nrow)
    call dagger_gpu(Gr(bl2,bl1)%val,GA%val)
    call matmul_gpu(hh, one, work2%val, GA%val, zero, TRS%val)
    call destroyAll(work2)
    call destroyAll(GA)
    if (bl2.gt.bl1+1) call destroyAll(Gr(bl2,bl1))

    call get_tun_mask(ESH, bl2, tun_proj, tun_mask)
    TUN = abs(real(trace_gpu(TRS%val, tun_mask)))
    
    call deleteGPU(tun_mask) 
    call log_deallocate(tun_mask)

    call destroyAll(TRS)

  end subroutine calculate_single_transmission_N_contacts_dp
  
! Based on projection indices build a logical mask just on contact block 
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
       !$acc data present(ESH) 
       do ii = 1, nbl-1
          istart = istart + ESH(ii,ii)%nrow
       end do
       !$acc end data 
       iend = istart + ESH(nbl,nbl)%nrow + 1  

       ! select the indices in tun_proj 
       
       !$acc kernels copyin(tun_proj, tun_mask)   
       do ii = 1, size(tun_proj%indexes)
          ind = tun_proj%indexes(ii)
          if (ind > istart .and. ind < iend) then
             tun_mask(ind - istart) = .true. 
          end if
       end do
       !$acc end kernels
       !$acc exit data delete(tun_proj)
    else
       tun_mask = .true.
       !$acc enter data copyin(tun_mask)
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
       !$acc data present(ESH) 
       do ii = 1, nbl-1
          istart = istart + ESH(ii,ii)%nrow
       end do
       iend = istart + ESH(nbl,nbl)%nrow + 1  
       !$acc end data 

       ! select the indices in tun_proj

       !$acc kernels copyin(tun_proj, tun_mask)   
       do ii = 1, size(tun_proj%indexes)
          ind = tun_proj%indexes(ii)
          if (ind > istart .and. ind < iend) then
             tun_mask(ind - istart) = .true. 
          end if
       end do
       !$acc end kernels
       !$acc exit data delete(tun_proj)
    else
       tun_mask = .true.
       !$acc enter data copyin(tun_mask)
    end if

  end subroutine get_tun_mask_dp

!---------------------------------------------------

  subroutine init_tridiag_blk(Matrix,S)
    type(z_DNS), dimension(:,:) :: Matrix,S

    integer :: nbl, j

    nbl = SIZE(Matrix,1)

    call create(Matrix(1,1),S(1,1)%nrow,S(1,1)%ncol)
    Matrix(1,1)%val=(0.0_dp,0.0_dp)
    do j=2,nbl-1
       call create(Matrix(j-1,j),S(j-1,j)%nrow,S(j-1,j)%ncol)
       Matrix(j-1,j)%val=(0.0_dp,0.0_dp)
      
       call create(Matrix(j,j),S(j,j)%nrow,S(j,j)%ncol)
       Matrix(j,j)%val=(0.0_dp,0.0_dp)
      
       call create(Matrix(j,j-1),S(j,j-1)%nrow,S(j,j-1)%ncol)
       Matrix(j,j-1)%val=(0.0_dp,0.0_dp)
    end do
    if (nbl.gt.1) then
       call create(Matrix(nbl,nbl),S(nbl,nbl)%nrow,S(nbl,nbl)%ncol)
       Matrix(nbl,nbl)%val=(0.0_dp,0.0_dp)

       call create(Matrix(nbl-1,nbl),S(nbl-1,nbl)%nrow,S(nbl-1,nbl)%ncol)
       Matrix(nbl-1,nbl)%val=(0.0_dp,0.0_dp)

       call create(Matrix(nbl,nbl-1),S(nbl,nbl-1)%nrow,S(nbl,nbl-1)%ncol)
       Matrix(nbl,nbl-1)%val=(0.0_dp,0.0_dp)
    endif

  end subroutine init_tridiag_blk

!---------------------------------------------------

  subroutine destroy_tridiag_blk(M)
    type(z_DNS), dimension(:,:), allocatable :: M 

    integer :: i, nbl

    if (.not.allocated(M)) return

    nbl=size(M,1)

    do i=1,nbl
      if (allocated(M(i,i)%val)) then
        call destroyAll(M(i,i))
      end if  
    end do
    do i=2,nbl
      if (allocated(M(i-1,i)%val)) then
        call destroyAll(M(i-1,i))
      end if
      if (allocated(M(i,i-1)%val)) then
        call destroyAll(M(i,i-1))
      end if  
    end do

  end subroutine destroy_tridiag_blk

!---------------------------------------------------

  subroutine allocate_blk_dns(blkM,nbl)
    type(z_DNS), dimension(:,:), allocatable :: blkM
    integer :: nbl, ierr

    allocate(blkM(nbl,nbl),stat=ierr)
    if (ierr.ne.0) stop 'ALLOCATION ERROR: could not allocate block-Matrix'

  end subroutine allocate_blk_dns

  !---------------------------------------------------

  subroutine deallocate_blk_dns(blkM)
    type(z_DNS), dimension(:,:), allocatable :: blkM
    integer :: ierr

    deallocate(blkM,stat=ierr)
    if (ierr.ne.0) stop 'DEALLOCATION ERROR: could not deallocate block-Matrix'

  end subroutine deallocate_blk_dns

end module iterative_gpu
