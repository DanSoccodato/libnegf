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


module iterative

  use ln_precision
  use ln_constants, only : pi, i_unit => j 
  use ln_allocation
  use mat_def
  use sparsekit_drv
  use inversions
  use elph
  use ln_structure, only : TStruct_Info
  use lib_param, only : MAXNCONT, Tnegf, intarray
  use mpi_globals, only : id, numprocs, id0
  use outmatrix, only : outmat_c, inmat_c, direct_out_c, direct_in_c
  use clock
#:if defined("GPU")
  use cudautils
  use cublas_v2
  use cusolverDn
  use iterative_gpu
#:else
  use iterative_cpu
#:endif
  !use transform

  implicit none
  private

  public :: calculate_transmissions
  public :: calculate_transmissions_and_dos

  public :: calculate_Gr
  public :: calculate_Gn_neq_components

  public :: calculate_Gr_outer
  public :: calculate_Gn_outer

  public :: iterative_meir_wingreen
  public :: iterative_layer_current
  public :: transmission_BP_corrected

  public :: destroy_all_blk

  logical, parameter :: debug=.false.
  ! These are here temporarily ...
  type(z_DNS), dimension(:), allocatable :: gsmr
  !type(z_DNS), dimension(:), allocatable :: gsml
  type(z_DNS), dimension(:,:), allocatable :: Gr
  type(z_DNS), dimension(:,:), allocatable :: ESH
  type(z_DNS), dimension(:,:), allocatable :: Gn


CONTAINS

  !****************************************************************************
  !
  ! Driver for computing Equilibrium Green Retarded (Outer blocks included)
  ! writing on memory
  !
  !****************************************************************************

  subroutine calculate_Gr(negf,E,SelfEneR,Tlc,Tcl,gsurfR,Grout,outer)

    !****************************************************************************
    !
    !Input
    !negf:    negf data container
    !E:        Energy point
    !SelfEneR: matrices array containing contacts Self Energy
    !Tlc:      matrices array containing contacts-device interaction blocks (ES-H)
    !Tcl:      matrices array containing device-contacts interaction blocks (ES-H)
    !gsurfR:   matrices array containing contacts surface green
    !outer:    optional parameter (0,1,2).
    !
    !Output:
    !Grout: Retarded Green's function (Device + Contacts overlap regions -> effective conductor)
    !   outer = 0  no outer parts are computed
    !   outer = 1  only D/C part is computed
    !   outer = 2  D/C and C/D parts are computed
    !              (needed for K-points calculations)
    !
    !*****************************************************************************

    type(Tnegf), intent(inout) :: negf
    complex(dp), intent(in) :: E
    type(z_DNS), dimension(:), intent(in) :: SelfEneR
    type(z_DNS), dimension(:), intent(inout) :: Tlc, Tcl, gsurfR
    type(z_CSR), intent(out) :: Grout
    integer, intent(in) :: outer

    !Work
    type(z_CSR) :: ESH_tot, Ain
    integer :: i,ierr, nbl, ncont,ii,n

    nbl = negf%str%num_PLs
    ncont = negf%str%num_conts

    call build_ESH(negf, E, SelfEneR, ncont, nbl, ESH)

    !! Add interaction self energy contribution, if any
    if (allocated(negf%inter)) call negf%inter%add_sigma_r(ESH)

    call allocate_gsm(gsmr,nbl)
    call calculate_gsmr_blocks(negf,ESH,nbl,2,gsmr)

    call allocate_blk_dns(Gr,nbl)

    call calculate_Gr_tridiag_blocks(negf,ESH,gsmr,Gr,1)
    call calculate_Gr_tridiag_blocks(negf,ESH,gsmr,Gr,2,nbl)

#:if defined("GPU")
    call copy_trid_toHOST(Gr) 
#:endif
    
    call destroy_tridiag_blk(ESH)
    call deallocate_blk_dns(ESH)

    call destroy_gsm(gsmr)
    call deallocate_gsm(gsmr)


    !! Deliver Gr to interaction models if any
    if (allocated(negf%inter)) call negf%inter%set_Gr(Gr, negf%iE)
    !-----------------------------------------------------------

    call blk2csr(Gr,negf%str,negf%S,Grout)

    SELECT CASE (outer)
    CASE(0)
    CASE(1)
      call calculate_Gr_outer(Tlc,Tcl,gsurfR,negf%str,.FALSE.,Grout)
    CASE(2)
      call calculate_Gr_outer(Tlc,Tcl,gsurfR,negf%str,.TRUE.,Grout)
    end SELECT

    !Distruzione dell'array Gr
    call destroy_blk(Gr)
    DEALLOCATE(Gr)

  end subroutine calculate_Gr


  !****************************************************************************
  !
  ! Driver for computing G_n contributions due to all contacts MINUS reference:
  ! Reference is necessary when splitting into contour + real-axis integration
  !
  !   Sum   [f_j(E)-f_r(E)] Gr Gam_j Ga
  !   j!=r
  !
  ! NOTE: Setting reference to ncont+1 removes the reference part
  !****************************************************************************

  subroutine calculate_Gn_neq_components(negf,E,SelfEneR,Tlc,Tcl,gsurfR,frm,Glout,outblocks)

    !****************************************************************************
    !
    !Input
    !negf:    negf data container
    !E:        Energy point
    !SelfEneR: matrices array containing contacts Self Energy
    !Tlc:      matrices array containing contacts-device interaction blocks (ES-H)
    !Tcl:      matrices array containing device-contacts interaction blocks (ES-H)
    !gsurfR:   matrices array containing contacts surface green
    !outblocks: optional parameter (0,1,2).
    !
    !Output:
    !Gn: NE GF (Device + Contacts overlap regions -> effective conductor)
    !   outblocks = 0  no outer parts are computed
    !   outblocks = 1  only D/C part is computed
    !   outblocks = 2  D/C and C/D parts are computed
    !              (needed for K-points calculations)
    !
    !*****************************************************************************

    implicit none

    !In/Out
    type(Tnegf), intent(inout) :: negf
    type(z_DNS), dimension(:), intent(in)  :: SelfEneR, gsurfR, Tlc, Tcl
    real(dp), intent(in)  :: E
    real(dp), dimension(:), intent(in)  :: frm
    type(z_CSR), intent(inout), optional  :: Glout
    integer, intent(in), optional  :: outblocks

    !Work
    integer :: ref, iter
    complex(dp) :: Ec
    integer :: i,ierr,ncont,nbl, lbl, rbl
    integer, dimension(:), allocatable :: Gr_columns
    type(z_CSR) :: ESH_tot, Gl
    logical :: mask(MAXNCONT)


    nbl = negf%str%num_PLs
    ncont = negf%str%num_conts
    ref = negf%refcont
    associate (cblk=>negf%str%cblk, indblk=>negf%str%mat_PL_start)

    Ec = cmplx(E,0.0_dp,dp)

    call build_ESH(negf, Ec, SelfEneR, ncont, nbl, ESH)

    !! Add interaction self energy if any and initialize scba counter
    if (allocated(negf%inter)) then
      call negf%inter%add_sigma_r(ESH)
    end if

    call allocate_gsm(gsmr,nbl)
    call calculate_gsmr_blocks(negf,ESH,nbl,2,gsmr)

    call allocate_blk_dns(Gr,nbl)

    ! compute Gr(1,1)
    call calculate_Gr_tridiag_blocks(negf,ESH,gsmr,Gr,1)
    ! compute Gr(n,n), Gr(n-1,n), Gr(n, n-1);  n = 2 .. nbl
    call calculate_Gr_tridiag_blocks(negf,ESH,gsmr,Gr,2,nbl)

#:if defined("GPU")
      call copy_trid_toHOST(Gr) 
#:endif
    
    !Passing Gr to interaction that builds Sigma_r
    if (allocated(negf%inter)) then
      call negf%inter%set_Gr(Gr, negf%iE)
    end if
    !Computing device G_n
    call allocate_blk_dns(Gn,nbl)
    call init_tridiag_blk(Gn,ESH)

#:if defined("GPU")
      call copy_vdns_toGPU(SelfEneR) 
#:endif
    call calculate_Gn_tridiag_blocks(negf,ESH,SelfEneR,frm,ref,negf%str,gsmr,Gr,Gn)

#:if defined("GPU")
      call copy_trid_toHOST(Gn) 
#:endif

    call destroy_gsm(gsmr) 
    call deallocate_gsm(gsmr)

    !Passing G^n to interaction that builds Sigma^n
    if (allocated(negf%inter)) then
      call negf%inter%set_Gn(Gn, negf%iE)
    end if

    if (present(Glout)) then
      call blk2csr(Gn,negf%str,negf%S,Glout)
    end if
    end associate

    !Computing the 'outer' blocks (device/contact overlapping elements)
    if (present(Glout)) then
      SELECT CASE (outblocks)
      CASE(0)
      CASE(1)
        call calculate_Gn_outer(Tlc,gsurfR,negf%str,frm,ref,.false.,Glout)
      CASE(2)
        call calculate_Gn_outer(Tlc,gsurfR,negf%str,frm,ref,.true.,Glout)
      end SELECT
    end if

    call destroy_all_blk(negf)
    !call destroy_blk(Gn)
    !deallocate(Gn)
    !call destroy_blk(Gr)
    !deallocate(Gr)

#:if defined("GPU")
    call delete_vdns_fromGPU(SelfEneR)
#:endif
    
    !call destroy_tridiag_blk(ESH)
    !deallocate(ESH)

  end subroutine calculate_Gn_neq_components

  !---------------------------------------------------------------------
  !>
  !  Iterative algorithm implementing Meir Wingreen formula for a given
  !  electrode
  !  Note: self consistent born approximation is not accounted for here
  !  It is assumed that the el-ph container already includes the
  !  desired values. SCBA loop should be run outside
  !  Many operations from calls_neq_ph are repeated here, as it is
  !  assumed that A and Gn are not available at the time of the call
  !
  !
  !  I_i = Tr[\Sigma_{i}^{n} A - \Gamma_{i} G^{n}] =
  !      = Tr[\Gamma_{i}( f_{i} A - G^{n} )]
  !
  ! The subroutine assumes that Gr and Gn are available
  ! set negf%tDestroyGr,Gn = .false. when computing Gn
  !---------------------------------------------------------------------
  subroutine iterative_meir_wingreen(negf,E,SelfEneR,frm,curr_mat)
    type(Tnegf), intent(inout) :: negf
    real(dp) :: E
    type(z_DNS), dimension(:) :: SelfEneR
    real(dp), dimension(:) :: frm
    real(dp), dimension(:) :: curr_mat

    !Work
    complex(dp) :: Ec, tmp
    integer :: ii, ncont, lead, lead_blk
    type(z_DNS) :: work1, Gam, A

    ncont = negf%str%num_conts

    do ii = 1, ncont 
      lead = ii 
      lead_blk = negf%str%cblk(lead)
      call zspectral(SelfEneR(lead),SelfEneR(lead), 0, Gam)
      call prealloc_mult(Gam, Gn(lead_blk, lead_blk), work1)
        
      call zspectral(Gr(lead_blk, lead_blk), Gr(lead_blk, lead_blk), 0, A)
      tmp = cmplx(-frm(lead),0.0_dp, dp)
      call prealloc_mult(Gam, A, tmp, work1)
      curr_mat(ii) = -real(trace(work1))

      call destroy(work1)
      call destroy(Gam)
      call destroy(A)
    end do

    call destroy_all_blk(negf)

  end subroutine iterative_meir_wingreen

  !---------------------------------------------------------------------
  !  Calculate the layer current per unit energy
  !
  !    I_LL'(E) = Tr[(ES-H)_LL' * Gn_L'L(E)-(ES-H)_L'L * Gn_LL'(E)]
  !
  ! The subroutine assumes that Gn(ii, ii+1) blocks are available
  ! set negf%tDestroyGn = .false. when computing Gn
  !
  ! Note : ESH(ii,ii+1) are not changed by diagonal Sigma^r
  !        With non-diagonal sigma^r ESH must be recomputed  
  !
  subroutine iterative_layer_current(negf,E,curr_mat)
    type(Tnegf), intent(inout) :: negf
    real(dp) :: E
    real(dp), dimension(:) :: curr_mat

    integer :: nbl, ii
    type(z_DNS) :: work1
    type(z_CSR) :: ESH_tot
    complex(dp), parameter :: minusone = (-1.0_dp, 0.0_dp)
    complex(dp) :: Ec

    nbl = negf%str%num_PLs

    if (size(curr_mat) .ne. nbl-1) then
       stop 'ERROR: curr_mat with wrong size in iterative_layer_current'
    end if
    
    Ec=cmplx(E,0.0_dp,dp)
    call prealloc_sum(negf%H, negf%S, minusone, Ec, ESH_tot)
    call allocate_blk_dns(ESH,nbl)
    call zcsr2blk_sod(ESH_tot, ESH, negf%str%mat_PL_start)
    call destroy(ESH_tot)

    do ii = 1, nbl-1
      call prealloc_mult(ESH(ii,ii+1),Gn(ii+1,ii),work1)
      call prealloc_mult(ESH(ii+1,ii),Gn(ii,ii+1), minusone, work1)
      curr_mat(ii) = real(i_unit*trace(work1)) 
      call destroy(work1)
    end do

    call destroy_all_blk(negf)

  end subroutine iterative_layer_current


  !------------------------------------------------------------------------------!
  ! Transmission_BP_corrected
  !
  ! The routine implements the current-probes, i.e. it assumes that the BPs are
  ! elastic dephasing sources that fullfill current conservation at any energy
  ! For a nice discussion see papers by D. Ryndyk and
  ! M. Kilgour, D. Segal, The J. of Chem Phys 144, 124107 (2016)
  !
  ! LIMITATIONS:
  !
  ! This routine is limited to 2 contacts
  ! It uses dense matrices and assumes there is only 1 PL in the device
  !------------------------------------------------------------------------------!

  subroutine transmission_BP_corrected(negf,SelfEneR,tun_mat)

    type(Tnegf), intent(inout) :: negf
    type(z_DNS), dimension(:) :: SelfEneR
    real(dp), dimension(:) :: tun_mat

    integer :: nn, mm, cont, ni, nf
    integer :: NumOrbs, ncont
    integer, allocatable ::  pls(:), ple(:)
    real(dp), allocatable  :: Trans(:,:), W(:,:), GG(:,:), R(:)
    type(z_DNS) :: Gam1, Gam2, GreenR, GreenA, Tmp1, Tmp2


    NumOrbs=negf%str%central_dim
    ncont = negf%str%num_conts

    allocate(Trans(NumOrbs+ncont,NumOrbs+ncont))
    Trans=0.0_dp

    call create(Gam1, NumOrbs, NumOrbs)
    call create(Gam2, NumOrbs, NumOrbs)

    call blk2dns(Gr,negf%str,GreenR)
    call zdagger(GreenR, GreenA)

    allocate(pls(ncont))
    allocate(ple(ncont))

    do nn = 1, ncont
      pls(nn) = negf%str%mat_PL_start(negf%str%cblk(nn))
      ple(nn) = negf%str%mat_PL_end(negf%str%cblk(nn))
    end do

    do nn = 1, NumOrbs+ncont

      if (nn > NumOrbs) then
         cont = nn - NumOrbs
         call zspectral(SelfEneR(cont),SelfEneR(cont),0,Tmp1)
         Gam1%val(pls(cont):ple(cont),pls(cont):ple(cont)) = Tmp1%val
         call destroy(Tmp1)
      else
         Gam1%val = 0.0_dp
         Gam1%val(nn,nn)=negf%bp_deph%coupling(nn)
      end if

      do mm = 1, NumOrbs+ncont

        if (mm > NumOrbs) then
           cont = mm - NumOrbs
           call zspectral(SelfEneR(cont),SelfEneR(cont),0,Tmp2)
           Gam2%val(pls(cont):ple(cont),pls(cont):ple(cont)) = Tmp2%val
           call destroy(Tmp2)
        else
           Gam2%val = 0.0_dp
           Gam2%val(mm,mm)=negf%bp_deph%coupling(mm)
        end if

        ! Compute coherent transmission: Tr[Gam1 Gr Gam2 Ga]
        ! The inefficient quadruple loop has been substituted with
        ! M*M multiplications exploting OMP parallelism
        call prealloc_mult(Gam1,GreenR,Tmp1)
        call prealloc_mult(Tmp1,Gam2,Tmp2)
        call destroy(Tmp1)
        call prealloc_mult(Tmp2,GreenA,Tmp1)
        call destroy(Tmp2)

        Trans(nn,mm) = real(trace(Tmp1))

        call destroy(Tmp1)

      end do
    end do

    call destroy(Gam1, Gam2, GreenR, GreenA)
    deallocate(pls, ple)

    allocate(W(NumOrbs,NumOrbs))
    allocate(R(NumOrbs+ncont))

    do nn = 1, NumOrbs+ncont
       R(nn) = 1.0_dp
       do mm= 1, NumOrbs+ncont
          if (mm.ne.nn) then
             R(nn) = R(nn)-Trans(mm,nn)
          end if
       end do
    end do

    do nn = 1, NumOrbs
      do mm = 1, NumOrbs
         W(nn,mm) = -Trans(nn,mm)
      end do
      W(nn,nn) = 1.0_dp - R(nn)
    end do

    deallocate(R)
    allocate(GG(NumOrbs,NumOrbs))

    call inverse(GG,W,NumOrbs)

    deallocate(W)

    allocate(R(NumOrbs))

    do nn = 1, size(negf%ni)
      ni = negf%ni(nn)
      nf = negf%nf(nn)
      R =  matmul(Trans(NumOrbs+ni,:),GG)
      tun_mat(nn) = Trans(NumOrbs+ni,NumOrbs+nf) + dot_product(R, Trans(:,NumOrbs+nf))
    end do

    deallocate(Trans,GG)
    deallocate(R)

  end subroutine transmission_BP_corrected

  !------------------------------------------------------------------------------!
  !DAR end
  !------------------------------------------------------------------------------!

  !**********************************************************************
  subroutine init_tridiag_blk(Matrix,S)
    type(z_DNS), dimension(:,:) :: Matrix,S

    integer :: nbl, j

    nbl = SIZE(Matrix,1)

#:if defined("GPU")
    call createAll(Matrix(1,1),S(1,1)%nrow,S(1,1)%ncol)
    Matrix(1,1)%val=(0.0_dp,0.0_dp)
    do j=2,nbl-1
       call createAll(Matrix(j-1,j),S(j-1,j)%nrow,S(j-1,j)%ncol)
       Matrix(j-1,j)%val=(0.0_dp,0.0_dp)
       call createAll(Matrix(j,j),S(j,j)%nrow,S(j,j)%ncol)
       Matrix(j,j)%val=(0.0_dp,0.0_dp)
       call createAll(Matrix(j,j-1),S(j,j-1)%nrow,S(j,j-1)%ncol)
       Matrix(j,j-1)%val=(0.0_dp,0.0_dp)
    end do
    if (nbl.gt.1) then
       call createAll(Matrix(nbl,nbl),S(nbl,nbl)%nrow,S(nbl,nbl)%ncol)
       Matrix(nbl,nbl)%val=(0.0_dp,0.0_dp)
       call createAll(Matrix(nbl-1,nbl),S(nbl-1,nbl)%nrow,S(nbl-1,nbl)%ncol)
       Matrix(nbl-1,nbl)%val=(0.0_dp,0.0_dp)
       call createAll(Matrix(nbl,nbl-1),S(nbl,nbl-1)%nrow,S(nbl,nbl-1)%ncol)
       Matrix(nbl,nbl-1)%val=(0.0_dp,0.0_dp)
    endif

#:else    
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
#:endif
  end subroutine init_tridiag_blk

  !---------------------------------------------------------------------
  subroutine destroy_all_blk(negf)
    type(Tnegf), intent(in) :: negf
    integer :: i

#:if defined("GPU")
    if (negf%tDestroyESH) then
      call destroy_tridiag_blk(ESH)
      if (allocated(ESH)) deallocate(ESH)
    end if  
    if (negf%tDestroyGr) then
      call destroy_blk(Gr)
      if (allocated(Gr)) deallocate(Gr)
    else
      call delete_trid_fromGPU(Gr)
    end if   
    if (negf%tDestroyGn) then
      call destroy_blk(Gn)
      if (allocated(Gn)) deallocate(Gn)
    else
      call delete_trid_fromGPU(Gn)
    end if   
#:else
    if (negf%tDestroyESH) then
      call destroy_tridiag_blk(ESH)
      if (allocated(ESH)) deallocate(ESH)
    end if  
    if (negf%tDestroyGr) then
      call destroy_blk(Gr)
      if (allocated(Gr)) deallocate(Gr)
    end if   
    if (negf%tDestroyGn) then
      call destroy_blk(Gn)
      if (allocated(Gn)) deallocate(Gn)
    end if   
#:endif
  end subroutine destroy_all_blk

  !**********************************************************************
  subroutine destroy_gsm(gsm)
    type(z_DNS), dimension(:) :: gsm
    integer :: i, i1, nbl

    nbl=size(gsm,1)

#:if defined("GPU")
    do i=1,nbl
       if (allocated(gsm(i)%val)) then          
          call destroyAll(gsm(i))
       end if
    end do
#:else
    do i=1,nbl
       if (allocated(gsm(i)%val)) then          
          call destroy(gsm(i))
       end if
    end do
#:endif
  end subroutine destroy_gsm

  !**********************************************************************
  subroutine destroy_blk(M)
    type(z_DNS), dimension(:,:), allocatable :: M
    integer :: i, i1, nbl

    if (.not.allocated(M)) return

    nbl=size(M,1)
#:if defined("GPU")
    do i=1,nbl
       do i1=1,nbl
          if (ALLOCATED(M(i1,i)%val)) THEN
             call destroyAll(M(i1,i))
          end if
       end do
    end do
#:else    
    do i=1,nbl
       do i1=1,nbl
          if (ALLOCATED(M(i1,i)%val)) THEN
             call destroy(M(i1,i))
          end if
       end do
    end do
#:endif
  end subroutine destroy_blk

  !**********************************************************************
  subroutine destroy_tridiag_blk(M)
    type(z_DNS), dimension(:,:), allocatable :: M 

    integer :: i, nbl

    if (.not.allocated(M)) return

    nbl=size(M,1)

    do i=1,nbl
      if (allocated(M(i,i)%val)) then
        call destroy(M(i,i))
      end if  
    end do
    do i=2,nbl
      if (allocated(M(i-1,i)%val)) then
        call destroy(M(i-1,i))
      end if
      if (allocated(M(i,i-1)%val)) then
        call destroy(M(i,i-1))
      end if  
    end do

  end subroutine destroy_tridiag_blk



  !***********************************************************************
  !
  !  g_small left (gsml) calculation - write on memory
  !
  !***********************************************************************

 ! subroutine calculate_gsml_blocks(ESH,sbl,ebl)

 !   !***********************************************************************
 !   !Input:
 !   !ESH: sparse matrices array ESH(nbl,nbl)
 !   !
 !   !nbl (number of layers), indblk(nbl+1) global variables needed
 !   !
 !   !Output:
 !   !sparse matrices array global variable gsml(nbl) is available in memory
 !   !single blocks are allocated internally, array Gr(nbl,nbl)
 !   !must be allocated externally
 !   !***********************************************************************


 !   implicit none

 !   !In/Out
 !   type(z_DNS), dimension(:,:), intent(in) :: ESH
 !   integer, intent(in) :: sbl,ebl                       ! start block, end block

 !   !Work
 !   type(z_DNS) :: work1, work2
 !   integer :: nrow
 !   integer :: i, nbl
 !   !type(z_DNS) :: INV(sbl,sbl)

 !   if (sbl.gt.ebl) return

 !   nbl = size(ESH,1)

 !   if (nbl.eq.1) return

 !   nrow=ESH(sbl,sbl)%nrow

 !   call create(gsml(sbl),nrow,nrow)

 !   call compGreen(gsml(sbl),ESH(sbl,sbl),nrow)


 !   do i=sbl+1,ebl

 !     nrow=ESH(i,i)%nrow

 !     call prealloc_mult(ESH(i,i-1),gsml(i-1),(-1.0_dp, 0.0_dp),work1)

 !     call prealloc_mult(work1,ESH(i-1,i),work2)

 !     call destroy(work1)

 !     call prealloc_sum(ESH(i,i),work2,work1)

 !     call destroy(work2)

 !     call create(gsml(i),work1%nrow,work1%nrow)

 !     call compGreen(gsml(i),work1,work1%nrow)

 !     call destroy(work1)

 !   end do

 !   if (debug) then
 !     WRITE(*,*) '********************'
 !     WRITE(*,*) 'calculate_gsml done'
 !     WRITE(*,*) '********************'
 !   endif

 ! end subroutine calculate_gsml_blocks


  !**************************************************************************
  !
  !  Calculate Green Retarded column "n" - writing on memory
  !
  !**************************************************************************
 ! subroutine calculate_Gr_column_blocks(ESH,n,indblk)

 !   !***********************************************************************
 !   !Input:
 !   !ESH: sparse matrices array ESH(nbl,nbl)
 !   !n: n umber of column to be calculated
 !   !
 !   !global variables needed: nbl (number of layers), indblk(nbl+1),
 !   !Gr diagonal, subadiagonal and superdiagonal, gsmr(:) for
 !   !downgoing and gsml(:) for upgoing
 !   !
 !   !Output:
 !   !sparse matrices array global variable Gr(:,n) is available in
 !   !memory - single blocks are allocated internally, array Gr(nbl,nbl)
 !   !must be allocated externally
 !   !***********************************************************************

 !   implicit none

 !   !In/Out
 !   type(z_DNS), dimension(:,:), intent(in) :: ESH
 !   integer, intent(in) :: n
 !   integer, dimension(:), intent(in) :: indblk

 !   !Work
 !   integer :: i,nrow,ncol,nbl
 !   type(z_DNS) :: work1
 !   real(dp) :: max

 !   nbl = size(ESH,1)

 !   if (n.GT.nbl) THEN
 !     STOP 'Error in calculate_Grcol : n is greater than nbl'
 !   endif

 !   !***************************************
 !   !  Downgoing (j>=n+2 && n<nbl-1)
 !   !
 !   !   G_j,n = -gR_jj T_j,j-1 G_j-1,n
 !   !
 !   !***************************************
 !   if (n.LT.(nbl-1)) THEN

 !     do i=n+2,nbl

 !       max=MAXVAL(ABS(Gr(i-1,n)%val))
 !       if (max.GT.EPS) THEN
 !         call prealloc_mult(gsmr(i),ESH(i,i-1),(-1.0_dp, 0.0_dp),work1)
 !         call prealloc_mult(work1,Gr(i-1,n),Gr(i,n))
 !         call destroy(work1)
 !       else
 !         ! WHEN BLOCK IS SMALLER THAN EPS IT IS NOT CREATED
 !         exit
 !       end if

 !     end do

 !   endif
 !   !*************************************
 !   !   Up-going (j<=n-2 && n>2)
 !   !
 !   !   G_j,n = -gL_jj T_j,j+1 G_j+1,n
 !   !
 !   !*************************************

 !   if (n.GT.2) THEN

 !     do i=n-2,1,(-1)

 !       max=MAXVAL(ABS(Gr(i+1,n)%val))

 !       if (max.GT.EPS) THEN
 !         call prealloc_mult(gsml(i),ESH(i,i+1),(-1.0_dp, 0.0_dp),work1)
 !         call prealloc_mult(work1,Gr(i+1,n),Gr(i,n))
 !         call destroy(work1)
 !       else
 !         ! WHEN BLOCK IS SMALLER THAN EPS IT IS NOT CREATED
 !         exit
 !       endif

 !     end do

 !   endif

 !   if (debug) then
 !     WRITE(*,*) '******************************'
 !     WRITE(*,*) 'calculate_Grcol done column',n
 !     WRITE(*,*) '******************************'
 !   endif

 ! end subroutine calculate_Gr_column_blocks

  !****************************************************************************
  !
  !  Calculate Green Retarded - writing on memory (optimized on mask)
  !
  !****************************************************************************

  subroutine Gr_blk2csr(P,nbl,indblk,A)

    !****************************************************************************
    !Input:
    !P: CSR matrix containing masking pattern
    !
    !global variable needed: nbl, indblk(nbl+1), Gr(:,:)
    !
    !Output:
    !A: sparse matrix containing Green Retarded of device (allocated internally)
    !****************************************************************************

    implicit none

    !In/Out
    integer :: nbl
    integer, dimension(:), pointer :: indblk
    type(z_CSR) :: A, P, GrCsr

    !Work
    integer :: i, j, i1, ix, iy, x, y, col, oldx

    !create A with same pattern of P
    call create(A,P%nrow,P%ncol,P%nnz)
    A%rowpnt(:)=P%rowpnt(:)
    A%colind(:)=P%colind(:)
    A%nzval = (0.0_dp,0.0_dp)

    !If only one block is present, concatenation is not needed
    !and it's implemented in a more trivial way
    if (nbl.EQ.1) THEN

      call create(GrCsr,Gr(1,1)%nrow,Gr(1,1)%ncol,Gr(1,1)%nrow*Gr(1,1)%ncol)
      call dns2csr(Gr(1,1),GrCsr)

      call mask(GrCsr,P,A)
      call destroy(GrCsr)

    ELSE

      !Cycle upon all rows
      x = 1
      do i = 1, A%nrow
        !Choose which block (row) we're dealing with
        oldx = x

        !Check if row is in same block of previous or in next block. Not needed
        !(and not allowed not to exceed indblk index boundaries) if we're in the last block
        if (oldx.EQ.nbl) THEN
          x = oldx
        ELSE
          do ix = oldx, oldx+1
            if ( (i.GE.indblk(ix)).AND.(i.LT.indblk(ix+1)) ) x = ix
          end do
        endif

        !Offset: i1 is the index for separate blocks
        i1 = i - indblk(x) + 1
        !Cycle upon columns
        do j = A%rowpnt(i), A%rowpnt(i+1) -1
          !Choose which block column we're dealing with
          y = 0
          if (x.EQ.1) THEN
            if ( (A%colind(j).GE.indblk(x)).AND.(A%colind(j).LT.indblk(x + 1)) ) then
              y = 1
            ELSEif ( (A%colind(j).GE.indblk(x + 1)).AND.(A%colind(j).LT.indblk(x + 2)) ) then
              y = 2
            endif
          elseif (x.eq.nbl) then
            if ( (A%colind(j).GE.indblk(x)).AND.(A%colind(j).LT.indblk(x + 1)) ) then
              y = nbl
            ELSEif ( (A%colind(j).GE.indblk(x - 1)).AND.(A%colind(j).LT.indblk(x)) ) then
              y = nbl - 1
            endif
          ELSE
            do iy = x-1, x+1
              if ( (A%colind(j).GE.indblk(iy)).AND.(A%colind(j).LT.indblk(iy + 1)) ) y = iy
            end do
          endif
          if (y.eq.0) then
            write(*,*)
            write(*,*) 'ERROR in Gr_blk2csr: probably wrong PL size',x
            write(*,*) 'row',i,A%colind(j)
            write(*,*) 'block indeces:',indblk(1:nbl)
            stop
          endif

          col = A%colind(j) - indblk(y) + 1

          A%nzval(j) = Gr(x,y)%val(i1,col)

        end do

      end do

    endif

    !if (debug) call writePeakInfo(6)
    if (debug) then
      WRITE(*,*) '**********************'
      WRITE(*,*) 'calculate_GreenR done'
      WRITE(*,*) '**********************'
    endif

  end subroutine Gr_blk2csr

  ! blk-sparse to dense converision
  subroutine blk2dns(G,str,Gdns)
    type(z_DNS), dimension(:,:), intent(in) :: G
    type(Tstruct_info), intent(in) :: str
    type(z_DNS) :: Gdns

    integer :: n, m, ii, jj, pl_start1, pl_end1, pl_start2, pl_end2, nbl, nrows

    nbl = str%num_PLs
    nrows = str%central_dim
    call create(Gdns,nrows,nrows)

    do n = 1, nbl
      do m= 1, nbl
         pl_start1 = str%mat_PL_start(n)
         pl_end1 = str%mat_PL_end(n)
         pl_start2 = str%mat_PL_start(m)
         pl_end2 = str%mat_PL_end(m)

         do ii = 1, pl_end1 - pl_start1 + 1
           do jj = 1, pl_end2 - pl_start2 + 1
               Gdns%val(pl_start1 + ii - 1,pl_start2 + jj - 1) = G(n,m)%val(ii,jj)
           end do
         end do
      end do
    end do

  end subroutine blk2dns


  !Concatenation for every contact in G_n. Performs a sum on elements, not a replacement
  !Similar to calculate_GreenR2, except for sum of elements
  !Note: to backup old version zconcat calls (and Glsub deallocations) must be
  !      uncommented and all this part removed
  !If only one block is present, concatenation is not needed and it's implemented in a
  !more trivial way
  subroutine blk2csr(G,struct,P,Gcsr)

    type(z_DNS), dimension(:,:) :: G
    type(Tstruct_info), intent(in) :: struct
    type(z_CSR) :: Gcsr
    type(z_CSR) :: P, G_sp

    integer :: nbl, oldx, row, col, iy, ix, x, y, ii, jj, nrows

    nbl = struct%num_PLs
    nrows = struct%mat_PL_end(nbl)

    !create Gcsr with same pattern of P
    call create(Gcsr,P%nrow,P%ncol,P%nnz)
    Gcsr%rowpnt = P%rowpnt
    Gcsr%colind = P%colind
    Gcsr%nzval = (0.0_dp, 0.0_dp)

    associate(indblk=>struct%mat_PL_start)
    !Cycle upon all rows
    x = 1
    do ii = 1, nrows
      !Search block x containing row ii
      oldx = x
      if (oldx.EQ.nbl) THEN
        x = oldx
      ELSE
        do ix = oldx, oldx+1
          if ( (ii.GE.indblk(ix)).AND.(ii.LT.indblk(ix+1)) ) x = ix
        end do
      endif

      !Offset: row is the index for separate blocks
      row = ii - indblk(x) + 1

      !Cycle upon columns of Gcsr (which has been ALREADY MASKED by S)
      do jj = Gcsr%rowpnt(ii), Gcsr%rowpnt(ii+1) -1
        if (Gcsr%colind(jj).gt.nrows) CYCLE
        !Choose which block column we're dealing with
        y = 0
        if (x.eq.1) then
          if ( (Gcsr%colind(jj).GE.indblk(x)).AND.(Gcsr%colind(jj).LT.indblk(x + 1)) ) then
            y = 1
          ELSEif ( (Gcsr%colind(jj).GE.indblk(x + 1)).AND.(Gcsr%colind(jj).LT.indblk(x + 2)) ) then
            y = 2
          endif
        elseif (x.eq.nbl) then
          if ( (Gcsr%colind(jj).GE.indblk(x)).AND.(Gcsr%colind(jj).LT.indblk(x + 1)) ) then
            y = nbl
          ELSEif ( (Gcsr%colind(jj).GE.indblk(x - 1)).AND.(Gcsr%colind(jj).LT.indblk(x)) ) then
            y = nbl - 1
          endif
        else
          do iy = x-1, x+1
            if ( (Gcsr%colind(jj).GE.indblk(iy)).AND.(Gcsr%colind(jj).LT.indblk(iy + 1)) ) y = iy
          end do
        endif

        if (y.EQ.0) THEN
          write(*,*)
          write(*,*) 'ERROR in blk2csr: probably wrong PL size', x
          write(*,*) 'row',ii,nrows,Gcsr%colind(jj)
          write(*,*) 'block indeces:',indblk(1:nbl)
          STOP
        endif

        col = Gcsr%colind(jj) - indblk(y) + 1

        if (allocated(G(x,y)%val)) THEN
          Gcsr%nzval(jj) = Gcsr%nzval(jj) + G(x,y)%val(row,col)
        endif

      end do

    end do
    end associate

  end subroutine blk2csr


  !****************************************************************************
  !
  !  Calculate Green Retarded in the
  !  contacts regions, where overlap with device orbitals is non-zero
  !  (only the upper part, needed in both K=0 and K points calculations,
  !   or both upper and lower parts)
  !  writing on memory
  !
  !****************************************************************************

  subroutine calculate_Gr_outer(Tlc,Tcl,gsurfR,struct,lower,Aout)

    !****************************************************************************
    !Input:
    !Tlc: sparse arraymatrices array containing contact-device interacting blocks
    !     ESH-H
    !Tlc: sparse arraymatrices array containing device-contact interacting blocks
    !     ESH-H
    !gsurfR: sparse matrices array containing contacts surface green
    !lower: if .true., also lower parts are calculated and concatenated
    !
    !global variables needed: nbl, indblk(nbl+1), ncont, cblk(ncont),
    !cindblk(ncont), Gr in the diagonal block related to the interaction layer
    !
    !Output:
    !Aout: sparse matrix containing density matrix in the region
    !      corresponding to non-zero overlap
    !
    !****************************************************************************

    implicit none

    !In/Out
    type(z_DNS), dimension(:) :: Tlc,Tcl,gsurfR
    logical :: lower
    type(Tstruct_info), intent(in) :: struct
    type(z_CSR) :: Aout


    !Work
    type(z_DNS) :: work1, Grcl, Grlc
    type(z_CSR) :: GrCSR, TCSR
    integer :: i,cb,nrow_tot,i1,j1
    integer :: ncont, nbl

    ncont = struct%num_conts
    nbl = struct%num_PLs
    nrow_tot = struct%total_dim

    if (.not.allocated(Aout%nzval)) THEN
      call create(Aout,nrow_tot,nrow_tot,0)
      Aout%rowpnt(:)=1
    endif

    do i=1,ncont

      !Numero di blocco del contatto
      cb=struct%cblk(i)
      call prealloc_mult(Gr(cb,cb),Tlc(i),(-1.0_dp, 0.0_dp),work1)
      call prealloc_mult(work1,gsurfR(i),Grlc)

      call destroy(work1)

      j1=nzdrop(Grlc,EPS)
      call create(GrCSR,Grlc%nrow,Grlc%ncol,j1)
      call dns2csr(Grlc,GrCSR)
      call destroy(Grlc)
      j1=nzdrop(Tlc(i),EPS)
      call create(TCSR,Tlc(i)%nrow,Tlc(i)%ncol,j1)
      call dns2csr(Tlc(i),TCSR)
      call zmask_realloc(GrCSR,TCSR)
      call destroy(TCSR)

      !Concatenazione di Asub nella posizione corrispondente
      i1=struct%mat_PL_start(cb)
      j1=struct%mat_B_start(i)

      call concat(Aout,GrCSR,i1,j1)

      call destroy(GrCSR)

      if (lower) THEN

        call prealloc_mult(gsurfR(i),Tcl(i),(-1.0_dp, 0.0_dp), work1)
        call prealloc_mult(work1, Gr(cb,cb), Grcl)

        call destroy(work1)

        j1=nzdrop(Grcl,EPS)
        call create(GrCSR,Grcl%nrow,Grcl%ncol,j1)
        call dns2csr(Grcl,GrCSR)
        call destroy(Grcl)
        j1=nzdrop(Tcl(i),EPS)
        call create(TCSR,Tcl(i)%nrow,Tcl(i)%ncol,j1)
        call dns2csr(Tcl(i),TCSR)
        call zmask_realloc(GrCSR,TCSR)
        call destroy(TCSR)

        i1 = struct%mat_B_start(i)-struct%central_dim+struct%mat_PL_start(nbl+1)-1
        j1 = struct%mat_PL_start(cb)

        call concat(Aout,GrCSR,i1,j1)

        call destroy(GrCSR)

      endif

    end do

    if (debug) then
      WRITE(*,*) '********************'
      WRITE(*,*) 'Outer_GreenR done'
      WRITE(*,*) '********************'
    endif

  end subroutine calculate_Gr_outer


  subroutine calculate_Gn_outer(Tlc,gsurfR,struct,frm,ref,lower,Gn_out)
    type(z_DNS), dimension(:), intent(in) :: Tlc, gsurfR
    real(dp), dimension(:), intent(in) :: frm
    type(Tstruct_info), intent(in) :: struct
    integer, intent(in) :: ref
    logical, intent(in) :: lower
    type(z_CSR) :: Gn_out

    !Work
    type(z_DNS) :: gsurfA, work1, work2, work3, Gn_lc
    type(z_CSR) :: GnCSR, TCSR
    integer :: k,cbk,i1,j1,nrow_tot
    integer :: ncont, nbl
    complex(dp) :: frmdiff

    ncont = struct%num_conts
    nbl = struct%num_PLs

    if (.not.allocated(Gn_out%nzval)) THEN
      nrow_tot = struct%total_dim
      call create(Gn_out,nrow_tot,nrow_tot,0)
      Gn_out%rowpnt(:)=1
    endif

    do k=1,ncont
      ! Sigma_n(cbk,ck) = Sigma_r(cbk,ck) = 0  => Tlc^r = Tlc^a = Tlc = T(cbk,ck)
      ! Gn(cbk,ck) =  -Gn(cbk,cbk)*Tlc*gA(cb) -(fk-fr)*Gr(cbk,cbk)*Tlc*j(gsurfR-gsurfA)]

      cbk=struct%cblk(k)
      call zdagger(gsurfR(k),gsurfA)
      call prealloc_mult(Tlc(k),gsurfA,work2)
      call destroy(gsurfA)

      call prealloc_mult(Gn(cbk,cbk),work2,(-1.0_dp,0.0_dp),work3)
      call destroy(work2)

      !Checks that Fermi levels are sufficiently different and contact is not reference
      if ((ABS(frm(k)-frm(ref)).GT.EPS).AND.(k.NE.ref)) THEN

        frmdiff = cmplx(frm(ref)-frm(k),0.0_dp,dp)
        call zspectral(gsurfR(k),gsurfR(k),0,work1)

        call prealloc_mult(Tlc(k),work1,work2)
        call destroy(work1)

        call prealloc_mult(Gr(cbk,cbk),work2,frmdiff,work1)
        call destroy(work2)
      else
        call create(work1, Gr(cbk,cbk)%nrow, Tlc(k)%ncol)
        work1%val=(0.0_dp, 0.0_dp)
      end if

      call prealloc_sum(work3,work1,Gn_lc)
      call destroy(work1)
      call destroy(work3)

      call mask(Gn_lc,Tlc(k))
      i1=nzdrop(Gn_lc,EPS)

      if (i1.gt.0) THEN
        call create(GnCSR,Gn_lc%nrow,Gn_lc%ncol,i1)
        call dns2csr(Gn_lc,GnCSR)
        call destroy(Gn_lc)

        ! GnCSR is concatenated (added) to Gn_out
        i1=struct%mat_PL_start(cbk)
        j1=struct%mat_B_start(k)-struct%central_dim+struct%mat_PL_start(nbl+1)-1
        call concat(Gn_out,GnCSR,i1,j1)

        ! lower Gn(c,cb) outer part is computed via Gn(c,cb) = Gn(cb,c)+
        if (lower) THEN
          call zdagger(GnCSR,TCSR)
          call concat(Gn_out,TCSR,j1,i1)
          call destroy(TCSR)
        endif

        call destroy(GnCSR)
      else
        call destroy(Gn_lc)
      end if

    end do


  end subroutine calculate_Gn_outer

  !---------------------------------------------------

  subroutine calculate_transmissions(negf,H,S,Ec,SelfEneR,ni,nf,str,tun_proj,tun_mat)
    Type(Tnegf), intent(in) :: negf
    Type(z_CSR) :: H
    Type(z_CSR) :: S
    Complex(dp) :: Ec
    Type(z_DNS), Dimension(MAXNCONT) :: SelfEneR
    Integer :: ni(:)
    Integer :: nf(:)
    Type(z_CSR) :: ESH_tot
    Type(TStruct_Info) :: str
    type(intarray), intent(in) :: tun_proj
    Real(dp), Dimension(:) :: tun_mat

    ! Local variables
    Real(dp) :: tun
    Integer :: nbl,ncont,ibl
    Integer :: i, ierr, icpl, nit, nft, nt, nt1

    nbl = str%num_PLs
    ncont = str%num_conts

    if (ncont == 1) then
      tun_mat = 0.0_dp
      return
    end if

    !With this method H and S have to be passed through negf. negf is also necessary for GPU handlers
    call build_ESH(negf, Ec, SelfEneR, ncont, nbl, ESH)

    nit=ni(1)
    nft=nf(1)
    ! find the contact with smaller block index
    if (str%cblk(nit).lt.str%cblk(nft)) then
      nt = str%cblk(nit)
    else
      nt = str%cblk(nft)
    endif

    ! Fall here when there are 2 contacts for fast transmission
    if (ncont == 2 .and. size(ni) == 1 .and. nt == 1) then
      call allocate_gsm(gsmr,nbl)
      call calculate_gsmr_blocks(negf,ESH,nbl,2,gsmr,.false.)
      call allocate_blk_dns(Gr,nbl)
      call calculate_Gr_tridiag_blocks(negf,ESH,gsmr,Gr,1)
#:if defined("GPU")
       call copy_vdns_toGPU(SelfEneR)
#:endif
      call calculate_single_transmission_2_contacts(negf,nit,nft,ESH,SelfEneR,str%cblk,tun_proj,Gr,tun)
      tun_mat(1) = tun

    else
      ! MULTITERMINAL case
      call allocate_gsm(gsmr,nbl)
      call calculate_gsmr_blocks(negf,ESH,nbl,2,gsmr)

#:if defined("GPU")
       call copy_vdns_toGPU(SelfEneR)
#:endif
      do icpl = 1, size(ni)

        !Computation of transmission(s) between contacts ni(:) -> nf(:)
        nit=ni(icpl)
        nft=nf(icpl)

        ! find the largest contact block between the two terminals
        if (str%cblk(nit).gt.str%cblk(nft)) then
          nt1 = str%cblk(nit)
        else
          nt1 = str%cblk(nft)
        endif

        if (icpl == 1) then
          ! Iterative calculation of Gr down to nt
          nt = nt1
          call allocate_blk_dns(Gr,nbl)
          call calculate_Gr_tridiag_blocks(negf,ESH,gsmr,Gr,1)
          if (nt.gt.1) then
            call calculate_Gr_tridiag_blocks(negf,ESH,gsmr,Gr,2,nt)
          end if
        else
          ! When more contacts are present sometimes we can re-use previous GF
          ! if nt1 > nt extend the Gr calculation
          if (nt1 .gt. nt) then
            call calculate_Gr_tridiag_blocks(negf,ESH,gsmr,Gr,nt+1,nt1)
            nt = nt1
          endif
        end if

        call calculate_single_transmission_N_contacts(negf,nit,nft,ESH,SelfEneR,str%cblk,tun_proj,gsmr,Gr,tun)

        tun_mat(icpl) = tun

      end do
    end if

#:if defined("GPU")
    call delete_vdns_fromGPU(SelfEneR)
#:endif
    !Deallocate energy-dependent matrices
    call destroy_gsm(gsmr)
    call deallocate_gsm(gsmr)

    call destroy_blk(Gr)
    call deallocate_blk_dns(Gr)

    call destroy_tridiag_blk(ESH)
    call deallocate_blk_dns(ESH)

  end subroutine calculate_transmissions

  !---------------------------------------------------!
  !Subroutine for transmission and dos calculation    !
  !---------------------------------------------------!

  subroutine calculate_transmissions_and_dos(negf,H,S,Ec,SelfEneR,Gs,ni,nf,str,tun_proj,TUN_MAT,dos_proj,LEDOS)
    Type(Tnegf), intent(in) :: negf
    Type(z_CSR), intent(in) :: H
    Type(z_CSR), intent(in) :: S
    Complex(dp), intent(in) :: Ec
    Type(z_DNS), Dimension(MAXNCONT), intent(in) :: SelfEneR, Gs
    Integer, intent(in) :: ni(:)
    Integer, intent(in) :: nf(:)
    Type(TStruct_Info), intent(in) :: str
    type(intarray), intent(in) :: tun_proj
    type(intarray), dimension(:), intent(in) :: dos_proj
    Real(dp), Dimension(:), intent(inout) :: TUN_MAT
    Real(dp), Dimension(:), intent(inout) :: LEDOS

    ! Local variables
    Type(z_CSR) :: ESH_tot, GrCSR
    Type(z_DNS), Dimension(:,:), allocatable :: ESH
    Type(r_CSR) :: Grm                          ! Green Retarded nella molecola
    real(dp), dimension(:), allocatable :: diag
    Real(dp) :: tun
    Complex(dp) :: zc
    Integer :: nbl,ncont, ierr
    Integer :: nit, nft, icpl
    Integer :: iLDOS, i2, i
    Character(1) :: Im


    nbl = str%num_PLs
    ncont = str%num_conts
    Im = 'I'

    call build_ESH(negf, Ec, SelfEneR, ncont, nbl, ESH)

    call allocate_gsm(gsmr,nbl)
    call calculate_gsmr_blocks(negf,ESH,nbl,2,gsmr)

    call allocate_blk_dns(Gr,nbl)
    ! call create
    call calculate_Gr_tridiag_blocks(negf,ESH,gsmr,Gr,1)
    call calculate_Gr_tridiag_blocks(negf,ESH,gsmr,Gr,2,nbl)
    !Computation of transmission(s) between contacts ni(:) -> nf(:)
#:if defined("GPU")
    call copy_vdns_toGPU(SelfEneR)
#:endif
    do icpl=1,size(ni)

      nit=ni(icpl)
      nft=nf(icpl)

      select case(ncont)
      case(1)
        tun = 0.0_dp
      case(2)
        call calculate_single_transmission_2_contacts(negf,nit,nft,ESH,SelfEneR,str%cblk,tun_proj,Gr,tun)
      case default
        call calculate_single_transmission_N_contacts(negf,nit,nft,ESH,SelfEneR,str%cblk,tun_proj,gsmr,Gr,tun)
      end select

      TUN_MAT(icpl) = tun

    end do

    !Deallocate energy-dependent matrices
    call destroy_gsm(gsmr)
    call deallocate_gsm(gsmr)
    call destroy_tridiag_blk(ESH)
    call deallocate_blk_dns(ESH)

    ! Destroy only off-diagonal blocks
#:if defined("GPU")    
    call copy_trid_toHOST(Gr) 
    do i=2,nbl
      call destroyAll(Gr(i-1,i))
      call destroyAll(Gr(i,i-1))
    end do
#:else
    do i=2,nbl
       call destroy(Gr(i-1,i))
       call destroy(Gr(i,i-1))
    end do
#:endif

    call create(Grm,str%mat_PL_start(nbl+1)-1,str%mat_PL_start(nbl+1)-1,0)

    Grm%rowpnt(:)=1

#:if defined("GPU")
    call delete_vdns_fromGPU(SelfEneR)
    do i=1,nbl
       call create(GrCSR,Gr(i,i)%nrow,Gr(i,i)%ncol,Gr(i,i)%nrow*Gr(i,i)%ncol)
       call dns2csr(Gr(i,i),GrCSR)
       !Concatena direttamente la parte immaginaria per il calcolo della doS
       zc=(-1.0_dp,0.0_dp)/pi

       call concat(Grm,zc,GrCSR,Im,str%mat_PL_start(i),str%mat_PL_start(i))
       call destroyAll(Gr(i,i))
       call destroy(GrCSR)
    end do
#:else
    do i=1,nbl
       call create(GrCSR,Gr(i,i)%nrow,Gr(i,i)%ncol,Gr(i,i)%nrow*Gr(i,i)%ncol)
       call dns2csr(Gr(i,i),GrCSR)
       !Concatena direttamente la parte immaginaria per il calcolo della doS
       zc=(-1.0_dp,0.0_dp)/pi

       call concat(Grm,zc,GrCSR,Im,str%mat_PL_start(i),str%mat_PL_start(i))
       call destroy(Gr(i,i))
       call destroy(GrCSR)
    end do
#:endif
    call deallocate_blk_dns(Gr)

    !Compute LDOS on the specified intervals
    if (size(dos_proj).gt.0) then
      call log_allocate(diag, Grm%nrow)
      call getdiag(Grm,diag)
      if (size(dos_proj) == size(diag) .and. dos_proj(1)%indexes(1) == 0) then
          LEDOS(:) = LEDOS(:) + diag(:)
       else
        do iLDOS=1,size(dos_proj)
          do i = 1, size(dos_proj(iLDOS)%indexes)
            i2 = dos_proj(iLDOS)%indexes(i)
            if (i2 .le. str%central_dim) then
              LEDOS(iLDOS) = LEDOS(iLDOS) + diag(i2)
            end if
          end do
        end do
      endif
      call log_deallocate(diag)
    endif

    call destroy(Grm)

  end subroutine calculate_transmissions_and_dos

  !---------------------------------------------------

  subroutine build_ESH(negf, Ec, SelfEneR, ncont, nbl, ESH)
    type(Tnegf), intent(in) :: negf
    complex(dp), intent(in) :: Ec      
    type(z_DNS), dimension(:), intent(in)  :: SelfEneR
    integer, intent(in) :: ncont
    integer, intent(in) :: nbl
    type(z_DNS), dimension(:,:), allocatable :: ESH

    Type(z_CSR) :: ESH_tot
    integer :: i

    associate (cblk=>negf%str%cblk)
    ! Take CSR H,S and build ES-H in dense blocks
      call prealloc_sum(negf%H,negf%S,(-1.0_dp, 0.0_dp),Ec,ESH_tot)

      call allocate_blk_dns(ESH, nbl)

      call zcsr2blk_sod(ESH_tot,ESH, negf%str%mat_PL_start)

      call destroy(ESH_tot)

      do i=1,ncont
         ESH(cblk(i),cblk(i))%val = ESH(cblk(i),cblk(i))%val-SelfEneR(i)%val
      end do
     end associate

   end subroutine build_ESH
  !---------------------------------------------------


  subroutine allocate_gsm(gsm,nbl)
    type(z_DNS), dimension(:), allocatable :: gsm
    integer :: nbl, ierr

    if (.not.allocated(gsm)) then
      allocate(gsm(nbl),stat=ierr)
    end if   
    if (ierr.ne.0) stop 'ALLOCATION ERROR: could not allocate gsm'

  end subroutine allocate_gsm

  !---------------------------------------------------

  subroutine allocate_blk_dns(blkM,nbl)
    type(z_DNS), dimension(:,:), allocatable :: blkM
    integer :: nbl, ierr

    allocate(blkM(nbl,nbl),stat=ierr)
    if (ierr.ne.0) stop 'ALLOCATION ERROR: could not allocate block-Matrix'

  end subroutine allocate_blk_dns

  !---------------------------------------------------

  subroutine deallocate_gsm(gsm)
    type(z_DNS), dimension(:), allocatable :: gsm
    integer :: ierr

    deallocate(gsm,stat=ierr)
    if (ierr.ne.0) stop 'DEALLOCATION ERROR: could not deallocate gsmr'

  end subroutine deallocate_gsm

  !---------------------------------------------------

  subroutine deallocate_blk_dns(blkM)
    type(z_DNS), dimension(:,:), allocatable :: blkM
    integer :: ierr

    deallocate(blkM,stat=ierr)
    if (ierr.ne.0) stop 'DEALLOCATION ERROR: could not deallocate block-Matrix'

  end subroutine deallocate_blk_dns

end module iterative


  !****************************************************************************
  !
  ! Calculate G_p=iG> contributions due to el-ph
  ! Writing on memory
  !
  !****************************************************************************
!!$  subroutine calculate_Gp_ph(negf,ESH,iter,Gp)
!!$
!!$    type(Tnegf) :: negf
!!$    type(z_DNS), dimension(:,:) :: ESH, Gp
!!$    integer :: iter
!!$
!!$    Type(z_DNS), dimension(:,:), allocatable :: Sigma_ph_p
!!$    Type(z_DNS) :: Ga, work1, work2
!!$    integer :: n, k, nbl, nrow, ierr
!!$
!!$    nbl = negf%str%num_PLs
!!$    ALLOCATE(Sigma_ph_p(nbl,nbl),stat=ierr)
!!$    if (ierr.NE.0) STOP 'ALLOCATION ERROR: could not allocate Sigma_ph_n'
!!$
!!$
!!$    do n = 1, nbl
!!$
!!$      nrow = ESH(n,n)%nrow
!!$
!!$      call create(Sigma_ph_p(n,n), nrow, nrow)
!!$
!!$      Sigma_ph_p(n,n)%val = (0.0_dp, 0.0_dp)
!!$      if (iter .gt. 0) then
!!$         call read_blkmat(Sigma_ph_p(n,n),negf%scratch_path,'Sigma_ph_p_',n,n,negf%iE)
!!$      else
!!$         call write_blkmat(Sigma_ph_p(n,n),negf%scratch_path,'Sigma_ph_p_',n,n,negf%iE)
!!$      endif
!!$
!!$    end do
!!$
!!$    do n = 1, nbl
!!$
!!$      do k = 1, nbl
!!$
!!$         if (Gr(n,k)%nrow.gt.0) then
!!$            call zdagger(Gr(n,k),Ga)
!!$            call prealloc_mult(Gr(n,k), Sigma_ph_p(k,k), work1)
!!$            call prealloc_mult(work1, Ga, work2)
!!$            Gp(n,n)%val = Gp(n,n)%val + work2%val
!!$            call destroy(work1, work2, Ga)
!!$         endif
!!$
!!$      end do
!!$
!!$    end do
!!$
!!$    do n = 1, nbl
!!$      call destroy(Sigma_ph_p(n,n))
!!$    end do
!!$
!!$    DEALLOCATE(Sigma_ph_p)
!!$
!!$  end subroutine calculate_Gp_ph
