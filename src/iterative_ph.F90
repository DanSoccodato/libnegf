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


MODULE iterative_ph

  USE ln_precision
  USE ln_constants, only : pi, zero, minusone
  USE ln_allocation
  USE mat_def
  USE sparsekit_drv
  USE inversions
  USE ln_structure, only : TStruct_Info
  USE lib_param, only : MAXNCONT, Tnegf, intarray
  USE outmatrix, only : outmat_c, inmat_c, direct_out_c, direct_in_c 
  USE clock
  !USE transform

  IMPLICIT NONE
  private

  public :: calls_Dr_ph
  public :: calls_Dn_ph
  public :: create_scratch
  public :: destroy_scratch

  LOGICAL, PARAMETER :: debug=.false. 

  TYPE(z_DNS), DIMENSION(:), ALLOCATABLE :: gsmr
  TYPE(z_DNS), DIMENSION(:), ALLOCATABLE :: gsml
  TYPE(z_DNS), DIMENSION(:,:), ALLOCATABLE :: Gr

  TYPE(z_DNS3), DIMENSION(:,:), ALLOCATABLE, SAVE :: DDn
  TYPE(z_DNS3), DIMENSION(:,:), ALLOCATABLE, SAVE :: DDp
  TYPE(z_DNS3), DIMENSION(:,:), ALLOCATABLE, SAVE :: DDr
  TYPE(z_DNS3), DIMENSION(:,:), ALLOCATABLE, SAVE :: Pin
  TYPE(z_DNS3), DIMENSION(:,:), ALLOCATABLE, SAVE :: Pip
  TYPE(z_DNS3), DIMENSION(:,:), ALLOCATABLE, SAVE :: Pir
  LOGICAL, PARAMETER :: memory = .true.

CONTAINS

  !****************************************************************************
  ! 
  ! Driver for computing Equilibrium Green Retarded + Sigma_c + Sigma_ph 
  !
  !****************************************************************************

  SUBROUTINE calls_Dr_ph(pnegf,E,SelfEneR,struct,comp_col,A)

    !****************************************************************************
    !
    !Input
    !pnegf:    negf data container
    !E:        Energy point
    !SelfEneR: matrices array containing contacts Self Energy
    !struct:   structure container to get nPls, nConts, indblk
    !comp_col: compute column GF
    !A:        optional, gives back the whole Gr masked by the overlap  
    !
    !*****************************************************************************

    IMPLICIT NONE
 
    !In/Out
    TYPE(Tnegf), intent(in) :: pnegf
    COMPLEX(dp), intent(in) :: E
    TYPE(z_DNS), DIMENSION(:), intent(in) :: SelfEneR
    TYPE(Tstruct_info), intent(in) :: struct
    logical, intent(in) :: comp_col 
    TYPE(z_CSR), intent(out), optional :: A

    !Work
    TYPE(z_DNS), DIMENSION(:,:), ALLOCATABLE :: ESH
    TYPE(z_CSR) :: ESH_csr 
    INTEGER :: i, m, cb, ierr, nbl, rbl, lbl, ncont, iE, np
    INTEGER, DIMENSION(:), POINTER :: cblk, indblk

    nbl = struct%num_PLs
    ncont = struct%num_conts
    cblk => struct%cblk
    indblk => struct%mat_PL_start
    iE = pnegf%iE
    np = pnegf%local_en_points

    ! Take CSR H,S and build ES-H in dense blocks
    CALL prealloc_sum(pnegf%H,pnegf%S,minusone,E,ESH_csr)
    call allocate_blk_dns(ESH,nbl)
    CALL csr2blkdns(ESH_csr,ESH,indblk)
    CALL destroy(ESH_csr)
    
    ! ADDING CONTACT GF
    DO i=1,ncont
       cb = cblk(i)
       ESH(cb,cb)%val = ESH(cb,cb)%val-SelfEneR(i)%val
    ENDDO

    ! ADDING PH SELF ENERGY
    if (pnegf%phph%include_phph .and. pnegf%phph%scba_iter.gt.0) then 
      call add_sigma_ph_r(pnegf, ESH)
    end if
 
    ! BLOCK-ITERATIVE ALGORITHM (Build Gr tridiagonal blocks)
    call allocate_blk_dns(Gr,nbl)
    call allocate_gsm(gsmr,nbl)
    rbl = minval(cblk(1:ncont)) + 1  
    lbl = maxval(cblk(1:ncont)) - 1

    ! -------------------------------------------------------------
    ! MESSY PART TO COMPUTE THE COLUMNS of Gr AT THE CONTACTS
    ! Needed for:  Gn = Gr(i,cb) Gamma(cb,cb) Ga(cb,j) 
    ! NOTE the behaviour of Make_Gr_mem:
    !      (ESH,i)   computes block i,i by inversion
    !      (ESH,i,i) computes block i,i by iteration
    ! -------------------------------------------------------------
     
    IF (comp_col) THEN
      call allocate_gsm(gsml,nbl)

      ! Fix to a bug when there are 2PLs
      ! later Make_Gr tries to compute Gr(1,1) but needs gsmr(2,2)
      ! 
      IF (nbl.eq.2) then
        CALL calculate_gsmr_blocks(ESH,nbl,rbl-1)
        CALL calculate_gsml_blocks(ESH,1,lbl+1)    
      ELSE
        CALL calculate_gsmr_blocks(ESH,nbl,rbl)
        CALL calculate_gsml_blocks(ESH,1,lbl)    
      ENDIF

      ! 1. rbl>lbl  => lbl+1=rbl-1 => compute first Gr(rbl-1,rbl-1)
      ! 2. rbl<lbl  => lbl=rbl-2 has been computed
      ! Make_Gr does not compute if sbl>nbl or sbl<1
      CALL calculate_Gr_tridiag_blocks(ESH,rbl-1)
      CALL calculate_Gr_tridiag_blocks(ESH,rbl,nbl)
      CALL calculate_Gr_tridiag_blocks(ESH,rbl-2,1)
 
      ! build contact column green's for later use 
      DO i=1,ncont
         CALL calculate_Gr_column_blocks(ESH,cblk(i),indblk)
      ENDDO

      call destroy_gsm(gsml)
      call deallocate_gsm(gsml)

    ELSE
      
      CALL calculate_gsmr_blocks(ESH,nbl,2)
      CALL calculate_Gr_tridiag_blocks(ESH,1)
      CALL calculate_Gr_tridiag_blocks(ESH,2,nbl)

    ENDIF
 
    CALL destroy_ESH(ESH)
    CALL deallocate_blk_dns(ESH)
    call destroy_gsm(gsmr)
    call deallocate_gsm(gsmr)

    ! SAVE ON FILES/MEMORY (for phph).........................
    IF (pnegf%phph%include_phph) THEN
      DO i = 1, nbl
         !print*,'G_r ',minval(abs(Gr(i,i)%val)), maxval(abs(Gr(i,i)%val))
         call write_blkmat(Gr(i,i),pnegf%scratch_path,'G_r_',i,i,iE,np)
      ENDDO
      DO i = 1, nbl-1
         call write_blkmat(Gr(i,i+1),pnegf%scratch_path,'G_r_',i,i+1,iE,np)
         call write_blkmat(Gr(i+1,i),pnegf%scratch_path,'G_r_',i+1,i,iE,np)
      ENDDO
    END IF
   
    ! SAVES COLUMN Gr (needed for Gn =  Gr Gamma Ga) 
    IF (comp_col) THEN
      DO m=1,ncont
        cb = cblk(m)
        DO i = 1, cb-2
          call write_blkmat(Gr(i,cb),pnegf%scratch_path,'G_r_',i,cb,iE,np)
        ENDDO
        DO i = cb+2, nbl
          call write_blkmat(Gr(i,cb),pnegf%scratch_path,'G_r_',i,cb,iE,np)
        ENDDO
      ENDDO
    END IF
   
    !..........................................................
    if (present(A)) then
       call blk2csr(Gr,struct,pnegf%S,A)
    endif

    CALL destroy_blk(Gr)
    DEALLOCATE(Gr)

  END SUBROUTINE calls_Dr_ph

  !****************************************************************************
  !
  ! Driver for computing G_n = (+/-)iG< contributions including interactions
  !
  !    Sum   f_j(E) Gr Gam_j Ga +   Gr Sigma_ph< Ga
  !     j
  !
  ! NOTE: 
  !
  !****************************************************************************
  SUBROUTINE calls_Dn_ph(pnegf,E,SelfEneR,nB,struct,Gout)
    
    !****************************************************************************
    !
    !Input
    !pnegf:    negf data container
    !E:        Energy point
    !SelfEneR: matrices array containing contacts Self Energy
    !nB:       bose distribution of each contact
    !struct:   structure container to get nPls, nConts, indblk,
    !A:        optional, gives back the whole Gr masked by the overlap  
    !
    !*****************************************************************************

    IMPLICIT NONE

    !In/Out
    TYPE(Tnegf), intent(in) :: pnegf
    REAL(dp), intent(in)  :: E
    TYPE(z_DNS), DIMENSION(:), intent(in)  :: SelfEneR
    real(dp), DIMENSION(:), intent(in)  :: nB
    TYPE(Tstruct_info), intent(in)  :: struct
    TYPE(z_CSR), intent(inout), optional  :: Gout

    !Work
    COMPLEX(dp) :: Ec
    INTEGER :: i, k, ierr,ncont,nbl, lbl, rbl, ref, iE, np
    INTEGER, DIMENSION(:), POINTER :: cblk, indblk
    TYPE(z_DNS), DIMENSION(:,:), ALLOCATABLE :: ESH
    TYPE(z_DNS), DIMENSION(:,:), ALLOCATABLE :: Gn
    TYPE(z_CSR) :: ESH_csr, Gl
    LOGICAL :: mask(MAXNCONT)

    nbl = struct%num_PLs
    ncont = struct%num_conts
    indblk => struct%mat_PL_start
    cblk => struct%cblk
    ref = pnegf%refcont
    iE = pnegf%iE
    np = pnegf%local_en_points

    Ec=cmplx(E, 0.0_dp, dp)


    DO i=1,ncont
       ESH(cblk(i),cblk(i))%val = ESH(cblk(i),cblk(i))%val-SelfEneR(i)%val
    ENDDO

    call allocate_blk_dns(Gn,nbl)

    ! COMPUTE: Gn =  Sum_i [ Gr Gamma_i Ga nB_i ] 
    call init_blkmat(Gn,ESH)
    CALL calculate_Gn_tridiag_blocks(ESH,SelfEneR,nB,ref,struct,Gn)

    ! ADD THE SCATTERING SELF-ENERGY:  Gn = Gr Sigma_n Ga
    IF (pnegf%phph%include_phph .and. pnegf%phph%scba_iter.gt.0) THEN
       call calculate_Gn_tridiag_elph_contributions(pnegf,ESH,Gn)
    ENDIF

    ! SAVE tridiagonal blocks of G_n
    IF (pnegf%phph%include_phph) THEN
      DO i = 1, nbl
        !print*,'(G_n) G_n',minval(abs(Gn(i,i)%val)), maxval(abs(Gn(i,i)%val))
        call write_blkmat(Gn(i,i),pnegf%scratch_path,'G_n_',i,i,iE,np)
      ENDDO
      DO i = 1, nbl-1
        call write_blkmat(Gr(i,i+1),pnegf%scratch_path,'G_n_',i,i+1,iE,np)
        call write_blkmat(Gr(i+1,i),pnegf%scratch_path,'G_n_',i+1,i,iE,np)
      ENDDO
    ENDIF

    if (present(Gout)) then
      call blk2csr(Gn,struct,pnegf%S,Gout)
    endif
 
    CALL destroy_blk(Gn)
    DEALLOCATE(Gn)

    CALL destroy_blk(Gr)
    DEALLOCATE(Gr)

    CALL destroy_ESH(ESH)
    DEALLOCATE(ESH)

  END SUBROUTINE calls_Dn_ph

  !***********************************************************************
  !
  !  FEW utility subrutines to allocate/deallocate block-dense matrices
  ! 
  !**********************************************************************
  SUBROUTINE init_blkmat(Matrix,S)
    type(z_DNS), dimension(:,:) :: Matrix,S

    integer :: nbl, j

    nbl = SIZE(Matrix,1)

    call create(Matrix(1,1),S(1,1)%nrow,S(1,1)%ncol)      
    Matrix(1,1)%val=zero
    DO j=2,nbl-1
       call create(Matrix(j-1,j),S(j-1,j)%nrow,S(j-1,j)%ncol)      
       Matrix(j-1,j)%val=zero
       call create(Matrix(j,j),S(j,j)%nrow,S(j,j)%ncol)      
       Matrix(j,j)%val=zero
       call create(Matrix(j,j-1),S(j,j-1)%nrow,S(j,j-1)%ncol)      
       Matrix(j,j-1)%val=zero
    ENDDO
    IF (nbl.gt.1) then
       call create(Matrix(nbl,nbl),S(nbl,nbl)%nrow,S(nbl,nbl)%ncol)      
       Matrix(nbl,nbl)%val=zero
    ENDIF

  END SUBROUTINE init_blkmat  


  !**********************************************************************
  SUBROUTINE destroy_gsm(gsm)
    type(z_DNS), DIMENSION(:) :: gsm
    integer :: i, i1, nbl

    nbl=size(gsm,1)

    do i=1,nbl
      if (allocated(gsm(i)%val)) call destroy(gsm(i))
    enddo

  END SUBROUTINE destroy_gsm

  !**********************************************************************
  SUBROUTINE destroy_blk(M)
    type(z_DNS), DIMENSION(:,:) :: M
    integer :: i, i1, nbl

    nbl=size(M,1)

    DO i=1,nbl
       DO i1=1,nbl
          IF (ALLOCATED(M(i1,i)%val)) THEN
              !print*,'kill Gr',i1,i
              CALL destroy(M(i1,i))
          ENDIF
       ENDDO
    ENDDO

  END SUBROUTINE destroy_blk

  !**********************************************************************
  SUBROUTINE destroy_ESH(ESH)

    integer :: i, nbl
    type(z_DNS), dimension(:,:) :: ESH

    nbl=size(ESH,1)

    DO i=1,nbl
       CALL destroy(ESH(i,i))
    ENDDO
    DO i=2,nbl
       CALL destroy(ESH(i-1,i))
       CALL destroy(ESH(i,i-1))
    ENDDO

  END SUBROUTINE destroy_ESH

  !**********************************************************************
  !
  !  Divides a sparse matrix A_csr into an array of dense matrices, 
  !  A(nbl,nbl) 
  !
  !**********************************************************************

  SUBROUTINE csr2blkdns(A_csr,A,indblk)

    !**********************************************************************
    !Input:
    !ESH_csr: sparse matrix ES-H related to device
    !
    !Output:
    !ESH(nbl,nbl): dense matrix array -> single matrices allocated 
    !              internally, array ESH(nbl,nbl) allocated externally
    !**********************************************************************

    IMPLICIT NONE 

    INTEGER :: i
    TYPE(z_CSR) :: A_csr
    INTEGER :: nbl
    TYPE(z_DNS), DIMENSION(:,:) :: A
    INTEGER, DIMENSION(:) :: indblk

    nbl = size(A,1)

    DO i=1,nbl
       CALL extract(A_csr,indblk(i),indblk(i+1)-1,indblk(i),indblk(i+1)-1,A(i,i))
    ENDDO

    DO i=2,nbl
       CALL extract(A_csr,indblk(i-1),indblk(i)-1,indblk(i),indblk(i+1)-1,A(i-1,i))
       CALL extract(A_csr,indblk(i),indblk(i+1)-1,indblk(i-1),indblk(i)-1,A(i,i-1))
    ENDDO

  END SUBROUTINE csr2blkdns

  !***********************************************************************
  !
  !  Reloads the retarded elph self-energy and add it to ES-H
  !
  !***********************************************************************
  SUBROUTINE add_sigma_ph_r(pnegf, ESH)
     TYPE(Tnegf), intent(in) :: pnegf
     TYPE(z_DNS), DIMENSION(:,:), intent(inout) :: ESH

     TYPE(z_DNS), DIMENSION(:,:), ALLOCATABLE :: Sigma_ph_r
     INTEGER :: n, nbl, nrow, ierr

     nbl = pnegf%str%num_PLs
    
     ALLOCATE(Sigma_ph_r(nbl,nbl),stat=ierr)
     IF (ierr.NE.0) STOP 'ALLOCATION ERROR: could not allocate Sigma_ph_r'


     DO n = 1, nbl
        call create(Sigma_ph_r(n,n), ESH(n,n)%nrow, ESH(n,n)%nrow)
        Sigma_ph_r(n,n)%val = zero
        call read_blkmat(Sigma_ph_r(n,n),pnegf%scratch_path,'Sigma_ph_r_',n,n,pnegf%iE)
        ESH(n,n)%val = ESH(n,n)%val - Sigma_ph_r(n,n)%val
        call destroy(Sigma_ph_r(n,n))
     END DO

     if (.not. pnegf%phph%diagonal) THEN
       DO n = 1, nbl-1       
          call create(Sigma_ph_r(n,n+1), ESH(n,n+1)%nrow, ESH(n,n+1)%ncol)
          Sigma_ph_r(n,n+1)%val = zero
          call read_blkmat(Sigma_ph_r(n,n+1),pnegf%scratch_path,'Sigma_ph_r_',n,n+1,pnegf%iE)
          ESH(n,n+1)%val = ESH(n,n+1)%val - Sigma_ph_r(n,n+1)%val
          call destroy(Sigma_ph_r(n,n+1))
          
          call create(Sigma_ph_r(n+1,n), ESH(n+1,n)%nrow, ESH(n+1,n)%ncol)
          Sigma_ph_r(n+1,n)%val = zero
          call read_blkmat(Sigma_ph_r(n+1,n),pnegf%scratch_path,'Sigma_ph_r_',n+1,n,pnegf%iE)
          ESH(n+1,n)%val = ESH(n+1,n)%val - Sigma_ph_r(n+1,n)%val
          call destroy(Sigma_ph_r(n+1,n))
       ENDDO
     ENDIF

     DEALLOCATE(Sigma_ph_r)

  END SUBROUTINE add_sigma_ph_r 


  !***********************************************************************
  !
  !  g_small right (gsmr) calculation - write on memory
  !
  !***********************************************************************

  SUBROUTINE calculate_gsmr_blocks(ESH,sbl,ebl)

    !***********************************************************************
    !Input:
    !ESH: sparse matrices array ESH(nbl,nbl)
    !
    !nbl (number of layers)  global needed   (indblk not anymore)
    !
    !Output:
    !sparse matrices array global variable gsmr(nbl) is available in memory
    !single blocks are allocated internally, array Gr(nbl,nbl) 
    !must be allocated externally
    !***********************************************************************

    IMPLICIT NONE

    !In/Out
    TYPE(z_DNS), DIMENSION(:,:), intent(in) :: ESH
    integer, intent(in) :: sbl,ebl                 ! start block, end block

    !Work
    !TYPE(z_DNS), DIMENSION(:,:), ALLOCATABLE :: INV
    TYPE(z_DNS) :: work1, work2
    INTEGER :: nrow, M, N
    INTEGER :: i, nbl

    if (sbl.lt.ebl) return

    nbl = size(ESH,1)

    if (nbl.gt.1) then

       nrow=ESH(sbl,sbl)%nrow

       call create(gsmr(sbl),nrow,nrow)

       call compGreen(gsmr(sbl),ESH(sbl,sbl),nrow)

    endif

    DO i=sbl-1,ebl,-1

       CALL prealloc_mult(ESH(i,i+1),gsmr(i+1),minusone,work1)
       
       CALL prealloc_mult(work1,ESH(i+1,i),work2)
       
       CALL destroy(work1)

       CALL prealloc_sum(ESH(i,i),work2,work1)

       CALL destroy(work2)

       CALL create(gsmr(i),work1%nrow,work1%nrow)
       
       CALL compGreen(gsmr(i),work1,work1%nrow)
       
       CALL destroy(work1)

    ENDDO


  END SUBROUTINE calculate_gsmr_blocks

 


  !***********************************************************************
  !
  !  g_small left (gsml) calculation - write on memory
  !
  !***********************************************************************

  SUBROUTINE calculate_gsml_blocks(ESH,sbl,ebl)

    !***********************************************************************
    !Input:
    !ESH: sparse matrices array ESH(nbl,nbl)
    !
    !nbl (number of layers), indblk(nbl+1) global variables needed
    !
    !Output:
    !sparse matrices array global variable gsml(nbl) is available in memory
    !single blocks are allocated internally, array Gr(nbl,nbl) 
    !must be allocated externally
    !***********************************************************************


    IMPLICIT NONE

    !In/Out
    TYPE(z_DNS), DIMENSION(:,:), intent(in) :: ESH
    integer, intent(in) :: sbl,ebl                       ! start block, end block

    !Work
    TYPE(z_DNS) :: work1, work2
    INTEGER :: nrow
    INTEGER :: i, nbl
   !TYPE(z_DNS) :: INV(sbl,sbl)

    if (sbl.gt.ebl) return

    !***
    !gsml(sbl)
    !***
    nbl = size(ESH,1)
    nrow=ESH(sbl,sbl)%nrow  

    CALL create(gsml(sbl),nrow,nrow)

    call compGreen(gsml(sbl),ESH(sbl,sbl),nrow)


    DO i=sbl+1,ebl

       nrow=ESH(i,i)%nrow   !indblk(i+1)-indblk(i)

       CALL prealloc_mult(ESH(i,i-1),gsml(i-1),minusone,work1)

       CALL prealloc_mult(work1,ESH(i-1,i),work2)

       CALL destroy(work1)

       CALL prealloc_sum(ESH(i,i),work2,work1)

       CALL destroy(work2)

       call create(gsml(i),work1%nrow,work1%nrow)

       call compGreen(gsml(i),work1,work1%nrow)

       CALL destroy(work1)

    ENDDO

    if (debug) then
       WRITE(*,*) '********************'
       WRITE(*,*) 'Make_gsml_mem done'
       WRITE(*,*) '********************'
    endif

  END SUBROUTINE calculate_gsml_blocks





  !***********************************************************************
  !
  !  Diagonal, Subdiagonal, Superdiagonal blocks of Green Retarded 
  !  Gr(nbl,nbl) - writing on memory
  !
  !*********************************************************************** 

  SUBROUTINE calculate_Gr_tridiag_blocks(ESH,sbl,ebl)

    !***********************************************************************
    !Input:
    !ESH: dense matrices array ESH(nbl,nbl)
    !sbl, ebl : block indexes
    ! If only sbl is specified, it calculates Gr(sbl, sbl)
    ! If sbl > ebl, it calculates Gr(ebl:sbl, ebl:sbl), Gr(ebl:sbl + 1, ebl:sbl),
    !    Gr(ebl:sbl, ebl:sbl + 1) (need gsml)          
    ! If sbl < ebl, it calculates Gr(sbl:ebl, sbl:ebl), Gr(sbl:ebl - 1, sbl:ebl),
    !    Gr(sbl:ebl, sbl:ebl - 1) (need gsmr)          
    !
    ! 
    !Output:
    !sparse matrices array global variable Gr(nbl,nbl) is available in 
    !memory - single blocks are allocated internally, array Gr(nbl,nbl) 
    !must be allocated externally
    !***********************************************************************

    IMPLICIT NONE 

    !In/Out
    TYPE(z_DNS), DIMENSION(:,:) :: ESH
    INTEGER :: sbl
    INTEGER, optional :: ebl

    !Work
    INTEGER :: i,nrow,nbl
    TYPE(z_DNS) :: work1, work2, work3

    nbl = size(ESH,1)
    
    if (sbl.gt.nbl) return 
    if (sbl.lt.1) return 

    if (.not.present(ebl)) then
       if (nbl.eq.1) then
          nrow = ESH(sbl,sbl)%nrow     
          CALL create(Gr(sbl,sbl),nrow,nrow)
          CALL compGreen(Gr(sbl,sbl),ESH(sbl,sbl),nrow)
       else
          nrow = ESH(sbl,sbl)%nrow     
          call create(work1,nrow,nrow)
          work1%val = ESH(sbl,sbl)%val
          if (sbl+1.le.nbl) then
            CALL prealloc_mult(ESH(sbl,sbl+1),gsmr(sbl+1),work2)
            CALL prealloc_mult(work2,ESH(sbl+1,sbl),work3)
            CALL destroy(work2)
            CALL prealloc_sum(work1,work3,minusone,work2)
            CALL destroy(work3)
            work1%val = work2%val
            CALL destroy(work2)
          endif 
          if (sbl-1.ge.1) then
            CALL prealloc_mult(ESH(sbl,sbl-1),gsml(sbl-1),work2)
            CALL prealloc_mult(work2,ESH(sbl-1,sbl),work3)
            CALL destroy(work2)
            CALL prealloc_sum(work1,work3,minusone,work2)
            CALL destroy(work3)
            work1%val = work2%val
            CALL destroy(work2)
          endif  
            
          CALL create(Gr(sbl,sbl),nrow,nrow)
          CALL compGreen(Gr(sbl,sbl),work1,nrow)
          CALL destroy(work1)
       endif   
       return
    endif


    !***
    !Diagonal, Subdiagonal and Superdiagonal blocks
    !***
    IF ((ebl.ge.sbl).and.(ebl.gt.1).and.(sbl.gt.1)) THEN
       DO i=sbl,ebl,1
          CALL prealloc_mult(gsmr(i),ESH(i,i-1),work1)
          CALL prealloc_mult(work1,Gr(i-1,i-1),minusone,Gr(i,i-1))
          CALL destroy(work1)
          
          CALL prealloc_mult(ESH(i-1,i),gsmr(i),work2)
          CALL prealloc_mult(Gr(i-1,i-1),work2,minusone,Gr(i-1,i))

          CALL prealloc_mult(Gr(i,i-1),work2,minusone,work1)
          CALL destroy(work2)
          
          CALL prealloc_sum(gsmr(i),work1,Gr(i,i))
          CALL destroy(work1) 
       ENDDO
    ELSE
       DO i=sbl,ebl,-1
          CALL prealloc_mult(gsml(i),ESH(i,i+1),work1)
          CALL prealloc_mult(work1,Gr(i+1,i+1),minusone,Gr(i,i+1))
          CALL destroy(work1)
          
          CALL prealloc_mult(ESH(i+1,i),gsml(i),work2)
          CALL prealloc_mult(Gr(i+1,i+1),work2,minuesone,Gr(i+1,i))

          CALL prealloc_mult(Gr(i,i+1),work2,minusone,work1)
          CALL destroy(work2)
          
          CALL prealloc_sum(gsml(i),work1,Gr(i,i))
          CALL destroy(work1) 
       ENDDO
    ENDIF 

  END SUBROUTINE calculate_Gr_tridiag_blocks

  !**************************************************************************
  !
  !  Calculate Green Retarded column "n" - writing on memory 
  !
  !**************************************************************************

  SUBROUTINE calculate_Gr_column_blocks(ESH,n,indblk)

    !***********************************************************************
    !Input:
    !ESH: sparse matrices array ESH(nbl,nbl)
    !n: n umber of column to be calculated
    !
    !global variables needed:  
    !Gr diagonal, subadiagonal and superdiagonal, 
    !gsmr(:) for downgoing and gsml(:) for upgoing 
    !
    !Output:
    !sparse matrices array global variable Gr(:,n) is available in 
    !memory - single blocks are allocated internally, array Gr(nbl,nbl) 
    !must be allocated externally
    !*********************************************************************** 

    IMPLICIT NONE

    !In/Out
    TYPE(z_DNS), DIMENSION(:,:), intent(in) :: ESH
    INTEGER, intent(in) :: n
    INTEGER, DIMENSION(:), intent(in) :: indblk 

    !Work
    INTEGER :: i,nrow,ncol,nbl
    TYPE(z_DNS) :: work1
    REAL(dp) :: max

    nbl = size(ESH,1)

    IF (n.GT.nbl) THEN
       STOP 'Error in Make_Grcol_mem : n is greater than nbl'
    ENDIF

    !***************************************
    !  Downgoing (j>=n+2 && n<nbl-1)
    !
    !   G_j,n = -gR_jj T_j,j-1 G_j-1,n
    !
    !***************************************
    ncol=indblk(n+1)-indblk(n)

    IF (n.LT.(nbl-1)) THEN

       DO i=n+2,nbl

          nrow=indblk(i+1)-indblk(i)

          max=MAXVAL(ABS(Gr(i-1,n)%val))

          IF (max.GT.EPS) THEN
             CALL prealloc_mult(gsmr(i),ESH(i,i-1),minusone,work1)
             CALL prealloc_mult(work1,Gr(i-1,n),Gr(i,n))
             CALL destroy(work1)
          ENDIF 

       ENDDO
    
    ENDIF
    !*************************************
    !   Up-going (j<=n-2 && n>2)
    !
    !   G_j,n = -gL_jj T_j,j+1 G_j+1,n
    !
    !*************************************

    IF (n.GT.2) THEN

       DO i=n-2,1,(-1)
          nrow=indblk(i+1)-indblk(i)

          max=MAXVAL(ABS(Gr(i+1,n)%val))

          IF (max.GT.EPS) THEN
             CALL prealloc_mult(gsml(i),ESH(i,i+1),minusone,work1)
             CALL prealloc_mult(work1,Gr(i+1,n),Gr(i,n))
             CALL destroy(work1)
          ENDIF
       ENDDO

    ENDIF

    if (debug) then
       WRITE(*,*) '********************'
       WRITE(*,*) 'Make_Grcol_mem done column',n
       WRITE(*,*) '********************'
    endif

  END SUBROUTINE calculate_Gr_column_blocks

  !****************************************************************************
  !
  !  Calculate Green Retarded - writing on memory (optimized on mask)
  !
  !****************************************************************************

  SUBROUTINE Gr_blk2csr(P,nbl,indblk,A)

    !****************************************************************************
    !Input:
    !P: CSR matrix containing masking pattern
    !
    !global variable needed: nbl, indblk(nbl+1), Gr(:,:) 
    !
    !Output:
    !A: sparse matrix containing Green Retarded of device (allocated internally)
    !****************************************************************************

    IMPLICIT NONE

    !In/Out
    INTEGER :: nbl
    INTEGER, DIMENSION(:), POINTER :: indblk
    TYPE(z_CSR) :: A, P, GrCsr

    !Work
    INTEGER :: i, j, i1, ix, iy, x, y, col, oldx

    !create A with same pattern of P
    CALL create(A,P%nrow,P%ncol,P%nnz)
    A%rowpnt(:)=P%rowpnt(:)
    A%colind(:)=P%colind(:)
    A%nzval = zero 

    !If only one block is present, concatenation is not needed 
    !and it's implemented in a more trivial way
    IF (nbl.EQ.1) THEN

       call create(GrCsr,Gr(1,1)%nrow,Gr(1,1)%ncol,Gr(1,1)%nrow*Gr(1,1)%ncol)
       call dns2csr(Gr(1,1),GrCsr)

       call mask(GrCsr,P,A)
       call destroy(GrCsr)

    ELSE  

       !Cycle upon all rows
       x = 1    
       DO i = 1, A%nrow
          !Choose which block (row) we're dealing with
          oldx = x

          !Check if row is in same block of previous or in next block. Not needed 
          !(and not allowed not to exceed indblk index boundaries) if we're in the last block
          IF (oldx.EQ.nbl) THEN 
             x = oldx
          ELSE
             DO ix = oldx, oldx+1
                IF ( (i.GE.indblk(ix)).AND.(i.LT.indblk(ix+1)) ) x = ix
             ENDDO
          ENDIF

          !Offset: i1 is the index for separate blocks
          i1 = i - indblk(x) + 1
          !Cycle upon columns 
          DO j = A%rowpnt(i), A%rowpnt(i+1) -1
             !Choose which block column we're dealing with
             y = 0
             IF (x.EQ.1) THEN
                IF ( (A%colind(j).GE.indblk(x)).AND.(A%colind(j).LT.indblk(x + 1)) ) then 
                   y = 1
                ELSEIF ( (A%colind(j).GE.indblk(x + 1)).AND.(A%colind(j).LT.indblk(x + 2)) ) then 
                   y = 2
                ENDIF
             elseif (x.eq.nbl) then
                IF ( (A%colind(j).GE.indblk(x)).AND.(A%colind(j).LT.indblk(x + 1)) ) then 
                   y = nbl
                ELSEIF ( (A%colind(j).GE.indblk(x - 1)).AND.(A%colind(j).LT.indblk(x)) ) then 
                   y = nbl - 1
                ENDIF
             ELSE
                DO iy = x-1, x+1
                   if ( (A%colind(j).GE.indblk(iy)).AND.(A%colind(j).LT.indblk(iy + 1)) ) y = iy
                ENDDO
             ENDIF
             IF (y.eq.0) then
                write(*,*)     
                write(*,*) 'ERROR in blk2csr: probably wrong PL size',x
                write(*,*) 'row',i,A%colind(j)
                write(*,*) 'block indeces:',indblk(1:nbl)
                stop
             ENDIF

             col = A%colind(j) - indblk(y) + 1

             A%nzval(j) = Gr(x,y)%val(i1,col) 

          ENDDO

       ENDDO

    ENDIF

    !if (debug) call writePeakInfo(6)    
    if (debug) then
       WRITE(*,*) '**********************'
       WRITE(*,*) 'Make_GreenR_mem done'
       WRITE(*,*) '**********************'
    endif

  END SUBROUTINE Gr_blk2csr


  !****************************************************************************
  !
  ! Calculate G_n contributions for all contacts (except reference)
  ! Writing on memory
  !
  !****************************************************************************

  SUBROUTINE calculate_Gn_tridiag_blocks(ESH,SelfEneR,nX,ref,struct,Gn)

    !******************************************************************************
    !Input:  
    !ESH(nbl,nbl): sparse matrices array ES-H
    !SelfEneR(ncont): sparse matrices array with contacts Self Energy 
    !frm(ncont): Fermi diistribution value for each contact 
    !ref:  reference contact 
    !
    !global variables needed: nbl, indblk(nbl+1), ncont, cblk(ncont), Gr(:,:) 
    !diagonal, subdiagonal, overdiagonal and in colums Gr(:,cb) where cb are
    !layers interacting with all contacts but collector
    !
    !Output:
    !Gl: sparse matrix containing G  contacts contribution
    !*******************************************************************************  

    IMPLICIT NONE

    !In/Out
    TYPE(z_DNS), DIMENSION(:,:), intent(in) :: ESH
    TYPE(z_DNS), DIMENSION(:,:), intent(inout) :: Gn
    TYPE(z_DNS), DIMENSION(:), intent(in) :: SelfEneR
    TYPE(Tstruct_info), intent(in) :: struct
    REAL(dp), DIMENSION(:), intent(in) :: nX !bosons or fermions 
    INTEGER, intent(in) :: ref
 
    !Work
    Type(z_DNS) :: Gam
    TYPE(z_DNS) :: work1,Ga
    INTEGER :: i,j,cb
    INTEGER :: ncont, nbl
    INTEGER, DIMENSION(:), POINTER :: cblk
    COMPLEX(dp) :: nXdiff

    ncont = struct%num_conts    
    nbl = struct%num_PLs
    cblk => struct%cblk
    
    !*******************************************
    ! Contact Iteration
    !*******************************************
    DO j=1,ncont
  
       ! NOTE: this soubroutine uses prealloc_mult that performs 
       ! C = C + A*B 
       IF (j.NE.ref .AND. ABS(nX(j)-nX(ref)).GT.EPS) THEN
 
          cb=cblk(j) ! block corresponding to contact j

          CALL zspectral(SelfEneR(j),SelfEneR(j),0,Gam)

          nXdiff = cmplx(nX(j)-nX(ref),0.0_dp,dp)
          ! Computation of Gl(1,1) = Gr(1,cb) Gam(cb) Ga(cb,1)
          if (allocated(Gr(1,cb)%val)) then
            CALL prealloc_mult(Gr(1,cb),Gam,work1)    
            CALL zdagger(Gr(1,cb),Ga)
            CALL prealloc_mult(work1,Ga,nXdiff,Gn(1,1))
            CALL destroy(work1, Ga)
          else
            Gn(1,1)%val= zero
          endif      

          ! Computation of all tridiagonal blocks
          DO i=2,nbl
          
             ! Computation of Gl(i,j) = Gr(i,cb) Gam(cb) Ga(cb,j)
             ! Both Gr(i,cb) and Gr(j,cb) must be non-zero
             if (Gr(i-1,cb)%nrow.gt.0 .and. Gr(i,cb)%nrow.gt.0) then
                CALL prealloc_mult(Gr(i-1,cb),Gam,work1)
                CALL zdagger(Gr(i,cb),Ga)
                CALL prealloc_mult(work1,Ga,nXdiff,Gn(i-1,i))
                CALL destroy(work1)
                
                CALL prealloc_mult(Gr(i,cb),Gam,work1)
                CALL prealloc_mult(work1,Ga,nXdiff,Gn(i,i))

                CALL destroy(work1, Ga)
               
                CALL prealloc_mult(Gr(i,cb),Gam,work1)
                CALL zdagger(Gr(i-1,cb),Ga)
                CALL prealloc_mult(work1,Ga,nXdiff,Gn(i,i-1))

                CALL destroy(work1, Ga)
             else
                Gn(i-1,i)%val= zero
                Gn(i,i-1)%val= zero
             endif      

          ENDDO
                  
          call destroy(Gam)

      ENDIF 

    ENDDO

  END SUBROUTINE calculate_Gn_tridiag_blocks


  !****************************************************************************
  !
  ! Calculate G_n contributions due to el-ph
  ! Writing on memory
  ! The subroutine assumes a block-tridiagonal structure for Sigma_ph_n
  ! Computes tridiagonal structure for Gn(i,j): 
  !           --       --
  ! Gn(i,j) = >        >       Gr(i,k) Sigma_ph(k,k+s) Ga(k+s,j)
  !           --       --
  !          k=1,nbl s=-1,0,1
  !****************************************************************************
  SUBROUTINE calculate_Gn_tridiag_elph_contributions(pnegf,ESH,Gn)

    TYPE(Tnegf), intent(in) :: pnegf
    TYPE(z_DNS), DIMENSION(:,:), intent(in) :: ESH
    TYPE(z_DNS), DIMENSION(:,:), intent(inout) :: Gn

    Type(z_DNS), DIMENSION(:,:), ALLOCATABLE :: Sigma_ph_n
    Type(z_DNS) :: Ga, work1, work2
    INTEGER :: n, k, s, nbl, nrow, ncol, ierr

    nbl = pnegf%str%num_PLs
    ALLOCATE(Sigma_ph_n(nbl,nbl),stat=ierr)
    IF (ierr.NE.0) STOP 'ALLOCATION ERROR: could not allocate Sigma_ph_n'

    DO n = 1, nbl
      nrow = ESH(n,n)%nrow
      call create(Sigma_ph_n(n,n), nrow, nrow)
      Sigma_ph_n(n,n)%val = zero
      call read_blkmat(Sigma_ph_n(n,n),pnegf%scratch_path,'Sigma_ph_n_',n,n,pnegf%iE)
    END DO
    DO n = 1, nbl-1 
      nrow = ESH(n,n+1)%nrow
      ncol = ESH(n,n+1)%ncol
      call create(Sigma_ph_n(n,n+1), nrow, ncol)
      call create(Sigma_ph_n(n+1,n), ncol, nrow)
      Sigma_ph_n(n,n+1)%val = zero
      Sigma_ph_n(n+1,n)%val = zero
      call read_blkmat(Sigma_ph_n(n,n+1),pnegf%scratch_path,'Sigma_ph_n_',n,n+1,pnegf%iE)
      call read_blkmat(Sigma_ph_n(n+1,n),pnegf%scratch_path,'Sigma_ph_n_',n+1,n,pnegf%iE)
    END DO

    DO n = 1, nbl
      DO k = 1, nbl
        DO s = -1,+1 
          if (k+s.ge.1 .and. k+s.le.nbl) then
            if (Gr(n,k+s)%nrow.gt.0) then
              CALL prealloc_mult(Gr(n,k+s), Sigma_ph_n(k,k+s), work1)
              CALL zdagger(Gr(n,k+s),Ga)
              CALL prealloc_mult(work1, Ga, work2)
              Gn(n,n)%val = Gn(n,n)%val + work2%val
              call destroy(work2,Ga)
            endif
          endif
          if (k+s.ge.1 .and. k+s.le.nbl .and. n+1.le.nbl) then
            if (Gr(n+1,k+s)%nrow.gt.0) then
               CALL zdagger(Gr(n+1,k+s),Ga)
               CALL prealloc_mult(work1, Ga, work2)
               Gn(n,n+1)%val = Gn(n,n+1)%val + work2%val
               call destroy(work1,work2,Ga)
            endif
          endif
        END DO
      END DO    
      Gn(n+1,n)%val =  Gn(n+1,n)%val + conjg(transpose(Gn(n,n+1)%val))
    END DO
    
    DO n = 1, nbl-1
      CALL destroy(Sigma_ph_n(n,n))
      CALL destroy(Sigma_ph_n(n,n+1))
      CALL destroy(Sigma_ph_n(n+1,n))
    END DO
    CALL destroy(Sigma_ph_n(nbl,nbl))

    DEALLOCATE(Sigma_ph_n)

  END SUBROUTINE calculate_Gn_tridiag_elph_contributions
 

  !****************************************************************************
  !Subroutine used to convert blk-dense to csr format, masked by the matrix P
  !If only one block is present, concatenation is not needed and it's implemented in a
  !more trivial way
  !****************************************************************************
  SUBROUTINE blk2csr(G,struct,P,Gcsr)

    TYPE(z_DNS), DIMENSION(:,:) :: G
    TYPE(Tstruct_info), intent(in) :: struct
    TYPE(z_CSR), intent(in) :: P
    TYPE(z_CSR), intent(out) :: Gcsr

    INTEGER, DIMENSION(:), POINTER :: indblk, cblk
    INTEGER :: nbl, oldx, row, col, iy, ix, x, y, ii, jj, nrows

    nbl = struct%num_PLs
    cblk => struct%cblk
    indblk => struct%mat_PL_start
    nrows = struct%mat_PL_end(nbl)

    !create Gcsr with same pattern of P
    CALL create(Gcsr,P%nrow,P%ncol,P%nnz)
    Gcsr%rowpnt = P%rowpnt
    Gcsr%colind = P%colind
    Gcsr%nzval = zero

    !Cycle upon all rows
    x = 1
    DO ii = 1, nrows
       !Search block x containing row ii
       oldx = x
       IF (oldx.EQ.nbl) THEN 
          x = oldx
       ELSE
          DO ix = oldx, oldx+1
             IF ( (ii.GE.indblk(ix)).AND.(ii.LT.indblk(ix+1)) ) x = ix
          ENDDO
       ENDIF

       !Offset: row is the index for separate blocks
       row = ii - indblk(x) + 1

       !Cycle upon columns of Gcsr (which has been ALREADY MASKED by S)
       DO jj = Gcsr%rowpnt(ii), Gcsr%rowpnt(ii+1) -1
          IF (Gcsr%colind(jj).gt.nrows) CYCLE
          !Choose which block column we're dealing with
          y = 0
          IF (x.eq.1) then
             IF ( (Gcsr%colind(jj).GE.indblk(x)).AND.(Gcsr%colind(jj).LT.indblk(x + 1)) ) then 
                y = 1
             ELSEIF ( (Gcsr%colind(jj).GE.indblk(x + 1)).AND.(Gcsr%colind(jj).LT.indblk(x + 2)) ) then 
                y = 2
             ENDIF
          elseif (x.eq.nbl) then
             IF ( (Gcsr%colind(jj).GE.indblk(x)).AND.(Gcsr%colind(jj).LT.indblk(x + 1)) ) then 
                y = nbl
             ELSEIF ( (Gcsr%colind(jj).GE.indblk(x - 1)).AND.(Gcsr%colind(jj).LT.indblk(x)) ) then 
                y = nbl - 1
             ENDIF
          else
             DO iy = x-1, x+1
                if ( (Gcsr%colind(jj).GE.indblk(iy)).AND.(Gcsr%colind(jj).LT.indblk(iy + 1)) ) y = iy
             ENDDO
          ENDIF

          IF (y.EQ.0) THEN
             write(*,*)     
             write(*,*) 'ERROR in blk2csr: probably wrong PL size', x
             write(*,*) 'row',ii,nrows,Gcsr%colind(jj)
             write(*,*) 'block indeces:',indblk(1:nbl)
             STOP
          ENDIF
          
          col = Gcsr%colind(jj) - indblk(y) + 1

          IF (allocated(G(x,y)%val)) THEN
              Gcsr%nzval(jj) = Gcsr%nzval(jj) + G(x,y)%val(row,col)
          ENDIF

       ENDDO

    ENDDO


  END SUBROUTINE blk2csr

 !---------------------------------------------------------
 !---------------------------------------------------------
 SUBROUTINE create_scratch(nbl, npoints)
   integer :: nbl, npoints

   integer :: i,j,k,err
   
   call destroy_scratch(nbl, npoints)
   
   ALLOCATE(DDn(nbl,nbl),stat=err)
   ALLOCATE(DDr(nbl,nbl),stat=err)
   ALLOCATE(DDp(nbl,nbl),stat=err)
   ALLOCATE(Pin(nbl,nbl),stat=err)
   ALLOCATE(Pir(nbl,nbl),stat=err)
   ALLOCATE(Pip(nbl,nbl),stat=err)
  
   ! Initialize everything to 0 
   do j=1,nbl
     do i=1,nbl  
         DDn(i,j)%nrow=0
         DDn(i,j)%ncol=0
         DDn(i,j)%npoints=0
         DDr(i,j)%nrow=0
         DDr(i,j)%ncol=0
         DDr(i,j)%npoints=0
         DDp(i,j)%nrow=0
         DDp(i,j)%ncol=0
         DDp(i,j)%npoints=0
         Pin(i,j)%nrow=0
         Pin(i,j)%ncol=0
         Pin(i,j)%npoints=0
         Pir(i,j)%nrow=0
         Pir(i,j)%ncol=0
         Pir(i,j)%npoints=0
         Pip(i,j)%nrow=0
         Pip(i,j)%ncol=0
         Pip(i,j)%npoints=0
     enddo
   enddo
   
   if(err.ne.0) then
      STOP 'ERROR: Cannot allocate GG'
   endif     
   print*,'Created memory scratch',nbl,'x',nbl,'x',npoints

 END SUBROUTINE create_scratch
 !---------------------------------------------------------
 !---------------------------------------------------------

 SUBROUTINE destroy_scratch(nbl, npoints)
   integer :: nbl, npoints
   
   integer :: i, j, iE,  err
   
   err = 0
   
   
   do i = 1, nbl
     do j = 1, nbl
         
        if (allocated(DDn)) then
            if (allocated(DDn(i,j)%val)) call destroy(DDn(i,j))
        endif       
        if (allocated(DDr)) then
            if (allocated(DDr(i,j)%val)) call destroy(DDr(i,j))
        endif       
        if (allocated(DDp)) then
            if (allocated(DDp(i,j)%val)) call destroy(DDp(i,j))
        endif       
        if (allocated(Pin)) then
            if (allocated(Pin(i,j)%val)) call destroy(Pin(i,j))     
        endif       
        if (allocated(Pir)) then
            if (allocated(Pir(i,j)%val)) call destroy(Pir(i,j))
        endif       
        if (allocated(Pip)) then
            if (allocated(Pip(i,j)%val)) call destroy(Pip(i,j))     
        endif       
       
     end do
   end do
   
   if (allocated(DDn)) DEALLOCATE(DDn,stat=err)
   if (allocated(DDr)) DEALLOCATE(DDr,stat=err)
   if (allocated(DDp)) DEALLOCATE(DDp,stat=err)
   if (allocated(Pin)) DEALLOCATE(Pin,stat=err)
   if (allocated(Pir)) DEALLOCATE(Pir,stat=err)
   if (allocated(Pip)) DEALLOCATE(Pip,stat=err)
   
   if(err.ne.0) then
      STOP 'ERROR: Cannot deallocate GG'
   endif

 END SUBROUTINE destroy_scratch
 !---------------------------------------------------------
 !---------------------------------------------------------

 ! READ Matrices
 SUBROUTINE read_blkmat(Matrix, path, name, i, j, iE)
    TYPE(z_DNS), intent(inout) :: Matrix
    CHARACTER(*), intent(in) :: path
    CHARACTER(*), intent(in) :: name
    INTEGER, intent(in) :: i, j, iE

    CHARACTER(64) :: filename
    character(4) :: ofblki, ofblkj
    character(10) :: ofpnt
    logical :: lex
    integer :: i1,i2
    complex(dp) :: mat_el

    if (memory) then
       select case(trim(name))
       case('G_n_')     
          Matrix%val = DDn(i,j)%val(:,:,iE) 
       case('G_r_')                   
          Matrix%val = DDr(i,j)%val(:,:,iE)
       case('G_p_')                   
          Matrix%val = DDp(i,j)%val(:,:,iE)
       case('Sigma_ph_r_')          
          Matrix%val = Pir(i,j)%val(:,:,iE) 
       case('Sigma_ph_n_')          
          Matrix%val = Pin(i,j)%val(:,:,iE)
       case('Sigma_ph_p_')          
          Matrix%val = Pip(i,j)%val(:,:,iE)
       case default 
         stop 'internal error: read_blkmat does not correspond'
       end select
       return
    endif        

    Matrix%val = (0.0_dp, 0.0_dp)

    if (i.le.9999) write(ofblki,'(i4.4)') i
    if (i.gt.9999) stop 'ERROR: too many blks (> 9999)'
    if (j.le.9999) write(ofblkj,'(i4.4)') j 
    if (j.gt.9999) stop 'ERROR: too many blks (> 9999)'

    if (iE.le.99999) write(ofpnt,'(i5.5)') iE

    filename = trim(name)//trim(ofblki)//'_'//trim(ofblkj)//'_'//trim(ofpnt)//'.dat'
    
    inquire(file=trim(path)//trim(filename),EXIST=lex)
    if (.not.lex) then
       RETURN
       !WRITE(*,*) 'ERROR: FILE '//trim(filename)//' DOES NOT EXIST'
       !STOP
    endif   

    open(9091,file=trim(path)//trim(filename), access='STREAM')
    
    call inmat_c(9091,.false.,Matrix%val,Matrix%nrow,Matrix%ncol)

    !open(9001,file=trim(path)//trim(filename), access='DIRECT', recl=4)
    !call direct_in_c(9001,Matrix%val,Matrix%nrow)

    close(9091)

  END SUBROUTINE read_blkmat

  ! WRITE Matrices
  SUBROUTINE write_blkmat(Matrix, path, name, i, j, iE, npoints)
    TYPE(z_DNS), intent(in) :: Matrix
    CHARACTER(*), intent(in) :: path
    CHARACTER(*), intent(in) :: name
    INTEGER, intent(in) :: i, j, iE, npoints

    CHARACTER(64) :: filename
    character(4) :: ofblki, ofblkj
    character(10) :: ofpnt
    integer :: m,n,nrow,ncol

    if (memory) then
              
       nrow=Matrix%nrow
       ncol=Matrix%ncol

       select case(trim(name))
       case('G_n_')     
          if (.not.allocated(DDn(i,j)%val)) then
              call create(DDn(i,j),nrow,ncol,npoints)
          endif    
          DDn(i,j)%val(:,:,iE) = Matrix%val
       case('G_r_')     
          if (.not.allocated(DDr(i,j)%val)) then
              call create(DDr(i,j),nrow,ncol,npoints)
          endif    
          DDr(i,j)%val(:,:,iE) = Matrix%val
       case('G_p_')     
          if (.not.allocated(DDp(i,j)%val)) then
              call create(DDp(i,j),nrow,ncol,npoints)
          endif    
          DDp(i,j)%val(:,:,iE) = Matrix%val 
       case('Sigma_ph_r_')     
          if (.not.allocated(Pir(i,j)%val)) then
              call create(Pir(i,j),nrow,ncol,npoints)
          endif    
          Pir(i,j)%val(:,:,iE) = Matrix%val
       case('Sigma_ph_n_')     
          if (.not.allocated(Pin(i,j)%val)) then
              call create(Pin(i,j),nrow,ncol,npoints)
          endif    
          Pin(i,j)%val(:,:,iE) = Matrix%val
       case('Sigma_ph_p_')     
          if (.not.allocated(Pip(i,j)%val)) then
              call create(Pip(i,j),nrow,ncol,npoints)
          endif    
          Pip(i,j)%val(:,:,iE) = Matrix%val
       case default 
         stop 'internal error: write_blkmat does not correspond'
       end select
       return
    endif        

    if (i.le.9999) write(ofblki,'(i4.4)') i
    if (i.gt.9999) stop 'ERROR: too many blks (> 9999)'
    if (j.le.9999) write(ofblkj,'(i4.4)') j
    if (j.gt.9999) stop 'ERROR: too many blks (> 9999)'

    if (iE.le.99999) write(ofpnt,'(i5.5)') iE

    filename = trim(name)//trim(ofblki)//'_'//trim(ofblkj)//'_'//trim(ofpnt)//'.dat'

    open(9001,file=trim(path)//trim(filename), access='STREAM', status='REPLACE')

    call outmat_c(9001,.false.,Matrix%val,Matrix%nrow,Matrix%ncol) !,1.0d-36)
    
    !open(9001,file=trim(path)//trim(filename), status='REPLACE', access='DIRECT', recl=4)
    !call direct_out_c(9001,Matrix%val,Matrix%nrow)

    close(9001)

  END SUBROUTINE write_blkmat



  !---------------------------------------------------


  subroutine allocate_gsm(gsm,nbl)
    type(z_DNS), dimension(:), allocatable :: gsm 
    integer :: nbl, ierr

    allocate(gsm(nbl),stat=ierr)
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

  !---------------------------------------------------


END MODULE iterative_ph


