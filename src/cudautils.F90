module cudautils
  use ln_precision
  use openacc
  use cublas_v2
  use cusolverDn
  use mat_def
  use, intrinsic :: ieee_arithmetic
  implicit none
  private
  public :: createGPU
  public :: copyToGPU
  public :: copyFromGPU
  public :: deleteGPU
  public :: createAll
  public :: destroyAll

  public :: deleteMatGPU
  public :: createMatGPU

  public :: copy_trid_toGPU
  public :: copy_trid_toHOST
  public :: copy_vdns_toGPU
  public :: copy_vdns_toHOST
  public :: delete_vdns_fromGPU

  public :: init_trid_toGPU
  public :: create_trid_toGPU
  public :: update_trid_toGPU
  public :: delete_trid_fromGPU

  public :: copy_mat_gpu
  public :: matmul_gpu
  public :: add_gpu
  public :: add_cublas
  public :: sum_gpu
  public :: inverse_gpu
  public :: init_mat_gpu
  public :: spectral_gpu
  public :: dagger_gpu
  public :: trace_gpu

  public :: checksum

  interface copyToGPU
     module procedure copyToGPU_dns_sp
     module procedure copyToGPU_dns_dp
     module procedure copyToGPU_comp_vec_sp 
     module procedure copyToGPU_comp_vec_dp
     module procedure copyToGPU_comp_mat_sp
     module procedure copyToGPU_comp_mat_dp
     module procedure copyToGPU_real_sp
     module procedure copyToGPU_real_dp
     module procedure copyToGPU_comp_sp
     module procedure copyToGPU_comp_dp
     module procedure copyToGPU_logic_vec
  end interface copyToGPU

  interface createGPU
     module procedure createGPU_sp
     module procedure createGPU_dp
  end interface createGPU

  interface copyFromGPU
     module procedure copyFromGPU_dns_sp
     module procedure copyFromGPU_dns_dp
     module procedure copyFromGPU_comp_mat_sp
     module procedure copyFromGPU_comp_mat_dp
     module procedure copyFromGPU_real_sp
     module procedure copyFromGPU_real_dp
  end interface copyFromGPU

  interface deleteGPU
     module procedure deleteGPU_dns_sp
     module procedure deleteGPU_dns_dp
     module procedure deleteGPU_logic_vec
  end interface deleteGPU

  interface createAll
     module procedure createAll_sp
     module procedure createAll_dp
  end interface createAll

  interface destroyAll
     module procedure destroyAll_sp
     module procedure destroyAll_dp
  end interface destroyAll
  ! ------------------------------------------------------------------

  interface createMatGPU
     module procedure createMatGPU_sp
     module procedure createMatGPU_dp
  end interface createMatGPU

  interface deleteMatGPU
     module procedure deleteMatGPU_sp
     module procedure deleteMatGPU_dp
  end interface deleteMatGPU
  ! ------------------------------------------------------------------

  interface copy_trid_toGPU
     module procedure copy_trid_toGPU_sp
     module procedure copy_trid_toGPU_dp
  end interface copy_trid_toGPU

  interface copy_trid_toHOST
     module procedure copy_trid_toHOST_sp
     module procedure copy_trid_toHOST_dp
  end interface copy_trid_toHOST

  interface delete_vdns_fromGPU
     module procedure delete_vdns_fromGPU_sp
     module procedure delete_vdns_fromGPU_dp
  end interface delete_vdns_fromGPU

  interface copy_vdns_toGPU
     module procedure copy_vdns_toGPU_sp
     module procedure copy_vdns_toGPU_dp
  end interface copy_vdns_toGPU

  interface copy_vdns_toHOST
     module procedure copy_vdns_toHOST_sp
     module procedure copy_vdns_toHOST_dp
  end interface copy_vdns_toHOST
 



  interface create_trid_toGPU
     module procedure create_trid_toGPU_dp
  end interface create_trid_toGPU
  
  interface init_trid_toGPU
     module procedure init_trid_toGPU_dp
  end interface init_trid_toGPU
  
  interface update_trid_toGPU
     module procedure update_trid_toGPU_dp
  end interface update_trid_toGPU
  
  interface delete_trid_fromGPU
     module procedure delete_trid_fromGPU_dp
  end interface delete_trid_fromGPU
  ! ------------------------------------------------------------------

  interface copy_mat_gpu
     module procedure copy_mat_gpu_sp
     module procedure copy_mat_gpu_dp
  end interface copy_mat_gpu

  interface matmul_gpu
     module procedure matmul_gpu_sp
     module procedure matmul_gpu_dp
  end interface matmul_gpu

  interface add_gpu
     module procedure add_gpu_sp, ccadd_gpu_sp
     module procedure add_gpu_dp, ccadd_gpu_dp
  end interface add_gpu

  interface add_cublas
     module procedure add_cublas_sp
     module procedure add_cublas_dp
  end interface add_cublas

  interface sum_gpu
     module procedure sum_gpu_sp
     module procedure sum_gpu_dp
  end interface sum_gpu

  interface inverse_gpu
     module procedure inverse_gpu_sp
     module procedure inverse_gpu_dp
  end interface inverse_gpu

  interface init_mat_gpu
     module procedure init_mat_gpu_sp
     module procedure init_mat_gpu_dp
  end interface init_mat_gpu

  interface spectral_gpu
     module procedure spectral_gpu_sp
     module procedure spectral_gpu_dp
  end interface spectral_gpu

  interface dagger_gpu
     module procedure dagger_gpu_sp
     module procedure dagger_gpu_dp
  end interface dagger_gpu

  interface trace_gpu
     module procedure trace_gpu_sp
     module procedure trace_gpu_dp
  end interface trace_gpu

contains

  subroutine createGPU_sp(A)
    type(c_DNS), intent(in) :: A     
    !$acc enter data create(A,A%val)
  end subroutine createGPU_sp

  subroutine createGPU_dp(A)
    type(z_DNS), intent(in) :: A

    !$acc enter data create(A,A%val)
  end subroutine createGPU_dp
  ! ------------------------------------------------------------------

  subroutine copyToGPU_dns_sp(A)
    type(c_DNS), intent(in) :: A     
    !$acc enter data copyin(A,A%val)
  end subroutine copyToGPU_dns_sp

  subroutine copyToGPU_dns_dp(A)
    type(z_DNS), intent(in) :: A     
    !$acc enter data copyin(A,A%val)
  end subroutine copyToGPU_dns_dp

  subroutine copyToGPU_real_sp(r)
    real(sp), intent(in) :: r     
    !$acc enter data copyin(r)
  end subroutine copyToGPU_real_sp
  
  subroutine copyToGPU_real_dp(r)
    real(dp), intent(in) :: r     
    !$acc enter data copyin(r)
  end subroutine copyToGPU_real_dp
  
  subroutine copyToGPU_comp_sp(r)
    complex(sp), intent(in) :: r     
    !$acc enter data copyin(r)
  end subroutine copyToGPU_comp_sp
  
  subroutine copyToGPU_comp_dp(r)
    complex(dp), intent(in) :: r     
    !$acc enter data copyin(r)
  end subroutine copyToGPU_comp_dp
  
  subroutine copyToGPU_comp_vec_sp(V)
    complex(sp), dimension(:), intent(in) :: V     
    !$acc enter data copyin(V)
  end subroutine copyToGPU_comp_vec_sp 

  subroutine copyToGPU_comp_vec_dp(V)
    complex(dp), dimension(:), intent(in) :: V     
    !$acc enter data copyin(V)
  end subroutine copyToGPU_comp_vec_dp

  subroutine copyToGPU_comp_mat_sp(Mat)
    complex(sp), intent(in) :: Mat(:,:)      
    !$acc enter data copyin(Mat)
  end subroutine copyToGPU_comp_mat_sp

  subroutine copyToGPU_comp_mat_dp(Mat)
    complex(dp), intent(in) :: Mat(:,:)      
    !$acc enter data copyin(Mat)
  end subroutine copyToGPU_comp_mat_dp

  subroutine copyToGPU_logic_vec(V)
    logical, dimension(:), intent(in) :: V
    !$acc enter data copyin(V)
  end subroutine copyToGPU_logic_vec

  ! ------------------------------------------------------------------

  subroutine copyFromGPU_dns_sp(A)
    type(c_DNS), intent(inout) :: A     
    !$acc update host(A%val)
  end subroutine copyFromGPU_dns_sp

  subroutine copyFromGPU_dns_dp(A)
    type(z_DNS), intent(inout) :: A     
    !$acc update host(A%val)
  end subroutine copyFromGPU_dns_dp

  subroutine copyFromGPU_comp_mat_sp(Mat)
    complex(sp), intent(out) :: Mat(:,:)      
    !$acc update host(Mat)
  end subroutine copyFromGPU_comp_mat_sp

  subroutine copyFromGPU_comp_mat_dp(Mat)
    complex(dp), intent(out) :: Mat(:,:)      
    !$acc update host(Mat)
  end subroutine copyFromGPU_comp_mat_dp
  
  subroutine copyFromGPU_real_sp(r)
    real(sp), intent(inout) :: r     
    !$acc update host(r)
  end subroutine copyFromGPU_real_sp

   subroutine copyFromGPU_real_dp(r)
    real(dp), intent(inout) :: r    
    !$acc update host(r)
  end subroutine copyFromGPU_real_dp
  ! ------------------------------------------------------------------

  subroutine deleteGPU_dns_sp(A)
    type(c_DNS), intent(inout) :: A     
    !$acc exit data delete(A%val,A)
  end subroutine deleteGPU_dns_sp

  subroutine deleteGPU_dns_dp(A)
    type(z_DNS), intent(inout) :: A     
    !$acc exit data delete(A%val,A)
  end subroutine deleteGPU_dns_dp

   subroutine deleteGPU_logic_vec(V)
    logical, dimension(:), intent(inout) :: V     
    !$acc exit data delete(V)
  end subroutine deleteGPU_logic_vec
  ! ------------------------------------------------------------------

  subroutine createAll_sp(Mat, nrow, ncol)
    type(c_DNS) :: Mat 
    integer, intent(in) :: nrow, ncol

    call create(Mat,nrow,ncol)
    call createGPU(Mat)

  end subroutine createAll_sp

  subroutine createAll_dp(Mat, nrow, ncol)
    type(z_DNS) :: Mat 
    integer, intent(in) :: nrow, ncol

    call create(Mat,nrow,ncol)
    call createGPU(Mat)

  end subroutine createAll_dp
  ! ------------------------------------------------------------------

  subroutine destroyAll_sp(Mat)
    type(c_DNS) :: Mat 
    call deleteGPU(Mat)
    call destroy(Mat)
  end subroutine destroyAll_sp

  subroutine destroyAll_dp(Mat)
    type(z_DNS) :: Mat 
    call deleteGPU(Mat)
    call destroy(Mat)
  end subroutine destroyAll_dp
  ! ------------------------------------------------------------------

  subroutine createMatGPU_sp(Mat)
    complex(sp), intent(in) :: Mat(:,:)      
    !$acc enter data create(Mat)
  end subroutine createMatGPU_sp

  subroutine createMatGPU_dp(Mat)
    complex(dp), intent(in) :: Mat(:,:)      
    !$acc enter data create(Mat)
  end subroutine createMatGPU_dp
  ! ------------------------------------------------------------------

  subroutine deleteMatGPU_sp(Mat)
    complex(sp), intent(in) :: Mat(:,:)      
    !$acc exit data delete(Mat)
  end subroutine deleteMatGPU_sp

  subroutine deleteMatGPU_dp(Mat)
    complex(dp), intent(in) :: Mat(:,:)      
    !$acc exit data delete(Mat)
  end subroutine deleteMatGPU_dp
  ! ------------------------------------------------------------------

  ! C = alpha*A*B + beta*C
  subroutine matmul_gpu_sp(hcublas,alpha,A,B,beta,C,dagger)
    type(cublasHandle), intent(in) :: hcublas
    complex(sp), intent(in) :: alpha      
    complex(sp), intent(in) :: A(:,:)      
    complex(sp), intent(in) :: B(:,:)      
    complex(sp), intent(in) :: beta     
    complex(sp), intent(inout) :: C(:,:)      
    character(*), intent(in), optional :: dagger 

    integer :: istat, m,n,k
    
    if (.not.present(dagger)) then
       m = size(A,1)
       n = size(C,2)
       k = size(B,1)
       !$acc data present(A, B, C)
       !$acc host_data use_device(A, B, C)
       istat = cublasCgemm(hcublas, CUBLAS_OP_N, CUBLAS_OP_N, m, n, k, &
            & alpha, A, m, B, k, beta, C, m)
       !$acc end host_data
       !$acc end data
     endif
     if (dagger == 'dag_1st') then
       m = size(A,2)
       n = size(B,2)
       k = size(A,1)
       !$acc data present(A, B, C)
       !$acc host_data use_device(A, B, C)
       istat = cublasCgemm(hcublas, CUBLAS_OP_C, CUBLAS_OP_N, m, n, k, &
            & alpha, A, k, B, k, beta, C, m)
       !$acc end host_data
       !$acc end data
    endif
    if (dagger == 'dag_2nd') then
       m = size(A,1)
       n = size(B,1)
       k = size(A,2)
       !$acc data present(A, B, C)
       !$acc host_data use_device(A, B, C)
       istat = cublasCgemm(hcublas, CUBLAS_OP_N, CUBLAS_OP_C, m, n, k, &
            & alpha, A, m, B, n, beta, C, m)
       !$acc end host_data
       !$acc end data
    endif

  end subroutine matmul_gpu_sp
  ! ------------------------------------------------------------------

  subroutine matmul_gpu_dp(hcublas,alpha,A,B,beta,C,dagger)
    type(cublasHandle), intent(in) :: hcublas
    complex(dp), intent(in) :: alpha      
    complex(dp), intent(in) :: A(:,:)      
    complex(dp), intent(in) :: B(:,:)      
    complex(dp), intent(in) :: beta    
    complex(dp), intent(inout) :: C(:,:)      
    character(*), intent(in), optional :: dagger 

    integer :: istat, m,n,k
    
    if (.not.present(dagger)) then
       m = size(A,1)
       n = size(C,2)
       k = size(B,1)
       call checksum(hcublas,A,'A')
       call checksum(hcublas,B,'B')
       !$acc data present(A, B, C)
       !$acc host_data use_device(A, B, C)
       istat = cublasZgemm(hcublas, CUBLAS_OP_N, CUBLAS_OP_N, m, n, k, &
            & alpha, A, m, B, k, beta, C, m)
       !$acc end host_data
       !$acc end data
     endif
     if (dagger == 'dag_1st') then
       m = size(A,2)
       n = size(B,2)
       k = size(A,1)
       call checksum(hcublas,A,'A')
       call checksum(hcublas,B,'B')
       !$acc data present(A, B, C)
       !$acc host_data use_device(A, B, C)
       istat = cublasZgemm(hcublas, CUBLAS_OP_C, CUBLAS_OP_N, m, n, k, &
            & alpha, A, k, B, k, beta, C, m)
       !$acc end host_data
       !$acc end data
    endif
    if (dagger == 'dag_2nd') then
       m = size(A,1)
       n = size(B,1)
       k = size(A,2)
       call checksum(hcublas,A,'A')
       call checksum(hcublas,B,'B')
       !$acc data present(A, B, C)
       !$acc host_data use_device(A, B, C)
       istat = cublasZgemm(hcublas, CUBLAS_OP_N, CUBLAS_OP_C, m, n, k, &
            & alpha, A, m, B, n, beta, C, m)
       !$acc end host_data
       !$acc end data
    endif

  end subroutine matmul_gpu_dp
  ! ------------------------------------------------------------------

  subroutine sum_gpu_sp(hcublas, Mat, summ)
    type(cublasHandle), intent(in) :: hcublas
    complex(sp), intent(in) :: Mat(:,:)
    real(sp), intent(out) :: summ
    !
    integer :: istat, n

    n=size(Mat)

    !$acc data present(Mat)
    !$acc host_data use_device(Mat)
    summ = 0.0_sp
    istat=cublasScasum(hcublas, n, Mat, 1, summ)
    !$acc end host_data
    !$acc end data
  end subroutine sum_gpu_sp
  ! ------------------------------------------------------------------

  subroutine sum_gpu_dp(hcublas, Mat, summ)
    type(cublasHandle), intent(in) :: hcublas
    complex(dp), intent(in) :: Mat(:,:)
    real(dp), intent(out) :: summ
    !
    integer :: istat, n

    n=size(Mat)
    !$acc data present(Mat)
    !$acc host_data use_device(Mat)
    summ = 0.0_dp
    istat=cublasDzasum(hcublas, n, Mat, 1, summ)
    !$acc end host_data
    !$acc end data
  end subroutine sum_gpu_dp
  ! ------------------------------------------------------------------

  subroutine copy_mat_gpu_sp(hcublas, A, B)
    type(cublasHandle), intent(in) :: hcublas
    complex(sp), intent(in) :: A(:,:)
    complex(sp), intent(out) :: B(:,:)

    integer :: istat, n 

    if (size(A) /= size(B)) then
       stop 'ERROR in copy_mat size mismatch'
    end if

    n=size(A)
    !$acc data present(A, B)
    !$acc host_data use_device(A,B)
    istat=cublasCcopy(hcublas, n, A, 1, B, 1)
    !$acc end host_data
    !$acc end data
  end subroutine copy_mat_gpu_sp
  ! ------------------------------------------------------------------

  subroutine copy_mat_gpu_dp(hcublas, A, B)
    type(cublasHandle), intent(in) :: hcublas
    complex(dp), intent(in) :: A(:,:)
    complex(dp), intent(out) :: B(:,:)

    integer :: istat, n 

    if (size(A) /= size(B)) then
       stop 'ERROR in copy_mat size mismatch'
    end if

    n=size(A)
    !$acc data present(A, B)
    !$acc host_data use_device(A,B)
    istat=cublasZcopy(hcublas, n, A, 1, B, 1)
    !$acc end host_data
    !$acc end data
  end subroutine copy_mat_gpu_dp
  ! ------------------------------------------------------------------

  subroutine inverse_gpu_sp(hcublas, hcusolver, A, Ainv, err)
    type(cublasHandle), intent(in) :: hcublas
    type(cusolverDnHandle), intent(in) :: hcusolver
    complex(sp), intent(in) :: A(:,:) 
    complex(sp), intent(inout) :: Ainv(:,:) 
    integer, intent(out) :: err

    complex(sp), allocatable :: LU(:,:) 
    complex(sp), allocatable, dimension(:) :: work 
    integer, allocatable, dimension(:) :: pivot
    integer :: n, ii, jj, lwork, istat, err1, err2
    !$acc declare present(A, Ainv)

    n= size(A,1)

    !$acc data present(A)
    istat=cusolverDnCgetrf_bufferSize(hcusolver,n,n,A,n,lwork)
    !$acc end data
    allocate(work(lwork))
    allocate(LU(n,n))
    allocate(pivot(n))

    !$acc kernels present(Ainv(n,n))
    Ainv = (0.0, 0.0)
    do jj = 1, n
       Ainv(jj,jj) = (1.0, 0.0)
    end do
    !$acc end kernels

    !$acc data create(LU, work, pivot, err1, err2) 
    !$acc host_data use_device(A, Ainv, LU, work, pivot, err1, err2)
    istat=cublasCcopy(hcublas, n*n, A, 1, LU, 1)
    istat=cusolverDnCgetrf(hcusolver,n,n,LU,n,work,pivot,err1)
    istat=cusolverDnCgetrs(hcusolver,CUBLAS_OP_N,n,n,LU,n,pivot,Ainv,n,err2)
    !$acc end host_data
    !$acc end data         
    err = err1

  end subroutine inverse_gpu_sp
  ! ------------------------------------------------------------------

  subroutine inverse_gpu_dp(hcublas, hcusolver, A, Ainv, err)
    type(cublasHandle), intent(in) :: hcublas
    type(cusolverDnHandle), intent(in) :: hcusolver
    complex(dp), intent(in) :: A(:,:) 
    complex(dp), intent(inout) :: Ainv(:,:) 
    integer, intent(out) :: err

    complex(dp), allocatable :: LU(:,:) 
    complex(dp), allocatable, dimension(:) :: work 
    integer, allocatable, dimension(:) :: pivot
    integer :: n, ii, jj, lwork, istat, err1, err2
    !$acc declare present(A, Ainv)

    n= size(A,1)

    !$acc data present(A)
    istat=cusolverDnZgetrf_bufferSize(hcusolver,n,n,A,n,lwork)
    !$acc end data
    allocate(work(lwork))
    allocate(LU(n,n))
    allocate(pivot(n))

    !$acc kernels present(Ainv(n,n))
    Ainv = (0.0, 0.0)
    do jj = 1, n
       Ainv(jj,jj) = (1.0, 0.0)
    end do
    !$acc end kernels

    !$acc data create(LU, work, pivot, err1, err2) 
    !$acc host_data use_device(A, Ainv, LU, work, pivot, err1, err2)
    istat=cublasZcopy(hcublas, n*n, A, 1, LU, 1)
    istat=cusolverDnZgetrf(hcusolver,n,n,LU,n,work,pivot,err1)
    istat=cusolverDnZgetrs(hcusolver,CUBLAS_OP_N,n,n,LU,n,pivot,Ainv,n,err2)
    !$acc end host_data
    !$acc end data         
    err = err1

  end subroutine inverse_gpu_dp
  ! ------------------------------------------------------------------

  subroutine init_mat_gpu_sp(A)
    type(c_DNS ) :: A      

    integer :: jj

    !$acc kernels present(A)
    A%val = (0.0, 0.0)
    do jj = 1, size(A%val,1) 
       A%val(jj,jj) = (1.0, 0.0)
    end do
    !$acc end kernels

  end subroutine init_mat_gpu_sp
  ! ------------------------------------------------------------------

  subroutine init_mat_gpu_dp(A)
    complex(dp), dimension(:,:) :: A      

    integer :: jj

    !$acc kernels present(A)
    A = (0.0, 0.0)
    do jj = 1, size(A,1) 
       A(jj,jj) = (1.0, 0.0)
    end do
    !$acc end kernels

  end subroutine init_mat_gpu_dp
  ! ------------------------------------------------------------------

  subroutine add_gpu_sp(A, B, C)
    type(c_DNS ), intent(in) :: A      
    type(c_DNS ), intent(in) :: B      
    type(c_DNS ), intent(inout) :: C      

    integer :: ii, jj, nrow, ncol
    nrow = size(A%val,1)
    ncol = size(A%val,2)

    !$acc kernels present(A, B, C)
    do jj = 1, ncol 
       do ii = 1, nrow 
          C%val(ii,jj) = A%val(ii,jj) + B%val(ii,jj)
       end do
    end do
    !$acc end kernels

  end subroutine add_gpu_sp
  ! ------------------------------------------------------------------

  subroutine add_gpu_dp(A, B, C)
    type(z_DNS), intent(in) :: A      
    type(z_DNS), intent(in) :: B      
    type(z_DNS), intent(inout) :: C      

    integer :: ii, jj, nrow, ncol
    nrow = size(A%val,1)
    ncol = size(A%val,2)

    !$acc kernels present(A, B, C)
    do jj = 1, ncol 
       do ii = 1, nrow 
          C%val(ii,jj) = A%val(ii,jj) + B%val(ii,jj)
       end do
    end do
    !$acc end kernels
  end subroutine add_gpu_dp
  ! ------------------------------------------------------------------

  subroutine ccadd_gpu_sp(alpha, A, beta,B, C)
    !C = alpha*A + beta*B      
    complex(sp), dimension(:,:), intent(in) :: A      
    complex(sp), dimension(:,:), intent(in) :: B 
    complex(sp) :: alpha, beta
    complex(sp), dimension(:,:), intent(inout) :: C      

    integer :: ii, jj, nrow, ncol
    nrow = size(A,1)
    ncol = size(A,2)

    !$acc kernels present(A, B, C)
    do jj = 1, ncol 
       do ii = 1, nrow 
          C(ii,jj) = alpha*A(ii,jj) + beta*B(ii,jj)
       end do
    end do
    !$acc end kernels
  end subroutine ccadd_gpu_sp
  ! ------------------------------------------------------------------

  subroutine ccadd_gpu_dp(alpha, A, beta,B, C)
    !C = alpha*A + beta*B      
    complex(dp), dimension(:,:), intent(in) :: A      
    complex(dp), dimension(:,:), intent(in) :: B 
    complex(dp) :: alpha, beta
    complex(dp), dimension(:,:), intent(inout) :: C      

    integer :: ii, jj, nrow, ncol
    nrow = size(A,1)
    ncol = size(A,2)

    !$acc kernels present(A, B, C)
    do jj = 1, ncol 
       do ii = 1, nrow 
          C(ii,jj) = alpha*A(ii,jj) + beta*B(ii,jj)
       end do
    end do
    !$acc end kernels
  end subroutine ccadd_gpu_dp
  ! ------------------------------------------------------------------


  subroutine spectral_gpu_sp(Gr, A)
    complex(sp), dimension(:,:), intent(in) :: Gr            
    complex(sp), dimension(:,:), intent(inout) :: A      

    complex(sp) :: ii
    integer :: kk, jj, nrow, ncol
    nrow = size(A,1)
    ncol = size(A,2)

    ii = (0.0_sp,1.0_sp) 
    !$acc kernels present(Gr, A)
    do jj = 1, ncol 
       do kk = 1, nrow 
          A(kk,jj) = ii*(Gr(kk,jj) - conjg(Gr(jj,kk)))
       end do
    end do
    !$acc end kernels
  end subroutine spectral_gpu_sp
  ! ------------------------------------------------------------------

  subroutine spectral_gpu_dp(Gr, A)
    complex(dp), dimension(:,:), intent(in) :: Gr            
    complex(dp), dimension(:,:), intent(inout) :: A      

    complex(dp) :: ii
    integer :: kk, jj, nrow, ncol
    nrow = size(A,1)
    ncol = size(A,2)

    ii = (0.0_dp,1.0_dp) 
    !$acc kernels present(Gr, A)
    do jj = 1, ncol 
       do kk = 1, nrow 
          A(kk,jj) = ii*(Gr(kk,jj) - conjg(Gr(jj,kk)))
       end do
    end do
    !$acc end kernels
  end subroutine spectral_gpu_dp
  ! ------------------------------------------------------------------

  subroutine dagger_gpu_sp(Gr, A)
    complex(sp), dimension(:,:), intent(in) :: Gr            
    complex(sp), dimension(:,:), intent(inout) :: A      

    integer :: kk, jj, nrow, ncol
    nrow = size(A,1)
    ncol = size(A,2)

    !$acc kernels present(Gr, A)
    do jj = 1, ncol 
       do kk = 1, nrow 
          A(kk,jj) = conjg(Gr(jj,kk))
       end do
    end do
    !$acc end kernels
  end subroutine dagger_gpu_sp
  ! ------------------------------------------------------------------

  subroutine dagger_gpu_dp(Gr, A)
    complex(dp), dimension(:,:), intent(in) :: Gr            
    complex(dp), dimension(:,:), intent(inout) :: A      

    integer :: kk, jj, nrow, ncol
    nrow = size(A,1)
    ncol = size(A,2)

    !$acc kernels present(Gr, A)
    do jj = 1, ncol 
       do kk = 1, nrow 
          A(kk,jj) = conjg(Gr(jj,kk))
       end do
    end do
    !$acc end kernels
  end subroutine dagger_gpu_dp
  ! ------------------------------------------------------------------

  subroutine add_cublas_sp(hcublas, C, A, alpha, B)
    !C = A + alpha B
    type(cublasHandle), intent(in) :: hcublas
    complex(sp), intent(in) :: A(:,:)      
    complex(sp), intent(in) :: B(:,:)      
    complex(sp), intent(inout) :: C(:,:)
    complex(sp), intent(in) :: alpha
    integer :: istat, n
    !complex(sp), parameter :: one = (1.0_sp,0.0_sp)
    if (size(A) /= size(B) .or. size(A) /= size(C)) then
       stop 'ERROR in add_cublas size mismatch'
    end if

    n = size(A)

    !$acc data present(A, B, C)
    !$acc host_data use_device(A, B, C)
    istat=cublasCcopy(hcublas, n, A, 1, C, 1)
    istat=cublasCaxpy(hcublas, n, alpha, B, 1, C, 1)
    !$acc end host_data
    !$acc end data

  end subroutine add_cublas_sp
  ! ------------------------------------------------------------------

  subroutine add_cublas_dp(hcublas, C, A, alpha, B)
    !C = A + alpha B
    type(cublasHandle), intent(in) :: hcublas
    complex(dp), intent(in) :: A(:,:)      
    complex(dp), intent(in) :: B(:,:)      
    complex(dp), intent(inout) :: C(:,:)
    complex(dp), intent(in) :: alpha
    integer :: istat, n
    !complex(sp), parameter :: one = (1.0_sp,0.0_sp)
    if (size(A) /= size(B) .or. size(A) /= size(C)) then
       stop 'ERROR in add_cublas size mismatch'
    end if

    n = size(A)

    !$acc data present(A, B, C)
    !$acc host_data use_device(A, B, C)
    istat=cublasZcopy(hcublas, n, A, 1, C, 1)
    istat=cublasZaxpy(hcublas, n, alpha, B, 1, C, 1)
    !$acc end host_data
    !$acc end data
  end subroutine add_cublas_dp
  ! ------------------------------------------------------------------

  subroutine copy_trid_toGPU_sp(M)
    type(c_DNS), dimension(:,:), intent(in) :: M
    integer :: ii, nbl

    nbl = size(M,1)
    call copyToGPU(M(1,1))
    do ii=2,nbl
       call copyToGPU(M(ii,ii))
       call copyToGPU(M(ii-1,ii))
       call copyToGPU(M(ii,ii-1))
    end do
  end subroutine copy_trid_toGPU_sp
  ! ------------------------------------------------------------------

  subroutine copy_trid_toGPU_dp(M)
    type(z_DNS), dimension(:,:), intent(in) :: M
    integer :: ii, nbl

    nbl = size(M,1)
    do ii=2,nbl
       call copyToGPU(M(ii,ii-1))
       call copyToGPU(M(ii-1,ii))
       call copyToGPU(M(ii,ii))
    end do
    call copyToGPU(M(1,1))
  end subroutine copy_trid_toGPU_dp
  ! ------------------------------------------------------------------

  subroutine delete_trid_fromGPU_dp(M)
    type(z_DNS), dimension(:,:), intent(in) :: M
    integer :: ii, nbl

    nbl = size(M,1)
    call deleteGPU(M(1,1))
    do ii=2,nbl
       call deleteGPU(M(ii,ii))
       call deleteGPU(M(ii-1,ii))
       call deleteGPU(M(ii,ii-1))
    end do
  end subroutine delete_trid_fromGPU_dp
  ! ------------------------------------------------------------------

  subroutine copy_trid_toHOST_sp(M)
    type(c_DNS), dimension(:,:), intent(inout) :: M
    integer :: ii, nbl

    nbl = size(M,1)
    call copyFromGPU(M(1,1))
    do ii=2,nbl
       call copyFromGPU(M(ii,ii))
       call copyFromGPU(M(ii-1,ii))
       call copyFromGPU(M(ii,ii-1))
    end do
  end subroutine copy_trid_toHOST_sp
  ! ------------------------------------------------------------------

  subroutine copy_trid_toHOST_dp(M)
    type(z_DNS), dimension(:,:), intent(inout) :: M
    integer :: ii, nbl

    nbl = size(M,1)
    call copyFromGPU(M(1,1))
    do ii=2,nbl
       call copyFromGPU(M(ii,ii))
       call copyFromGPU(M(ii-1,ii))
       call copyFromGPU(M(ii,ii-1))
    end do
  end subroutine copy_trid_toHOST_dp
  ! ------------------------------------------------------------------
  
  subroutine create_trid_toGPU_dp(nbl)
    integer, intent(in) :: nbl

    type(z_DNS), dimension(nbl,nbl) :: ESH
    integer :: ii

    call createGPU(ESH(1,1))
    do ii=2,nbl
       call createGPU(ESH(ii,ii))
       call createGPU(ESH(ii-1,ii))
       call createGPU(ESH(ii,ii-1))
    end do
  end subroutine create_trid_toGPU_dp
  ! ------------------------------------------------------------------
  
  subroutine init_trid_toGPU_dp(M)
    type(z_DNS), dimension(:,:), intent(inout) :: M
    integer :: ii, nbl

    nbl = size(M,1)
    call init_mat_gpu(M(1,1)%val)
    do ii=2,nbl
       call init_mat_gpu(M(ii,ii)%val)
       call init_mat_gpu(M(ii-1,ii)%val)
       call init_mat_gpu(M(ii,ii-1)%val)
    end do
  end subroutine init_trid_toGPU_dp
  ! ------------------------------------------------------------------
  
  subroutine update_trid_toGPU_dp(M)
    type(z_DNS), dimension(:,:), intent(inout) :: M
    integer :: ii, nbl

    nbl = size(M,1)
    !$acc update device(M(1,1)%val) 
    do ii=2,nbl
       !$acc update device(M(ii,ii)%val) 
       !$acc update device(M(ii-1,ii)%val) 
       !$acc update device(M(ii,ii-1)%val) 
    end do
  end subroutine update_trid_toGPU_dp
  ! ------------------------------------------------------------------

  subroutine copy_vdns_toGPU_sp(V)
    type(c_DNS), dimension(:), intent(in) :: V
    integer :: ii, nbl

    nbl = size(V)
    do ii=1,nbl
       if(allocated(V(ii)%val)) then
          call copyToGPU(V(ii))
       endif
    end do
  end subroutine copy_vdns_toGPU_sp
  ! ------------------------------------------------------------------

  subroutine copy_vdns_toGPU_dp(V)
    type(z_DNS), dimension(:), intent(in) :: V
    integer :: ii, nbl

    nbl = size(V)
    do ii=1,nbl
       if(allocated(V(ii)%val)) then
          call copyToGPU(V(ii))
       endif   
    end do
  end subroutine copy_vdns_toGPU_dp
  ! ------------------------------------------------------------------

  subroutine copy_vdns_toHOST_sp(V)
    type(c_DNS), dimension(:), intent(in) :: V
    integer :: ii, nbl

    nbl = size(V)
    do ii=1,nbl
       call copyFromGPU(V(ii))
    end do
  end subroutine copy_vdns_toHOST_sp
  ! ------------------------------------------------------------------

  subroutine copy_vdns_toHOST_dp(V)
    type(z_DNS), dimension(:), intent(in) :: V
    integer :: ii, nbl

    nbl = size(V)
    do ii=1,nbl
       call copyFromGPU(V(ii))
    end do
  end subroutine copy_vdns_toHOST_dp
  ! ------------------------------------------------------------------

  subroutine delete_vdns_fromGPU_sp(V)
    type(c_DNS), dimension(:) :: V
    integer :: ii, nbl

    nbl = size(V)
    do ii=1,nbl
       if(allocated(V(ii)%val)) then
          call deleteGPU(V(ii))
       endif
    end do
  end subroutine delete_vdns_fromGPU_sp
  ! ------------------------------------------------------------------

  subroutine delete_vdns_fromGPU_dp(V)
    type(z_DNS), dimension(:) :: V
    integer :: ii, nbl

    nbl = size(V)
    do ii=1,nbl
       if(allocated(V(ii)%val)) then
          call deleteGPU(V(ii))
       endif
    end do
  end subroutine delete_vdns_fromGPU_dp
  ! ------------------------------------------------------------------

  function trace_gpu_sp(mat, mask) result(trace)
    complex(sp), dimension(:,:), intent(in) :: mat
    logical, intent(in), optional :: mask(:)
    complex(sp) :: trace

    integer :: i, nrow

    nrow = size(mat,1)
    trace = (0_sp,0_sp)
    if (present(mask)) then
       if (size(mask) /= nrow) then
          stop 'Error in ztrace_csr: size(mask) /= nrow'
       end if
       !$acc kernels present(mat, mask) copy(trace)
       do i = 1,nrow
          if (mask(i)) then
             trace = trace + mat(i,i)
          end if
       end do
       !$acc end kernels
    else
       !$acc kernels present(mat) copy(trace)
       do i = 1,nrow
          trace = trace + mat(i,i)
       end do
       !$acc end kernels
    end if
  end function trace_gpu_sp
  ! ------------------------------------------------------------------

  function trace_gpu_dp(mat, mask) result(trace)
    complex(dp), dimension(:,:), intent(in) :: mat
    logical, intent(in), optional :: mask(:)
    complex(dp) :: trace

    integer :: i, nrow

    nrow = size(mat,1)
    trace = (0_dp,0_dp)
    if (present(mask)) then
       if (size(mask) /= nrow) then
          stop 'Error in ztrace_csr: size(mask) /= nrow'
       end if
       !$acc kernels present(mat, mask) copy(trace)
       do i = 1,nrow
          if (mask(i)) then
             trace = trace + mat(i,i)
          end if
       end do
       !$acc end kernels
    else
       !$acc kernels present(mat) copy(trace)
       do i = 1,nrow
          trace = trace + mat(i,i)
       end do
       !$acc end kernels
    end if
  end function trace_gpu_dp

  ! ----------------------------------------------------------------
  subroutine checksum(hcublas, A, nome)
    type(cublasHandle), intent(in) :: hcublas
    complex(dp), intent(in) :: A(:,:)      
    character(*), intent(in) :: nome

    real(dp) :: summ

    call sum_gpu(hcublas, A, summ)

    if (ieee_is_nan(summ)) then
       write(*,*) 'GPU:   ',trim(nome),summ
    end if

  end subroutine checksum

end module cudautils
