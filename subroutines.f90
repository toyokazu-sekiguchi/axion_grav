#include "mkl_vsl.fi"

module subroutines
#ifdef _OPENMP
    use omp_lib
#endif
  use mkl_vsl_type
  use mkl_vsl
  use utils
  implicit none
  
  integer n_streams ! should be the same as num_threads in settings.f90
  type(vsl_stream_state) stream0
  type(vsl_stream_state),allocatable :: streams(:)
!8!  integer,parameter :: rng=VSL_BRNG_MCG59 ! type of random number generator
  integer(8),parameter :: rng=VSL_BRNG_MCG59 ! type of random number generator

  integer(4) seed(2)
  
contains

  subroutine check_syminv(n,mat,imat,tol)
    integer,intent(in) :: n
    double precision,intent(in) :: mat(n*(n+1)/2),imat(n*(n+1)/2)
    double precision,intent(in) :: tol
    double precision d,dum
    integer i,j,k,ij1,ij2
    
    d=0.d0
    do i=1,n
       do j=1,n
          dum=0.d0
          do k=1,n
             ij1=max(i,k)
             ij1=ij1*(ij1-1)/2+min(i,k)
             ij2=max(j,k)
             ij2=ij2*(ij2-1)/2+min(j,k)
             dum=dum+mat(ij1)*imat(ij2)
          end do
          if(i==j) dum=dum-1
          d=max(d,abs(dum))
       end do
    end do
    
    if(d>tol) write(*,'(A,E,A,E)') 'WARNING: residue in check_syminv is significantly large: res=',d,'> tol=',tol
    
  end subroutine check_syminv
  
  subroutine init_rand_gen(n,str)
#ifdef MPI
    include "mpif.h"
#endif
    integer,intent(in) :: n ! number of threads
    character(*),intent(in) :: str
    integer info,i
    integer(8) m
    integer(8),parameter :: m0=2147483648 ! 2^31
    
    n_streams=n
    
    if(str/='') then ! seed is given
       read(str,*,iostat=info) seed(1:2)
       if(info/=0) call MpiStop('error: cannot read seed; we need two integer*4')
    else
       if(MpiRank==0) then
          call system_clock(count=m)
          m=mod(m,m0*m0) ! m should be less than m0*m0 to get two integer*4
          seed(2)=int(m/m0,kind=4)
          seed(1)=int(m-m0*seed(2),kind=4)
       end if
#ifdef MPI
       call MPI_BCAST(seed(1:2),2,MPI_INTEGER8,0,MPI_COMM_WORLD,info)
#endif
    end if
    
    info=VSLNewStreamEx(stream0,rng,2,seed)
    if(info/=VSL_STATUS_OK) call MpiStop('error: something wrong with VSLNewStreamEx')
    
    if(MpiRank==0) then
       write(*,'(" random seed = ",2I)',iostat=info) seed(1:2)
       if(str/='') then
          write(*,'(A)') '  random seed is given'
       else
          write(*,'(A)') '  random seed is taken from system clock'
       end if
    end if
    
    allocate(streams(n_streams),stat=info)
    if(info/=0) call MpiStop('error: cannot allocate streams')
    do i=1,n_streams
       info=VSLCopyStream(streams(i),stream0)
       if(info/=VSL_STATUS_OK) call MpiStop('error: something wrong with VSLCopyStream')
    end do
    
  end subroutine init_rand_gen

  subroutine fin_rand_gen
    integer info,i
    
    do i=1,n_streams
       info=VSLDeleteStream(streams(i))
       if(info/=VSL_STATUS_OK) call MpiStop('error: something wrong with VSLDeleteStream')
    end do
    deallocate(streams,stat=info)
    
  end subroutine fin_rand_gen

  subroutine skip_rand_stream0(n,mult)
    integer(8),intent(in) :: n
    integer,intent(in) :: mult
    integer info
    
    info=VSLSkipAheadStream(stream0,n*mult)
    if(info/=0) call MpiStop('error: something wrong in VSLSkipAheadStream')
    
  end subroutine skip_rand_stream0
  
  subroutine skip_rand_stream(i,n,mult)
    integer,intent(in) :: i,mult
    integer(8),intent(in) :: n
    integer info
    
    if(.not.(i>=1.and.i<=n_streams)) call MpiStop('error: wrong i_stream')
    info=VSLSkipAheadStream(streams(i),n*mult)
    if(info/=0) call MpiStop('error: something wrong in VSLSkipAheadStream')
    
  end subroutine skip_rand_stream
  
  subroutine put_rand_stream0
    integer info,i
    
    do i=1,n_streams
       info=VSLCopyStreamState(streams(i),stream0)
       if(info/=0) call MpiStop('error: something wrong in VSLCopyStreamState')
    end do
    
  end subroutine put_rand_stream0
  
  subroutine gen_rand_gauss(i,x,n)
    integer,intent(in) :: i,n
    double precision,intent(out) :: x(n)
    integer,parameter :: meth=VSL_RNG_METHOD_GAUSSIAN_BOXMULLER
    integer info
    
    info=VDRngGaussian(meth,streams(i),n,x,0d0,1d0)
    if(info/=0) call MpiStop('error: something wrong in VDRngGaussian')
    
  end subroutine gen_rand_gauss
  
  subroutine gen_rand_uniform(i,x,n)
    integer,intent(in) :: i,n
    double precision,intent(out) :: x(n)
    integer,parameter :: meth=VSL_RNG_METHOD_UNIFORM_STD ! this suffices
    integer info
    
    info=VDRngUniform(meth,streams(i),n,x,0d0,1d0)
    if(info/=0) call MpiStop('error: something wrong in VDRngUniform')
    
  end subroutine gen_rand_uniform
  
  function get_openmp_nskip(i,n,niter)
    ! compute the number of skipped iterations when SCHEDULE(STATIC) is specified in OMP DO
    integer,intent(in) :: i,n ! i is thread num, n is num_threads, niter is number of iterations
    integer(8),intent(in) :: niter
    integer(8) get_openmp_nskip
    
    get_openmp_nskip=niter/n
    get_openmp_nskip=get_openmp_nskip*(i-1)
    get_openmp_nskip=get_openmp_nskip+min(i-1,mod(niter,n))
    
  end function get_openmp_nskip
  
end module subroutines
