module spectrum
  use settings
  use potential,only : pi,twopi,fourpi
  implicit none
  
  integer,parameter :: num_largeq=500
  integer num_cov
  
  double precision,allocatable :: mat_spec(:,:),largeq_spec(:,:),dq_spec(:),wav_spec(:), &
       power_spec(:,:),power_cov(:,:),mat_deconv(:),conv_spec(:),num_conv_spec(:),crho_spec(:,:)
  
  integer(8) ngrid_spec,ngrid_tot_spec,ngo2_spec,alloc_local_spec,local_0_spec,local_0_offset_spec
  type(C_PTR) plan_forward_spec,plan_backward_spec,cdata_spec1,cdata_spec2
#if DIM==2
  complex(C_DOUBLE_COMPLEX),pointer :: data_spec1(:,:),data_spec2(:,:)
#elif DIM==3
  complex(C_DOUBLE_COMPLEX),pointer :: data_spec1(:,:,:),data_spec2(:,:,:)
#endif
  double precision dspa_spec,dwav_spec,dv_spa_spec,dv_wav_spec
  
contains
  
  subroutine init_spectrum
    integer i,j,iq,info,ij,deg
    double precision q,maxq,minq,dq,li,si,lj,sj,u,v
    integer,parameter :: num_int=200
    
    num_cov=num_bin*(num_bin+1)/2
    call MPI_BARRIER(MPI_COMM_WORLD,info)
    if(MpiRank==0) then
       write(*,*)
       write(*,'(A)') ' #precomputing matrix for spectrum estimation ...'
       write(*,'(A,I)') '  number of bins in spectrum:',num_bin
       write(*,'(A,I)') '  number of bins in cov mat:',num_cov
    end if
    
    allocate(mat_spec(0:num_largeq,num_cov),largeq_spec(0:num_largeq,num_cov),dq_spec(num_cov), &
         wav_spec(0:num_bin),power_spec(num_bin,n_chk-n_id+1:n_chk),power_cov(num_cov,n_chk-n_id+1:n_chk), &
         crho_spec(2,n_chk-n_id+1:n_chk),stat=info)
    if(info/=0) stop 'error: cannot allocate mat_spec, ..., power_spec'
    
    deg=2**degrade_spec
    ngrid_spec=ngrid/deg
    ngrid_tot_spec=ngrid_spec**ndim
    ngo2_spec=ngrid_spec/2
    dspa_spec=dspa*deg
    dwav_spec=twopi/dspa_spec/ngrid_spec
    dv_spa_spec=dspa_spec**ndim
    dv_wav_spec=(dwav_spec/twopi)**ndim
    
    if(MpiRank==0) then
       write(*,'(A,I)') '  degrade_spec =',degrade_spec
       write(*,'(A,I)') '  ngrid_spec =',ngrid_spec
       write(*,'(A,I)') '  ngrid_tot_spec =',ngrid_tot_spec
       write(*,'(A,E)') '  dspa_spec =',dspa_spec
       write(*,'(A,E)') '  dwav_spec =',dwav_spec
       write(*,'(A,E)') '  dv_spa_spec =',dv_spa_spec
       write(*,'(A,E)') '  dv_wav_spec =',dv_wav_spec
    end if
    
    call init_fftw_spec
    
    q=pi/dspa_spec ! maximum q
    dq=q/num_bin
    if(MpiRank==0) then
       write(*,'(A,E)') '  bin size of wav_spec (delta_k*tau_cr/4pi): ',dq/fourpi
       write(*,'(A,E)') '  bin size of wav_spec (delta_k*tau_fin/4pi): ',dq/fourpi*tau_fin
    end if
    wav_spec(0)=0d0
    do iq=1,num_bin
       wav_spec(iq)=wav_spec(iq-1)+dq
    end do
    
    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(j,lj,sj,i,li,si,ij,maxq,minq,dq,iq,q) SCHEDULE(DYNAMIC,1)
    do ij=1,num_cov
       i=nint(sqrt(2.d0*ij))
       j=ij-i*(i-1)/2
       
       lj=wav_spec(j)
       sj=wav_spec(j-1)
       li=wav_spec(i)
       si=wav_spec(i-1)
       
       maxq=li+lj
       minq=max(si,sj)-min(li,lj)
       if(minq<0) minq=0.d0
       dq=maxq-minq
       dq=dq/num_largeq
       
       dq_spec(ij)=dq
       do iq=0,num_largeq
          q=minq+dq*iq
          largeq_spec(iq,ij)=q
          call integrate_p_check(num_int,li,si,lj,sj,q,mat_spec(iq,ij))
       end do
       
    end do
    !$OMP END PARALLEL DO
    
    call MPI_BARRIER(MPI_COMM_WORLD,info)
    if(MpiRank==0) write(*,'(A)') ' *... done'
    
  end subroutine init_spectrum
  
  subroutine fin_spectrum
    integer info
    
    call get_epsilon
    
    call fin_fftw_spec
    
    deallocate(mat_spec,largeq_spec,dq_spec,wav_spec, &
         power_spec,power_cov,crho_spec,stat=info)
    
  end subroutine fin_spectrum
  
  subroutine integrate_p_check(num_p_int,li,si,lj,sj,q,v)
    integer,intent(in) :: num_p_int
    double precision,intent(in) :: li,si,lj,sj,q
    double precision,intent(out) :: v
    double precision p,maxp,minp,dp,x,maxmu,minmu,mu1,mu2,mutol
#if DIM==2
    double precision maxth,minth
#endif
    integer ip
    double precision,parameter :: tol=1.d-5
    
    maxp=min(lj,q+li)
    minp=max(q-li,si-q)
    minp=max(minp,sj)
    dp=(maxp-minp)/num_p_int
    
    mutol=1.d0+tol
    
    x=0.d0
    if(dp>0) then
       
       do ip=1,num_p_int
          p=dp*(ip-0.5d0)+minp
          
          if(q==0.d0) then
             maxmu=1.d0
             minmu=-1.d0
          else
             maxmu=(li*li-p*p-q*q)/2/p/q
             minmu=(si*si-p*p-q*q)/2/p/q
             
             if(maxmu<-mutol.or.minmu>mutol) &
                  call MpiStop('error: something is wrong in integrate_p_check') 
             
             maxmu=min(1.d0,maxmu)
             minmu=max(-1.d0,minmu)
          end if
          
#if DIM==2
          x=x+p*dp*(acos(maxmu)-acos(minmu))
#elif DIM==3
          x=x+p*p*dp*(maxmu-minmu)
#endif
          
       end do
       
    end if
    
#if DIM==2
    v=x*2
#elif DIM==3
    v=x*twopi
#endif
    
  end subroutine integrate_p_check
  
  subroutine init_fftw_spec
#include <fftw3-mpi.f03>
    integer info
    
#if DIM==2
    alloc_local_spec=FFTW_MPI_LOCAL_SIZE_2D(ngrid_spec, &
         ngrid_spec,MPI_COMM_WORLD,local_0_spec,local_0_offset_spec)
#elif DIM==3
    alloc_local_spec=FFTW_MPI_LOCAL_SIZE_3D(ngrid_spec,ngrid_spec, &
         ngrid_spec,MPI_COMM_WORLD,local_0_spec,local_0_offset_spec)
#endif
    write(*,'(A,5I)') ' information of data_spec:',MpiRank,ngrid_spec,local_0_spec,local_0_offset_spec,alloc_local_spec
    
    cdata_spec1=FFTW_ALLOC_COMPLEX(alloc_local_spec)
    cdata_spec2=FFTW_ALLOC_COMPLEX(alloc_local_spec)
    if(.not.c_associated(cdata_spec1)) call MpiStop('error: cdata_spec1 is null')
    if(.not.c_associated(cdata_spec2)) call MpiStop('error: cdata_spec2 is null')
    
#if DIM==2
    call c_f_pointer(cdata_spec1,data_spec1,[ngrid_spec,local_0_spec])
    call c_f_pointer(cdata_spec2,data_spec2,[ngrid_spec,local_0_spec])
#elif DIM==3
    call c_f_pointer(cdata_spec1,data_spec1,[ngrid_spec,ngrid_spec,local_0_spec])
    call c_f_pointer(cdata_spec2,data_spec2,[ngrid_spec,ngrid_spec,local_0_spec])
#endif
    if(.not.(associated(data_spec1).and.associated(data_spec2))) &
         call MpiStop('error: data_spec1 or data_spec2 is not allocated')
    
    call FFTW_PLAN_WITH_NTHREADS(int(num_threads,kind=4))
#if DIM==2
    plan_forward_spec=FFTW_MPI_PLAN_DFT_2D(ngrid_spec,ngrid_spec, &
         data_spec1,data_spec2,MPI_COMM_WORLD,FFTW_FORWARD,FFTW_ESTIMATE)
    plan_backward_spec=FFTW_MPI_PLAN_DFT_2D(ngrid_spec,ngrid_spec, &
         data_spec2,data_spec1,MPI_COMM_WORLD,FFTW_BACKWARD,FFTW_ESTIMATE)
#elif DIM==3
    plan_forward_spec=FFTW_MPI_PLAN_DFT_3D(ngrid_spec,ngrid_spec,ngrid_spec, &
         data_spec1,data_spec2,MPI_COMM_WORLD,FFTW_FORWARD,FFTW_ESTIMATE)
    plan_backward_spec=FFTW_MPI_PLAN_DFT_3D(ngrid_spec,ngrid_spec,ngrid_spec, &
         data_spec2,data_spec1,MPI_COMM_WORLD,FFTW_BACKWARD,FFTW_ESTIMATE)
#endif
    if(.not.c_associated(plan_forward_spec)) call MpiStop('error: plan_forward_spec is null')
    if(.not.c_associated(plan_backward_spec)) call MpiStop('error: plan_backward_spec is null')
    
    call memory_check_mkl()
    
  end subroutine init_fftw_spec
  
  subroutine fin_fftw_spec
#include <fftw3-mpi.f03>
    
    call FFTW_DESTROY_PLAN(plan_forward_spec)
    call FFTW_DESTROY_PLAN(plan_backward_spec)
    
    data_spec1=>NULL()
    data_spec2=>NULL()
    call FFTW_FREE(cdata_spec1)
    call FFTW_FREE(cdata_spec2)
    
  end subroutine fin_fftw_spec
  
  subroutine get_comwav_spec(ixyz,q) ! should be modified to include degrading?
    integer(8),intent(in) :: ixyz(ndim)
    double precision,intent(out) :: q
    
    q=sum(dble((mod(ixyz(1:ndim)+ngo2_spec-1,ngrid_spec)-ngo2_spec)**2))
    q=sqrt(q)*dwav_spec
    
  end subroutine get_comwav_spec
  
  subroutine spectrum_estimation(idx)
    integer,intent(in) :: idx
    integer,parameter :: max_mc=100,min_mc=10
    double precision,parameter :: tol=1.d-1
    integer info,i,j,ij
    double precision dq,t2,rho,drho
    character(len=charlen) str,file_spec,file_cov
    
    if(MpiRank==0) write(*,'(A)') ' *estimating spectrum ...'
    
#if DIM==2
    call MpiStop('NOT SURE t2=tau is appropriate')
    t2=tau_id(idx) ! this factor counts for tau^3 from expansion and tau^-2 from derivative respect to conformal time instead of proper one
#elif DIM==3
    t2=tau_id(idx)**2 ! this factor counts for tau^4 from expansion and tau^-2 from derivative respect to conformal time instead of proper one
#endif
    
    allocate(mat_deconv(num_cov),stat=info)
    if(info/=0) call MpiStop('error: cannot allocate mat_deconv')
    call compute_matrix_deconv(idx)
    
    allocate(conv_spec(num_bin),num_conv_spec(num_bin),stat=info)
    if(info/=0) call MpiStop('error: cannot allocate conv_spec,num_conv_spec')
    call compute_conv_spec(idx)
    
    do i=1,num_bin
       conv_spec(i)=conv_spec(i)*t2
       do j=1,i
          ij=i*(i-1)/2+j
          mat_deconv(ij)=mat_deconv(ij)/num_conv_spec(i)/num_conv_spec(j)
       end do
    end do
    
    call dspmv('U',num_bin,1.d0,mat_deconv,conv_spec,1,0.d0,power_spec(1:num_bin,idx),1)
    
    if(MpiRank==0) then
       write(str,*) idx
       file_spec=trim(rootname)//'_spectrum_'//trim(adjustl(str))//'.dat'
       open(100,file=trim(file_spec),action='write',iostat=info)
       do i=1,num_bin
          write(100,'(4E20.10e3)',iostat=info) (wav_spec(i)+wav_spec(i-1))/2, &
               power_spec(i,idx),conv_spec(i)/num_conv_spec(i),num_conv_spec(i)
       end do
       close(100)
    end if
    deallocate(conv_spec,num_conv_spec,stat=info)
    
    call estimate_covmat(idx,min_mc,max_mc,tol)
    
    deallocate(mat_deconv,stat=info)
    
    if(MpiRank==0) then
       write(str,*) idx
       file_cov=trim(rootname)//'_covariance_'//trim(adjustl(str))//'.dat'
       open(100,file=trim(file_cov),action='write',iostat=info)
       if(info/=0) call MpiStop('cannot open covariance.dat')
       do i=1,num_bin
          do j=1,i
             ij=i*(i-1)/2+j
             write(100,'(2I,E20.10e3)',iostat=info) i,j,power_cov(ij,idx)
          end do
       end do
       close(100)
    end if
    
    ! compute conformal energy density and its error
    call compute_crho(idx,crho_spec(1,idx),crho_spec(2,idx))
    if(MpiRank==0) write(*,'(A,3E16.4)') '  tau, crho, delta crho: ',tau_id(idx),crho_spec(1:2,idx)
    
    call MPI_BARRIER(MPI_COMM_WORLD,info)
    if(MpiRank==0) write(*,'(A)') ' *... done'
    
  end subroutine spectrum_estimation
  
  subroutine compute_crho(idx,rho,drho)
    integer,intent(in) :: idx
    double precision,intent(out) :: rho,drho
    integer i,j,ij
    double precision q,dq,q2dq,ddrho
    
    rho=0.d0
    drho=0.d0
    do i=1,num_bin
       ddrho=0.d0
       do j=1,i
          ij=i*(i-1)/2+j
          q=(wav_spec(j)+wav_spec(j-1))/2
          dq=wav_spec(j)-wav_spec(j-1)
          q2dq=q*q*dq/pi/pi ! factor 2 counted for i<->j
          if(j==i) q2dq=q2dq/2
          ddrho=ddrho+power_cov(ij,idx)*q2dq
       end do
       q=(wav_spec(i)+wav_spec(i-1))/2
       dq=wav_spec(i)-wav_spec(i-1)
       q2dq=q*q*dq/twopi/pi
       drho=drho+ddrho*q2dq
       rho=rho+power_spec(i,idx)*q2dq
    end do
    drho=sqrt(drho)
    
  end subroutine compute_crho
  
  subroutine compute_matrix_deconv(idx)
    use subroutines,only : check_syminv
#include <fftw3-mpi.f03>
    integer,intent(in) :: idx
    complex(kind=C_DOUBLE_COMPLEX) xc
    integer(8) ix,iy,iz
    double precision q,u,v
    integer info,imax2,imin2,iqx2,iqy2,iqz2,iq2,ij
    double precision,allocatable :: mat(:),mat0(:),tat(:),tat0(:)
    double precision,parameter :: tol=1.d-10
    
    call MPI_BARRIER(MPI_COMM_WORLD,info)
    if(MpiRank==0) write(*,'(A)') '  computing deconvolution matrix ...'
    
    xc=cmplx(dv_spa_spec)!,0d0,kind=C_DOUBLE_COMPLEX)
    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(iz,iy,ix) COLLAPSE(2)
    do iz=1,local_0_spec
       do iy=1,ngrid_spec
#if DIM==2
          data_spec2(iy,iz)=xc
#elif DIM==3
          !dir$ivdep
          do ix=1,ngrid_spec
             data_spec2(ix,iy,iz)=xc
          end do
#endif
       end do
    end do
    !$OMP END PARALLEL DO
    
    call masking(idx)
    
    call FFTW_MPI_EXECUTE_DFT(plan_backward_spec,data_spec2,data_spec1)
    if(MpiRank==0) write(*,*) '   backwward FFT done'
    
    allocate(mat(num_cov),mat0(num_cov),tat(num_cov),tat0(num_cov),stat=info)
    if(info/=0) call MpiStop('error: cannot allocate mat,mat0,...')
    
    mat0(1:num_cov)=0.d0
    tat0(1:num_cov)=0.d0
    
    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ij,imax2,imin2,ix,iy,iz,iqx2,iqy2,iqz2,iq2,q,v,u) REDUCTION(+:mat0,tat0) SCHEDULE(DYNAMIC,1)
    do ij=1,num_cov
       
!!$ debug
       imax2=ceiling((largeq_spec(num_largeq,ij)/dwav_spec)**2)-1
       imin2=ceiling((largeq_spec(0,ij)/dwav_spec)**2)
!!$ debug
       
       do iz=1,local_0_spec
          iqz2=(mod(iz+local_0_offset_spec+ngo2_spec-1,ngrid_spec)-ngo2_spec)**2
          if(iqz2>imax2) cycle
          
          do iy=1,ngrid_spec
             iqy2=(mod(iy+ngo2_spec-1,ngrid_spec)-ngo2_spec)**2
             iq2=iqz2+iqy2
#if DIM==2
!!$ debug
             if(.not.(iq2<=imax2.and.iq2=>imin2)) cycle
             q=sqrt(dble(iq2))*dwav_spec
!!$ debug
             v=abs(data_spec1(iy,iz))**2
             u=get_mat_spec(ij,q)
             mat0(ij)=mat0(ij)+v*u
             tat0(ij)=tat0(ij)+u
#elif DIM==3
             if(iq2>imax2) cycle
             do ix=1,ngrid_spec
                iqx2=(mod(ix+ngo2_spec-1,ngrid_spec)-ngo2_spec)**2
                iq2=iqx2+iqy2+iqz2
!!$ debug
                if(.not.(iq2<=imax2.and.iq2>=imin2)) cycle
                q=sqrt(dble(iq2))*dwav_spec
!!$ debug
                u=get_mat_spec(ij,q)
                v=abs(data_spec1(ix,iy,iz))**2
                mat0(ij)=mat0(ij)+v*u
                tat0(ij)=tat0(ij)+u
             end do
#endif
          end do
       end do
       
    end do
    !$OMP END PARALLEL DO
    write(*,*) '   MpiRank',MpiRank,'done k summation'
    
    call MPI_ALLREDUCE(mat0,mat,num_cov,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
    call MPI_ALLREDUCE(tat0,tat,num_cov,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
    mat(1:num_cov)=mat(1:num_cov)/tat(1:num_cov)*dv_wav_spec*dv_wav_spec
    
    ! inversion of matrix
    mat0=mat
    call dpptrf('U',num_bin,mat0,info)
    call dpptri('U',num_bin,mat0,info)
    call check_syminv(num_bin,mat,mat0,tol)
    mat_deconv(1:num_cov)=mat0(1:num_cov)
    
    deallocate(mat,mat0,tat0,tat,stat=info)
    
  end subroutine compute_matrix_deconv
  
  function get_mat_spec(ij,q)
    integer,intent(in) :: ij
    double precision,intent(in) :: q
    double precision get_mat_spec,d
    integer iq
    
    d=(q-largeq_spec(0,ij))/dq_spec(ij)
    iq=floor(d) ! round down
!!$ debug
    if(iq<num_largeq.and.iq>=0) then
       d=d-iq
       get_mat_spec=(1-d)*mat_spec(iq,ij)+d*mat_spec(iq+1,ij)
    else
       print '(I,2E20.12e2,I,2E20.12e2)',MpiRank,q,d,iq,largeq_spec(0,ij),largeq_spec(num_largeq,ij)
       get_mat_spec=0d0
    endif
!!$ debug
    
  end function get_mat_spec
  
  subroutine compute_conv_spec(idx)
#include <fftw3-mpi.f03>
    integer,intent(in) :: idx
    double precision,allocatable :: pow(:),pow0(:),num_spec(:),num_spec0(:)
    integer(8) ix,iy,iz,ixyz(ndim)
    double precision q,dq
    integer i,info
    
    call MPI_BARRIER(MPI_COMM_WORLD,info)
    if(MpiRank==0) write(*,'(A)') '  ... done and computing masked spectrum ...'
    dq=wav_spec(1)-wav_spec(0) ! bin spacing of wave number in power spectrum estimation
    
    allocate(pow(num_bin),pow0(num_bin),num_spec(num_bin),num_spec0(num_bin),stat=info)
    if(info/=0) call MpiStop('error: cannot allocate pow,...')

!!$ deallocate data_spec1 to saving memory usage
!!$    call FFTW_DESTROY_PLAN(plan_forward_spec)
!!$    call FFTW_DESTROY_PLAN(plan_backward_spec)
!!$    data_spec1=>NULL()
!!$    call FFTW_FREE(cdata_spec1)
!!$ deallocate data_spec1 to saving memory usage
    
    call data2_to_data_spec_degrade ! phidot
    
!!$ allocate data_spec1 to saving memory usage
!!$    cdata_spec1=FFTW_ALLOC_COMPLEX(alloc_local_spec)
!!$    if(.not.c_associated(cdata_spec1)) call MpiStop('error: cdata_spec1 is null')
!!$#if DIM==2
!!$    call c_f_pointer(cdata_spec1,data_spec1,[ngrid_spec,local_0_spec])
!!$#elif DIM==3
!!$    call c_f_pointer(cdata_spec1,data_spec1,[ngrid_spec,ngrid_spec,local_0_spec])
!!$#endif
!!$    if(.not.associated(data_spec1)) call MpiStop('error: data_spec1 is not allocated')
!!$    call FFTW_PLAN_WITH_NTHREADS(num_threads)
!!$#if DIM==2
!!$    plan_forward_spec=FFTW_MPI_PLAN_DFT_2D(ngrid_spec,ngrid_spec, &
!!$         data_spec1,data_spec2,MPI_COMM_WORLD,FFTW_FORWARD,FFTW_ESTIMATE)
!!$    plan_backward_spec=FFTW_MPI_PLAN_DFT_2D(ngrid_spec,ngrid_spec, &
!!$         data_spec2,data_spec1,MPI_COMM_WORLD,FFTW_BACKWARD,FFTW_ESTIMATE)
!!$#elif DIM==3
!!$    plan_forward_spec=FFTW_MPI_PLAN_DFT_3D(ngrid_spec,ngrid_spec,ngrid_spec, &
!!$         data_spec1,data_spec2,MPI_COMM_WORLD,FFTW_FORWARD,FFTW_ESTIMATE)
!!$    plan_backward_spec=FFTW_MPI_PLAN_DFT_3D(ngrid_spec,ngrid_spec,ngrid_spec, &
!!$         data_spec2,data_spec1,MPI_COMM_WORLD,FFTW_BACKWARD,FFTW_ESTIMATE)
!!$#endif
!!$    if(.not.c_associated(plan_forward_spec)) call MpiStop('error: plan_forward_spec is null')
!!$    if(.not.c_associated(plan_backward_spec)) call MpiStop('error: plan_backward_spec is null')
!!$ allocate data_spec1 to saving memory usage
    
    call masking(idx)
    call FFTW_MPI_EXECUTE_DFT(plan_backward_spec,data_spec2,data_spec1)
    
    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(iz,iy) COLLAPSE(2)
    do iz=1,local_0_spec
       do iy=1,ngrid_spec
#if DIM==2
          data_spec1(iy,iz)=data_spec1(iy,iz)*dv_spa_spec
#elif DIM==3
          data_spec1(1:ngrid_spec,iy,iz)=data_spec1(1:ngrid_spec,iy,iz)*dv_spa_spec
#endif
       end do
    end do
    !$OMP END PARALLEL DO
    
    pow(1:num_bin)=0.d0
    num_spec(1:num_bin)=0.d0
    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(iz,iy,ix,ixyz,q,i) REDUCTION(+:pow,num_spec) COLLAPSE(2)
    do iz=1,local_0_spec
       do iy=1,ngrid_spec
#if DIM==2
          ixyz(1:ndim)=(/iy,iz+local_0_offset_spec/)
          call get_comwav_spec(ixyz,q)
          i=nint(q/dq+0.5d0)
          if(i>=1.and.i<=num_bin) then
             pow(i)=pow(i)+abs(data_spec1(iy,iz))**2*dv_wav_spec
             num_spec(i)=num_spec(i)+1.d0
          end if
#elif DIM==3
          do ix=1,ngrid_spec
             ixyz(1:ndim)=(/ix,iy,iz+local_0_offset_spec/)
             call get_comwav_spec(ixyz,q)
             i=nint(q/dq+0.5d0)
             if(i>=1.and.i<=num_bin) then
                pow(i)=pow(i)+abs(data_spec1(ix,iy,iz))**2*dv_wav_spec
                num_spec(i)=num_spec(i)+1.d0
             end if
          end do
#endif
       end do
    end do
    !$OMP END PARALLEL DO
    
    call MPI_ALLREDUCE(pow,pow0,num_bin,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
    call MPI_ALLREDUCE(num_spec,num_spec0,num_bin,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
    
    conv_spec(1:num_bin)=pow0(1:num_bin)
    num_conv_spec(1:num_bin)=num_spec0(1:num_bin)
    
    deallocate(pow,pow0,num_spec,num_spec0,stat=info)
    if(MpiRank==0) write(*,'(A)') '  ... done'
    
  end subroutine compute_conv_spec
  
  subroutine estimate_covmat(idx,min_mc,max_mc,tol)
    use subroutines
#include <fftw3-mpi.f03>
    integer,intent(in) :: idx,min_mc,max_mc
    double precision,intent(in) :: tol
    integer info,i,j,k,ij,ii,iter,ithread,mult,nxx
    double precision q,dq,diff,x,u,v
    
    integer(8) ix,iy,iz,ixyz(ndim),nskip
    complex(C_DOUBLE_COMPLEX) xc
    double precision,allocatable :: xx(:,:),pow(:),pow0(:),power_cov_tmp(:,:),sigma(:)
    
    call MPI_BARRIER(MPI_COMM_WORLD,info)
    if(MpiRank==0) write(*,'(A)') ' *estimating covariance matrix using MC (Gaussianity assumed)'

    nxx=mult_comp*ngrid_spec
    allocate(xx(nxx,num_threads),pow(num_bin),pow0(num_bin),power_cov_tmp(num_cov,max_mc),sigma(num_bin),stat=info)
    if(info/=0) call MpiStop('error: cannot allocate xx, ...')
    
    dq=wav_spec(1)-wav_spec(0) ! bin spacing of wave number in power spectrum estimation
    mult=mult_gauss*mult_comp
    sigma(1:num_bin)=sqrt(power_spec(1:num_bin,idx)/2)/ngrid_tot_spec ! amplitude
    
    power_cov(1:num_cov,idx)=0d0
    do iter=1,max_mc
       
       nskip=local_0_offset_spec*ngrid_spec**(ndim-1)
       call skip_rand_stream0(nskip,mult)
       call put_rand_stream0
       call skip_rand_stream0(ngrid_tot_spec-nskip,mult)
       
       !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ithread,nskip,iz,iy,ix,ixyz,q,i,xc)
#ifdef _OPENMP
       ithread=OMP_GET_THREAD_NUM()+1
#else
       ithread=1
#endif
       nskip=get_openmp_nskip(ithread,num_threads,local_0_spec)
       nskip=nskip*ngrid_spec**(ndim-1)
       call skip_rand_stream(ithread,nskip,mult)
       !$OMP DO SCHEDULE(STATIC) ! static and no-collapse for reproducibility
       do iz=1,local_0_spec
#if DIM==2
          call gen_rand_gauss(ithread,xx(1:nxx,ithread),nxx)
          do iy=1,ngrid_spec
             ixyz(1:ndim)=(/iy,iz+local_0_offset_spec/)
             call get_comwav_spec(ixyz,q)
             i=nint(q/dq+0.5d0)
             if(i<=num_bin) then
                xc=cmplx(xx(mult_comp*iy-1,ithread),xx(mult_comp*iy,ithread),kind=C_DOUBLE_COMPLEX)
                data_spec1(iy,iz)=xc*sigma(i)
             else
                data_spec1(iy,iz)=cmplx(0d0,kind=C_DOUBLE_COMPLEX)
             end if
          end do
#elif DIM==3
          do iy=1,ngrid_spec
             call gen_rand_gauss(ithread,xx(1:nxx,ithread),nxx)
             do ix=1,ngrid_spec
                ixyz(1:ndim)=(/ix,iy,iz+local_0_offset_spec/)
                call get_comwav_spec(ixyz,q)
                i=nint(q/dq+0.5d0)
                if(i<=num_bin) then
                   xc=cmplx(xx(mult_comp*ix-1,ithread),xx(mult_comp*ix,ithread),kind=C_DOUBLE_COMPLEX)
                   data_spec1(ix,iy,iz)=xc*sigma(i)
                else
                   data_spec1(ix,iy,iz)=cmplx(0d0,kind=C_DOUBLE_COMPLEX)
                end if
             end do
          end do
#endif
       end do
       !$OMP END DO
       !$OMP END PARALLEL
       call put_rand_stream0       
       
       call FFTW_MPI_EXECUTE_DFT(plan_forward_spec,data_spec1,data_spec2)
       if(MpiRank==0) write(*,'(A)') '   forward FFT done'
       
       call masking(idx)
       call FFTW_MPI_EXECUTE_DFT(plan_backward_spec,data_spec2,data_spec1)
       if(MpiRank==0) write(*,'(A)') '   backward FFT done'
       
       pow(1:num_bin)=0.d0
       !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(iz,iy,ix,ixyz,q,i) REDUCTION(+:pow) COLLAPSE(2)
       do iz=1,local_0_spec
          do iy=1,ngrid_spec
#if DIM==2
             ixyz(1:ndim)=(/iy,iz+local_0_offset_spec/)
             call get_comwav_spec(ixyz,q)
             i=nint(q/dq+0.5d0)
             if(i<=num_bin) pow(i)=pow(i)+abs(data_spec1(iy,iz))**2
#elif DIM==3
             do ix=1,ngrid_spec
                ixyz(1:ndim)=(/ix,iy,iz+local_0_offset_spec/)
                call get_comwav_spec(ixyz,q)
                i=nint(q/dq+0.5d0)
                if(i<=num_bin) pow(i)=pow(i)+abs(data_spec1(ix,iy,iz))**2
             end do
#endif
          end do
       end do
       !$OMP END PARALLEL DO
       
       call MPI_ALLREDUCE(pow,pow0,num_bin,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
       if(MpiRank==0) write(*,'(A)') '   spectrum estimation done'
       
       pow(1:num_bin)=-power_spec(1:num_bin,idx)
       call dspmv('U',num_bin,1.d0,mat_deconv,pow0,1,1.d0,pow,1)
       
       do i=1,num_bin
          do j=1,i
             ij=i*(i-1)/2+j
             u=pow(i)*pow(j)
             power_cov_tmp(ij,iter)=power_cov(ij,idx)
             do k=1,iter-1
                power_cov_tmp(ij,k)=power_cov_tmp(ij,k)+u
             end do
             power_cov(ij,idx)=power_cov(ij,idx)+u
          end do
       end do
       
       call MPI_BARRIER(MPI_COMM_WORLD,info)
       if(iter<min_mc) cycle
       
       if(MpiRank==0) write(*,'(A)') '   diagnosing convergence ... '
       
       ! convergence diagnostic with jackknife method
       ! only diagonal part
       diff=0d0
       do i=1,num_bin
          ij=i*(i+1)/2
          u=0d0
          v=0d0
          do k=1,iter
             x=power_cov_tmp(ij,k)/(iter-1)
             u=u+x
             v=v+x*x
          end do
          u=u/iter
          v=v/iter
          v=sqrt(v-u*u)
          x=power_cov(ij,idx)/iter
          diff=max(diff,v/x)
       end do
       
       call MPI_BARRIER(MPI_COMM_WORLD,info)
       if(MpiRank==0) write(*,'(A,I,E20.10e3)') '  iter, diff:',iter,diff
       
       if(diff<tol) exit
       
    end do
    
    deallocate(xx,pow,pow0,power_cov_tmp,sigma,stat=info)
    
    if(MpiRank==0) then
       if(iter<=max_mc) then
          write(*,'(A)') '  MC converged'
       else
          write(*,'(A)') '  WARNING: MC failed to converge'
       end if
    end if
    
    power_cov(1:num_cov,idx)=power_cov(1:num_cov,idx)/iter    
    
  end subroutine estimate_covmat
  
  subroutine data2_to_data_spec_degrade
    integer,parameter :: tg=0
    integer info,deg,MpiRankD,MpiRankA,count,ssend(MPI_STATUS_SIZE),srecv(MPI_STATUS_SIZE), &
         irank,irecv0,isend0
    integer(8),allocatable :: min_1t(:),max_1t(:),local_1t(:),local_1t_offset(:)
    integer,allocatable :: irecv(:),isend(:)
#if DIM==2
    complex(kind(0d0)),allocatable :: data_phi_degrade(:,:),data_phi_tmp(:)
#elif DIM==3
    complex(kind(0d0)),allocatable :: data_phi_degrade(:,:,:),data_phi_tmp(:,:)
#endif
    complex(kind(0d0)) xc
    integer(8) ix,iy,iz,ix0,iy0,iz0,j,k,min_1,max_1
    
    call MPI_BARRIER(MPI_COMM_WORLD,info)
    if(MpiRank==0) write(*,*) '  degrading and distributing phi/phidot from data2 to data_spec ...'
    deg=2**degrade_spec
    min_1=local_0_offset/deg+1
    max_1=(local_0_offset+local_0-1)/deg+1
    write(*,'(A,4I)') '  information of data_phi_degrade:',MpiRank,deg,min_1,max_1
    
#if DIM==2
    allocate(data_phi_degrade(ngrid_spec,min_1:max_1),data_phi_tmp(ngrid_spec),stat=info)
#elif DIM==3
    allocate(data_phi_degrade(ngrid_spec,ngrid_spec,min_1:max_1),data_phi_tmp(ngrid_spec,ngrid_spec),stat=info)
#endif
    if(info/=0) call MpiStop('error: cannot allocate data_phi_degrade, data_phi_tmp')
    
    ! degrade
    xc=cmplx(0d0)
    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(iz,iy,ix,iz0,iy0,ix0) COLLAPSE(ndim)
    do iz=min_1,max_1
       do iy=1,ngrid_spec
#if DIM==2
          data_phi_degrade(iy,iz)=xc
          do iz0=max(1,(iz-1)*deg+1-local_0_offset)+1,min(local_0,iz*deg-local_0_offset)+1
             do iy0=(iy-1)*deg+1,iy*deg
                data_phi_degrade(iy,iz)=data_phi_degrade(iy,iz)+data2(iy0,iz0)
             end do
          end do
          data_phi_degrade(iy,iz)=data_phi_degrade(iy,iz)/deg**ndim
#elif DIM==3
          do ix=1,ngrid_spec
             data_phi_degrade(ix,iy,iz)=xc
             do iz0=max(1,(iz-1)*deg+1-local_0_offset)+1,min(local_0,iz*deg-local_0_offset)+1
                do iy0=(iy-1)*deg+1,iy*deg
                   do ix0=(ix-1)*deg+1,ix*deg
                      data_phi_degrade(ix,iy,iz)=data_phi_degrade(ix,iy,iz)+data2(ix0,iy0,iz0)
                   end do
                end do
             end do
             data_phi_degrade(ix,iy,iz)=data_phi_degrade(ix,iy,iz)/deg**ndim
          end do
#endif
       end do
    end do
    !$OMP END PARALLEL DO
    
    count=ngrid_spec**(ndim-1)
    if(mod(local_0+local_0_offset,deg)/=0) then
       MpiRankA=modulo(MpiRank+1,MpiSize)
#if DIM==2
       call MPI_ISEND(data_phi_degrade(1,max_1),count,MPI_DOUBLE_COMPLEX,MpiRankA,tg,MPI_COMM_WORLD,isend0,info)
#elif DIM==3
       call MPI_ISEND(data_phi_degrade(1,1,max_1),count,MPI_DOUBLE_COMPLEX,MpiRankA,tg,MPI_COMM_WORLD,isend0,info)
#endif
    end if
    if(mod(local_0_offset,deg)/=0) then
       MpiRankD=modulo(MpiRank-1,MpiSize)
       call MPI_IRECV(data_phi_tmp,count,MPI_DOUBLE_COMPLEX,MpiRankD,tg,MPI_COMM_WORLD,irecv0,info)
    end if
    
    if(mod(local_0+local_0_offset,deg)/=0) call MPI_WAIT(isend0,ssend,info)
    if(mod(local_0_offset,deg)/=0) then
       call MPI_WAIT(irecv0,ssend,info)
       !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(iy0)
       do iy0=1,ngrid_spec
#if DIM==2
          data_phi_degrade(iy0,min_1)=data_phi_degrade(iy0,min_1)+data_phi_tmp(iy0)
#elif DIM==3
          data_phi_degrade(1:ngrid_spec,iy0,min_1)=data_phi_degrade(1:ngrid_spec,iy0,min_1)+data_phi_tmp(1:ngrid_spec,iy0)
#endif
       end do
       !$OMP END PARALLEL DO
    end if
    
    deallocate(data_phi_tmp,stat=info)
    
    ! data_phi_degrade to data_spec
    allocate(local_1t(0:MpiSize-1),local_1t_offset(0:MpiSize-1),min_1t(0:MpiSize-1),max_1t(0:MpiSize-1), &
         irecv(0:MpiSize-1),isend(0:MpiSize-1),stat=info)
    if(info/=0) call MpiStop('error: cannot allocate local_1t, ...')
    call MPI_ALLGATHER(min_1,1,MPI_INTEGER8,min_1t,1,MPI_INTEGER8,MPI_COMM_WORLD,info)
    call MPI_ALLGATHER(max_1,1,MPI_INTEGER8,max_1t,1,MPI_INTEGER8,MPI_COMM_WORLD,info)
    call MPI_ALLGATHER(local_0_spec,1,MPI_INTEGER8,local_1t,1,MPI_INTEGER8,MPI_COMM_WORLD,info)
    call MPI_ALLGATHER(local_0_offset_spec,1,MPI_INTEGER8,local_1t_offset,1,MPI_INTEGER8,MPI_COMM_WORLD,info)
    
    do irank=0,MpiSize-1
       
       if(irank==MpiRank) cycle
       
       ! send
       j=max(min_1,local_1t_offset(irank)+1)
       k=min(max_1,local_1t_offset(irank)+local_1t(irank))
       if(k>=j) then
          count=(k-j+1)*ngrid_spec**(ndim-1)
#if DIM==2
          call MPI_ISEND(data_phi_degrade(1,j),count,MPI_DOUBLE_COMPLEX,irank,tg,MPI_COMM_WORLD,isend(irank),info)
#elif DIM==3
          call MPI_ISEND(data_phi_degrade(1,1,j),count,MPI_DOUBLE_COMPLEX,irank,tg,MPI_COMM_WORLD,isend(irank),info)
#endif
       end if
       
       ! receive
       j=max(min_1t(irank),local_0_offset_spec+1)
       k=min(max_1t(irank),local_0_offset_spec+local_0_spec)
       if(k>=j) then
          iz=j-local_0_offset_spec
          count=(k-j+1)*ngrid_spec**(ndim-1)
#if DIM==2
          call MPI_IRECV(data_spec2(1,iz),count,MPI_DOUBLE_COMPLEX,irank,tg,MPI_COMM_WORLD,irecv(irank),info)
#elif DIM==3
          call MPI_IRECV(data_spec2(1,1,iz),count,MPI_DOUBLE_COMPLEX,irank,tg,MPI_COMM_WORLD,irecv(irank),info)
#endif
       end if
       
    end do
    
    ! within node
    j=max(min_1,local_0_offset_spec+1)
    k=min(max_1,local_0_offset_spec+local_0_spec)
    if(k>=j) then
       !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(iz,iy) COLLAPSE(2)
       do iz=j-local_0_offset_spec,k-local_0_offset_spec
          do iy=1,ngrid_spec
#if DIM==2
             data_spec2(iy,iz)=data_phi_degrade(iy,iz+local_0_offset_spec)
#elif DIM==3
             data_spec2(1:ngrid_spec,iy,iz)=data_phi_degrade(1:ngrid_spec,iy,iz+local_0_offset_spec)
#endif
          end do
       end do
       !$OMP END PARALLEL DO
    end if
    
    do irank=0,MpiSize-1
       if(irank==MpiRank) cycle
       
       j=max(min_1,local_1t_offset(irank)+1)
       k=min(max_1,local_1t_offset(irank)+local_1t(irank))
       if(k>=j) call MPI_WAIT(isend(irank),ssend,info)
       
       j=max(min_1t(irank),local_0_offset_spec+1)
       k=min(max_1t(irank),local_0_offset_spec+local_0_spec)
       if(k>=j) call MPI_WAIT(irecv(irank),srecv,info)
    end do
    
    deallocate(data_phi_degrade,min_1t,max_1t,local_1t,local_1t_offset,irecv,isend,stat=info)
    
    call MPI_BARRIER(MPI_COMM_WORLD,info)
    if(MpiRank==0) write(*,*) '  ... done'
    
  end subroutine data2_to_data_spec_degrade
    
  subroutine masking(idx)
    integer,intent(in) :: idx
    integer(8) ix0,iy0,iz0
    integer info,irank,min_1t,max_1t,deg
    character(len=charlen) file_mask,str,str0
    
    write(str,*) idx
    str=adjustl(str)
    deg=2**degrade_spec
    
    do irank=0,MpiSize-1
       
       min_1t=arr_local_0_offset(irank)+1
       max_1t=arr_local_0_offset(irank)+arr_local_0(irank)
       
       min_1t=(min_1t-1)/deg+1
       max_1t=(max_1t-1)/deg+1
       
       if(local_0_offset_spec+1<=max_1t.and.local_0_offset_spec+local_0_spec>=min_1t) then
          
          write(str0,*) irank
          str0=adjustl(str0)
          file_mask=trim(rootname)//'_mask_'//trim(str)//'-'//trim(str0)//'.dat'
          
          open(100,file=file_mask,status='old',action='read',iostat=info)
          if(info/=0) call MpiStop('error: cannot open mask.dat')
          do 
#if DIM==2
             read(100,*,iostat=info) iy0,iz0
#elif DIM==3
             read(100,*,iostat=info) ix0,iy0,iz0
#endif
             if(info>1) call MpiStop('error: cannot read mask.dat')
             if(info<0) exit
             
             ix0=(ix0-1)/deg+1
             iy0=(iy0-1)/deg+1
             iz0=(iz0-1)/deg+1
             iz0=iz0-local_0_offset_spec
#if DIM==2
             if(iz0>=1.and.iz0<=local_0_spec) data_spec2(iy0,iz0)=0.d0
#elif DIM==3
             if(iz0>=1.and.iz0<=local_0_spec) data_spec2(ix0,iy0,iz0)=0.d0
#endif
          end do
          close(100)
          
       end if
       
    end do
    
    call MPI_BARRIER(MPI_COMM_WORLD,info)
    if(MpiRank==0) write(*,*) '   masking done'
    
  end subroutine masking
  
  subroutine get_epsilon
    integer idx,info,i,j,ij
    double precision eps,deps,num,den,d,q,tau,dum1,dum2,f
    character(len=charlen) str,file_spec,file_cov,file_eps
    
    if(MpiRank==0) then
       write(*,*)
       write(*,'(A)') ' #epsilon'
    end if
    
    do idx=n_chk-n_id+1,n_chk
       
       write(str,*) idx
       file_spec=trim(rootname)//'_spectrum_'//trim(adjustl(str))//'.dat'
       open(100,file=trim(file_spec),action='read',status='old',iostat=info)
       if(info/=0) call MpiStop('error: cannot open spectrum.dat')
       do i=1,num_bin
          read(100,*,iostat=info) dum1, &
               power_spec(i,idx),dum2
          if(info/=0) call MpiStop('error: cannot read spectrum.dat')
       end do
       close(100)
       
       file_cov=trim(rootname)//'_covariance_'//trim(adjustl(str))//'.dat'
       open(100,file=trim(file_cov),action='read',status='old',iostat=info)
       if(info/=0) call MpiStop('error: cannot open covariance.dat')
       do i=1,num_bin
          do j=1,i
             ij=i*(i-1)/2+j
             read(100,*,iostat=info) dum1,dum2,power_cov(ij,idx)
             if(info/=0) call MpiStop('error: cannot read spectrum.dat')
          end do
       end do
       close(100)
       
!!$       call compute_crho(idx,crho_spec(1,idx),crho_spec(2,idx))
!!$       if(MpiRank==0) write(200,'(I,3E)',iostat=info) idx,tau_id(idx),crho_spec(1,idx),sqrt(crho_spec(2,idx))
       
    end do
    
    if(MpiRank==0) then
       file_eps=trim(rootname)//'_eps.txt'
       open(200,file=trim(file_eps),action='write',status='replace',iostat=info)
       if(info/=0) call MpiStop('error: cannot open eps.txt')
    end if
    
    do idx=n_chk-n_id+1,n_chk-n_mod_eps,n_mod_eps
       tau=(tau_id(idx+n_mod_eps)+tau_id(idx))/2
       f=tau/2/twopi
       num=0.d0
       den=0.d0
       do i=1,num_bin
          d=power_spec(i,idx+n_mod_eps)-power_spec(i,idx) ! differential power spectrum
          q=(wav_spec(i-1)+wav_spec(i))/2
          !!$num=num+d
          !!$den=den+d/q
          num=num+d*q*q
          den=den+d*q
       end do
       eps=num/den*f
       
       deps=0d0
       do i=1,num_bin
          dum1=(wav_spec(i-1)+wav_spec(i))/2
          do j=1,i
             ij=i*(i-1)/2+j
             dum2=(wav_spec(j-1)+wav_spec(j))/2
             !!$d=(1-eps/dum1)*(1-eps/dum2)*(power_cov(ij,idx+n_mod_eps)+power_cov(ij,idx)) ! this would be too conservative???
             d=(1-eps/dum1/f)*(1-eps/dum2/f)*(power_cov(ij,idx+n_mod_eps)+power_cov(ij,idx))*dum1**2*dum2**2
             if(j/=i) d=d*2
             deps=deps+d
          end do
       end do
       deps=sqrt(deps/den/den)*f
       
       !!if(idx==n_chk-n_id+1) call compute_crho(idx,crho_spec(1,idx),crho_spec(2,idx))
       !!call compute_crho(idx+n_mod_eps,crho_spec(1,idx+n_mod_eps),crho_spec(2,idx+n_mod_eps))
       if(MpiRank==0) &
            write(200,'(2I,7E20.10e3)') idx,idx+n_mod_eps,tau_id(idx),tau_id(idx+n_mod_eps), &
            tau,eps,deps,crho_spec(1,idx+n_mod_eps)-crho_spec(1,idx),sqrt(crho_spec(2,idx+n_mod_eps)**2+crho_spec(2,idx)**2)
    end do
    
    if(MpiRank==0) close(200)
    
  end subroutine get_epsilon

end module spectrum
