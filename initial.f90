module initial
  use settings
  use subroutines
  implicit none
  
  type(C_PTR) plan_forward,cdata1,cdata2
  complex(C_DOUBLE_COMPLEX),pointer :: ldata1(:),ldata2(:)
  integer(8) alloc_local,alloc_local_mod
  
  double precision,parameter :: srj_tol=5.d-3
  integer,parameter :: srj_levels=6
  integer,parameter :: min_srj_levels=6,max_srj_levels=15
  type srj_param
     integer levels
     double precision,allocatable :: omega(:),beta(:)
     integer,allocatable :: q(:)
  end type srj_param
  type(srj_param) srjp
  
contains
  subroutine set_initial_condition(chk)
    logical,intent(out) :: chk
    integer info,idx
    double precision tau
    integer,external :: access
    character(len=charlen) chk_phi,str
    
    if(MpiRank==0) then
       write(*,*) 
       write(*,'(A)') ' #set initial condition'
    end if
    
    chk=.false.
    
    if(access(file_resume,'r')==0) then ! resume file exists
       
       open(100,file=file_resume,action='read',status='old',iostat=info)
       read(100,*,iostat=info) idx,tau
       close(100)
       
       if(info==0.and.idx>=0.and.idx<=n_chk) then
          
          write(str,*) idx
          str=adjustl(str)
          if(tmpdir/='') then
             chk_phi=trim(tmpdir)//trim(tmprn)//'_phi_'//trim(str)//'.chk'
          else
             chk_phi=trim(rootname)//'_phi_'//trim(str)//'.chk'
          end if
          
          if(access(chk_phi,'r')==0) then
             if(MpiRank==0) write(*,'(A)') ' *** chk files EXIST!; restart from previous run ***'
             chk=.true.
             call allocate_data(chk)
             return
          end if
          
       end if
    end if
    
    if(MpiRank==0) write(*,'(A)') ' *generating new initial configuration ...'
    
    if(initial_loop) then
       call get_initial_loop
    else
       call get_initial_thermal(chk)
    end if
    
    if(MpiRank==0) then
       open(100,file=file_resume,action='write',status='replace',iostat=info)
       close(100)
       write(*,'(A)') '  *param.chk is created'
    end if
    
    call MPI_BARRIER(MPI_COMM_WORLD,info)
    if(MpiRank==0) write(*,'(A)') ' *... done'
    
  end subroutine set_initial_condition
  
  subroutine get_initial_thermal(chk)
#include <fftw3-mpi.f03>
    logical,intent(in) :: chk
    integer(8) ix,iy,iz,ixyz(ndim),nskip
    double precision,allocatable :: xx(:)
    integer info,idx,ip,mult,ithread,nxx,iter
    integer,external :: access
!!$!!$!!$
    integer(8) iz0,iz1
    character(len=1024) fn
!!$!!$!!$
    
    if(MpiRank==0) write(*,*) ' #starting from thermal distribution before phase transition'
!!$!!$!!$ debug
    call allocate_data(chk)
!!$!!$!!$ debug
    
    mult=mult_gauss*mult_comp
    nxx=ngrid*mult_comp
    
    do ip=1,mult_phi
       if(ip==1) then
          if(MpiRank==0) write(*,*) '  phi:'
       else if(ip==mult_phi) then
          if(MpiRank==0) write(*,*) '  phidot:'
       end if
       
       nskip=local_0_offset*ngrid_slice
       call skip_rand_stream0(nskip,mult)
       call put_rand_stream0
       call skip_rand_stream0(ngrid_tot-nskip,mult)
       
       !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ithread,nskip,iz,iy,ix,ixyz,xx)
       allocate(xx(nxx),stat=info)
       if(info/=0) call MpiStop('error: cannot allocate xx')
#ifdef _OPENMP
       ithread=OMP_GET_THREAD_NUM()+1
#else
       ithread=1
#endif
       nskip=get_openmp_nskip(ithread,num_threads,local_0)
       nskip=nskip*ngrid_slice
       call skip_rand_stream(ithread,nskip,mult)
       !$OMP DO SCHEDULE(STATIC) ! static and no-collapse for reproducibility
       do iz=1+local_0_offset,local_0+local_0_offset
#if DIM==2
          call gen_rand_gauss(ithread,xx,nxx)
          xx(1:nxx)=xx(1:nxx)*dv_wav
          do iy=1,ngrid
             ixyz(1:ndim)=(/iy,iz/)
             data2(iy,iz-local_0_offset+1)=cmplx(xx(mult_comp*iy-1),xx(mult_comp*iy))
             call initial_phi_phidot(ip,ixyz,data2(iy,iz-local_0_offset+1))
          end do
#elif DIM==3
          do iy=1,ngrid
             call gen_rand_gauss(ithread,xx,nxx)
             xx(1:nxx)=xx(1:nxx)*dv_wav
             do ix=1,ngrid
                ixyz(1:ndim)=(/ix,iy,iz/)
                data2(ix,iy,iz-local_0_offset+1)=cmplx(xx(mult_comp*ix-1),xx(mult_comp*ix))
                call initial_phi_phidot(ip,ixyz,data2(ix,iy,iz-local_0_offset+1))
             end do
          end do
#endif
       end do
       !$OMP END DO
       deallocate(xx,stat=info)
       !$OMP END PARALLEL
       
       call put_rand_stream0
       
       call MPI_BARRIER(MPI_COMM_WORLD,info)
       if(MpiRank==0) write(*,*) '   generation in Fourier space done'
       
       if(ip==1) call init_fftw_initial
#if DIM==2
       call FFTW_MPI_EXECUTE_DFT(plan_forward,data2(:,min_0),data2(:,min_0))
#else if DIM==3
       call FFTW_MPI_EXECUTE_DFT(plan_forward,data2(:,:,min_0),data2(:,:,min_0))
#endif
       if(MpiRank==0) write(*,*) '   forward FFT done'
       if(ip==mult_phi) call fin_fftw_initial
       
       if(ip==1) then
          
          call allocate_buffer(mult_phi)
!!$
          write(*,'(A,I,E)') '    temporal buffer size:',MpiRank,byte_dc*alloc_local/1024/1024/1024
#if USETMP==1
          write(fn,*) mpirank
          fn=trim(rootname)//'_phi_'//trim(adjustl(fn))//'.dat.tmp_'
          open(199,file=trim(fn),access='sequential',form='unformatted',action='readwrite',iostat=info)
#else
          open(199,access='sequential',form='unformatted',status='scratch',action='readwrite',iostat=info)
#endif
          if(info/=0) call MpiStop('error: cannot open temporal output file')

!!$
!!$          write(*,'(A,I,E)') '    temporal buffer size:',MpiRank,byte_dc*alloc_local/1024/1024/1024
!!$          open(199,access='sequential',form='unformatted',status='scratch',action='readwrite',iostat=info)
!!$          if(info/=0) call MpiStop('error: cannot open temporal output file')
          iz0=ngrid_slice*(1+local_0)
          iter=0
          do iz=1+ngrid_slice,iz0,count_dc_max
             iz1=min(iz+count_dc_max-1,iz0)
             iter=iter+1
             write(*,'(A,3I)') '    iterations:',MpiRank,iter,(iz0-ngrid_slice-1)/count_dc_max+1
             write(199,iostat=info) ldata2(iz:iz1)
             if(info/=0) call MpiStop('error: cannot write in temporal output file')
          end do
!!$          iter=0
!!$          do iz=1+ngrid_slice,alloc_local+ngrid_slice,count_dc_max
!!$             iter=iter+1
!!$             write(199,iostat=info) ldata2(iz:min(iz+count_dc_max-1,alloc_local+ngrid_slice))
!!$             if(info/=0) call MpiStop('error: cannot write in temporal output file')
!!$             write(*,'(A,3I)') '    iterations:',MpiRank,iter,(alloc_local-1)/count_dc_max+1
!!$          end do
          
          if(MpiRank==0) write(*,*) '   phi done'
          
       else if(ip==mult_phi) then
          
          rewind(199)
          
          iz0=ngrid_slice*(1+local_0)
          iter=0
          do iz=1+ngrid_slice,iz0,count_dc_max
             iz1=min(iz+count_dc_max-1,iz0)
             iter=iter+1
             write(*,'(A,3I)') '    iterations:',MpiRank,iter,(iz0-ngrid_slice-1)/count_dc_max+1
             read(199,iostat=info) ldata1(iz:iz1)
             if(info/=0) call MpiStop('error: cannot read from temporal output file')
          end do
!!$          iter=0
!!$          do iz=1+ngrid_slice,alloc_local+ngrid_slice,count_dc_max
!!$             iter=iter+1
!!$             read(199,iostat=info) ldata1(iz:min(iz+count_dc_max-1,alloc_local+ngrid_slice))
!!$             if(info/=0) call MpiStop('error: cannot read from temporal output file')
!!$             write(*,'(A,3I)') '    iterations:',MpiRank,iter,(alloc_local-1)/count_dc_max+1
!!$          end do
!!$          close(199)
#if USETMP==1
          close(199,status='delete')
#else
          close(199)
#endif
          
          call deallocate_buffer
          if(MpiRank==0) write(*,*) '   phidot done'
          
       end if
       
    end do
    
#if DIM==3
    call check_initial_var
#endif
    
  end subroutine get_initial_thermal
  
  subroutine allocate_data(chk)
#include <fftw3-mpi.f03>
    logical,intent(in) :: chk ! if true, no initialization is needed and data1 is allocated; otherwise only data2 is allocated here
    integer info
    
#if DIM==2
    alloc_local=FFTW_MPI_LOCAL_SIZE_2D(ngrid, &
         ngrid,MPI_COMM_WORLD,local_0,local_0_offset)
#elif DIM==3
    alloc_local=FFTW_MPI_LOCAL_SIZE_3D(ngrid,ngrid, &
         ngrid,MPI_COMM_WORLD,local_0,local_0_offset)
#endif

    allocate(arr_local_0(0:MpiSize-1),arr_local_0_offset(0:MpiSize-1),stat=info)
    if(info/=0) stop 'error: cannot allocate arr_local_0, arr_local_0_offset'
    call MPI_ALLGATHER(local_0,1,MPI_INTEGER8,arr_local_0,1,MPI_INTEGER8,MPI_COMM_WORLD,info)
    call MPI_ALLGATHER(local_0_offset,1,MPI_INTEGER8,arr_local_0_offset,1,MPI_INTEGER8,MPI_COMM_WORLD,info)
    
    alloc_local_mod=max(alloc_local+ngrid_slice,(local_0+2)*ngrid_slice)
    if(feedback>1) write(*,'(A,5I)') '  information of allocated data:',MpiRank,local_0,local_0_offset,alloc_local,alloc_local_mod
    min_0=2
    max_0=local_0+1

    if(chk) then
       cdata1=FFTW_ALLOC_COMPLEX(alloc_local_mod)
       if(.not.c_associated(cdata1)) call MpiStop('error: cdata1 is null')
#if DIM==2
       call c_f_pointer(cdata1,data1,[ngrid,local_0+2])
#elif DIM==3
       call c_f_pointer(cdata1,data1,[ngrid,ngrid,local_0+2])
#endif
       if(.not.associated(data1)) call MpiStop('error: data1 are not associated')
    end if
    
    cdata2=FFTW_ALLOC_COMPLEX(alloc_local_mod)
    if(.not.c_associated(cdata2)) call MpiStop('error: cdata2 is null')
#if DIM==2
    call c_f_pointer(cdata2,data2,[ngrid,local_0+2])
#elif DIM==3
    call c_f_pointer(cdata2,data2,[ngrid,ngrid,local_0+2])
#endif
    if(.not.associated(data2)) call MpiStop('error: data2 are not associated')
    
    if(feedback>1) call memory_check_mkl()
    
  end subroutine allocate_data
  
  subroutine deallocate_data
#include <fftw3-mpi.f03>
    integer info
    
    deallocate(arr_local_0,arr_local_0_offset,stat=info)
    data1=>NULL()
    data2=>NULL()
    call FFTW_FREE(cdata1)
    call FFTW_FREE(cdata2)
    
  end subroutine deallocate_data
  
  subroutine allocate_buffer(i)
#include <fftw3-mpi.f03>
    integer,intent(in) :: i

    if(i==1) then
       call c_f_pointer(cdata1,ldata1,[alloc_local_mod])
       if(.not.associated(ldata1)) call MpiStop('error: ldata2 is not associated')
    else if(i==mult_phi) then
       call c_f_pointer(cdata2,ldata2,[alloc_local_mod])
       if(.not.associated(ldata2)) call MpiStop('error: ldata2 is not associated')
    end if
    if(feedback>1) call memory_check_mkl()

  end subroutine allocate_buffer
  
  subroutine deallocate_buffer
#include <fftw3-mpi.f03>
    
    ldata1=>NULL()
    ldata2=>NULL()
    
    if(feedback>1) call memory_check_mkl()
    
  end subroutine deallocate_buffer
  
  subroutine init_fftw_initial
#include <fftw3-mpi.f03>
    
    if(fft_threaded) call FFTW_PLAN_WITH_NTHREADS(int(num_threads,kind=4))
#if DIM==2
    plan_forward=FFTW_MPI_PLAN_DFT_2D(ngrid,ngrid, & ! in-place
         data2(:,min_0),data2(:,min_0),MPI_COMM_WORLD,FFTW_FORWARD,FFTW_ESTIMATE)
#elif DIM==3
    plan_forward=FFTW_MPI_PLAN_DFT_3D(ngrid,ngrid,ngrid, & ! in-place
         data2(:,:,min_0),data2(:,:,min_0),MPI_COMM_WORLD,FFTW_FORWARD,FFTW_ESTIMATE)
#endif
    if(.not.c_associated(plan_forward)) call MpiStop('error: plan_forward is null')
    if(feedback>1) call memory_check_mkl()
    
  end subroutine init_fftw_initial
  
  subroutine fin_fftw_initial
#include <fftw3-mpi.f03>
    
    call FFTW_DESTROY_PLAN(plan_forward)
    
    ! allocate data1
    cdata1=FFTW_ALLOC_COMPLEX(alloc_local_mod)
    if(.not.c_associated(cdata1)) call MpiStop('error: cdata1 is null')
#if DIM==2
    call c_f_pointer(cdata1,data1,[ngrid,local_0+2])
#elif DIM==3
    call c_f_pointer(cdata1,data1,[ngrid,ngrid,local_0+2])
#endif
    if(.not.associated(data1)) call MpiStop('error: data1 are not associated')
    
!!$    call c_f_pointer(cdata1,ldata1,[alloc_local_mod])
!!$    if(.not.associated(ldata1)) call MpiStop('error: ldata1 are not associated')
    call allocate_buffer(1)
    
    if(feedback>1) call memory_check_mkl()
    
  end subroutine fin_fftw_initial
  
  subroutine get_comwav_initial(ixyz,q) ! should be modified to include degrading?
    integer(8),intent(in) :: ixyz(ndim)
    double precision,intent(out) :: q
    
    q=sum(dble((mod(ixyz(1:ndim)+ngo2-1,ngrid)-ngo2)**2))
    q=sqrt(q)*dwav
    
  end subroutine get_comwav_initial
  
  subroutine initial_phi_phidot(ip,idx,xc)
    use potential
    integer,intent(in) :: ip
    complex(C_DOUBLE_COMPLEX),intent(inout) :: xc
    integer(8),intent(in) :: idx(ndim)
    double precision thrm,osc,wav2,var,fbe
    
    call get_comwav_initial(idx,wav2)
    wav2=(wav2/tau_ini/sigma_cr)**2 ! physical momentum (in sigma)
    thrm=thrm_sigma/tau_ini ! initial temperature (in sigma)
    osc=sqrt(msquared(thrm)+wav2) ! (in sigma)
    
    fbe=exp(-osc/thrm)
    fbe=fbe/(1-fbe)
    
    var=fbe*osc/dv_wav/sigma_cr ! variance of phidot (in sigma/tau_cr)
    xc=xc*sqrt(var/2)/(osc*sigma_cr)**(2-ip)
    
  end subroutine initial_phi_phidot
  
#if DIM==3
  subroutine check_initial_var
    use potential,only : pi,thrm_sigma
    integer(8) ix,iy,iz
    integer info
    integer,parameter :: ncov=mult_comp*(mult_comp+1)/2
    double precision cov(ncov,mult_phi),cov0(ncov,mult_phi),x(mult_comp),v
    
    cov(1:ncov,1:mult_phi)=0.d0
    
    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(iz,iy,ix,x) REDUCTION(+:cov) COLLAPSE(3)
    do iz=min_0,max_0
       do iy=1,ngrid
          do ix=1,ngrid
             
             x(1)=dble(data1(ix,iy,iz))
             x(mult_comp)=aimag(data1(ix,iy,iz))
             cov(1:ncov,1)=cov(1:ncov,1)+(/x(1)*x(1),x(1)*x(mult_phi),x(mult_phi)*x(mult_phi)/)*2
             
             x(1)=dble(data2(ix,iy,iz))
             x(mult_comp)=aimag(data2(ix,iy,iz))
             cov(1:ncov,mult_phi)=cov(1:ncov,mult_phi)+(/x(1)*x(1),x(1)*x(mult_phi),x(mult_phi)*x(mult_phi)/)*2
             
          end do
       end do
    end do
    !$OMP END PARALLEL DO
    
    call MPI_ALLREDUCE(cov,cov0,ncov*mult_phi,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
    
    cov0(1:ncov,1:mult_phi)=cov0(1:ncov,1:mult_phi)/ngrid_tot
    
    v=thrm_sigma**2/12 ! <phi phi>
    if(MpiRank==0) then
       write(*,'(A,E)') '   check: var of phi/sigma = ',v
       write(*,'(A,3E)')'         <phi phi>=',cov0(1:ncov,1)
    end if
    v=(thrm_sigma**2*pi*sigma_cr)**2/30 !<phidot phidot>
    if(MpiRank==0) then
       write(*,'(A,E)') '   check: var of phidot/sigma*tau_cr = ',v
       write(*,'(A,3E)')'         <phidot phidot>=',cov0(1:ncov,mult_phi)
    end if
    
  end subroutine check_initial_var
#endif
  
  subroutine get_initial_loop
    logical first_time
    integer(8) iz
    integer info,iter,i,ip
    double precision,parameter :: tol=5d-3
    double precision err,err0,omega,atol,d
    integer,parameter :: max_iter=1000
    double precision,parameter :: crot=0.1d0 ! small rotation angle (+inward displacement) for differentiation
!!$!!$!!$
    integer(8) iz0,iz1,iz2
    !!$!!$!!$logical ex
    character(len=1024) fn
!!$!!$!!$
    if(MpiRank==0) write(*,*) ' #starting from loop configuration'
    call allocate_data(.true.) ! both of data1 and data2 are allocated
    
    do ip=1,2
       
       if(ip==1) then
          d=0d0
       else if(ip==2) then
          d=crot
       end if

#ifdef MPI
       call MPI_BARRIER(MPI_COMM_WORLD,info)
#endif
       
       if(MpiRank==0) then
          write(*,*);write(*,'(A)') ' *over-relaxation:'
       end if
       first_time=.true.
       call init_srj(srj_levels)
       iter=1
       err0=1d0
       
       do
          
          call constraint_loop(d,first_time)
          
          if(MpiRank==0) write(*,'(A,I,F)') '           iter, err :',iter,err0
          if(err0<srj_tol.or.iter>=max_iter) exit
          
          do i=1,srjp%levels-1
             if(srjp%q(i)+.5d0<srjp%beta(i)*iter) exit
          end do
          srjp%q(i)=srjp%q(i)+1
          omega=srjp%omega(i)
          if(MpiRank==0) then
             write(*,'(A,I)') '                 level =',i
             write(*,'(A,F)') '                 omega =',omega
          end if
          
          atol=min(tol,sqrt(err0))
          call over_relaxation(tau_ini,omega,atol,err)
          
#ifdef MPI
          call MPI_ALLREDUCE(err,err0,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,info)
#else
          err0=err
#endif
          iter=iter+1
          
       end do
       
       call fin_srj
       if(MpiRank==0) then
          if(err0<srj_tol) then
             write(*,'(A,F)') ' over-relaxation converged with err:',err0
          else
             write(*,'(A,F)') ' WARNING: over-relaxation failed to converged with err:',err0
          end if
       endif
       
       if(ip==1) then
          
          call allocate_buffer(mult_phi)
          
          write(*,'(A,I,E)') '    temporal buffer size:',MpiRank,byte_dc*alloc_local/1024/1024/1024
#if USETMP==1
          write(fn,*) mpirank
          fn=trim(rootname)//'_phi_'//trim(adjustl(fn))//'.dat.tmp_'
          !!$!!$!!$print *,trim(fn)
          open(199,file=trim(fn),access='sequential',form='unformatted',action='readwrite',iostat=info)
#else
          open(199,access='sequential',form='unformatted',status='scratch',action='readwrite',iostat=info)
#endif
          if(info/=0) call MpiStop('error: cannot open temporal output file')
!!$          inquire(199,exist=ex,name=fn)
!!$          print '(I,L,X,A)',mpirank,ex,trim(fn)


!!$!!$!!$
!!$          iz1=1+ngrid_slice
!!$          iz2=alloc_local+ngrid_slice
!!$          print '(3I,8E20.8)',mpirank,iz1,iz2,dble(ldata2(1)),aimag(ldata2(1)),dble(ldata2(iz1)),aimag(ldata2(iz1)),&
!!$               dble(ldata2(iz2)),aimag(ldata2(iz2)),dble(ldata2(alloc_local_mod)),aimag(ldata2(alloc_local_mod))
!!$          print '(3I,8E20.8)',mpirank,ngrid_slice+1,ngrid_slice*max_0,dble(data2(1,1,1)),aimag(data2(1,1,1)),dble(data2(1,1,min_0)),aimag(data2(1,1,min_0)),&
!!$               dble(data2(ngrid,ngrid,max_0)),aimag(data2(ngrid,ngrid,max_0)),dble(data2(ngrid,ngrid,max_0+1)),aimag(data2(ngrid,ngrid,max_0+1))
!!$          
!!$          iz0=alloc_local+ngrid_slice
!!$          iter=0
!!$          do iz=1,alloc_local,count_dc_max
!!$             iz1=iz+ngrid_slice
!!$             iz2=iz1+count_dc_max-1
!!$             if(.not.(iz2<=iz0)) iz2=iz0
!!$             !!$print '(3I,4E20.8)',mpirank,iz1,iz2,dble(ldata2(iz1)),aimag(ldata2(iz1)),dble(ldata2(iz2)),aimag(ldata2(iz2))
!!$             iter=iter+1
!!$             write(*,'(A,3I)') '    iterations:',MpiRank,iter,(alloc_local-1)/count_dc_max+1
!!$          end do
!!$          call MPI_BARRIER(MPI_COMM_WORLD,info)
!!$          call sleep(15)
!!$          print *
!!$          call MPI_BARRIER(MPI_COMM_WORLD,info)
!!$!!$!!$
!!$          iz1=1+ngrid_slice
!!$          iz2=ngrid_slice*(1+local_0)
!!$          print '(4I)',mpirank,iz1,iz2,alloc_local_mod
!!$          write(199,iostat=info) ldata2(iz1:iz2)
!!$          if(info/=0) call MpiStop('error: cannot write in temporal output file')
!!$!!$!!$          
          iz0=ngrid_slice*(1+local_0)
          iter=0
          do iz=1+ngrid_slice,iz0,count_dc_max
             iz1=min(iz+count_dc_max-1,iz0)
             iter=iter+1
             write(199,iostat=info) ldata2(iz:iz1)
             if(info/=0) call MpiStop('error: cannot write in temporal output file')
             write(*,'(A,3I)') '    iterations:',MpiRank,iter,(iz0-ngrid_slice-1)/count_dc_max+1
          end do
!!$!!$!!$
!!$          iter=0
!!$          do iz=1+ngrid_slice,alloc_local+ngrid_slice,count_dc_max
!!$             iter=iter+1
!!$             write(199,iostat=info) ldata2(iz:min(iz+count_dc_max-1,alloc_local+ngrid_slice))
!!$             if(info/=0) call MpiStop('error: cannot write in temporal output file')
!!$             write(*,'(A,3I)') '    iterations:',MpiRank,iter,(alloc_local-1)/count_dc_max+1
!!$          end do
!!$!!$!!$
          
          if(MpiRank==0) write(*,*) '   phi done'
          
       else if(ip==mult_phi) then
          
          rewind(199)
          
          call allocate_buffer(1)
          
!!$!!$!!$
!!$          iz1=1+ngrid_slice
!!$          iz2=ngrid_slice*(1+local_0)
!!$          print '(4I)',mpirank,iz1,iz2,alloc_local_mod
!!$          read(199,iostat=info) ldata1(iz1:iz2)
!!$          if(info/=0) call MpiStop('error: cannot read from temporal output file')
!!$          ldata2(iz1:iz2)=(ldata2(iz1:iz2)-ldata1(iz1:iz2))*sqrt(loop_velocity**2+loop_velocity_in**2)/crot/(loop_radius*ngo2*dspa)
!!$!!$!!$
          iz0=ngrid_slice*(1+local_0)
          iter=0
          do iz=1+ngrid_slice,iz0,count_dc_max
             iz1=min(iz+count_dc_max-1,iz0)
             iter=iter+1
             read(199,iostat=info) ldata1(iz:iz1)
             if(info/=0) call MpiStop('error: cannot read from temporal output file')
             write(*,'(A,3I)') '    iterations:',MpiRank,iter,(iz0-ngrid_slice-1)/count_dc_max+1
             
             ldata2(iz:iz1)=(ldata2(iz:iz1)-ldata1(iz:iz1))*sqrt(loop_velocity**2+loop_velocity_in**2)/crot/(loop_radius*ngo2*dspa)
          end do
!!$!!$!!$
!!$          iter=0
!!$          do iz=1+ngrid_slice,alloc_local+ngrid_slice,count_dc_max
!!$             iter=iter+1
!!$             read(199,iostat=info) ldata1(iz:min(iz+count_dc_max-1,alloc_local+ngrid_slice))
!!$             if(info/=0) call MpiStop('error: cannot read from temporal output file')
!!$             write(*,'(A,3I)') '    iterations:',MpiRank,iter,(alloc_local-1)/count_dc_max+1
!!$             
!!$             ! conformal time derivative
!!$             ldata2(iz:min(iz+count_dc_max-1,alloc_local+ngrid_slice))=(ldata2(iz:min(iz+count_dc_max-1,alloc_local+ngrid_slice)) &
!!$                  -ldata1(iz:min(iz+count_dc_max-1,alloc_local+ngrid_slice)))*sqrt(loop_velocity**2+loop_velocity_in**2)/crot/(loop_radius*ngo2*dspa)
!!$          end do
!!$!!$!!$

#if USETMP==1
          close(199,status='delete')
#else
          close(199)
#endif
          
          call deallocate_buffer
          if(MpiRank==0) write(*,*) '   phidot done'
          
       end if
       
    end do
    
  end subroutine get_initial_loop
  
  subroutine over_relaxation(tau,omega,tol,err)
    use potential
    double precision,intent(in) :: tau,omega,tol
    double precision,intent(out) :: err
    integer(8) ix,iy,iz,ix1,ix2,iy1,iy2,iz1,iz2
    double precision,parameter :: thrm=0d0 ! assuming zero temperature for loop
    double precision m2,ma2,fact,eps,derr
    complex(C_DOUBLE_COMPLEX) dphi,phi,phi0
    
    integer info,MpiRankD,MpiRankA,isend,irecv, &
         ssend(MPI_STATUS_SIZE),srecv(MPI_STATUS_SIZE)
    integer,parameter :: tag=0
    
    integer count,iter
    double precision,parameter :: eps0=2d-1 ! stepsize for local search
    integer,parameter :: max_iter=20
    
    MpiRankD=modulo(MpiRank-1,MpiSize)
    MpiRankA=modulo(MpiRank+1,MpiSize)
    
    count=ngrid_slice
    
#if DIM==2
    data1(1:ngrid,min_0:max_0)=data2(1:ngrid,min_0:max_0)
#elif DIM==3
    data1(1:ngrid,1:ngrid,min_0:max_0)=data2(1:ngrid,1:ngrid,min_0:max_0)
#endif
    
    ! send lower margin & receive upper margin
#if DIM==2
    call MPI_ISEND(data2(1,min_0),count,MPI_DOUBLE_COMPLEX,MpiRankD,tag,MPI_COMM_WORLD,isend,info)
    call MPI_IRECV(data1(1,max_0+1),count,MPI_DOUBLE_COMPLEX,MpiRankA,tag,MPI_COMM_WORLD,irecv,info)
#elif DIM==3
    call MPI_ISEND(data2(1,1,min_0),count,MPI_DOUBLE_COMPLEX,MpiRankD,tag,MPI_COMM_WORLD,isend,info)
    call MPI_IRECV(data1(1,1,max_0+1),count,MPI_DOUBLE_COMPLEX,MpiRankA,tag,MPI_COMM_WORLD,irecv,info)
#endif
    call MPI_WAIT(isend,ssend,info)
    call MPI_WAIT(irecv,srecv,info)
    
    ! send upper margin & receive lower margin
#if DIM==2
    call MPI_ISEND(data2(1,max_0),count,MPI_DOUBLE_COMPLEX,MpiRankA,tag,MPI_COMM_WORLD,isend,info)
    call MPI_IRECV(data1(1,1),count,MPI_DOUBLE_COMPLEX,MpiRankD,tag,MPI_COMM_WORLD,irecv,info)
#elif DIM==3
    call MPI_ISEND(data2(1,1,max_0),count,MPI_DOUBLE_COMPLEX,MpiRankA,tag,MPI_COMM_WORLD,isend,info)
    call MPI_IRECV(data1(1,1,1),count,MPI_DOUBLE_COMPLEX,MpiRankD,tag,MPI_COMM_WORLD,irecv,info)
#endif
    call MPI_WAIT(isend,ssend,info)
    call MPI_WAIT(irecv,srecv,info)
    
    !!$thrm=thrm_sigma/tau ! temperature (in sigma)
    m2=msquared(thrm)
    ma2=maxion2(thrm)
    fact=1/dspa**2/tau**2/sigma_cr**2
    
    eps=eps0/(fact+abs(msquared(0d0)))
    ! step size in searching for local solution
    ! should be comparable to the curvature of the minimized function
    
    err=0d0
    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ix,iy,iz,ix1,iy1,iz1,ix2,iy2,iz2,phi0,phi,dphi,derr) COLLAPSE(2) REDUCTION(max:err)
    do iz=min_0,max_0
       do iy=1,ngrid
          iz1=iz-1
          iz2=iz+1
          iy1=modulo(iy-2,ngrid)+1
          iy2=modulo(iy,ngrid)+1
#if DIM==2
          
          phi0=data1(iy1,iz)+data1(iy2,iz)+data1(iy,iz1)+data1(iy,iz2)
          phi=phi0/4 ! initial value for local search
          
          do iter=1,max_iter
             dphi=phi0-4*phi
             dphi=dphi*fact-dVdphi(phi,m2,ma2)
             
             if(abs(dphi)<tol) exit
             phi=phi+eps*dphi ! update local phi
          end do
          
          derr=abs(phi-data2(iy,iz))
          err=max(err,derr)
          data2(iy,iz)=(1-omega)*data2(iy,iz)+omega*phi ! update data2
          
#elif DIM==3
          
          do ix=1,ngrid
             ix1=modulo(ix-2,ngrid)+1
             ix2=modulo(ix,ngrid)+1
             phi0=data1(ix1,iy,iz)+data1(ix2,iy,iz)+data1(ix,iy1,iz)+data1(ix,iy2,iz) &
                  +data1(ix,iy,iz1)+data1(ix,iy,iz2) ! initial value for local search
             phi=phi0/6
             
             do iter=1,max_iter
                dphi=phi0-6*phi
                dphi=dphi*fact-dVdphi(phi,m2,ma2)
                
                if(abs(dphi)<tol) exit
                phi=phi+eps*dphi ! update local phi
             end do
             
             derr=abs(phi-data2(ix,iy,iz))
             err=max(err,derr)
             data2(ix,iy,iz)=(1-omega)*data2(ix,iy,iz)+omega*phi ! update data2
             
          end do
#endif
       end do
    end do
    !$OMP END PARALLEL DO
    
  end subroutine over_relaxation
  
  subroutine constraint_loop(d,first_time)
    double precision,intent(in) :: d
    logical,intent(inout) :: first_time
    integer(8) ix,iy,iz,ix0,iy0,iz0
    integer nt,it
    double precision r0,c,s,x0,y0,z0,x,y,z,t,dr,dw
    double precision, parameter :: veps2=1d-6 ! avoiding divide by zero
    
    r0=loop_radius*ngo2*dspa ! radius of loop
    
    dw=d*loop_velocity/sqrt(loop_velocity**2+loop_velocity**2+veps2)
    dr=d*loop_velocity_in/sqrt(loop_velocity**2+loop_velocity**2+veps2)
    
    r0=r0*(1-dr)
    c=cos(dw)
    s=sin(dw)
    
    if(first_time) then
       
       !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(iz,iy,ix,x0,y0,z0,x,y,z) COLLAPSE(2)
       do iz=1+local_0_offset,local_0+local_0_offset
          do iy=1,ngrid
             y0=(iy-ngo2)*dspa
             z0=(iz-ngo2)*dspa
             y= y0*c+z0*s
             z=-y0*s+z0*c
#if DIM==2
             data2(iy,iz-local_0_offset+1)=loop_config(0d0,y,z,r0,width_cr)
#elif DIM==3
             do ix=1,ngrid
                x0=(ix-ngo2)*dspa
                x=x0
                data2(ix,iy,iz-local_0_offset+1)=loop_config(x,y,z,r0,width_cr)
             end do
#endif
          end do
       end do
       !$OMP END PARALLELDO
       
       first_time=.false.
       
    else ! constraint
       
#if DIM==2
       y0=c
       z0=s
       iy0=nint(y0/dspa)+ngo2
       iz0=nint(z0/dspa)+ngo2
       do iz=max(iz0-1,local_0_offset+1),min(iz0+1,local_0+local_0_offset)
          do iy=modulo(iy0-2,ngrid)+1,modulo(iy0,ngrid)+1
             
             y0=(iy-ngo2)*dspa
             z0=(iz-ngo2)*dspa
             y= y0*c+z0*s
             z=-y0*s+z0*c
             data2(iy,iz-local_0_offset+1)=loop_config(0d0,y,z,r0,width_cr)
             
          end do
       end do
       
#elif DIM==3
       
       nt=nint(twopi*r0/dspa) ! probably large enough...
       !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(it,t,x0,y0,z0,ix0,iy0,iz0,iz,iy,ix,x,y,z)
       do it=1,nt
          t=(twopi*it)/nt
          x0=r0*cos(t)
          y0=c*r0*sin(t)
          z0=s*r0*sin(t)
          ix0=nint(x0/dspa)+ngo2
          iy0=nint(y0/dspa)+ngo2
          iz0=nint(z0/dspa)+ngo2
          
          do iz=max(iz0-1,local_0_offset+1),min(iz0+1,local_0+local_0_offset)
             do iy=modulo(iy0-2,ngrid)+1,modulo(iy0,ngrid)+1
                do ix=modulo(ix0-2,ngrid)+1,modulo(ix0,ngrid)+1
                   
                   x0=(ix-ngo2)*dspa
                   y0=(iy-ngo2)*dspa
                   z0=(iz-ngo2)*dspa
                   x=x0
                   y= y0*c+z0*s
                   z=-y0*s+z0*c
                   data2(ix,iy,iz-local_0_offset+1)=loop_config(x,y,z,r0,width_cr)
                end do
             end do
          end do
          
       end do
       !$OMP END PARALLEL DO
       
#endif
       
    end if
    
  end subroutine constraint_loop
  
  function loop_config(x,y,z,r0,width)
    double precision,intent(in) :: x,y,z,r0,width
    complex(C_DOUBLE_COMPLEX) loop_config
    double precision dr,d
    double precision, parameter :: eps=1d-10 ! avoid divide by zero
    double precision, parameter :: c=0.6d0 ! Hagmann & Sikivie 
    
    ! string contribution
    dr=sqrt(x**2+y**2)-r0
    d=sqrt(dr**2+z**2)
    if(d<eps*width) then
       loop_config=cmplx(dr,z)/width*c
    else
       loop_config=tanh(d/width*c)*cmplx(dr,z)/d
    end if
    
    dr=sqrt(x**2+y**2)+r0
    d=sqrt(dr**2+z**2)
    if(d<eps*width) then
       loop_config=loop_config*cmplx(dr,-z)/width*c
    else
       loop_config=loop_config*tanh(d/width*c)*cmplx(dr,-z)/d
    end if
    
  end function loop_config
  
  subroutine init_srj(i)
    integer,intent(in) :: i
    integer info,n
    character(256) fn,str
    double precision omega(max_srj_levels),beta(max_srj_levels)
    
    if(.not.(i>=min_srj_levels.and.i<=max_srj_levels)) &
         call MpiSTop('error: srj levels should be in [min_srj_levels,max_srj_levels]')
    
    write(str,*) i
    fn='param_srj/srj_'//trim(adjustl(str))//'.txt'
    if(MpiRank==0) then
       write(*,'("   SRJ scheme with",3I," levels")') i
       write(*,'("   optimal parameters are read from ",A)') trim(fn)
    end if
    open(101,file=fn,action='read',iostat=info)
    if(info/=0) call MpiStop('error: cannot open param_srj/srj_***.txt')
    
    omega(1:max_srj_levels)=0.d0
    beta(1:max_srj_levels)=0.d0
    do
       read(101,*,iostat=info,err=100) str
       read(str(4:),*,iostat=info,err=100) n
       read(101,*,iostat=info,err=100) omega(1:i)
       read(101,*,iostat=info,err=100) beta(1:i)
       if(n==ngrid) exit
    end do
    
    close(101)
    
    srjp%levels=i
    allocate(srjp%omega(i),srjp%beta(i),srjp%q(i),stat=info)
    srjp%omega(1:i)=omega(1:i)
    srjp%beta(1:i)=beta(1:i)
    srjp%q(1:i)=0
    
    return
    
100 call MpiStop('error: cannot read srj/srj_***.txt')
    
  end subroutine init_srj
  
  subroutine fin_srj
    integer info
    deallocate(srjp%omega,srjp%beta,srjp%q,stat=info)
  end subroutine fin_srj
  
end module initial
