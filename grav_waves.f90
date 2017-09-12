module grav_waves
  use settings
  use potential,only : pi,twopi,fourpi
  implicit none
  
  integer,parameter :: num_bin_grav=25
  
  integer(8),parameter :: num_stress=ndim*(ndim+1)/2,num_tt=num_stress-1-ndim
  double precision,allocatable :: power_grav(:,:),wave_grav(:)
  
#if DIM==3
  complex(C_DOUBLE_COMPLEX),allocatable :: arr_dphi(:,:,:,:),arr_grav(:,:,:,:,:)
#endif
  
  integer(8) ngrid_grav,ngrid_tot_grav,ngo2_grav,ngrid_slice_grav,alloc_local_grav,local_0_grav,local_0_offset_grav
  
  type(C_PTR) plan_backward_grav,cdata_grav1,cdata_grav2
#if DIM==3
  complex(C_DOUBLE_COMPLEX),pointer :: data_grav1(:,:,:),data_grav2(:,:,:)
#endif
  double precision dspa_grav,dwav_grav,dv_spa_grav,dv_wav_grav
  
  logical no_degrade
  
contains
  
  subroutine init_grav
    integer deg,info,iq
    double precision q,dq
    logical no_degrade_tmp
    
    if(MpiRank==0) write(*,*) ' *initialization for gravitational wave calculaiton'
#if DIM==2
    call MpiStop('error: grav_wave.f90 only support DIM==3')
#endif
    
    deg=2**degrade_spec
    ngrid_grav=ngrid/deg
    ngrid_tot_grav=ngrid_grav**ndim
    ngrid_slice_grav=ngrid_tot_grav/ngrid_grav
    ngo2_grav=ngrid_grav/2
    dspa_grav=dspa*deg
    dwav_grav=twopi/dspa_grav/ngrid_grav
    dv_spa_grav=dspa_grav**ndim
    dv_wav_grav=(dwav_grav/twopi)**ndim
    
    call init_fftw_grav
    
    no_degrade_tmp=(degrade_spec==0.and.local_0_grav==local_0.and.local_0_offset_grav==local_0_offset)
    call MPI_ALLREDUCE(no_degrade_tmp,no_degrade,1,MPI_LOGICAL,MPI_LAND,MPI_COMM_WORLD,info)
    if(MPIRank==0) then
       if(no_degrade) then
          write(*,*) '  no degradation in grav wave calculation'
       else
          write(*,*) '  degradation in grav wave calculation'
          write(*,*) '  degrade:',degrade_spec
       end if
    end if
    
    allocate(power_grav(num_bin_grav,n_chk),wave_grav(0:num_bin_grav),stat=info)
    if(info/=0) call MpiStop('error: cannot allocate power_grav, wave_grav')
    q=pi/dspa_grav ! maximum q
    dq=q/num_bin_grav
    if(MpiRank==0) then
       write(*,'(A,E)') '  bin size of wav_grav (delta_k*tau_cr/4pi): ',dq/fourpi
       write(*,'(A,E)') '  bin size of wav_grav (delta_k*tau_fin/4pi): ',dq/fourpi*tau_fin
    end if
    wave_grav(0)=0d0
    do iq=1,num_bin_grav
       wave_grav(iq)=wave_grav(iq-1)+dq
    end do
    
#if DIM==3
    allocate(arr_grav(2,num_tt,ngrid_grav,ngrid_grav,local_0_grav),arr_dphi(ndim,ngrid_grav,ngrid_grav,local_0_grav),stat=info)
#endif    
    if(info/=0) call MpiStop('error: cannot allocate arr_grav or arr_dphi')
    arr_grav=0d0
    
    call MPI_BARRIER(MPI_COMM_WORLD,info)
    if(MpiRank==0) write(*,*) ' *done'
    
  end subroutine init_grav
  
  subroutine fin_grav
    integer info
    
    call fin_fftw_grav
    deallocate(power_grav,wave_grav,arr_grav,arr_dphi,stat=info)
    
  end subroutine fin_grav
  
  subroutine init_fftw_grav
#include <fftw3-mpi.f03>
    integer info
    
#if DIM==3
    alloc_local_grav=FFTW_MPI_LOCAL_SIZE_3D(ngrid_grav,ngrid_grav, &
         ngrid_grav,MPI_COMM_WORLD,local_0_grav,local_0_offset_grav)
#endif
    
    write(*,'(A,5I)') ' information of data_grav:',MpiRank,ngrid_grav,local_0_grav,local_0_offset_grav,alloc_local_grav
    
    cdata_grav1=FFTW_ALLOC_COMPLEX(alloc_local_grav)
    cdata_grav2=FFTW_ALLOC_COMPLEX(alloc_local_grav)
    if(.not.c_associated(cdata_grav1)) call MpiStop('error: cdata_grav1 is null')
    if(.not.c_associated(cdata_grav2)) call MpiStop('error: cdata_grav2 is null')
    
#if DIM==3
    call c_f_pointer(cdata_grav1,data_grav1,[ngrid_grav,ngrid_grav,local_0_grav])
    call c_f_pointer(cdata_grav2,data_grav2,[ngrid_grav,ngrid_grav,local_0_grav])
#endif
    if(.not.(associated(data_grav1).and.associated(data_grav2))) &
         call MpiStop('error: data_grav1 or data_grav2 is not allocated')
    
    call FFTW_PLAN_WITH_NTHREADS(int(num_threads,kind=4))
#if DIM==3
    plan_backward_grav=FFTW_MPI_PLAN_DFT_3D(ngrid_grav,ngrid_grav,ngrid_grav, &
         data_grav2,data_grav1,MPI_COMM_WORLD,FFTW_BACKWARD,FFTW_ESTIMATE)
#endif
    if(.not.c_associated(plan_backward_grav)) call MpiStop('error: plan_backward_grav is null')
    
    call memory_check_mkl()
    
  end subroutine init_fftw_grav
  
  subroutine fin_fftw_grav
#include <fftw3-mpi.f03>
    
    call FFTW_DESTROY_PLAN(plan_backward_grav)
    
    data_grav1=>NULL()
    data_grav2=>NULL()
    call FFTW_FREE(cdata_grav1)
    call FFTW_FREE(cdata_grav2)
    
  end subroutine fin_fftw_grav
  
!!$  subroutine get_comwav_grav(ixyz,q,qhat)
!!$    integer(8),intent(in) :: ixyz(ndim)
!!$    double precision,intent(out) :: q
!!$    double precision,intent(out),optional :: qhat(ndim)
!!$    integer iq(ndim)
!!$    
!!$    iq(1:ndim)=mod(ixyz(1:ndim)+ngo2_grav-1,ngrid_grav)-ngo2_grav
!!$    q=sqrt(dble(sum(iq(1:ndim)**2)))
!!$    if(present(qhat)) then
!!$       if(q>0d0) then
!!$          qhat(1:ndim)=iq(1:ndim)/q
!!$       else
!!$          qhat(1:ndim)=0d0
!!$       end if
!!$    end if
!!$    q=q*dwav_grav
!!$    
!!$  end subroutine get_comwav_grav
  
  subroutine get_comwav_grav(ixyz,q,pol,i,j)
    integer(8),intent(in) :: ixyz(ndim)
    double precision,intent(out) :: q
    integer,intent(in),optional :: i,j
    complex(C_DOUBLE_COMPLEX),intent(out),optional :: pol ! multiplicity-weighted ij component of polarization tensor with positive helicity(++)
    double precision qhat(ndim),c,s,r
    complex(C_DOUBLE_COMPLEX) p(ndim)
    
    qhat(1:ndim)=dble(mod(ixyz(1:ndim)+ngo2_grav-1,ngrid_grav)-ngo2_grav)
    q=sqrt(sum(qhat(1:ndim)**2))
    
    if(present(pol)) then
       if(q>0d0) then
          qhat(1:ndim)=qhat(1:ndim)/q
          r=sqrt(qhat(1)**2+qhat(2)**2)
          if(r>0d0) then
             c=qhat(1)/r
             s=qhat(2)/r
          else
             c=1d0
             s=0d0
          end if
          p(1)=cmplx(qhat(3)*c,-s)
          p(2)=cmplx(qhat(3)*s,c)
          p(3)=cmplx(-r,0d0)
          pol=p(i)*p(j)
          if(i==j) pol=pol/2
       else
          pol=0d0
       end if
    end if
    
    q=q*dwav_grav
    
  end subroutine get_comwav_grav
  
  subroutine grav_wave_integration(tau,dtau)
#include <fftw3-mpi.f03>
    double precision,intent(in) :: tau,dtau
    double precision dw,q,c,s
    integer(8) ix,iy,iz,ixyz(ndim)
    integer info,i,j
    complex(C_DOUBLE_COMPLEX) polp,poln
    
    if(MpiRank==0) write(*,'(A)') '   integrating gravitational wave ...'
    
    if(no_degrade) then
       call data1_to_arr_dphi
    else
       call data1_to_arr_dphi_degrade
    end if
    
    dw=tau*dtau*dv_spa_grav*sigma_cr**2
    do i=1,ndim
       do j=1,i
          
          !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(iz,iy,ix) COLLAPSE(ndim)
          do iz=1,local_0_grav
             do iy=1,ngrid_grav
#if DIM==3
                do ix=1,ngrid_grav
                   data_grav2(ix,iy,iz)=arr_dphi(i,ix,iy,iz)*conjg(arr_dphi(j,ix,iy,iz))+conjg(arr_dphi(i,ix,iy,iz))*arr_dphi(j,ix,iy,iz)
                end do
#endif
             end do
          end do
          !$END OMP PARALLEL DO
          
          call FFTW_MPI_EXECUTE_DFT(plan_backward_grav,data_grav2,data_grav1)
          
          !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(iz,iy,ix,ixyz,q,c,s,polp,poln) COLLAPSE(ndim)
          do iz=1,local_0_grav
             do iy=1,ngrid_grav
#if DIM==3
                do ix=1,ngrid_grav
                   ixyz(1:ndim)=(/ix,iy,iz+local_0_offset_grav/)
                   call get_comwav_grav(ixyz,q,polp,i,j)
                   q=q*tau
                   c=cos(q)
                   s=sin(q)
                   polp=polp*dw
                   poln=conjg(polp)
                   arr_grav(1,1,ix,iy,iz)=arr_grav(1,1,ix,iy,iz)+c*polp*data_grav1(ix,iy,iz)
                   arr_grav(2,1,ix,iy,iz)=arr_grav(1,2,ix,iy,iz)+s*polp*data_grav1(ix,iy,iz)
                   arr_grav(1,2,ix,iy,iz)=arr_grav(2,1,ix,iy,iz)+c*poln*data_grav1(ix,iy,iz)
                   arr_grav(2,2,ix,iy,iz)=arr_grav(2,2,ix,iy,iz)+s*poln*data_grav1(ix,iy,iz)
                end do
#endif
             end do
          end do
          !$OMP END PARALLEL DO
          
       end do
    end do
    
    call MPI_BARRIER(MPI_COMM_WORLD,info)
    if(MpiRank==0) write(*,'(A)') '   ... done'
    
  end subroutine grav_wave_integration
  
  subroutine grav_spec_estimation(idx)
    use potential
    integer,intent(in) :: idx
    integer(8) ix,iy,iz,ixyz(ndim)
    integer iq,info
    double precision q,numer(num_bin_grav),denom(num_bin_grav),numer0(num_bin_grav),denom0(num_bin_grav),dq,freq_cr_in_Hz,fact
    character(len=charlen) str,file_grav
    
    dq=wave_grav(1)-wave_grav(0)
    freq_cr_in_Hz=Mstar_in_Hz/tau_cr/a_now/twopi
    fact=1/tau_cr**6/BigH_in_Mstar**2/a_now**4/3/twopi**2*dv_wav_grav
    
    numer0(1:num_bin_grav)=0d0
    denom0(1:num_bin_grav)=0d0
    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(iz,iy,ix,ixyz,q,iq) REDUCTION(+:numer,denom) COLLAPSE(ndim)
    do iz=1,local_0_grav
       do iy=1,ngrid_grav
#if DIM==3
          do ix=1,ngrid_grav
             
             ixyz(1:ndim)=(/ix,iy,iz+local_0_offset_grav/)
             call get_comwav_grav(ixyz,q)
             iq=nint(q/dq+0.5d0)
             if(q>0d0.and.iq>=1.and.iq<=num_bin_grav) then
                
                denom0(iq)=denom0(iq)+1
                numer0(iq)=numer0(iq)+abs(arr_grav(1,1,ix,iy,iz))**2+abs(arr_grav(2,1,ix,iy,iz))**2 &
                     +abs(arr_grav(1,2,ix,iy,iz))**2+abs(arr_grav(2,2,ix,iy,iz))**2
                
             end if
             
          end do
#endif
       end do
    end do
    !$END PARALLEL DO
    
    call MPI_ALLREDUCE(denom0,denom,num_bin_grav,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
    call MPI_ALLREDUCE(numer0,numer,num_bin_grav,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
    
    do iq=1,num_bin_grav
       q=(wave_grav(iq)+wave_grav(iq-1))/2
       if(MpiRank==0) print '(I,3E20.4)',iq,q,numer(iq),denom(iq)
       if(denom(iq)>0d0) then
          power_grav(iq,idx)=numer(iq)/denom(iq)*q**3*fact
       else
          power_grav(iq,idx)=0d0
       end if
    end do
    
    if(MpiRank==0) then
       write(str,*) idx
       file_grav=trim(rootname)//'_gravspec_'//trim(adjustl(str))//'.dat'
       open(100,file=trim(file_grav),action='write',iostat=info)
       do iq=1,num_bin_grav
          q=(wave_grav(iq)+wave_grav(iq-1))/2
          write(100,'(5E20.10e3)',iostat=info) tau_id(idx),q,q*freq_cr_in_Hz, &
               power_grav(iq,idx),denom(iq)
       end do
       close(100)
    end if
    
    call MPI_BARRIER(MPI_COMM_WORLD,info)
    
  end subroutine grav_spec_estimation
  
  subroutine data1_to_arr_dphi
    integer(8) min_1,max_1,ix0,iy0,iz0,ix1,iy1,iz1,ix2,iy2,iz2
    double precision d
    integer info
    
    call MPI_BARRIER(MPI_COMM_WORLD,info)
    if(MpiRank==0) write(*,*) '  undegraded round phi from data1 to arr_dphi ...'
    
    d=dspa*2
    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(iz0,iy0,ix0,iz1,iy1,ix1,iz2,iy2,ix2) COLLAPSE(2)
    do iz0=min_0,max_0
       do iy0=1,ngrid
          iz1=iz0-1
          iz2=iz0+1
          iy1=modulo(iy0-2,ngrid)+1
          iy2=modulo(iy0,ngrid)+1
#if DIM==3
          do ix0=1,ngrid
             ix1=modulo(ix0-2,ngrid)+1
             ix2=modulo(ix0,ngrid)+1
             arr_dphi(1,ix0,iy0,iz1)=data1(ix2,iy0,iz0)-data1(ix1,iy0,iz0)
             arr_dphi(2,ix0,iy0,iz1)=data1(ix0,iy2,iz0)-data1(ix0,iy1,iz0)
             arr_dphi(3,ix0,iy0,iz1)=data1(ix0,iy0,iz2)-data1(ix0,iy0,iz1)
             arr_dphi(1:ndim,ix0,iy0,iz1)=arr_dphi(1:ndim,ix0,iy0,iz1)/d
          end do
#endif
       end do
    end do
    !$OMP END PARALLEL DO
    
    call MPI_BARRIER(MPI_COMM_WORLD,info)
    if(MpiRank==0) write(*,*) '  ... done'
    
  end subroutine data1_to_arr_dphi

  subroutine data1_to_arr_dphi_degrade
    integer,parameter :: tg=0
    integer info,deg,MpiRankD,MpiRankA,count,ssend(MPI_STATUS_SIZE),srecv(MPI_STATUS_SIZE),irank,irecv0,isend0
    integer(8),allocatable :: min_1t(:),max_1t(:),local_1t(:),local_1t_offset(:)
    integer,allocatable :: irecv(:),isend(:)
#if DIM==3
    complex(kind(0d0)),allocatable :: data_dphi_degrade(:,:,:,:),data_dphi_tmp(:,:,:)
#endif
    complex(kind(0d0)) xc
    integer(8) ix,iy,iz,ix0,iy0,iz0,j,k,min_1,max_1,ix1,ix2,iy1,iy2,iz1,iz2
    double precision d
    
    call MPI_BARRIER(MPI_COMM_WORLD,info)
    if(MpiRank==0) write(*,*) '  degrading and distributing round phi from data1 to arr_dphi ...'
    deg=2**degrade_spec
    min_1=local_0_offset/deg+1
    max_1=(local_0_offset+local_0-1)/deg+1
    !!write(*,'(A,4I)') '  information of data_dphi_degrade:',MpiRank,deg,min_1,max_1
    
#if DIM==3
    allocate(data_dphi_degrade(ndim,ngrid_grav,ngrid_grav,min_1:max_1),data_dphi_tmp(ndim,ngrid_grav,ngrid_grav),stat=info)
#endif
    if(info/=0) call MpiStop('error: cannot allocate data_dphi_degrade, data_dphi_tmp')
    
    ! degraded spatial derivative
    xc=cmplx(0d0)
    d=deg**ndim*dspa*2
    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(iz,iy,ix,iz0,iy0,ix0,iz2,iy2,ix2,iz1,iy1,ix1) COLLAPSE(ndim)
    do iz=min_1,max_1
       do iy=1,ngrid_grav
#if DIM==3
          do ix=1,ngrid_grav
             data_dphi_degrade(1:ndim,ix,iy,iz)=xc
             do iz0=max(1,(iz-1)*deg+1-local_0_offset)+1,min(local_0,iz*deg-local_0_offset)+1
                iz1=iz0-1
                iz2=iz0+1
                do iy0=(iy-1)*deg+1,iy*deg
                   iy1=modulo(iy0-2,ngrid)+1
                   iy2=modulo(iy0,ngrid)+1
                   do ix0=(ix-1)*deg+1,ix*deg
                      ix1=modulo(ix0-2,ngrid)+1
                      ix2=modulo(ix0,ngrid)+1
                      data_dphi_degrade(1,ix,iy,iz)=data_dphi_degrade(1,ix,iy,iz)+data1(ix2,iy0,iz0)-data1(ix1,iy0,iz0)
                      data_dphi_degrade(2,ix,iy,iz)=data_dphi_degrade(2,ix,iy,iz)+data1(ix0,iy2,iz0)-data1(ix0,iy1,iz0)
                      data_dphi_degrade(3,ix,iy,iz)=data_dphi_degrade(3,ix,iy,iz)+data1(ix0,iy0,iz2)-data1(ix0,iy0,iz1)
                   end do
                end do
             end do
             data_dphi_degrade(1:ndim,ix,iy,iz)=data_dphi_degrade(1:ndim,ix,iy,iz)/d
          end do
#endif
       end do
    end do
    !$OMP END PARALLEL DO
    
    count=ngrid_slice_grav*ndim
    if(mod(local_0+local_0_offset,deg)/=0) then
       MpiRankA=modulo(MpiRank+1,MpiSize)
#if DIM==3
       call MPI_ISEND(data_dphi_degrade(1,1,1,max_1),count,MPI_DOUBLE_COMPLEX,MpiRankA,tg,MPI_COMM_WORLD,isend0,info)
#endif
    end if
    if(mod(local_0_offset,deg)/=0) then
       MpiRankD=modulo(MpiRank-1,MpiSize)
       call MPI_IRECV(data_dphi_tmp,count,MPI_DOUBLE_COMPLEX,MpiRankD,tg,MPI_COMM_WORLD,irecv0,info)
    end if
    
    if(mod(local_0+local_0_offset,deg)/=0) call MPI_WAIT(isend0,ssend,info)
    if(mod(local_0_offset,deg)/=0) then
       call MPI_WAIT(irecv0,ssend,info)
       !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(iy0)
       do iy0=1,ngrid_grav
#if DIM==3
          data_dphi_degrade(1:ndim,1:ngrid_grav,iy0,min_1)=data_dphi_degrade(1:ndim,1:ngrid_grav,iy0,min_1)+data_dphi_tmp(1:ndim,1:ngrid_grav,iy0)
#endif
       end do
       !$OMP END PARALLEL DO
    end if
    
    deallocate(data_dphi_tmp,stat=info)
    
    ! data_dphi_degrade to arr_dphi
    allocate(local_1t(0:MpiSize-1),local_1t_offset(0:MpiSize-1),min_1t(0:MpiSize-1),max_1t(0:MpiSize-1), &
         irecv(0:MpiSize-1),isend(0:MpiSize-1),stat=info)
    if(info/=0) call MpiStop('error: cannot allocate local_1t, ...')
    call MPI_ALLGATHER(min_1,1,MPI_INTEGER8,min_1t,1,MPI_INTEGER8,MPI_COMM_WORLD,info)
    call MPI_ALLGATHER(max_1,1,MPI_INTEGER8,max_1t,1,MPI_INTEGER8,MPI_COMM_WORLD,info)
    call MPI_ALLGATHER(local_0_grav,1,MPI_INTEGER8,local_1t,1,MPI_INTEGER8,MPI_COMM_WORLD,info)
    call MPI_ALLGATHER(local_0_offset_grav,1,MPI_INTEGER8,local_1t_offset,1,MPI_INTEGER8,MPI_COMM_WORLD,info)
    
    do irank=0,MpiSize-1
       
       if(irank==MpiRank) cycle
       
       ! send
       j=max(min_1,local_1t_offset(irank)+1)
       k=min(max_1,local_1t_offset(irank)+local_1t(irank))
       if(k>=j) then
          count=(k-j+1)*ngrid_slice_grav*ndim
#if DIM==3
          call MPI_ISEND(data_dphi_degrade(1,1,1,j),count,MPI_DOUBLE_COMPLEX,irank,tg,MPI_COMM_WORLD,isend(irank),info)
#endif
       end if
       
       ! receive
       j=max(min_1t(irank),local_0_offset_grav+1)
       k=min(max_1t(irank),local_0_offset_grav+local_0_grav)
       if(k>=j) then
          iz=j-local_0_offset_grav
          count=(k-j+1)*ngrid_slice_grav*ndim
#if DIM==3
          call MPI_IRECV(arr_dphi(1,1,1,iz),count,MPI_DOUBLE_COMPLEX,irank,tg,MPI_COMM_WORLD,irecv(irank),info)
#endif
       end if
       
    end do
    
    ! within node
    j=max(min_1,local_0_offset_grav+1)
    k=min(max_1,local_0_offset_grav+local_0_grav)
    if(k>=j) then
       !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(iz,iy) COLLAPSE(2)
       do iz=j-local_0_offset_grav,k-local_0_offset_grav
          do iy=1,ngrid_grav
#if DIM==3
             arr_dphi(1:ndim,1:ngrid_grav,iy,iz)=data_dphi_degrade(1:ndim,1:ngrid_grav,iy,iz+local_0_offset_grav)
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
       
       j=max(min_1t(irank),local_0_offset_grav+1)
       k=min(max_1t(irank),local_0_offset_grav+local_0_grav)
       if(k>=j) call MPI_WAIT(irecv(irank),srecv,info)
    end do
    
    deallocate(data_dphi_degrade,min_1t,max_1t,local_1t,local_1t_offset,irecv,isend,stat=info)
    
    call MPI_BARRIER(MPI_COMM_WORLD,info)
    if(MpiRank==0) write(*,*) '  ... done'
    
  end subroutine data1_to_arr_dphi_degrade
  
end module grav_waves
