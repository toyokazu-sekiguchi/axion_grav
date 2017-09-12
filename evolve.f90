module evolve
  use settings
  use subroutines
  implicit none
  
  double precision,allocatable :: array_sw(:,:),array_crho(:,:)
  
contains
  
  subroutine get_dtau(tau1,tau2,dtau,n)
    use potential,only : thrm_sigma
    double precision,intent(in) :: tau1,tau2
    double precision,intent(out) :: dtau
    integer,intent(out) :: n
    double precision omega,x
    
    omega=thrm_sigma*tau2*sigma_cr ! R*mass (in 1/tau_cr) at tau_fin >> tau_cr
    omega=sqrt((dwav*ngrid)**2+omega**2) ! maximum comoving wave number 
    dtau=eps_dtau/omega
    n=ceiling((tau2-tau1)/dtau)
    dtau=(tau2-tau1)/n
    if(MpiRank==0) write(*,'(A,E,I)') '   dtau & number of timesteps: ',dtau,n
    
  end subroutine get_dtau
  
  subroutine evolution(chk)
    use potential
    logical,intent(in) :: chk
    integer idx,info
    double precision tau_start,tau_end,dtau_chk,tau0
    double precision,parameter :: tau_at_pt=1.d0
    logical flag_defects
    
    if(MpiRank==0) then
       write(*,*)
       write(*,'(A)') ' #evolution'
    end if
    call init_evolve
    
#if DIM==2
    allocate(array_sw(4,0:n_chk),stat=info)
#elif DIM==3
    allocate(array_sw(11,0:n_chk),array_crho(3,0:n_chk),stat=info)
#endif
    if(info/=0) call MpiStop('error: cannot allocate array_sw (and array_crho for DIM==3)')
    
    tau0=max(tau_ini,tau_at_pt)
    dtau_chk=(tau_fin-tau0)/n_chk
    if(tau_ini+dtau_chk<tau_at_pt) call MpiStop(' error: too large n_chk; tau at idx=1 is smaller than tau_at_pt')
    
    if(spec_est.or.grav_spec) then
       tau_id(0)=tau_ini
       do idx=1,n_chk
          tau_id(idx)=tau0+dtau_chk*idx
       end do
    end if
    
    if(chk) then
       call read_chk(tau_end,idx,flag_defects)
    else
       tau_end=tau_ini
       idx=0
       flag_defects=.true.
       
       call notice_walltime
       if((tsim_now-tsim_begin)*(1+1d0/max(idx,1))>tsim_max) then
          if(MpiRank==0) write(*,'(A,I3)') ' *update chk files: idx = ',idx
          call write_chk(tau_end,idx,flag_defects)
       end if
    end if
    
    if(idx==0) then
       open(300,file=trim(rootname)//'_string_wall_parameters.txt',position='rewind',iostat=info)
    else
       open(300,file=trim(rootname)//'_string_wall_parameters.txt',position='append',iostat=info)
    end if
#if DIM==3
    if(idx==0) then
       open(400,file=trim(rootname)//'_crho.txt',position='rewind',iostat=info)
    else
       open(400,file=trim(rootname)//'_crho.txt',position='append',iostat=info)
    end if
#endif
    
    if(idx+n_id>n_chk.and.flag_defects) then
       call defects_and_spectrum(idx,tau_end)
#if DIM==2
       write(300,'(I,4E20.10e3)') idx,array_sw(1:4,idx)
#elif DIM==3
       write(300,'(I,11E20.10e3)') idx,array_sw(1:11,idx)
       write(400,'(I,3E20.10e3)') idx,array_crho(1:3,idx)
#endif
       flag_defects=.false.
       call write_chk(tau_end,idx,flag_defects)
    end if
    if(info/=0) call MpiStop('error: cannot open file string_wall_parameters.txt')
    
    do while(idx<n_chk)
       
       tau_start=tau_end
       idx=idx+1
       tau_end=tau0+dtau_chk*idx
       
       if(MpiRank==0) then
          write(*,*)
          write(*,'(A,I3)') ' #idx = ',idx
          write(*,'(A,E,A,E)') '  tau_start = ',tau_start,', tau_end = ',tau_end
          write(*,'(A,E)')     '  zeta_start = ',zeta*tau_start**2
          write(*,'(A,F,E)') '   tau_start, m_axion/H:',tau_start,sqrt(maxion2(thrm_sigma/tau_start))*tau_start**2*sigma_cr
          write(*,'(A,F,E)') '   tau_end, m_axion/H:',tau_end,sqrt(maxion2(thrm_sigma/tau_end))*tau_end**2*sigma_cr
       end if
       
       call symplectic_integration(tau_start,tau_end)
       
       call notice_walltime
       if((tsim_now-tsim_begin)*(1+1d0/idx)>tsim_max) then
          if(MpiRank==0) write(*,'(A,I3)') ' *update chk files: idx = ',idx
          flag_defects=.true.
          call write_chk(tau_end,idx,flag_defects)
       end if
       
       if(idx+n_id>n_chk) then
          call defects_and_spectrum(idx,tau_end)
#if DIM==2
          write(300,'(I,4E20.10e3)') idx,array_sw(1:4,idx)
#elif DIM==3
          write(300,'(I,11E20.10e3)') idx,array_sw(1:11,idx)
          write(400,'(I,3E20.10e3)') idx,array_crho(1:3,idx)
#endif
          if(tsim_now-tsim_begin>tsim_max) then
             flag_defects=.false.
             call write_chk(tau_end,idx,flag_defects)
          end if
       end if
       
    end do
    
    deallocate(array_sw,stat=info)
#if DIM==3
    deallocate(array_crho,stat=info)
#endif
    
    flag_defects=.false.
    call write_chk(tau_end,idx,flag_defects)
    close(300)
#if DIM==3
    close(400)
#endif
    
    call fin_evolve
    
  end subroutine evolution
  
  subroutine defects_and_spectrum(idx,tau)
#if DIM==2
!!$!!$    use defects_2d
#elif DIM==3
    use defects_3d
#endif
    use spectrum
    use grav_waves
    use potential
    integer,intent(in) :: idx
    double precision,intent(in) :: tau
#if DIM==2
    integer num_strings
    double precision xi,largea
#elif DIM==3
    double precision xi(3),largea(3),vel(3)
#endif
    
    call update_margin(mult_phi) ! phidot
    
    ! separate strings
#if DIM==2
    xi=0d0
    largea=0d0
!!$!!$    call id_defects_2d(idx,thickness_cr/tau/dspa,num_strings,largea)
!!$!!$    xi=dble(num_strings)
!!$!!$    xi=(tau/dspa)**2*xi/ngrid_tot/4 ! string parameter
!!$!!$    largea=(tau/dspa)*largea/ngrid_tot/2 ! wall parameter
    if(MpiRank==0) then
       write(*,'(A,F,E)')'  tau, string parameter: ',tau,xi
       write(*,'(A,F,E)') '  tau, wall parameter: ',tau,largea
       write(*,'(A,F,E)') '  tau, wall-string energy ratio:', &
            tau,eratio_wall_string(xi,thrm_sigma/tau,tau**2*sigma_cr/2)
    end if
    array_sw(1:4,idx)=(/tau,tau**2*zeta,xi,largea/)
#elif DIM==3
    call id_defects(idx,tau,xi,largea,vel) ! xi is the string length in units of dspa; largea is the wall area in units of dspa^2
    xi(1:3)=(tau/dspa)**2*xi(1:3)/ngrid_tot/4 ! string parameter
    largea(1:3)=(tau/dspa)*largea(1:3)/ngrid_tot/2 ! wall parameter
    if(MpiRank==0) then
       if(id_strings) then 
          write(*,'(A,F,3E)')'  tau, string parameter: ',tau,xi(1:3)
          write(*,'(A,F,3E)')'  tau, string velocity parameter: ',tau,vel(1:3)
       end if
       if(id_walls) write(*,'(A,F,3E)') '  tau, wall parameter: ',tau,largea(1:3)
       if(id_strings.and.id_walls) write(*,'(A,F,E)') '  tau, wall-string energy ratio:', &
            tau,eratio_wall_string(xi(1),thrm_sigma/tau,tau**2*sigma_cr/2)
    end if
    array_sw(1:11,idx)=(/tau,tau**2*zeta,xi(1:3),largea(1:3),vel(1:3)/)
#endif
    if(spec_est) call spectrum_estimation(idx)
    if(grav_spec) call grav_spec_estimation(idx)
    
#if DIM==3
    array_crho(1:3,idx)=0d0
    !if(spec_est.and.group_strings) call analyze_energy(idx)
    if(spec_est) call analyze_energy(idx)
#endif
    
  end subroutine defects_and_spectrum
  
#if DIM==3
  subroutine analyze_energy(idx)
    use defects_3d
    use spectrum
    integer,intent(in) :: idx
    double precision tau,crho0,crho1,crho2,xi0,gam,v2,xi1
    
    tau=array_sw(1,idx)
    xi0=array_sw(5,idx)
    xi1=array_sw(3,idx)-array_sw(5,idx)
    gam=array_sw(11,idx)
    v2=array_sw(10,idx)
    
    crho0=xi0*twopi*log(0.5d0*tau**2/sqrt(xi0)/width_cr)*4*gam
    crho1=xi1*twopi*log(0.5d0*tau**2/sqrt(xi0)/width_cr)*4*gam
    crho2=crho_spec(1,idx)
    if(MpiRank==0) write(*,'(A,4E)') ' conformal comoving energy densities: ',tau,crho0,crho1,crho2
    
    array_crho(1:3,idx)=(/crho0,crho1,crho2/)
    
  end subroutine analyze_energy
#endif
  
  subroutine symplectic_integration(t1,t2)
    double precision,intent(in) :: t1,t2 ! not proper, but conformal time
    double precision t,dt ! not proper, but conformal time
    integer it,nt
    
    call get_dtau(t1,t2,dt,nt)
    t=t1
    
    do it=1,nt
       if(it==1) call evolve_phi(t,dt/2)
       call evolve_phidot(t,dt)
       if(it==nt) then
          call evolve_phi(t,dt/2)
       else
          call evolve_phi(t,dt)
       end if
       if(MpiRank==0) write(*,'(A,I)') '  iterations: ',it
    end do
    
  end subroutine symplectic_integration
  
  subroutine evolve_phi(t,dt)
    double precision,intent(inout) :: t ! not proper, but conformal time
    double precision,intent(in) :: dt ! not proper, but conformal time
    integer(8) iz,iy,ix
    
    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(iz,iy,ix) COLLAPSE(2)
    do iz=min_0,max_0
       do iy=1,ngrid
#if DIM==2
          data1(iy,iz)=data1(iy,iz)+dt*data2(iy,iz)
#elif DIM==3
          do ix=1,ngrid
             data1(ix,iy,iz)=data1(ix,iy,iz)+dt*data2(ix,iy,iz) 
          end do
#endif
       end do
    end do
    !$OMP END PARALLEL DO
    
    call update_margin(1)
    t=t+dt
    
  end subroutine evolve_phi
  
  subroutine evolve_phidot(t,dt)
    use potential 
    use grav_waves
    double precision,intent(in) :: t,dt ! not proper, but conformal time
    integer(8) ix,iy,iz,ix1,iy1,iz1,ix2,iy2,iz2
    double precision hubble,thrm,fact0,fact1,fact2,fact3,m2,ma2
    complex(kind(0d0)) dlap,ddv
    
    hubble=1/t ! conformal Hubble (in 1/tau_cr)
    if(initial_loop) then
       thrm=0d0
    else
       thrm=thrm_sigma/t ! temperature (in sigma)
    end if
    m2=msquared(thrm)
    ma2=maxion2(thrm)
    fact0=1+dt*hubble
    fact1=1-dt*hubble
    fact2=dt/dspa**2
    fact3=-dt*t**2*sigma_cr**2
    fact1=fact1/fact0
    fact2=fact2/fact0
    fact3=fact3/fact0
    
    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(iz,iz1,iz2,iy,iy1,iy2,ix,ix1,ix2,dlap,ddv) COLLAPSE(2)
    do iz=min_0,max_0
       do iy=1,ngrid
          iz1=iz-1
          iz2=iz+1
          iy1=modulo(iy-2,ngrid)+1
          iy2=modulo(iy,ngrid)+1
#if DIM==2
          dlap=data1(iy1,iz)+data1(iy2,iz) &
               +data1(iy,iz1)+data1(iy,iz2) &
               -4*data1(iy,iz)
          ddv=dVdphi(data1(iy,iz),m2,ma2)
          data2(iy,iz)=fact1*data2(iy,iz)+fact2*dlap+fact3*ddv
#elif DIM==3
          do ix=1,ngrid
             ix1=modulo(ix-2,ngrid)+1
             ix2=modulo(ix,ngrid)+1
             
             dlap=data1(ix1,iy,iz)+data1(ix2,iy,iz) &
                  +data1(ix,iy1,iz)+data1(ix,iy2,iz) &
                  +data1(ix,iy,iz1)+data1(ix,iy,iz2) &
                  -6*data1(ix,iy,iz)
             ddv=dVdphi(data1(ix,iy,iz),m2,ma2)
             data2(ix,iy,iz)=fact1*data2(ix,iy,iz)+fact2*dlap+fact3*ddv
          end do
#endif
       end do
    end do
    !$OMP END PARALLEL DO
    
    if(grav_spec) call grav_wave_integration(t,dt)
    
  end subroutine evolve_phidot
  
  subroutine update_margin(ip0)
    integer,intent(in) :: ip0
    integer info,MpiRankD,MpiRankA,ssend(MPI_STATUS_SIZE),srecv(MPI_STATUS_SIZE),isend,irecv,count,ip
    integer,parameter :: tag=0
    
    MpiRankD=modulo(MpiRank-1,MpiSize)
    MpiRankA=modulo(MpiRank+1,MpiSize)
    
    ip=modulo(ip0-1,mult_phi)+1
    count=ngrid_slice
    if(count/ngrid/=ngrid**(ndim-2)) call MpiStop('error: may be count is too large')

#if DIM==2
    if(ip==1) then
       call MPI_ISEND(data1(1,min_0),count,MPI_DOUBLE_COMPLEX,MpiRankD,tag,MPI_COMM_WORLD,isend,info)
       call MPI_IRECV(data1(1,max_0+1),count,MPI_DOUBLE_COMPLEX,MpiRankA,tag,MPI_COMM_WORLD,irecv,info)
    else if(ip==mult_phi) then
       call MPI_ISEND(data2(1,min_0),count,MPI_DOUBLE_COMPLEX,MpiRankD,tag,MPI_COMM_WORLD,isend,info)
       call MPI_IRECV(data2(1,max_0+1),count,MPI_DOUBLE_COMPLEX,MpiRankA,tag,MPI_COMM_WORLD,irecv,info)
    end if
#elif DIM==3
    if(ip==1) then
       call MPI_ISEND(data1(1,1,min_0),count,MPI_DOUBLE_COMPLEX,MpiRankD,tag,MPI_COMM_WORLD,isend,info)
       call MPI_IRECV(data1(1,1,max_0+1),count,MPI_DOUBLE_COMPLEX,MpiRankA,tag,MPI_COMM_WORLD,irecv,info)
    else if(ip==mult_phi) then
       call MPI_ISEND(data2(1,1,min_0),count,MPI_DOUBLE_COMPLEX,MpiRankD,tag,MPI_COMM_WORLD,isend,info)
       call MPI_IRECV(data2(1,1,max_0+1),count,MPI_DOUBLE_COMPLEX,MpiRankA,tag,MPI_COMM_WORLD,irecv,info)
    end if
#endif
    call MPI_WAIT(isend,ssend,info)
    call MPI_WAIT(irecv,srecv,info)
    
#if DIM==2
    if(ip==1) then
       call MPI_ISEND(data1(1,max_0),count,MPI_DOUBLE_COMPLEX,MpiRankA,tag,MPI_COMM_WORLD,isend,info)
       call MPI_IRECV(data1(1,min_0-1),count,MPI_DOUBLE_COMPLEX,MpiRankD,tag,MPI_COMM_WORLD,irecv,info)
    else if(ip==mult_phi) then
       call MPI_ISEND(data2(1,max_0),count,MPI_DOUBLE_COMPLEX,MpiRankA,tag,MPI_COMM_WORLD,isend,info)
       call MPI_IRECV(data2(1,min_0-1),count,MPI_DOUBLE_COMPLEX,MpiRankD,tag,MPI_COMM_WORLD,irecv,info)
    end if
#elif DIM==3
    if(ip==1) then
       call MPI_ISEND(data1(1,1,max_0),count,MPI_DOUBLE_COMPLEX,MpiRankA,tag,MPI_COMM_WORLD,isend,info)
       call MPI_IRECV(data1(1,1,min_0-1),count,MPI_DOUBLE_COMPLEX,MpiRankD,tag,MPI_COMM_WORLD,irecv,info)
    else if(ip==mult_phi) then
       call MPI_ISEND(data2(1,1,max_0),count,MPI_DOUBLE_COMPLEX,MpiRankA,tag,MPI_COMM_WORLD,isend,info)
       call MPI_IRECV(data2(1,1,min_0-1),count,MPI_DOUBLE_COMPLEX,MpiRankD,tag,MPI_COMM_WORLD,irecv,info)
    end if
#endif
    call MPI_WAIT(isend,ssend,info)
    call MPI_WAIT(irecv,srecv,info)
    
  end subroutine update_margin
  
  subroutine init_evolve
  end subroutine init_evolve
  
  subroutine fin_evolve
  end subroutine fin_evolve
  
  subroutine write_chk(tau,idx,flag_defects)
    logical,intent(in) :: flag_defects
    double precision,intent(in) :: tau
    integer idx,info,iu
    character(len=charlen) str,chk_phi
    
    if(.not.flag_defects) then ! just modify param.chk
       if(MpiRank==0) then
          open(100,file=file_resume,action='write',status='old',iostat=info)
          write(100,'(I,E,L)',iostat=info) idx,tau,flag_defects
          close(100)
          write(*,'(A)') '  *param.chk is updated'
       end if
       return
    end if
    
    call io_chk(.false.,idx)
    
    if(MpiRank==0) then
       open(100,file=file_resume,action='write',status='old',iostat=info)
       write(100,'(I,E,L)',iostat=info) idx,tau,flag_defects
       close(100)
       write(*,'(A)') '  *param.chk is updated; deleting all previous chk files'
    end if
    
    do iu=0,idx-1
       write(str,*) iu
       str=adjustl(str)
       if(tmpdir/='') then
          chk_phi=trim(tmpdir)//trim(tmprn)//'_phi_'//trim(str)//'.chk'
       else
          chk_phi=trim(rootname)//'_phi_'//trim(str)//'.chk'
       end if
       call MPI_FILE_DELETE(chk_phi,MPI_INFO_NULL,info)
    end do
    
  end subroutine write_chk
  
  subroutine read_chk(tau,idx,flag_defects)
    integer idx
    logical,intent(out) :: flag_defects
    double precision,intent(out) :: tau
    integer info
    
    open(100,file=file_resume,action='read',status='old',iostat=info)
    if(info==0) read(100,*,iostat=info) idx,tau,flag_defects
    if(info/=0) call MpiStop('error: cannot read param.chk')
    close(100)
    
    call io_chk(.true.,idx)
    call update_margin(1)
    
  end subroutine read_chk
  
  subroutine io_chk(flag,idx)
    logical,intent(in) :: flag
    integer,intent(in) :: idx
    integer info,amode,count,count0,iu,i,inc0,inc,iter,niter
    integer(8) iz
    
    integer(kind=MPI_OFFSET_KIND) disp
    character(len=charlen) str,chk_phi
    integer sreqs(MPI_STATUS_SIZE)
    
    write(str,*) idx
    str=adjustl(str)
    if(tmpdir/='') then
       chk_phi=trim(tmpdir)//trim(tmprn)//'_phi_'//trim(str)//'.chk'
    else
       chk_phi=trim(rootname)//'_phi_'//trim(str)//'.chk'
    end if
    
    if(flag) then ! read
       amode=MPI_MODE_RDONLY
       str='read'
    else ! write
       amode=MPI_MODE_WRONLY+MPI_MODE_CREATE
       str='write'
    end if
    
    if(MpiRank==0) then
       write(*,'(A)') '  calling io_chk to '//trim(str)//' phi and phidot ...'
       write(*,'(A)') '   file: '//trim(chk_phi)
    end if
    
    call MPI_FILE_OPEN(MPI_COMM_WORLD,trim(chk_phi),amode,MPI_INFO_NULL,iu,info)
    if(info/=MPI_SUCCESS) call MpiStop('error: MPI_FILE_OPEN failed')
    
    count0=ngrid_slice
    call MPI_FILE_SET_VIEW(iu,0_MPI_OFFSET_KIND,MPI_DOUBLE_COMPLEX,MPI_DOUBLE_COMPLEX,'native',MPI_INFO_NULL,info)
    if(info/=0) call MpiStop('error: MPI_FILE_SET_VIEW failed')
    
    inc0=local_0/ceiling((dble(count0)*local_0)/count_dc_max) ! number of z-slices input/output at one time
    if(.not.inc0>0) call MpiStop('error: need further decomposition along y?')
    
    niter=(local_0-1)/inc0+1

    do i=1,mult_phi
       do iz=min_0,max_0,inc0
          inc=min(inc0,local_0-iz+min_0)
          count=count0*inc
          disp=int(count0,kind=MPI_OFFSET_KIND)
          disp=disp*(ngrid*(i-1)+local_0_offset+iz-min_0)
          write(*,'(A,I,A,I,A,I)') '   MPIRANK: ',MpiRank,' disp: ',disp,' count: ',count
          iter=(iz-min_0)/inc0+1+niter*(i-1)
          
          if(flag) then
#if DIM==2
             if(i==1) then
                call MPI_FILE_READ_AT(iu,disp,data1(1,iz),count,MPI_DOUBLE_COMPLEX,sreqs,info)
             else if(i==mult_phi) then
                call MPI_FILE_READ_AT(iu,disp,data2(1,iz),count,MPI_DOUBLE_COMPLEX,sreqs,info)
             end if
#elif DIM==3
             if(i==1) then
                call MPI_FILE_READ_AT(iu,disp,data1(1,1,iz),count,MPI_DOUBLE_COMPLEX,sreqs,info)
             else if(i==mult_phi) then
                call MPI_FILE_READ_AT(iu,disp,data2(1,1,iz),count,MPI_DOUBLE_COMPLEX,sreqs,info)
             end if
#endif
          else
#if DIM==2
             if(i==1) then
                call MPI_FILE_WRITE_AT(iu,disp,data1(1,iz),count,MPI_DOUBLE_COMPLEX,sreqs,info)
             else if(i==mult_phi) then
                call MPI_FILE_WRITE_AT(iu,disp,data2(1,iz),count,MPI_DOUBLE_COMPLEX,sreqs,info)
             end if
#elif DIM==3
             if(i==1) then
                call MPI_FILE_WRITE_AT(iu,disp,data1(1,1,iz),count,MPI_DOUBLE_COMPLEX,sreqs,info)
             else if(i==mult_phi) then
                call MPI_FILE_WRITE_AT(iu,disp,data2(1,1,iz),count,MPI_DOUBLE_COMPLEX,sreqs,info)
             end if
#endif
          end if
          if(info/=MPI_SUCCESS) call MpiStop('error: MPI_FILE_READ/WRITE_AT failed',info)
          write(*,'(A,I,A,I,A,I,A,I)') '    MpiRank ',MpiRank,', i=',i,', iterations',iter,'/',niter*mult_phi
       end do
    end do
    
    call MPI_FILE_CLOSE(iu,info)
    if(info/=MPI_SUCCESS) call MpiStop('error: MPI_FILE_CLOSE failed')
    
    if(MpiRank==0) write(*,'(A)') '  ... done'
    
  end subroutine io_chk
  
end module evolve
