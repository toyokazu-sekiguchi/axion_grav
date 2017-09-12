#ifndef MPI
#error 'MPI is required'
#endif

#define THREADED_MPI MPI_THREAD_MULTIPLE
#define THREADED_FFT .true.
  
module settings
#ifdef _OPENMP
  use omp_lib
#endif
  use utils
  use potential,only : pi,twopi
  use,intrinsic :: ISO_C_BINDING
  implicit none
  
  integer,parameter :: feedback=2
  integer,parameter :: requested_mpi_thread=THREADED_MPI
  integer provided_mpi_thread
  logical,parameter :: fft_threaded=THREADED_FFT
  
#if DIM==2
  integer,parameter :: ndim=2
#elif DIM==3
  integer,parameter :: ndim=3
#endif

  
  double precision,parameter :: byte_dc=16d0, &
       hard_limit_cdata=2*byte_dc, & ! memory limit of cdata (GB) ie 2^31 double complex
       hard_limit_io_buffer=2d0 ! limit of I/O buffering should not be larger than 2GB
  integer,parameter :: mult_gauss=2,mult_comp=2,mult_phi=2, &
       count_dc_max=int(hard_limit_io_buffer/byte_dc*1024*1024*1024) ! maximum number of double complex that can be I/O at one time, ie 2^27
  
  character(len=charlen) rootname,rootdir,tmprn,tmpdir,file_resume
  integer n_chk,n_id,num_bin,num_threads,n_mod_eps
  logical spec_est,dry_run,out_strings,out_walls,id_strings,id_walls,group_strings,group_walls,grav_spec
  double precision limit_memory ! available memory in GB; read at runtime
  
  !=========================
  ! for MPI FFT
  !=========================
  integer(8) ngrid,ngrid_tot,ngo2,ngrid_slice,local_0,local_0_offset,min_0,max_0
  integer(8),allocatable :: arr_local_0(:),arr_local_0_offset(:)
  integer :: degrade_spec=0
#if DIM==2
  complex(C_DOUBLE_COMPLEX),pointer :: data1(:,:),data2(:,:)
#elif DIM==3
  complex(C_DOUBLE_COMPLEX),pointer :: data1(:,:,:),data2(:,:,:)
#endif
  
  !=========================
  ! simulation parameters
  !=========================
  double precision sigma,g_rel,sigma_cr,zeta,width_cr,thickness_cr,tau_cr,a_now
  double precision tau_ini,tau_fin,boxsize,dspa,dwav,dv_spa,dv_wav,eps_dtau
  real tsim_begin,tsim_now,tsim_max ! walltime 
  
  double precision,allocatable :: tau_id(:)
  
!!$ initial loop
  logical initial_loop
  double precision loop_radius, & ! loop size in units of half of box size
       loop_velocity, & ! loop velocity at the point of largest torque; 
       loop_velocity_in ! inward loop velocity
!!$ initial_loop
  
contains
  
  subroutine initialization(ParamFile,str)
    use inifile
    use potential
    use subroutines
    character(*),intent(in) :: ParamFile,str
    character(len=charlen) ss
    integer i,info
    double precision b
    logical bad
    
#ifdef _OPENMP
    call Mpi_Init_Thread(requested_mpi_thread,provided_mpi_thread,info)
#else
    call Mpi_Init(info)
#endif
    if(info/=MPI_SUCCESS) stop 'error: cannot initialize MPI'
    call MPI_COMM_RANK(MPI_COMM_WORLD,MpiRank,info)
    if(info/=MPI_SUCCESS) call MpiStop('error: cannot get MpiRank')
    call MPI_COMM_SIZE(MPI_COMM_WORLD,MpiSize,info)
    if(info/=MPI_SUCCESS) call MpiStop('error: cannot get MpiSize')
#ifdef _OPENMP
    if(MpiRank==0) write(*,'(A,I,A,I)') ' requested_mpi_thread:',requested_mpi_thread, &
         ',  provided_mpi_thread:',provided_mpi_thread
#endif
    if(MpiRank==0) then
       write(*,*)
       write(*,*) '#initialization'
    end if
    
    call memory_check_mkl(0)
    
    ! output several sentinel values used for computation
    if(C_INTPTR_T/=8) call MpiStop('error: C_INTPTR_T/=8: this is not supported')
    if(MpiRank==0) then
       write(*,'(A)') ' *check some sentinel parameters'
       write(*,'(A,I)') '   C_INTPTR_T = ',C_INTPTR_T
       write(*,'(A,I)') '   C_DOUBLE_COMPLEX = ',C_DOUBLE_COMPLEX
       write(*,'(A,I)') '   MPI_OFFSET_KIND = ',MPI_OFFSET_KIND
    end if
    
    if(MpiRank==0) write(*,'(A,I)') ' number of tasks = ',MpiSize
    
    call Ini_Open(ParamFile,1,bad,.false.)
    Ini_fail_on_not_found=.false.
    
#ifdef _OPENMP
    i=OMP_GET_MAX_THREADS()
    num_threads=Ini_Read_Int('num_threads',i)
    if(.not.(num_threads>0.and.num_threads<=i)) num_threads=i
    call OMP_SET_NUM_THREADS(int(num_threads,kind=4))
    if(MpiRank==0) write(*,'(A,I3)') ' openmp parallelized with num_threads = ',num_threads
#else
    num_threads=1
    if(MpiRank==0) write(*,'(A)') ' no openmp parallelization'
#endif
    
    dry_run=Ini_Read_Logical('dry_run')
    rootdir=Ini_Read_String('rootdir')
    if(rootdir(len_trim(rootdir):len_trim(rootdir))/='/') rootdir=trim(rootdir)//'/'
    tmprn=Ini_Read_String('rootname')
    rootname=trim(rootdir)//trim(tmprn)
    if(str/='') then
       tmprn=trim(tmprn)//'_'//trim(str)
       rootname=trim(rootname)//'_'//trim(str)
    end if
    if(MpiRank==0) write(*,'(A)') ' rootname = '//trim(rootname)
    
    tmpdir=Ini_Read_String('tmpdir')
    if(tmpdir/='') then
       if(tmpdir(len_trim(tmpdir):len_trim(tmpdir))/='/') tmpdir=trim(tmpdir)//'/'
    else
       tmpdir=rootdir
    end if
    if(MpiRank==0.and.tmpdir/='') write(*,'(A)') ' temporary output is in '//trim(tmpdir)
    
    ss=Ini_Read_String('ss')
    call init_rand_gen(num_threads,ss)
    
    ngrid=Ini_Read_Int('ngrid')
    ngo2=ngrid/2
    if(ngo2*2/=ngrid) call MpiStop('error: ngrid should be even')
    ngrid_slice=ngrid**(ndim-1)
    if(ngrid_slice/ngrid/=ngrid**(ndim-2)) call MpiStop('error: womething wrong with ngrid_slice')
    ngrid_tot=ngrid*ngrid_slice
    if(ngrid_tot/ngrid/=ngrid_slice) call MpiStop('error: womething wrong with ngrid_tot')
    
    n_chk=Ini_Read_Int('n_chk')
    if(.not.n_chk>=1) call MpiStop('error: n_chk should be 1 or larger')
    num_bin=Ini_Read_Int('num_bin',50)
!!$ initial_loop
    initial_loop=Ini_Read_Logical('initial_loop',.false.) ! default is thermal
    if(initial_loop) then
       loop_radius=Ini_Read_Double('loop_radius')
       loop_velocity=Ini_Read_Double('loop_velocity')
       loop_velocity_in=Ini_Read_Double('loop_velocity_in')
    end if
!!$ initial_loop
    id_strings=Ini_Read_Logical('id_strings')
    if(id_strings) then
#if DIM==2
       group_strings=.false.
#elif DIM==3
       group_strings=Ini_Read_Logical('group_strings')
#endif
       out_strings=Ini_Read_Logical('out_strings')
    else
       out_strings=.false.
    end if
    id_walls=Ini_Read_Logical('id_walls')
    if(id_walls) then
#if DIM==2
       group_walls=.false.
#elif DIM==3
       group_walls=Ini_Read_Logical('group_walls')
#endif
       out_walls=Ini_Read_Logical('out_walls')
    else
       out_walls=.false.
    end if

    spec_est=Ini_Read_Logical('spec_est')
    grav_spec=Ini_Read_Logical('grav_spec')

    n_id=0
    if(id_strings.or.id_walls.or.spec_est.or.grav_spec) then
       n_id=Ini_Read_Int('n_id')
       if(spec_est) n_mod_eps=Ini_Read_Int('n_mod_eps')
    end if
    
    if(spec_est.or.grav_spec) then
       allocate(tau_id(0:n_chk),stat=info)
       if(info/=0) call MpiStop('error: cannot allocate tau_id')
       degrade_spec=Ini_Read_Int('degrade_spec')
    end if
    
    if(MpiRank==0) then
       write(*,'(A,I)') ' ndim = ',ndim
       write(*,'(A,I)') ' ngrid = ',ngrid
       write(*,'(A,I)') ' ngrid_tot = ',ngrid_tot
       write(*,'(A,I)') ' ngrid_slice = ',ngrid_slice
       write(*,'(A,I)') ' n_chk = ',n_chk
       
       if(initial_loop) then
          write(*,'(A)') ' initial configuration: loop'
          write(*,'(A,F)') '  loop radius in units of box size: ',loop_radius
          write(*,'(A,F)') '  loop rotation velocity: ',loop_velocity
          write(*,'(A,F)') '  loop inward velocity: ',loop_velocity_in
       else
          write(*,'(A)') ' initial configuration: random (thermal)'
       end if
       
       write(*,'(A,L)') ' identify strings: ', id_strings
       if(id_strings) then
          write(*,'(A,L)') '  grouping strings: ', group_strings
          write(*,'(A,L)') '  output string positions: ', out_strings
       end if
       write(*,'(A,L)') ' identify walls: ', id_walls
       if(id_walls) then
          write(*,'(A,L)') '  grouping walls: ', group_walls
          write(*,'(A,L)') '  output wall positions: ', out_walls
       end if
       write(*,'(A,L)') ' estimation of spectrum = ',spec_est
       if(id_strings.or.id_walls.or.spec_est) then
          write(*,'(A,I)') ' n_id = ',n_id
          if(spec_est) then
             write(*,'(A,I)') ' n_mod_eps = ',n_mod_eps
             write(*,'(A,I)') ' degrade_spec = ',degrade_spec
          end if
       end if
    end if
    
    b=byte_dc*mult_phi
    if(spec_est) then
       b=b+2*byte_dc/2**(degrade_spec*ndim) ! 2 counts for data_spec1 and data_spec2
       b=b+max(byte_dc/2**(degrade_spec*ndim),1d0) ! data_phi_degrade or exist_defect
    else if(id_strings.or.id_walls) then
       b=b+1d0 ! exist_defect
    end if
    b=b*ngrid_tot/1024/1024/1024
    if(MpiRank==0) write(*,'(A,F)') &
         ' approx required total memory (GB) (without margin) = ',b
    b=b/MpiSize
    if(MpiRank==0) write(*,'(A,F)') &
         ' approx required memory per task (GB) (without margin) = ',b
    if(b>=mult_phi*hard_limit_cdata.and..not.dry_run) call MpiStop('error: to large cdata; increase the number of tasks')
    b=b+2*byte_dc*mult_phi*ngrid_slice/1024/1024/1024
    if(MpiRank==0) write(*,'(A,F)') &
         ' approx required memory per task (GB) (with margin but wihout sys) = ',b
    if(MpiRank==0) write(*,'(A,F)') &
         ' available memory per task (GB) = ',limit_memory
    if(b>=limit_memory.and..not.dry_run) call MpiStop('error: not enough memory')
    
    if(tmpdir/='') then
       file_resume=trim(tmpdir)//trim(tmprn)//'_param.chk'
    else
       file_resume=trim(rootname)//'_param.chk'
    end if
    if(MpiRank==0) write(*,'(A)') ' resume file:'//trim(file_resume)
    
    sigma=Ini_Read_Double('sigma')
    g_rel=Ini_Read_Double('g_rel')
    kappa=Ini_Read_Double('kappa')
!!$
!!$    sigma_cr = 1/pi/sqrt(g_rel/90.d0)/sigma/thrm_sigma**2 ! sigma*tau_cr
    a_now=sigma*thrm_sigma/Tcmb_in_Mstar*(g_rel/g_rel_nowS)**(1d0/3) ! present scale factor with a normalization a_cr=1
    tau_cr = 1/pi/sqrt(g_rel/90.d0)/sigma**2/thrm_sigma**2 ! tau_cr in units of (reduced planck mass)^-1
    sigma_cr = tau_cr*sigma ! sigma*tau_cr
!!$
    zeta=sqrt(22.5d0/g_rel)/sigma/pi
    width_cr = 1/sqrt(2d0)/sigma_cr
    thickness_cr = 1/sqrt(2d0)/(sqrt(maxion2(1.d-3*kappa))*sigma_cr)
    if(MpiRank==0) then
       write(*,'(A,F)') ' sigma/M_* = ',sigma
       write(*,'(A,F)') ' g_rel = ',g_rel
       write(*,'(A,F)') ' kappa = ',kappa
       write(*,'(A,F)') ' zeta_s = ',zeta
       write(*,'(A,F)') ' zeta_w = ',zeta*c0_axion
       write(*,'(A,3F)') ' maxion parameters n, cT, c0:',n_axion_IILM,cT_axion_IILM,c0_axion
    end if
    
    tau_ini=Ini_Read_Double('tau_ini')
    tau_fin=Ini_Read_Double('tau_fin')
    eps_dtau=Ini_Read_Double('eps_dtau')
    boxsize=Ini_Read_Double('boxsize')
    dspa=boxsize*tau_fin/ngrid
    dwav=twopi/dspa/ngrid
    dv_spa=dspa**ndim
    dv_wav=(dwav/twopi)**ndim
    
    if(MpiRank==0) then
       write(*,'(A,F)') ' tau_ini/tau_cr = ',tau_ini
       if(.not.initial_loop) call check_tau_ini
       write(*,'(A,F)') ' tau_fin/tau_cr = ',tau_fin
       call check_tau_fin
       write(*,'(A,F)') ' eps_dtau = ',eps_dtau
       write(*,'(A,F)') ' box size (in units of tau_fin) = ',boxsize
       write(*,'(A,F)') ' delta_lattice (in units of tau_cr) = ',dspa
    end if

#ifdef MPI
    call MPI_BARRIER(MPI_COMM_WORLD,info)
#endif

    
    call Ini_Close
    
  end subroutine initialization
  
  subroutine check_tau_fin
    double precision width,thickness
    
    width=width_cr/tau_fin ! comoving size of strings width at tau_fin
    write(*,'(A,E)') '  dspa/width at tau_fin (should be < 1??)= ',dspa/width
    if(dspa>width) write(*,'(A)') '   WARNING??: resolution is not enough; tau_fin should be shorter'
    
    thickness=thickness_cr/tau_fin
    write(*,'(A,E)') '  dspa/thickness at tau_fin (should be < 1)= ',dspa/thickness
    
  end subroutine check_tau_fin
  
  subroutine check_tau_ini
    use potential
    double precision thrm,wav,mass
    
    thrm=thrm_sigma/tau_ini ! initial temperature (in sigma)
    wav=dwav/tau_ini/sigma_cr ! minimum nonzero physical momentum (in sigma)
    mass=sqrt(msquared(thrm))
    write(*,'(A,E)') '  dwav/R_ini/mass_ini (should be < 1/a few) = ',wav/mass
    if(wav>mass/3) then
       write(*,'(A)') '   WARNING: boxsize not large enough or initial mass too small;'
       write(*,'(A)') '   WARNING: initial configuration may not be accurate enough'
    end if
    write(*,'(A,E)') '  mass_ini/T_ini (should be < O(0.1)) = ',mass/thrm
    if(mass>thrm/10) then
       write(*,'(A)') '   WARNING: too large initial mass;'
       write(*,'(A)') '   WARNING: tau_ini may be not enough close to tau_cr'
    end if
    wav=wav*ngrid 
    wav=sqrt(wav**2+mass**2) ! ~ maximum physical momentum
    write(*,'(A,E)') '  wav_max/R_ini/T_ini (should be > O(10)) = ',wav/thrm
    if(wav<10*thrm) then
       write(*,'(A)') '   WARNING: not enough resolution;'
    end if
    
  end subroutine check_tau_ini
  
  subroutine finalization
    use subroutines
    integer info
    character(len=charlen) file_chk
    
    if(MpiRank==0) then
       write(*,*)
       write(*,'(A)') ' #finalising simulation...'
       
       if(tmpdir=='') then
          file_chk=trim(rootname)//'_p*.chk'
          call system('rm -rf '//trim(file_chk))
       end if
       
    end if
    
    if(spec_est.or.grav_spec) deallocate(tau_id,stat=info)
    
    call fin_rand_gen
    call memory_check_mkl()
    call memory_check_mkl(-1)
    
    call MPI_BARRIER(MPI_COMM_WORLD,info)
    call notice_walltime
    call MPI_FINALIZE(info)
    
    if(MpiRank==0) write(*,'(A)') ' ... done'
    
  end subroutine finalization
  
  subroutine notice_walltime
#include <mkl.fi>
    tsim_now=second()/3600.d0
    if(MpiRank==0) write(*,'(A,2E)') '   wall time [hour] from start, ellapse time [hour]:',tsim_now-tsim_begin,tsim_max
  end subroutine notice_walltime
  
  subroutine memory_check_mkl(i)
#include <mkl.fi>
    integer,intent(in),optional :: i
    integer(8) byte
    integer,save :: iter=0
    integer info
    
    call MPI_BARRIER(MPI_COMM_WORLD,info)
    
    if(present(i)) then
       if(i==0) then ! start
          byte=mkl_peak_mem_usage(MKL_PEAK_MEM_ENABLE)
          if(byte<0) call MpiStop('error: MKL_PEAK_MEM_ENABEL failed')
       else if(i==-1) then
          byte=mkl_peak_mem_usage(MKL_PEAK_MEM_DISABLE)
          if(byte<0) call MpiStop('error: MKL_PEAK_MEM_DISABLEL failed')
       else if(i==1) then
          byte=mkl_peak_mem_usage(MKL_PEAK_MEM_RESET)
          if(byte<0) call MpiStop('error: MKL_PEAK_MEM_RESET failed')
       end if
    else
       iter=iter+1
       byte=mkl_peak_mem_usage(MKL_PEAK_MEM)
       if(byte<0) call MpiStop('error: MKL_PEAK_MEM failed')
       print '(A,2I,E)',' - peak memory usage by MKL [GByte]:',MpiRank,iter,dble(byte)/1024/1024/1024
    end if
    
    call MPI_BARRIER(MPI_COMM_WORLD,info)
    if(MpiRank==0) write(*,*)
    
  end subroutine memory_check_mkl
  
end module settings
