program lattice_axion
  use utils
  use settings
  use initial
  use evolve
  use spectrum
  use grav_waves
  implicit none
#include <fftw3-mpi.f03>
#include <mkl.fi>
  integer,parameter :: iarg_max=7
  character(len=charlen) paramfile,str,str0(iarg_max)
  integer iarg,i,info
  logical chk
  
  iarg=iargc()
  if(iarg>iarg_max) stop 'error: too many arguments'
  if(.not.iarg>=1) stop 'error: parameter file is absent'
  do i=1,iarg_max
     call getarg(i,str0(i))
     str0(i)=adjustl(str0(i))
  end do
  paramfile=trim(str0(1))
  str=''
  tsim_max=0d0
  limit_memory=0d0
  if(iarg>=2) then
     do i=2,iarg,2
        if(trim(str0(i))=='-n') str=str0(i+1)
        if(trim(str0(i))=='-t') read(str0(i+1),*) tsim_max
        if(trim(str0(i))=='-b') read(str0(i+1),*) limit_memory
     end do
  end if
  tsim_begin=second()/3600.d0
  
  call initialization(paramfile,str)
  
  if(dry_run) then
     
     if(MpiRank==0) write(*,'(A)') ' *** dry run just for checking parameters ***'
     
  else ! simulation
     
     if(provided_mpi_thread>MPI_THREAD_FUNNELED) then
        info=FFTW_INIT_THREADS()
        if(info==0) call MpiStop(' error something wroing in FFTW_INIT_THREADS')
        if(MpiRank==0) write(*,*) ' threaded FFT is enabled'
     end if
     call FFTW_MPI_INIT()
     
     call set_initial_condition(chk)

     if(spec_est) call init_spectrum
     if(grav_spec) call init_grav
     
     call evolution(chk)
     
     if(spec_est) call fin_spectrum
     if(grav_spec) call fin_grav
     
     call FFTW_MPI_CLEANUP()
     
  end if
  
  call finalization
  
end program lattice_axion
