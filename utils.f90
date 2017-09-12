module utils
  implicit none
#ifdef MPI
  include "mpif.h"
#endif

#ifdef SINGLE
  integer,parameter :: lpc=kind(1.0)
#else
  integer,parameter :: lpc=kind(1d0)
#endif
#ifdef MPI
#ifdef SINGLE
  integer,parameter :: Mpi_lpc=MPI_REAL
#else
  integer,parameter :: Mpi_lpc=MPI_DOUBLE_PRECISION
#endif
#endif

  integer,parameter :: charlen=1024
  integer MpiRank,MpiSize

contains

!!$  function GetMpiRank()
!!$    integer GetMpiRank
!!$#ifdef MPI
!!$    integer info
!!$    call MPI_COMM_RANK(MPI_COMM_WORLD,GetMpiRank,info)
!!$    if(info/=MPI_SUCCESS) call MpiStop('MPI failed')
!!$#else
!!$    GetMpiRank=0
!!$#endif
!!$
!!$  end function GetMpiRank
    
!!$  function IsMasterMpi()
!!$    logical IsMasterMpi
!!$
!!$    IsMasterMpi=GetMpiRank()==0
!!$
!!$  end function IsMasterMpi

  subroutine MpiStop(str,i)
    character(*),intent(in),optional :: str
!!$
    integer,intent(in),optional :: i
!!$
#ifdef MPI
    integer info,MpiRank
!!$
    character(len=charlen) str0
    integer n
!!$
#endif
    
    if(present(str)) write(*,*) trim(str)

#ifdef MPI
    call MPI_COMM_RANK(MPI_COMM_WORLD,MpiRank,info)
    write(*,*) 'MpiStop: ', MpiRank
!!$
    if(present(i)) then
       call MPI_ERROR_STRING(i,str0,n,info)
       write(*,'(A)') trim(str0)
    end if
!!$
    call MPI_ABORT(MPI_COMM_WORLD,i)
#endif
!!$    i=1
    stop

  end subroutine MpiStop

!!$  subroutine MpiStat(MpiID,MpiSize)
!!$    implicit none
!!$    integer MpiID,MpiSize
!!$#ifdef MPI
!!$    integer info
!!$    call MPI_COMM_RANK(MPI_COMM_WORLD,MpiID,info)
!!$    if(info/=MPI_SUCCESS) call MpiStop('MpiStat: MPI rank')
!!$    call MPI_COMM_SIZE(MPI_COMM_WORLD,MpiSize,info)
!!$#else
!!$    MpiID=0
!!$    MpiSize=1
!!$#endif
!!$  end subroutine MpiStat
  
!!$  subroutine MpiQuietWait
!!$    ! Set MPI thread to sleep, e.g. so can run openmp on cpu instead
!!$#ifdef MPI
!!$    integer flag,info,status(MPI_STATUS_SIZE)
!!$    integer i,MpiID,MpiSize
!!$    
!!$    call MPiStat(MpiID,MpiSize)
!!$    if(MpiID/=0) then
!!$       do 
!!$          call MPI_IPROBE(0,0,MPI_COMM_WORLD,flag,MPI_STATUS_IGNORE,info)
!!$          if(flag/=0) then
!!$             call MPI_RECV(i,1,MPI_INTEGER,0,0,MPI_COMM_WORLD,status,info)
!!$             exit
!!$          end if
!!$          call sleep(1)
!!$       end do
!!$    end if
!!$#endif
!!$  end subroutine MpiQuietWait
!!$
!!$  subroutine MpiWakeQuietWait
!!$#ifdef MPI
!!$    integer j,MpiID,MpiSize,info,r
!!$    
!!$    call MpiStat(MpiID,MpiSize)
!!$    if (MpiID==0) then
!!$       do j=1, MpiSize-1
!!$          call MPI_ISSEND(MpiID,1,MPI_INTEGER,j,0,MPI_COMM_WORLD,r,info)
!!$       end do
!!$    end if
!!$#endif
!!$   end subroutine MpiWakeQuietWait

end module utils
