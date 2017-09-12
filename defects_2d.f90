module defects_2d
  use settings
  implicit none

  double precision,parameter :: string_phdiff=pi, &
       removal_phdiff=pi/4
!!$  double precision,parameter :: def_th_dist=2.d0 ! this might need to be tuned...

  type mask_list
     integer ixyz(ndim)
     type(mask_list),pointer :: next
  end type mask_list
  
  type string_list
     integer ixyz(ndim),sgn
     double precision pos(ndim),vel(ndim),gof
     type(string_list),pointer :: next
  end type string_list
  
  type wall_list
     integer ixyz(ndim),sgn
     double precision pos(ndim),vel,gof,wa
     type(wall_list),pointer :: next
  end type wall_list
  
contains

  subroutine id_defects_2d(idx,wid_dw,ns,wa)
    double precision,intent(in) :: wid_dw ! wall thickness in units of lattice spacing
    integer,intent(in) :: idx
    integer,intent(out) :: ns ! number of strings
    double precision,intent(out) :: wa ! wall area 
    
    integer iy1,iz1,iy2,iz2,info,MpiRank,MpiSize,i
    double precision ph(0:1,0:1),phdiff,ph0
    character(len=charlen) file_string,file_wall,file_mask,str
    type(string_list),pointer :: first_string,tmp_string
    type(wall_list),pointer :: first_wall,tmp_wall
    type(mask_list),pointer :: first_mask,tmp_mask
    logical exist_string,exist_wall
    
    call MpiStat(MpiRank,MpiSize)
    
    if(MpiRank==0) write(*,'(A)') ' *identifying strings and walls ...'
    
    if(id_strings) nullify(first_string)
    if(id_walls) nullify(first_wall)
    if(spec_est) nullify(first_mask)
    
    ns=0
    do iz1=min_1,max_1
       iz2=iz1+1
       do iy1=1,ngrid
          iy2=iy1+1
          if(iy1==ngrid) iy2=1
          
          ph(0,0)=aimag(log(data_phi(iy1,iz1)))
          ph(1,0)=aimag(log(data_phi(iy2,iz1)))
          ph(0,1)=aimag(log(data_phi(iy1,iz2)))
          ph(1,1)=aimag(log(data_phi(iy2,iz2)))
          
          call id_with_phase_2d(ph,exist_string,exist_wall,phdiff,ph0)
          
          if(id_walls.and.exist_wall) then ! wall exists
             allocate(tmp_wall)
             tmp_wall%ixyz(1:ndim)=(/iy1,iz1/)
             call position_wall(wid_dw,ph,tmp_wall%pos,tmp_wall%vel,tmp_wall%sgn,tmp_wall%gof,tmp_wall%wa)
             tmp_wall%next=>first_wall
             first_wall=>tmp_wall
          end if
          
          if(id_strings.and.exist_string) then ! string exists
             call rotate_phase_2d(ph,ph0)
             allocate(tmp_string)
             tmp_string%ixyz(1:ndim)=(/iy1,iz1/)
             call position_string(ph,tmp_string%pos,tmp_string%vel,tmp_string%sgn,tmp_string%gof)
             tmp_string%next=>first_string
             first_string=>tmp_string
             ns=ns+1 
          end if
          
          if(spec_est.and.phdiff>removal_phdiff) then ! mask
             allocate(tmp_mask)
             tmp_mask%ixyz(1:ndim)=(/iy1,iz1/)
             tmp_mask%next=>first_mask
             first_mask=>tmp_mask
          end if
          
       end do
    end do
    
    if(id_walls) call add_wall_area(first_wall,wa)
    
    write(str,*) idx
    str=adjustl(str)
    if(out_strings) file_string=trim(rootname)//'_string_'//trim(str)//'.dat'
    if(out_walls) file_wall=trim(rootname)//'_wall_'//trim(str)//'.dat'
    if(spec_est) file_mask=trim(rootname)//'_mask_'//trim(adjustl(str))//'.dat'
    
    if(out_strings.or.out_walls) then
       if(MpiRank==0) write(*,'(A)') '  sequencially outputting string, wall (and mask)'
    end if
    
    do i=0,MpiSize-1 ! sequencially output string and mask; this may be inefficient ...
       if(i==MpiRank) then
          
          if(id_strings) then
             if(out_strings) then
                if(MpiRank==0) then
                   open(100,file=trim(file_string),action='write',status='replace',iostat=info)
                   if(info/=0) call MpiStop(' error: cannot create string.dat')
                else
                   open(100,file=trim(file_string),action='write',status='old',position='append',iostat=info)
                   if(info/=0) call MpiStop(' error: cannot open string.dat')
                end if
             end if
             tmp_string=>first_string
             do while(associated(tmp_string))
                if(out_strings) &
                     write(100,'(2I12,4E20.8,I3,E20.8)',iostat=info) tmp_string%ixyz(1:ndim), &
                     tmp_string%ixyz(1:ndim)+tmp_string%pos(1:ndim),tmp_string%vel(1:ndim),tmp_string%sgn,tmp_string%gof
                tmp_string=>tmp_string%next
                deallocate(first_string)
                first_string=>tmp_string
             end do
             if(out_strings) close(100)
          end if
          
          if(id_walls) then
             if(out_walls) then
                if(MpiRank==0) then
                   open(100,file=trim(file_wall),action='write',status='replace',iostat=info)
                   if(info/=0) call MpiStop(' error: cannot create wall.dat')
                else
                   open(100,file=trim(file_wall),action='write',status='old',position='append',iostat=info)
                   if(info/=0) call MpiStop(' error: cannot open wall.dat')
                end if
             end if
             tmp_wall=>first_wall
             do while(associated(tmp_wall))
                if(out_walls) &
                     write(100,'(2I12,3E20.8,I3,2E20.8)',iostat=info) tmp_wall%ixyz(1:ndim), &
                     tmp_wall%ixyz(1:ndim)+tmp_wall%pos(1:ndim),tmp_wall%vel,tmp_wall%sgn,tmp_wall%gof,tmp_wall%wa
                tmp_wall=>tmp_wall%next
                deallocate(first_wall)
                first_wall=>tmp_wall
             end do
             if(out_walls) close(100)
          end if
          
          if(spec_est) then
             if(MpiRank==0) then
                open(100,file=trim(file_mask),action='write',status='replace',iostat=info)
                if(info/=0) call MpiStop(' error: cannot create mask.dat')
             else
                open(100,file=trim(file_mask),action='write',status='old',position='append',iostat=info)
                if(info/=0) call MpiStop(' error: cannot open mask.dat')
             end if
             tmp_mask=>first_mask
             do while(associated(tmp_mask))
                write(100,'(2I)',iostat=info) tmp_mask%ixyz(1:ndim)
                tmp_mask=>tmp_mask%next
                deallocate(first_mask)
                first_mask=>tmp_mask
             end do
             close(100)
          end if
          
       end if
       
#ifdef MPI
       call MPI_BARRIER(MPI_COMM_WORLD,info)
#endif
       
    end do
    
    if(id_strings) deallocate(tmp_string)
    if(id_walls) deallocate(tmp_wall)
    if(spec_est) deallocate(tmp_mask)
    
    if(MpiRank==0) write(*,'(A)') ' *...done'
    
  end subroutine id_defects_2d
  
  subroutine id_with_phase_2d(ph,exist_string,exist_wall,phase_diff,ph0)
    use potential,only : twopi
    double precision,intent(in) :: ph(0:1,0:1)
    logical,intent(out) :: exist_string,exist_wall
    double precision,intent(out) :: phase_diff,ph0
    
    integer,parameter :: n=2**ndim
    integer i,j,ix,iy,iz
    double precision t,pht(n),tmax
    
    pht(1:2)=ph(0:1,0)
    pht(3:4)=ph(0:1,1)
    
    ! sort
    do i=n-1,1,-1
       do j=1,i
          if(pht(j)>pht(j+1)) then
             t=pht(j)
             pht(j)=pht(j+1)
             pht(j+1)=t
          end if
       end do
    end do
    
    exist_wall=.false.
    exist_string=.false.
    
    t=pht(n)-twopi
    tmax=0.d0
    do i=1,n
       t=pht(i)-t
       if(t>tmax) then
          tmax=t
          ph0=pht(i)-t/2 ! center of the largest phase separation
          if(i>1) exist_wall=.true.! largest phase separation does not contain theta=+-pi ==> DW exist
       end if
       t=pht(i)
    end do
    if(tmax<=string_phdiff) exist_string=.true.
    
    phase_diff=twopi-tmax
    
  end subroutine id_with_phase_2d

  subroutine rotate_phase_2d(ph,ph0)
    use potential,only : pi,twopi
    double precision,intent(inout) :: ph(0:1,0:1)
    double precision,intent(in) :: ph0
    integer iy,iz
    double precision t
    
    ! in the end ph is rotated so as to avoid the discontinuity at -pi==pi if strings exist
    do iz=0,1
       do iy=0,1
          t=ph(iy,iz)-ph0
          if(t>pi) then
             t=t-twopi
          else if(t<-pi) then
             t=t+twopi
          end if
          ph(iy,iz)=t
       end do
    end do
    
  end subroutine rotate_phase_2d
  
  function calc_chi2_string_2d(ph,y0,z0,u0,v0,sgn)
    ! Given string position (y0,z0), velocity (u0,v0) and its direction sgn,
    ! this function gives the value [sum_i |exp(phase_i)-exp(phase(y0,z0,y0,v0,sgn,th0;point_i))|^2]
    ! minimized over th0
    double precision,intent(in) :: ph(0:1,0:1),y0,z0,u0,v0
    integer,intent(in) :: sgn
    integer,parameter :: n=2**ndim
    integer iy,iz
    double precision calc_chi2_string_2d,beta,gamma,y,z,pht,cr,si
    
    beta=u0**2+v0**2
    gamma=1/sqrt(1-beta)
    beta=sqrt(beta)
    
    cr=0.d0
    si=0.d0
    do iy=0,1
       do iz=0,1
          y=(u0*(iy-y0)+v0*(iz-z0))*gamma/beta
          z=(-v0*(iy-y0)+u0*(iz-z0))/beta
          pht=sgn*aimag(log(cmplx(y,z)))
          pht=pht-ph(iy,iz)
          cr=cr+cos(pht)
          si=si+sin(pht)
       end do
    end do
    calc_chi2_string_2d=2.d0*(n-sqrt(cr**2+si**2))
    
  end function calc_chi2_string_2d
  
  subroutine position_string(ph,pos,vel,sgn,gof)
    use potential,only : twopi
    double precision,intent(in) :: ph(0:1,0:1)
    double precision,intent(out) :: pos(ndim),vel(ndim),gof
    integer,intent(out) :: sgn
    ! pos and vel are respectively the position and velocity of string
    ! sgn the direction of the string; +1 for clockwise and -1 for counter-clockwise 
    integer,parameter :: n=2**ndim
    double precision y0,z0,u0,v0,b0,vp0,chi2,chi2min 
    integer isgn,iy0,iz0,ib0,ivp0
    integer,parameter :: nyz0=10,nb0=10,nvp0=20
    double precision,parameter :: maxg=10.d0
    
    chi2min=2.d0*n
    pos(1:ndim)=0.d0
    vel(1:ndim)=0.d0
    sgn=0
    do isgn=-1,1,2
       do iy0=1,nyz0
          y0=(iy0-0.5d0)/nyz0
          do iz0=1,nyz0
             z0=(iz0-0.5d0)/nyz0
             do ib0=0,nb0
                b0=1.d0+(maxg-1)*ib0/nb0
                b0=sqrt(1-1/b0**2)
                do ivp0=0,nvp0-1
                   vp0=twopi*ivp0/nvp0
                   u0=b0*cos(vp0)
                   v0=b0*sin(vp0)
                   
                   chi2=calc_chi2_string_2d(ph,y0,z0,u0,v0,isgn)
                   if(chi2<chi2min) then
                      pos(1:2)=(/y0,z0/)
                      vel(1:2)=(/u0,v0/)
                      sgn=isgn
                      chi2min=chi2
                   end if

                end do
             end do
          end do
       end do
    end do
    
    if(sgn==0) call MpiStop ('error: something wrong in position_string')
    
    ! goodness of fit
    gof=exp(-chi2min/2)

  end subroutine position_string
  
  function calc_chi2_wall_2d(wid_dw,ph,y0,z0,vel,sgn)
    double precision,intent(in) :: wid_dw,ph(0:1,0:1),y0,z0,vel
    integer,intent(in) :: sgn
    double precision calc_chi2_wall_2d
    integer iy,iz
    double precision vec(1:ndim),pht,gm,dy,dz,d
    
    vec(1:ndim)=(/y0,z0/)-0.5d0
    vec(1:ndim)=sgn*vec(1:ndim)/sqrt(sum(vec(1:ndim)**2)) ! direction of wall
    gm=wid_dw/sqrt(1.d0-vel**2) ! gamma factor times the axion mass in units of lattice separation
    
    calc_chi2_wall_2d=0.d0
    do iy=0,1
       dy=iy-0.5d0
       do iz=0,1
          dz=iz-0.5d0
          d=dy*vec(1)+dz*vec(2)
          pht=4*atan(exp(gm*d)) ! sine-Gordon solution
          calc_chi2_wall_2d=calc_chi2_wall_2d+((pht-ph(iy,iz))/pi)**2
       end do
    end do
    
  end function calc_chi2_wall_2d
  
  subroutine position_wall(wid_dw,ph,pos,vel,sgn,gof,wa)
    ! conpute the position and velocity of wall by fitting with sine-Gordon solution
    use potential,only : pi,twopi
    double precision,intent(in) :: wid_dw
    double precision,intent(in) :: ph(0:1,0:1)
    double precision,intent(out) :: pos(ndim),vel,gof,wa
    integer,intent(out) :: sgn
    ! pos is point on wall closest to the origin of the cell
    ! sgn is the direction of wall along the wall direction (from the center of cell to the wall)
    ! and vel is the velocity along the unit vector of pos
    integer iy0,iz0,ib0,isgn
    double precision y0,z0,b0,chi2min,chi2
    integer,parameter :: nyz0=9,nb0=19 ! nyz0 should be odd to avoid /0
    integer,parameter :: n=2**ndim
    double precision,parameter :: maxg=5.d0
    
    chi2min=4.d0*n
    pos(1:2)=0.d0
    vel=0.d0
    sgn=0
    do isgn=-1,1,2
       do iy0=0,nyz0
          y0=dble(iy0)/nyz0
          do iz0=0,nyz0
             z0=dble(iz0)/nyz0

             do ib0=0,nb0 ! degeneracies between b0<->-b0
                b0=1.d0+(maxg-1)*abs(ib0)/nb0
                b0=sqrt(1-1/b0**2)*sign(1,ib0)
                
                chi2=calc_chi2_wall_2d(wid_dw,ph,y0,z0,b0,isgn)
                if(chi2<chi2min) then
                   pos(1:2)=(/y0,z0/)
                   vel=b0
                   sgn=isgn
                   chi2min=chi2
                end if
             end do

          end do
       end do
    end do
    
    if(sgn==0) call MpiStop ('error: something wrong in position_wall')
    
    ! goodness of fit
    gof=exp(-chi2min/2)
    
    ! wall area
    y0=abs(pos(1)-0.5d0)
    z0=abs(pos(2)-0.5d0)
    wa=abs(y0**2+z0**2-(y0+z0)/2)*sqrt(y0**2+z0**2)/y0/z0
    
  end subroutine position_wall

!!$  subroutine grouping_strings(first_string,wid_st,hor_size,num_strings)
!!$    type(string_list),pointer :: first_string
!!$    double precision,intent(in) :: wid_st,hor_size
!!$    integer,intent(out) :: num_strings
!!$
!!$    integer MpiRank,MpiSize,num,info
!!$    double precision th_dist
!!$    double precision,parameter :: pow_st=2.d0
!!$
!!$    call MpiStat(MpiRank,MpiSize)
!!$    
!!$    th_dist=(wid_st**pow_st*hor_size)**(1/(pow_st+1)) ! weighted geometric mean of string width and string separation
!!$    th_dist=max(def_th_dist,th_dist)
!!$    if(MpiRank==0) then
!!$       write(*,'(" size of string width in lattice units:",F)') wid_st
!!$       write(*,'(" size of string separation in lattice units:",F)') hor_size
!!$       write(*,'(" threshold for string grouping in lattice units:",F)') th_dist
!!$    end if
!!$
!!$    ! get local id in each MPI process
!!$    call grouping_local(first_string,th_dist,num)
!!$
!!$    ! get connections between different MPI processes
!!$    call grouping_global(first_string,th_dist,num)
!!$    num_strings=num
!!$
!!$  end subroutine grouping_strings
!!$
!!$  subroutine grouping_global(first_string,th_dist,num_strings)
!!$    type(string_list),pointer :: first_string
!!$    double precision,intent(in) :: th_dist
!!$    integer,intent(inout) :: num_strings
!!$
!!$    type(string_list),pointer :: tmp_string
!!$    integer MpiRank,MpiSize,MpiRank_out,MpiRank_in, &
!!$         num_out,num_in,info,num_conn,i,j,n,local_id,global_id, &
!!$         max_num_conn,max_num_strings
!!$    double precision dist
!!$    double precision,allocatable :: pos_out(:,:),pos_in(:,:)
!!$    integer,allocatable :: id_out(:),id_in(:),id_conn(:,:), &
!!$         array_num_conn(:),array_num_strings(:), &
!!$         connection(:,:,:),local_global(:,:)
!!$    logical new_group
!!$    ! buffter parameter for multiple-connection
!!$    integer,parameter :: mult_buff=100
!!$    integer isend,irecv,ssend(MPI_STATUS_SIZE),srecv(MPI_STATUS_SIZE)
!!$    
!!$    call MpiStat(MpiRank,MpiSize)
!!$    MpiRank_out=mod(MpiRank+1,MpiSize)
!!$    MpiRank_in=mod(MpiRank-1+MpiSize,MpiSize)
!!$    
!!$    ! get connections between different MPI processes
!!$    num_out=0
!!$    tmp_string=>first_string
!!$    do while(associated(tmp_string))
!!$       if(tmp_string%ixyz(3)+tmp_string%pos(3)>=max_1-th_dist) num_out=num_out+1
!!$       tmp_string=>tmp_string%next
!!$    end do
!!$    
!!$    call MPI_ISEND(num_out,1,MPI_INTEGER,MpiRank_out,0,MPI_COMM_WORLD,isend,info)
!!$    call MPI_IRECV(num_in,1,MPI_INTEGER,MpiRank_in,0,MPI_COMM_WORLD,irecv,info)
!!$    call MPI_WAIT(isend,ssend,info)
!!$    call MPI_WAIT(irecv,srecv,info)
!!$
!!$    if(num_out/=0) then
!!$       allocate(pos_out(3,num_out),id_out(num_out),stat=info)
!!$       if(info/=0) call MpiStop('error: cannot allocate pos_out or id_out')
!!$    end if
!!$    if(num_in/=0) then
!!$       allocate(pos_in(3,num_in),id_in(num_in),id_conn(2,num_in*mult_buff),stat=info)
!!$       if(info/=0) call MpiStop('error: cannot allocate pos_in, id_in or id_conn')
!!$    end if
!!$    
!!$    if(num_out/=0) then
!!$       i=0
!!$       tmp_string=>first_string
!!$       do while(associated(tmp_string))
!!$          if(tmp_string%ixyz(3)+tmp_string%pos(3)>=max_1-th_dist) then
!!$             i=i+1
!!$             if(i>num_out) call MpiStop('error: num_out is odd')
!!$             pos_out(1:3,i)=tmp_string%ixyz(1:3)+tmp_string%pos(1:3)
!!$             id_out(i)=tmp_string%loc_id
!!$          end if
!!$          tmp_string=>tmp_string%next
!!$       end do
!!$       if(i/=num_out) call MpiStop('error: something wrong with num_out')
!!$       call MPI_ISEND(pos_out,3*num_out,MPI_DOUBLE_PRECISION,MpiRank_out,0,MPI_COMM_WORLD,isend,info)
!!$    end if
!!$    if(num_in/=0) &
!!$         call MPI_IRECV(pos_in,3*num_in,MPI_DOUBLE_PRECISION,MpiRank_in,0,MPI_COMM_WORLD,irecv,info)
!!$    if(num_out/=0) call MPI_WAIT(isend,ssend,info)
!!$    if(num_in/=0) call MPI_WAIT(irecv,srecv,info)
!!$    
!!$    if(num_out/=0) call MPI_ISEND(id_out,num_out,MPI_INTEGER,MpiRank_out,0,MPI_COMM_WORLD,isend,info)
!!$    if(num_in/=0) call MPI_IRECV(id_in,num_in,MPI_INTEGER,MpiRank_in,0,MPI_COMM_WORLD,irecv,info)
!!$    if(num_out/=0) call MPI_WAIT(isend,ssend,info)
!!$    if(num_in/=0) call MPI_WAIT(irecv,srecv,info)
!!$    
!!$    deallocate(pos_out,id_out,stat=info)
!!$    
!!$    if(num_in/=0) then
!!$       id_conn(1:2,1:num_in*mult_buff)=0
!!$       tmp_string=>first_string
!!$       num_conn=0
!!$       do while(associated(tmp_string))
!!$          if(tmp_string%ixyz(3)+tmp_string%pos(3)<=min_1+th_dist) then
!!$             do i=1,num_in
!!$                call lat_dist_periodic(tmp_string%ixyz(1:3)+tmp_string%pos(1:3),pos_in(1:3,i),dist)
!!$                if(dist<=th_dist) then
!!$                   j=1
!!$                   do while(j<=num_conn)
!!$                      if(id_conn(1,j)==id_in(i).and.id_conn(2,j)==tmp_string%loc_id) exit ! known connection
!!$                      j=j+1
!!$                   end do
!!$                   if(j>num_conn) then ! new connection
!!$                      if(j/=num_conn+1) call MpiStop ('error: something is wrong in grouping_global')
!!$                      num_conn=num_conn+1
!!$                      if(num_conn>num_in*mult_buff) call MpiStop ('error: need larger mult_buff')
!!$                      id_conn(1:2,num_conn)=(/id_in(i),tmp_string%loc_id/)
!!$                   end if
!!$                end if
!!$             end do
!!$          end if
!!$          tmp_string=>tmp_string%next
!!$       end do
!!$    end if
!!$
!!$    deallocate(pos_in,id_in,stat=info)
!!$
!!$    ! asign global id
!!$    allocate(array_num_conn(0:MpiSize-1),array_num_strings(0:MpiSize-1),stat=info)
!!$#ifdef MPI
!!$    call MPI_ALLGATHER(num_conn,1,MPI_INTEGER,array_num_conn,1,MPI_INTEGER,MPI_COMM_WORLD,info)
!!$    call MPI_ALLGATHER(num_strings,1,MPI_INTEGER,array_num_strings,1,MPI_INTEGER,MPI_COMM_WORLD,info)
!!$#else
!!$    array_num_conn(0)=num_conn
!!$    array_num_strings(0)=num_strings
!!$#endif
!!$    
!!$    max_num_conn=maxval(array_num_conn)
!!$    max_num_strings=maxval(array_num_strings)
!!$    allocate(connection(2,max_num_conn,0:MpiSize-1),stat=info)
!!$    connection(:,:,:)=0
!!$    
!!$#ifdef MPI
!!$    call MPI_ALLGATHER(id_conn,2*max_num_conn,MPI_INTEGER, &
!!$         connection,2*max_num_conn,MPI_INTEGER,MPI_COMM_WORLD,info)
!!$#else
!!$    connection(1:2,1:max_num_connection,0)=id_conn(1:2,1:max_num_connection)
!!$#endif
!!$    deallocate(id_conn,stat=info)
!!$
!!$    allocate(local_global(max_num_strings,0:MpiSize-1),stat=info)
!!$    local_global(:,:)=0
!!$    
!!$    global_id=0
!!$    do
!!$
!!$       new_group=.false.
!!$       do i=0,MpiSize-1
!!$          do local_id=1,array_num_strings(i)
!!$             if(local_global(local_id,i)==0) then
!!$                new_group=.true.
!!$                global_id=global_id+1
!!$                local_global(local_id,i)=global_id
!!$                exit
!!$             end if
!!$          end do
!!$          if(new_group) exit
!!$       end do
!!$       if(.not.new_group) exit
!!$       
!!$       do 
!!$          
!!$          new_group=.false.
!!$
!!$          do i=0,MpiSize-1
!!$             MpiRank_out=mod(i+1,MpiSize)
!!$             MpiRank_in=mod(i-1+MpiSize,MpiSize)
!!$          
!!$             do local_id=1,array_num_strings(i)
!!$                if(local_global(local_id,i)==global_id) then
!!$                   ! check connections with lower block
!!$                   do j=1,array_num_conn(i)
!!$                      if(connection(2,j,i)==local_id) then
!!$                         if(local_global(connection(1,j,i),MpiRank_in)==0) then
!!$                            new_group=.true.
!!$                            local_global(connection(1,j,i),MpiRank_in)=global_id
!!$                         else if(local_global(connection(1,j,i),MpiRank_in)/=global_id) then
!!$                            call MpiStop ('error: incorrect asignment of global id with upper block')
!!$                         end if
!!$                      end if
!!$                   end do
!!$                   ! check connections with upper block
!!$                   do j=1,array_num_conn(MpiRank_out)
!!$                      if(connection(1,j,MpiRank_out)==local_id) then
!!$                         if(local_global(connection(2,j,MpiRank_out),MpiRank_out)==0) then
!!$                            new_group=.true.
!!$                            local_global(connection(2,j,MpiRank_out),MpiRank_out)=global_id
!!$                         else if(local_global(connection(2,j,MpiRank_out),MpiRank_out)/=global_id) then
!!$                            call MpiStop ('error: incorrect asignment of global id with lower block')
!!$                         end if
!!$                      end if
!!$                   end do
!!$                end if
!!$             end do
!!$          end do
!!$
!!$          if(.not.new_group) exit
!!$       end do
!!$       
!!$    end do
!!$    
!!$    deallocate(array_num_conn,array_num_strings,connection,stat=info)
!!$    
!!$    tmp_string=>first_string
!!$    do while(associated(tmp_string))
!!$       tmp_string%glo_id=local_global(tmp_string%loc_id,MpiRank)
!!$       tmp_string=>tmp_string%next
!!$    end do
!!$    num_strings=global_id
!!$    
!!$    deallocate(local_global,stat=info)
!!$    
!!$#ifdef MPI
!!$    call MPI_BARRIER(MPI_COMM_WORLD,info)
!!$#endif
!!$
!!$  end subroutine grouping_global
!!$
!!$  subroutine grouping_local(first_string,th_dist,local_num)
!!$    ! get local id in each MPI process
!!$    type(string_list),pointer :: first_string
!!$    double precision,intent(in) :: th_dist
!!$    integer,intent(out) :: local_num
!!$
!!$    double precision dist
!!$    integer MpiRank,MpiSize,local_id,info,n,j,i
!!$    type(string_list),pointer :: tmp_string,scan_string
!!$    logical new_group
!!$    integer,allocatable :: connection(:,:),tmp_local(:)
!!$
!!$    call MpiStat(MpiRank,MpiSize)
!!$    
!!$    local_id=0
!!$    do 
!!$       new_group=.true.
!!$       
!!$       tmp_string=>first_string
!!$       do while(associated(tmp_string))
!!$
!!$          if(new_group) then
!!$             if(tmp_string%loc_id==0) then ! new group
!!$                new_group=.false.
!!$                local_id=local_id+1
!!$                tmp_string%loc_id=local_id
!!$             end if
!!$          endif
!!$             
!!$          if(.not.new_group) then
!!$
!!$             if(tmp_string%loc_id==0) then ! inquery with current group
!!$                   
!!$                scan_string=>tmp_string%next
!!$                do while(associated(scan_string))
!!$                   if(scan_string%loc_id==local_id) then
!!$                      call lat_dist_periodic(tmp_string%ixyz(1:3)+tmp_string%pos(1:3), &
!!$                           scan_string%ixyz(1:3)+scan_string%pos(1:3),dist)
!!$                      if(dist<=th_dist) tmp_string%loc_id=local_id
!!$                   endif
!!$                   scan_string=>scan_string%next
!!$                end do
!!$
!!$             else if(tmp_string%loc_id==local_id) then
!!$
!!$                scan_string=>tmp_string%next
!!$                do while(associated(scan_string))
!!$                   if(scan_string%loc_id==0) then
!!$                      call lat_dist_periodic(tmp_string%ixyz(1:3)+tmp_string%pos(1:3), &
!!$                           scan_string%ixyz(1:3)+scan_string%pos(1:3),dist)
!!$                      if(dist<=th_dist) scan_string%loc_id=local_id
!!$                   endif
!!$                   scan_string=>scan_string%next
!!$                end do
!!$                
!!$             end if
!!$             
!!$          end if
!!$             
!!$          tmp_string=>tmp_string%next
!!$       end do
!!$
!!$       if(new_group) exit ! there is no more ungrouped members
!!$          
!!$    end do
!!$
!!$    ! check connections of groups
!!$    allocate(connection(2,2*local_id),stat=info)
!!$    n=0
!!$    tmp_string=>first_string
!!$    do while(associated(tmp_string))
!!$       scan_string=>tmp_string%next
!!$       do while(associated(scan_string))
!!$          
!!$          if(tmp_string%loc_id/=scan_string%loc_id) then
!!$             call lat_dist_periodic(tmp_string%ixyz(1:3)+tmp_string%pos(1:3), &
!!$                  scan_string%ixyz(1:3)+scan_string%pos(1:3),dist)
!!$             if(dist<=th_dist) then
!!$                new_group=.true.
!!$                do j=1,n
!!$                   if(connection(1,j)==max(tmp_string%loc_id,scan_string%loc_id) .and. &
!!$                        connection(2,j)==min(tmp_string%loc_id,scan_string%loc_id)) then
!!$                      new_group=.false.
!!$                      exit
!!$                   end if
!!$                end do
!!$                if(new_group) then
!!$                   n=n+1
!!$                   if(n>2*local_id) call MpiStop('something is wrong in connection')
!!$                   connection(1:2,n)=(/max(tmp_string%loc_id,scan_string%loc_id),&
!!$                        min(tmp_string%loc_id,scan_string%loc_id)/)
!!$                end if
!!$             end if
!!$          endif
!!$          
!!$          scan_string=>scan_string%next
!!$       end do
!!$       tmp_string=>tmp_string%next
!!$    end do
!!$
!!$    allocate(tmp_local(local_id),stat=info)
!!$    do i=1,local_id
!!$       tmp_local(i)=i
!!$    end do
!!$    do
!!$       new_group=.true.
!!$       do j=1,n
!!$          if(tmp_local(connection(1,j))/=tmp_local(connection(2,j))) then
!!$             new_group=.false.
!!$             tmp_local(connection(1:2,j))=minval(tmp_local(connection(1:2,j)))
!!$          end if
!!$       end do
!!$       if(new_group) exit
!!$    end do
!!$    
!!$    deallocate(connection,stat=info)
!!$    
!!$    n=local_id
!!$    local_id=0
!!$    do while(any(tmp_local(1:n)>local_id))
!!$       local_id=local_id+1
!!$       j=n
!!$       do i=1,n
!!$          if(tmp_local(i)>=local_id.and.tmp_local(i)<j) &
!!$               j=tmp_local(i)
!!$       end do
!!$       do i=1,n
!!$          if(tmp_local(i)==j) tmp_local(i)=local_id
!!$       end do
!!$    end do
!!$
!!$    tmp_string=>first_string
!!$    do while(associated(tmp_string))
!!$       tmp_string%loc_id=tmp_local(tmp_string%loc_id)
!!$       tmp_string=>tmp_string%next
!!$    end do
!!$    
!!$    deallocate(tmp_local,stat=info)
!!$    
!!$    !! check
!!$    tmp_string=>first_string
!!$    do while(associated(tmp_string))
!!$       scan_string=>tmp_string%next
!!$       do while(associated(scan_string))
!!$          call lat_dist_periodic(tmp_string%ixyz(1:3)+tmp_string%pos(1:3), &
!!$               scan_string%ixyz(1:3)+scan_string%pos(1:3),dist)
!!$          if(dist<=th_dist.and.tmp_string%loc_id/=scan_string%loc_id) call MpiStop('error: omitted connection')
!!$          scan_string=>scan_string%next
!!$       end do
!!$       tmp_string=>tmp_string%next
!!$    end do
!!$    !! check
!!$
!!$    local_num=local_id
!!$    
!!$  end subroutine grouping_local

  subroutine lat_dist_periodic(p1,p2,dist)
    double precision,intent(in) :: p1(3),p2(3)
    double precision,intent(out) :: dist
    double precision dp(3),s2
    integer i

    dp(1:3)=abs(p1(1:3)-p2(1:3))
    s2=0.d0
    do i=1,3
       s2=s2+min(dp(i),ngrid-dp(i))**2
    end do
    dist=sqrt(s2)
  end subroutine lat_dist_periodic

!!$  subroutine sort_points(first,wid_st)
!!$    type(string_list),pointer :: first
!!$    double precision,intent(in) :: wid_st
!!$    type(string_list),pointer :: pool,scan,tmp,tmp_first,closest
!!$    integer i,max_loc_id,n,m,MpiRank,MpiSize,info,num
!!$    double precision dist,th_dist,dist_min
!!$    logical new_point
!!$    double precision,parameter :: boost_th=2.d0
!!$
!!$    call MpiStat(MpiRank,MpiSize)
!!$    th_dist=max(def_th_dist,wid_st)*boost_th
!!$    
!!$    max_loc_id=0
!!$    scan=>first
!!$    num=0
!!$    do while(associated(scan))
!!$       if(scan%loc_id>max_loc_id) max_loc_id=scan%loc_id
!!$       scan=>scan%next
!!$       num=num+1
!!$    end do
!!$    
!!$    pool=>first
!!$    nullify(first)
!!$    
!!$    do i=1,max_loc_id
!!$
!!$       ! extract points with loc_id==i from pool
!!$       nullify(tmp_first)
!!$
!!$       n=0
!!$       do 
!!$          new_point=.false.
!!$          
!!$          if(.not.associated(pool)) exit
!!$
!!$          if(pool%loc_id==i) then
!!$             if(associated(tmp_first)) then
!!$                tmp%next=>pool
!!$                tmp=>tmp%next
!!$             else
!!$                tmp_first=>pool
!!$                tmp=>tmp_first
!!$             end if
!!$             pool=>pool%next
!!$             new_point=.true.
!!$          else
!!$             scan=>pool
!!$             do while(associated(scan%next))
!!$                if(scan%next%loc_id==i) then
!!$                   if(associated(tmp_first)) then
!!$                      tmp%next=>scan%next
!!$                      tmp=>tmp%next
!!$                   else
!!$                      tmp_first=>scan%next
!!$                      tmp=>tmp_first
!!$                   end if
!!$                   scan%next=>scan%next%next
!!$                   new_point=.true.
!!$                   exit
!!$                end if
!!$                scan=>scan%next
!!$             end do
!!$          end if 
!!$         
!!$          if(.not.new_point) exit
!!$
!!$          tmp%order=-1
!!$          n=n+1
!!$       end do
!!$       nullify(tmp%next)
!!$
!!$       ! reordering with greedy algorithm
!!$       tmp=>tmp_first
!!$       do m=1,n
!!$          if(.not.associated(tmp)) call MpiStop('error: tmp should be associated')
!!$          tmp%order=m
!!$          if(m==n) exit
!!$          
!!$          ! find the closest point
!!$          nullify(closest)
!!$          dist_min=1.d30
!!$          scan=>tmp
!!$          do while(associated(scan%next))
!!$             call lat_dist_periodic(tmp%pos(1:3)+tmp%ixyz(1:3), &
!!$                  scan%next%pos(1:3)+scan%next%ixyz(1:3),dist)
!!$             if(.not.associated(closest).or.dist<dist_min) then
!!$                closest=>scan
!!$                dist_min=dist
!!$             end if
!!$             scan=>scan%next
!!$          end do
!!$          
!!$          if(.not.associated(closest,tmp)) then
!!$             scan=>tmp%next
!!$             tmp%next=>closest%next
!!$             closest%next=>closest%next%next
!!$             tmp%next%next=>scan
!!$          end if
!!$          tmp=>tmp%next
!!$       end do
!!$       nullify(tmp%next)
!!$
!!$       ! reordering
!!$       
!!$
!!$       ! accumulate
!!$       if(associated(first)) then
!!$          tmp=>first
!!$          do while(associated(tmp%next))
!!$             tmp=>tmp%next
!!$          end do
!!$          tmp%next=>tmp_first
!!$       else
!!$          first=>tmp_first
!!$       end if
!!$
!!$    end do
!!$    
!!$!! test
!!$    tmp=>first
!!$    n=0
!!$    do while(associated(tmp))
!!$       n=n+1
!!$       write(MpiRank+200,'(4I,3E)') tmp%loc_id,tmp%order,n,num,tmp%ixyz(1:3)+tmp%pos(1:3)
!!$       if(associated(tmp%next)) then
!!$          if(tmp%loc_id/=tmp%next%loc_id) then
!!$             write(MpiRank+200,*);write(MpiRank+200,*)
!!$          end if
!!$       end if
!!$       tmp=>tmp%next
!!$    end do
!!$!! test
!!$
!!$#ifdef MPI
!!$    call MPI_BARRIER(MPI_COMM_WORLD,info)
!!$#endif
!!$  end subroutine sort_points

!!$  subroutine add_string_length(first,num_strings,hor_size,sl)
!!$    use potential,only : twopi
!!$    type(string_list),pointer :: first
!!$    integer,intent(in) :: num_strings
!!$    double precision,intent(in) :: hor_size
!!$    double precision,intent(out) :: sl(3)
!!$    double precision dsl(num_strings),ddsl(num_strings)
!!$    type(string_list),pointer :: tmp 
!!$    integer i,info
!!$
!!$    ddsl(1:num_strings)=0.d0
!!$    tmp=>first
!!$    do while(associated(tmp))
!!$       ddsl(tmp%glo_id)=ddsl(tmp%glo_id)+tmp%sl
!!$       tmp=>tmp%next
!!$    end do
!!$    
!!$#ifdef MPI
!!$    call MPI_ALLREDUCE(ddsl,dsl,num_strings,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
!!$#else
!!$    dsl(1:num_strings)=ddsl(1:num_strings)
!!$#endif
!!$
!!$    sl(1:3)=0.d0
!!$    do i=1,num_strings
!!$       sl(1)=sl(1)+dsl(i)
!!$       if(dsl(i)>hor_size) sl(2)=sl(2)+dsl(i)
!!$       if(dsl(i)>twopi*hor_size) sl(3)=sl(3)+dsl(i)
!!$    end do
!!$    
!!$  end subroutine add_string_length
  
  subroutine add_wall_area(first,wa)
    type(wall_list),pointer :: first
    double precision,intent(out) :: wa
    double precision dwa,ddwa
    type(wall_list),pointer :: tmp 
    integer info

    ddwa=0.d0
    tmp=>first
    do while(associated(tmp))
       ddwa=ddwa+tmp%wa
       tmp=>tmp%next
    end do
    
#ifdef MPI
    call MPI_ALLREDUCE(ddwa,dwa,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
#else
    dwa=ddwa
#endif
    
    wa=dwa

  end subroutine add_wall_area

end module defects_2d
