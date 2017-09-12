module defects_3d
  use settings
  use subroutines
  implicit none
  
  double precision,parameter :: string_phdiff=pi,removal_phdiff=7*pi/4 !3*pi/2
  
  type mask_list
     integer ixyz(ndim)
  end type mask_list
  
  type string_list
     integer loc_id,glo_id
     double precision pos(ndim),dir(ndim),gof,sl,vel(ndim)
  end type string_list
  
  type wall_list
     integer loc_id,glo_id
     double precision pos(ndim),dir(ndim),gof,wa
  end type wall_list
  
  integer(1),allocatable :: exist_defect(:,:,:)
  integer(1),parameter :: flag_string=1,flag_wall=2,flag_mask=4, &
       flag_defect=flag_string+flag_wall

  integer,parameter :: max_iter_string_pos=1000,num_simp_string_pos=4, &
       mult_string_pos=num_simp_string_pos*(max_iter_string_pos+1)
  
contains
  
  subroutine id_defects(idx,tau,sl,wa,svel)
    integer,intent(in) :: idx
    double precision,intent(in) :: tau
    double precision,intent(out) :: sl(3),wa(3),svel(3)
    
    integer,parameter :: n=2**ndim,n2=n/2
    integer,parameter :: &
         icx1(n2)=(/1,3,5,7/),icx2(n2)=(/2,4,6,8/), &
         icy1(n2)=(/1,2,5,6/),icy2(n2)=(/3,4,7,8/), &
         icz1(n2)=(/1,2,3,4/),icz2(n2)=(/5,6,7,8/)
    
    integer ix1,iy1,iz1,ix2,iy2,iz2,info,i,ixyz(n,ndim),ithread,ns,nw,nm, &
         nstot,nwtot,nmtot,nsvel,nsveltot,ixyz0(ndim),num_strings,num_walls
    double precision ph(n)
    character(len=charlen) file_string,file_wall,file_mask,str,str0
    type(string_list),allocatable :: arr_string(:)
    type(wall_list),allocatable :: arr_wall(:)
    type(mask_list),allocatable :: arr_mask(:)
    integer,allocatable :: arr_ns(:),arr_nw(:),arr_nm(:),arr_nstot(:)
    integer(8) nskip
    
    sl=0.d0
    wa=0.d0
    svel=0.d0
    
    if(.not.(id_strings.or.id_walls.or.spec_est)) return
    
    if(MpiRank==0) write(*,'(A)') ' *identifying strings and walls ...'
    
    allocate(exist_defect(ngrid,ngrid,local_0_offset+1:local_0_offset+local_0),stat=info)
    if(info/=0) call MpiStop('error: cannot allocate exist_defect')
    
    if(id_strings) allocate(arr_ns(num_threads),stat=info)
    if(id_walls) allocate(arr_nw(num_threads),stat=info)
    if(spec_est) allocate(arr_nm(num_threads),stat=info)
    
    !$OMP PARALLEL DEFAULT(SHARED) &
    !$OMP PRIVATE(iz1,iz2,iy1,iy2,ix1,ix2,ixyz,ph,ns,nw,nm,ithread)
#ifdef _OPENMP
    ithread=OMP_GET_THREAD_NUM()+1
#else
    ithread=1
#endif
    ns=0
    nw=0
    nm=0
    !$OMP DO SCHEDULE(STATIC) ! should be static
    do iz1=1+local_0_offset,local_0+local_0_offset
       iz2=iz1+1
       ixyz(icz1(1:n2),3)=iz1
       ixyz(icz2(1:n2),3)=iz2
       do iy1=1,ngrid
          iy2=mod(iy1,ngrid)+1
          ixyz(icy1(1:n2),2)=iy1
          ixyz(icy2(1:n2),2)=iy2
          do ix1=1,ngrid
             ix2=mod(ix1,ngrid)+1
             ixyz(icx1(1:n2),1)=ix1
             ixyz(icx2(1:n2),1)=ix2
             
             do i=1,n
                ph(i)=aimag(log(data1(ixyz(i,1),ixyz(i,2),ixyz(i,3)-local_0_offset+1)))
             end do
             call id_with_phase(ph,exist_defect(ix1,iy1,iz1),ns,nw,nm)
             
          end do
       end do
    end do
    !$OMP END DO
    if(id_strings) arr_ns(ithread)=ns
    if(id_walls) arr_nw(ithread)=nw
    if(spec_est) arr_nm(ithread)=nm
    !$OMP END PARALLEL
    
    nstot=0
    if(id_strings) then
       nstot=sum(arr_ns(1:num_threads))
       if(nstot>0) then
          allocate(arr_string(nstot),stat=info)
       else
          allocate(arr_string(1),stat=info)
       end if
       if(info/=0) call MpiStop('error: cannot allocate arr_string(nstot)')
    end if
    nwtot=0
    if(id_walls) then
       nwtot=sum(arr_nw(1:num_threads))
       if(nwtot>0) then
          allocate(arr_wall(nwtot),stat=info)
       else
          allocate(arr_wall(1),stat=info)
       end if
       if(info/=0) call MpiStop('error: cannot allocate arr_wall(nwtot)')
    end if
    nmtot=0
    if(spec_est) then
       nmtot=sum(arr_nm(1:num_threads))
       if(nmtot>0) then
          allocate(arr_mask(nmtot),stat=info)
       else
          allocate(arr_mask(1),stat=info)
       end if
       if(info/=0) call MpiStop('error: cannot allocate arr_mask(nmtot)')
    end if
    write(*,'(A,I4,A,3I)') '  MpiRank ',MpiRank,' has found string, wall and mask points: ',nstot,nwtot,nmtot
    
    if(id_strings) then
       allocate(arr_nstot(0:MpiSize-1),stat=info)
       if(info/=0) call MpiStop('error: cannot allocate arr_nstot')
       call MPI_ALLGATHER(nstot,1,MPI_INTEGER,arr_nstot,1,MPI_INTEGER,MPI_COMM_WORLD,info)
       
       nskip=0
       if(MpiRank/=0) nskip=sum(arr_nstot(0:MpiRank-1))
       call skip_rand_stream0(nskip,mult_string_pos)
       call put_rand_stream0
       nskip=sum(arr_nstot(MpiRank:MpiSize-1))
       call skip_rand_stream0(nskip,mult_string_pos)

       deallocate(arr_nstot,stat=info)
    end if

    !$OMP PARALLEL DEFAULT(SHARED) &
    !$OMP PRIVATE(iz1,iz2,iy1,iy2,ix1,ix2,ixyz,ph,ns,nw,nm,ithread,ixyz0,nskip)
#ifdef _OPENMP
    ithread=OMP_GET_THREAD_NUM()+1
#else
    ithread=1
#endif
    if(nstot>0) then
       ns=sum(arr_ns(1:ithread-1))
       nskip=ns
       call skip_rand_stream(ithread,nskip,mult_string_pos)
    end if
    if(nwtot>0) nw=sum(arr_nw(1:ithread-1))
    if(nmtot>0) nm=sum(arr_nm(1:ithread-1))
    !$OMP DO SCHEDULE(STATIC) ! should be static
    do iz1=local_0_offset+1,local_0_offset+local_0
       iz2=iz1+1
       ixyz(icz1(1:n2),3)=iz1
       ixyz(icz2(1:n2),3)=iz2
       do iy1=1,ngrid
          iy2=mod(iy1,ngrid)+1
          ixyz(icy1(1:n2),2)=iy1
          ixyz(icy2(1:n2),2)=iy2
          do ix1=1,ngrid
             ixyz0(1:ndim)=(/ix1,iy1,iz1/)
             
             if(IAND(exist_defect(ix1,iy1,iz1),flag_defect)/=0) then ! defects exist
                
                ix2=mod(ix1,ngrid)+1
                ixyz(icx1(1:n2),1)=ix1
                ixyz(icx2(1:n2),1)=ix2
                do i=1,n
                   ph(i)=aimag(log(data1(ixyz(i,1),ixyz(i,2),ixyz(i,3)-local_0_offset+1)))
                end do
                
                if(IAND(exist_defect(ix1,iy1,iz1),flag_wall)/=0) then ! walls exist
                   nw=nw+1
                   call position_wall(ph,arr_wall(nw)%pos,arr_wall(nw)%dir,arr_wall(nw)%gof,arr_wall(nw)%wa)
                   arr_wall(nw)%pos(1:ndim)=arr_wall(nw)%pos(1:ndim)+ixyz0(1:ndim)-0.5d0
                   arr_wall(nw)%loc_id=0
                   arr_wall(nw)%glo_id=0
                end if
                
                if(IAND(exist_defect(ix1,iy1,iz1),flag_string)/=0) then ! strings exist
                   ns=ns+1
                   call position_string(ph,arr_string(ns)%pos,arr_string(ns)%dir,arr_string(ns)%gof,arr_string(ns)%sl)
#if VEL==1
                   call velocity_string_YY(ixyz0,arr_string(ns)%pos,arr_string(ns)%vel)
#else
                   call velocity_string(ixyz0,arr_string(ns)%pos,arr_string(ns)%vel)
#endif
                   arr_string(ns)%pos(1:ndim)=arr_string(ns)%pos(1:ndim)+ixyz0(1:ndim)-0.5d0
                   arr_string(ns)%loc_id=0
                   arr_string(ns)%glo_id=0
                end if
                
             end if
             
             if(IAND(exist_defect(ix1,iy1,iz1),flag_mask)/=0) then ! mask
                nm=nm+1
                arr_mask(nm)%ixyz(1:ndim)=ixyz0(1:ndim)
             end if
             
          end do
       end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL
    
    if(id_strings) call put_rand_stream0
    write(*,'(A,I4,A)') '  MpiRank ',MpiRank,': done position identification'
    
    deallocate(exist_defect,stat=info)
    
    if(group_strings) call grouping_strings(nstot,arr_string,width_cr/dspa/tau,tau/dspa,num_strings)
    if(group_walls) call grouping_walls(nwtot,arr_wall,thickness_cr/dspa/tau,tau/dspa,num_walls)
    
    if(nstot>0) deallocate(arr_ns,stat=info)
    if(nwtot>0) deallocate(arr_nw,stat=info)
    if(nmtot>0) deallocate(arr_nm,stat=info)
    
    if(out_strings.or.out_walls.or.spec_est) then
       write(str0,*) MpiRank
       str0=adjustl(str0)
       write(str,*) idx
       str=adjustl(str)
       if(MpiRank==0) &
            write(*,'(A)') '  outputting string, wall (and mask)'
       
       if(out_strings) then
          file_string=trim(rootname)//'_string_'//trim(str)//'-'//trim(str0)//'.dat'
          open(100,file=trim(file_string),action='write',status='replace',iostat=info)
          if(info/=0) call MpiStop(' error: cannot create string.dat')
          do ns=1,nstot
             write(100,'(3F10.3,8F8.3,2I4)',iostat=info) &
                  arr_string(ns)%pos(1:3),arr_string(ns)%dir(1:3),arr_string(ns)%vel(1:3),arr_string(ns)%gof, &
                  arr_string(ns)%sl,arr_string(ns)%loc_id,arr_string(ns)%glo_id
          end do
          close(100)
       end if
       
       if(out_walls) then
          file_wall=trim(rootname)//'_wall_'//trim(str)//'-'//trim(str0)//'.dat'
          open(100,file=trim(file_wall),action='write',status='replace',iostat=info)
          if(info/=0) call MpiStop(' error: cannot create wall.dat')
          do nw=1,nwtot
             write(100,'(3F10.3,5F8.3,2I4)',iostat=info) &
                  arr_wall(nw)%pos(1:3),arr_wall(nw)%dir(1:3),arr_wall(nw)%gof,arr_wall(nw)%wa, &
                  arr_wall(ns)%loc_id,arr_wall(ns)%glo_id
          end do
          close(100)
       end if
       
       if(spec_est) then
          file_mask=trim(rootname)//'_mask_'//trim(adjustl(str))//'-'//trim(str0)//'.dat'
          open(100,file=trim(file_mask),action='write',status='replace',iostat=info)
          if(info/=0) call MpiStop(' error: cannot create mask.dat')
          do nm=1,nmtot
             write(100,'(3I6)',iostat=info) arr_mask(nm)%ixyz(1:3)
          end do
          close(100)
       end if
    end if
    
    if(id_strings) then
       call add_string_length(nstot,arr_string,num_strings,tau/dspa,sl,svel)
       deallocate(arr_string)
    end if
    if(id_walls) then
       call add_wall_area(nwtot,arr_wall,num_walls,tau/dspa,wa)
       deallocate(arr_wall)
    endif
    if(spec_est) deallocate(arr_mask)

#ifdef MPI
    call MPI_BARRIER(MPI_COMM_WORLD,info)
    if(MpiRank==0) then
       if(out_strings) then
          file_string=trim(rootname)//'_string_'//trim(str)
          call execute_command_line('cat '//trim(file_string)//'-*.dat >'//trim(file_string)//'.dat; rm -rf '//trim(file_string)//'-*.dat')
       end if
       if(out_walls) then
          file_wall=trim(rootname)//'_wall_'//trim(str)
          call execute_command_line('cat '//trim(file_wall)//'-*.dat >'//trim(file_wall)//'.dat; rm -rf '//trim(file_wall)//'-*.dat')
       end if
!!$       if(spec_est) then
!!$          file_mask=trim(rootname)//'_mask_'//trim(str)
!!$          call execute_command_line('cat '//trim(file_mask)//'-*.dat >'//trim(file_mask)//'.dat; rm -rf '//trim(file_mask)//'-*.dat')
!!$       end if
    end if
    call MPI_BARRIER(MPI_COMM_WORLD,info)
#endif
    
    if(MpiRank==0) write(*,'(A)') ' *...done'
    
  end subroutine id_defects
  
  subroutine id_with_phase(ph,flag,ns,nw,nm)
    use potential,only : twopi
    integer,parameter :: n=2**ndim
    double precision,intent(in) :: ph(n)
    integer(1),intent(out) :: flag
    integer,intent(inout) :: ns,nw,nm
    integer i,j
    double precision pht(n),tmax,t
    
    pht(1:n)=ph(1:n)
    
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
    
    flag=0
    t=pht(1)-pht(n)+twopi
    tmax=maxval(pht(2:n)-pht(1:n-1))
    if(id_walls.and.tmax>t) then ! largest separation does not include +-pi
       flag=flag+flag_wall
       nw=nw+1
    else
       tmax=max(t,tmax)
    end if
    if(id_strings.and.tmax<string_phdiff) then ! string exists
       flag=flag+flag_string
       ns=ns+1
    end if
    if(spec_est.and.tmax<removal_phdiff) then
       flag=flag+flag_mask
       nm=nm+1
    end if
    
  end subroutine id_with_phase
  
  subroutine rotate_phase(ph) ! set the center of the largest phase sparation to +-pi
    use potential,only : pi,twopi
    integer,parameter :: n=2**ndim
    double precision,intent(inout) :: ph(n)
    integer i,j
    double precision t,ph0,pht(0:n),tmax
    
    pht(1:n)=ph(1:n)
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
    pht(0)=pht(n)-twopi
    
    tmax=0.d0
    do j=1,n
       t=pht(j)-pht(j-1)
       if(tmax<t) then
          tmax=t
          i=j
       end if
    end do
    ph0=(pht(i)+pht(i-1))/2
    
    ph(1:n)=modulo(ph(1:n)-ph0+pi,twopi)-pi
    
  end subroutine rotate_phase
  
  function calc_chi2_string(ph,p)
    use potential,only : pi,twopi
    integer,parameter :: n=2**ndim
    double precision,intent(in) :: ph(n),p(4)
    double precision calc_chi2_string,x,y,z,du,dv, &
         cp,sp,cm,sm,a,b,d,mu,phi,u0,v0
    integer ix,iy,iz,i
    
    mu=p(1)
    phi=p(2)
    u0=p(3)*cos(p(4))
    v0=p(3)*sin(p(4))
    
    cm=mu
    sm=sqrt(1.d0-mu*mu)
    cp=cos(phi)
    sp=sin(phi)
    
    a=0.d0
    b=0.d0
    do iz=0,1
       z=iz-0.5d0
       do iy=0,1
          y=iy-0.5d0
          do ix=0,1
             x=ix-0.5d0
             i=(iz*2+iy)*2+ix+1
             
             du=-u0+cm*(x*cp+y*sp)-z*sm
             dv=-v0-x*sp+y*cp
             
             d=ph(i)-atan2(dv,du)
             if(d>=pi) then
                d=d-twopi
             else if(d<-pi) then
                d=d+twopi
             end if
             
             a=a+d**2
             b=b+d
             
          end do
       end do
    end do
    a=a/n
    b=b/n
    calc_chi2_string=a-b**2
    
  end function calc_chi2_string
  
  subroutine position_string(ph,pos,dir,gof,sl)
    use potential,only : pi
    integer,parameter :: n=2**ndim
    double precision,intent(inout) :: ph(n)
    double precision,intent(out) :: pos(ndim),dir(ndim),gof,sl
    ! pos is point on string closest to the origin of the cell
    ! dir is the direction of the string
    double precision mu,phi,u0,v0,chi2
    
    double precision,parameter :: chi2max=8*pi*pi,tol=0.05d0
    double precision,parameter :: &
         max_mu=1.d0,min_mu=-1.d0,max_phi=pi,min_phi=-pi, &
         max_uv=sqrt(3.d0)/2,min_uv=0.d0,max_rot=pi,min_rot=-pi
    double precision x(mult_string_pos)
    
    double precision,parameter :: &
         max_simp(num_simp_string_pos)=(/max_mu,max_phi,max_uv,max_rot/), &
         min_simp(num_simp_string_pos)=(/min_mu,min_phi,min_uv,min_rot/)
    
    double precision p_simp(num_simp_string_pos,0:num_simp_string_pos),chi2_simp(0:num_simp_string_pos), &
         c_simp(num_simp_string_pos),t_simp(num_simp_string_pos),s_simp(num_simp_string_pos),chi2s
    integer i,j,i0,n0,ithread
    
    !!$call rotate_phase(ph)
    
    ithread=1
#ifdef _OPENMP
    if(OMP_IN_PARALLEL()) ithread=OMP_GET_THREAD_NUM()+1
#endif
    
    call gen_rand_uniform(ithread,x,mult_string_pos)
    ! initial points for simplex method
    n0=-1
    do i=0,max_iter_string_pos
       j=i*num_simp_string_pos
       t_simp(1)=max_mu*x(j+1)+min_mu*(1-x(j+1))
       t_simp(2)=max_phi*x(j+2)+min_phi*(1-x(j+2))
       t_simp(3)=max_uv*x(j+3)+min_uv*(1-x(j+3))
       t_simp(4)=max_rot*x(j+4)+min_rot*(1-x(j+4))
       
       chi2=calc_chi2_string(ph,t_simp)
       if(n0<=0) then
          n0=n0+1
          chi2_simp(n0)=chi2
          p_simp(1:num_simp_string_pos,n0)=t_simp(1:num_simp_string_pos)
          cycle
       end if
       
       if(n0<num_simp_string_pos.or.chi2<=chi2_simp(n0)) then
          n0=min(n0+1,num_simp_string_pos)
          
          i0=n0
          do while(chi2<=chi2_simp(i0-1))
             if(i0<=num_simp_string_pos) then
                chi2_simp(i0)=chi2_simp(i0-1)
                p_simp(1:num_simp_string_pos,i0)=p_simp(1:num_simp_string_pos,i0-1)
             end if
             i0=i0-1
             if(i0==0) exit
          end do
          
          chi2_simp(i0)=chi2
          p_simp(1:num_simp_string_pos,i0)=t_simp(1:num_simp_string_pos)
       end if
       
    end do
    
    ! simplex method
    i=0
    do while(chi2_simp(0)>tol.and.i<max_iter_string_pos)!.and.chi2_simp(num_simp_string_pos)/chi2_simp(0)>1+tol/10) ! not very much tested; does this work well?
       i=i+1
       if(i>max_iter_string_pos) exit!call MpiStop('error: simplex method does not converge')
       
       c_simp(1:num_simp_string_pos)=0.d0
       do i0=0,num_simp_string_pos-1
          c_simp(1:num_simp_string_pos)=c_simp(1:num_simp_string_pos)+p_simp(1:num_simp_string_pos,i0)
       end do
       c_simp(1:num_simp_string_pos)=c_simp(1:num_simp_string_pos)/num_simp_string_pos
       
       ! reflection
       t_simp(1:num_simp_string_pos)=2*c_simp(1:num_simp_string_pos)-p_simp(1:num_simp_string_pos,num_simp_string_pos)
       if(any(t_simp>max_simp).or.any(t_simp<min_simp)) then
          chi2=chi2max
       else
          chi2=calc_chi2_string(ph,t_simp)
       end if
       
       if(chi2<chi2_simp(num_simp_string_pos)) then
          
          ! expansion
          s_simp(1:num_simp_string_pos)=3*c_simp(1:num_simp_string_pos)-4*p_simp(1:num_simp_string_pos,num_simp_string_pos)
          if(any(s_simp>max_simp).or.any(s_simp<min_simp)) then
             chi2s=chi2max
          else
             chi2s=calc_chi2_string(ph,s_simp)
          end if
          
          if(chi2<=chi2s) then
             chi2_simp(num_simp_string_pos)=chi2
             p_simp(1:num_simp_string_pos,num_simp_string_pos)=t_simp(1:num_simp_string_pos)
          else
             chi2_simp(num_simp_string_pos)=chi2s
             p_simp(1:num_simp_string_pos,num_simp_string_pos)=s_simp(1:num_simp_string_pos)
          end if
          
       else
          
          s_simp(1:num_simp_string_pos)=1.5d0*c_simp(1:num_simp_string_pos)-0.5d0*p_simp(1:num_simp_string_pos,num_simp_string_pos)
          if(any(s_simp>max_simp).or.any(s_simp<min_simp)) then
             chi2s=chi2max
          else
             chi2s=calc_chi2_string(ph,s_simp)
          end if
          
          if(chi2s<=chi2_simp(num_simp_string_pos)) then
             ! contraction
             chi2_simp(num_simp_string_pos)=chi2s
             p_simp(1:num_simp_string_pos,num_simp_string_pos)=s_simp(1:num_simp_string_pos)
          else
             ! simplex contraction
             do i0=1,num_simp_string_pos
                p_simp(1:num_simp_string_pos,i0)=0.5d0*(p_simp(1:num_simp_string_pos,i0)+p_simp(1:num_simp_string_pos,0))
                chi2_simp(i0)=calc_chi2_string(ph,p_simp(1:num_simp_string_pos,i0))
             end do
          end if
          
       end if

       ! sort
       do i0=num_simp_string_pos-1,0,-1
          do n0=0,i0
             if(chi2_simp(n0)>chi2_simp(n0+1)) then
                
                t_simp(1:num_simp_string_pos)=p_simp(1:num_simp_string_pos,n0)
                p_simp(1:num_simp_string_pos,n0)=p_simp(1:num_simp_string_pos,n0+1)
                p_simp(1:num_simp_string_pos,n0+1)=t_simp(1:num_simp_string_pos)
                
                chi2=chi2_simp(n0)
                chi2_simp(n0)=chi2_simp(n0+1)
                chi2_simp(n0+1)=chi2

             end if
          end do
       end do
              
    end do
    
    mu=p_simp(1,0)
    phi=p_simp(2,0)
    u0=p_simp(3,0)*cos(p_simp(4,0))
    v0=p_simp(3,0)*sin(p_simp(4,0))
    chi2=chi2_simp(0)
    
    ! once mu, phi, u0, v0 are determined, one can calculate string position and direction
    dir(1)=sqrt(1.d0-mu*mu)*cos(phi)
    dir(2)=sqrt(1.d0-mu*mu)*sin(phi)
    dir(3)=mu
    
    pos(1)=u0*mu*cos(phi)-v0*sin(phi)
    pos(2)=u0*mu*sin(phi)+v0*cos(phi)
    pos(3)=-u0*sqrt(1-mu*mu)
    
    ! goodness of fit
    gof=exp(-chi2/2)
    
    ! string length
    sl=string_length(pos,dir)
    
  end subroutine position_string
  
  subroutine velocity_string(ixyz,pos,vel)
    use subroutines,only : check_syminv
    ! slightly different from Yamaguchi & Yokoyama (2003)
    integer,intent(in) :: ixyz(ndim)
    double precision,intent(in) :: pos(ndim)
    double precision,intent(out) :: vel(ndim)
    integer,parameter :: np=2**ndim,nmat=(ndim+1)*(ndim+2)/2
    integer ix,iy,iz,i,j,k,ii,jj,kk,info
    complex(kind(0d0)) w(ndim+1),wdot(ndim+1), & ! zeroth and first derivatives
         w0(ndim+1),wdot0(ndim+1)
    double precision d(ndim+1),mat(nmat),mat0(nmat),u(ndim),v(ndim),det,a,b,c,x,y,vel0(ndim),vel1(ndim)
    double precision,parameter :: tol=1.d-5

    w(1:ndim+1)=cmplx(0.d0)
    wdot(1:ndim+1)=cmplx(0.d0)
    mat(1:nmat)=0.d0
    
    d(1)=1.d0
    do k=0,1
       iz=ixyz(3)+k
       d(4)=k-pos(3)-0.5d0
       
       do j=0,1
          iy=ixyz(2)+j
          if(iy>ngrid) iy=iy-ngrid
          d(3)=j-pos(2)-0.5d0
          
          do i=0,1
             ix=ixyz(1)+i
             if(ix>ngrid) ix=ix-ngrid
             d(2)=i-pos(1)-0.5d0
             
             w(1:ndim+1)=w(1:ndim+1)+d(1:ndim+1)*data1(ix,iy,iz-local_0_offset+1)
             wdot(1:ndim+1)=wdot(1:ndim+1)+d(1:ndim+1)*data2(ix,iy,iz-local_0_offset+1)
             call dspr('U',ndim+1,1.d0,d,1,mat)
             
          end do
       end do
    end do
    
    w(1:ndim+1)=w(1:ndim+1)/np
    wdot(1:ndim+1)=wdot(1:ndim+1)/np
    mat(1:nmat)=mat(1:nmat)/np
    
    ! inversion of matrix
    mat0(1:nmat)=mat(1:nmat)
    call dpptrf('U',ndim+1,mat0,info)
    call dpptri('U',ndim+1,mat0,info)
    call check_syminv(ndim+1,mat,mat0,tol)
    
    w0(1:ndim+1)=cmplx(0.d0)
    wdot0(1:ndim+1)=cmplx(0.d0)
    do ii=1,ndim+1
       do jj=1,ndim+1
          kk=max(ii,jj)
          kk=kk*(kk-1)/2+min(ii,jj)
          w0(ii)=w0(ii)+mat0(kk)*w(jj)
          wdot0(ii)=wdot0(ii)+mat0(kk)*wdot(jj)
       end do
    end do
    
    x=dble(wdot0(1)) ! Re[phidot]
    y=aimag(wdot0(1)) ! Im[phidot]
    u(1:ndim)=dble(w0(2:ndim+1)) ! Re[del phi]
    v(1:ndim)=aimag(w0(2:ndim+1)) ! Im[del phi]
    a=sum(u(1:ndim)**2)
    b=sum(u(1:ndim)*v(1:ndim)) ! Yamaguchi and Yokoyama (2003) assumes this vanishes
    c=sum(v(1:ndim)**2)
    det=a*c-b**2

!!$ Yamaguchi & Yokoyama 2002
    vel1(1:ndim)=x*v(1:ndim)-y*u(1:ndim)
    vel1(1:ndim)=vel1(1:ndim)/sqrt(det)*dspa
!!$ Yamaguchi & Yokoyama 2002
    vel0(1:ndim)=(x*c-y*b)*u(1:ndim)+(-x*b+y*a)*v(1:ndim)
    vel0(1:ndim)=-vel0(1:ndim)/det*dspa

#if VEL==1
    vel(1:ndim)=vel1(1:ndim)
#else
    vel(1:ndim)=vel0(1:ndim)
#endif

!!$ test
!!$    write(MpiRank+1000,'(3E12.4)') sqrt(sum(vel1(1:ndim)**2)),sqrt(sum(vel0(1:ndim)**2)),(sum(vel1(1:ndim)**2)-sum(vel0(1:ndim)**2))/sqrt(sum(vel0(1:ndim)**2)*sum(vel1(1:ndim)**2))
!!$ test

  end subroutine velocity_string
  
  subroutine velocity_string_YY0(ixyz,pos,vel)
    integer,intent(in) :: ixyz(ndim)
    double precision,intent(in) :: pos(ndim)
    double precision,intent(out) :: vel(ndim)
    
    integer,parameter :: num_pla=3,max_npen=2*num_pla
    double precision,parameter :: tol=1d-5
    integer info,ipla,iuvw(ndim),ixyz0(ndim),i,j,npen,k,l1,l2
    complex(kind(0d0)) vphi(0:2),dvphi(0:ndim)
    double precision mat(1:6),invmat(1:6),w(0:2),d(0:2,2),a11,a12,a22,det,du,dv,vpen(ndim,max_npen),weight
    
    ! matrix used for least chi^2 fit to phi and partial phi
    mat=0d0
    do i=0,1
       do j=0,1
          w(0:2)=(/1d0,i-0.5d0,j-0.5d0/)
          call dspr('U',3,1d0,w(0:2),1,mat)
       end do
    end do
    invmat(1:6)=mat(1:6)
    call dpptrf('U',3,invmat,info)
    call dpptri('U',3,invmat,info)
    call check_syminv(3,invmat,mat,tol)

    npen=0
    do ipla=1,num_pla ! each of three cells
       iuvw(1:ndim)=cshift((/1,2,3/),SHIFT=ipla-1)

       do k=0,1

          vphi(0:2)=0d0
          do i=0,1
             do j=0,1
                
                ixyz0(1:ndim)=ixyz(1:ndim)
                ixyz0(iuvw(1))=ixyz0(iuvw(1))+i
                ixyz0(iuvw(2))=ixyz0(iuvw(2))+j
                ixyz0(iuvw(3))=ixyz0(iuvw(3))+k
                if(ixyz0(1)>ngrid) ixyz0(1)=ixyz0(1)-ngrid
                if(ixyz0(2)>ngrid) ixyz0(2)=ixyz0(2)-ngrid
                
                vphi(0:2)=vphi(0:2)+(/1d0,i-0.5d0,j-0.5d0/)*data1(ixyz0(1),ixyz0(2),ixyz0(3)-local_0_offset+1)
                
             end do
          end do
          
          w(0:2)=dble(vphi(0:2))
          call dspmv('U',3,1d0,invmat,w(0:2),1,0d0,d(0:2,1),1)
          w(0:2)=aimag(vphi(0:2))
          call dspmv('U',3,1d0,invmat,w(0:2),1,0d0,d(0:2,2),1)
          
          ! compute the penetration point
          det=d(1,1)*d(2,2)-d(2,1)*d(1,2)
          du=-( d(2,2)*d(0,1)-d(2,1)*d(0,2))/det
          dv=-(-d(1,2)*d(0,1)+d(1,1)*d(0,2))/det
          
          if(abs(du)<0.5d0.and.abs(dv)<0.5d0) then ! string actually penetrates the plaquette
             ! string velocity
             npen=npen+1
             
             dvphi(0:ndim)=0
             do i=0,1
                do j=0,1
                   
                   ! bilinear interpolation
                   weight=(0.5d0+sign(du,i-0.5d0))*(0.5d0+sign(dv,j-0.5d0))
                   
                   ixyz0(1:ndim)=ixyz(1:ndim)
                   ixyz0(iuvw(1))=ixyz0(iuvw(1))+i
                   ixyz0(iuvw(2))=ixyz0(iuvw(2))+j
                   ixyz0(iuvw(3))=ixyz0(iuvw(3))+k
                   if(ixyz0(1)>ngrid) ixyz0(1)=ixyz0(1)-ngrid
                   if(ixyz0(2)>ngrid) ixyz0(2)=ixyz0(2)-ngrid
                   
                   ! dt
                   dvphi(0)=dvphi(0)+weight*data2(ixyz0(1),ixyz0(2),ixyz0(3)-local_0_offset+1)
                   ! dx
                   l1=ixyz0(1)+1
                   if(l1>ngrid) l1=l1-ngrid
                   l2=ixyz0(1)-1
                   if(l2<1) l2=l2+ngrid
                   dvphi(1)=dvphi(1)+weight*(data1(l1,ixyz0(2),ixyz0(3)-local_0_offset+1) &
                        -data1(l2,ixyz0(2),ixyz0(3)-local_0_offset+1))/2
                   ! dy
                   l1=ixyz0(2)+1
                   if(l1>ngrid) l1=l1-ngrid
                   l2=ixyz0(2)-1
                   if(l2<1) l2=l2+ngrid
                   dvphi(2)=dvphi(2)+weight*(data1(ixyz0(1),l1,ixyz0(3)-local_0_offset+1) &
                        -data1(ixyz0(1),l2,ixyz0(3)-local_0_offset+1))/2
                   ! dz
                   l1=ixyz0(3)+1
                   if(l1>local_0_offset+local_0+1) l1=ixyz0(3)
                   l2=ixyz0(3)-1
                   dvphi(3)=dvphi(3)+weight*(data1(ixyz0(1),ixyz0(2),l1-local_0_offset+1) &
                        -data1(ixyz0(1),ixyz0(2),l2-local_0_offset+1))/(l1-l2)
                end do
             end do
             
             a11=sum(dble(dvphi(1:ndim))**2)
             a12=sum(dble(dvphi(1:ndim))*aimag(dvphi(1:3)))
             a22=sum(aimag(dvphi(1:ndim))**2)
             det=a11*a22-a12**2
             vpen(1:ndim,npen)=dble(dvphi(0))*(a22*dble(dvphi(1:ndim))-a12*aimag(dvphi(1:ndim))) &
                  +aimag(dvphi(0))*(a11*aimag(dvphi(1:ndim))-a12*dble(dvphi(1:ndim)))
             vpen(1:ndim,npen)=-vpen(1:ndim,npen)/det*dspa
             
             !! debug
             write(*,'(A,7I5,6E10.2)') '#YY',MpiRank,ixyz(1:ndim),npen,ipla,k,du,dv,vpen(1:ndim,npen),sum(vpen(1:ndim,npen)**2)
             !! debug
             
          end if
       end do
    end do
    
    if(npen/=0) then
       vel(1:ndim)=0
       do i=1,npen
          vel(1:ndim)=vel(1:ndim)+vpen(1:ndim,npen)
       end do
       vel(1:ndim)=vel(1:ndim)/npen
    else
       !! debug
       vel(1:ndim)=1.d-10
       !!call MpiStop('error: something wrong in velocity_string_YY')
       !! debug
    end if
    
  end subroutine velocity_string_YY0

  subroutine velocity_string_YY(ixyz,pos,vel)
    integer,intent(in) :: ixyz(ndim)
    double precision,intent(in) :: pos(ndim)
    double precision,intent(out) :: vel(ndim)
    
    integer,parameter :: num_pla=3,max_npen=2*num_pla
    double precision,parameter :: tol=1d-5
    integer info,ipla,iuvw(ndim),ixyz0(ndim),i,j,npen,k,m,n,l1,l2
    complex(kind(0d0)) vphi(0:2),dvphi(0:ndim)
    double precision ph,mat(1:6),invmat(1:6),w(0:2),d(0:2,2),a11,a12,a22,det,du,dv,vpen(ndim,max_npen),weight
    logical quad(4)
    
    ! matrix used for least chi^2 fit to phi and partial phi
    mat=0d0
    do i=0,1
       do j=0,1
          w(0:2)=(/1d0,i-0.5d0,j-0.5d0/)
          call dspr('U',3,1d0,w(0:2),1,mat)
       end do
    end do
    invmat(1:6)=mat(1:6)
    call dpptrf('U',3,invmat,info)
    call dpptri('U',3,invmat,info)
    call check_syminv(3,invmat,mat,tol)
    
    npen=0
    do ipla=1,num_pla ! each of three cells
       iuvw(1:ndim)=cshift((/1,2,3/),SHIFT=ipla-1)
       
       do k=0,1
          
          ! check string if string penetrates plaquette
          n=0
          quad(1:4)=.false.
          do i=0,1
             do j=0,1
                
                ixyz0(1:ndim)=ixyz(1:ndim)
                ixyz0(iuvw(1))=ixyz0(iuvw(1))+i
                ixyz0(iuvw(2))=ixyz0(iuvw(2))+j
                ixyz0(iuvw(3))=ixyz0(iuvw(3))+k
                if(ixyz0(1)>ngrid) ixyz0(1)=ixyz0(1)-ngrid
                if(ixyz0(2)>ngrid) ixyz0(2)=ixyz0(2)-ngrid
                
                ph=aimag(log(data1(ixyz0(1),ixyz0(2),ixyz0(3)-local_0_offset+1)))
                m=ceiling(2*ph/pi)+2
                if(m<1.or.m>4) call MpiStop('error: wrong quadrant')
                if(quad(m)==.false.) then
                   quad(m)=.true.
                   n=n+1
                end if
                
             end do
          end do
          
          if(n<3) cycle ! string doesn't exist
          
          npen=npen+1
          ! compute the string position
          vphi(0:2)=0d0
          do i=0,1
             do j=0,1
                
                ixyz0(1:ndim)=ixyz(1:ndim)
                ixyz0(iuvw(1))=ixyz0(iuvw(1))+i
                ixyz0(iuvw(2))=ixyz0(iuvw(2))+j
                ixyz0(iuvw(3))=ixyz0(iuvw(3))+k
                if(ixyz0(1)>ngrid) ixyz0(1)=ixyz0(1)-ngrid
                if(ixyz0(2)>ngrid) ixyz0(2)=ixyz0(2)-ngrid
                
                vphi(0:2)=vphi(0:2)+(/1d0,i-0.5d0,j-0.5d0/)*data1(ixyz0(1),ixyz0(2),ixyz0(3)-local_0_offset+1)
                
             end do
          end do
          
          w(0:2)=dble(vphi(0:2))
          call dspmv('U',3,1d0,invmat,w(0:2),1,0d0,d(0:2,1),1)
          w(0:2)=aimag(vphi(0:2))
          call dspmv('U',3,1d0,invmat,w(0:2),1,0d0,d(0:2,2),1)
          
          ! compute the penetration point
          det=d(1,1)*d(2,2)-d(2,1)*d(1,2)
          du=-( d(2,2)*d(0,1)-d(2,1)*d(0,2))/det
          dv=-(-d(1,2)*d(0,1)+d(1,1)*d(0,2))/det
          ! force to fit in the plaquette
          du=min(du,0.5d0)
          du=max(du,-0.5d0)
          dv=min(dv,0.5d0)
          dv=max(dv,-0.5d0)

          ! string velocity
          dvphi(0:ndim)=0
          do i=0,1
             do j=0,1
                
                ! bilinear interpolation
                weight=(0.5d0+sign(du,i-0.5d0))*(0.5d0+sign(dv,j-0.5d0))
                
                ixyz0(1:ndim)=ixyz(1:ndim)
                ixyz0(iuvw(1))=ixyz0(iuvw(1))+i
                ixyz0(iuvw(2))=ixyz0(iuvw(2))+j
                ixyz0(iuvw(3))=ixyz0(iuvw(3))+k
                if(ixyz0(1)>ngrid) ixyz0(1)=ixyz0(1)-ngrid
                if(ixyz0(2)>ngrid) ixyz0(2)=ixyz0(2)-ngrid
                
                ! dt
                dvphi(0)=dvphi(0)+weight*data2(ixyz0(1),ixyz0(2),ixyz0(3)-local_0_offset+1)
                ! dx
                l1=ixyz0(1)+1
                if(l1>ngrid) l1=l1-ngrid
                l2=ixyz0(1)-1
                if(l2<1) l2=l2+ngrid
                dvphi(1)=dvphi(1)+weight*(data1(l1,ixyz0(2),ixyz0(3)-local_0_offset+1) &
                     -data1(l2,ixyz0(2),ixyz0(3)-local_0_offset+1))/2
                ! dy
                l1=ixyz0(2)+1
                if(l1>ngrid) l1=l1-ngrid
                l2=ixyz0(2)-1
                if(l2<1) l2=l2+ngrid
                dvphi(2)=dvphi(2)+weight*(data1(ixyz0(1),l1,ixyz0(3)-local_0_offset+1) &
                     -data1(ixyz0(1),l2,ixyz0(3)-local_0_offset+1))/2
                ! dz
                l1=ixyz0(3)+1
                if(l1>local_0_offset+local_0+1) l1=ixyz0(3)
                l2=ixyz0(3)-1
                dvphi(3)=dvphi(3)+weight*(data1(ixyz0(1),ixyz0(2),l1-local_0_offset+1) &
                     -data1(ixyz0(1),ixyz0(2),l2-local_0_offset+1))/(l1-l2)
             end do
          end do
          
          a11=sum(dble(dvphi(1:ndim))**2)
          a12=sum(dble(dvphi(1:ndim))*aimag(dvphi(1:3)))
          a22=sum(aimag(dvphi(1:ndim))**2)
          det=a11*a22-a12**2
          vpen(1:ndim,npen)=dble(dvphi(0))*(a22*dble(dvphi(1:ndim))-a12*aimag(dvphi(1:ndim))) &
               +aimag(dvphi(0))*(a11*aimag(dvphi(1:ndim))-a12*dble(dvphi(1:ndim)))
          vpen(1:ndim,npen)=-vpen(1:ndim,npen)/det*dspa
          
          !! debug
          !! write(*,'(A,7I5,6E10.2)') '#YY',MpiRank,ixyz(1:ndim),npen,ipla,k,du,dv,vpen(1:ndim,npen),sum(vpen(1:ndim,npen)**2)
          !! debug
          
       end do
    end do
    
    if(npen/=0) then
       vel(1:ndim)=0
       do i=1,npen
          vel(1:ndim)=vel(1:ndim)+vpen(1:ndim,npen)
       end do
       vel(1:ndim)=vel(1:ndim)/npen
    else
       !! debug
       vel(1:ndim)=1.d-10
       !!call MpiStop('error: something wrong in velocity_string_YY')
       !! debug
    end if
    
  end subroutine velocity_string_YY
  
!!$  function string_length(pos,dir)
!!$    double precision,intent(in) :: pos(ndim),dir(ndim)
!!$    double precision sect(ndim,2*ndim),w,string_length,f1,f2,g1,g2
!!$    integer iu,iv,iw,i,j,k,m
!!$    
!!$    iu=1
!!$    iv=2
!!$    iw=3
!!$    k=0
!!$    do i=1,ndim
!!$       do j=0,1
!!$          w=j-0.5d0
!!$          w=w-pos(iw)
!!$          f1=(0.5d0-pos(iu))*dir(iw)-w*dir(iu)
!!$          f2=(-0.5d0-pos(iu))*dir(iw)-w*dir(iu)
!!$          g1=(0.5d0-pos(iv))*dir(iw)-w*dir(iv)
!!$          g2=(-0.5d0-pos(iv))*dir(iw)-w*dir(iv)
!!$          if(f1*f2<0.and.g1*g2<0) then
!!$             k=k+1
!!$             sect(iu,k)=(f1+f2)/(f1-f2)/2
!!$             sect(iv,k)=(g1+g2)/(g1-g2)/2
!!$             sect(iw,k)=w
!!$          end if
!!$       end do
!!$       ! rotation
!!$       m=iu
!!$       iu=iv
!!$       iv=iw
!!$       iw=m
!!$    end do
!!$    
!!$    if(k==0) then
!!$       string_length=0.d0
!!$    else if(k==2) then
!!$       string_length=sqrt(sum((sect(1:ndim,1)-sect(1:ndim,2))**2))
!!$    else
!!$       do i=1,k
!!$          write(*,'(I,3F)') MpiRank,sect(1:3,i)
!!$       end do
!!$       call MpiStop('error: something is wrong in string_length')
!!$    end if
!!$    
!!$  end function string_length
  function string_length(pos,dir)
    double precision,intent(in) :: pos(ndim),dir(ndim)
    double precision w,string_length,t,tp,tn
    double precision, parameter :: atmax=0.867d0, & ! slightly larger than sqrt(3)/2
         atmax0=0.86602540378378444d0 ! sqrt(3)/2
    integer i,j
    
    tp=atmax
    tn=-atmax
    do i=1,ndim
       do j=0,1
          w=j-0.5d0
          w=w-pos(i)
          if(abs(w)>atmax*abs(dir(i))) cycle
          t=w/dir(i)
          if(t>0) then
             if(t<tp) tp=t
          else
             if(t>tn) tn=t
          end if
       end do
    end do
    
    if(tp<atmax0.and.tn>-atmax0) then
       string_length=tp-tn
    else
       string_length=0d0
    end if
    
  end function string_length
    
  subroutine position_wall(ph,pos,dir,gof,wa)
    use potential,only : pi,twopi
    integer,parameter :: n0=2**ndim
    double precision,intent(in) :: ph(n0)
    double precision,intent(out) :: pos(ndim),dir(ndim),gof,wa
    ! pos is point on wall closest to the origin of the cell
    ! dir is the direction normal to the wall surface
    integer,parameter :: n=12
    integer ix,iy,iz,i,j,m,imu,iphi,num
    double precision chi2,chi2min,pp(ndim,n),u,v,w,dir0(ndim), &
         mu,phi,d,mut,phit,wt,disp(ndim),th(n),uv(2,n),dir1(ndim),dir2(ndim)
    integer,parameter :: num_mu=20,num_phi=20
    
    m=0
    do iz=0,1
       do iy=0,1
          ix=0
          i=(iz*2+iy)*2+ix+1
          u=sin(ph(i))
          ix=1
          i=(iz*2+iy)*2+ix+1
          v=sin(ph(i))
          if(u*v<0d0) then
             m=m+1
             pp(1,m)=(u+v)/(u-v)/2
             pp(2,m)=iy-0.5d0
             pp(3,m)=iz-0.5d0
          end if
       end do
    end do
    do iy=0,1
       do ix=0,1
          iz=0
          i=(iz*2+iy)*2+ix+1
          u=sin(ph(i))
          iz=1
          i=(iz*2+iy)*2+ix+1
          v=sin(ph(i))
          if(u*v<0d0) then
             m=m+1
             pp(1,m)=ix-0.5d0
             pp(2,m)=iy-0.5d0
             pp(3,m)=(u+v)/(u-v)/2
          end if
       end do
    end do
    do iz=0,1
       do ix=0,1
          iy=0
          i=(iz*2+iy)*2+ix+1
          u=sin(ph(i))
          iy=1
          i=(iz*2+iy)*2+ix+1
          v=sin(ph(i))
          if(u*v<0d0) then
             m=m+1
             pp(1,m)=ix-0.5d0
             pp(2,m)=(u+v)/(u-v)/2
             pp(3,m)=iz-0.5d0
          end if
       end do
    end do
    
    if(m<3) then
       write(*,'(A,I)') ' error: # number of segments with opposite signs of ImPhi on edges= ',m
       do ix=0,1
          do iy=0,1
             do iz=0,1
                i=(iz*2+iy)*2+ix+1
                write(*,'(3I,2E)') ix,iy,iz,ph(i),sin(ph(i))
             end do
          end do
       end do
       call MpiStop('error: something wrong in position_wall')
    end if
    
    chi2min=6d0
    do imu=1,num_mu
       mu=(-num_mu+2*imu-1.d0)/num_mu
       do iphi=1,num_phi
          phi=pi*(2*iphi-1.d0-num_phi)/num_phi
          
          dir0(1)=sqrt(1-mu*mu)*cos(phi)
          dir0(2)=sqrt(1-mu*mu)*sin(phi)
          dir0(3)=mu
          
          w=0d0
          chi2=0d0
          do i=1,m
             d=sum(pp(1:ndim,i)*dir0(1:ndim))
             w=w+d
             chi2=chi2+d*d
          end do
          if(w<0d0) then
             w=0d0
          else
             w=w/m
          end if
          chi2=chi2-m*w*w
          
          if(chi2<chi2min) then
             chi2min=chi2
             mut=mu
             phit=phi
             wt=w
          end if
       end do
    end do
    
    ! determines the direction
    dir0(1)=sqrt(1-mut*mut)*cos(phit)
    dir0(2)=sqrt(1-mut*mut)*sin(phit)
    dir0(3)=mut
    num=0
    do iz=0,1
       disp(3)=iz-0.5d0
       do iy=0,1
          disp(2)=iy-0.5d0
          do ix=0,1
             disp(1)=ix-0.5d0
             d=sum(disp(1:ndim)*dir0(1:ndim))-wt
             i=(iz*2+iy)*2+ix+1
             if(d*cos(ph(i))>0d0) num=num+1 ! direction is correct
          end do
       end do
    end do
    if(num<=4) then ! flip the direction
       dir0(1:ndim)=-dir0(1:ndim)
       wt=-wt
    end if
    
    dir(1:ndim)=dir0(1:ndim)
    pos(1:ndim)=wt*dir0(1:ndim)
    
    ! goodness of fit
    gof=exp(-chi2min/2)
        
    ! wall area
    disp(1:ndim)=pos(1:ndim)
    dir1(1:ndim)=(/mut*cos(phit),mut*sin(phit),-sqrt(1-mut*mut)/)
    dir2(1:ndim)=(/-sin(phit),cos(phit),0d0/)
    do i=1,m
       pp(1:ndim,i)=pp(1:ndim,i)-disp(1:ndim)
       u=sum(dir1(1:ndim)*pp(1:ndim,i))
       v=sum(dir2(1:ndim)*pp(1:ndim,i))
       th(i)=atan2(v,u)
       uv(1:2,i)=(/u,v/)
    end do
        
    ! sort
    do i=m-1,1,-1
       do j=1,i
          if(th(j)>th(j+1)) then
             w=th(j)
             th(j)=th(j+1)
             th(j+1)=w
             u=uv(1,j)
             v=uv(2,j)
             uv(1:2,j)=uv(1:2,j+1)
             uv(1:2,j+1)=(/u,v/)
          end if
       end do
    end do
        
    wa=0d0
    w=th(m)-twopi
    u=sqrt(uv(1,m)**2+uv(2,m)**2)
    do i=1,m
       v=sqrt(uv(1,i)**2+uv(2,i)**2)
       w=th(i)-w
       wa=wa+0.5d0*u*v*sin(w)
       u=v
       w=th(i)
    end do
    
  end subroutine position_wall
  
  subroutine grouping_strings(nstot,arr_string,wid_st,hor_size,num_strings)
    integer,intent(in) :: nstot
    type(string_list) :: arr_string(nstot)
    double precision,intent(in) :: wid_st,hor_size
    integer,intent(out) :: num_strings
    
    integer num
    double precision th_dist
    double precision,parameter :: pow_st=2.d0,def_th_dist_st=2.d0 ! this might need to be tuned ...
    
    th_dist=(wid_st**pow_st*hor_size)**(1/(pow_st+1)) ! weighted geometric mean of string width and string separation
    th_dist=max(def_th_dist_st,th_dist)
    if(MpiRank==0) then
       write(*,'(" size of string width in lattice units:",F)') wid_st
       write(*,'(" size of string separation in lattice units:",F)') hor_size
       write(*,'(" threshold for string grouping in lattice units:",F)') th_dist
    end if
    
    ! get local id in each MPI process
    call grouping_strings_local(nstot,arr_string,th_dist,num)
    
    ! get connections between different MPI processes
    call grouping_strings_global(nstot,arr_string,th_dist,num)
    num_strings=num
    
  end subroutine grouping_strings
  
  subroutine grouping_strings_global(nstot,arr_string,th_dist,global_num)
    integer,intent(in) :: nstot
    type(string_list) :: arr_string(nstot)
    double precision,intent(in) :: th_dist
    integer,intent(inout) :: global_num ! note: local_num at the input
    
    integer MpiRankA,MpiRankD,num_out,num_in,info,i,j,k,n,ng,ngp,loc_id_offset, &
         isend,irecv,ssend(MPI_STATUS_SIZE),srecv(MPI_STATUS_SIZE)
    double precision dist
    double precision,allocatable :: pos_out(:,:),pos_in(:,:)
    integer,allocatable :: id_out(:),id_in(:),arr_ng(:),global_ids(:)
    logical,allocatable :: conn(:),conn0(:)
    
    if(MpiRank==0) write(*,'(A)') '  grouping strings globally...'
    
    MpiRankA=modulo(MpiRank+1,MpiSize)
    MpiRankD=modulo(MpiRank-1,MpiSize)
    
    allocate(arr_ng(0:MpiSize-1),stat=info)
    if(info/=0) call MpiStop('error: cannot allocate arr_ng')
#ifdef MPI
    call MPI_ALLGATHER(global_num,1,MPI_INTEGER,arr_ng,1,MPI_INTEGER,MPI_COMM_WORLD,info)
#else
    arr_ng(0)=global_num
#endif
    ng=sum(arr_ng(0:MpiSize-1))
    ngp=ng*(ng-1)/2
    loc_id_offset=sum(arr_ng(0:MpiRank-1))
    deallocate(arr_ng,stat=info)
    
    ! get connections between different MPI processes
    num_out=0
    do i=1,nstot
       if(arr_string(i)%pos(ndim)>=local_0_offset+local_0-th_dist) num_out=num_out+1
    end do
    
    call MPI_ISEND(num_out,1,MPI_INTEGER,MpiRankA,0,MPI_COMM_WORLD,isend,info)
    call MPI_IRECV(num_in,1,MPI_INTEGER,MpiRankD,0,MPI_COMM_WORLD,irecv,info)
    call MPI_WAIT(isend,ssend,info)
    call MPI_WAIT(irecv,srecv,info)
    
    if(num_out/=0) then
       allocate(pos_out(ndim,num_out),id_out(num_out),stat=info)
       if(info/=0) call MpiStop('error: cannot allocate pos_out or id_out')
    end if
    if(num_in/=0) then
       allocate(pos_in(ndim,num_in),id_in(num_in),stat=info)
       if(info/=0) call MpiStop('error: cannot allocate pos_in, id_in')
    end if
    
    if(num_out/=0) then
       n=0
       do i=1,nstot
          if(.not.arr_string(i)%pos(ndim)>=local_0_offset+local_0-th_dist) cycle
          
          n=n+1
          if(n>num_out) call MpiStop('error: num_out is odd')
          pos_out(1:ndim,n)=arr_string(i)%pos(1:ndim)
          id_out(n)=arr_string(i)%loc_id+loc_id_offset
       end do
       if(n/=num_out) call MpiStop('error: something wrong with num_out')
       call MPI_ISEND(pos_out,ndim*num_out,MPI_DOUBLE_PRECISION,MpiRankA,0,MPI_COMM_WORLD,isend,info)
    end if
    if(num_in/=0) call MPI_IRECV(pos_in,ndim*num_in,MPI_DOUBLE_PRECISION,MpiRankD,0,MPI_COMM_WORLD,irecv,info)
    if(num_out/=0) call MPI_WAIT(isend,ssend,info)
    if(num_in/=0) call MPI_WAIT(irecv,srecv,info)
    
    if(num_out/=0) call MPI_ISEND(id_out,num_out,MPI_INTEGER,MpiRankA,0,MPI_COMM_WORLD,isend,info)
    if(num_in/=0) call MPI_IRECV(id_in,num_in,MPI_INTEGER,MpiRankD,0,MPI_COMM_WORLD,irecv,info)
    if(num_out/=0) call MPI_WAIT(isend,ssend,info)
    if(num_in/=0) call MPI_WAIT(irecv,srecv,info)
    
    if(num_out/=0) deallocate(pos_out,id_out,stat=info)
    write(*,'(A,I,A)') '   MpiRank ',MpiRank,' preparation of connection: done'
    
    allocate(conn0(ngp),stat=info)
    if(info/=0) call MpiStop('error: cannot allocate conn0')
    conn0(1:ngp)=.false.
    if(num_in/=0) then
       !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,k,n,dist) SCHEDULE(DYNAMIC)
       do i=1,nstot
          if(.not.arr_string(i)%pos(ndim)<=local_0_offset+1+th_dist) cycle
          
          n=arr_string(i)%loc_id+loc_id_offset
          
          do j=1,num_in
             
             if(MpiRankD>MpiRank) then
                k=(id_in(j)-1)*(id_in(j)-2)/2+n
             else
                k=(n-1)*(n-2)/2+id_in(j)
             end if
             
             if(conn0(k)) cycle ! already connected
             
             call lat_dist_periodic(arr_string(i)%pos(1:ndim),pos_in(1:ndim,j),dist)
             if(dist<=th_dist) then
                !$OMP CRITICAL
                conn0(k)=.true.
                !$OMP END CRITICAL
             end if
             
          end do
          
       end do
       !$OMP END PARALLEL DO
       deallocate(pos_in,id_in,stat=info)
    end if
    write(*,'(A,I,A)') '   MpiRank ',MpiRank,' connecting local groups: done'
    
    allocate(conn(ngp),stat=info)
    if(info/=0) call MpiStop('error: cannot allocate conn')
#ifdef MPI
    call MPI_ALLREDUCE(conn0,conn,ngp,MPI_LOGICAL,MPI_LOR,MPI_COMM_WORLD,info)
#else
    conn(1:ngp)=conn0(1:ngp)
#endif
    deallocate(conn0,stat=info)
    
    allocate(global_ids(ngp),stat=info)
    if(info/=0) call MpiStop('error: cannot allocate global_ids')
    global_ids(1:ngp)=0
    n=0
    do
       
       do i=1,ng
          if(global_ids(i)/=0) cycle
          n=n+1
          global_ids(i)=-n ! new groups have negative index
          exit
       end do
       
       if(i>ng) exit ! no more group
       
       do while(i<=ng)
          k=i
          i=k+1
          
          if(global_ids(k)/=-n) cycle
          
          do j=1,ng
             if(global_ids(j)/=0) cycle
             
             if(k==j) call MpiStop('error: something is wrong in global_ids')
             if(conn((max(k,j)-1)*(max(k,j)-2)/2+min(k,j))) then 
                global_ids(j)=-n
                i=1
             end if
          end do
          
          global_ids(k)=n
       end do
       
       ! check
       do i=1,ng
          if(global_ids(i)/=0) cycle
          do j=1,ng
             if(global_ids(j)/=n) cycle
             if(conn((max(i,j)-1)*(max(i,j)-2)/2+min(i,j))) call MpiStop('error: found non-grouped groups')
          end do
       end do
       
    end do
    global_num=n
    write(*,'(A,I,A)') '   MpiRank ',MpiRank,' assignment of global ids: done'
    
    deallocate(conn,stat=info)
    
    do i=1,nstot
       arr_string(i)%glo_id=global_ids(arr_string(i)%loc_id+loc_id_offset)
    end do
    
    deallocate(global_ids,stat=info)
    
    call MPI_BARRIER(MPI_COMM_WORLD,info)
    if(MpiRank==0) write(*,'(A)') '  ...done'
    
  end subroutine grouping_strings_global
  
  subroutine grouping_strings_local(nstot,arr_string,th_dist,local_num)
    ! get local id in each MPI process
    integer,intent(in) :: nstot
    type(string_list) :: arr_string(nstot)
    double precision,intent(in) :: th_dist
    integer,intent(out) :: local_num
    
    double precision y(nstot),ya,yd,dist
    integer idx(nstot),info,it,n,i,j,k,it2,ng,ngp
    integer,allocatable :: min_ns(:),max_ns(:),arr_ng(:),local_ids(:)
    type(string_list) tmp_string(nstot)
    logical,allocatable :: conn(:)
    
    call MPI_BARRIER(MPI_COMM_WORLD,info)
    if(MpiRank==0) write(*,'(A)') '  grouping strings locally...'
    
    ! sort by y-coodinate
    do i=1,nstot
       y(i)=arr_string(i)%pos(ndim-1)
       idx(i)=i
    end do
    call dlasrt2('I',nstot,y,idx,info) ! quicksort; parallelized?
    do i=1,nstot
       tmp_string(i)%pos(1:ndim)=arr_string(idx(i))%pos(1:ndim)
       tmp_string(i)%loc_id=0
    end do
    write(*,'(A,I,A)') '   MpiRank ',MpiRank,' quicksort in y: done'
    
    ! domain decomposition
    allocate(min_ns(num_threads),max_ns(num_threads),arr_ng(num_threads),stat=info)
    if(info/=0) call MpiStop('error: cannot allocate min_ns, max_ns or arr_ng')
    i=nstot/num_threads
    j=nstot-i*num_threads
    min_ns(1)=1
    do it=1,num_threads
       max_ns(it)=min_ns(it)+i
       if(it>j) max_ns(it)=max_ns(it)-1
       if(it==num_threads) exit
       min_ns(it+1)=max_ns(it)+1
    end do
    if(max_ns(num_threads)/=nstot) call MpiStop('error: something is wrong in computing min_ns & max_ns') !check
    write(*,'(A,I,A)') '   MpiRank ',MpiRank,' domain-decomposition in y: done'
    
    !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(it,n,i,j,k,dist)
#ifdef _OPENMP
    it=OMP_GET_THREAD_NUM()+1
#else
    it=1
#endif
    n=0
    do
       
       do i=min_ns(it),max_ns(it)
          if(tmp_string(i)%loc_id/=0) cycle
          n=n+1
          tmp_string(i)%loc_id=-n ! new points have negative loc_id
          exit
       end do
       
       if(i>max_ns(it)) exit ! no more new group
       
       do while(i<=max_ns(it))
          k=i
          i=k+1
          
          if(tmp_string(k)%loc_id/=-n) cycle
          
          do j=min_ns(it),max_ns(it)
             if(tmp_string(j)%loc_id/=0) cycle
             
             call lat_dist_periodic(tmp_string(k)%pos,tmp_string(j)%pos,dist)
             if(dist<=th_dist) then
                tmp_string(j)%loc_id=-n ! new point
                i=min_ns(it)
             end if
          end do
          
          tmp_string(k)%loc_id=n ! done; old points have positive loc_id
       end do
       
       ! check
       do i=min_ns(it),max_ns(it)
          if(tmp_string(i)%loc_id/=0) cycle
          do j=min_ns(it),max_ns(it)
             if(tmp_string(j)%loc_id/=n) cycle
             call lat_dist_periodic(tmp_string(i)%pos,tmp_string(j)%pos,dist)
             if(dist<=th_dist) call MpiStop('error: found non-grouped points')
          end do
       end do
       
    end do
    arr_ng(it)=n
    !$OMP END PARALLEL
    write(*,'(A,I,A)') '   MpiRank ',MpiRank,' pre-grouping: done'
    
    do it=1,num_threads
       n=sum(arr_ng(1:it-1)) !offset
       do i=min_ns(it),max_ns(it)
          tmp_string(i)%loc_id=tmp_string(i)%loc_id+n
       end do
    end do
    
    ng=sum(arr_ng)
    deallocate(arr_ng,stat=info)
    
    ngp=ng*(ng-1)/2
    allocate(conn(ngp),stat=info)
    if(info/=0) call MpiStop('error: cannot allocate conn')
    conn(1:ngp)=.false.
    
    !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(it,it2,i,j,ya,yd,dist,k)
#ifdef _OPENMP
    it=OMP_GET_THREAD_NUM()+1
#else
    it=1
#endif
    it2=modulo(it,num_threads)+1
    
    ya=tmp_string(max_ns(it))%pos(ndim-1)+th_dist
    yd=tmp_string(min_ns(it2))%pos(ndim-1)-th_dist
    do i=max_ns(it),min_ns(it),-1
       
       if(.not.tmp_string(i)%pos(ndim-1)>yd) exit
       
       do j=min_ns(it2),max_ns(it2)
          
          if(.not.tmp_string(j)%pos(ndim-1)<ya) exit

          if(it2>it) then
             k=(tmp_string(j)%loc_id-1)*(tmp_string(j)%loc_id-2)/2+tmp_string(i)%loc_id
          else
             k=(tmp_string(i)%loc_id-1)*(tmp_string(i)%loc_id-2)/2+tmp_string(j)%loc_id
          end if
          if(conn(k)) cycle
          
          call lat_dist_periodic(tmp_string(i)%pos,tmp_string(j)%pos,dist)
          if(dist<=th_dist) conn(k)=.true.
       end do
    end do
    !$OMP END PARALLEL
    deallocate(min_ns,max_ns,stat=info)
    write(*,'(A,I,A)') '   MpiRank ',MpiRank,' connection of pregroups: done'
    
    allocate(local_ids(ng),stat=info)
    if(info/=0) call MpiStop('error: cannot allocate local_ids')
    local_ids(1:ng)=0
    
    n=0
    do 
       
       do i=1,ng
          if(local_ids(i)/=0) cycle
          n=n+1
          local_ids(i)=-n ! new group has negative index
          exit
       end do
       
       if(i>ng) exit ! no more new group
       
       do while(i<=ng)
          k=i
          i=k+1
          
          if(local_ids(k)/=-n) cycle
          
          do j=1,ng
             if(local_ids(j)/=0) cycle
             
             if(k==j) call MpiStop('error: something is wrong in local_ids')
             if(conn((max(k,j)-1)*(max(k,j)-2)/2+min(k,j))) then 
                local_ids(j)=-n
                i=1
             end if
          end do
          
          local_ids(k)=n
       end do
       
       ! check
       do i=1,ng
          if(local_ids(i)/=0) cycle
          do j=1,ng
             if(local_ids(j)/=n) cycle
             if(conn((max(i,j)-1)*(max(i,j)-2)/2+min(i,j))) call MpiStop('error: found non-grouped pregroups')
          end do
       end do

    end do
    local_num=n
    write(*,'(A,I,A)') '   MpiRank ',MpiRank,' assignment of local ids: done'
    
    deallocate(conn,stat=info)
    
    do i=1,nstot
       arr_string(idx(i))%loc_id=local_ids(tmp_string(i)%loc_id)
    end do
    
    deallocate(local_ids,stat=info)

    call MPI_BARRIER(MPI_COMM_WORLD,info)
    if(MpiRank==0) write(*,'(A)') '  ...done'
    
  end subroutine grouping_strings_local
  
  subroutine lat_dist_periodic(p1,p2,dist)
    double precision,intent(in) :: p1(ndim),p2(ndim)
    double precision,intent(out) :: dist
    double precision dp(ndim)
    
    dp(1:ndim)=abs(p1(1:ndim)-p2(1:ndim))
    dist=sqrt(sum(min(dp(1:ndim),ngrid-dp(1:ndim))**2))
  end subroutine lat_dist_periodic
  
  subroutine add_string_length(nstot,arr_string,num_strings,hor_size,sl,svel)
    use potential,only : twopi
    integer,intent(in) :: nstot
    type(string_list) arr_string(nstot)
    integer,intent(in) :: num_strings
    double precision,intent(in) :: hor_size
    
    double precision,intent(out) :: sl(3),svel(3)
    double precision,allocatable :: dsl(:),ddsl(:)
    integer i,info,n,nn
    double precision ddsvel(3),dsvel(3),v2
    
    sl(1:3)=0.d0
    
    if(.not.group_strings) then
       
       allocate(dsl(1),ddsl(1),stat=info)
       ddsl(1)=0.d0
       
       do i=1,nstot
          ddsl(1)=ddsl(1)+arr_string(i)%sl
       end do
       
       call MPI_ALLREDUCE(ddsl,dsl,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
       sl(1)=dsl(1)
       deallocate(dsl,ddsl,stat=info)
       
    else
       
       if(.not.num_strings==0) then
       
          allocate(dsl(num_strings),ddsl(num_strings),stat=info)
          ddsl(1:num_strings)=0.d0
          do i=1,nstot
             ddsl(arr_string(i)%glo_id)=ddsl(arr_string(i)%glo_id)+arr_string(i)%sl
          end do
          
          call MPI_ALLREDUCE(ddsl,dsl,num_strings,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
          sl(1:3)=0.d0
          do i=1,num_strings
             sl(1)=sl(1)+dsl(i)
             if(dsl(i)>2*hor_size) sl(2)=sl(2)+dsl(i)
             if(dsl(i)>twopi*hor_size) sl(3)=sl(3)+dsl(i)
          end do
          deallocate(dsl,ddsl,stat=info)
          
       end if
       
    end if
    
    if(.not.group_strings) sl(2:3)=sl(1)
    
    n=0
    ddsvel(1:3)=0.d0
    do i=1,nstot
       v2=sum(arr_string(i)%vel(1:ndim)**2)
!======================================================================================================================================================================================
       if(v2>=1.d0) cycle ! omit these
       !!$if(v2>1) v2=1.d0
!======================================================================================================================================================================================
       
       n=n+1
       ddsvel(1)=ddsvel(1)+sqrt(v2)
       ddsvel(2)=ddsvel(2)+v2
!======================================================================================================================================================================================
       ddsvel(3)=ddsvel(3)+1/sqrt(1-v2)
!======================================================================================================================================================================================
    end do
    
    call MPI_ALLREDUCE(ddsvel,dsvel,3,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
    call MPI_ALLREDUCE(n,nn,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
    
    if(nn==0) then
       if(MpiRank==0) write(*,*) ' velocity cannot be computed because there are no strings'
       svel(1:3)=0.d0
    else
       svel(1:3)=dsvel(1:3)/nn
    end if
    
  end subroutine add_string_length
  
  subroutine grouping_walls(nwtot,arr_wall,wid_dw,hor_size,num_walls)
    integer,intent(in) :: nwtot
    type(wall_list) :: arr_wall(nwtot)
    double precision,intent(in) :: wid_dw,hor_size
    integer,intent(out) :: num_walls
    
    integer num
    double precision th_dist
    double precision,parameter :: pow_dw=2.d0,def_th_dist_dw=2.d0 ! this might need to be tuned ...
    
    th_dist=(wid_dw**pow_dw*hor_size)**(1/(pow_dw+1)) ! weighted geometric mean of wall width and wall separation
    th_dist=max(def_th_dist_dw,th_dist)
    if(MpiRank==0) then
       write(*,'(" size of wall width in lattice units:",F)') wid_dw
       write(*,'(" size of wall separation in lattice units:",F)') hor_size
       write(*,'(" threshold for wall grouping in lattice units:",F)') th_dist
    end if
    
    ! get local id in each MPI process
    call grouping_walls_local(nwtot,arr_wall,th_dist,num)
    
    ! get connections between different MPI processes
    call grouping_walls_global(nwtot,arr_wall,th_dist,num)
    num_walls=num
    
  end subroutine grouping_walls
  
  subroutine grouping_walls_global(nwtot,arr_wall,th_dist,global_num)
    integer,intent(in) :: nwtot
    type(wall_list) :: arr_wall(nwtot)
    double precision,intent(in) :: th_dist
    integer,intent(inout) :: global_num ! note: local_num at the input
    
    integer MpiRankA,MpiRankD,num_out,num_in,info,i,j,k,n,ng,ngp,loc_id_offset, &
         isend,irecv,ssend(MPI_STATUS_SIZE),srecv(MPI_STATUS_SIZE)
    double precision dist
    double precision,allocatable :: pos_out(:,:),pos_in(:,:)
    integer,allocatable :: id_out(:),id_in(:),arr_ng(:),global_ids(:)
    logical,allocatable :: conn(:),conn0(:)
    
    if(MpiRank==0) write(*,'(A)') '  grouping walls globally...'
    
    MpiRankA=modulo(MpiRank+1,MpiSize)
    MpiRankD=modulo(MpiRank-1,MpiSize)
    
    allocate(arr_ng(0:MpiSize-1),stat=info)
    if(info/=0) call MpiStop('error: cannot allocate arr_ng')
#ifdef MPI
    call MPI_ALLGATHER(global_num,1,MPI_INTEGER,arr_ng,1,MPI_INTEGER,MPI_COMM_WORLD,info)
#else
    arr_ng(0)=global_num
#endif
    ng=sum(arr_ng(0:MpiSize-1))
    ngp=ng*(ng-1)/2
    loc_id_offset=sum(arr_ng(0:MpiRank-1))
    deallocate(arr_ng,stat=info)
    
    ! get connections between different MPI processes
    num_out=0
    do i=1,nwtot
       if(arr_wall(i)%pos(ndim)>=local_0_offset+local_0-th_dist) num_out=num_out+1
    end do
    
    call MPI_ISEND(num_out,1,MPI_INTEGER,MpiRankA,0,MPI_COMM_WORLD,isend,info)
    call MPI_IRECV(num_in,1,MPI_INTEGER,MpiRankD,0,MPI_COMM_WORLD,irecv,info)
    call MPI_WAIT(isend,ssend,info)
    call MPI_WAIT(irecv,srecv,info)
    
    if(num_out/=0) then
       allocate(pos_out(ndim,num_out),id_out(num_out),stat=info)
       if(info/=0) call MpiStop('error: cannot allocate pos_out or id_out')
    end if
    if(num_in/=0) then
       allocate(pos_in(ndim,num_in),id_in(num_in),stat=info)
       if(info/=0) call MpiStop('error: cannot allocate pos_in, id_in')
    end if
    
    if(num_out/=0) then
       n=0
       do i=1,nwtot
          if(.not.arr_wall(i)%pos(ndim)>=local_0_offset+local_0-th_dist) cycle
          
          n=n+1
          if(n>num_out) call MpiStop('error: num_out is odd')
          pos_out(1:ndim,n)=arr_wall(i)%pos(1:ndim)
          id_out(n)=arr_wall(i)%loc_id+loc_id_offset
       end do
       if(n/=num_out) call MpiStop('error: something wrong with num_out')
       call MPI_ISEND(pos_out,ndim*num_out,MPI_DOUBLE_PRECISION,MpiRankA,0,MPI_COMM_WORLD,isend,info)
    end if
    if(num_in/=0) call MPI_IRECV(pos_in,ndim*num_in,MPI_DOUBLE_PRECISION,MpiRankD,0,MPI_COMM_WORLD,irecv,info)
    if(num_out/=0) call MPI_WAIT(isend,ssend,info)
    if(num_in/=0) call MPI_WAIT(irecv,srecv,info)
    
    if(num_out/=0) call MPI_ISEND(id_out,num_out,MPI_INTEGER,MpiRankA,0,MPI_COMM_WORLD,isend,info)
    if(num_in/=0) call MPI_IRECV(id_in,num_in,MPI_INTEGER,MpiRankD,0,MPI_COMM_WORLD,irecv,info)
    if(num_out/=0) call MPI_WAIT(isend,ssend,info)
    if(num_in/=0) call MPI_WAIT(irecv,srecv,info)
    
    if(num_out/=0) deallocate(pos_out,id_out,stat=info)
    write(*,'(A,I,A)') '   MpiRank ',MpiRank,' preparation for connection: done'
    
    allocate(conn0(ngp),stat=info)
    if(info/=0) call MpiStop('error: cannot allocate conn0')
    conn0(1:ngp)=.false.
    if(num_in/=0) then
       !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,k,n,dist) SCHEDULE(DYNAMIC)
       do i=1,nwtot
          if(.not.arr_wall(i)%pos(ndim)<=local_0_offset+1+th_dist) cycle
          
          n=arr_wall(i)%loc_id+loc_id_offset
          
          do j=1,num_in
             
             if(MpiRankD>MpiRank) then
                k=(id_in(j)-1)*(id_in(j)-2)/2+n
             else
                k=(n-1)*(n-2)/2+id_in(j)
             end if
             
             if(conn0(k)) cycle ! already connected
             
             call lat_dist_periodic(arr_wall(i)%pos(1:ndim),pos_in(1:ndim,j),dist)
             if(dist<=th_dist) then
                !$OMP CRITICAL
                conn0(k)=.true.
                !$OMP END CRITICAL
             end if
             
          end do
          
       end do
       !$OMP END PARALLEL DO
       deallocate(pos_in,id_in,stat=info)
    end if
    write(*,'(A,I,A)') '   MpiRank ',MpiRank,' connecting local groups: done'
    
    allocate(conn(ngp),stat=info)
    if(info/=0) call MpiStop('error: cannot allocate conn')
#ifdef MPI
    call MPI_ALLREDUCE(conn0,conn,ngp,MPI_LOGICAL,MPI_LOR,MPI_COMM_WORLD,info)
#else
    conn(1:ngp)=conn0(1:ngp)
#endif
    deallocate(conn0,stat=info)
    
    allocate(global_ids(ngp),stat=info)
    if(info/=0) call MpiStop('error: cannot allocate global_ids')
    global_ids(1:ngp)=0
    n=0
    do
       
       do i=1,ng
          if(global_ids(i)/=0) cycle
          n=n+1
          global_ids(i)=-n ! new groups have negative index
          exit
       end do
       
       if(i>ng) exit ! no more group
       
       do while(i<=ng)
          k=i
          i=k+1
          
          if(global_ids(k)/=-n) cycle
          
          do j=1,ng
             if(global_ids(j)/=0) cycle
             
             if(k==j) call MpiStop('error: something is wrong in global_ids')
             if(conn((max(k,j)-1)*(max(k,j)-2)/2+min(k,j))) then 
                global_ids(j)=-n
                i=1
             end if
          end do
          
          global_ids(k)=n
       end do
       
       ! check
       do i=1,ng
          if(global_ids(i)/=0) cycle
          do j=1,ng
             if(global_ids(j)/=n) cycle
             if(conn((max(i,j)-1)*(max(i,j)-2)/2+min(i,j))) call MpiStop('error: found non-grouped groups')
          end do
       end do
       
    end do
    global_num=n
    write(*,'(A,I,A)') '   MpiRank ',MpiRank,' assignment of global ids: done'
    
    deallocate(conn,stat=info)
    
    do i=1,nwtot
       arr_wall(i)%glo_id=global_ids(arr_wall(i)%loc_id+loc_id_offset)
    end do
    
    deallocate(global_ids,stat=info)
    
    call MPI_BARRIER(MPI_COMM_WORLD,info)
    if(MpiRank==0) write(*,'(A)') '  ...done'
    
  end subroutine grouping_walls_global
  
  subroutine grouping_walls_local(nwtot,arr_wall,th_dist,local_num)
    ! get local id in each MPI process
    integer,intent(in) :: nwtot
    type(wall_list) :: arr_wall(nwtot)
    double precision,intent(in) :: th_dist
    integer,intent(out) :: local_num
    
    double precision y(nwtot),ya,yd,dist
    integer idx(nwtot),info,it,n,i,j,k,it2,ng,ngp
    integer,allocatable :: min_nw(:),max_nw(:),arr_ng(:),local_ids(:)
    type(wall_list) tmp_wall(nwtot)
    logical,allocatable :: conn(:)
    
    call MPI_BARRIER(MPI_COMM_WORLD,info)
    if(MpiRank==0) write(*,'(A)') '  grouping walls locally...'
    
    ! sort by y-coodinate
    do i=1,nwtot
       y(i)=arr_wall(i)%pos(ndim-1)
       idx(i)=i
    end do
    call dlasrt2('I',nwtot,y,idx,info) ! quicksort; parallelized?
    do i=1,nwtot
       tmp_wall(i)%pos(1:ndim)=arr_wall(idx(i))%pos(1:ndim)
       tmp_wall(i)%loc_id=0
    end do
    write(*,'(A,I,A)') '   MpiRank ',MpiRank,' quickosrt in y: done'

    ! domain decomposition
    allocate(min_nw(num_threads),max_nw(num_threads),arr_ng(num_threads),stat=info)
    if(info/=0) call MpiStop('error: cannot allocate min_nw, max_nw or arr_ng')
    i=nwtot/num_threads
    j=nwtot-i*num_threads
    min_nw(1)=1
    do it=1,num_threads
       max_nw(it)=min_nw(it)+i
       if(it>j) max_nw(it)=max_nw(it)-1
       if(it==num_threads) exit
       min_nw(it+1)=max_nw(it)+1
    end do
    if(max_nw(num_threads)/=nwtot) call MpiStop('error: something is wrong in computing min_nw & max_nw') !check
    write(*,'(A,I,A)') '   MpiRank ',MpiRank,' domain-decomposition in y: done'
    
    !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(it,n,i,j,k,dist)
#ifdef _OPENMP
    it=OMP_GET_THREAD_NUM()+1
#else
    it=1
#endif
    n=0
    do
       
       do i=min_nw(it),max_nw(it)
          if(tmp_wall(i)%loc_id/=0) cycle
          n=n+1
          tmp_wall(i)%loc_id=-n ! new points have negative loc_id
          exit
       end do
       
       if(i>max_nw(it)) exit ! no more new group

       do while(i<=max_nw(it))
          k=i
          i=k+1
          
          if(tmp_wall(k)%loc_id/=-n) cycle
          
          do j=min_nw(it),max_nw(it)
             if(tmp_wall(j)%loc_id/=0) cycle
             
             call lat_dist_periodic(tmp_wall(k)%pos,tmp_wall(j)%pos,dist)
             if(dist<=th_dist) then
                tmp_wall(j)%loc_id=-n ! new point
                i=min_nw(it)
             end if
          end do
          
          tmp_wall(k)%loc_id=n ! done; old points have positive loc_id
       end do
          
       ! check
       do i=min_nw(it),max_nw(it)
          if(tmp_wall(i)%loc_id/=0) cycle
          do j=min_nw(it),max_nw(it)
             if(tmp_wall(j)%loc_id/=n) cycle
             call lat_dist_periodic(tmp_wall(i)%pos,tmp_wall(j)%pos,dist)
             if(dist<=th_dist) call MpiStop('error: found non-grouped points')
          end do
       end do
       
    end do
    arr_ng(it)=n
    !$OMP END PARALLEL
    write(*,'(A,I,A)') '   MpiRank ',MpiRank,' pregrouping: done'
    
    do it=1,num_threads
       n=sum(arr_ng(1:it-1)) !offset
       do i=min_nw(it),max_nw(it)
          tmp_wall(i)%loc_id=tmp_wall(i)%loc_id+n
       end do
    end do
    
    ng=sum(arr_ng)
    deallocate(arr_ng,stat=info)
    
    ngp=ng*(ng-1)/2
    allocate(conn(ngp),stat=info)
    if(info/=0) call MpiStop('error: cannot allocate conn')
    conn(1:ngp)=.false.
    
    !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(it,it2,i,j,ya,yd,dist,k)
#ifdef _OPENMP
    it=OMP_GET_THREAD_NUM()+1
#else
    it=1
#endif
    it2=modulo(it,num_threads)+1
    
    ya=tmp_wall(max_nw(it))%pos(ndim-1)+th_dist
    yd=tmp_wall(min_nw(it2))%pos(ndim-1)-th_dist
    do i=max_nw(it),min_nw(it),-1
       
       if(.not.tmp_wall(i)%pos(ndim-1)>yd) exit
       
       do j=min_nw(it2),max_nw(it2)
          
          if(.not.tmp_wall(j)%pos(ndim-1)<ya) exit
          
          if(it2>it) then
             k=(tmp_wall(j)%loc_id-1)*(tmp_wall(j)%loc_id-2)/2+tmp_wall(i)%loc_id
          else
             k=(tmp_wall(i)%loc_id-1)*(tmp_wall(i)%loc_id-2)/2+tmp_wall(j)%loc_id
          end if
          if(conn(k)) cycle
          
          call lat_dist_periodic(tmp_wall(i)%pos,tmp_wall(j)%pos,dist)
          if(dist<=th_dist) conn(k)=.true.
       end do
    end do
    !$OMP END PARALLEL
    deallocate(min_nw,max_nw,stat=info)
    write(*,'(A,I,A)') '   MpiRank ',MpiRank,' connecting pregroups: done'
    
    allocate(local_ids(ng),stat=info)
    if(info/=0) call MpiStop('error: cannot allocate local_ids')
    local_ids(1:ng)=0
    
    n=0
    do 
       
       do i=1,ng
          if(local_ids(i)/=0) cycle
          n=n+1
          local_ids(i)=-n ! new group has negative index
          exit
       end do
       
       if(i>ng) exit ! no more new group
       
       do while(i<=ng)
          k=i
          i=k+1
          
          if(local_ids(k)/=-n) cycle
          
          do j=1,ng
             if(local_ids(j)/=0) cycle
             
             if(k==j) call MpiStop('error: something is wrong in local_ids')
             if(conn((max(k,j)-1)*(max(k,j)-2)/2+min(k,j))) then 
                local_ids(j)=-n
                i=1
             end if
          end do
          
          local_ids(k)=n
       end do
       
       ! check
       do i=1,ng
          if(local_ids(i)/=0) cycle
          do j=1,ng
             if(local_ids(j)/=n) cycle
             if(conn((max(i,j)-1)*(max(i,j)-2)/2+min(i,j))) call MpiStop('error: found non-grouped pregroups')
          end do
       end do

    end do
    local_num=n
    write(*,'(A,I,A)') '   MpiRank ',MpiRank,' assignment of local identities: done'
    
    deallocate(conn,stat=info)
    
    do i=1,nwtot
       arr_wall(idx(i))%loc_id=local_ids(tmp_wall(i)%loc_id)
    end do
    
    deallocate(local_ids,stat=info)
    
    call MPI_BARRIER(MPI_COMM_WORLD,info)
    if(MpiRank==0) write(*,'(A)') '  ...done'
    
  end subroutine grouping_walls_local
  
  subroutine add_wall_area(nwtot,arr_wall,num_walls,hor_size,wa)
    use potential,only : pi,fourpi
    integer,intent(in) :: nwtot
    type(wall_list) arr_wall(nwtot)
    integer,intent(in) :: num_walls
    double precision,intent(in) :: hor_size
    double precision,intent(out) :: wa(3)

    double precision,allocatable :: dwa(:),ddwa(:)
    integer i,info
    
    wa(1:3)=0.d0
    
    if(.not.group_walls) then
       
       allocate(dwa(1),ddwa(1),stat=info)
       ddwa(1)=0.d0
       do i=1,nwtot
          ddwa(1)=ddwa(1)+arr_wall(i)%wa
       end do
       
#ifdef MPI
       call MPI_ALLREDUCE(ddwa,dwa,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
#else
       dwa(1)=ddwa(1)
#endif
       wa(1)=dwa(1)
       deallocate(dwa,ddwa,stat=info)
       
    else
       
       if(.not.num_walls==0) then
          
          allocate(dwa(num_walls),ddwa(num_walls),stat=info)
          ddwa(1:num_walls)=0.d0
          do i=1,nwtot
             ddwa(arr_wall(i)%glo_id)=ddwa(arr_wall(i)%glo_id)+arr_wall(i)%wa
          end do
          
#ifdef MPI
          call MPI_ALLREDUCE(ddwa,dwa,num_walls,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
#else
          dwa(1:num_walls)=ddwa(1:num_walls)
#endif
          
          wa(1:3)=0.d0
          do i=1,num_walls
             wa(1)=wa(1)+dwa(i)
             if(dwa(i)>pi*hor_size**2) wa(2)=wa(2)+dwa(i)
             if(dwa(i)>fourpi*hor_size**2) wa(3)=wa(3)+dwa(i)
          end do
          deallocate(dwa,ddwa,stat=info)
       
       end if
       
    end if
    
    if(.not.group_walls) wa(2:3)=wa(1)
    
  end subroutine add_wall_area
  
end module defects_3d
