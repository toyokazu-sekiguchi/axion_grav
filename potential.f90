module potential
  ! module for scalar potential
  implicit none

  !=========================
  ! mathematical constants
  !=========================
  double precision,parameter :: pi=3.14159265359d0,twopi=2*pi,fourpi=4*pi

  !=========================
  ! physical parameters
  !=========================
  double precision,parameter :: Tcmb_in_K=2.726d0, &
       Mstar_in_GeV=2.435d18, & ! reduced planck mass in GeV
       GeV_in_K=1.16045250061657d13, & ! GeV in Kelvin
       GeV_in_Hz=1.5192488833521d24, & ! GeV in Hz
       BigH=2.1332d-42, & ! H0/h in GeV
       g_rel_now=3.3626d0, & ! relativistic degrees of freedom in terms of energy
       g_rel_nowS=3.9091d0 ! rel. deg. of freedom in terms of entropy
  
  double precision,parameter :: Tcmb_in_Mstar=Tcmb_in_K/GeV_in_K/Mstar_in_GeV, &
       BigH_in_Mstar=BigH/Mstar_in_GeV, &
       Mstar_in_Hz=Mstar_in_GeV*GeV_in_Hz
  
  double precision kappa
  double precision,parameter :: n_axion_IILM=6.68d0, &
       cT_axion_IILM=6.26d0,c0_axion=1.d0
  
  ! T_cr/sigma = sqrt6
  double precision,parameter :: thrm_sigma=2.449489742783178d0
  
contains
  
  function msquared(thrm)
    ! mass^2/sigma^2
    double precision,intent(in) :: thrm ! T/sigma
    double precision msquared
    
    msquared=(thrm**2-thrm_sigma**2)/3
    
  end function msquared

  function dVdphi(phi,m2,ma2)
    ! complex conjugate of dV/dphi/sigma^3;
    ! effective potential is the same as in arXiv:1202.5851 
    ! -> slightly different; for differentiability at Phi=0, f_a=sigma/N_DW in axion potential is replaced with Phi/N_DW
    complex(kind(0d0)),intent(in) :: phi
    double precision,intent(in) :: m2,ma2
    complex(kind(0d0)) dVdphi
    double precision a,b,s2 !!,Nth

    a=dble(phi)
    b=aimag(phi)
    s2=a*a+b*b
    dVdphi=(2*s2+m2)*phi
    
    ! add axion potential
    dVdphi=dVdphi-ma2/2 ! this is the potential adopted in Saikawa+ arXiv:1202.5851
    
  end function dVdphi

  function maxion2(thrm)
    ! squared mass of axion in unit of sigma
    double precision,intent(in) :: thrm
    double precision maxion2,thrm12 !!$,mm1,mm2
    
    thrm12=kappa*(cT_axion_IILM/c0_axion)**(1d0/n_axion_IILM)
    if(thrm>thrm12) then
       maxion2=cT_axion_IILM*(kappa/thrm)**n_axion_IILM
    else
       maxion2=c0_axion
    endif
    maxion2=maxion2*kappa**4
    
    !!$mm1=cT_axion_IILM*(kappa/thrm)**n_axion_IILM
    !!$mm2=c0_axion
    !!$maxion2=min(mm1,mm2)*kappa**4
    
  end function maxion2
  
  function eratio_wall_string(xi,thrm,tsigma)
    ! energy density ratio of walls to strings
    double precision,intent(in) :: xi,thrm,tsigma ! tsigma is t*sigma
    double precision eratio_wall_string
    double precision,parameter :: c=9.23d0/pi

    eratio_wall_string=tsigma*c*sqrt(maxion2(thrm))/log(tsigma/sqrt(xi))
  end function eratio_wall_string
  
end module potential
