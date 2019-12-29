SUBROUTINE aemfwd1d_genesis(mode,alt,nlyr,mvec,calc_data)
!f2py depend(nlyr) mvec
!f2py intent(in)   mode,nlyr,js,alt, mvec
!f2py intent (out) calc_data
!f2py threadsafe
!----------------------------------------------------------------------
!  Uses AEM1D_FD to compute  the frequency-domain layered earth 
!  H field for a dipole of unit moment and current.
!.
!                             INPUT
!                             -----
!        JS = station reference
!      FREQ - array of NFRQ frequencies
!     TXCLN - angle in radians that TX dipole makes with vertical (climb = +)
!     TXA90 - true for vertical co-planar briadside array
!     NSTAT - number of stations in survey line.
!       ALT - array of transmitter altitudes
!       ZRX - vertical offset of each receiver from transmitter  (below = +)
!       XRX - in-line horizontal offset of RX J;  (behind = +)
!       YRX - transverse horizontal offset of RX J (port = +)thk
!      NLYR - number of layers
!       RES - layer resistivities
!      REPS - array of relative dislectric constants
!      RMUX - mu(i) / mu(0)
!       THK - array of layer thicknesses
!     CALF, CTAU, CFREQ are the layered earth Cole-Cole parameters.
!
!------------------------------------------------------------------------
!  SIGN CONVENTION:
!  The normal layered earth field coordinate system used in thiscalc_data
!  subroutine has X (JC=1) positive along the flight path, Y (JC=2)
!  positive to starboard, and Z (JCnlyr)=3) positive down.
!------------------------------------------------------------------------
! LAST CHANGE:   19 Aug 2016   VR  
!------------------------------------------------------------------------
!   
    IMPLICIT NONE
    
    INTEGER, INTENT(in) :: nlyr, mode
    INTEGER, PARAMETER :: isw=4, iunits=6, step=1, nsx=1, npuls=1,   &
    &                     nppf=3, nchnl=11, nrx=1, kppm=0, cmp= 3,   &
    &                     ider=4, gstrp=0 , astrp=0, nrxf=1
! CMP=13 or 2?!         
! CMP = 11 => invert on horizontal in-line component only
!             = 13 => invert on vertical component only
!             = 2  => joint inversion on vertical & horizontal in-line components
!             = 3  => invert on all 3 components
!             = 4  => invert on total field

    REAL(KIND=8), PARAMETER :: pi = 3.141592654d0, twopi = 6.283185307d0
 
    REAL(KIND=8), DIMENSION(1)   :: swx
    REAL(KIND=8), DIMENSION(1,3) :: swy
    REAL(KIND=8), DIMENSION(3) ::   prm_td
    
    ! instrument geometry and constants
    ! not used :offtym = 0.d0, txfreq = 90.d0, 
    REAL(KIND=8), DIMENSION(nrxf) ::   zrx = 43.d0, xrx = 90.d0, yrx = 0.d0,           &
    &                            txcln = 0.087266462599716d0    
    REAL(KIND=8) :: txampl =52394.d0,  txfreq = 90.d0, pulse =  0.005555555555556d0
    ! txarea =  90.d0       
    ! txampl = 391.d0
    ! pulse = 0.5d0/txfreq
    ! txcln = 5 deg*pi/180
    REAL(KIND=8) , PARAMETER :: T0_min = 1.d-7
    ! windows
    REAL(KIND=8), DIMENSION(nchnl) ::  &
    &  topn =  (/ 0.01d0, 0.017d0, 0.035d0, 0.069d0, 0.122d0,   &
    &       0.191d0, 0.295d0, 0.434d0, 0.660d0, 1.007d0,    & 
    &       1.510d0 /),                                     &
    &  tcls =  (/ 0.017d0, 0.035d0, 0.069d0, 0.122d0, 0.191d0, &
    &          0.295d0, 0.434d0, 0.660d0, 1.007d0, 1.510d0, &
    &         2.205d0 /)
    ! units  punit = 'fT ' 
    REAl*8, PARAMETER :: ppfac = 1.d6, bffac = 1.d6 
    
    INTEGER, PARAMETER :: mxtym =200
    REAL(KIND=8) :: t0 , extent, tbase , qtym , tq
    REAL(KIND=8), ALLOCATABLE :: trp(:), tmp(:)

    INTEGER :: ntypls , ntyrp , kk
    
    REAL(KIND=8), INTENT(in) :: alt
    REAL(KIND=8), DIMENSION(nchnl) :: topns, tclss
    REAL(KIND=8), DIMENSION(nchnl,3) :: td_out
    REAL(KIND=8), DIMENSION(3) :: norm 
    ! model 
      
    REAL(KIND=8),  INTENT(in), DIMENSION((7*nlyr)) :: mvec      
    
    REAL(KIND=8), ALLOCATABLE :: mcurrent(:)
    
    REAL(KIND=8),  INTENT(out), DIMENSION(33) :: calc_data
    REAL(KIND=8) , DIMENSION(nlyr) :: res, reps , ctau , cfreq , calf , rmu
    REAL(KIND=8) , DIMENSION(nlyr-1) :: thk

    !SAVE 
    
 
! Set system-specific parameters for CGG Genesis system
!  only caled once: current
!  Sets up interpolation times for FD -> TD tjsransform which use the
!  exact 6 points per decade frequency-domain data plus 6 per decade
!  interpolated values.  These are based on a 12 point per decade
!  cosine filter derived from the Niels Christensen routine FILCOA
!  with  OMEGA = .3 PI and shift 0.
!  Calulates:REAL(KIND=8) swx(nsx) , swy(nsx,3)
!        TRP - array of time values for FD -> TD transformations
!      NTYRP - number of values in TRP
!     EXTENT - the latest time for which time-domain output is required.
!      PULSE - time length of one signal pulse
!     NTYPLS - number of TRP values in 1 PULSE
  
      topns = 1.d-3*topn
      tclss = 1.d-3*tcls
   
      swx(1) = 0.d0
      swy(1,1) = txampl

       ALLOCATE (tmp(mxtym))
       tmp = 0.d0
 
       qtym = dexp(dlog(10.D0)/12.D0)
       extent = 2.0d0*npuls*pulse
 

       t0 = max(minval(topns) - swx(nsx),T0_min)
       tbase = 1.D0/twopi
       DO kk = 1 , mxtym
          IF ( tbase<t0 ) EXIT
          tbase = tbase/qtym
       ENDDO
      
       tq = tbase
       tmp(1) = tq
       DO kk = 2 , mxtym
         ntyrp = kk
         tq = tq*qtym
         tmp(kk) = tq
         IF ( tmp(kk)<pulse ) ntypls = kk + 2
         IF ( tmp(kk)>extent ) EXIT
       ENDDO
 
       ALLOCATE (trp(ntyrp))
       trp(1:ntyrp) = tmp(1:ntyrp)
       DEALLOCATE (tmp)
            

       
       CALL setup_genesis(xrx,yrx,zrx,txcln,prm_td,kppm,ppfac,norm)
       
       ALLOCATE (mcurrent(7*nlyr))
       mcurrent = mvec
       
! unpack parameter vector
      CALL untrans_mvec(mode,nlyr,mcurrent)
      CALL unpack_mvec(nlyr,res,reps,rmu,calf,ctau,cfreq,thk,mcurrent)
!      WRITE(*,*) nlyr,mode
!      WRITE(*,'(a4,3x,3G12.4)') 'RES ',res  !,reps,rmu,calf,ctau,cfreq,thk

      CALL aem1d_td(step,ider,nsx,swx,swy,npuls,pulse,ntypls,     &
     &                     ntyrp,trp,nchnl,topns,tclss,txcln,alt, &
     &                     zrx,xrx,yrx,                           &
     &                     nlyr,res,reps,rmu,thk,calf,ctau,cfreq, &
     &                     gstrp,astrp,td_out)
!------------------------------------------------------
!      WRITE(*,'(/A/11G12.4/G12.4)') 'in-line ',td_out(1:nchnl,1)*norm(1)
!      WRITE(*,'(/A/11G12.4/G12.4)') 'vertical ',td_out(1:nchnl,3)*norm(3)
!      WRITE(*,'(/A/11G12.4/G12.4)') 'cross-line ',td_out(1:nchnl,2)*norm(2)
 
      
     ! in-line component (fT)
         calc_data(1:nchnl)            = td_out(1:nchnl,1)*norm(1)
     ! vertical component (fT) 
         calc_data(1*nchnl+1:2*nchnl)  = td_out(1:nchnl,3)*norm(3)
      !WRITE(*,'(/A/11G12.4)') ' result is: ',calc_data
         calc_data(2*nchnl+1:3*nchnl)  = td_out(1:nchnl,2)*norm(2)

        DEALLOCATE (mcurrent)
        
END SUBROUTINE aemfwd1d_genesis
    
     
SUBROUTINE aemjac1d_genesis(mode,alt,nlyr,mvec,calc_data,npara,isactive,jacobian)
!f2py depend(nlyr)  mvec,isactive
!f2py depend(npara) jacobian
!f2py intent(in)    nlyr,js,alt,npara,isactive, mvec
!f2py intent (out)  calc_data, jacobian
!f2py threadsafe
!----------------------------------------------------------------------
!  Uses AEM1D_FD to compute  the frequency-domain layered earth 
!  H field for a dipole of unit moment and current.
!.
!                             INPUT 
!                             -----
!        JS = station reference
!      FREQ - array of NFRQ frequencies
!     TXCLN - angle in radians that TX dipole makes with vertical (climb = +)
!     TXA90 - true for vertical co-planar briadside array
!     NSTAT - number of stations in survey line.
!       ALT - array of transmitter altitudes
!       ZRX - vertical offset of each receiver from transmitter  (below = +)
!       XRX - in-line horizontal offset of RX J;  (behind = +)
!       YRX - transverse horizontal offset of RX J (port = +)
!      NLYR - number of layers
!      RES - layer resistivities (nlyr)
!      REPS - array of relative dislectric constants (nlyr)
!      RMUX - mu(i) / mu(0)   (nlyr)
!      CALF, CTAU, CFREQ are the layered earth Cole-Cole parameters. (nlyr)
!      THK - array of layer thicknesses (nlyr-1)       
!
!------------------------------------------------------------------------
!  SIGN CONVENTION:
!  The normal layered earth field coordinate system used in this
!  subroutine has X (JC=1) positive along the flight path, Y (JC=2)
!  positive to starboard, and Z (JC=3) positive down.
!------------------------------------------------------------------------
! LAST CHANGE:   15 Jul 2017   VR  
!------------------------------------------------------------------------

    IMPLICIT NONE
    
    INTEGER, INTENT(in) :: nlyr, npara, mode
    
    INTEGER, PARAMETER :: isw=4, iunits=6, step=1, nsx=1, npuls=1,   &
    &                     nppf=3, nchnl=11, nrx=1, kppm=0, cmp= 3,   &
    &                     ider=4, gstrp=0 , astrp=0, nrxf=1
! CMP=13 or 2?!         
! CMP = 11 => invert on horizontal in-line component only
!             = 13 => invert on vertical component only
!             = 2  => joint inversion on vertical & horizontal in-line components
!             = 3  => invert on all 3 components
!             = 4  => invert on total field

    REAL(KIND=8), PARAMETER :: pi = 3.141592654d0, twopi = 6.283185307d0
 
    REAL(KIND=8), DIMENSION(1)   :: swx
    REAL(KIND=8), DIMENSION(1,3) :: swy
    REAL(KIND=8), DIMENSION(3) ::   prm_td
    
    ! instrument geometry and constants
    ! not used :offtym = 0.d0, txfreq = 90.d0, 
    REAL(KIND=8), DIMENSION(nrxf) ::   zrx = 43.d0, xrx = 90.d0, yrx = 0.d0,           &
    &                            txcln = 0.087266462599716d0    
    REAL(KIND=8) :: txampl =52394.d0,  txfreq = 90.d0, pulse =  0.005555555555556d0
    ! txarea =  90.d0       
    ! txampl = 391.d0
    ! pulse = 0.5d0/txfreq
    ! txcln = 5 deg*pi/180
    REAL(KIND=8) , PARAMETER :: T0_min = 1.d-7
    ! windows
    REAL(KIND=8), DIMENSION(nchnl) ::  &
    &  topn =  (/ 0.01d0, 0.017d0, 0.035d0, 0.069d0, 0.122d0,   &
    &       0.191d0, 0.295d0, 0.434d0, 0.660d0, 1.007d0,    & 
    &       1.510d0 /),                                     &
    &  tcls =  (/ 0.017d0, 0.035d0, 0.069d0, 0.122d0, 0.191d0, &
    &          0.295d0, 0.434d0, 0.660d0, 1.007d0, 1.510d0, &
    &         2.205d0 /)
    ! units  punit = 'fT ' 
    REAl*8, PARAMETER :: ppfac = 1.d6, bffac = 1.d6 
     
    INTEGER, PARAMETER :: mxtym =200
    REAL(KIND=8) :: t0 , extent, tbase , qtym , tq
    REAL(KIND=8), ALLOCATABLE :: trp(:), tmp(:)

    INTEGER :: ntypls , ntyrp ,  ii, jj, kk
    
    REAL(KIND=8)  :: alt
    REAL(KIND=8), DIMENSION(nchnl) :: topns, tclss
    REAL(KIND=8), DIMENSION(nchnl,3) :: td_out
    REAL(KIND=8), DIMENSION(3) :: norm 
   
    REAL(KIND=8), PARAMETER :: deltap = 1.e-3
    REAL(KIND=8), DIMENSION(7)  :: dp 
    ! model 

    INTEGER, INTENT(in), DIMENSION((7*nlyr)) :: isactive      
    REAL(KIND=8),  INTENT(in), DIMENSION((7*nlyr)) :: mvec      
    REAL(KIND=8),  INTENT(out), DIMENSION(33,npara) :: jacobian
    
    REAL(KIND=8), ALLOCATABLE :: mcurrent(:)
        
    REAL(KIND=8),  DIMENSION(33) :: calc_data, calc_data1    

    REAL(KIND=8) , DIMENSION(nlyr) :: res, reps , ctau , cfreq , calf , rmu
    REAL(KIND=8) , DIMENSION(nlyr-1) :: thk

   
    !SAVE 
    
    
!   Set system-specific parameters for CGG Genesis system
!  only caled once: 
!  Sets up interpolation times for FD -> TD transform which use the
!  exact 6 points per decade frequency-domain data plus 6 per decade
!  interpolated values.  These are based on a 12 point per decade
!  cosine filter derived from the Niels Christensen routine FILCOA
!  with  OMEGA = .3 PI and shift 0.
!  Calulates:REAL(KIND=8) swx(nsx) , swy(nsx,3)
!        TRP - array of time values for FD -> TD SUBROUTINE aem1d_tdtransformations
!      NTYRP - number of values in TRP
!     EXTENT - the latest time for which time-domain output is required.
!      PULSE - time length of one signal pulse
!     NTYPLS - number of TRP values in 1 PULSE
    
      topns = 1.d-3*topn
      tclss = 1.d-3*tcls
   
      swx(1) = 0.d0
      swy(1,1) = txampl

       ALLOCATE (tmp(mxtym))
       tmp = 0.d0
 
       qtym = dexp(dlog(10.D0)/12.D0)
       extent = 2.0d0*npuls*pulse
 
 
       t0 = max(minval(topns) - swx(nsx),T0_min)
       tbase = 1.D0/twopi
       DO kk = 1 , mxtym
          IF ( tbase<t0 ) EXIT
          tbase = tbase/qtym
       ENDDO
      
       tq = tbase
       tmp(1) = tq
       DO kk = 2 , mxtym
         ntyrp = kk
         tq = tq*qtym
         tmp(kk) = tq
         IF ( tmp(kk)<pulse ) ntypls = kk + 2
         IF ( tmp(kk)>extent ) EXIT
       ENDDO
 
 
       ALLOCATE (trp(ntyrp))
       trp(1:ntyrp) = tmp(1:ntyrp)
       
       DEALLOCATE (tmp)


      CALL setup_genesis(xrx,yrx,zrx,txcln,prm_td,kppm,ppfac,norm)
      
      ALLOCATE (mcurrent(7*nlyr))
      mcurrent = mvec
       
      CALL untrans_mvec(mode,nlyr,mcurrent)
      CALL unpack_mvec(nlyr,res,reps,rmu,calf,ctau,cfreq,thk,mcurrent)
      CALL aem1d_td(step,ider,nsx,swx,swy,npuls,pulse,ntypls,        &
     &                     ntyrp,trp,nchnl,topns,tclss,txcln,alt,    &
     &                     zrx,xrx,yrx,                              &
     &                     nlyr,res,reps,rmu,thk,calf,ctau,cfreq,    &
     &                     gstrp,astrp,td_out)
     
     ! in-line component (fT)
      calc_data(1:nchnl)            = td_out(1:nchnl,1)*norm(1)
     ! vertical component (fT) 
      calc_data(nchnl+1:2*nchnl)    = td_out(1:nchnl,3)*norm(3) 
     ! cross-line component (fT) 
      calc_data(2*nchnl+1:3*nchnl)  = td_out(1:nchnl,2)*norm(2)
     
      jj=0
      DO ii =1, size(isactive)
    
         IF (isactive(ii) /= 1) CYCLE
         
         jj=jj+1
         ! write(*,*) jj
         mcurrent = mvec
         CALL trans_mvec(mode,nlyr,mcurrent, dp)
         ! perturb mvec for paramter ii, to be sored in column jj
         mcurrent(ii) = mcurrent(ii)+dp(ii) 
         ! write(*,*) mcurrent(ii), ii,jj
         ! calculate forward model 
         CALL untrans_mvec(mode,nlyr,mcurrent)
         CALL unpack_mvec(nlyr,res,reps,rmu,calf,ctau,cfreq,thk,mcurrent)
         CALL aem1d_td(step,ider,nsx,swx,swy,npuls,pulse,ntypls,  &
     &                     ntyrp,trp,nchnl,topns,tclss,txcln,alt,    &
     &                     zrx,xrx,yrx,                           &
     &                     nlyr,res,reps,rmu,thk,calf,ctau,cfreq, &
     &                     gstrp,astrp,calc_data)
     
     ! in-line component (fT)
         calc_data1(1:nchnl)            = td_out(1:nchnl,1)*norm(1) 
     ! vertical component (fT) 
         calc_data1(1*nchnl+1:2*nchnl)  = td_out(1:nchnl,3)*norm(3)
     ! vertical component (fT) 
         calc_data1(2*nchnl+1:3*nchnl)  = td_out(1:nchnl,2)*norm(2)

    ! store divided difference into jacobian  
         jacobian(:,jj)  = (calc_data1 -calc_data)/(deltap) 
      
    ENDDO
    DEALLOCATE (trp) 
    DEALLOCATE (mcurrent)

    
    
END SUBROUTINE aemjac1d_genesis
    
SUBROUTINE setup_genesis(xrx,yrx,zrx,txcln,prm_td,kppm,ppfac,norm)
!--------------------------------------------------------- 
! For time-domain, PRM_TD is the 3 component Tx-Rx dc coupling factor
! per unit dipole moment, expressed in NANOTESLAS per unit amp
! Multiply it by dI/dt and get nanovolts per setup_genesism^2 which is the
! same as nT/s.  Multiply it by current and get nT.
!
!  SIGN CONVENTION:
!  The normal layered earth field coordinate system used in this
!  subroutine has X (JC=1) positive along the flight path, Y (JC=2)
!  positive to starboard, and Z (JC=3) positive down.
! 
!  Computes all fields in nT or nT/s.  In order to match field data for
!   inversion, it computes factors (NORM) to convert computed data into pT, pT/s,
!   fT, fT/s, pct, ppt, ppm, ppb as required.
!-----------------------------------------------------------------------------
!                             INPUT (time domain)
!                             -------------------
!      TXCLN0 - angle in degrees that the transmitting dipole makes with vertical
!              (climb = positive for VMD transmitter)
!      TXAREA - transmitter area in sq. metres
!      ZRX0, XRX0 & YRX0 are the initial vertical, in-line and transverse offsets
!                        of the receiver relative to transmitter
!                        below = + ;  behind = + ;  port = +
!
!      KPPM = 0   => No PPM normalisation (automatic if ISW = 4) ?????
!      KPPM = 1   => All components are normalised to in-line primary field
!      KPPM = 3   => All components are normalised to vertical primary field
!      KPPM = 123 => Vertical component normalised to the vertical primary &
!      KPPM = 4   => All components are normalised to total primary field
!
!      For KPPM > 0:
!        PUNIT can have values pct, ppt, ppm or ppb; eg,
!           parts per hundred (percent)
!           parts per thousand
!           parts per million
!           parts per billion
!
!      PPFAC = 1e2, 1e3, 1e6 or 1e9 is the conversion factor to achieve this
!!      PRM_TD(I) - peak primary dB/dt in nT/s if STEP = 0 or
!               peak primary B (in nT if STEP = 1)
!               I = 1, 2 & 3 are the in-line, transverse & vertical components.
!----------------------------------------------------------------------------
!                                OUTPUT (time domain)
!                                --------------------
!     PRM_TD = primary field coupling factor for B in nT (unit TX moment)
!     PRM_TD(1) = in-line B
!     PRM_TD(2) = transverse B
!     PRM_TD(3) = vertical B for station
!     NORM      = Factor to convert time-domain response in nT or nT/s into relevant units.

!     BFAC = 1.0E9 * MU / (4 * PI) for time domain  ! NANOTESLAS

      IMPLICIT NONE
      INTEGER :: kppm
      INTEGER, PARAMETER :: nrxf = 1
      REAL(KIND=8), PARAMETER ::  pi = 3.141592654d0
      REAL(KIND=8), ALLOCATABLE :: tmp(:)
      REAL(KIND=8)  :: ppfac, norm(3)
      REAL(KIND=8) :: sntx , cstx , xbd , ybd , zbd , rbrd , rsq , rsq1 , bfac ,   &
     &     fac , faczx ,  theta , inline , vert , trans , prm_td(3)
      REAL(KIND=8) , DIMENSION(nrxf) :: zrx, xrx, yrx, txcln 
 

      bfac = 100.d0
      theta = txcln(1)*pi/180.d0
      sntx = dsin(theta)
      cstx = dcos(theta)
 
      xbd = -xrx(1) !  XRX is defined as positive behind the TX.
      ybd = yrx(1)
      rbrd = dsqrt(xbd**2+ybd**2)
      zbd = zrx(1)
      rsq = zbd**2 + rbrd**2
      rsq1 = dsqrt(rsq)
      fac = bfac/rsq1**5
      faczx = 3.d0*fac*xbd*zbd
      vert = fac*cstx*(3.*zbd**2-rsq) + sntx*faczx
      inline = fac*sntx*(3.*xbd**2-rsq) + cstx*faczx
      trans = 3.d0*fac*ybd*((cstx*zbd)+(sntx*xbd))
      
      prm_td(1:3) = 0.d0
      prm_td(1) = inline
      prm_td(2) = trans
      prm_td(3) = vert
      
      
       ALLOCATE (tmp(4))
        
       tmp(1:3) = abs(prm_td(1:3))
       tmp(4)   = sqrt(tmp(1)**2+tmp(2)**2+tmp(3)**2)
        
       norm(1:3) = ppfac     
       if (kppm==1) norm(1:3) = ppfac/tmp(1)
       if (kppm==3) norm(1:3)= ppfac/tmp(3)
       if (kppm==4) norm(1:3)= ppfac/tmp(4)
       if (kppm==123) norm(1:3) = ppfac/tmp(1:3)
             
       DEALLOCATE (tmp)
       
      
END SUBROUTINE setup_genesis
       
