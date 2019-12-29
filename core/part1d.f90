SUBROUTINE unpack_mvec(nlyr,res,reps,rmu,calf,ctau,cfreq,thk,mvec)
!f2py intent(in) nlyr, mvec 
!f2py depend(nlyr)  mvec,res,reps,rmu,calf,ctau,cfreq,thk
!f2py intent (out) res,reps,rmu,calf,ctau,cfreq,thk
!------------------------------------------------------------------------
! Unpacks the parameter vector to physical properties
!
!      NLYR - number of layers
!      RES - layer resistivities (nlyr)
!      REPS - array of relative dislectric constants (nlyr)
!      RMUX - mu(i) / mu(0)   (nlyr)
!      CALF, CTAU, CFREQ are the layered earth Cole-Cole parameters. (nlyr)
!      THK - array of layer thicknesses (nlyr-1)
!
!------------------------------------------------------------------------
! LAST CHANGE:   1 Sep 2016   VR  
!------------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER :: nlyr
      REAL(KIND=8) mvec(7*nlyr)
      REAL(KIND=8) thk(nlyr-1),res(nlyr),reps(nlyr),rmu(nlyr),calf(nlyr),ctau(nlyr),& 
     & cfreq(nlyr)
     
       res(1:nlyr)      = mvec(0*nlyr+1:1*nlyr)  
       rmu(1:nlyr)      = mvec(1*nlyr+1:2*nlyr)
       reps(1:nlyr)     = mvec(2*nlyr+1:3*nlyr)
       calf(1:nlyr)     = mvec(3*nlyr+1:4*nlyr)  
       ctau(1:nlyr)     = mvec(4*nlyr+1:5*nlyr)    
       cfreq(1:nlyr)    = mvec(5*nlyr+1:6*nlyr)
       thk(1:nlyr-1)    = mvec(6*nlyr+1:7*nlyr-1)

END SUBROUTINE unpack_mvec

SUBROUTINE pack_mvec(nlyr,res,reps,rmu,calf,ctau,cfreq,thk,mvec)
!f2py intent(out) mvec 
!f2py intent (in) nlyr,res,reps,rmu,calf,ctau,cfreq,thk
!f2py depend(nlyr)  mvec,res,reps,rmu,calf,ctau,cfreq,thk
!------------------------------------------------------------------------
! packs the parameter vector to physical properties
!
!      NLYR - number of layers
!      RES - layer resistivities (nlyr)
!      REPS - array of relative dislectric constants (nlyr)
!      RMUX - mu(i) / mu(0)   (nlyr)
!      CALF, CTAU, CFREQ are the layered earth Cole-Cole parameters. (nlyr)
!      THK - array of layer thicknesses (nlyr-1)
!
!------------------------------------------------------------------------
! LAST CHANGE:   1 Sep 2016   VR  
!------------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER :: nlyr 
      REAL(KIND=8) mvec(7*nlyr)
      REAL(KIND=8) thk(nlyr-1),res(nlyr),reps(nlyr),rmu(nlyr),calf(nlyr),ctau(nlyr),& 
     & cfreq(nlyr)
     
       mvec(0*nlyr+1:1*nlyr) = res(1:nlyr)
       mvec(1*nlyr+1:2*nlyr) = rmu(1:nlyr)   
       mvec(2*nlyr+1:3*nlyr) = reps(1:nlyr) 
       mvec(3*nlyr+1:4*nlyr) = calf(1:nlyr)
       mvec(4*nlyr+1:5*nlyr) = ctau(1:nlyr)
       mvec(5*nlyr+1:6*nlyr) = cfreq(1:nlyr)
       mvec(6*nlyr+1:7*nlyr-1) = thk(1:nlyr-1)

END SUBROUTINE pack_mvec
      
      

SUBROUTINE trans_mvec(mode,nlyr,mvec,dp)
!f2py intent(inout) mvec 
!f2py intent(out)   dp 
!f2py intent (in)   nlyr, mode
!f2py depend(nlyr)  mvec
!------------------------------------------------------------------------, 
! transforms the parameter vector to pre-defined functions for inversion
! generates appropriate DD pertutbations
!
!      NLYR - number of layers
!
!------------------------------------------------------------------------
! LAST CHANGE:   14 Aug 2018   VR  
!------------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER :: nlyr,  mode  
      REAL(KIND=8) mvec(7*nlyr),dp(7)
      REAL(KIND=8)  :: dp0=1.d-4
    
      dp = dp0
     
      SELECT CASE (abs(mode))
         CASE (1)
! all log Hoenig 2002 
            mvec(0*nlyr+1:1*nlyr) = dlog10(mvec(0*nlyr+1:1*nlyr))
!            mvec(1*nlyr+1:2*nlyr) =    
!            mvec(2*nlyr+1:3*nlyr) =  
            mvec(3*nlyr+1:4*nlyr) = dlog10(mvec(3*nlyr+1:4*nlyr))
            mvec(4*nlyr+1:5*nlyr) = dlog10(mvec(4*nlyr+1:5*nlyr))
            mvec(5*nlyr+1:6*nlyr) = dlog10(mvec(5*nlyr+1:6*nlyr))
            mvec(6*nlyr+1:7*nlyr) = dlog10(mvec(6*nlyr+1:7*nlyr))

        CASE (2)
! 10^mu/(10^mu + 1), INVERSE normalized chargeability Ghorbani 2007 
            mvec(0*nlyr+1:1*nlyr) = dlog10(mvec(0*nlyr+1:1*nlyr))
!            mvec(1*nlyr+1:2*nlyr) =    
!            mvec(2*nlyr+1:3*nlyr) =  
            mvec(3*nlyr+1:4*nlyr) = dlog10(mvec(3*nlyr+1:4*nlyr)/(1D0-mvec(3*nlyr+1:4*nlyr)))
            mvec(4*nlyr+1:5*nlyr) = dlog10(mvec(4*nlyr+1:5*nlyr))
            mvec(5*nlyr+1:6*nlyr) = dlog10(mvec(5*nlyr+1:6*nlyr))
            mvec(6*nlyr+1:7*nlyr) = dlog10(mvec(6*nlyr+1:7*nlyr))
    END SELECT
    
END SUBROUTINE trans_mvec

SUBROUTINE untrans_mvec(mode,nlyr,mvec)
!f2py intent(inout) mvec 
!f2py intent (in)   nlyr,mode
!f2py depend(nlyr)  mvec
!------------------------------------------------------------------------
! back-transforms the parameter vector to physical properties
!
!      NLYR - number of layers
!------------------------------------------------------------------------
! LAST CHANGE:   14 Aug 2018   VR  
!------------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER :: nlyr, mode 
      REAL(KIND=8) mvec(7*nlyr)
  
     
      
      SELECT CASE (mode)
         CASE (1)
! all log Hoenig 2002 
            mvec(0*nlyr+1:1*nlyr) = 10.d0**(mvec(0*nlyr+1:1*nlyr))
!            mvec(1*nlyr+1:2*nlyr) =    
!            mvec(2*nlyr+1:3*nlyr) =  
            mvec(3*nlyr+1:4*nlyr) = 10.d0**(mvec(3*nlyr+1:4*nlyr))
            mvec(4*nlyr+1:5*nlyr) = 10.d0**(mvec(4*nlyr+1:5*nlyr))
            mvec(5*nlyr+1:6*nlyr) = 10.d0**(mvec(5*nlyr+1:6*nlyr))
            mvec(6*nlyr+1:7*nlyr) = 10.d0**(mvec(6*nlyr+1:7*nlyr))

        CASE (2)
! 10^mu/(10^mu + 1), INVERSE normalized chargeability Ghorbani 2007 
            mvec(0*nlyr+1:1*nlyr) = 10.d0**(mvec(0*nlyr+1:1*nlyr))
!            mvec(1*nlyr+1:2*nlyr) =    
!            mvec(2*nlyr+1:3*nlyr) =  
            mvec(3*nlyr+1:4*nlyr) = 10.d0**mvec(3*nlyr+1:4*nlyr)/ (10.d0**mvec(3*nlyr+1:4*nlyr)+1.d0)  
            mvec(4*nlyr+1:4*nlyr) = 10.d0**(mvec(4*nlyr+1:5*nlyr))
            mvec(5*nlyr+1:6*nlyr) = 10.d0**(mvec(5*nlyr+1:6*nlyr))
            mvec(6*nlyr+1:7*nlyr) = 10.d0**(mvec(6*nlyr+1:7*nlyr))
    END SELECT

 
END SUBROUTINE untrans_mvec
