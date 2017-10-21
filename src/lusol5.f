!***********************************************************************
!
!     File lusol5.f
!
!     lu5prt
!
! 23 Jun 2004: First version of lu5prt.  Removing write statements
!              to make LUSOL better suited to use as a DLL.
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine lu5prt( mode, luparm, parmlu, msg, length )

      implicit           none
      integer            mode        ! Input  
      integer            luparm(30)  ! Input
      double precision   parmlu(30)  ! Input
      character          msg*(*)     ! Output
      integer            length      ! Output



      character          mnkey*1, kPiv(0:3)*2


      integer            lPiv
      logical            keepLU, TCP, TPP, TRP, TSP
      double precision   Lmax, Ltol

      double precision   zero         ,  one
      parameter        ( zero = 0.0d+0,  one = 1.0d+0 )

!     Grab relevant input parameters.

      nelem0 = nelem
      nout   = luparm(1)
      lprint = luparm(2)
      lPiv   = luparm(6)

      Ltol   = parmlu(1)    ! Limit on size of Lij

      kPiv(0)= 'PP'
      kPiv(1)= 'RP'
      kPiv(2)= 'CP'
      kPiv(3)= 'SP'
