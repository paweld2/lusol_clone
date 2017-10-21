************************************************************************
*
*     File  testlusol.f
*
*     testlusol (main)   Aset     
*
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*     This program solves Ax = b using LUSOL.
*
*     06 Apr 2002: First version, for SIAM Optimization 2002.
*                  Michael Saunders, Dept of MS&E, Stanford University.
*                  Derived from testma48.f.
*
*                  Data is read in Harwell-Boeing format.
*     23 Apr 1998: Iterative refinement added.
*     20 Mar 2000: ||A|| = 1.  dnrmx (infinity norm) used everywhere.
*                  All printed norms are absolute, not relative.
*-----------------------------------------------------------------------

      program          testlusol
      implicit	       none

!     Externals      

      integer          idamax
      double precision dnrmx
      external         dnrmx , idamax

*     Storage for data in Harwell-Boeing format.

      character        title*72, key*8
      integer          maxm, maxn, maxnz
      parameter	      (maxn   =  25000,
     $                 maxm   =   maxn,
     $                 maxnz  = 500000)
      integer          pointr( maxn+1)
      integer          rowind( maxnz )
      double precision values( maxnz ), rhsval( maxm ), xexact( maxn )

*     Storage for LUSOL.

      integer          lena
      parameter	      (lena = 20000000)

      double precision parmlu(30),  a(lena), w(maxn)
      integer          luparm(30),
     &                 indc(lena),  indr(lena),
     &                 lenc(maxn),  lenr(maxm),
     &                 ip(maxm),    iq(maxn),
     &                 iploc(maxn), iqloc(maxm),
     &                 ipinv(maxm), iqinv(maxn),
     &                 locc(maxn),  locr(maxm)

*     Local storage.

      character*8      ProbName
      integer          ijmax , inform, itn   ,
     &                 iname , iperm , iprint, ispecs, lenL  , lenU,
     &                 m     , maxitn, mode  , matfil, msgfil,
     &                 n     , nelem , nnzero, Scale , TPiv
      double precision Aijmax, factol, scltol,
     &                 bnorm , dxnorm, enorm , rnorm , rnorm0, 
     &                 snorm , xnorm , xenorm,
     &                 time1 , time2 , timeF , timeR , timeS
      double precision rscale(maxm), cscale(maxn),
     $                 rhs(maxm), w1(maxm), w2(maxm)  
      double precision r(maxn), x(maxn), dx(maxn)

      double precision zero         ,  one
      parameter       (zero = 0.0d+0,  one = 1.0d+0)

*     ------------------------------------------------------------------
*     Define files.
*     ------------------------------------------------------------------
      msgfil = 6
      iname  = 9
      ispecs = 10
      matfil = 12
      iprint = 20
      iperm  = 0    ! 30 if we want to save the permutation

      read (iname ,'(a)') ProbName
      write(msgfil, 1000) ProbName
      write(iprint, 1000) ProbName

      read (ispecs, 1100) Scale
      read (ispecs, 1110) factol
      read (ispecs, 1100) TPiv
      write(msgfil, 1200) Scale, factol, TPiv
      write(iprint, 1200) Scale, factol, TPiv

*     Initialize timer.

      call m1cpu ( 1, time1 )

*     ------------------------------------------------------------------
*     Read a matrix A in Harwell-Boeing format.
*     A is loaded into the column list (pointr, rowind, values).
*     ------------------------------------------------------------------
      call HBread( matfil, msgfil, title , key   ,
     $             maxm  , maxn  , maxnz ,
     $             m     , n     , nnzero,
     $             pointr, rowind, values, rhsval, xexact )

      write(iprint, 1300) title
      write(msgfil, 2000) n, nnzero
      write(iprint, 2000) n, nnzero

      if (m .gt. maxm  .or.  n .gt. maxn  .or.  nnzero .gt. maxnz) then
         write(msgfil, 1320)
         write(iprint, 1320)
         go to 900
      end if

*     ------------------------------------------------------------------
*     If requested, scale the matrix.
*     For this test, we just want an A with ||A|| = 1,
*     so we don't get round to unscaling.
*     ------------------------------------------------------------------
      if (Scale .gt. 0) then
         scltol = 0.95
         call gmscal( 'Set    ', m     , n     , nnzero,
     $        values   , rowind, pointr, 
     $        iprint   , scltol, 
     $        rscale   , cscale, w1    , w2 )

         call gmscal( 'Scale  ', m     , n     , nnzero,
     $        values   , rowind, pointr, 
     $        iprint   , scltol, 
     $        rscale   , cscale, w1    , w2 )

         ! TEST: Make sure ||A|| = 1 so the delta values are sensible.

         ijmax  = idamax( nnzero, values, 1 )
         Aijmax = abs( values(ijmax) )
         call dscal ( nnzero, (one/Aijmax), values, 1 )
      end if

      write(iprint, '(a)') ' ' ! to separate lu1fac output

*     ------------------------------------------------------------------
*     Set xexact, then rhsval = A*xexact
*     (over-riding any known HB rhs and solution).
*     ------------------------------------------------------------------
      call dload ( n, one , xexact, 1 )
      call dload ( n, zero, rhsval, 1 )
      call Aprod1( n     , n     , nnzero,
     $             pointr, rowind, values,
     $             xexact, rhsval )
      bnorm  = dnrmx ( n, rhsval, 1 )

      !------------------------------------------------------------------
      ! Set parameters for LUSOL's lu1fac.
      !------------------------------------------------------------------
      luparm(1) = iprint   ! File number for printed messages
      luparm(2) = 10       ! Print level. >= 0 to get singularity info.
                           !              >=10 to get more LU statistics.
                           !              >=50 to get info on each pivot.
      luparm(3) = 5        ! maxcol
      luparm(6) = TPiv     ! Threshold Pivoting: 0 = TPP, 1 = TRP, 2 = TCP
      luparm(8) = 1        ! keepLU
      parmlu(1) = factol   ! Ltol1:  max |Lij| during Factor
      parmlu(2) = factol   ! Ltol2:  max |Lij| during Update 
      parmlu(3) = 3.0d-13  ! small:  drop tolerance
      parmlu(4) = 3.7d-11  ! Utol1:  absolute tol for small Uii
      parmlu(5) = 3.7d-11  ! Utol2:  relative tol for small Uii
      parmlu(6) = 3.0d+0   ! Uspace: 
      parmlu(7) = 0.3d+0   ! dens1
      parmlu(8) = 0.5d+0   ! dens2

*     ------------------------------------------------------------------
*     Load A into (a, indc, indr).
*     ------------------------------------------------------------------
      call Aset  ( n     , nnzero,
     $             pointr, rowind, values,
     $             a     , indc  , indr   )

      nelem  = nnzero

      !------------------------------------------------------------------
      ! Factor  A = L U.
      !------------------------------------------------------------------
      call m1cpu ( 0, time1 )
      call lu1fac( m    , n    , nelem, lena , luparm, parmlu,
     $             a    , indc , indr , ip   , iq    ,
     $             lenc , lenr , locc , locr ,
     $             iploc, iqloc, ipinv, iqinv, w     , inform )
      call m1cpu ( 0, time2 )
      timeF  = time2 - time1

      lenL   = luparm(21)
      lenU   = luparm(22)
      write(msgfil, 2200) lenL, lenU, lenL+lenU
      write(iprint, 2200) lenL, lenU, lenL+lenU

      if (inform .gt. 1) then
         write(msgfil, 1500) inform
         write(iprint, 1500) inform

         timeS  = zero
         write(msgfil, 2500) timeF, timeS
         write(iprint, 2500) timeF, timeS
         go to 900
      end if

*     ------------------------------------------------------------------
*     SOLVE  A rhs(new) = rhs.
*     ------------------------------------------------------------------
      call dcopy ( n, rhsval, 1, rhs, 1 )
      mode   = 5

      call m1cpu ( 0, time1 )
      call lu6sol( mode, m, n, rhs, x,
     $             lena, luparm, parmlu,
     $             a, indc, indr, ip, iq,
     $             lenc, lenr, locc, locr,
     $             inform )
      call m1cpu ( 0, time2 )
      timeS  = time2 - time1

      write(msgfil, 2500) timeF, timeS
      write(iprint, 2500) timeF, timeS

*     ------------------------------------------------------------------
*     Set r = b - Ax.
*     Find norms of r and error in x.
*     ------------------------------------------------------------------
      snorm  = zero
      xnorm  = dnrmx ( n, x, 1 )
      xenorm = dnrmx ( n, xexact, 1 )

      call dcopy ( n, rhsval  , 1, r, 1 )
      call dcopy ( n, x       , 1, w, 1 )
      call dscal ( n,     (- one), w, 1 )
      call Aprod1( n     , n     , nnzero,
     $             pointr, rowind, values,
     $             w     , r     )
      rnorm  = dnrmx ( n, r, 1 )

      call dcopy ( n, x, 1, w, 1 )
      call daxpy ( n, (- one), xexact, 1, w, 1 )
      enorm  = dnrmx ( n, w, 1 )

      write(msgfil, 3000) snorm, xnorm, rnorm, enorm
      write(iprint, 3000) snorm, xnorm, rnorm, enorm

*     ------------------------------------------------------------------
*     ITERATIVE REFINEMENT.
*     ------------------------------------------------------------------
      call m1cpu ( 0, time1)
      write(msgfil, 3200)
      write(iprint, 3200)
      rnorm0 = 1.0d+30
      maxitn = 10

      do itn = 1, maxitn
*        ---------------------------------------------------------------
*        Solve A dx = r.
*        Set      x = x + dx.
*        ---------------------------------------------------------------
         call lu6sol( mode, m, n, r, dx,
     $                lena, luparm, parmlu,
     $                a, indc, indr, ip, iq,
     $                lenc, lenr, locc, locr,
     $                inform )
         call daxpy ( n, one, dx, 1, x, 1 )
         dxnorm = dnrmx ( n, dx, 1 )

*        ---------------------------------------------------------------
*        Set r = b - Ax.
*        Find max residual.
*        ---------------------------------------------------------------
         call dcopy ( n, rhsval , 1, r, 1 )
         call dcopy ( n, x      , 1, w, 1 )
         call dscal ( n, (- one),    w, 1 )
         call Aprod1( n     , n     , nnzero,
     $                pointr, rowind, values,
     $                w     , r     )
         rnorm  = dnrmx ( n, r, 1 )

         call dcopy ( n, x, 1, w, 1 )
         call daxpy ( n, (- one), xexact, 1, w, 1 )
         enorm  = dnrmx ( n, w, 1 )

         write(msgfil, 3300) itn, dxnorm, rnorm, enorm
         write(iprint, 3300) itn, dxnorm, rnorm, enorm

         if (rnorm .ge. 0.5d+0 * rnorm0) go to 800
         rnorm0 = rnorm
      end do

  800 call m1cpu ( 0, time2 )
      timeR  = time2 - time1
      write(msgfil, 3500) timeR
      write(iprint, 3500) timeR

  900 stop

 1000 format(/ ' ----------------------------------------------------'
     &       / ' Problem  ', a
     &       / ' ----------------------------------------------------')
 1100 format(i10)
 1110 format(e10.1)
 1200 format(    i10  , 2x, 'Scale:  0=Noscale  1=Scale'
     $     / 0p, f10.4, 2x, 'factol: Max Lij (> 1.0)'
     $     / 0p,   i10, 2x, 'TPiv:   0=TPP      1=TRP      2=TCP')
 1300 format(/ ' ', a)
 1320 format(/ ' The problem is too big for maxm, maxn or maxnz')
 1500 format(/ ' lu1fac error.   inform =', i5)
 2000 format(/
     $       ' A size  ', i8,   3x, ' A nonz  ', i8)
 2200 format(1p,
     $       ' L nonz  ', i8,   3x, ' U nonz  ', i8,   3x,
     $       ' L+U     ', i8)
 2500 format(' Factor  ', f8.2, 3x,
     $       ' Solve   ', f8.2)
 3000 format(1p,
     $     / ' snorm   ', e8.1, 3x, ' xnorm   ', e8.1
     $     / ' Residual', e8.1, 3x, ' Error   ', e8.1)
 3200 format(/ ' Refine    dxnorm     rnorm     enorm')
 3300 format(1p, i7, 3e10.1)
 3500 format(' Refine  ', f8.2) 

      end ! program testlusol

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine Aset  ( n     , nnzero,
     $                   pointr, rowind, values,
     $                   a     , irn   , jcn   )

      implicit		 double precision (a-h,o-z)
      integer            pointr( n+1)
      integer            rowind( nnzero ), irn( nnzero ), jcn( nnzero )
      double precision   values( nnzero ), a( nnzero )

*     ------------------------------------------------------------------
*     Aset   takes a matrix A in the column list
*        (pointr, rowind, values)
*     and defines a matrix A as a list (Aij, i, j) of ne entries in
*        (a, irn, jcn)
*
*     03 Jun 1996  First version.
*     ------------------------------------------------------------------

      parameter         (zero = 0.0d+0)

      ne     = 0

      do j = 1, n
         do l = pointr(j), pointr(j+1) - 1
            ne      = ne + 1
            a  (ne) = values(l)
            irn(ne) = rowind(l)
            jcn(ne) = j
         end do
      end do

*     end of Aset
      end

