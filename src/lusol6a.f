!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!     File  lusol6a.f
!
!     lu6sol   lu6L     lu6Lt     lu6U     Lu6Ut   lu6LD
!     lu6chk
!
! 26 Apr 2002: lu6 routines put into a separate file.
! 15 Dec 2002: lu6sol modularized via lu6L, lu6Lt, lu6U, lu6Ut.
!              lu6LD implemented to allow solves with LDL' or L|D|L'.
! 23 Apr 2004: lu6chk modified.  TRP can judge singularity better
!              by comparing all diagonals to DUmax.
! 27 Jun 2004: lu6chk.  Allow write only if nout .gt. 0 .
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine lu6sol( mode, m, n, v, w,
     &                   lena, luparm, parmlu,
     &                   a, indc, indr, ip, iq,
     &                   lenc, lenr, locc, locr,
     &                   inform )

      implicit
     &     none
      integer
     &     luparm(30), mode, m, n, lena, inform,
     &     indc(lena), indr(lena), ip(m), iq(n),
     &     lenc(n), lenr(m), locc(n), locr(m)
      double precision
     &     parmlu(30), a(lena), v(m), w(n)

!-----------------------------------------------------------------------
!     lu6sol  uses the factorization  A = L U  as follows:
!
!     mode
!      1    v  solves   L v = v(input).   w  is not touched.
!      2    v  solves   L'v = v(input).   w  is not touched.
!      3    w  solves   U w = v.          v  is not altered.
!      4    v  solves   U'v = w.          w  is destroyed.
!      5    w  solves   A w = v.          v  is altered as in 1.
!      6    v  solves   A'v = w.          w  is destroyed.
!
!     If mode = 3,4,5,6, v and w must not be the same arrays.
!
!     If lu1fac has just been used to factorize a symmetric matrix A
!     (which must be definite or quasi-definite), the factors A = L U
!     may be regarded as A = LDL', where D = diag(U).  In such cases,
!
!     mode
!      7    v  solves   A v = L D L'v = v(input).   w  is not touched.
!      8    v  solves       L |D| L'v = v(input).   w  is not touched.
!
!     ip(*), iq(*)      hold row and column numbers in pivotal order.
!     lenc(k)           is the length of the k-th column of initial L.
!     lenr(i)           is the length of the i-th row of U.
!     locc(*)           is not used.
!     locr(i)           is the start  of the i-th row of U.
!
!     U is assumed to be in upper-trapezoidal form (nrank by n).
!     The first entry for each row is the diagonal element
!     (according to the permutations  ip, iq).  It is stored at
!     location locr(i) in a(*), indr(*).
!
!     On exit, inform = 0 except as follows.
!     If mode = 3,4,5,6 and if U (and hence A) is singular, then
!     inform = 1 if there is a nonzero residual in solving the system
!     involving U.  parmlu(20) returns the norm of the residual.
!
!       July 1987: Early version.
!     09 May 1988: f77 version.
!     27 Apr 2000: Abolished the dreaded "computed go to".
!                  But hard to change other "go to"s to "if then else".
!     15 Dec 2002: lu6L, lu6Lt, lu6U, lu6Ut added to modularize lu6sol.
!-----------------------------------------------------------------------

      if      (mode .eq. 1) then             ! Solve  L v(new) = v.
         call lu6L  (
     &        inform, m, n, v,
     &        lena, luparm, parmlu,
     &        a, indc, indr, lenc )

      else if (mode .eq. 2) then             ! Solve  L'v(new) = v.
         call lu6Lt (
     &        inform, m, n, v,
     &        lena, luparm, parmlu,
     &        a, indc, indr, lenc )

      else if (mode .eq. 3) then             ! Solve  U w = v.
         call lu6U  (
     &        inform, m, n, v, w,
     &        lena, luparm, parmlu,
     &        a, indr, ip, iq, lenr, locr )

      else if (mode .eq. 4) then             ! Solve  U'v = w.
         call lu6Ut (
     &        inform, m, n, v, w,
     &        lena, luparm, parmlu,
     &        a, indr, ip, iq, lenr, locr )

      else if (mode .eq. 5) then             ! Solve  A w      = v
         call lu6L  (                        ! via    L v(new) = v
     &        inform, m, n, v,
     &        lena, luparm, parmlu,
     &        a, indc, indr, lenc )
         call lu6U  (                        ! and    U w = v(new).
     &        inform, m, n, v, w,
     &        lena, luparm, parmlu,
     &        a, indr, ip, iq, lenr, locr )

      else if (mode .eq. 6) then             ! Solve  A'v = w
         call lu6Ut (                        ! via    U'v = w
     &        inform, m, n, v, w,
     &        lena, luparm, parmlu,
     &        a, indr, ip, iq, lenr, locr )
         call lu6Lt (                        ! and    L'v(new) = v.
     &        inform, m, n, v,
     &        lena, luparm, parmlu,
     &        a, indc, indr, lenc )

      else if (mode .eq. 7) then
         call lu6LD (                        ! Solve  LDv(bar) = v
     &        inform, 1, m, n, v,
     &        lena, luparm, parmlu,
     &        a, indc, indr, lenc, locr )
         call lu6Lt (                        ! and    L'v(new) = v(bar).
     &        inform, m, n, v,
     &        lena, luparm, parmlu,
     &        a, indc, indr, lenc )

      else if (mode .eq. 8) then
         call lu6LD (                        ! Solve  L|D|v(bar) = v
     &        inform, 2, m, n, v,
     &        lena, luparm, parmlu,
     &        a, indc, indr, lenc, locr )
         call lu6Lt (                        ! and    L'v(new) = v(bar).
     &        inform, m, n, v,
     &        lena, luparm, parmlu,
     &        a, indc, indr, lenc )
      end if

      end ! subroutine lu6sol

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine lu6L  (
     &     inform, m, n, v,
     &     lena, luparm, parmlu,
     &     a, indc, indr, lenc )

      implicit
     &     none
      integer
     &     inform, m, n, lena, luparm(30),
     &     indc(lena), indr(lena), lenc(n)
      double precision
     &     parmlu(30), a(lena), v(m)

!     ------------------------------------------------------------------
!     lu6L   solves   L v = v(input).
!
!     15 Dec 2002: First version derived from lu6sol.
!     15 Dec 2002: Current version.
!     ------------------------------------------------------------------

      integer
     &     i, ipiv, j, k, l, l1, ldummy, len, lenL, lenL0, numL, numL0
      double precision
     &     small, vpiv

      numL0  = luparm(20)
      lenL0  = luparm(21)
      lenL   = luparm(23)
      small  = parmlu(3)
      inform = 0
      l1     = lena + 1

      do k = 1, numL0
         len   = lenc(k)
         l     = l1
         l1    = l1 - len
         ipiv  = indr(l1)
         vpiv  = v(ipiv)

         if (abs( vpiv ) .gt. small) then
            !***** This loop could be coded specially.
            do ldummy = 1, len
               l    = l - 1
               j    = indc(l)
               v(j) = v(j)  +  a(l) * vpiv
            end do
         end if
      end do

      l      = lena - lenL0 + 1
      numL   = lenL - lenL0

      !***** This loop could be coded specially.

      do ldummy = 1, numL
         l      = l - 1
         i      = indr(l)
         if (abs( v(i) ) .gt. small) then
            j    = indc(l)
            v(j) = v(j)  +  a(l) * v(i)
         end if
      end do

!     Exit.

      luparm(10) = inform

      end ! subroutine lu6L

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine lu6Lt (
     &     inform, m, n, v,
     &     lena, luparm, parmlu,
     &     a, indc, indr, lenc )

      implicit
     &     none
      integer
     &     inform, m, n, lena, luparm(30),
     &     indc(lena), indr(lena), lenc(n)
      double precision
     &     parmlu(30), a(lena), v(m)

!     ------------------------------------------------------------------
!     lu6Lt  solves   L'v = v(input).
!
!     15 Dec 2002: First version derived from lu6sol.
!     15 Dec 2002: Current version.
!     ------------------------------------------------------------------

      integer
     &     i, ipiv, j, k, l, l1, l2, len, lenL, lenL0, numL0
      double precision
     &     small, sum

!     ------------------------------------------------------------------
      double precision   zero
      parameter        ( zero = 0.0d+0 )
!     ------------------------------------------------------------------

      numL0  = luparm(20)
      lenL0  = luparm(21)
      lenL   = luparm(23)
      small  = parmlu(3)
      inform = 0
      l1     = lena - lenL + 1
      l2     = lena - lenL0

      !***** This loop could be coded specially.
      do l = l1, l2
         j     = indc(l)
         if (abs( v(j) ) .gt. small) then
            i     = indr(l)
            v(i)  = v(i)  +  a(l) * v(j)
         end if
      end do

      do k = numL0, 1, -1
         len   = lenc(k)
         sum   = zero
         l1    = l2 + 1
         l2    = l2 + len

         !***** This loop could be coded specially.
         do l = l1, l2
            j     = indc(l)
            sum   = sum  +  a(l) * v(j)
         end do

         ipiv    = indr(l1)
         v(ipiv) = v(ipiv) + sum
      end do

!     Exit.

      luparm(10) = inform

      end ! subroutine lu6Lt

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine lu6U  (
     &     inform, m, n, v, w,
     &     lena, luparm, parmlu,
     &     a, indr, ip, iq, lenr, locr )

      implicit
     &     none
      integer
     &     inform, m, n, lena, luparm(30),
     &     indr(lena), ip(m), iq(n), lenr(m), locr(m)
      double precision
     &     parmlu(30), a(lena), v(m), w(n)

!     ------------------------------------------------------------------
!     lu6U   solves   U w = v.          v  is not altered.
!
!     15 Dec 2002: First version derived from lu6sol.
!     15 Dec 2002: Current version.
!     ------------------------------------------------------------------

      integer
     &     i, j, k, klast, l, l1, l2, l3, nrank, nrank1
      double precision
     &     resid, small, t

!     ------------------------------------------------------------------
      double precision   zero
      parameter        ( zero = 0.0d+0 )
!     ------------------------------------------------------------------

      nrank  = luparm(16)
      small  = parmlu(3)
      inform = 0
      nrank1 = nrank + 1
      resid  = zero

!     Find the first nonzero in v(1:nrank), counting backwards.

      do klast = nrank, 1, -1
         i      = ip(klast)
         if (abs( v(i) ) .gt. small) go to 320
      end do

  320 do k = klast + 1, n
         j     = iq(k)
         w(j)  = zero
      end do

!     Do the back-substitution, using rows 1:klast of U.

      do k  = klast, 1, -1
         i      = ip(k)
         t      = v(i)
         l1     = locr(i)
         l2     = l1 + 1
         l3     = l1 + lenr(i) - 1

         !***** This loop could be coded specially.
         do l = l2, l3
            j     = indr(l)
            t     = t  -  a(l) * w(j)
         end do

         j      = iq(k)
         if (abs( t ) .le. small) then
            w(j)  = zero
         else
            w(j)  = t / a(l1)
         end if
      end do

!     Compute residual for overdetermined systems.

      do k = nrank1, m
         i     = ip(k)
         resid = resid  +  abs( v(i) )
      end do

!     Exit.

      if (resid .gt. zero) inform = 1
      luparm(10) = inform
      parmlu(20) = resid

      end ! subroutine lu6U

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine lu6Ut (
     &     inform, m, n, v, w,
     &     lena, luparm, parmlu,
     &     a, indr, ip, iq, lenr, locr )

      implicit
     &     none
      integer
     &     inform, m, n, lena, luparm(30),
     &     indr(lena), ip(m), iq(n), lenr(m), locr(m)
      double precision
     &     parmlu(30), a(lena), v(m), w(n)

!     ------------------------------------------------------------------
!     lu6Ut  solves   U'v = w.          w  is destroyed.
!
!     15 Dec 2002: First version derived from lu6sol.
!     15 Dec 2002: Current version.
!     ------------------------------------------------------------------

      integer
     &     i, j, k, l, l1, l2, nrank, nrank1
      double precision
     &     resid, small, t

!     ------------------------------------------------------------------
      double precision   zero
      parameter        ( zero = 0.0d+0 )
!     ------------------------------------------------------------------

      nrank  = luparm(16)
      small  = parmlu(3)
      inform = 0
      nrank1 = nrank + 1
      resid  = zero

      do k = nrank1, m
         i     = ip(k)
         v(i)  = zero
      end do

!     Do the forward-substitution, skipping columns of U(transpose)
!     when the associated element of w(*) is negligible.

      do 480 k = 1, nrank
         i      = ip(k)
         j      = iq(k)
         t      = w(j)
         if (abs( t ) .le. small) then
            v(i) = zero
            go to 480
         end if

         l1     = locr(i)
         t      = t / a(l1)
         v(i)   = t
         l2     = l1 + lenr(i) - 1
         l1     = l1 + 1

         !***** This loop could be coded specially.
         do l = l1, l2
            j     = indr(l)
            w(j)  = w(j)  -  t * a(l)
         end do
  480 continue

!     Compute residual for overdetermined systems.

      do k = nrank1, n
         j     = iq(k)
         resid = resid  +  abs( w(j) )
      end do

!     Exit.

      if (resid .gt. zero) inform = 1
      luparm(10) = inform
      parmlu(20) = resid

      end ! subroutine lu6Ut

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine lu6LD (
     &     inform, mode, m, n, v,
     &     lena, luparm, parmlu,
     &     a, indc, indr, lenc, locr )

      implicit
     &     none
      integer
     &     luparm(30), inform, mode, m, n, lena,
     &     indc(lena), indr(lena), lenc(n), locr(m)
      double precision
     &     parmlu(30), a(lena), v(m)

!-----------------------------------------------------------------------
!     lu6LD  assumes lu1fac has computed factors A = LU of a
!     symmetric definite or quasi-definite matrix A,
!     using Threshold Symmetric Pivoting (TSP),   luparm(6) = 3,
!     or    Threshold Diagonal  Pivoting (TDP),   luparm(6) = 4.
!     It also assumes that no updates have been performed.
!     In such cases,  U = D L', where D = diag(U).
!     lu6LDL returns v as follows:
!
!     mode
!      1    v  solves   L D v = v(input).
!      2    v  solves   L|D|v = v(input).
!
!     15 Dec 2002: First version of lu6LD.
!     15 Dec 2002: Current version.
!-----------------------------------------------------------------------

      ! Solve L D v(new) = v  or  L|D|v(new) = v, depending on mode.
      ! The code for L is the same as in lu6L,
      ! but when a nonzero entry of v arises, we divide by
      ! the corresponding entry of D or |D|.

      integer
     &     ipiv, j, k, l, l1, ldummy, len, numL0
      double precision
     &     diag, small, vpiv

      numL0  = luparm(20)
      small  = parmlu(3)
      inform = 0
      l1     = lena + 1

      do k = 1, numL0
         len   = lenc(k)
         l     = l1
         l1    = l1 - len
         ipiv  = indr(l1)
         vpiv  = v(ipiv)

         if (abs( vpiv ) .gt. small) then
            !***** This loop could be coded specially.
            do ldummy = 1, len
               l    = l - 1
               j    = indc(l)
               v(j) = v(j)  +  a(l) * vpiv
            end do

            ! Find diag = U(ipiv,ipiv) and divide by diag or |diag|.

            l    = locr(ipiv)
            diag = A(l)
            if (mode .eq. 2) diag = abs( diag )
            v(ipiv) = vpiv / diag
         end if
      end do

      end ! subroutine lu6LD

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine lu6chk( mode, m, n, w,
     &                   lena, luparm, parmlu,
     &                   a, indc, indr, ip, iq,
     &                   lenc, lenr, locc, locr,
     &                   inform )

      implicit
     &     none
      integer
     &     mode, m, n, lena, inform,
     &     luparm(30), indc(lena), indr(lena), ip(m), iq(n),
     &     lenc(n), lenr(m), locc(n), locr(m)
      double precision
     &     parmlu(30), a(lena), w(n)

!     ------------------------------------------------------------------
!     lu6chk  looks at the LU factorization  A = L*U.
!
!     If mode = 1, lu6chk is being called by lu1fac.
!     (Other modes not yet implemented.)
!     The important input parameters are
!
!                    lprint = luparm(2)
!                             luparm(6) = 1 if TRP
!                    keepLU = luparm(8)
!                    Utol1  = parmlu(4)
!                    Utol2  = parmlu(5)
!
!     and the significant output parameters are
!
!                    inform = luparm(10)
!                    nsing  = luparm(11)
!                    jsing  = luparm(12)
!                    jumin  = luparm(19)
!                    Lmax   = parmlu(11)
!                    Umax   = parmlu(12)
!                    DUmax  = parmlu(13)
!                    DUmin  = parmlu(14)
!                    and      w(*).
!
!     Lmax  and Umax  return the largest elements in L and U.
!     DUmax and DUmin return the largest and smallest diagonals of U
!                     (excluding diagonals that are exactly zero).
!
!     In general, w(j) is set to the maximum absolute element in
!     the j-th column of U.  However, if the corresponding diagonal
!     of U is small in absolute terms or relative to w(j)
!     (as judged by the parameters Utol1, Utol2 respectively),
!     then w(j) is changed to - w(j).
!
!     Thus, if w(j) is not positive, the j-th column of A
!     appears to be dependent on the other columns of A.
!     The number of such columns, and the position of the last one,
!     are returned as nsing and jsing.
!
!     Note that nrank is assumed to be set already, and is not altered.
!     Typically, nsing will satisfy      nrank + nsing = n,  but if
!     Utol1 and Utol2 are rather large,  nsing > n - nrank   may occur.
!
!     If keepLU = 0,
!     Lmax  and Umax  are already set by lu1fac.
!     The diagonals of U are in the top of A.
!     Only Utol1 is used in the singularity test to set w(*).
!
!     inform = 0  if  A  appears to have full column rank  (nsing = 0).
!     inform = 1  otherwise  (nsing .gt. 0).
!
!     00 Jul 1987: Early version.
!     09 May 1988: f77 version.
!     11 Mar 2001: Allow for keepLU = 0.
!     17 Nov 2001: Briefer output for singular factors.
!     05 May 2002: Comma needed in format 1100 (via Kenneth Holmstrom).
!     06 May 2002: With keepLU = 0, diags of U are in natural order.
!                  They were not being extracted correctly.
!     23 Apr 2004: TRP can judge singularity better by comparing
!                  all diagonals to DUmax.
!     27 Jun 2004: (PEG) Allow write only if nout .gt. 0 .
!     ------------------------------------------------------------------

      character
     &     mnkey
      logical
     &     keepLU, TRP
      integer
     &     i, j, jsing, jumin, k, l, l1, l2, ldiagU, lenL, lprint,
     &     ndefic, nout, nrank, nsing
      double precision
     &     aij, diag, DUmax, DUmin, Lmax, Umax, Utol1, Utol2

      double precision   zero
      parameter        ( zero = 0.0d+0 )

      nout   = luparm(1)
      lprint = luparm(2)
      TRP    = luparm(6) .eq. 1  ! Threshold Rook Pivoting
      keepLU = luparm(8) .ne. 0
      nrank  = luparm(16)
      lenL   = luparm(23)
      Utol1  = parmlu(4)
      Utol2  = parmlu(5)

      inform = 0
      Lmax   = zero
      Umax   = zero
      nsing  = 0
      jsing  = 0
      jumin  = 0
      DUmax  = zero
      DUmin  = 1.0d+30

      do j = 1, n
         w(j) = zero
      end do


      if (keepLU) then
         !--------------------------------------------------------------
         ! Find  Lmax.
         !--------------------------------------------------------------
         do l = lena + 1 - lenL, lena
            Lmax  = max( Lmax, abs(a(l)) )
         end do

         !--------------------------------------------------------------
         ! Find Umax and set w(j) = maximum element in j-th column of U.
         !--------------------------------------------------------------
         do k = 1, nrank
            i     = ip(k)
            l1    = locr(i)
            l2    = l1 + lenr(i) - 1

            do l = l1, l2
               j     = indr(l)
               aij   = abs( a(l) )
               w(j)  = max( w(j), aij )
               Umax  = max( Umax, aij )
            end do
         end do

         parmlu(11) = Lmax
         parmlu(12) = Umax

         !--------------------------------------------------------------
         ! Find DUmax and DUmin, the extreme diagonals of U.
         !--------------------------------------------------------------
         do k = 1, nrank
            j      = iq(k)
            i      = ip(k)
            l1     = locr(i)
            diag   = abs( a(l1) )
            DUmax  = max( DUmax, diag )
            if (DUmin .gt. diag) then
               DUmin  =   diag
               jumin  =   j
            end if
         end do

      else
         !--------------------------------------------------------------
         ! keepLU = 0.
         ! Only diag(U) is stored.  Set w(*) accordingly.
         ! Find DUmax and DUmin, the extreme diagonals of U.
         !--------------------------------------------------------------
         ldiagU = lena - n

         do k = 1, nrank
            j      = iq(k)
          !!diag   = abs( a(ldiagU + k) ) ! 06 May 2002: Diags
            diag   = abs( a(ldiagU + j) ) ! are in natural order
            w(j)   = diag
            DUmax  = max( DUmax, diag )
            if (DUmin .gt. diag) then
               DUmin  =   diag
               jumin  =   j
            end if
         end do
      end if


      !--------------------------------------------------------------
      ! Negate w(j) if the corresponding diagonal of U is
      ! too small in absolute terms or relative to the other elements
      ! in the same column of  U.
      !
      ! 23 Apr 2004: TRP ensures that diags are NOT small relative to
      !              other elements in their own column.
      !              Much better, we can compare all diags to DUmax.
      !--------------------------------------------------------------
      if (mode .eq. 1  .and.  TRP) then
         Utol1 = max( Utol1, Utol2*DUmax )
      end if

      if (keepLU) then
         do k = 1, n
            j     = iq(k)
            if (k .gt. nrank) then
               diag   = zero
            else
               i      = ip(k)
               l1     = locr(i)
               diag   = abs( a(l1) )
            end if

            if (diag .le. Utol1  .or.  diag .le. Utol2*w(j)) then
               nsing  =   nsing + 1
               jsing  =   j
               w(j)   = - w(j)
            end if
         end do

      else ! keepLU = 0

         do k = 1, n
            j      = iq(k)
            diag   = w(j)

            if (diag .le. Utol1) then
               nsing  =   nsing + 1
               jsing  =   j
               w(j)   = - w(j)
            end if
         end do
      end if


      !-----------------------------------------------------------------
      ! Set output parameters.
      !-----------------------------------------------------------------
      if (jumin .eq. 0) DUmin = zero
      luparm(11) = nsing
      luparm(12) = jsing
      luparm(19) = jumin
      parmlu(13) = DUmax
      parmlu(14) = DUmin

      if (nsing .gt. 0) then  ! The matrix has been judged singular.
         inform = 1
         ndefic = n - nrank
         if (nout .gt. 0  .and.  lprint .ge. 0) then
            if (m .gt. n) then
               mnkey  = '>'
            else if (m .eq. n) then
               mnkey  = '='
            else
               mnkey  = '<'
            end if
            write(nout, 1100) mnkey, nrank, ndefic, nsing
         end if
      end if

!     Exit.

      luparm(10) = inform
      return

 1100 format(' Singular(m', a, 'n)',
     &       '  rank', i9, '  n-rank', i8, '  nsing', i9)

      end ! subroutine lu6chk

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
