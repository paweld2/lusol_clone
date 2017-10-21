************************************************************************
*
*     File  lusol7a.f
*
*     lu7add   lu7cyc   lu7elm   lu7for   lu7rnk   lu7zap
*
*     Utilities for LUSOL's update routines.
*     lu7for is the most important -- the forward sweep.
*
* 01 May 2002: Derived from LUSOL's original lu7a.f file.
* 01 May 2002: Current version of lusol7a.f.
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine lu7add( m, n, jadd, v,
     $                   lena, luparm, parmlu,
     $                   lenL, lenU, lrow, nrank,
     $                   a, indr, ip, lenr, locr,
     $                   inform, klast, vnorm )

      implicit           double precision (a-h,o-z)
      integer            luparm(30)
      double precision   parmlu(30), a(lena), v(m)
      integer            indr(lena), ip(m), lenr(m)
      integer            locr(m)

*     ------------------------------------------------------------------
*     lu7add  inserts the first nrank elements of the vector v(*)
*     as column  jadd  of  U.  We assume that  U  does not yet have any
*     entries in this column.
*     Elements no larger than  parmlu(3)  are treated as zero.
*     klast  will be set so that the last row to be affected
*     (in pivotal order) is row  ip(klast).
*
*     09 May 1988: First f77 version.
*     ------------------------------------------------------------------

      parameter        ( zero = 0.0d+0 )

      small  = parmlu(3)
      vnorm  = zero
      klast  = 0

      do 200 k  = 1, nrank
         i      = ip(k)
         if (abs( v(i) ) .le. small) go to 200
         klast  = k
         vnorm  = vnorm  +  abs( v(i) )
         leni   = lenr(i)

*        Compress row file if necessary.

         minfre = leni + 1
         nfree  = lena - lenL - lrow
         if (nfree .lt. minfre) then
            call lu1rec( m, .true., luparm, lrow, lena,
     $                   a, indr, lenr, locr )
            nfree  = lena - lenL - lrow
            if (nfree .lt. minfre) go to 970
         end if

*        Move row  i  to the end of the row file,
*        unless it is already there.
*        No need to move if there is a gap already.

         if (leni .eq. 0) locr(i) = lrow + 1
         lr1    = locr(i)
         lr2    = lr1 + leni - 1
         if (lr2    .eq.   lrow) go to 150
         if (indr(lr2+1) .eq. 0) go to 180
         locr(i) = lrow + 1

         do 140 l = lr1, lr2
            lrow       = lrow + 1
            a(lrow)    = a(l)
            j          = indr(l)
            indr(l)    = 0
            indr(lrow) = j
  140    continue

  150    lr2     = lrow
         lrow    = lrow + 1

*        Add the element of  v.

  180    lr2       = lr2 + 1
         a(lr2)    = v(i)
         indr(lr2) = jadd
         lenr(i)   = leni + 1
         lenU      = lenU + 1
  200 continue

*     Normal exit.

      inform = 0
      go to 990

*     Not enough storage.

  970 inform = 7

  990 return

      end ! subroutine lu7add

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine lu7cyc( kfirst, klast, ip )

      integer            ip(klast)

*     ------------------------------------------------------------------
*     lu7cyc performs a cyclic permutation on the row or column ordering
*     stored in ip, moving entry kfirst down to klast.
*     If kfirst .ge. klast, lu7cyc should not be called.
*     Sometimes klast = 0 and nothing should happen.
*
*     09 May 1988: First f77 version.
*     ------------------------------------------------------------------

      if (kfirst .lt. klast) then
         ifirst = ip(kfirst)

         do 100 k = kfirst, klast - 1
            ip(k) = ip(k + 1)
  100    continue

         ip(klast) = ifirst
      end if

      end ! subroutine lu7cyc

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine lu7elm( m, n, jelm, v,
     $                   lena, luparm, parmlu,
     $                   lenL, lenU, lrow, nrank,
     $                   a, indc, indr, ip, iq, lenr, locc, locr,
     $                   inform, diag )

      implicit           double precision (a-h,o-z)
      integer            luparm(30)
      double precision   parmlu(30), a(lena), v(m)
      integer            indc(lena), indr(lena), ip(m), iq(n), lenr(m)
      integer            locc(n), locr(m)

*     ------------------------------------------------------------------
*     lu7elm  eliminates the subdiagonal elements of a vector  v(*),
*     where  L*v = y  for some vector y.
*     If  jelm > 0,  y  has just become column  jelm  of the matrix  A.
*     lu7elm  should not be called unless  m  is greater than  nrank.
*
*     inform = 0 if y contained no subdiagonal nonzeros to eliminate.
*     inform = 1 if y contained at least one nontrivial subdiagonal.
*     inform = 7 if there is insufficient storage.
*
*     09 May 1988: First f77 version.
*                  No longer calls lu7for at end.  lu8rpc, lu8mod do so.
*     ------------------------------------------------------------------

      parameter        ( zero = 0.0d+0 )

      small  = parmlu(3)
      nrank1 = nrank + 1
      diag   = zero

*     Compress row file if necessary.

      minfre = m - nrank
      nfree  = lena - lenL - lrow
      if (nfree .ge. minfre) go to 100
      call lu1rec( m, .true., luparm, lrow, lena, a, indr, lenr, locr )
      nfree  = lena - lenL - lrow
      if (nfree .lt. minfre) go to 970

*     Pack the subdiagonals of  v  into  L,  and find the largest.

  100 vmax   = zero
      kmax   = 0
      l      = lena - lenL + 1

      do 200 k = nrank1, m
         i       = ip(k)
         vi      = abs( v(i) )
         if (vi .le. small) go to 200
         l       = l - 1
         a(l)    = v(i)
         indc(l) = i
         if (vmax .ge. vi ) go to 200
         vmax    = vi
         kmax    = k
         lmax    = l
  200 continue

      if (kmax .eq. 0) go to 900

*     ------------------------------------------------------------------
*     Remove  vmax  by overwriting it with the last packed  v(i).
*     Then set the multipliers in  L  for the other elements.
*     ------------------------------------------------------------------

      imax       = ip(kmax)
      vmax       = a(lmax)
      a(lmax)    = a(l)
      indc(lmax) = indc(l)
      l1         = l + 1
      l2         = lena - lenL
      lenL       = lenL + (l2 - l)

      do 300 l = l1, l2
         a(l)    = - a(l) / vmax
         indr(l) =   imax
  300 continue

*     Move the row containing vmax to pivotal position nrank + 1.

      ip(kmax  ) = ip(nrank1)
      ip(nrank1) = imax
      diag       = vmax

*     ------------------------------------------------------------------
*     If jelm is positive, insert  vmax  into a new row of  U.
*     This is now the only subdiagonal element.
*     ------------------------------------------------------------------

      if (jelm .gt. 0) then
         lrow       = lrow + 1
         locr(imax) = lrow
         lenr(imax) = 1
         a(lrow)    = vmax
         indr(lrow) = jelm
      end if

      inform = 1
      go to 990

*     No elements to eliminate.

  900 inform = 0
      go to 990

*     Not enough storage.

  970 inform = 7

  990 return

      end ! subroutine lu7elm

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine lu7for( m, n, kfirst, klast,
     $                   lena, luparm, parmlu,
     $                   lenL, lenU, lrow,
     $                   a, indc, indr, ip, iq, lenr, locc, locr,
     $                   inform, diag )

      implicit           double precision (a-h,o-z)
      integer            luparm(30)
      double precision   parmlu(30), a(lena)
      integer            indc(lena), indr(lena), ip(m), iq(n), lenr(m)
      integer            locc(n), locr(m)

*     ------------------------------------------------------------------
*     lu7for  (forward sweep) updates the LU factorization  A = L*U
*     when row  iw = ip(klast)  of  U  is eliminated by a forward
*     sweep of stabilized row operations, leaving  ip * U * iq  upper
*     triangular.
*
*     The row permutation  ip  is updated to preserve stability and/or
*     sparsity.  The column permutation  iq  is not altered.
*
*     kfirst  is such that row  ip(kfirst)  is the first row involved
*     in eliminating row  iw.  (Hence,  kfirst  marks the first nonzero
*     in row  iw  in pivotal order.)  If  kfirst  is unknown it may be
*     input as  1.
*
*     klast   is such that row  ip(klast)  is the row being eliminated.
*     klast   is not altered.
*
*     lu7for  should be called only if  kfirst .le. klast.
*     If  kfirst = klast,  there are no nonzeros to eliminate, but the
*     diagonal element of row  ip(klast)  may need to be moved to the
*     front of the row.
*
*     On entry,  locc(*)  must be zero.
*
*     On exit:
*     inform = 0  if row iw has a nonzero diagonal (could be small).
*     inform = 1  if row iw has no diagonal.
*     inform = 7  if there is not enough storage to finish the update.
*
*     On a successful exit (inform le 1),  locc(*)  will again be zero.
*
*        Jan 1985: Final f66 version.
*     09 May 1988: First f77 version.
*     ------------------------------------------------------------------

      parameter        ( zero = 0.0d+0 )

      double precision   Ltol
      logical            swappd

      Ltol   = parmlu(2)
      small  = parmlu(3)
      uspace = parmlu(6)
      kbegin = kfirst
      swappd = .false.

*     We come back here from below if a row interchange is performed.

  100 iw     = ip(klast)
      lenw   = lenr(iw)
      if (lenw   .eq.   0  ) go to 910
      lw1    = locr(iw)
      lw2    = lw1 + lenw - 1
      jfirst = iq(kbegin)
      if (kbegin .ge. klast) go to 700

*     Make sure there is room at the end of the row file
*     in case row  iw  is moved there and fills in completely.

      minfre = n + 1
      nfree  = lena - lenL - lrow
      if (nfree .lt. minfre) then
         call lu1rec( m, .true., luparm, lrow, lena,
     $                a, indr, lenr, locr )
         lw1    = locr(iw)
         lw2    = lw1 + lenw - 1
         nfree  = lena - lenL - lrow
         if (nfree .lt. minfre) go to 970
      end if

*     Set markers on row  iw.

      do 120 l = lw1, lw2
         j       = indr(l)
         locc(j) = l
  120 continue


*     ==================================================================
*     Main elimination loop.
*     ==================================================================
      kstart = kbegin
      kstop  = min( klast, n )

      do 500 k  = kstart, kstop
         jfirst = iq(k)
         lfirst = locc(jfirst)
         if (lfirst .eq. 0) go to 490

*        Row  iw  has its first element in column  jfirst.

         wj     = a(lfirst)
         if (k .eq. klast) go to 490

*        ---------------------------------------------------------------
*        We are about to use the first element of row  iv
*               to eliminate the first element of row  iw.
*        However, we may wish to interchange the rows instead,
*        to preserve stability and/or sparsity.
*        ---------------------------------------------------------------
         iv     = ip(k)
         lenv   = lenr(iv)
         lv1    = locr(iv)
         vj     = zero
         if (lenv      .eq.   0   ) go to 150
         if (indr(lv1) .ne. jfirst) go to 150
         vj     = a(lv1)
         if (            swappd               ) go to 200
         if (Ltol * abs( wj )  .lt.  abs( vj )) go to 200
         if (Ltol * abs( vj )  .lt.  abs( wj )) go to 150
         if (            lenv  .le.  lenw     ) go to 200

*        ---------------------------------------------------------------
*        Interchange rows  iv  and  iw.
*        ---------------------------------------------------------------
  150    ip(klast) = iv
         ip(k)     = iw
         kbegin    = k
         swappd    = .true.
         go to 600

*        ---------------------------------------------------------------
*        Delete the eliminated element from row  iw
*        by overwriting it with the last element.
*        ---------------------------------------------------------------
  200    a(lfirst)    = a(lw2)
         jlast        = indr(lw2)
         indr(lfirst) = jlast
         indr(lw2)    = 0
         locc(jlast)  = lfirst
         locc(jfirst) = 0
         lenw         = lenw - 1
         lenU         = lenU - 1
         if (lrow .eq. lw2) lrow = lrow - 1
         lw2          = lw2  - 1

*        ---------------------------------------------------------------
*        Form the multiplier and store it in the  L  file.
*        ---------------------------------------------------------------
         if (abs( wj ) .le. small) go to 490
         amult   = - wj / vj
         l       = lena - lenL
         a(l)    = amult
         indr(l) = iv
         indc(l) = iw
         lenL    = lenL + 1

*        ---------------------------------------------------------------
*        Add the appropriate multiple of row  iv  to row  iw.
*        We use two different inner loops.  The first one is for the
*        case where row  iw  is not at the end of storage.
*        ---------------------------------------------------------------
         if (lenv .eq. 1) go to 490
         lv2    = lv1 + 1
         lv3    = lv1 + lenv - 1
         if (lw2 .eq. lrow) go to 400

*        ...............................................................
*        This inner loop will be interrupted only if
*        fill-in occurs enough to bump into the next row.
*        ...............................................................
         do 350 lv = lv2, lv3
            jv     = indr(lv)
            lw     = locc(jv)
            if (lw .gt. 0) then

*              No fill-in.

               a(lw)  = a(lw)  +  amult * a(lv)
               if (abs( a(lw) ) .le. small) then

*                 Delete small computed element.

                  a(lw)     = a(lw2)
                  j         = indr(lw2)
                  indr(lw)  = j
                  indr(lw2) = 0
                  locc(j)   = lw
                  locc(jv)  = 0
                  lenU      = lenU - 1
                  lenw      = lenw - 1
                  lw2       = lw2  - 1
               end if
            else

*              Row  iw  doesn't have an element in column  jv  yet
*              so there is a fill-in.

               if (indr(lw2+1) .ne. 0) go to 360
               lenU      = lenU + 1
               lenw      = lenw + 1
               lw2       = lw2  + 1
               a(lw2)    = amult * a(lv)
               indr(lw2) = jv
               locc(jv)  = lw2
            end if
  350    continue

         go to 490

*        Fill-in interrupted the previous loop.
*        Move row  iw  to the end of the row file.

  360    lv2      = lv
         locr(iw) = lrow + 1

         do 370 l = lw1, lw2
            lrow       = lrow + 1
            a(lrow)    = a(l)
            j          = indr(l)
            indr(l)    = 0
            indr(lrow) = j
            locc(j)    = lrow
  370    continue

         lw1    = locr(iw)
         lw2    = lrow

*        ...............................................................
*        Inner loop with row  iw  at the end of storage.
*        ...............................................................
  400    do 450 lv = lv2, lv3
            jv     = indr(lv)
            lw     = locc(jv)
            if (lw .gt. 0) then

*              No fill-in.

               a(lw)  = a(lw)  +  amult * a(lv)
               if (abs( a(lw) ) .le. small) then

*                 Delete small computed element.

                  a(lw)     = a(lw2)
                  j         = indr(lw2)
                  indr(lw)  = j
                  indr(lw2) = 0
                  locc(j)   = lw
                  locc(jv)  = 0
                  lenU      = lenU - 1
                  lenw      = lenw - 1
                  lw2       = lw2  - 1
               end if
            else

*              Row  iw  doesn't have an element in column  jv  yet
*              so there is a fill-in.

               lenU      = lenU + 1
               lenw      = lenw + 1
               lw2       = lw2  + 1
               a(lw2)    = amult * a(lv)
               indr(lw2) = jv
               locc(jv)  = lw2
            end if
  450    continue

         lrow   = lw2

*        The  k-th  element of row  iw  has been processed.
*        Reset  swappd  before looking at the next element.

  490    swappd = .false.
  500 continue

*     ==================================================================
*     End of main elimination loop.
*     ==================================================================

*     Cancel markers on row  iw.

  600 lenr(iw) = lenw
      if (lenw .eq. 0) go to 910
      do 620 l = lw1, lw2
         j       = indr(l)
         locc(j) = 0
  620 continue

*     Move the diagonal element to the front of row  iw.
*     At this stage,  lenw gt 0  and  klast le n.

  700 do 720 l = lw1, lw2
         ldiag = l
         if (indr(l) .eq. jfirst) go to 730
  720 continue
      go to 910

  730 diag        = a(ldiag)
      a(ldiag)    = a(lw1)
      a(lw1)      = diag
      indr(ldiag) = indr(lw1)
      indr(lw1)   = jfirst

*     If an interchange is needed, repeat from the beginning with the
*     new row  iw,  knowing that the opposite interchange cannot occur.

      if (swappd) go to 100
      inform = 0
      go to 950

*     Singular.

  910 diag   = zero
      inform = 1

*     Force a compression if the file for  U  is much longer than the
*     no. of nonzeros in  U  (i.e. if  lrow  is much bigger than  lenU).
*     This should prevent memory fragmentation when there is far more
*     memory than necessary  (i.e. when  lena  is huge).

  950 limit  = uspace * lenU + m + n + 1000
      if (lrow .gt. limit) then
         call lu1rec( m, .true., luparm, lrow, lena,
     $                a, indr, lenr, locr )
      end if
      go to 990

*     Not enough storage.

  970 inform = 7

*     Exit.

  990 return

      end ! subroutine lu7for

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine lu7rnk( m, n, jsing,
     $                   lena, luparm, parmlu,
     $                   lenL, lenU, lrow, nrank,
     $                   a, indc, indr, ip, iq, lenr, locc, locr,
     $                   inform, diag )

      implicit           double precision (a-h,o-z)
      integer            luparm(30)
      double precision   parmlu(30), a(lena)
      integer            indc(lena), indr(lena), ip(m), iq(n), lenr(m)
      integer            locc(n), locr(m)

*     ------------------------------------------------------------------
*     lu7rnk (check rank) assumes U is currently nrank by n
*     and determines if row nrank contains an acceptable pivot.
*     If not, the row is deleted and nrank is decreased by 1.
*
*     jsing is an input parameter (not altered).  If jsing is positive,
*     column jsing has already been judged dependent.  A substitute
*     (if any) must be some other column.
*
*     -- Jul 1987: First version.
*     09 May 1988: First f77 version.
*     ------------------------------------------------------------------

      parameter        ( zero = 0.0d+0 )

      Utol1    = parmlu(4)
      diag     = zero

*     Find Umax, the largest element in row nrank.

      iw       = ip(nrank)
      lenw     = lenr(iw)
      if (lenw .eq. 0) go to 400
      l1       = locr(iw)
      l2       = l1 + lenw - 1
      Umax     = zero
      lmax     = l1

      do 100 l = l1, l2
         if (Umax .lt. abs( a(l) )) then
             Umax   =  abs( a(l) )
             lmax   =  l
         end if
  100 continue

*     Find which column that guy is in (in pivotal order).
*     Interchange him with column nrank, then move him to be
*     the new diagonal at the front of row nrank.

      diag   = a(lmax)
      jmax   = indr(lmax)

      do 300 kmax = nrank, n
         if (iq(kmax) .eq. jmax) go to 320
  300 continue

  320 iq(kmax)  = iq(nrank)
      iq(nrank) = jmax
      a(lmax)   = a(l1)
      a(l1)     = diag
      indr(lmax)= indr(l1)
      indr(l1)  = jmax

*     See if the new diagonal is big enough.

      if (Umax .le. Utol1) go to 400
      if (jmax .eq. jsing) go to 400

*     ------------------------------------------------------------------
*     The rank stays the same.
*     ------------------------------------------------------------------
      inform = 0
      return

*     ------------------------------------------------------------------
*     The rank decreases by one.
*     ------------------------------------------------------------------
  400 inform = -1
      nrank  = nrank - 1
      if (lenw .gt. 0) then

*        Delete row nrank from U.

         lenU     = lenU - lenw
         lenr(iw) = 0
         do 420 l = l1, l2
            indr(l) = 0
  420    continue

         if (l2 .eq. lrow) then

*           This row was at the end of the data structure.
*           We have to reset lrow.
*           Preceding rows might already have been deleted, so we
*           have to be prepared to go all the way back to 1.

            do 450 l = 1, l2
               if (indr(lrow) .gt. 0) go to 900
               lrow  = lrow - 1
  450       continue
         end if
      end if

  900 return

      end ! subroutine lu7rnk

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine lu7zap( m, n, jzap, kzap,
     $                   lena, lenU, lrow, nrank,
     $                   a, indr, ip, iq, lenr, locr )

      implicit           double precision (a-h,o-z)
      double precision   a(lena)
      integer            indr(lena), ip(m), iq(n), lenr(m)
      integer            locr(m)

*     ------------------------------------------------------------------
*     lu7zap  eliminates all nonzeros in column  jzap  of  U.
*     It also sets  kzap  to the position of  jzap  in pivotal order.
*     Thus, on exit we have  iq(kzap) = jzap.
*
*     -- Jul 1987: nrank added.
*     10 May 1988: First f77 version.
*     ------------------------------------------------------------------

      do 100 k  = 1, nrank
         i      = ip(k)
         leni   = lenr(i)
         if (leni .eq. 0) go to 90
         lr1    = locr(i)
         lr2    = lr1 + leni - 1
         do 50 l = lr1, lr2
            if (indr(l) .eq. jzap) go to 60
   50    continue
         go to 90

*        Delete the old element.

   60    a(l)      = a(lr2)
         indr(l)   = indr(lr2)
         indr(lr2) = 0
         lenr(i)   = leni - 1
         lenU      = lenU - 1

*        Stop if we know there are no more rows containing  jzap.

   90    kzap   = k
         if (iq(k) .eq. jzap) go to 800
  100 continue

*     nrank must be smaller than n because we haven't found kzap yet.

      do 200 k = nrank+1, n
         kzap  = k
         if (iq(k) .eq. jzap) go to 800
  200 continue

*     See if we zapped the last element in the file.

  800 if (lrow .gt. 0) then
         if (indr(lrow) .eq. 0) lrow = lrow - 1
      end if

      end ! subroutine lu7zap
