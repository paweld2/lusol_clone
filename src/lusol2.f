*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*
*     File  lusol2.f
*
*     Hbuild   Hchange  Hdelete  Hdown    Hinsert  Hup
*
*     Heap-management routines for LUSOL's lu1fac.
*     May be useful for other applications.
*
* 11 Feb 2002: MATLAB version derived from "Algorithms" by R. Sedgewick.
* 03 Mar 2002: F77    version derived from MATLAB version.
* 07 May 2002: Safeguard input parameters k, N, Nk.
*              We don't want them to be output!
* 19 Dec 2004: Hdelete: Nin is new input parameter for length of Hj, Ha.
* 19 Dec 2004: Current version of lusol2.f.
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!     For LUSOL, the heap structure involves three arrays of length N.
!     N        is the current number of entries in the heap.
!     Ha(1:N)  contains the values that the heap is partially sorting.
!              For LUSOL they are double precision values -- the largest
!              element in each remaining column of the updated matrix.
!              The biggest entry is in Ha(1), the top of the heap.
!     Hj(1:N)  contains column numbers j.
!              Ha(k) is the biggest entry in column j = Hj(k).
!     Hk(1:N)  contains indices within the heap.  It is the
!              inverse of Hj(1:N), so  k = Hk(j)  <=>  j = Hj(k).
!              Column j is entry k in the heap.
!     hops     is the number of heap operations,
!              i.e., the number of times an entry is moved
!              (the number of "hops" up or down the heap).
!     Together, Hj and Hk let us find values inside the heap
!     whenever we want to change one of the values in Ha.
!     For other applications, Ha may need to be some other data type,
!     like the keys that sort routines operate on.
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine Hbuild( Ha, Hj, Hk, N, Nk, hops )

      implicit
     &     none
      integer
     &     N, Nk, hops, Hj(N), Hk(Nk)
      double precision
     &     Ha(N)

*     ==================================================================
*     Hbuild initializes the heap by inserting each element of Ha.
*     Input:  Ha, Hj.
*     Output: Ha, Hj, Hk, hops.
*
*     01 May 2002: Use k for new length of heap, not k-1 for old length.
*     05 May 2002: Use kk in call to stop loop variable k being altered.
*                  (Actually Hinsert no longer alters that parameter.)
*     07 May 2002: ftnchek wants us to protect Nk, Ha(k), Hj(k) too.
*     07 May 2002: Current version of Hbuild.
*     ==================================================================

      integer
     &     h, jv, k, kk, Nkk
      double precision
     &     v

      Nkk  = Nk
      hops = 0
      do k = 1, N
         kk    = k
         v     = Ha(k)
         jv    = Hj(k)
         call Hinsert( Ha, Hj, Hk, kk, Nkk, v, jv, h )
         hops  = hops + h
      end do

      end ! subroutine Hbuild

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine Hchange( Ha, Hj, Hk, N, Nk, k, v, jv, hops )

      implicit
     &     none
      integer
     &     N, Nk, k, jv, hops, Hj(N), Hk(Nk)
      double precision
     &     v, Ha(N)

*     ==================================================================
*     Hchange changes Ha(k) to v in heap of length N.
*
*     01 May 2002: Need Nk for length of Hk.
*     07 May 2002: Protect input parameters N, Nk, k.
*     07 May 2002: Current version of Hchange.
*     ==================================================================

      integer
     &     kx, Nx, Nkx
      double precision
     &     v1

      Nx     = N
      Nkx    = Nk
      kx     = k
      v1     = Ha(k)
      Ha(k)  = v
      Hj(k)  = jv
      Hk(jv) = k
      if (v1 .lt. v) then
         call Hup   ( Ha, Hj, Hk, Nx, Nkx, kx, hops )
      else
         call Hdown ( Ha, Hj, Hk, Nx, Nkx, kx, hops )
      end if

      end ! subroutine Hchange

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine Hdelete( Ha, Hj, Hk, Nin, N, Nk, k, hops )

      implicit
     &     none
      integer
     &     N, Nin, Nk, k, hops, Hj(Nin), Hk(Nk)
      double precision
     &     Ha(Nin)

*     ==================================================================
*     Hdelete deletes Ha(k) from heap of length N.
*
*     03 Apr 2002: Current version of Hdelete.
*     01 May 2002: Need Nk for length of Hk.
*     07 May 2002: Protect input parameters N, Nk, k.
*     19 Dec 2004: Nin is new input parameter for length of Hj, Ha.
*     19 Dec 2004: Current version of Hdelete.
*     ==================================================================

      integer
     &     jv, kx, Nkx, Nx
      double precision
     &     v

      kx    = k
      Nkx   = Nk
      Nx    = N
      v     = Ha(N)
      jv    = Hj(N)
      N     = N - 1
      hops  = 0
      if (k .le. N) then
         call Hchange( Ha, Hj, Hk, Nx, Nkx, kx, v, jv, hops )
      end if

      end ! subroutine Hdelete

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine Hdown ( Ha, Hj, Hk, N, Nk, kk, hops )

      implicit
     &     none
      integer
     &     N, Nk, kk, hops, Hj(N), Hk(Nk)
      double precision
     &     Ha(N)

*     ==================================================================
*     Hdown  updates heap by moving down tree from node k.
*
*     01 May 2002: Need Nk for length of Hk.
*     05 May 2002: Change input paramter k to kk to stop k being output.
*     05 May 2002: Current version of Hdown.
*     ==================================================================

      integer
     &     j, jj, jv, k, N2
      double precision
     &     v

      k     = kk
      hops  = 0
      v     = Ha(k)
      jv    = Hj(k)
      N2    = N/2

!     while 1
  100    if (k .gt. N2   ) go to 200   ! break
         hops   = hops + 1
         j      = k+k
         if (j .lt. N) then
            if (Ha(j) .lt. Ha(j+1)) j = j+1
         end if
         if (v .ge. Ha(j)) go to 200   ! break
         Ha(k)  = Ha(j)
         jj     = Hj(j)
         Hj(k)  = jj
         Hk(jj) =  k
         k      =  j
         go to 100
!     end while

  200 Ha(k)  =  v
      Hj(k)  = jv
      Hk(jv) =  k

      end ! subroutine Hdown

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine Hinsert( Ha, Hj, Hk, N, Nk, v, jv, hops )

      implicit
     &     none
      integer
     &     N, Nk, jv, hops, Hj(N), Hk(Nk)
      double precision
     &     v, Ha(N)

*     ==================================================================
*     Hinsert inserts (v,jv) into heap of length N-1
*     to make heap of length N.
*
*     03 Apr 2002: First version of Hinsert.
*     01 May 2002: Require N to be final length, not old length.
*                  Need Nk for length of Hk.
*     07 May 2002: Protect input parameters N, Nk.
*     07 May 2002: Current version of Hinsert.
*     ==================================================================

      integer
     &     kk, Nkk, Nnew

      Nnew     = N
      Nkk      = Nk
      kk       = Nnew
      Ha(Nnew) =  v
      Hj(Nnew) = jv
      Hk(jv)   = Nnew
      call Hup   ( Ha, Hj, Hk, Nnew, Nkk, kk, hops )

      end ! subroutine Hinsert

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine Hup   ( Ha, Hj, Hk, N, Nk, kk, hops )

      implicit
     &     none
      integer
     &     N, Nk, kk, hops, Hj(N), Hk(Nk)
      double precision
     &     Ha(N)

*     ==================================================================
*     Hup updates heap by moving up tree from node k.
*
*     01 May 2002: Need Nk for length of Hk.
*     05 May 2002: Change input paramter k to kk to stop k being output.
*     05 May 2002: Current version of Hup.
*     ==================================================================

      integer
     &     j, jv, k, k2
      double precision
     &     v

      k     = kk
      hops  = 0
      v     = Ha(k)
      jv    = Hj(k)
!     while 1
  100    if (k .lt.  2    ) go to 200   ! break
         k2    = k/2
         if (v .lt. Ha(k2)) go to 200   ! break
         hops  = hops + 1
         Ha(k) = Ha(k2)
         j     = Hj(k2)
         Hj(k) =  j
         Hk(j) =  k
         k     = k2
         go to 100
!     end while

  200 Ha(k)  =  v
      Hj(k)  = jv
      Hk(jv) =  k

      end ! subroutine Hup
