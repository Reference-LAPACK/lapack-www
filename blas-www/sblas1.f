      subroutine sblas()
      i = isamax(n,sx,incx)
      r = sasum(n,sx,incx)
      call saxpy(n,sa,sx,incx,sy,incy)
      call scopy(n,sx,incx,sy,incy)
      r = sdot(n,sx,incx,sy,incy)
      r = smach(job)
      r = snrm2 ( n, sx, incx)
      call srot (n,sx,incx,sy,incy,c,s)
      call srotg(sa,sb,c,s)
      call sscal(n,sa,sx,incx)
      call sswap (n,sx,incx,sy,incy)
      stop
      end
