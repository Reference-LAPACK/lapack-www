      subroutine dblas()
      d = dasum(n,dx,incx)
      call daxpy(n,da,dx,incx,dy,incy)
      call  dcopy(n,dx,incx,dy,incy)
      d = ddot(n,dx,incx,dy,incy)
      d = dmach(job)
      d = dnrm2 ( n, dx, incx)
      call  drot (n,dx,incx,dy,incy,c,s)
      call drotg(da,db,c,s)
      call  dscal(n,da,dx,incx)
      call  dswap (n,dx,incx,dy,incy)
      i = idamax(n,dx,incx)
      stop
      end
