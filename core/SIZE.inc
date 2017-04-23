c
c SIZE INTERNALS
c
      integer lelt
      parameter(lelt=lelg/lpmin + 4)

      integer lelv
      parameter(lelv=lelt)

      integer ly1,lz1
      parameter(ly1=lx1)
      parameter(lz1=1 + (ldim-2)*(lx1-1))

      integer lyd,lzd
      parameter(lyd=lxd)
      parameter(lzd=1 + (ldim-2)*(lxd-1))

      integer ly2,lz2
      parameter(ly2=lx2)
      parameter(lz2=1 + (ldim-2)*(lx2-1))

      integer ly1m,lz1m
      parameter(ly1m=lx1m)
      parameter(lz1m=1 + (ldim-2)*(lx1m-1))

      ! Averaging
      integer ay1,az1
      parameter (ay1=ax1)
      parameter (az1=1 + (ldim-2)*(ax1-1))

      integer ay2,az2
      parameter(ay2=ax2)
      parameter(az2=1 + (ldim-2)*(ax2-1))

      ! Adjoint
      integer lpelt,lpelv
      integer lpy1,lpz1
      parameter(lpelt=lelt)
      parameter(lpelv=lpelt)
      parameter(lpy1=lpx1)
      parameter(lpz1=1 + (ldim-2)*(lpx1-1))

      integer lpy2,lpz2
      parameter(lpy2=lpx2)
      parameter(lpz2=1 + (ldim-2)*(lpx2-1))

      ! MHD
      integer lbelt,lbelv
      integer lby1,lbz1
      parameter(lbelt=lelt)
      parameter(lbelv=lbelt)
      parameter(lby1=lbx1)
      parameter(lbz1=1 + (ldim-2)*(lbx1-1))

      integer lby2,lbz2
      parameter(lby2=lbx2)
      parameter(lbz2=1 + (ldim-2)*(lbx2-1))

      integer lxz
      parameter (lxz=lx1*lz1)

      integer ldimt1,ldimt3
      parameter (ldimt1=ldimt+1)
      parameter (ldimt3=ldimt+3)

      integer lzl
      parameter (lzl=3 + 2*(ldim-3))

      integer lx3,ly3,lz3
      parameter (lx3=lx1)
      parameter (ly3=ly1)
      parameter (lz3=lz1)

      integer lctmp0,lctmp1
      parameter (lctmp0 =2*lx1*ly1*lz1*lelt)
      parameter (lctmp1 =4*lx1*ly1*lz1*lelt)

      integer lvec
      parameter (lvec=1) 

      integer maxmor
      parameter (maxmor = lelt)

      integer nio
      common/IOFLAG/ nio  ! for logfile verbosity control

      integer lxo
      parameter(lxo=lx1) ! max output grid size (lxo>=lx1)

      integer lxs,lys,lzs
      parameter (lxs=1,lys=lxs,lzs=(lxs-1)*(ldim-2)+1) !New Pressure Preconditioner
      integer lfdm
      parameter (lfdm=0)  ! == 1 for fast diagonalization method

      integer lcvy1,lcvz1
      parameter(lcvy1=lcvx1)
      parameter(lcvz1=1 + (ldim-2)*(lcvx1-1))

      integer lpart
      parameter (lpart=lhis)

      integer lp, maxmbr
      parameter (lp=lpmax)
      parameter (maxmbr=lelt*6)

      integer lcvelt
      parameter (lcvelt=lelt)   ! cvode

      ! nek-nek
      integer nsessmax,nmaxl,nfldmax,nmaxcom
      parameter (nsessmax=2)
      parameter (nmaxl=lx1*lz1*lelt)
      parameter (nfldmax=4+ldimt)
      parameter (nmaxcom=nmaxl*lpratio)

      integer loglevel,optlevel
      common /lolevels/ loglevel,optlevel

      integer       nelv,nelt,nx1,ny1,nz1,nx2,ny2,nz2
     $ ,nx3,ny3,nz3,ndim,nfield,npert,nid,nxd,nyd,nzd
      common/dimn/ nelv,nelt,nx1,ny1,nz1,nx2,ny2,nz2
     $ ,nx3,ny3,nz3,ndim,nfield,npert,nid,nxd,nyd,nzd
