!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
!!
!!   Author: Prabal Negi
!!   Email : negi@mech.kth.se 
!!   Description: High pass filtering for stabilization of SEM
!!   Last Modified: 30/01/2017 
!!
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
!!----------------------------------------------------------------------  
!     Main interface for high-pass filter
      subroutine MAKE_HPF

      implicit none

      include 'SIZE'
      include 'SOLN'
      include 'INPUT'         ! param(110),(111),(112)
      include 'TSTEP'         ! ifield
      include 'MASS'          ! BM1
      include 'HPF'

      integer nxyz
      parameter (nxyz=lx1*ly1*lz1)
      integer n

      integer lm,lm2
      parameter (lm=40)
      parameter (lm2=lm*lm)
      real hpf_filter(lm2)

      integer hpf_kut
      real hpf_kai
      logical hpf_ifboyd

      integer nel

      integer icalld
      save icalld
      data icalld /0/

!     prabal
      real hpfcx(lx1,ly1,lz1,lelt), hpfcy(lx1,ly1,lz1,lelt)
      real hpfcz(lx1,ly1,lz1,lelt)

!---------------------------------------- 

      hpf_kut = int(param(110))+1
      hpf_kai = param(111)
      if (param(112).eq.0.) then
        hpf_ifboyd = .false.
      else
        hpf_ifboyd = .true.
      endif

      nel = nelv
      n = nxyz*nel

      if (hpf_kai.eq.0) then
        if (ifield.eq.1) then
          call opzero(hpfx,hpfy,hpfz)
        else
          call rzero(hpft(1,1,1,1,ifield-1),n)
        endif

        return
      endif

      if (icalld.eq.0) then

        if (hpf_kai.gt.0) then
          if (nio.eq.0) then
            write(6,*) 'Positive filtering is Numerically Unstable.'
            write(6,*) 'Remove check in hpf.f if this was intentional.'
          endif
          call exitt   
        endif

!       prabal
        call err_est_init

!       Create the filter transfer function
        call hpf_trns_fcn(hpf_filter,hpf_kut)

!       Build the matrix to apply the filter function
!       to an input field
        call build_hpf_mat(hpf_op,hpf_filter,hpf_ifboyd)

!       Only initialize once    
        icalld=icalld+1 
      endif

      if (ifield.eq.1) then
!       Apply the filter
!       to velocity fields
        call build_hpf_fld(hpfx,vx,hpf_op,nx1,nz1)
        call build_hpf_fld(hpfy,vy,hpf_op,nx1,nz1)
        if (if3d) call build_hpf_fld(hpfz,vz,hpf_op,nx1,nz1)

!       Multiply by filter weight (kai)
        call cmult(hpfx,hpf_kai,n)    
        call cmult(hpfy,hpf_kai,n)    
        if (if3d) call cmult(hpfz,hpf_kai,n)    

!       prabal
        call get_spectra(hpfcx,hpfx)
        call get_spectra(hpfcy,hpfy)
        if (if3d) call get_spectra(hpfcz,hpfz)

!       Multiply by Mass matrix 
!       and add to forcing term 
        call opadd2col (bfx,bfy,bfz,hpfx,hpfy,hpfz,bm1)

!       prabal. Testing
        if (mod(ISTEP,IOSTEP).eq.0) then
          call outpost(hpfcx,hpfcy,hpfcz,pr,t,'hpf')
        endif      

      else

!       Apply filter to temp/passive scalar fields      
        call build_hpf_fld(hpft(1,1,1,1,ifield-1),t(1,1,1,1,ifield-1),
     $       hpf_op,nx1,nz1)

!       Multiply by filter weight (kai)
        call cmult(hpft(1,1,1,1,ifield-1),hpf_kai,n)    

!       Multiply by Mass matrix    
!       and add to source term
        call addcol3(bq(1,1,1,1,ifield-1),hpft(1,1,1,1,ifield-1),bm1,n)

      endif

      return
      end subroutine MAKE_HPF

!----------------------------------------------------------------------

      subroutine build_hpf_mat(op_mat,f_filter,ifboyd)

!     Builds the operator for high pass filtering
!     Transformation matrix from nodal to modal space.
!     Applies f_filter to the the legendre coefficients
!     Transforms back to nodal space
!     Operation: V * f_filter * V^(-1)
!     Where V is the transformation matrix from modal to nodal space

      implicit none

      include 'SIZE'

      logical IFBOYD 
      integer n
      parameter (n=lx1*lx1)
      integer lm, lm2
      parameter (lm=40)
      parameter (lm2=lm*lm)

      real f_filter(lm2)
      real op_mat(lx1,lx1)

      real ref_xmap(lm2)
      real wk_xmap(lm2)

      real wk1(lm2),wk2(lm2)
      real indr(lm),ipiv(lm),indc(lm)

      real rmult(lm)
      integer ierr

      integer i,j

      call spec_coeff_init(ref_xmap,ifboyd)
      
      call copy(wk_xmap,ref_xmap,lm2)
      call copy(wk1,wk_XMAP,lm2)

      call gaujordf  (wk1,lx1,lx1,indr,indc,ipiv,ierr,rmult)  ! xmap inverse

      call mxm  (f_filter,lx1,wk1,lx1,wk2,lx1)        !          -1
      call mxm  (wk_xmap,lx1,wk2,lx1,op_mat,lx1)      !     V D V

      return
      end subroutine build_hpf_mat

c---------------------------------------------------------------------- 

      subroutine build_hpf_fld(v,u,f,nx,nz)

!     Appies the operator f to field u
!     using tensor operations
!     v = f*u

      implicit none

      include 'SIZE'
      include 'TSTEP'         ! ifield

      integer nxyz 
      parameter (nxyz=lx1*ly1*lz1)

      real w1(nxyz),w2(nxyz)              ! work arrays
      real v(nxyz,lelt),u(nxyz,lelt)      ! output and input flds
c
      integer nx,nz

      real f(nx,nx),ft(nx,nx)             ! operator f and its transpose
c
      integer e,i,j,k
      integer nel

      nel = nelv

      call copy(v,u,nxyz*nel)
c
      call transpose(ft,nx,f,nx)
c
      if (ndim.eq.3) then
        do e=1,nel
c         Filter
          call copy(w2,v(1,e),nxyz)
          call mxm(f,nx,w2,nx,w1,nx*nx)
          i=1
          j=1
          do k=1,nx
             call mxm(w1(i),nx,ft,nx,w2(j),nx)
             i = i+nx*nx
             j = j+nx*nx
          enddo
          call mxm (w2,nx*nx,ft,nx,w1,nx)

          call sub3(w2,u(1,e),w1,nxyz)
          call copy(v(1,e),w2,nxyz)
!          call copy(v(1,e),w1,nxyz)  

        enddo
      else
        do e=1,nel
c         Filter
          call copy(w1,v(1,e),nxyz)
          call mxm(f ,nx,w1,nx,w2,nx)
          call mxm(w2,nx,ft,nx,w1,nx)

          call sub3(w2,u(1,e),w1,nxyz)
          call copy(v(1,e),w2,nxyz)
!          call copy(v(1,e),w1,nxyz) 
        enddo
      endif
c
      return
      end subroutine build_hpf_fld

c---------------------------------------------------------------------- 

      subroutine hpf_trns_fcn(diag,kut)

      implicit none

      include 'SIZE'
!      include 'INPUT'
      include 'PARALLEL'

      real diag(lx1*lx1)
      real intv(lx1*lx1)
      integer nx,k0,kut,kk,k

      real amp

c     Set up transfer function
c
      nx = lx1
      call ident   (diag,nx)
      call rzero   (intv,nx*nx) 
c
      kut=kut                                     ! kut=additional modes
      k0 = nx-kut
      do k=k0+1,nx
        kk = k+nx*(k-1)
        amp = ((k-k0)*(k-k0)+0.)/(kut*kut+0.)     ! Normalized amplitude. quadratic growth
        diag(kk) = 1.-amp
        intv(kk) = amp 
      enddo

!     Output normalized transfer function
      k0 = lx1+1
      if (nio.eq.0) then
        write(6,6) 'HPF :',((intv(k)), k=1,lx1*lx1,k0)
   6    format(a8,16f9.6,6(/,8x,16f9.6))
      endif

      return
      end subroutine hpf_trns_fcn

!---------------------------------------------------------------------- 

      subroutine spec_coeff_init(ref_xmap,ifboyd)
!     Initialise spectral coefficients
!     For legendre transform

      implicit none

      include 'SIZE'
      include 'WZ'

      integer lm, lm2
      parameter (lm=40)
      parameter (lm2=lm*lm)

!     local variables
      integer i, j, k, n, nx, kj
!     Legendre polynomial
      real plegx(lm)
      real z
      real ref_xmap(lm2)
      real pht(lm2)

!     Change of basis
      logical IFBOYD
!---------------------------------------- 

      nx = LX1
      kj = 0
      n  = nx-1
      do j=1,nx
        z = ZGM1(j,1)
        call legendre_poly(plegx,z,n)
        kj = kj+1
        pht(kj) = plegx(1)
        kj = kj+1
        pht(kj) = plegx(2)

        if (IFBOYD) then        ! change basis to preserve element boundary values
          do k=3,nx
             kj = kj+1
             pht(kj) = plegx(k)-plegx(k-2)
          enddo
        else                    ! legendre basis    
          do k=3,nx
             kj = kj+1
             pht(kj) = plegx(k)
          enddo         
        endif
      enddo

      call transpose (ref_xmap,nx,pht,nx)

      return
      end subroutine spec_coeff_init

!=======================================================================

