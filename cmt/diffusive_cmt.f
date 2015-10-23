      subroutine compute_cip
      include 'SIZE'
      include 'CMTDATA'
      include 'GEOM'
      include 'MASS'
      include 'DG'
      integer   heresize
      parameter (lfq=lx1*lz1*2*ldim*lelt)
      parameter (heresize=nqq*3*lfq+toteq*2*ldim*lfq-2*lfq)
      common /CMTSURFLX/ vminus(lx1*lz1,2*ldim,lelt),
     >                    vplus(lx1*lz1,2*ldim,lelt),
     >                   notyet(heresize)
      real    notyet

      integer e,f
      real    hk,the_area

      nface=2*ndim
      nxyz=nx1*ny1*nz1
      nxz=nx1*nz1
      nfc =nxz*nface*nelt

! Lol I'm in way to big of a hurry to make this less ridiculously
! wasteful and code up a new gs_op. It won't be happening more than once
! until we have a moving mesh
      do e=1,nelt
         call copy(notyet,bm1(1,1,1,e),nxyz)
         call col2(notyet,jacmi(1,e),nxyz)
         vol=vlsum(notyet,nxyz)
         do f=1,nface
            do i=1,nxz
               vminus(i,f,e)=vol
            enddo
         enddo
      enddo

      call face_state_commo(vminus,vplus,nfc,1,flux_hndl)

      do e=1,nelt
         do f=1,nface
            the_area=vlsum(area(1,1,f,e),nxz)
            if (vplus(1,f,e) .gt. 0.0) then
               hk=min(vminus(1,f,e),vplus(1,f,e))
            else
               hk=vminus(1,f,e)
            endif
            hk=hk/the_area
!           cip(f,e)=cip_adhoc/hk*nx1**2
            cip(f,e)=cip_adhoc/hk*real(nx1**2)
         enddo
      enddo

      return
      end subroutine compute_cip

!-----------------------------------------------------------------------

      subroutine tauij_gdu(gijklu,visco,gvar,gvart,ngvar,du,dut,nu,npt,
     >                     xpflag,eq,jflux,kdir,e,mu,lambda,kond,cvg)
      include 'SIZE'
! subroutine for computing flux of a conserved variable by higher-order
! differential operators. This one computes viscous fluxes for everything
! except gas density.
! ngvar      number of primitive variables to form Gijkl
! nu         number of conserved variables
! eq         index i; LHS equation
! jflux      index j; flux direction
! kdir       index k; direction of derivative or jump in U
      integer  ngvar,nu,npt,eq,jflux,kdir,e
      real   gvar(ngvar,npt),gvart(npt,lelv,*),
     >                      du(nu,   npt),dut  (npt,nu)
      real  mu(npt),lambda(npt),kond(npt),cvg(npt)
      real visco(npt) ! you know, you should probably just
                         ! pass mu+lambda and mu-k/cv when eq=5
                         ! so you don't have to recompute them
                         ! so many times
! variables making up Gjkil terms, viscous stress tensor and total energy
! equation, compressible Navier-Stokes equations
!     gvar(1,:) or gvart(:,e,1) rho
!     gvar(2,:) or gvart(:,e,2) u
!     gvar(3,:) or gvart(:,e,3) v
!     gvar(4,:) or gvart(:,e,4) w
!     gvar(5,:) or gvart(:,e,5) p
!     gvar(6,:) or gvart(:,e,6) T
! derivatives or jumps, conserved variables, compressible Navier-Stokes
! equations
!     du(1,:) or dut(:,1) rho
!     du(2,:) or dut(:,2) rho u
!     du(3,:) or dut(:,3) rho v
!     du(4,:) or dut(:,4) rho w
!     du(5,:) or dut(:,5) rho E
      real gijklu(npt) !incremented. never give this exclusive intent
      logical   xpflag ! flag determining if variable is
                                    ! innermost for both G and dU.
      integer eijk3(3,3,3)
!     data eijk2 / 0, -1, 1, 0/
      data eijk3
     >/0,0,0,0,0,-1,0,1,0,0,0,1,0,0,0,-1,0,0,0,-1,0,1,0,0,0,0,0/
      parameter (itemp=6)

      call rzero(visco,npt)

      if (eq .lt. 5) then

         if (jflux .eq. eq-1) then
            m=kdir
            call copy(visco,lambda,npt)
            if (kdir .eq. jflux) then
               call add2s2(visco,mu,2.0,npt)
            endif
         else
            call copy(visco,mu,npt)
            if (kdir .eq. jflux) then
               m=eq-1
            else
               m=jflux
            endif
         endif

         m=m+1 ! skip density

         if (xpflag) then ! ipt innermost
            call invcol2(visco,gvart(1,e,1),npt)
            call addcol3(gijklu,visco,dut(1,m),npt)
            call subcol4(gijklu,visco,dut(1,1),gvart(1,e,m),npt)
         else
            do ipt=1,npt
               visco(ipt)=visco(ipt)/gvar(1,ipt)
            enddo
            do ipt=1,npt
               gijklu(ipt)=gijklu(ipt)-visco(ipt)*gvar(m,ipt)*du(1,ipt)
            enddo
            do ipt=1,npt
               gijklu(ipt)=gijklu(ipt)+visco(ipt)*du(m,ipt)
            enddo
         endif ! xpflag

      else ! energy equation is very different. and could use a rewrite

         if (jflux .eq. kdir) then
            kp1=kdir+1

            l=1 ! sigh
            if (xpflag) then ! ipt innermost
               call copy(visco,gvart(1,e,2),npt)
               call vsq(visco,npt)
               do m=3,ndim+1
                  call addcol3(visco,gvart(1,e,m),gvart(1,e,m),npt)
               enddo
! visco now contains uiui=2KE. only happens here
               call invcol2(visco,gvart(1,e,1),npt)
               call addcol4(gijklu,visco,kond,dut(1,l),npt)
! gdu-=(mu-k/cv)uiui*drho
               call subcol4(gijklu,visco,mu,dut(1,l),npt)
               call cmult(visco,0.5,npt)
               do ipt=1,npt
                  visco(ipt)=visco(ipt)+cvg(ipt)*gvart(ipt,e,itemp)/
     >                                           gvart(ipt,e,1) ! sigh
               enddo
! visco now contains E =(cvg*T+0.5*(ux**2+uy**2+uz**2))/rho
! gdu-=k/cv*E*drho/rho
               call subcol4(gijklu,visco,kond,dut(1,l),npt)

               call add3(visco,mu,lambda,npt)
               call invcol2(visco,gvart(1,e,1),npt)
! visco now contains mu+lambda. Much better
! but I've completely given up on doing this through math.f
               do ipt=1,npt
               gijklu(ipt)=gijklu(ipt)-visco(ipt)*gvart(ipt,e,kp1)**2*
     >                                            dut(ipt,l)
               enddo
            else ! xpflag
               do ipt=1,npt ! someone else can make this loop more clever
                  twoke=gvar(2,ipt)**2+gvar(3,ipt)**2+gvar(4,ipt)**2
                  gdu=(mu(ipt)-kond(ipt))*twoke
                  energy=cvg(ipt)*gvar(itemp,ipt)+0.5*twoke
                  gdu=gdu+(mu(ipt)+lambda(ipt))*gvar(kp1,ipt)**2
                  gijklu(ipt)=gijklu(ipt)-(kond(ipt)*energy-gdu)*
     >                                     du(l,ipt)/gvar(1,ipt)
               enddo
            endif

            call sub3(visco,mu,kond,npt) ! form mu-K/cv
            if (xpflag) then ! ipt innermost
               call invcol2(visco,gvart(1,e,1),npt)
               do l=2,ldim+1 ! both gvar and du are indexed by l
                  call addcol4(gijklu,visco,gvart(1,e,l),dut(1,l),npt)
               enddo
            else
               do ipt=1,npt
                  visco(ipt)=visco(ipt)/gvar(1,ipt)
                  gdu=0.0
                  do l=2,ldim+1 ! both gvar and du are indexed by l
                     gdu=gdu+gvar(l,ipt)*du(l,ipt)
                  enddo
                  gijklu(ipt)=gijklu(ipt)+gdu*visco(ipt)
               enddo
            endif
            call add3(visco,mu,lambda,npt)
            l=jflux+1
            if (xpflag) then ! ipt innermost
               call invcol2(visco,gvart(1,e,1),npt)
               call addcol4(gijklu,visco,gvart(1,e,l),dut(1,l),npt)
            else
               do ipt=1,npt
               gijklu(ipt)=gijklu(ipt)+visco(ipt)*gvar(l,ipt)*du(l,ipt)/
     >                                            gvar(1,ipt)
               enddo
            endif

            l=5
            call copy(visco,kond,npt)
            if (xpflag) then
               call invcol2(visco,gvart(1,e,1),npt)
               call addcol3(gijklu,visco,dut(1,l),npt)
            else
               do ipt=1,npt
                  visco(ipt)=visco(ipt)/gvar(1,ipt)
               enddo
               do ipt=1,npt
                  gijklu(ipt)=gijklu(ipt)+visco(ipt)*du(l,ipt)
               enddo
            endif

         else ! dU is off-diagonal

            call add3(visco,mu,lambda,npt)
            jp1=jflux+1
            kp1=kdir+1

            if (xpflag) then ! ipt innermost
               call invcol2(visco,gvart(1,e,1),npt)
! just clobber visco relentlessly
               call col2(visco,gvart(1,e,jp1),npt)
               call subcol4(gijklu,visco,gvart(1,e,kp1),dut(1,1),npt)
            else
               do ipt=1,npt
                  visco(ipt)=visco(ipt)/gvar(1,ipt)
               enddo
               do ipt=1,npt
                  gijklu(ipt)=gijklu(ipt)-
     >                  visco(ipt)*gvar(jp1,ipt)*gvar(kp1,ipt)*du(1,ipt)
               enddo
            endif

            do l=2,ldim+1
               lm1=l-1
               if (eijk3(jflux,kdir,lm1) .ne. 0) cycle

               if (lm1 .eq. kdir) then
                  call copy(visco,lambda,npt)
                  m=jflux+1
               else
                  call copy(visco,mu,npt)
                  m=kdir+1
               endif

               if (xpflag) then
                  call invcol2(visco,gvart(1,e,1),npt)
                  call addcol4(gijklu,visco,dut(1,l),gvart(1,e,m),npt)
               else
                  do ipt=1,npt
                     visco(ipt)=visco(ipt)/gvar(1,ipt)
                  enddo
                  do ipt=1,npt
                     gijklu(ipt)=gijklu(ipt)+visco(ipt)*gvar(m,ipt)*
     >                                                    du(l,ipt)
                  enddo
               endif

            enddo ! l

         endif ! diagonal?

      endif ! energy equation

      return
      end subroutine tauij_gdu

!-----------------------------------------------------------------------

      subroutine compute_transport_props
! get vdiff props (viscosity in imu, second viscosity in ilam, and
! thermal conductivity in iknd; second viscosity is usually -2/3
! visc, but we refuse to assume Stokes' hypothesis for the user)
! via nekasn
      include 'SIZE'
      include 'PARALLEL'
      include 'NEKUSE'
      include 'SOLN'

      integer   e

      do e=1,nelt
         ieg=lglel(e)
         do k=1,nz1
         do j=1,ny1
         do i=1,nx1
            call nekasgn(i,j,k,e)
            call uservp(i,j,k,ieg)
            vdiff(i,1,1,e,imu) = mu
            vdiff(i,1,1,e,ilam) = lambda
            vdiff(i,1,1,e,iknd) = udiff
         enddo
         enddo
         enddo
      enddo
      return
      end subroutine compute_transport_props

! shovel all this into usr files
!     prrl = 1.0/prlam

!     if(viscmodel.eq.1)then ! constant properties
!        do i=1,npts
!           vdiff(i,1,1,e,imu) = muref !defined in <casename>.usr/.rea
!           vdiff(i,1,1,e,ilam) = coeflambda*muref ! stokes hypothesis
!           cpg=vtrans(i,1,1,e,icv)+rgasref
!           vdiff(i,1,1,e,iknd)  = cpg*vdiff(i,1,1,e,imu)*prrl/
!    >                                  vtrans(i,1,1,e,icv) ! use gamma
!        enddo
!     else if (viscmodel.eq.2)then ! sutherland's law
!        s2 = suthCoef
!        s1 = reftemp
!        s12 = 1.0 + s1/s2 
!        do i=1,npts
!           rat = sqrt(t(i,1,1,e,1)/s2)*s12/(1.0 
!    >           +s1/t(i,1,1,e,1))

!           vdiff(i,1,1,e,imu) = muref*rat
!           cpg=vtrans(i,1,1,e,icv)+rgasref
!           vdiff(i,1,1,e,ilam) = coeflambda*vdiff(i,1,1,e,imu) ! Stokes hypothesis
!           vdiff(i,1,1,e,iknd)  = cpg*vdiff(i,1,1,e,imu)*prrl/
!    >                                 vtrans(i,1,1,e,icv)
!        enddo
!     else 
!        if(nio.eq.0) write(6,*)'Transport model unknown.'
!    >                         , ' viscmodel = ', viscmodel
!     endif

