C> @file drive2_cmt.f mid-level initialization drivers. Not long for this world.
c-----------------------------------------------------------------------
      subroutine nek_cmt_init
      include 'SIZE'
      include 'TOTAL'
      include 'DG'
      include 'CMTDATA'
      if (nio.eq.0) write(6,*)'Set up CMT-Nek'    
      if (toteq.ne.5) then
         if (nio.eq.0) write(6,*)'toteq is low ! toteq = ',toteq
         if (nio.eq.0) write(6,*) 'Reset toteq in SIZE to 5'
         call exitt
      endif
      if (ifrestart) then
         ifheat = .true. ! almost certainly incorrect
      endif
      call setup_cmt_commo
      
c     call setup_cmt_param
      return
      end

!-----------------------------------------------------------------------

      subroutine izero8(a,n)
      integer*8 a(1)
      do i=1,n
         a(i)=0
      enddo
      return
      end
!-----------------------------------------------------------------------
      subroutine setup_cmt_param
      INCLUDE 'SIZE'
      INCLUDE 'INPUT'
      INCLUDE 'CMTDATA'
      INCLUDE 'CMTBCDATA'

      real  MixtPerf_R_CpG, MixtPerf_T_DPR, MixtPerf_C_GRT
     >                 ,MixtPerf_Ho_CpTUVW,MixtPerf_Cp_CvR,MixtPerf_R_M
     >                 ,MixtPerf_G_CpR      
      external MixtPerf_R_CpG, MixtPerf_T_DPR, MixtPerf_C_GRT
     >                 ,MixtPerf_Ho_CpTUVW,MixtPerf_Cp_CvR,MixtPerf_R_M
     >                 ,MixtPerf_G_CpR      

      cip_adhoc=10.0
      cvgref     = param(104)
c     gmaref     = param(105)
      molmass    = param(106)
      muref      = param(107)
      coeflambda = param(108)
      suthcoef   = param(109)
      reftemp    = param(110)
      prlam      = param(111)
      pinfty     = param(112)
      rgasref    = MixtPerf_R_M(molmass,dum)
      cpgref     = MixtPerf_Cp_CvR(cvgref,rgasref)
      gmaref     = MixtPerf_G_CpR(cpgref,rgasref) 
! put these in rea file someday
      return
      end
c------------------------------------------------------------------------

      subroutine limiter
! EBDG Stuff. WHERE'S PHI????
      include 'SIZE'
      include 'TOTAL'
      include 'CMTDATA'
      include 'NEKUSE'
      parameter (lxyz=lx1*ly1*lz1)
      common /scrns/ scrent(lxyz,lelt)
      real scrent
      integer e,eg

      nxyz=lx1*ly1*lz1
      ntot=nxyz*nelt

      rgam=rgasref/(gmaref-1.0)
      do i=1,ntot
         rho=max(vtrans(i,1,1,1,irho),1.0e-10)
!        scrent(i,1)=rgam*log(pr(i,1,1,1)/(rho**gmaref))
         scrent(i,1)=log(pr(i,1,1,1)/(rho**gmaref))
      enddo
!     call dsop(scrent,'MIN',lx1,ly1,lz1)
      call copy(t(1,1,1,1,4),scrent,ntot)
! elemental entropy minimum
!     do e=1,nelt
!        se0(e)=vlmin(scrent(1,e),nxyz)
!     enddo
      do e=1,nelt
         do i=1,nxyz
         scrent(i,e)=pr(i,1,1,e)-exp(se0const)*
     >                         (vtrans(i,1,1,e,irho)**gmaref)
         enddo
      enddo
      call copy(t(1,1,1,1,5),scrent,ntot)
      do e=1,nelt
         tauebdg(e)=vlmin(scrent(1,e),nxyz)
         tauebdg(e)=min(tauebdg(e),0.0)
         call cfill(t(1,1,1,e,6),tauebdg(e),nxyz)
      enddo
      
      do e=1,nelt
      
         do m=1,toteq
            avstate(m)=vlsc2(bm1(1,1,1,e),u(1,1,1,m,e),nxyz)/volel(e)
         enddo
         rho=avstate(1)
         e_internal=avstate(5)-
     >              0.5*(avstate(2)**2+avstate(3)**2+avstate(4)**2)/
     >                   rho
         e_internal=e_internal/rho
         eg=gllel(e)
         call cmt_userEOS(1,1,1,eg) ! assigns elm avg to  pres and temp
         epsebdg=tauebdg(e)/
     >          (tauebdg(e)-(pres-exp(se0const)*rho**gmaref))
         do i=1,nxyz
            t(i,1,1,e,3)=epsebdg
         enddo

         do m=1,toteq
            do i=1,nxyz
               uold=u(i,1,1,m,e)
               u(i,1,1,m,e)=uold+epsebdg*(avstate(m)-uold)
            enddo
         enddo

      enddo
      
      return
      end
