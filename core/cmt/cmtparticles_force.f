c-----------------------------------------------------------------------
      subroutine usr_particles_forcing
c
c     calculate the rhs of particle equation
c
      include 'SIZE'
      include 'TOTAL'
      include 'CMTDATA'
      include 'CMTPART'
      include 'PERFECTGAS'

      real vel_diff,pmass,pmassf

      common /PARTRK3/ kv_stage_p
      real kv_stage_p(llpart,7)

      ! begin timer
      ptdum(11) = dnekclock()

      pi  = 4.0d+0*atan(1.0d+0)

c     scheme 1 ------------------------------------------------------
      if (time_integ.eq.0) then 

c     RK3 ------------------------------------------------------------
      elseif (time_integ.eq.1) then 

      call calc_substantial_derivative

      do i=1,n
c        setup values ------------------------------------------------
         pmass = rpart(jvol,i)*rpart(jrhop,i)
         pmassf= rpart(jvol,i)*rpart(jrho,i)
         if(part_force(4).ne.0) pmass = pmass + rpart(jcmiu,i)*pmassf ! am

         vel_diff = sqrt((rpart(ju0  ,i)-rpart(jv0  ,i))**2+
     >                   (rpart(ju0+1,i)-rpart(jv0+1,i))**2+
     >                   (rpart(ju0+2,i)-rpart(jv0+2,i))**2)
         rpart(ja,i)  = MixtPerf_C_GRT(gmaref,rgasref,rpart(jtempf,i))
         rpart(ja,i)  = vel_diff/rpart(ja,i) ! relative mach number
         rpart(jre,i) = rpart(jrho,i)*rpart(jdp,i)*vel_diff/mu_0 ! Re
            
c        momentum rhs ------------------------------------------------
         do j=0,ndim-1
            call usr_particles_f_user(i,j)
            call usr_particles_f_qs(i,j)
            call usr_particles_f_un(i,j)
            call usr_particles_f_iu(i,j)

            rdum = 0.
            rdum = rdum + rpart(jfusr+j,i)
            rdum = rdum + rpart(jfqs+j,i)
            rdum = rdum + rpart(jfun+j,i)
            rdum = rdum + rpart(jfiu+j,i)

           rpart(jf0+j,i) = rdum/pmass ! mass weighted force
         enddo

c        energy rhs --------------------------------------------------
         call usr_particles_q_uu(i)
         call usr_particles_q_qs(i)

         rdum = 0. 
         rdum = rdum + rpart(jquu,i)
         rdum = rdum + rpart(jqqs,i)

         pmass = rpart(jvol,i)*rpart(jrhop,i)
         rpart(jq0,i) = rdum/(pmass*cp_p)

      enddo

c     other ----------------------------------------------------------
      elseif (time_integ.eq.2) then 

      endif

      ! end timer
      pttime(11) = pttime(11) + dnekclock() - ptdum(11)

      return
      end
c-----------------------------------------------------------------------
      subroutine compute_forcing_post_part
c
c     post calculate forces due to factoring of equations
c
      include 'SIZE'
      include 'TOTAL'
      include 'CMTDATA'
      include 'CMTPART'

      real uvel(0:2), vvel(0:2), pmass, pmassf,S_qs

      common /PARTRK3/ kv_stage_p
      real kv_stage_p(llpart,7)

      common /myparts/ times(0:3),alpha(0:3),beta(0:3)

      ! begin timer
      ptdum(12) = dnekclock()

      pi  = 4.0d+0*atan(1.0d+0)

c     scheme 1 ------------------------------------------------------- 
      if (time_integ.eq.0) then 

c     rk3 ------------------------------------------------------------ 
      elseif (time_integ.eq.1) then
         do i=1,n
            pmass = rpart(jvol,i)*rpart(jrhop,i)
            pmassf= rpart(jvol,i)*rpart(jrho,i)
            if (part_force(4).ne.0) pmass =pmass + rpart(jcmiu,i)*pmassf

c           momentum forcing to fluid
            do j=0,ndim-1
               rdum = 0.

               rdvdt = rpart(jf0+j,i) ! note already divided by Mp + am
               ram_s = rdvdt*rpart(jcmiu,i)*pmassf

               rpart(jfiu+j,i) = pmass*rpart(jf0+j,i) - 
     >                                            ( rpart(jfusr+j,i) +
     >                                              rpart(jfqs+j,i)  +
     >                                              rpart(jfun+j,i)  )
     >                                                - ram_s

c              note that no coupled f_un in this formulation
               rdum = rdum + rpart(jfqs+j,i)
               rdum = rdum + rpart(jfiu+j,i)

               rpart(jf0+j,i) = rdum
            enddo

c           energy forcing to fluid (quasi-steady)
            rpart(jg0,i) = rpart(jv0  ,i)*rpart(jfqs  ,i) + !force work
     >                     rpart(jv0+1,i)*rpart(jfqs+1,i) +
     >                     rpart(jv0+2,i)*rpart(jfqs+2,i)
            rpart(jg0,i) = rpart(jg0,i) + 
     >                     rpart(ju0  ,i)*rpart(jfiu  ,i) + !iu
     >                     rpart(ju0+1,i)*rpart(jfiu+1,i) +
     >                     rpart(ju0+2,i)*rpart(jfiu+2,i)

            rdum = 0.
            rdum = rdum + rpart(jqqs,i)

            rpart(jq0,i) = rdum

         enddo
c     other ---------------------------------------------------------- 
      elseif (time_integ.eq.2) then 

      endif

      ! end timer
      pttime(12) = pttime(12) + dnekclock() - ptdum(12)

      return
      end
c-----------------------------------------------------------------------
      subroutine calc_substantial_derivative
c 
c     calculate rhs of NS, which is just pressure gradient. used for
c     undisturbed and invisicid unsteady force
c
c     no forcing included...should it be?
c
c     rhs_fluidp(i,j,k,e,l)
c   
c        l = 1     >    dP/dx
c        l = 2     >    dP/dy
c        l = 3     >    dP/dz
c        l = 4     >    -P* div(phi_p v), v is eulerian particle vel
c        l = 5     >    P * d phi_g/dx
c        l = 6     >    P * d phi_g/dy
c        l = 7     >    P * d phi_g/dz
c
      include 'SIZE'
      include 'INPUT'
      include 'GEOM'
      include 'SOLN'
      include 'CMTDATA'
      include 'CMTPART'

      integer e
      real ur(lx1,ly1,lz1),us(lx1,ly1,lz1),ut(lx1,ly1,lz1),
     >         prsm(lx1,ly1,lz1,lelt)

      ! begin timer
      ptdum(13) = dnekclock()

      nxyz=nx1*ny1*nz1
      nlxyze = lx1*ly1*lz1*lelt

      call rzero(rhs_fluidp,nlxyze*7)

c     compute grad pr
      do e=1,nelt
        call gradl_rst(ur(1,1,1),us(1,1,1),ut(1,1,1),
     >                                        pr(1,1,1,e),lx1,if3d)
        if(if3d) then ! 3d
            do i=1,nxyz
              rhs_fluidp(i,1,1,e,1) = 1.0d+0/JACM1(i,1,1,e)* !d/dx
     >             (ur(i,1,1)*RXM1(i,1,1,e) +
     >              us(i,1,1)*SXM1(i,1,1,e) +
     >              ut(i,1,1)*TXM1(i,1,1,e))
            enddo

            do i=1,nxyz
              rhs_fluidp(i,1,1,e,2) = 1.0d+0/JACM1(i,1,1,e)* !d/dy
     >             (ur(i,1,1)*RYM1(i,1,1,e) +
     >              us(i,1,1)*SYM1(i,1,1,e) +
     >              ut(i,1,1)*TYM1(i,1,1,e))
            enddo

            do i=1,nxyz
              rhs_fluidp(i,1,1,e,3) = 1.0d+0/JACM1(i,1,1,e)* !d/dz
     >             (ur(i,1,1)*RZM1(i,1,1,e) +
     >              us(i,1,1)*SZM1(i,1,1,e) +
     >              ut(i,1,1)*TZM1(i,1,1,e))
            enddo
        endif ! end 3d
      enddo

      ! div (phi_p * v)
      do e=1,nelt
        call gradl_rst(ur(1,1,1),us(1,1,1),ut(1,1,1), ! x dir
     >                                        ptw(1,1,1,e,6),lx1,if3d)
        if(if3d) then ! 3d
            do i=1,nxyz
              rhs_fluidp(i,1,1,e,4) = 1.0d+0/JACM1(i,1,1,e)* !d/dx
     >             (ur(i,1,1)*RXM1(i,1,1,e) +
     >              us(i,1,1)*SXM1(i,1,1,e) +
     >              ut(i,1,1)*TXM1(i,1,1,e))
            enddo
        endif ! end 3d

        call gradl_rst(ur(1,1,1),us(1,1,1),ut(1,1,1), ! y dir
     >                                        ptw(1,1,1,e,7),lx1,if3d)
        if(if3d) then ! 3d
            do i=1,nxyz
              rhs_fluidp(i,1,1,e,4) = rhs_fluidp(i,1,1,e,4) +
     >             1.0d+0/JACM1(i,1,1,e)* !d/dy
     >             (ur(i,1,1)*RYM1(i,1,1,e) +
     >              us(i,1,1)*SYM1(i,1,1,e) +
     >              ut(i,1,1)*TYM1(i,1,1,e))
            enddo
         endif

        call gradl_rst(ur(1,1,1),us(1,1,1),ut(1,1,1), ! z dir
     >                                        ptw(1,1,1,e,8),lx1,if3d)
        if(if3d) then ! 3d
            do i=1,nxyz
              rhs_fluidp(i,1,1,e,4) = rhs_fluidp(i,1,1,e,4) +
     >             1.0d+0/JACM1(i,1,1,e)* !d/dz
     >             (ur(i,1,1)*RZM1(i,1,1,e) +
     >              us(i,1,1)*SZM1(i,1,1,e) +
     >              ut(i,1,1)*TZM1(i,1,1,e))
            enddo
        endif
      enddo

      do e=1,nelt
      do k=1,nz1
      do j=1,ny1
      do i=1,nx1
         rhs_fluidp(i,j,k,e,4) = -pr(i,j,k,e)*rhs_fluidp(i,j,k,e,4)
      enddo
      enddo
      enddo
      enddo

c     compute grad phi_g
      do e=1,nelt
        call gradl_rst(ur(1,1,1),us(1,1,1),ut(1,1,1),
     >                                        phig(1,1,1,e),lx1,if3d)
        if(if3d) then ! 3d
            do i=1,nxyz
              rhs_fluidp(i,1,1,e,5) = 1.0d+0/JACM1(i,1,1,e)*
     >             (ur(i,1,1)*RXM1(i,1,1,e) +
     >              us(i,1,1)*SXM1(i,1,1,e) +
     >              ut(i,1,1)*TXM1(i,1,1,e))
            enddo

            do i=1,nxyz
              rhs_fluidp(i,1,1,e,6) = 1.0d+0/JACM1(i,1,1,e)*
     >             (ur(i,1,1)*RYM1(i,1,1,e) +
     >              us(i,1,1)*SYM1(i,1,1,e) +
     >              ut(i,1,1)*TYM1(i,1,1,e))
            enddo

            do i=1,nxyz
              rhs_fluidp(i,1,1,e,7) = 1.0d+0/JACM1(i,1,1,e)*
     >             (ur(i,1,1)*RZM1(i,1,1,e) +
     >              us(i,1,1)*SZM1(i,1,1,e) +
     >              ut(i,1,1)*TZM1(i,1,1,e))
            enddo
        endif ! end 3d
      enddo

      do e=1,nelt
      do k=1,nz1
      do j=1,ny1
      do i=1,nx1
         rhs_fluidp(i,j,k,e,5) = pr(i,j,k,e)*rhs_fluidp(i,j,k,e,5)
         rhs_fluidp(i,j,k,e,6) = pr(i,j,k,e)*rhs_fluidp(i,j,k,e,6)
         rhs_fluidp(i,j,k,e,7) = pr(i,j,k,e)*rhs_fluidp(i,j,k,e,7)
      enddo
      enddo
      enddo
      enddo

      ! end timer
      pttime(13) = pttime(13) + dnekclock() - ptdum(13)

      return
      end
c-----------------------------------------------------------------------
      subroutine usr_particles_q_uu(ii)
c
c
      include 'SIZE'
      include 'TOTAL'
      include 'CMTDATA'
      include 'CMTPART'

      if (part_force(5) .eq. 0) then
         rpart(jquu,ii) = rpart(jvol,ii)*0. ! none for now (del . q_f) = 0
      endif

      ! end timer

      return
      end
c-----------------------------------------------------------------------
      subroutine usr_particles_q_qs(ii)
c
c
      include 'SIZE'
      include 'TOTAL'
      include 'CMTDATA'
      include 'CMTPART'

      real nu_g,S_qs,rpra

      ! begin timer
      pi      = 4.0d+0*atan(1.0d+0)

      nu_g    = mu_0/rpart(jrho,ii) ! kinematic viscosity
      rpra    = nu_g/kappa_g
      S_qs    = 2.*pi*kappa_g*rpart(jdp,ii)

      if (part_force(6) .gt. 0) then
         S_qs    = S_qs*(1. + 0.3*sqrt(rpart(jre,ii))*rpra**(2./3.)) !re < 500
      elseif (part_force(6) .eq. 0) then
         S_qs    = 0.
      endif

      rpart(jqqs,ii) = S_qs*(rpart(jtempf,ii) - rpart(jtemp,ii))

      ! end timer

      return
      end
c-----------------------------------------------------------------------
      subroutine usr_particles_f_un(ii,jj)
c
c     extra body forces
c
      include 'SIZE'
      include 'TOTAL'
      include 'CMTDATA'
      include 'CMTPART'

      if (part_force(3) .ne. 0) then
         rpart(jfun+jj,ii) = -rpart(jvol,ii)*rpart(jDuDt+jj,ii)
      else
         rpart(jfun+jj,ii) = 0.
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine usr_particles_f_iu(ii,jj)
c
c     extra body forces
c
      include 'SIZE'
      include 'TOTAL'
      include 'CMTDATA'
      include 'CMTPART'

      real S_qs

      S_qs = 0.5
      if (part_force(4) .gt. 0) then
         S_qs = S_qs*(1. + 1.8*rpart(ja,ii)**2 + 7.6*rpart(ja,ii)**4)
         rphip = rpart(jvol1,ii)
         if (rphip .gt. 0.3) rphip = 0.3
         S_qs = S_qs*(1.+2.*rphip)/(1.-rphip)
      elseif (part_force(4) .eq. 0) then
         S_qs = 0.
      endif

      rpart(jcmiu,ii) = S_qs

c     need to fix, fake, no D(rho_f u)/Dt contribution
      rpart(jfiu+jj,ii) = -1.*0.*S_qs*rpart(jvol,ii)*
     >                    rpart(jDuDt+jj,ii)

      return
      end
c-----------------------------------------------------------------------
      subroutine usr_particles_f_qs(i,j)
c     calculate quasi steady force with drag corrections
      include 'SIZE'
      include 'TOTAL'
      include 'PERFECTGAS'
      include 'CMTDATA'
      include 'CMTPART'

      real cd,S_qs

c     cd_std = 24/re_p already taken into account below
      S_qs = rpart(jvol,i)*rpart(jrhop,i)/rpart(jtaup,i)

      if (part_force(2).gt.0) then
         rrep = rpart(jre,i)
         rrma = rpart(ja,i)
         rmacr= 0.6
         rcd_std = 1.+0.15*rrep**(0.687) + 
     >               0.42*(1.+42500./rrep**(1.16))**(-1)
         rcd_mcr = 1.+0.15*rrep**(0.684) + 
     >               0.513*(1. + 483./rrep**(0.669))**(-1)
         rcd1 = rcd_std + (rcd_mcr - rcd_std)*rrma/rmacr
         
         S_qs = S_qs*rcd1

         rphip = rpart(jvol1,i)
         if (rphip .gt. 0.3) rphip = 0.3
         rcd2 = (1. - 2.*rphip)/(1. - rphip)**3

         S_qs = S_qs*rcd2

      elseif (part_force(2).eq.0) then
         S_qs = 0.
      endif

      rpart(jcd+j,i)  = 0.
      rpart(jfqs+j,i) = S_qs*(rpart(ju0+j,i) - rpart(jv0+j,i))

      return
      end
c-----------------------------------------------------------------------
      subroutine usr_particles_f_user(ii,jj)
c
c     extra body forces
c
      include 'SIZE'
      include 'TOTAL'
      include 'CMTDATA'
      include 'CMTPART'

      real pmass,rxdum,rydum

      if (part_force(1) .eq. 0) then
         rpart(jfusr+jj,ii) = 0.
      else
         pmass = rpart(jvol,ii)*rpart(jrhop,ii)
         if (jj.eq.0) rpart(jfusr+jj,ii) = 0.0 !-9.8*pmass
         if (jj.eq.1) rpart(jfusr+jj,ii) = 0.0
         if (jj.eq.2) rpart(jfusr+jj,ii) = 0.0
      endif

      return
      end
c-----------------------------------------------------------------------
