c-----------------------------------------------------------------------
      subroutine plan4_acc_data_copyin_onetime()
c-----------------------------------------------------------------------
      include 'SIZE'
      include 'TOTAL'    

      parameter (lg=lx1*ly1*lz1*lelt)
      parameter (maxcg=900)

      common /tdarray/ diagt(maxcg),upper(maxcg)
      common /scrcg/ d(lg), scalar(2)
      common /scrcg2/ r(lg), w(lg), p(lg), z(lg)

!$acc enter data copyin(xm1,ym1,zm1)
!$acc enter data copyin(dxm1,dxtm1,w3m1)
!$acc enter data copyin(bm1,bm1lag,binvm1,bintm1)
!$acc enter data copyin(jacm1,jacmi)
!$acc enter data copyin(rxm1,sxm1,txm1)
!$acc enter data copyin(rym1,sym1,tym1)
!$acc enter data copyin(rzm1,szm1,tzm1)
!$acc enter data copyin(rxm2,sxm2,txm2)
!$acc enter data copyin(rym2,sym2,tym2)
!$acc enter data copyin(rzm2,szm2,tzm2)
!$acc enter data copyin(g1m1,g2m1,g3m1,g4m1,g5m1,g6m1)
!$acc enter data copyin(v1mask,v2mask,v3mask,pmask,tmask,omask)
!$acc enter data copyin(b1mask,b2mask,b3mask,bpmask)
!$acc enter data copyin(tmult,vmult)
!$acc enter data copyin(unx,uny,unz,area)
!$acc enter data copyin(cbc,bc)
!$acc enter data copyin(param,nelfld)

      return
      end
c-----------------------------------------------------------------------
      subroutine plan4_acc_data_copyout_onetime()
c-----------------------------------------------------------------------
      include 'SIZE'
      include 'TOTAL'

      parameter (lg=lx1*ly1*lz1*lelt)
      parameter (maxcg=900)

      common /tdarray/ diagt(maxcg),upper(maxcg)
      common /scrcg/ d(lg), scalar(2)
      common /scrcg2/ r(lg), w(lg), p(lg), z(lg)

!$acc exit data copyout(xm1,ym1,zm1)
!$acc exit data copyout(dxm1,dxtm1,w3m1)
!$acc exit data copyout(bm1,bm1lag,binvm1,bintm1)
!$acc exit data copyout(jacm1,jacmi)
!$acc exit data copyout(rxm1,sxm1,txm1)
!$acc exit data copyout(rym1,sym1,tym1)
!$acc exit data copyout(rzm1,szm1,tzm1)
!$acc exit data copyout(rxm2,sxm2,txm2)
!$acc exit data copyout(rym2,sym2,tym2)
!$acc exit data copyout(rzm2,szm2,tzm2)
!$acc exit data copyout(g1m1,g2m1,g3m1,g4m1,g5m1,g6m1)
!$acc exit data copyout(v1mask,v2mask,v3mask,pmask,tmask,omask)
!$acc exit data copyout(b1mask,b2mask,b3mask,bpmask)
!$acc exit data copyout(tmult,vmult)
!$acc exit data copyout(unx,uny,unz,area)
!$acc exit data copyout(cbc,bc)
!$acc enter data copyin(param,nelfld)

      return
      end
c-----------------------------------------------------------------------
      subroutine hsmg_acc_data_copyin_onetime()
c-----------------------------------------------------------------------
      include 'SIZE'
      include 'HSMG'

!$acc enter data copyin(mg_nx)
!$acc enter data copyin(mg_ny,mg_nz)
!$acc enter data copyin(mg_nh,mg_nhz)
!$acc enter data copyin(mg_gsh_schwarz_handle)
!$acc enter data copyin(mg_gsh_handle)
!$acc enter data copyin(mg_rstr_wt_index)
!$acc enter data copyin(mg_mask_index)
!$acc enter data copyin(mg_solve_index)
!$acc enter data copyin(mg_fast_s_index)
!$acc enter data copyin(mg_fast_d_index)
!$acc enter data copyin(mg_schwarz_wt_index)
!$acc enter data copyin(mg_g_index)

!$acc enter data copyin(mg_jh)
!$acc enter data copyin(mg_jht)
!$acc enter data copyin(mg_jhfc )
!$acc enter data copyin(mg_jhfct)
!$acc enter data copyin(mg_ah)
!$acc enter data copyin(mg_bh)
!$acc enter data copyin(mg_dh)
!$acc enter data copyin(mg_dht)
!$acc enter data copyin(mg_zh)
!$acc enter data copyin(mg_rstr_wt)
!$acc enter data copyin(mg_mask)
!$acc enter data copyin(mg_fast_s)
!$acc enter data copyin(mg_fast_d)
!$acc enter data copyin(mg_schwarz_wt)
!$acc enter data copyin(mg_solve_e)
!$acc enter data copyin(mg_solve_r)
!$acc enter data copyin(mg_h1)
!$acc enter data copyin(mg_h2)
!$acc enter data copyin(mg_b)
!$acc enter data copyin(mg_g)
!$acc enter data copyin(mg_work)
!$acc enter data copyin(mg_work2)
!$acc enter data copyin(mg_worke)
!$acc enter data copyin(mg_imask)
!$acc enter data copyin(mg_h1_n)
!$acc enter data copyin(p_mg_h1)
!$acc enter data copyin(p_mg_b)
!$acc enter data copyin(p_mg_msk)

      return
      end
c-----------------------------------------------------------------------
      subroutine plan4_acc_data_copyin()
c-----------------------------------------------------------------------
      include 'SIZE'
      include 'TOTAL'

      parameter (lg=lx1*ly1*lz1*lelt)
      parameter (maxcg=900)

      common /tdarray/ diagt(maxcg),upper(maxcg)
      common /scrcg/ d(lg), scalar(2)
      common /scrcg2/ r(lg), w(lg), p(lg), z(lg)

!$acc enter data copyin(vxlag,vylag,vzlag,tlag,vgradt1,vgradt2)
!$acc enter data copyin(abx1,aby1,abz1,abx2,aby2,abz2,vdiff_e)
!$acc enter data copyin(vtrans,vdiff,bfx,bfy,bfz,cflf,c_vx,fw)
!$acc enter data copyin(bmnv,bmass,bdivw,bx,by,bz,pm,bmx,bmy,bmz)
!$acc enter data copyin(vx,vy,vz,pr,t,vx_e,vy_e,vz_e)
!$acc enter data copyin(bbx1,bby1,bbz1,bbx2,bby2,bbz2,bxlag,bylag,bzlag)

!$acc enter data copyin(ab,bd)
!$acc enter data copyin(pr,pmlag,prlag,qtl,usrdiv)
!$acc enter data copyin(vxd,vyd,vzd)
!$acc enter data copyin(diagt,upper)
!$acc enter data copyin(d,scalar,r,w,p,z)

!$acc enter data create (ibc_acc)
!$acc enter data copyin (c_vx)

      return
      end
c-----------------------------------------------------------------------
      subroutine hmh_gmres_acc_data_copyin()
c-----------------------------------------------------------------------
#ifdef _OPENACC
      include 'SIZE'
      include 'HSMG'            ! Same array space as HSMG
      include 'DXYZ'
      include 'GEOM'
      include 'GMRES'
      include 'INPUT'
      include 'MASS'
      include 'SOLN'
      include 'TSTEP'
      include 'CTIMER'
      include 'PARALLEL'

      parameter (lt=lx1*ly1*lz1*lelt)
      common /scrmg/ e(2*lt),w(lt),r(lt)
      parameter (lwk=(lx1+2)*(ly1+2)*(lz1+2))
      common /hsmgw/ work(0:lwk-1),work2(0:lwk-1)
      common /ctmp0/ w1   (lx1,ly1,lz1,lelt)
     $             , w2   (lx1,ly1,lz1,lelt)

!$acc enter data copyin(work,work2)
!$acc enter data copyin(mg_mask,mg_imask,pmask)
!$acc enter data copyin(mg_jht,mg_jh,mg_rstr_wt,mg_schwarz_wt)
!$acc enter data copyin(mg_work,mg_fast_s,mg_fast_d)
!$acc enter data copyin(h_gmres,w_gmres,v_gmres,z_gmres)
!$acc enter data copyin(c_gmres,s_gmres,x_gmres,gamma_gmres)
!$acc enter data copyin(r_gmres)
!$acc enter data copyin(ml_gmres,mu_gmres)

!$ACC ENTER DATA CREATE(e,w,r)
!$ACC ENTER DATA CREATE(w1,w2)
!$ACC ENTER DATA CREATE(wk1,wk2)

#endif
      return
      end
c-----------------------------------------------------------------------
      subroutine hmh_gmres_acc_data_copyout()
c-----------------------------------------------------------------------
#ifdef _OPENACC
      include 'SIZE'
      include 'HSMG'            ! Same array space as HSMG
      include 'DXYZ'
      include 'GEOM'
      include 'GMRES'
      include 'INPUT'
      include 'MASS'
      include 'SOLN'
      include 'TSTEP'
      include 'CTIMER'
      include 'PARALLEL'

      parameter (lt=lx1*ly1*lz1*lelt)
      common /scrmg/ e(2*lt),w(lt),r(lt)
      parameter (lwk=(lx1+2)*(ly1+2)*(lz1+2))
      common /hsmgw/ work(0:lwk-1),work2(0:lwk-1)
      common /ctmp0/ w1   (lx1,ly1,lz1,lelt)
     $ ,             w2   (lx1,ly1,lz1,lelt)

!$ACC EXIT DATA DELETE(work,work2)
!$ACC EXIT DATA DELETE(mg_mask,mg_imask,pmask)
!$ACC EXIT DATA DELETE(mg_jht,mg_jh,mg_rstr_wt,mg_schwarz_wt)
!$ACC EXIT DATA DELETE(mg_work,mg_fast_s,mg_fast_d)
!$ACC EXIT DATA COPYOUT(h_gmres,w_gmres,v_gmres,z_gmres)
!$ACC EXIT DATA COPYOUT(c_gmres,s_gmres,x_gmres,gamma_gmres)
!$ACC EXIT DATA COPYOUT(r_gmres)
!$ACC EXIT DATA COPYOUT(ml_gmres,mu_gmres)
!$ACC EXIT DATA DELETE(w1,w2)
!$ACC EXIT DATA DELETE(e,w,r)

#endif
      return
      end
c-----------------------------------------------------------------------
      subroutine print_acc(a,n,text,mode)
c-----------------------------------------------------------------------
      include 'SIZE'
      integer n,i,mode
      real    a(n)
      character*4 text

!$ACC update host(a) if (mode.eq.1)  !copies the data from device to host
      if (nid.eq.0) then
!         do i=1,n
!            write(6,10) text,i,a(i)
!         enddo
! 10      format(a4,i300,1p1e15.4)
          write(6,102) text, (a(i), i=1,4)
      endif
  102    format(a4,2x,1p4e15.4)

      return
      end
c-----------------------------------------------------------------------
      subroutine print_acc_int(a,n,text,mode)
c-----------------------------------------------------------------------
      include 'SIZE'
      integer n,i,mode
      integer a(n)
      character*4 text

#ifdef _OPENACC
!$ACC update host(a) if (mode.eq.1)  !copies the data from device to !host
#endif
      if (nid.eq.0) write(6,*) text,(a(i),i=1,n) ! prints data on host

      return
      end
c-----------------------------------------------------------------------
      subroutine global_grad3(d,u,u1,u2,u3)
c-----------------------------------------------------------------------

      include 'SIZE'
      integer i,j,k,l,e
      real d (lx1,lx1)
      real u (lx1,ly1,lz1,lelt)
      real u1(lx1,ly1,lz1,lelt)
      real u2(lx1,ly1,lz1,lelt)
      real u3(lx1,ly1,lz1,lelt)
      real tmpu1,tmpu2,tmpu3

!$ACC DATA PRESENT (d(nx1,nx1))
!$ACC&     PRESENT (u(nx1,ny1,nz1,nelt))
!$ACC&     PRESENT (u1(nx1,ny1,nz1,nelt),u2(nx1,ny1,nz1,nelt))
!$ACC&     PRESENT (u3(nx1,ny1,nz1,nelt))

!$ACC PARALLEL LOOP COLLAPSE(4) GANG WORKER VECTOR
!$ACC&    PRIVATE(tmpu1,tmpu2,tmpu3)
!dir$ NOBLOCKING
      do e=1,nelt
         do k=1,nz1
            do j=1,ny1
               do i=1,nx1
                  tmpu1 = 0.0
                  tmpu2 = 0.0
                  tmpu3 = 0.0
!$ACC LOOP SEQ
                  do l=1,nx1
                     tmpu1 = tmpu1 + d(i,l)*u(l,j,k,e)
                     tmpu2 = tmpu2 + d(j,l)*u(i,l,k,e)
                     tmpu3 = tmpu3 + d(k,l)*u(i,j,l,e)
                  enddo
                  u1(i,j,k,e) = tmpu1
                  u2(i,j,k,e) = tmpu2
                  u3(i,j,k,e) = tmpu3
               enddo
            enddo
         enddo
      enddo
!$ACC END PARALLEL LOOP

!$ACC END DATA

      return
      end
c-----------------------------------------------------------------------
      subroutine global_div3(d,u1,u2,u3,v1,v2,v3)
c-----------------------------------------------------------------------

      include 'SIZE'
      integer i,j,k,l,e
      real d (lx1,ly1)
      real u1(lx1,ly1,lz1,lelt)
      real u2(lx1,ly1,lz1,lelt)
      real u3(lx1,ly1,lz1,lelt)
      real v1(lx1,ly1,lz1,lelt)
      real v2(lx1,ly1,lz1,lelt)
      real v3(lx1,ly1,lz1,lelt)
      real tmpu1,tmpu2,tmpu3

!$acc update host(d,u1,u2,u3,v1)
c     call outpost(u1,u2,u3,d,v1,'wgd')
c     stop

!$acc data present (d (nx1,nx1))
!$acc&     present (u1(nx1,ny1,nz1,nelv))
!$acc&     present (u2(nx1,ny1,nz1,nelv),u3(nx1,ny1,nz1,nelv))
!$acc$     present (v1(nx1,ny1,nz1,nelv),v2(nx1,ny1,nz1,nelv))
!$acc&     present (v3(nx1,ny1,nz1,nelv))

!$acc parallel loop collapse(4) gang worker vector
!$acc&    private(tmpu1,tmpu2,tmpu3)

cc!dir$ noblocking
c!$acc kernels
      do e=1,nelv
      do k=1,nz1
      do j=1,ny1
      do i=1,nx1
         tmpu1 = 0.0
         tmpu2 = 0.0
         tmpu3 = 0.0
!$acc    loop seq
         do l=1,nx1
            tmpu1 = tmpu1 + d(i,l)*u1(l,j,k,e)
            tmpu2 = tmpu2 + d(j,l)*u2(i,l,k,e)
            tmpu3 = tmpu3 + d(k,l)*u3(i,j,l,e)
         enddo
!$acc    end loop
         v1(i,j,k,e) = tmpu1
         v2(i,j,k,e) = tmpu2
         v3(i,j,k,e) = tmpu3
      enddo
      enddo
      enddo
      enddo
c!$acc end kernels
!$acc end parallel loop

!$acc end data
      end
c-----------------------------------------------------------------------
      subroutine cresvsp_acc (resv1,resv2,resv3,h1,h2)
      include 'SIZE'
      include 'TOTAL'

      real resv1(lx1*ly1*lz1*lelv)
     $   , resv2(lx1*ly1*lz1*lelv)
     $   , resv3(lx1*ly1*lz1*lelv)
     $   , h1   (lx1*ly1*lz1*lelv)
     $   , h2   (lx1*ly1*lz1*lelv)

      common /scruz/ ta1   (lx1*ly1*lz1*lelv)
     $ ,             ta2   (lx1*ly1*lz1*lelv)
     $ ,             ta3   (lx1*ly1*lz1*lelv)
     $ ,             ta4   (lx1*ly1*lz1*lelv)

      n = lx1*ly1*lz1*nelv
      intype = -1

      scale = -1./3.
      if (ifstrs) scale =  2./3.

      dtbd = bd(1)/dt

!$acc data present(resv1,resv2,resv3,h1,h2) create(ta1,ta2,ta3,ta4)
!$acc& present(vdiff,qtl,pr,vtrans,vx,vy,vz,bfx,bfy,bfz)

!$acc parallel loop 
      do i=1,n
         ta4(i)=scale*(vdiff(i,1,1,1,1)*qtl(i,1,1,1))+pr(i,1,1,1)
         h1 (i)=vdiff (i,1,1,1,1)
         h2 (i)=vtrans(i,1,1,1,1)*dtbd
      enddo
!$acc end parallel loop 

      call ophx_acc(resv1,resv2,resv3,vx,vy,vz,h1,h2)
c!$acc update host(resv1,resv2,resv3)

      do i=1,lx1*ly1*lz1
c        write (6,*) 'resv1=',resv1(i)
c        write (6,*) 'resv2=',resv2(i)
c        write (6,*) 'resv3=',resv3(i)
      enddo
c     stop

!$acc update host(ta4)
      do i=1,lx1*ly1*lz1*nelv
c        write (6,*) 'ta4=',ta4(i)
      enddo
c     stop

      call wgradm1_acc(ta1,ta2,ta3,ta4,nelv)

c!$acc update host(ta1,ta2,ta3)
      do i=1,lx1*ly1*lz1*nelv
c        write (6,*) 'ta1=',ta1(i)
c        write (6,*) 'ta2=',ta2(i)
c        write (6,*) 'ta3=',ta3(i)
      enddo
c     stop

c!$acc update host(resv1,resv2,resv3)
      do i=1,lx1*ly1*lz1
c        write (6,*) 'resv1=',resv1(i)
c        write (6,*) 'resv2=',resv2(i)
c        write (6,*) 'resv3=',resv3(i)
      enddo
c     stop

!$acc parallel loop
      do i=1,n
         resv1(i)=bfx(i,1,1,1)-resv1(i)-ta1(i)
         resv2(i)=bfy(i,1,1,1)-resv2(i)-ta2(i)
         resv3(i)=bfz(i,1,1,1)-resv3(i)-ta3(i)
      enddo
!$acc end parallel loop 

c!$acc update host(resv1,resv2,resv3)
      do i=1,lx1*ly1*lz1*nelv
c        write (6,*) 'resv1=',resv1(i)
c        write (6,*) 'resv2=',resv2(i)
c        write (6,*) 'resv3=',resv3(i)
      enddo
c     stop

!$acc end data

      return
      end
c-----------------------------------------------------------------------
      subroutine wgradm1_acc(ux,uy,uz,u,nel) ! weak form of grad 
c
c     Compute gradient -- mesh 1 to mesh 1 (vel. to vel.)
c
      include 'SIZE'
      include 'DXYZ'
      include 'GEOM'
      include 'INPUT'
      include 'TSTEP'
      include 'WZ'

      parameter (lxyz=lx1*ly1*lz1)
      real ux(lxyz,lelv),uy(lxyz,lelv),uz(lxyz,lelv),u(lx1,ly1,lz1,lelv)

      real ur(lx1,ly1,lz1),us(lx1,ly1,lz1),ut(lx1,ly1,lz1)

      integer e,q

!$acc data present(rxm1,sxm1,txm1,rym1,sym1,tym1,rzm1,szm1,tzm1)
!$acc& present(dxm1,w3m1,ux,uy,uz,u)

!$acc parallel loop gang private(ur,us,ut)
      do e=1,nel
!$acc    loop collapse(3) vector private(w1,w2,w3) 
         do k=1,lz1
         do j=1,ly1
         do i=1,lx1
            w1 = 0.0
            w2 = 0.0
            w3 = 0.0
!$acc       loop seq
            do q=1,lx1
               w1 = w1 + dxm1(i,q)*u(q,j,k,e) ! D u
               w2 = w2 + dxm1(j,q)*u(i,q,k,e) ! D u
               w3 = w3 + dxm1(k,q)*u(i,j,q,e) ! D u
            enddo
!$acc       end loop
            ur(i,j,k)=w1
            us(i,j,k)=w2
            ut(i,j,k)=w3
         enddo
         enddo
         enddo
!$acc    end loop

!$acc    loop vector
         do i=1,lxyz
            ux(i,e)=w3m1(i,1,1)*(ur(i,1,1)*rxm1(i,1,1,e)
     $                         + us(i,1,1)*sxm1(i,1,1,e)
     $                         + ut(i,1,1)*txm1(i,1,1,e))
            uy(i,e)=w3m1(i,1,1)*(ur(i,1,1)*rym1(i,1,1,e)
     $                         + us(i,1,1)*sym1(i,1,1,e)
     $                         + ut(i,1,1)*tym1(i,1,1,e))
            uz(i,e)=w3m1(i,1,1)*(ur(i,1,1)*rzm1(i,1,1,e)
     $                         + us(i,1,1)*szm1(i,1,1,e)
     $                         + ut(i,1,1)*tzm1(i,1,1,e))
         enddo
!$acc    end loop
      enddo
!$acc end parallel loop

!$acc end data

      return
      end
c-----------------------------------------------------------------------
      subroutine ophx_acc(out1,out2,out3,inp1,inp2,inp3,h1,h2)

c     OUT = (H1*A+H2*B) * INP  

      include 'SIZE'
c     include 'INPUT'
c     include 'SOLN'
      
      parameter (lt=lx1*ly1*lz1*lelt)

      real out1(lt),out2(lt),out3(lt),inp1(lt),inp2(lt),inp3(lt),
     $     h1(lt),h2(lt)

      imesh = 1

c     IF (IFSTRS) THEN
c        MATMOD = 0
c        CALL AXHMSF (OUT1,OUT2,OUT3,INP1,INP2,INP3,H1,H2,MATMOD)
c     ELSE

c     Later - we can make this fast


!$acc data present(inp1,inp2,inp3,h1,h2)
      call axhelm_acc(out1,inp1,h1,h2,imesh,1)
      call axhelm_acc(out2,inp2,h1,h2,imesh,2)
      call axhelm_acc(out3,inp3,h1,h2,imesh,3)
!$acc end data

c!$acc update host(out3)

      do i=1,lx1*ly1*lz1*nelv
c        write (6,*) 'out3=',out3(i)
      enddo
c     stop

      return
      end
c-----------------------------------------------------------------------
      subroutine makef_acc
c
c     Compute and add: (1) user specified forcing function (FX,FY,FZ)
c                      (2) driving force due to natural convection
c                      (3) convection term
c
c     !! NOTE: Do not change the arrays BFX, BFY, BFZ until the
c              current time step is completed.
c
      include 'SIZE'
      include 'SOLN'
      include 'MASS'
      include 'INPUT'
      include 'TSTEP'

      include 'MVGEOM'

      call chck('abb')
      call makeuf_acc    ! paul and som
      call chck('bbb')
      if (ifnav .and..not.ifchar) call advab_acc
      call chck('cbb')
      call makeabf_acc
      call chck('dbb')
      call makebdf_acc

      call chck('ebb')

      return
      end
C-----------------------------------------------------------------------
      subroutine makeuf_acc
C---------------------------------------------------------------------
C
C     Compute and add: (1) user specified forcing function (FX,FY,FZ)
C     for the gpu version of the code, FX=FY=FZ=0  
C
C----------------------------------------------------------------------
      include 'SIZE'
      include 'SOLN'
      include 'MASS'
      include 'TSTEP'

      TIME = TIME-DT

!$acc update host      (bfx,bfy,bfz)
      call nekuf       (bfx,bfy,bfz)
!$acc update device    (bfx,bfy,bfz)

      CALL opcolv_acc (BFX,BFY,BFZ,BM1)

      TIME = TIME+DT

      return
      END
c-----------------------------------------------------------------------
      subroutine oprzero_acc (a,b,c)
      include 'SIZE'
      real a(1),b(1),c(1)

      ntot1=nx1*ny1*nz1*nelv

      call rzero_acc(a,ntot1)
      call rzero_acc(b,ntot1)
      call rzero_acc(c,ntot1)

      return
      end
c------------------------------------------------------------------------
      subroutine opcolv_acc (a1,a2,a3,c)
      include 'SIZE'
      real a1(nx1*ny1*nz1*nelv),
     $     a2(nx1*ny1*nz1*nelv),
     $     a3(nx1*ny1*nz1*nelv),
     $      c(nx1*ny1*nz1*nelv)

      ntot1=nx1*ny1*nz1*nelv

      call col2_acc(a1,c,ntot1)
      call col2_acc(a2,c,ntot1)
      call col2_acc(a3,c,ntot1)

      return
      end

c------------------------------------------------------------------------
      subroutine advab_acc
C---------------------------------------------------------------
C
C     Eulerian scheme, add convection term to forcing function 
C     at current time step.
C
C---------------------------------------------------------------
      include 'SIZE'
      include 'SOLN'
      include 'MASS'
      include 'TSTEP'
      include 'DEALIAS'

      real ta1 (lx1*ly1*lz1*lelv)
     $   , ta2 (lx1*ly1*lz1*lelv)
     $   , ta3 (lx1*ly1*lz1*lelv)

!$acc data create(ta1,ta2,ta3) present(vx,vy,vz,bfx,bfy,bfz,bm1)

      n = lx1*ly1*lz1*nelv

c!$acc update host(vxd,vyd,vzd)
      call chck('a11')
      call convop_acc(ta1,vx)
      call chck('b11')
      call convop_acc(ta2,vy)
      call chck('c11')
      call convop_acc(ta3,vz)

c     call outpost(ta1,ta2,ta3,pr,t,'wta')

!$acc parallel loop
      do i=1,n
c        bfx(i,1,1,1)=bfx(i,1,1,1)-ta1(i)*bm1(i,1,1,1)
c        bfy(i,1,1,1)=bfy(i,1,1,1)-ta2(i)*bm1(i,1,1,1)
c        bfz(i,1,1,1)=bfz(i,1,1,1)-ta3(i)*bm1(i,1,1,1)
         bfx(i,1,1,1)=bfx(i,1,1,1)-ta1(i)
         bfy(i,1,1,1)=bfy(i,1,1,1)-ta2(i)
         bfz(i,1,1,1)=bfz(i,1,1,1)-ta3(i)
      enddo
!$acc end loop
!$acc end data 

      return
      end
c-----------------------------------------------------------------------
      subroutine get_jggl(jggl,dggl) ! acc copyin happens here
c-------------------------------------------------------------------
      include 'SIZE'
      real jggl(lxd,lx1),dggl(lxd,lx1)

      parameter (ldg=lxd**3,lwkd=4*lxd*lxd)
      common /dgrad/ d(ldg),dt(ldg),dg(ldg),dgt(ldg),jgl(ldg),jgt(ldg)
     $             , wkd(lwkd)
      real jgl,jgt

      call get_int_ptr (i,nx1,nxd)
      call copy(jggl,jgl(i),lx1*lxd)          ! jggl = jgl

      call get_dgl_ptr (ip,nxd,nxd)
      call mxm(dg(ip),nxd,jggl,nxd,dggl,nx1)  ! dggl = dg*jggl

      return
      end
c-----------------------------------------------------------------------
      subroutine makeabf_acc
C-----------------------------------------------------------------------
C
C     Sum up contributions to kth order extrapolation scheme.
C
C-----------------------------------------------------------------------
      include 'SIZE'
      include 'INPUT'
      include 'SOLN'
      include 'TSTEP'
      common /SCRUZ/ ta1 (LX1,LY1,LZ1,LELV)
     $ ,             ta2 (LX1,LY1,LZ1,LELV)
     $ ,             ta3 (LX1,LY1,LZ1,LELV)
      ntot1 = nx1*ny1*nz1*nelv

!$acc data create (ta1,ta2,ta3)

      ab0 = ab(1)
      ab1 = ab(2)
      ab2 = ab(3)

      call add3s2_acc (ta1,abx1,abx2,ab1,ab2,ntot1)
      call add3s2_acc (ta2,aby1,aby2,ab1,ab2,ntot1)
      call copy_acc (abx2,abx1,ntot1)
      call copy_acc (aby2,aby1,ntot1)
      call copy_acc (abx1,bfx,ntot1)
      call copy_acc (aby1,bfy,ntot1)
      call add2s1_acc (bfx,ta1,ab0,ntot1)
      call add2s1_acc (bfy,ta2,ab0,ntot1)

      call add3s2_acc (ta3,abz1,abz2,ab1,ab2,ntot1)
      call copy_acc (abz2,abz1,ntot1)
      call copy_acc (abz1,bfz,ntot1)
      call add2s1_acc (bfz,ta3,ab0,ntot1)

      if (.not.iflomach) then
         call col2_acc (bfx,vtrans,ntot1)
         call col2_acc (bfy,vtrans,ntot1)
         call col2_acc (bfz,vtrans,ntot1)
      endif

!$acc end data

      return
      end

c-----------------------------------------------------------------------
      subroutine makebdf_acc
c-----------------------------------------------------------------------
C     Add contributions to F from lagged BD terms.
      include 'SIZE'
      include 'SOLN'
      include 'MASS'
      include 'GEOM'
      include 'INPUT'
      include 'TSTEP'

      integer e

      real TA1(LX1,LY1,LZ1,LELV)
     $ ,   TA2(LX1,LY1,LZ1,LELV)
     $ ,   TA3(LX1,LY1,LZ1,LELV)
     $ ,   TB1(LX1,LY1,LZ1,LELV)
     $ ,   TB2(LX1,LY1,LZ1,LELV)
     $ ,   TB3(LX1,LY1,LZ1,LELV)
     $ ,   H2 (LX1,LY1,LZ1,LELV)

      NTOT1 = NX1*NY1*NZ1*NELV
      const = 1./DT
!$acc enter data create(ta1,ta2,ta3,tb1,tb2,tb3,h2)

c     INLINED:
c     call cmult2_acc(h2,vtrans(1,1,1,1,ifield),const,ntot1)
!$acc parallel loop collapse(4)
      do e = 1, nelv
      do k = 1, lz1
      do j = 1, ly1
      do i = 1, lx1
         h2(i,j,k,e)=vtrans(i,j,k,e,ifield)*const
      enddo
      enddo
      enddo
      enddo
!$acc end parallel

c     INLINED:
c     CALL opcolv3c_acc (tb1,tb2,tb3,vx,vy,vz,bm1,bd(2))

!$acc update host(vx,vy,vz,bm1)
c     call outpost(vx,vy,vz,bm1,t,'wso')

!$acc update host(bd)
c     write (6,*) 'bd2=',bd(2)
!$acc parallel loop collapse(4)
      do e = 1, nelv
      do k = 1, lz1
      do j = 1, ly1
      do i = 1, lx1
         tb1(i,j,k,e)=vx(i,j,k,e)*bm1(i,j,k,e)*bd(2)
         tb2(i,j,k,e)=vy(i,j,k,e)*bm1(i,j,k,e)*bd(2)
         tb3(i,j,k,e)=vz(i,j,k,e)*bm1(i,j,k,e)*bd(2)
      enddo
      enddo
      enddo
      enddo
!$acc end parallel

      do ilag=2,nbd

         if (ifgeom) then
c        INLINED:
c           CALL opcolv3c_acc(TA1,TA2,TA3,VXLAG (1,1,1,1,ILAG-1),
c    $                                VYLAG (1,1,1,1,ILAG-1),
c    $                                VZLAG (1,1,1,1,ILAG-1),
c    $                                BM1LAG(1,1,1,1,ILAG-1),bd(ilag+1))
!$acc parallel loop collapse(4)
            do e = 1, nelv
            do k = 1, lz1
            do j = 1, ly1
            do i = 1, lx1
               ta1(i,j,k,e) = vxlag(i,j,k,e,ilag-1)  * 
     $                        bm1lag(i,j,k,e,ilag-1) *
     $                        bd(ilag+1)
               ta2(i,j,k,e) = vylag(i,j,k,e,ilag-1)  * 
     $                        bm1lag(i,j,k,e,ilag-1) *
     $                        bd(ilag+1)
               ta3(i,j,k,e) = vzlag(i,j,k,e,ilag-1)  * 
     $                        bm1lag(i,j,k,e,ilag-1) *
     $                        bd(ilag+1)
            enddo
            enddo
            enddo
            enddo
!$acc end parallel
         else
c           INLINED:
c           CALL opcolv3c_acc(TA1,TA2,TA3,VXLAG (1,1,1,1,ILAG-1),
c    $                                VYLAG (1,1,1,1,ILAG-1),
c    $                                VZLAG (1,1,1,1,ILAG-1),
c    $                                BM1                   ,bd(ilag+1))
!$acc parallel loop collapse(4)
            do e = 1, nelv
            do k = 1, lz1
            do j = 1, ly1
            do i = 1, lx1
               ta1(i,j,k,e) = vxlag(i,j,k,e,ilag-1)  * 
c    $                        bm1(i,j,k,ilag-1+e) *
     $                        bm1(i,j,k,e) *
     $                        bd(ilag+1)
               ta2(i,j,k,e) = vylag(i,j,k,e,ilag-1)  * 
c    $                        bm1(i,j,k,ilag-1+e) *
     $                        bm1(i,j,k,e) *
     $                        bd(ilag+1)
               ta3(i,j,k,e) = vzlag(i,j,k,e,ilag-1)  * 
c    $                        bm1(i,j,k,ilag-1+e) *
     $                        bm1(i,j,k,e) *
     $                        bd(ilag+1)
            enddo
            enddo
            enddo
            enddo
!$acc end parallel
c!$acc update host(vxlag,vylag,vzlag)
c            write (6,*) 'vlag dev',istep
c            do i=1,lx1*ly1*lz1
c               write (6,*) 'vlag vxlag=',vxlag(i,1,1,1,ilag-1),ilag
c               write (6,*) 'vlag vylag=',vylag(i,1,1,1,ilag-1),ilag
c               write (6,*) 'vlag vzlag=',vzlag(i,1,1,1,ilag-1),ilag
c               write (6,*) 'vlag bm1=',bm1(i,1,1,1)
c               write (6,*) 'vlag bd=',bd(ilag+1)
c            enddo
         endif
c        INLINED:
c        call opadd2_acc(TB1,TB2,TB3,TA1,TA2,TA3)
c         =
c          CALL ADD2(TB1,TA1,NTOT1)
c          CALL ADD2(TB2,TA2,NTOT1)
c          IF(NDIM.EQ.3)CALL ADD2(TB3,TA3,NTOT1)
!$acc parallel loop collapse(4)
         do e = 1, nelv
         do k = 1, lz1
         do j = 1, ly1
         do i = 1, lx1
            tb1(i,j,k,e) = tb1(i,j,k,e) + ta1(i,j,k,e)
            tb2(i,j,k,e) = tb2(i,j,k,e) + ta2(i,j,k,e)
c           tb3(i,j,k,e) = tb3(i,j,k,e) + ta1(i,j,k,e)
            tb3(i,j,k,e) = tb3(i,j,k,e) + ta3(i,j,k,e)
         enddo
         enddo
         enddo
         enddo
!$acc end parallel
      enddo

c     INLINED:
c     call opadd2col_acc(BFX,BFY,BFZ,TB1,TB2,TB3,h2)
!$acc parallel loop collapse(4)
      do e = 1, nelv
      do k = 1, lz1
      do j = 1, ly1
      do i = 1, lx1
         bfx(i,j,k,e)=bfx(i,j,k,e)+tb1(i,j,k,e)*h2(i,j,k,e)
         bfy(i,j,k,e)=bfy(i,j,k,e)+tb2(i,j,k,e)*h2(i,j,k,e)
         bfz(i,j,k,e)=bfz(i,j,k,e)+tb3(i,j,k,e)*h2(i,j,k,e)
      enddo
      enddo
      enddo
      enddo
!$acc end parallel
!$acc exit data

      return
      end
c-----------------------------------------------------------------------
      subroutine sumab_vx_acc ! Extrapolate velocity
c
c     sum up AB/BDF contributions 
c
      include 'SIZE'
      include 'TOTAL'

      n = lx1*ly1*lz1*nelv

!$acc update device(ab)
!$acc data present(vx,vy,vz,vxlag,vylag,vzlag,vx_e,vy_e,vz_e)

      if (nab.eq.3) then

!$acc    parallel loop
         do i=1,n
            vx_e(i)=ab(1)*vx(i,1,1,1)
     $             +ab(2)*vxlag(i,1,1,1,1)
     $             +ab(3)*vxlag(i,1,1,1,2)
            vy_e(i)=ab(1)*vy(i,1,1,1)
     $             +ab(2)*vylag(i,1,1,1,1)
     $             +ab(3)*vylag(i,1,1,1,2)
            vz_e(i)=ab(1)*vz(i,1,1,1)
     $             +ab(2)*vzlag(i,1,1,1,1)
     $             +ab(3)*vzlag(i,1,1,1,2)
         enddo
!$acc    end parallel
      else
!$acc    parallel loop
         do i=1,n
            vx_e(i)=ab(1)*vx(i,1,1,1)
     $             +ab(2)*vxlag(i,1,1,1,1)
            vy_e(i)=ab(1)*vy(i,1,1,1)
     $             +ab(2)*vylag(i,1,1,1,1)
            vz_e(i)=ab(1)*vz(i,1,1,1)
     $             +ab(2)*vzlag(i,1,1,1,1)
         enddo
!$acc    end parallel
      endif

!$acc end data

      return
      end
c-----------------------------------------------------------------------
      subroutine chck(s3)

      include 'SIZE'
      character*3 s3

c     write(6,*) 'checker: ',s3

      return
      end
c-----------------------------------------------------------------------
c      subroutine opcolv3c_acc(a1,a2,a3,b1,b2,b3,c,d)
c-----------------------------------------------------------------------
c      include 'SIZE'
c      parameter (n = lx1*ly1*lz1*lelt)
c      real a1(n),a2(n),a3(n)
c      real b1(n),b2(n),b3(n)
c      real c (n)
c
c!$acc parallel loop present(a1,a2,a3,b1,b2,b3,c)
c      do i=1,n
c         a1(i)=b1(i)*c(i)*d
c         a2(i)=b2(i)*c(i)*d
c         a3(i)=b3(i)*c(i)*d
c      enddo
c!$acc end parallel 
c
c      return
c      end
c-----------------------------------------------------------------------
      subroutine opadd2_acc (a1,a2,a3,b1,b2,b3)
      include 'SIZE'
      REAL A1(lx1*ly1*lz1*lelt),
     $     A2(lx1*ly1*lz1*lelt),
     $     A3(lx1*ly1*lz1*lelt),
     $     B1(lx1*ly1*lz1*lelt),
     $     B2(lx1*ly1*lz1*lelt),
     $     B3(lx1*ly1*lz1*lelt)

      NTOT1=lx1*ly1*lz1*lelt

      CALL ADD2_ACC(A1,B1,NTOT1)
      CALL ADD2_ACC(A2,B2,NTOT1)
      IF(NDIM.EQ.3)CALL ADD2_ACC(A3,B3,NTOT1)

      return
      end
c-----------------------------------------------------------------------
      subroutine opadd2col_acc(a1,a2,a3,b1,b2,b3,c)
      include 'SIZE'
      REAL A1(lx1*ly1*lz1*lelt),
     $     A2(lx1*ly1*lz1*lelt),
     $     A3(lx1*ly1*lz1*lelt),
     $     B1(lx1*ly1*lz1*lelt),
     $     B2(lx1*ly1*lz1*lelt),
     $     B3(lx1*ly1*lz1*lelt)

      NTOT1=lx1*ly1*lz1*lelt

      IF (NDIM.EQ.3) THEN
         call add2col2_acc(a1,b1,c,ntot1)
         call add2col2_acc(a2,b2,c,ntot1)
         call add2col2_acc(a3,b3,c,ntot1)
      ELSE
         call add2col2_acc(a1,b1,c,ntot1)
         call add2col2_acc(a2,b2,c,ntot1)
      ENDIF

      return
      end
c-----------------------------------------------------------------------
      subroutine curl_acc(u1r, u1s, u1t, u2r, u2s, u2t, u3r, u3s, u3t,
     $                    rxm1,sxm1,txm1,rym1,sym1,tym1,rzm1,szm1,tzm1,
     $                    w1,  w2,  w3,  jacmi)
c-----------------------------------------------------------------------

      include 'SIZE'

      real u1r(lx1,ly1,lz1,lelt)
      real u1s(lx1,ly1,lz1,lelt)
      real u1t(lx1,ly1,lz1,lelt)
      real u2r(lx1,ly1,lz1,lelt)
      real u2s(lx1,ly1,lz1,lelt)
      real u2t(lx1,ly1,lz1,lelt)
      real u3r(lx1,ly1,lz1,lelt)
      real u3s(lx1,ly1,lz1,lelt)
      real u3t(lx1,ly1,lz1,lelt)
      real w1 (lx1,ly1,lz1,lelt)
      real w2 (lx1,ly1,lz1,lelt)
      real w3 (lx1,ly1,lz1,lelt)

      real jacmi(lx1*ly1*lz1,lelt)
      real rxm1(lx1,ly1,lz1,lelt)
      real sxm1(lx1,ly1,lz1,lelt)
      real txm1(lx1,ly1,lz1,lelt)
      real rym1(lx1,ly1,lz1,lelt)
      real sym1(lx1,ly1,lz1,lelt)
      real tym1(lx1,ly1,lz1,lelt)
      real rzm1(lx1,ly1,lz1,lelt)
      real szm1(lx1,ly1,lz1,lelt)
      real tzm1(lx1,ly1,lz1,lelt)

      integer i,j,k,l,e

      nxyz = lx1*ly1*lz1

!$ACC DATA PRESENT(u1r,u1s,u1t,u2r,u2s,u2t,u3r,u3s,u3t,w1,w2,w3)
!$ACC& PRESENT(jacmi,rxm1,sxm1,txm1,rym1,sym1,tym1,rzm1,szm1,tzm1)


!$ACC PARALLEL LOOP GANG WORKER VECTOR
         do i = 1,nxyz*nelv
            w1(i,1,1,1) = (u3r(i,1,1,1)*rym1(i,1,1,1)
     $                   + u3s(i,1,1,1)*sym1(i,1,1,1)
     $                   + u3t(i,1,1,1)*tym1(i,1,1,1)
     $                   - u2r(i,1,1,1)*rzm1(i,1,1,1)
     $                   - u2s(i,1,1,1)*szm1(i,1,1,1)
     $                   - u2t(i,1,1,1)*tzm1(i,1,1,1))*jacmi(i,1)

            w2(i,1,1,1)= (u1r(i,1,1,1)*rzm1(i,1,1,1)
     $                  + u1s(i,1,1,1)*szm1(i,1,1,1)
     $                  + u1t(i,1,1,1)*tzm1(i,1,1,1)
     $                  - u3r(i,1,1,1)*rxm1(i,1,1,1)
     $                  - u3s(i,1,1,1)*sxm1(i,1,1,1)
     $                  - u3t(i,1,1,1)*txm1(i,1,1,1))*jacmi(i,1)

            w3(i,1,1,1)= (u2r(i,1,1,1)*rxm1(i,1,1,1)
     $                  + u2s(i,1,1,1)*sxm1(i,1,1,1)
     $                  + u2t(i,1,1,1)*txm1(i,1,1,1)
     $                  - u1r(i,1,1,1)*rym1(i,1,1,1)
     $                  - u1s(i,1,1,1)*sym1(i,1,1,1)
     $                  - u1t(i,1,1,1)*tym1(i,1,1,1))*jacmi(i,1)
         enddo
!$ACC END PARALLEL LOOP

!$ACC END DATA

c     flop_a = flop_a + (18.*nxyz*lx1)*nelt



      return
      end
c---------------------------------------------------------------------------           
c-----------------------------------------------------------------------
      subroutine global_curl_grad3_acc
     $           (u1r,u1s,u1t,u2r,u2s,u2t,u3r,u3s,u3t,u1,u2,u3,d)
c-----------------------------------------------------------------------
      include 'SIZE'

      real u1r(lx1,ly1,lz1,lelt)
      real u1s(lx1,ly1,lz1,lelt)
      real u1t(lx1,ly1,lz1,lelt)
      real u2r(lx1,ly1,lz1,lelt)
      real u2s(lx1,ly1,lz1,lelt)
      real u2t(lx1,ly1,lz1,lelt)
      real u3r(lx1,ly1,lz1,lelt)
      real u3s(lx1,ly1,lz1,lelt)
      real u3t(lx1,ly1,lz1,lelt)

      real u1 (lx1,ly1,lz1,lelt)
      real u2 (lx1,ly1,lz1,lelt)
      real u3 (lx1,ly1,lz1,lelt)

      real d(lx1,ly1)

      real tmpr1,tmps1,tmpt1
      real tmpr2,tmps2,tmpt2
      real tmpr3,tmps3,tmpt3
      integer i,j,k,l,e

c     tempo = 0.0

c     do i=1,lx1
c        write (6,*) 'u3=',u3(i,1,1,1)
c        write (6,*) 'd=',d(1,i)
c        tempo=tempo+u3(i,1,1,1)*d(1,i)
c     enddo

c     write (6,*) 'u3r(1,1,1,1)=',tempo

c     call exitti('exit before mxm$',1)

!$acc update device(u3,d)

!$ACC DATA PRESENT(u1r,u1s,u1t,u2r,u2s,u2t,u3r,u3s,u3t)
!$ACC&    PRESENT(u1,u2,u3,d)

!$ACC PARALLEL LOOP COLLAPSE(4) GANG WORKER VECTOR
!$ACC&     private(tmpr1,tmpr2,tmpr3,tmps1,tmps2,tmps3,
!$ACC&             tmpt1,tmpt2,tmpt3)
!dir$ NOBLOCKING
      do e = 1,nelv
         do k = 1,nz1
            do j = 1,ny1
               do i = 1,nx1
                  tmpr1 = 0.0
                  tmpr2 = 0.0
                  tmpr3 = 0.0
                  tmps1 = 0.0
                  tmps2 = 0.0
                  tmps3 = 0.0
                  tmpt1 = 0.0
                  tmpt2 = 0.0
                  tmpt3 = 0.0
!$ACC LOOP SEQ
                  do l=1,nx1
                     tmpr1=tmpr1+d(i,l)*u1(l,j,k,e)

                     tmpr2=tmpr2+d(i,l)*u2(l,j,k,e)
                     tmpr3=tmpr3+d(i,l)*u3(l,j,k,e)

                     tmps1=tmps1+d(j,l)*u1(i,l,k,e)

                     tmps2=tmps2+d(j,l)*u2(i,l,k,e)
                     tmps3=tmps3+d(j,l)*u3(i,l,k,e)

                     tmpt1=tmpt1+d(k,l)*u1(i,j,l,e)

                     tmpt2=tmpt2+d(k,l)*u2(i,j,l,e)
                     tmpt3=tmpt3+d(k,l)*u3(i,j,l,e)
                  enddo
                  u1r(i,j,k,e) = tmpr1
                  u2r(i,j,k,e) = tmpr2
                  u3r(i,j,k,e) = tmpr3
                  u1s(i,j,k,e) = tmps1
                  u2s(i,j,k,e) = tmps2
                  u3s(i,j,k,e) = tmps3
                  u1t(i,j,k,e) = tmpt1
                  u2t(i,j,k,e) = tmpt2
                  u3t(i,j,k,e) = tmpt3
               enddo
            enddo
         enddo
      enddo
                                                                     
!$ACC END PARALLEL LOOP
!$acc update host(u3r)
!$ACC END DATA

c     write (6,*) 'u3r(1,1,1,1)=',u3r(1,1,1,1)
c     call exitti('exit after mxm$',1)

      return
      end
c-----------------------------------------------------------------------
      subroutine op_curl_acc(w1,w2,w3,u1,u2,u3)

      include 'SIZE'
      include 'TOTAL'

      real u1r(lx1,ly1,lz1,lelt)
      real u1s(lx1,ly1,lz1,lelt)
      real u1t(lx1,ly1,lz1,lelt)
      real u2r(lx1,ly1,lz1,lelt)
      real u2s(lx1,ly1,lz1,lelt)
      real u2t(lx1,ly1,lz1,lelt)
      real u3r(lx1,ly1,lz1,lelt)
      real u3s(lx1,ly1,lz1,lelt)
      real u3t(lx1,ly1,lz1,lelt)

      real w1(lx1*ly1*lz1*lelt),
     $     w2(lx1*ly1*lz1*lelt),
     $     w3(lx1*ly1*lz1*lelt),
     $     u1(lx1*ly1*lz1*lelt),
     $     u2(lx1*ly1*lz1*lelt),
     $     u3(lx1*ly1*lz1*lelt)
c

!$ACC DATA PRESENT(w1,w2,w3,u1,u2,u3)
!$ACC& PRESENT(jacmi,rxm1,sxm1,txm1,rym1,sym1,tym1,rzm1,szm1,tzm1)
!$ACC& PRESENT(dxm1,dxtm1)
!$ACC& CREATE(u1r,u1s,u1t,u2r,u2s,u2t,u3r,u3s,u3t)

      call global_curl_grad3_acc(u1r,u1s,u1t,
     $   u2r,u2s,u2t,u3r,u3s,u3t,u1,u2,u3,dxm1)

      call curl_acc (u1r,u1s,u1t,u2r,u2s,u2t,u3r,u3s,u3t,
     $               rxm1,sxm1,txm1,rym1,sym1,tym1,rzm1,szm1,tzm1,
     $               w1,w2,w3,jacmi)

!$acc update host(w1,w2,w3)
      ifielt = ifield
      ifield = 1
      call opcolv_acc (w1,w2,w3,bm1)
c!$acc update host(w1,w2,w3)
c      do i=1,lx1*ly1*lz1*nelv
c         write (6,*) 'dst1 ',w1(i)
c         write (6,*) 'dst1 ',w2(i)
c         write (6,*) 'dst1 ',w3(i)
c      enddo
      call dssum      (w1,nx1,ny1,nz1)
      call dssum      (w2,nx1,ny1,nz1)
      call dssum      (w3,nx1,ny1,nz1)
c!$acc update host(w1,w2,w3)
c      do i=1,lx1*ly1*lz1*nelv
c         write (6,*) 'dst2 ',w1(i)
c         write (6,*) 'dst2 ',w2(i)
c         write (6,*) 'dst2 ',w3(i)
c      enddo
c     stop
c!$acc update host(binvm1)
c      do i=1,lx1*ly1*lz1*nelv
c         write (6,*) 'binvm',binvm1(i,1,1,1)
c      enddo
c      stop
      call opcolv_acc (w1,w2,w3,binvm1)
!$acc update host(w1,w2,w3) !<- necessary
c      do i=1,lx1*ly1*lz1*nelv
c         write (6,*) 'opcolv',w1(i)
c         write (6,*) 'opcolv',w1(i)
c         write (6,*) 'opcolv',w1(i)
c      enddo
      ifield = ifielt

!$ACC END DATA

      return
      end

c-----------------------------------------------------------------------
      subroutine opadd2cm_acc (a1,a2,a3,b1,b2,b3,c)
      INCLUDE 'SIZE'
      REAL A1(lx1*ly1*lz1*lelt),
     $     A2(lx1*ly1*lz1*lelt),
     $     A3(lx1*ly1*lz1*lelt),
     $     B1(lx1*ly1*lz1*lelt),
     $     B2(lx1*ly1*lz1*lelt),
     $     B3(lx1*ly1*lz1*lelt),
     $     C
      NTOT1=NX1*NY1*NZ1*NELV
      if (ldim.eq.3) then
!$ACC PARALLEL LOOP PRESENT(a1,a2,a3,b1,b2,b3)
         do i=1,ntot1
            a1(i) = a1(i) + b1(i)*c
            a2(i) = a2(i) + b2(i)*c
            a3(i) = a3(i) + b3(i)*c
         enddo
!$ACC END PARALLEL
      else
!$ACC PARALLEL LOOP PRESENT(a1,a2,b1,b2)
         do i=1,ntot1
            a1(i) = a1(i) + b1(i)*c
            a2(i) = a2(i) + b2(i)*c
         enddo
!$ACC END PARALLEL
      endif
      return
      END
c-----------------------------------------------------------------------
      subroutine convop_acc(du,u)

      include 'SIZE'
      include 'TOTAL'
      include 'CTIMER'    ! Contains icalld

      real    du (lx1,ly1,lz1,lelt)
      real    u  (lx1,ly1,lz1,lelt)

      real cvx(lxd,lyd,lzd,nelv,3)

      real jggl(lxd,lx1),dggl(lxd,lx1)
      save jggl,dggl

      call chck('q11')

      if (icalld.eq.0) then
         tadvc=0.0
         call get_jggl(jggl,dggl) ! run on CPU and then acc copyin after
      call chck('q21')
      call chck('q31')
      endif

!$acc data copy(jggl,dggl)
      icalld=icalld+1
      nadvc=icalld
      etime1=dnekclock()

      call chck('q41')

      call convop_fst_3d_acc(du,u,vxd,vyd,vzd,bm1,jggl,dggl)

      call chck('q51')

      tadvc=tadvc+(dnekclock()-etime1)
!$acc end data

      return
      end
c-------------------------------------------------------------------
      subroutine convop_fst_3d_acc(du,u,c1,c2,c3,b,jj,dd)
c-------------------------------------------------------------------
      include 'SIZE'
      include 'DEALIAS'

c     apply convecting field c to scalar field u
c
c           T
c     du = J   ( C . grad Ju)


      real du(lx1,ly1,lz1,lelt),u(lx1,ly1,lz1,lelt),b(lx1,ly1,lz1,lelt)
      real c1(lxd*lyd*lzd,lelt),
     $     c2(lxd*lyd*lzd,lelt),
     $     c3(lxd*lyd*lzd,lelt)

      real jj(lxd,lx1)
      real dd(lxd,lx1)

      real wk1(lxd,ly1,lz1),wk2(lxd,ly1,lz1)
      real wk3(lxd,lyd,lz1),wk4(lxd,lyd,lz1)
      real wk5(lxd,lyd,lz1),wk6(lxd,lyd,lzd)

      integer e,p,q,r

      call chck('p11')

!$acc parallel loop present(du,u,c1,c2,c3,b,jj,dd) gang
!$acc&              private(wk1,wk2,wk3,wk4,wk5,wk6)

      do e=1,nelv

!$acc    loop collapse(2) vector
         do j=1,ly1*lz1
         do i=1,lxd
            w1 = 0.0
            w2 = 0.0
!$acc       loop seq
            do p=1,lx1
               w1 = w1 + dd(i,p) * u(p,j,1,e) ! Grad:  crs->fine
               w2 = w2 + jj(i,p) * u(p,j,1,e) ! Inter: crs->fine
            enddo
!$acc       end loop
            wk1(i,j,1)=w1
            wk2(i,j,1)=w2
         enddo
         enddo
!$acc    end loop

!$acc    loop collapse(3) vector
         do k=1,lz1
         do i=1,lxd
         do j=1,lyd
            w1 = 0.0
            w2 = 0.0
            w3 = 0.0
!$acc       loop seq
            do q=1,ly1
               w1 = w1 + jj(j,q) * wk1(i,q,k) ! JxD u
               w2 = w2 + dd(j,q) * wk2(i,q,k) ! DxJ u
               w3 = w3 + jj(j,q) * wk2(i,q,k) ! JxJ u
            enddo
!$acc       end loop
            wk3(i,j,k)=w1
            wk4(i,j,k)=w2
            wk5(i,j,k)=w3
         enddo
         enddo
         enddo
!$acc    end loop

!$acc    loop collapse(2) vector
         do i=1,lxd*lyd
         do k=1,lzd
            w1 = 0.0
            w2 = 0.0
            w3 = 0.0
!$acc       loop seq
            do r=1,lz1
               w1 = w1 + jj(k,r) * wk3(i,1,r) ! JxJxD u
               w2 = w2 + jj(k,r) * wk4(i,1,r) ! JxDxJ u
               w3 = w3 + dd(k,r) * wk5(i,1,r) ! DxJxJ u
            enddo
!$acc       end loop
            wk6(i,1,k)=w1*vxd(i,1,k,e)
     $                +w2*vyd(i,1,k,e)
     $                +w3*vzd(i,1,k,e) ! c1*ur+c2*us+c3*ut
         enddo
         enddo
!$acc    end loop

!! START COLLAPSING BACK with J'

!$acc    loop collapse(2) vector
         do i=1,lxd*lyd
         do k=1,lz1
            w1 = 0.0
!$acc       loop seq
            do r=1,lzd                        !  T
               w1 = w1 + jj(r,k) * wk6(i,1,r) ! J  x I x I w6
            enddo
!$acc       end loop
            wk5(i,1,k)=w1
         enddo
         enddo
!$acc    end loop

!$acc    loop collapse(3) vector
         do k=1,lz1
         do i=1,lxd
         do j=1,ly1
            w1 = 0.0
!$acc       loop seq
            do q=1,lyd
               w1 = w1 + jj(q,j) * wk5(i,q,k)
            enddo
!$acc       end loop
            wk2(i,j,k)=w1
         enddo
         enddo
         enddo
!$acc    end loop

!$acc    loop collapse(2) vector
         do j=1,ly1*lz1
         do i=1,lx1
            w1 = 0.0
!$acc       loop seq
            do p=1,lxd
               w1 = w1 + jj(p,i) * wk2(p,j,1)
            enddo
!$acc       end loop
c           du(i,j,1,e) = w1*b(i,j,1,e)
            du(i,j,1,e) = w1
         enddo
         enddo
!$acc    end loop
      enddo
!$acc end parallel loop

      return
      end

c-------------------------------------------------------------------
      subroutine plan4_acc_update_device
c-----------------------------------------------------------------------
      include 'SIZE'
      include 'TOTAL'    

      parameter (lg=lx1*ly1*lz1*lelt)
      parameter (maxcg=900)

      common /tdarray/ diagt(maxcg),upper(maxcg)

      common /scrcg/ d(lg), scalar(2)
      common /scrcg2/ r(lg), w(lg), p(lg), z(lg)

!$acc update device(vxlag,vylag,vzlag,tlag,vgradt1,vgradt2)
!$acc update device(abx1,aby1,abz1,abx2,aby2,abz2,vdiff_e)
!$acc update device(vtrans,vdiff,bfx,bfy,bfz,cflf,c_vx,fw)
!$acc update device(bmnv,bmass,bdivw,bx,by,bz,pm,bmx,bmy,bmz)
!$acc update device(vx,vy,vz,pr,t,vx_e,vy_e,vz_e)
!$acc update device(bbx1,bby1,bbz1,bbx2,bby2,bbz2,bxlag,bylag,bzlag)

!$acc update device(abx1,aby1,abz1,abx2,aby2,abz2)
!$acc update device(ab,bd)
!$acc update device(pr,pmlag,prlag,qtl,usrdiv)

!$acc update device(vxd,vyd,vzd)

!$acc update device(d,scalar,r,w,p,z)
!$acc update device(ibc_acc)
!$acc update device(c_vx)

      return
      end
c-----------------------------------------------------------------------
      subroutine hsmg_acc_update_device
      include 'SIZE'
      include 'HSMG'

!$acc update device(mg_nx)
!$acc update device(mg_ny,mg_nz)
!$acc update device(mg_nh,mg_nhz)
!$acc update device(mg_gsh_schwarz_handle)
!$acc update device(mg_gsh_handle)
!$acc update device(mg_rstr_wt_index)
!$acc update device(mg_mask_index)
!$acc update device(mg_solve_index)
!$acc update device(mg_fast_s_index)
!$acc update device(mg_fast_d_index)
!$acc update device(mg_schwarz_wt_index)
!$acc update device(mg_g_index)

!$acc update device(mg_jh)
!$acc update device(mg_jht)
!$acc update device(mg_jhfc )
!$acc update device(mg_jhfct)
!$acc update device(mg_ah)
!$acc update device(mg_bh)
!$acc update device(mg_dh)
!$acc update device(mg_dht)
!$acc update device(mg_zh)
!$acc update device(mg_rstr_wt)
!$acc update device(mg_mask)
!$acc update device(mg_fast_s)
!$acc update device(mg_fast_d)
!$acc update device(mg_schwarz_wt)
!$acc update device(mg_solve_e)
!$acc update device(mg_solve_r)
!$acc update device(mg_h1)
!$acc update device(mg_h2)
!$acc update device(mg_b)
!$acc update device(mg_g)
!$acc update device(mg_work)
!$acc update device(mg_work2)
!$acc update device(mg_worke)

!$acc update device(mg_imask)

!$acc update device(mg_h1_n)
!$acc update device(p_mg_h1)
!$acc update device(p_mg_b)
!$acc update device(p_mg_msk)

      return
      end
c-----------------------------------------------------------------------
      subroutine plan4_acc_update_host  
c-----------------------------------------------------------------------
      include 'SIZE'
      include 'TOTAL'    

      parameter (lg=lx1*ly1*lz1*lelt)
      parameter (maxcg=900)

      common /tdarray/ diagt(maxcg),upper(maxcg)

      common /scrcg/ d(lg), scalar(2)
      common /scrcg2/ r(lg), w(lg), p(lg), z(lg)

!$acc update host(vxlag,vylag,vzlag,tlag,vgradt1,vgradt2)
!$acc update host(abx1,aby1,abz1,abx2,aby2,abz2,vdiff_e)
!$acc update host(vtrans,vdiff,bfx,bfy,bfz,cflf,c_vx,fw)
!$acc update host(bmnv,bmass,bdivw,bx,by,bz,pm,bmx,bmy,bmz)
!$acc update host(vx,vy,vz,pr,t,vx_e,vy_e,vz_e)
!$acc update host(bbx1,bby1,bbz1,bbx2,bby2,bbz2,bxlag,bylag,bzlag)

!$acc update host(ab,bd)
!$acc update host(pr,pmlag,prlag,qtl,usrdiv)
!$acc update host(vxd,vyd,vzd)
!$acc update host(diagt,upper)
!$acc update host(d,scalar,r,w,p,z)
!$acc update host(ibc_acc)
!$acc update host(c_vx)

      return
      end
c-----------------------------------------------------------------------
c      subroutine convop_fst_3d_acc2(du,u,c1,c2,c3,b,jj,dd)
cc-------------------------------------------------------------------
c      include 'SIZE'
c
cc     apply convecting field c to scalar field u
cc
cc           T
cc     du = J   ( C . grad Ju)
c
c      real du(lx1,ly1,lz1,lelt),u(lx1,ly1,lz1,lelt),b(lx1,ly1,lz1,lelt)
c      real c1(lxd*lyd*lzd,lelt),
c     $     c2(lxd*lyd*lzd,lelt),
c     $     c3(lxd*lyd*lzd,lelt)
c
c      real jj(lxd,lx1)
c      real dd(lxd,lx1)
c
c      real wk1(lxd,ly1,lz1),wk2(lxd,ly1,lz1)
c      real wk3(lxd,lyd,lz1),wk4(lxd,lyd,lz1)
c      real wk5(lxd,lyd,lz1),wk6(lxd,lyd,lzd)
c
c      integer e,p,q,r
c
c      call chck('p11')
c
c!$acc data create(wk1,wk2,wk3,wk4,wk5,wk6)
c      do e=1,nelv
c!$acc    kernels
c         do j=1,ly1*lz1
c         do i=1,lxd
c            w1 = 0.0
c            w2 = 0.0
c!$acc       loop seq
c            do p=1,lx1
c               w1 = w1 + dd(i,p) * u(p,j,1,e) ! Grad:  crs->fine
c               w2 = w2 + jj(i,p) * u(p,j,1,e) ! Inter: crs->fine
c            enddo
c!$acc       end loop
c            wk1(i,j,1)=w1
c            wk2(i,j,1)=w2
c         enddo
c         enddo
c!$acc    end kernels
c
c!$acc    kernels
c         do k=1,lz1
c         do i=1,lxd
c         do j=1,lyd
c            w1 = 0.0
c            w2 = 0.0
c            w3 = 0.0
c!$acc       loop seq
c            do q=1,ly1
c               w1 = w1 + jj(j,q) * wk1(i,q,k) ! JxD u
c               w2 = w2 + dd(j,q) * wk2(i,q,k) ! DxJ u
c               w3 = w3 + jj(j,q) * wk2(i,q,k) ! JxJ u
c            enddo
c!$acc       end loop
c            wk3(i,j,k)=w1
c            wk4(i,j,k)=w2
c            wk5(i,j,k)=w3
c         enddo
c         enddo
c         enddo
c!$acc    end kernels
c
cc        l=0
c!$acc    kernels
c         do i=1,lxd*lyd
c         do k=1,lzd
c            w1 = 0.0
c            w2 = 0.0
c            w3 = 0.0
c!$acc       loop seq
c            do r=1,lz1
c               w1 = w1 + jj(k,r) * wk3(i,1,r) ! JxJxD u
c               w2 = w2 + jj(k,r) * wk4(i,1,r) ! JxDxJ u
c               w3 = w3 + dd(k,r) * wk5(i,1,r) ! DxJxJ u
cc              w3 = w3 + jj(k,r) * wk5(i,1,r) ! DxJxJ u
c            enddo
c!$acc       end loop
cc           l=l+1
cc           wk6(i,1,k)=w1*c(l,e,1)+w2*c(l,e,2)+w3*c(l,e,3) ! c1*ur+c2*us+c3*ut
c            wk6(i,1,k)=w1*vxd(i,1,k,e)
c     $                +w2*vyd(i,1,k,e)
c     $                +w3*vzd(i,1,k,e) ! c1*ur+c2*us+c3*ut
c         enddo
c         enddo
c!$acc    end kernels
c
c!$acc    update host(wk6)
c         do i=1,lx1*ly1*lz1
c            write (6,*) 'wk6=',wk6(i,1,1)
c         enddo
c         stop
c
c!! START COLLAPSING BACK with J'
c
c!$acc    loop collapse(2) vector
c         do i=1,lxd*lyd
c         do k=1,lz1
c            w1 = 0.0
c!$acc       loop seq
c            do r=1,lzd                        !  T
c               w1 = w1 + jj(r,k) * wk6(i,1,r) ! J  x I x I w6
c            enddo
c!$acc       end loop
c            wk5(i,1,k)=w1
c         enddo
c         enddo
c!$acc    end loop
c
c!$acc    loop collapse(3) vector
c         do k=1,lz1
c         do i=1,lxd
c         do j=1,ly1
c            w1 = 0.0
c!$acc       loop seq
c            do q=1,lyd
c               w1 = w1 + jj(q,j) * wk5(i,q,k)
c            enddo
c!$acc       end loop
c            wk2(i,j,k)=w1
c         enddo
c         enddo
c         enddo
c!$acc    end loop
c
c!$acc    loop collapse(2) vector
c         do j=1,ly1*lz1
c         do i=1,lx1
c            w1 = 0.0
c!$acc       loop seq
c            do p=1,lxd
c               w1 = w1 + jj(p,i) * wk2(p,j,1)
c            enddo
c!$acc       end loop
cc           du(i,j,1,e) = w1*b(i,j,1,e)
c            du(i,j,1,e) = w1
c         enddo
c         enddo
c!$acc    end loop
c      enddo
c!$acc end data
c
c      return
c      end
c-------------------------------------------------------------------
