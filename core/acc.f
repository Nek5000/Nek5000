c-----------------------------------------------------------------------
      subroutine plan4_acc_data_copyin()
c-----------------------------------------------------------------------
      include 'SIZE'
      include 'TOTAL'    

      common /scrns/ res1  (lx1,ly1,lz1,lelv)
     $ ,             res2  (lx1,ly1,lz1,lelv)
     $ ,             res3  (lx1,ly1,lz1,lelv)
     $ ,             dv1   (lx1,ly1,lz1,lelv)
     $ ,             dv2   (lx1,ly1,lz1,lelv)
     $ ,             dv3   (lx1,ly1,lz1,lelv)
     $ ,             respr (lx2,ly2,lz2,lelv)

      common /scrvh/ h1    (lx1,ly1,lz1,lelv)
     $ ,             h2    (lx1,ly1,lz1,lelv)

      REAL           dpr   (lx2,ly2,lz2,lelv)
      EQUIVALENCE   (dpr,dv1)


      integer icalld
      save    icalld
      data    icalld/0/
      if (icalld.eq.0) then ! ONE TIME ONLY

        icalld=1


ccc!$acc   enter data copyin (h1,h2,vtrans,respr,pr,vx,vy,vz) 
ccc!$acc   enter data copyin (res1,res2,res3,dv1,dv2,dv3)

!$acc   enter data copyin (v1mask,v2mask,v3mask,pmask,tmask,omask)
!$acc   enter data copyin (pmult,tmult,vmult)
!$acc   enter data copyin (dxm1,dxtm1,w3m1)
!$acc   enter data copyin (bm1,bm1lag,binvm1,bintm1)
!$acc   enter data copyin (jacm1,jacmi)
!$acc   enter data copyin (xm1,ym1,zm1)
!$acc   enter data copyin (unx,uny,unz,area)
!$acc   enter data copyin (rxm1,sxm1,txm1)
!$acc   enter data copyin (rym1,sym1,tym1)
!$acc   enter data copyin (rzm1,szm1,tzm1)
!$acc   enter data copyin (g1m1,g2m1,g3m1,g4m1,g5m1,g6m1)
!$acc   enter data copyin (cbc,bc)
!$acc   enter data copyin (abx1,aby1,abz1,abx2,aby2,abz2)
!$acc   enter data copyin (ab)
!$acc   enter data copyin (bfx,bfy,bfz)
!$acc   enter data copyin (vx,vy,vz,vx_e,vy_e,vz_e)
!$acc   enter data copyin (vtrans)
!$acc   enter data copyin (vxlag,vylag,vzlag)
!$acc   enter data copyin (bd)
!$acc   enter data copyin (pr,prlag,qtl,usrdiv)
!$acc   enter data copyin (vtrans,vdiff)
!$acc   enter data copyin (param,nelfld)

!$acc   enter data create (ibc_acc)

c!$acc   enter data copyin (tlag,vgradt1,vgradt2)
c!$acc   enter data copyin (vdiff_e)
c!$acc   enter data copyin (bq,t,usrdiv)


      endif

!$acc enter data copyin (c_vx)

      return
      end
c-----------------------------------------------------------------------
c     subroutine hmh_gmres_acc_data_copyin()
      subroutine hsmg_acc_data_copyin
      include 'SIZE'
      include 'HSMG'
      integer icalld
      save    icalld
      data    icalld/0/
      if (icalld.eq.0) then ! ONE TIME ONLY

!$ACC ENTER DATA COPYIN(mg_nx)
!$ACC ENTER DATA COPYIN(mg_ny,mg_nz)
!$ACC ENTER DATA COPYIN(mg_nh,mg_nhz)
!$ACC ENTER DATA COPYIN(mg_gsh_schwarz_handle)
!$ACC ENTER DATA COPYIN(mg_gsh_handle)
!$ACC ENTER DATA COPYIN(mg_rstr_wt_index)
!$ACC ENTER DATA COPYIN(mg_mask_index)
!$ACC ENTER DATA COPYIN(mg_solve_index)
!$ACC ENTER DATA COPYIN(mg_fast_s_index)
!$ACC ENTER DATA COPYIN(mg_fast_d_index)
!$ACC ENTER DATA COPYIN(mg_schwarz_wt_index)
!$ACC ENTER DATA COPYIN(mg_g_index)

!$ACC ENTER DATA COPYIN(mg_jh)
!$ACC ENTER DATA COPYIN(mg_jht)
!$ACC ENTER DATA COPYIN(mg_jhfc )
!$ACC ENTER DATA COPYIN(mg_jhfct)
!$ACC ENTER DATA COPYIN(mg_ah)
!$ACC ENTER DATA COPYIN(mg_bh)
!$ACC ENTER DATA COPYIN(mg_dh)
!$ACC ENTER DATA COPYIN(mg_dht)
!$ACC ENTER DATA COPYIN(mg_zh)
!$ACC ENTER DATA COPYIN(mg_rstr_wt)
!$ACC ENTER DATA COPYIN(mg_mask)
!$ACC ENTER DATA COPYIN(mg_fast_s)
!$ACC ENTER DATA COPYIN(mg_fast_d)
!$ACC ENTER DATA COPYIN(mg_schwarz_wt)
!$ACC ENTER DATA COPYIN(mg_solve_e)
!$ACC ENTER DATA COPYIN(mg_solve_r)
!$ACC ENTER DATA COPYIN(mg_h1)
!$ACC ENTER DATA COPYIN(mg_h2)
!$ACC ENTER DATA COPYIN(mg_b)
!$ACC ENTER DATA COPYIN(mg_g)
!$ACC ENTER DATA COPYIN(mg_work)
!$ACC ENTER DATA COPYIN(mg_work2)
!$ACC ENTER DATA COPYIN(mg_worke)

!$ACC ENTER DATA COPYIN(mg_imask)

!$ACC ENTER DATA COPYIN(mg_h1_n)
!$ACC ENTER DATA COPYIN(p_mg_h1)
!$ACC ENTER DATA COPYIN(p_mg_b)
!$ACC ENTER DATA COPYIN(p_mg_msk)
      endif

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

!$ACC ENTER DATA COPYIN(work,work2)
!$ACC ENTER DATA COPYIN(mg_mask,mg_imask,pmask)
!$ACC ENTER DATA COPYIN(mg_jht,mg_jh,mg_rstr_wt,mg_schwarz_wt)
!$ACC ENTER DATA COPYIN(mg_work,mg_fast_s,mg_fast_d)
!$ACC ENTER DATA COPYIN(h_gmres,w_gmres,v_gmres,z_gmres)
!$ACC ENTER DATA COPYIN(c_gmres,s_gmres,x_gmres,gamma_gmres)
!$ACC ENTER DATA COPYIN(r_gmres)
!$ACC ENTER DATA COPYIN(ml_gmres,mu_gmres)

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

!$ACC DATA PRESENT (d (nx1,nx1))
!$ACC&     PRESENT (u1(nx1,ny1,nz1,nelt))
!$ACC&     PRESENT (u2(nx1,ny1,nz1,nelt),u3(nx1,ny1,nz1,nelt))
!$ACC$     PRESENT (v1(nx1,ny1,nz1,nelt),v2(nx1,ny1,nz1,nelt))
!$ACC&     PRESENT (v3(nx1,ny1,nz1,nelt))
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
                     tmpu1 = tmpu1 + d(i,l)*u1(l,j,k,e)
                     tmpu2 = tmpu2 + d(j,l)*u2(i,l,k,e)
                     tmpu3 = tmpu3 + d(k,l)*u3(i,j,l,e)
                  enddo
                  v1(i,j,k,e) = tmpu1
                  v2(i,j,k,e) = tmpu2
                  v3(i,j,k,e) = tmpu3
               enddo
            enddo
         enddo
      enddo
!$ACC END PARALLEL LOOP
!$ACC END DATA
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
     $ ,             ta4   (lx1*ly1*lz1*lelv*ldimt1)

      n = lx1*ly1*lz1*nelv
      intype = -1

      scale = -1./3.
      if (ifstrs) scale =  2./3.

!$acc data present(resv1,resv2,resv3,h1,h2) create(ta1,ta2,ta3,ta4)
!$ac& present(vdiff,qtl,pr,vtrans,vx,vy,vz,bfx,bfy,bfz,bd)
      dtbd = bd(1)/dt

!$acc parallel loop 
      do i=1,n
         ta4(i)=vdiff (i,1,1,1,1)*qtl(i,1,1,1)+scale*pr(i,1,1,1)
         h1 (i)=vdiff (i,1,1,1,1)
         h2 (i)=vtrans(i,1,1,1,1)*dtbd
      enddo
!$acc end parallel loop 

      call wgradm1_acc  (ta1,ta2,ta3,ta4,nelv)
      call ophx_acc     (resv1,resv2,resv3,vx,vy,vz,h1,h2)

!$acc parallel loop
      do i=1,n
         resv1(i)=bfx(i,1,1,1)-resv1(i)-ta1(i)
         resv2(i)=bfy(i,1,1,1)-resv2(i)-ta2(i)
         resv3(i)=bfz(i,1,1,1)-resv3(i)-ta3(i)
      enddo
!$acc end parallel loop 

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
!$acc& present(dxm1,w3m1,ux,uy,uz,u) create(ur,us,ut)

!$acc parallel loop gang
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
      subroutine ophx_acc (out1,out2,out3,inp1,inp2,inp3,h1,h2)

c     OUT = (H1*A+H2*B) * INP  

      include 'SIZE'
      include 'INPUT'
      include 'SOLN'
      real out1(1),out2(1),out3(1),inp1(1),inp2(1),inp3(1),h1(1),h2(1)

      imesh = 1

c     IF (IFSTRS) THEN
c        MATMOD = 0
c        CALL AXHMSF (OUT1,OUT2,OUT3,INP1,INP2,INP3,H1,H2,MATMOD)
c     ELSE

c     Later - we can make this fast

      call axhelm_acc(out1,inp1,h1,h2,imesh,1)
      call axhelm_acc(out2,inp2,h1,h2,imesh,2)
      call axhelm_acc(out3,inp3,h1,h2,imesh,3)

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

      call makeuf_acc    ! paul and som
      call advab_acc
      call makeabf_acc
      call makebdf_acc

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
      real a1(1),a2(1),a3(1),c(1)

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

      real ta1 (lx1*ly1*lz1*lelv)
     $   , ta2 (lx1*ly1*lz1*lelv)
     $   , ta3 (lx1*ly1*lz1*lelv)


!$acc data create(ta1,ta2,ta3) present(vx,vy,vz,bfx,bfy,bfz,bm1)

      n = lx1*ly1*lz1*nelv

      call convop_acc  (ta1,vx)
      call convop_acc  (ta2,vy)
      call convop_acc  (ta3,vz)

!$acc parallel loop
      do i=1,n
         bfx(i,1,1,1)=bfx(i,1,1,1)-ta1(i)*bm1(i,1,1,1)
         bfy(i,1,1,1)=bfy(i,1,1,1)-ta2(i)*bm1(i,1,1,1)
         bfz(i,1,1,1)=bfz(i,1,1,1)-ta3(i)*bm1(i,1,1,1)
      enddo
!$acc end loop
!$acc end data 

      return
      end
c-----------------------------------------------------------------------
      subroutine convop_acc(du,u)

      include 'SIZE'
      include 'TOTAL'
      include 'CTIMER'    ! Contains icalld

      real    du (lx1,ly1,lz1,lelt)
      real    u  (lx1,ly1,lz1,lelt)

      real jggl(lxd,lx1),dggl(lxd,lx1)
      save jggl,dggl


      if (icalld.eq.0) then
         tadvc=0.0
         call get_jggl(jggl,dggl) ! run on CPU and then acc copyin after
!$acc    enter data copyin(jggl,dggl)
      endif

      icalld=icalld+1
      nadvc=icalld
      etime1=dnekclock()

      call convop_fst_3d_acc(du,u,c_vx,bm1,jggl,dggl)

      tadvc=tadvc+(dnekclock()-etime1)

      return
      END
c-------------------------------------------------------------------
      subroutine convop_fst_3d_acc(du,u,c,b,jj,dd)
c-------------------------------------------------------------------
      include 'SIZE'

c     apply convecting field c to scalar field u
c
c           T
c     du = J   ( C . grad Ju)


      real du(lx1,ly1,lz1,lelt)  , u(lx1,ly1,lz1,lelt)
      real  c(lxd*lyd*lzd,lelv,3), b(lx1,ly1,lz1,lelt)

      real jj(lxd,lx1)
      real dd(lxd,lx1)

      real wk1(lxd,ly1,lz1),wk2(lxd,ly1,lz1)
      real wk3(lxd,lyd,lz1),wk4(lxd,lyd,lz1)
      real wk5(lxd,lyd,lz1),wk6(lxd,lyd,lzd)

      integer e,p,q,r

!$acc parallel loop present(du,u,c,b,jj,dd) gang
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

         l=0
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
            l=l+1
            wk6(i,1,k)=w1*c(l,e,1)+w2*c(l,e,2)+w3*c(l,e,3) ! c1*ur+c2*us+c3*ut
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
            du(i,j,1,e) = w1*b(i,j,1,e)
         enddo
         enddo
!$acc    end loop

      enddo
!$acc end parallel loop

      return
      end

c-------------------------------------------------------------------
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

!$acc data create (ta1,ta2,ta3,ab)
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
      call add2s1 (bfz,ta3,ab0,ntot1)


cc    if(.not.iflomach)
      call col2_acc (bfx,vtrans,ntot1)
      call col2_acc (bfy,vtrans,ntot1)
      call col2_acc (bfz,vtrans,ntot1)

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
!$acc parallel loop
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
!$acc parallel loop
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
!$acc parallel loop
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
!$acc parallel
            do e = 1, nelv
            do k = 1, lz1
            do j = 1, ly1
            do i = 1, lx1
               ta1(i,j,k,e) = vxlag(i,j,k,e,ilag-1)  * 
     $                        bm1(i,j,k,ilag-1+e) *
     $                        bd(ilag+1)
               ta2(i,j,k,e) = vylag(i,j,k,e,ilag-1)  * 
     $                        bm1(i,j,k,ilag-1+e) *
     $                        bd(ilag+1)
               ta3(i,j,k,e) = vzlag(i,j,k,e,ilag-1)  * 
     $                        bm1(i,j,k,ilag-1+e) *
     $                        bd(ilag+1)
            enddo
            enddo
            enddo
            enddo
!$acc end parallel
         endif
c        INLINED:
c        call opadd2_acc(TB1,TB2,TB3,TA1,TA2,TA3)
c         =
c          CALL ADD2(TB1,TA1,NTOT1)
c          CALL ADD2(TB2,TA2,NTOT1)
c          IF(NDIM.EQ.3)CALL ADD2(TB3,TA3,NTOT1)
!$acc parallel
         do e = 1, nelv
         do k = 1, lz1
         do j = 1, ly1
         do i = 1, lx1
            tb1(i,j,k,e) = tb1(i,j,k,e) + ta1(i,j,k,e)
            tb2(i,j,k,e) = tb2(i,j,k,e) + ta2(i,j,k,e)
            tb3(i,j,k,e) = tb3(i,j,k,e) + ta1(i,j,k,e)
         enddo
         enddo
         enddo
         enddo
!$acc end parallel
      enddo

c     INLINED:
c     call opadd2col_acc(BFX,BFY,BFZ,TB1,TB2,TB3,h2)
!$acc parallel
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

!$acc    data present(ab,vx,vy,vz,vxlag,vylag,vzlag,vx_e,vy_e,vz_e)
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

      write(6,*) 'checker: ',s3

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
c      implicit none

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


!$ACC DATA PRESENT(u1r,u1s,u1t,u2r,u2s,u2t,u3r,u3s,u3t,w1,w2,w3)
!$ACC& PRESENT(jacmi,rxm1,sxm1,txm1,rym1,sym1,tym1,rzm1,szm1,tzm1)

!$ACC PARALLEL LOOP COLLAPSE(4) GANG WORKER VECTOR
         do e = 1,nelt
            do k = 1,lz1
               do j = 1,ly1
                  do i = 1,lx1
           
                   l = i + (j-1)*lx1 + (k-1)*lx1*ly1 
                   
                     w1(i,j,k,e)= (u3r(i,j,k,e)*rym1(i,j,k,e)
     $                    + u3s(i,j,k,e)*sym1(i,j,k,e)
     $                    + u3t(i,j,k,e)*tym1(i,j,k,e)
     $                    - u2r(i,j,k,e)*rzm1(i,j,k,e)
     $                    - u2s(i,j,k,e)*szm1(i,j,k,e)
     $                    - u2t(i,j,k,e)*tzm1(i,j,k,e))*jacmi(l,e)

                     w2(i,j,k,e)= (u1r(i,j,k,e)*rzm1(i,j,k,e)
     $                    + u1s(i,j,k,e)*szm1(i,j,k,e)
     $                    + u1t(i,j,k,e)*tzm1(i,j,k,e)
     $                    - u3r(i,j,k,e)*rxm1(i,j,k,e)
     $                    - u3s(i,j,k,e)*sxm1(i,j,k,e)
     $                    - u3t(i,j,k,e)*txm1(i,j,k,e))*jacmi(l,e)

                     w3(i,j,k,e)= (u2r(i,j,k,e)*rxm1(i,j,k,e)
     $                    + u2s(i,j,k,e)*sxm1(i,j,k,e)
     $                    + u2t(i,j,k,e)*txm1(i,j,k,e)
     $                    - u1r(i,j,k,e)*rym1(i,j,k,e)
     $                    - u1s(i,j,k,e)*sym1(i,j,k,e)
     $                    - u1t(i,j,k,e)*tym1(i,j,k,e))*jacmi(l,e)
                  enddo
               enddo
            enddo
         enddo
!$ACC END PARALLEL LOOP

!$ACC END DATA

c         write(6,*) iter, ' iter'


      flop_a = flop_a + (18.*nxyz*lx1)*nelt



      return
      end
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

!$ACC DATA PRESENT(u1r,u1s,u1t,u2r,u2s,u2t,u3r,u3s,u3t)
!$ACC&    PRESENT(u1,u2,u3,d)

!$ACC PARALLEL LOOP COLLAPSE(4) GANG WORKER VECTOR
!$ACC&     private(tmpr1,tmpr2,tmpr3,tmps1,tmps2,tmps3,
!$ACC&             tmpt1,tmpt2,tmpt3)
!dir$ NOBLOCKING
      do e = 1,nelt
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

!$ACC END DATA

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

!$ACC END DATA




ccc      if (ifavg.and..not.ifcyclic) then
ccc
ccc         ifielt = ifield
ccc         ifield = 1
ccc
ccc
ccc         call opcolv(w1,w2,w3,bm1)
ccc         call opdssum(w1,w2,w3)
ccc         call opcolv(w1,w2,w3,binvm1)
ccc
ccc!$ACC UPDATE DEVICE(w1,w2,w3)
ccc
ccc
ccc         ifield = ifielt
ccc      endif
ccc

      return
      end

c-----------------------------------------------------------------------
      subroutine crespsp_acc (respr)

C     Compute startresidual/right-hand-side in the pressure

      include 'SIZE'
      include 'TOTAL'

      real           respr (lx1,ly1,lz1,lelv)
c
      real           ta1   (lx1,ly1,lz1,lelv)
     $ ,             ta2   (lx1,ly1,lz1,lelv)
     $ ,             ta3   (lx1,ly1,lz1,lelv)
     $ ,             wa1   (lx1*ly1*lz1*lelv)
     $ ,             wa2   (lx1*ly1*lz1*lelv)
     $ ,             wa3   (lx1*ly1*lz1*lelv)
     $ ,             wa4   (lx1*ly1*lz1*lelv*ldimt1)
      real           w1    (lx1,ly1,lz1)
     $ ,             w2    (lx1,ly1,lz1)
     $ ,             w3    (lx1,ly1,lz1)

      character cb*3,c1*1
      integer e,f

!$acc routine(facind) seq
      
      nxyz1  = nx1*ny1*nz1
      ntot1  = nxyz1*nelv
      ntot2  = LX1*LY1*LZ1*LELV
      nfaces = 2*ndim

!$ACC ENTER DATA CREATE (ta1,ta2,ta3,wa1,wa2,wa3,wa4,w1,w2,w3)
c     -mu*curl(curl(v))
      call op_curl_acc (ta1,ta2,ta3,vx_e,vy_e,vz_e)

      call op_curl_acc (wa1,wa2,wa3,ta1,ta2,ta3)

      call rzero_acc(ta1,ntot1)
      call rzero_acc(ta2,ntot1)
      call rzero_acc(ta3,ntot1)
      ! No support for lowmach
      !call wgradm1_acc (ta1,ta2,ta3,qtl,nelt)
      call opcolv_acc  (wa1,wa2,wa3,bm1)

      scale = -4./3. 
      call opadd2cm_acc(wa1,wa2,wa3,ta1,ta2,ta3,scale)

c compute stress tensor for ifstrs formulation - variable viscosity Pn-Pn
      if (ifstrs) 
     $   call exitti('ifstrs not yet support on gpu$',nelv)

!$acc update device (vdiff,vtrans)
      call invcol3_acc  (wa4,vdiff,vtrans,ldimt1*ntot1)
        
      call opcolv_acc   (wa1,wa2,wa3,wa4)

c     add old pressure term because we solve for delta p 

      call invers2_acc (ta1,vtrans,ntot1)
      call rzero_acc   (ta2,ntot1)

!$ACC UPDATE HOST(pr)
      call bcdirsc (pr)
!$ACC UPDATE DEVICE(pr)

      call axhelm_acc  (respr,pr,ta1,ta2,1,1)
      call chsign_acc  (respr,ntot1)

c     add explicit (NONLINEAR) terms 
      n = nx1*ny1*nz1*nelv

!$acc parallel loop present(ta1,ta2,ta3,wa1,wa2,wa3)
      do i=1,n
         ta1(i,1,1,1) = bfx(i,1,1,1)/vtrans(i,1,1,1,1)-wa1(i)
         ta2(i,1,1,1) = bfy(i,1,1,1)/vtrans(i,1,1,1,1)-wa2(i)
         ta3(i,1,1,1) = bfz(i,1,1,1)/vtrans(i,1,1,1,1)-wa3(i)
      enddo
!$acc end parallel

      call dssum (ta1,nx1,ny1,nz1)
      call dssum (ta2,nx1,ny1,nz1)
      call dssum (ta3,nx1,ny1,nz1)

!$acc parallel loop present(ta1,ta2,ta3,binvm1)
      do i=1,n
         ta1(i,1,1,1) = ta1(i,1,1,1)*binvm1(i,1,1,1)
         ta2(i,1,1,1) = ta2(i,1,1,1)*binvm1(i,1,1,1)
         ta3(i,1,1,1) = ta3(i,1,1,1)*binvm1(i,1,1,1)
      enddo
!$acc end parallel

      dtbd = bd(1)/dt
c     call admcol3(respr,qtl,bm1,dtbd,ntot1)

!$acc update host(respr,ta1,ta2,ta3)
      call crespsp_face_update(respr,ta1,ta2,ta3,dtbd,w1,w2,w3)
!$acc update device(respr,ta1,ta2,ta3)           


!$acc parallel loop gang
!$acc&private(w1,w2,w3)
!$acc&present(ta1,ta2,ta3,w3m1,rxm1,rym1,rzm1
!$acc&         ,sxm1,sym1,szm1,txm1,tym1,tzm1
!$acc&         ,dxtm1,bm1,qtl,respr)
      do e=1,nelv
!$acc  loop vector
       do i=1,lx1*ly1*lz1
         w1(i,1,1) = (rxm1(i,1,1,e)*ta1(i,1,1,e) ! Jacobian
     $               +rym1(i,1,1,e)*ta2(i,1,1,e) ! included
     $               +rzm1(i,1,1,e)*ta3(i,1,1,e))*w3m1(i,1,1)
         w2(i,1,1) = (sxm1(i,1,1,e)*ta1(i,1,1,e)
     $               +sym1(i,1,1,e)*ta2(i,1,1,e)
     $               +szm1(i,1,1,e)*ta3(i,1,1,e))*w3m1(i,1,1)
         w3(i,1,1) = (txm1(i,1,1,e)*ta1(i,1,1,e)
     $               +tym1(i,1,1,e)*ta2(i,1,1,e)
     $               +tzm1(i,1,1,e)*ta3(i,1,1,e))*w3m1(i,1,1)
       enddo

!$acc  loop vector collapse(3)
       do k=1,nz1
       do j=1,ny1
       do i=1,nx1
          t1 = 0.0
!$acc     loop seq
          do l=1,nx1
             t1 = t1 + dxtm1(i,l)*w1(l,j,k) ! D^T
     $               + dxtm1(j,l)*w2(i,l,k)
     $               + dxtm1(k,l)*w3(i,j,l)
          enddo
          respr(i,j,k,e) = respr(i,j,k,e) + t1
     $                   + dtbd*bm1(i,j,k,e)*qtl(i,j,k,e)
       enddo
       enddo
       enddo

      enddo
!$acc end parallel


C     Orthogonalize to (1,1,...,1)T for all-Dirichlet case
      call ortho_acc (respr)

!$acc exit data

      return
      end
c----------------------------------------------------------------------
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
      subroutine chk2(s3,a)
      include 'SIZE'
      character*3 s3
      parameter (lt=lx1*ly1*lz1*lelt)
      real a(lt)

      n=nx1*ny1*nz1*nelt

      amx=glamax_acc(a,n)
      ams=glasum_acc(a,n)

      write(6,*) s3,   ' dev  max ',amx
      write(6,*) '   ','      sum ',ams

      return
      end
c-----------------------------------------------------------------------
      subroutine chk3(s3,a)
      include 'SIZE'
      character*3 s3

      parameter (lt=lx1*ly1*lz1*lelt)
      real a(lt)

      n=nx1*ny1*nz1*nelt

      amx=glamax(a,n)
      ams=glasum(a,n)

      write(6,*) s3,   ' host max ',amx
      write(6,*) '   ','      sum ',ams

      return
      end
c-----------------------------------------------------------------------
      subroutine chk2n(s4,a,n)
      include 'SIZE'
      character*4 s4

      real a(n)

      amx=glamax_acc(a,n)
      ams=glasum_acc(a,n)

      write(6,*) s4,     ' dev  max ',amx
      write(6,*) '    ', '      sum ',ams

      return
      end
c-----------------------------------------------------------------------
      subroutine chk3n(s4,a,n)
      include 'SIZE'
      character*4 s4

      real a(n)

      amx=glamax(a,n)
      ams=glasum(a,n)

      write(6,*) s4,     ' host max ',amx
      write(6,*) '    ', '      sum ',ams

      return
      end
c-----------------------------------------------------------------------
      subroutine chk4(s4,a,n1)
      include 'SIZE'

      real a(n1,n1,n1,nelv)
      integer e
      character*4 s4

      do e = 1,nelv
         max = glamax_acc(a(1,1,1,e),n1**3)
         if (e.eq.1) then
            write (6,*) s4,    ' e=',e,' max=',max
         else
            write (6,*) '    ',' e=',e,' max=',max
         endif
      enddo

      do e = 1,nelv
         sum = glasum_acc(a(1,1,1,e),n1**3)
         if (e.eq.1) then
            write (6,*) s4,    ' e=',e,' sum=',sum
         else
            write (6,*) '    ',' e=',e,' sum=',sum
         endif
      enddo

      return
      end
c-----------------------------------------------------------------------
