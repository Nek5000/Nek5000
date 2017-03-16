
      subroutine acc_copy_all_in()
      include 'SIZE'
      include 'HSMG'            ! Same array space as HSMG
      include 'GEOM'
      include 'INPUT'
      include 'MASS'
      include 'SOLN'
      include 'TSTEP'
      include 'CTIMER'
      include 'PARALLEL'


      common /scrhi/ h2inv (lx1,ly1,lz1,lelv)
      common /scrvh/ h1    (lx1,ly1,lz1,lelv),
     $     h2    (lx1,ly1,lz1,lelv)
      parameter (lt=lx1*ly1*lz1*lelt)
      common /scrmg/ e(2*lt),w(lt),r(lt)

      print*,'============in acc_copy_all_in============'
!$ACC ENTER DATA COPYIN(mg_imask)
!$ACC ENTER DATA CREATE(e,w,r)
!$ACC ENTER DATA COPYIN(mg_jht,mg_jh,mg_rstr_wt,mg_schwarz_wt)
!$ACC ENTER DATA COPYIN(mg_work,mg_mask,mg_fast_s,mg_fast_d)

      end
