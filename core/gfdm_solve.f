c-----------------------------------------------------------------------
      subroutine gfdm_pres_solv(z,r,ug,wg)
c
      include 'SIZE'
      include 'INPUT'
      include 'PARALLEL'
      include 'TSTEP'
      include 'ZPER'
c
      real z (lx2*ly2*lz2*lelv),r (lx2*ly2*lz2*lelv)
      real ug(lx2*ly2*lz2*lelv),wg(lx2*ly2*lz2*lelv)
c
      mfld = 1
      if (ifmhd.and.ifield.gt.1) mfld = 2
c
      mp   = ngfdm_p(pst2lex(1))
      ms   = ngfdm_p(pst2lex(2))
      mt   = ngfdm_p(pst2lex(3))
c
      msfp = msp(pst2lex(1),mfld)
      msfs = msp(pst2lex(2),mfld)
      msft = msp(pst2lex(3),mfld)
c
      ntot  = nx2*ny2*nz2*nelv
c
      m     = ntot/mp
      nwave = mcex
      mpt   = mcex/(ms*mt)
c
      if (ifmhd) then
         dtbdi = 1   ! scale by dt*backwd-diff coefficient
      else
         dtbdi = bd(1)/dt   ! scale by dt*backwd-diff coefficient
      endif
c
c     write(6,*) 'msfp:',msfp,msfs,msft,mp,ms,mt,sp(msfp)
c     call outmat(sp(msfp),mp,mp,'S  x  ',istep)
c     call outmat(sp(msfs),ms,ms,'S  y  ',ifield)
c     call outmat(sp(msft),mt,mt,'S  z  ',mt)
c     call exitt
c
      if (if3d) then
         call copy    (ug,r,ntot)
c        write(6,*) nid,'SWAP_IP 1:',ntot
         call swap_ip (ug,tpn2,ntot)
         call mxm     (ug,m,sp (msfp),mp,wg,mp)
         call cexr    (ug,wg,m,mp,part_out,part_in,msg_id,wdsize,nid,np)
c        write(6,*) nid,'SWAP_IP 2:',ntot
         call swap_ip (ug,ind23,mcex)
         call mxm     (spt(msfs),ms,ug,ms,wg,mpt*mt)
         call mxm     (wg,ms*mpt,sp (msft),mt,ug,mt)
         call col2s2  (ug,wavep,dtbdi,mcex)
         call mxm     (ug,ms*mpt,spt(msft) ,mt,wg,mt)
         call mxm     (sp (msfs),ms,wg,ms,ug,mpt*mt)
         call swapt_ip(ug,ind23,mcex)
         call cextr   (wg,m,mp,ug,part_out,part_in,msg_id,wdsize,nid,np)
         call mxm     (wg,m,spt(msfp),mp,z,mp)
         call swapt_ip(z,tpn2,ntot)
      else
         call copy    (ug,r,ntot)
         call swap_ip (ug,tpn2,ntot)
         call mxm     (ug,m,sp (msfp),mp,wg,mp)
         call cexr    (ug,wg,m,mp,part_out,part_in,msg_id,wdsize,nid,np)
         call swap_ip (ug,ind23,mcex)
         call mxm     (spt(msfs),ms,ug,ms,wg,mpt)
         call col2s2  (wg,wavep,dtbdi,mcex)
         call mxm     (sp (msfs),ms,wg,ms,ug,mpt)
         call swapt_ip(ug,ind23,mcex)
         call cextr   (wg,m,mp,ug,part_out,part_in,msg_id,wdsize,nid,np)
         call mxm     (wg,m,spt(msfp),mp,z,mp)
         call swapt_ip(z,tpn2,ntot)
      endif
c
      return
      end
c-----------------------------------------------------------------------
