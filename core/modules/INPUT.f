      module input_mod
      include 'SIZE'
c
c     Input parameters from preprocessors.
c
c     Note that in parallel implementations, we distinguish between
c     distributed data (LELT) and uniformly distributed data.
c
c     Input common block structure:
c
c     INPUT1:  REAL            INPUT5: REAL      with LELT entries
c     INPUT2:  INTEGER         INPUT6: INTEGER   with LELT entries
c     INPUT3:  LOGICAL         INPUT7: LOGICAL   with LELT entries
c     INPUT4:  CHARACTER       INPUT8: CHARACTER with LELT entries
c
      real param(200),rstim,vnekton
     $    ,cpfld(ldimt1,3)
     $    ,cpgrp(-5:10,ldimt1,3)
     $    ,qinteg(ldimt3,maxobj)
     $    ,uparam(20)
     $    ,atol(0:ldimt1)
     $    ,restol(0:ldimt1)
     $    ,fem_amg_param(15)
     $    ,crs_param(15)
     $    ,filterType
     $    ,connectivityTol

      common /input1/ param,rstim,vnekton,cpfld,cpgrp,qinteg,uparam,
     $                atol,restol,fem_amg_param,crs_param,
     $                filterType,connectivityTol

      integer matype(-5:10,ldimt1)
     $       ,nktonv,nhis,lochis(4,lhis+maxobj)
     $       ,ipscal,npscal,ipsco, ifldmhd
     $       ,nrst_rd,nrst_cm,irstv,irstt,irstim,nmember(maxobj),nobj
     $       ,ngeom,idpss(ldimt),fluid_partitioner,solid_partitioner
      common /input2/ matype,nktonv,nhis,lochis,ipscal,npscal,ipsco
     $               ,ifldmhd,nrst_rd,nrst_cm,irstv,irstt,irstim,nmember,nobj
     $               ,ngeom,idpss,fluid_partitioner,solid_partitioner

      logical         if3d,ifflow,ifheat,iftran,ifaxis,ifstrs,ifsplit
     $               ,ifmgrid
     $               ,ifadvc(ldimt1),ifdiff(ldimt1),ifdeal(ldimt1)
     $               ,iffilter(ldimt1),ifprojfld(0:ldimt1)
     $               ,iftmsh(0:ldimt1),ifdgfld(0:ldimt1),ifdg
     $               ,ifmvbd,ifchar,ifnonl(ldimt1)
     $               ,ifvarp(ldimt1),ifpsco(ldimt1),ifvps
     $               ,ifmodel,ifkeps,ifintq,ifcons
     $               ,ifxyo,ifpo,ifvo,ifto,iftgo,ifpso(ldimt1),iffmtin
     $               ,ifbo,ifanls,ifanl2,ifmhd,ifessr,ifpert,ifbase
     $               ,ifcvode,iflomach,ifexplvis,ifschclob,ifuservp
     $               ,ifcyclic,ifmoab,ifcoup, ifvcoup, ifusermv,ifreguo
     $               ,ifxyo_,ifaziv,ifneknek,ifneknekm,ifneknekc
     $               ,ifcvfld(ldimt1),ifdp0dt
     $               ,ifmpiio,ifcrrs,ifrich,ifvvisp
     $               ,ifbmap(ldimt1),ifjac0_abort

      common /input3/ if3d,ifflow,ifheat,iftran,ifaxis,ifstrs,ifsplit
     $               ,ifmgrid 
     $               ,ifadvc,ifdiff,ifdeal
     $               ,iffilter, ifprojfld
     $               ,iftmsh,ifdgfld,ifdg
     $               ,ifmvbd,ifchar,ifnonl
     $               ,ifvarp        ,ifpsco        ,ifvps
     $               ,ifmodel,ifkeps,ifintq,ifcons
     $               ,ifxyo,ifpo,ifvo,ifto,iftgo,ifpso        ,iffmtin
     $               ,ifbo,ifanls,ifanl2,ifmhd,ifessr,ifpert,ifbase
     $               ,ifcvode,iflomach,ifexplvis,ifschclob,ifuservp
     $               ,ifcyclic,ifmoab,ifcoup, ifvcoup, ifusermv,ifreguo
     $               ,ifxyo_,ifaziv,ifneknek,ifneknekm,ifneknekc
     $               ,ifcvfld,ifdp0dt
     $               ,ifmpiio,ifcrrs,ifrich,ifvvisp,ifbmap,ifjac0_abort

      logical         ifnav
      equivalence    (ifnav, ifadvc(1))

      character*1     hcode(11,lhis+maxobj)
      character*2     ocode(8)
      character*10    drivc(5)
      character*14    rstv,rstt
      character*40    textsw(100,2)
      character*132   initc(15)
      common /input4/ hcode,ocode,rstv,rstt,drivc,initc,textsw

      character*40    turbmod
      equivalence    (turbmod,textsw(1,1))

      character*132   reafle,fldfle,dmpfle,hisfle,schfle,orefle,nrefle
      common /cfiles/ reafle,fldfle,dmpfle,hisfle,schfle,orefle,nrefle

      character*132   session,path,re2fle,parfle,amgfile
      common /cfile2/ session,path,re2fle,parfle,amgfile

      integer cr_re2,fh_re2
      common /handles_re2/ cr_re2,fh_re2

      integer*8 re2off_b
      common /off_re2/ re2off_b
c
c proportional to LELT
c
      real, allocatable, target :: cb_input5(:)
      real, pointer :: xc(:,:),yc(:,:),zc(:,:)
     $               , bc(:,:,:,:)
     $               , curve(:,:,:)
     $               , cerror(:)

      integer, allocatable, target :: igroup(:),object(:,:,:)

      integer lbid
      parameter(lbid = 100)

      character*1, allocatable, target :: ccurve(:,:),cdof(:,:)
      character*3, allocatable, target :: cbc(:,:,:)
      character*3     cbc_bmap(lbid,ldimt1)
      integer         cbc_imap(lbid)
      integer         nbctype
      common /input8b/ cbc_bmap,cbc_imap,nbctype

      integer ieact(lelt),neact
      common /input9/ ieact,neact
c
c material set ids, BC set ids, materials (f=fluid, s=solid), bc types
c
      integer numsts
      parameter (numsts=50)
      
      integer numflu,numoth,numbcs 
     $       ,matindx(numsts),matids(numsts),imatie(lelt)
     $       ,ibcsts(numsts)
      common /inputmi/ numflu,numoth,numbcs,matindx,matids,imatie
     $                ,ibcsts
      
      integer bcf(numsts)
      common /inputmr/ bcf

      character*3 bctyps(numsts)
      common /inputmc/ bctyps

      integer out_mask(lelt)
      common /cbout_mask/ out_mask
c
c h refinement
c
      integer nhref, hrefcuts(lhref)
      common /hrefinei/ nhref, hrefcuts

      contains

      subroutine init
         implicit none
         integer ierr, ioff

         allocate(cb_input5(97*lelt + 30*(ldimt1+1)*lelt), stat=ierr)
         if (ierr.ne.0) call exitti('alloc cb_input5$',ierr)
         cb_input5 = 0

         ioff = 1
         xc(1:8,1:lelt) => cb_input5(ioff : ioff+8*lelt-1)
         ioff = ioff + 8*lelt
         yc(1:8,1:lelt) => cb_input5(ioff : ioff+8*lelt-1)
         ioff = ioff + 8*lelt
         zc(1:8,1:lelt) => cb_input5(ioff : ioff+8*lelt-1)
         ioff = ioff + 8*lelt
         bc(1:5,1:6,1:lelt,0:ldimt1) =>
     $      cb_input5(ioff : ioff+30*(ldimt1+1)*lelt-1)
         ioff = ioff + 30*(ldimt1+1)*lelt
         curve(1:6,1:12,1:lelt) => cb_input5(ioff : ioff+72*lelt-1)
         ioff = ioff + 72*lelt
         cerror(1:lelt) => cb_input5(ioff : ioff+lelt-1)

         allocate(igroup(lelt), stat=ierr)
         if (ierr.ne.0) call exitti('alloc igroup$',ierr)
         igroup = 0
         allocate(object(maxobj,maxmbr,2), stat=ierr)
         if (ierr.ne.0) call exitti('alloc object$',ierr)
         object = 0
         allocate(ccurve(12,lelt), stat=ierr)
         if (ierr.ne.0) call exitti('alloc ccurve$',ierr)
         ccurve = ''
         allocate(cdof(6,lelt), stat=ierr)
         if (ierr.ne.0) call exitti('alloc cdof$',ierr)
         cdof = ''
         allocate(cbc(6,lelt,0:ldimt1), stat=ierr)
         if (ierr.ne.0) call exitti('alloc cbc$',ierr)
         cbc = ''

      end subroutine init
      end module input_mod
