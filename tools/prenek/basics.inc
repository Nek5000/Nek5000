C     Utilities to be included in all subroutines in Preprocessor Program
C     BASICS.INC

      parameter (nelm=0 500 000,nxm=7,nym=7,nzm=7,maxfld=5)
      parameter (lx1=nxm,ly1=nym,lz1=nzm)

      PARAMETER ( MAXCRV=4,MAXIBC=4,MAXHIS=2,MPSCAL=4)
      PARAMETER ( MSOBJS=5, MFACES=200)

      LOGICAL     IFHEAT,IFFLOW,IFTRAN,IFNAV,IFAXIS,IFCEIL,IFERROR
     $,IFSTRS,IFADVC(11),IFTMSH(0:11),IFPLOT(NELM),IFDRAX,IFPSCO(9)
     $,IFSPLIT,IFMGRID,IFMODEL,IFKEPS,IF3BUT
     $,IFMVBD,IFCHAR,IFNEEDC,IFUNIV,IFEXPT,ifdssum,ifavgupt
     $,ifmid,ifcstd,ifmerge,ifconj_merge

      COMMON /LOGI/IFHEAT,IFFLOW,IFTRAN,IFAXIS,IFCEIL,IFERROR
     $,IFSTRS,IFADVC,IFTMSH,IFPLOT,IFDRAX,IFPSCO
     $,IFSPLIT,IFMGRID,IFMODEL,IFKEPS,IF3BUT
     $,IFMVBD,IFCHAR,IFNEEDC,IFUNIV,IFEXPT,ifdssum,ifavgupt
     $,ifmid,ifcstd,ifmerge,ifconj_merge

      COMMON /OBJECT/ NVOBJS,NSOBJS,NEOBJS,NPOBJS,NFACE(MSOBJS),
     $ILSURF(2,MSOBJS,MFACES)

      COMMON /OBJECC/ SOBJ(MSOBJS)
      CHARACTER SOBJ*20

      LOGICAL      IFXYO,IFXYH,IFPO,IFVO,IFTO,IFTGO,IFMOVB,IFPLAN
     $            ,IFPROJ,IFDMSH,IFZOOM,IFGNGM,IFCFBD,IFCLRB,IFMLTI
     $            ,IFFMTIN,ifptplt,ifauto,if_unstructure
      COMMON /LOGO/IFXYO,IFXYH,IFPO,IFVO,IFTO,IFTGO,IFMOVB,IFPLAN
     $            ,IFPROJ,IFDMSH,IFZOOM,IFGNGM,IFCFBD,IFCLRB,IFMLTI
     $            ,IFFMTIN,ifptplt,ifauto,if_unstructure

      common /nek2r/ curve(6,12,nelm),bc(5,6,nelm,0:maxfld)
     $         ,pi,theta,phi,vhobs(3),xhobs(3),yhobs(3)
     $         ,zpts(nxm),wght(nxm),dfdr(nxm,nxm),dfdrt(nxm,nxm)
     $         ,dgdr(nxm,nxm),dgdrt(nxm,nxm)
     $         ,points(3,3),pcond(11),prhocp(11)
     $         ,dmptim(1000),cprop(-5:10,10,3),cerror(6,nelm)
     $         ,veltyp

      common /nek2i/ nptpro,ncsegs,nrake
     $         ,ncurve,iellev,isrt(nelm)
     $         ,numapt(nelm),iletap,maxlet,ilevel,nsides,ncorns
     $         ,nedges,igrp
     $         ,lochis(4,maxhis)
     $         ,npscal,ncond
     $         ,ipsco,ipsch,idstep(1000)
     $         ,igroup(nelm),isobj
     $         ,ibc(6,nelm,0:maxfld) ! For periodic with nel>1e6
     $         ,matype(-5:10,0:maxfld)

      COMMON /CHAAR/LINE,SESION,ITEM,CHOICE,OLDVER,CMD(500),RSTV,RSTT
     $,X13,CPARAM(500),HCODE(11,MAXHIS),OCODE(20),CBC1,CBC3,DARROW
     $,PLANE,PLANE2,VPROP(-5:10,10,3),TEXTSW(100,2),S
     $,QUANTY,COMPON,SCALPT,VECTPT,XYATTR,SUATTR,PLFORM
     $,INBC(MAXIBC),CCURVE(12,NELM),CBC(6,NELM,0:MAXFLD),LETAPT(NELM)
     $,REQD,DRIVC(20),INITC(15),INITP(15),PSNAME(9),DERIV,DATE
     $,PLTMOD
     $,session_name
      character*14 session_name

      PARAMETER (LHI=50)
      COMMON /IBASIC/ NLINF,NLINP,NLINR,NHI,IEH(LHI)
      COMMON /BASICSI/ IN,NEL,NELF,NX,NY,NZ,NFLDS,NDIM,NELV,NELT
     $,NODES(MAXFLD),NODESB(MAXFLD),ICHOIC,NCHOIC,ICRV(NELM),nlevel
     $,NKTONV,LOCLIN,LPSDMP,NPSDMP

       COMMON /BASICSR/DXMIN,DYMIN
     $,CWRITE,YBS(100),YTS(100),XLMEN,XRMEN
     $,HEIGHT(nelm),VNEKTON,VNEKOLD

      COMMON /MATS/XPTS(NXM,NYM,NELM),YPTS(NXM,NYM,NELM)
     $,ZPTN(NXM,NYM,NELM),XXPTS(NXM),YYPTS(NYM)
     $,x(8,nelm),y(8,nelm),z(8,nelm)
     $,XCEN(NELM),YCEN(NELM),ZCEN(NELM),RCEN(NELM)

       COMMON /MATSR/ZDEPTH(NELM,6)
     $,SIDES(NELM,6,4),MASKEL(NELM,MAXFLD)
     $,zgml(nxm)

       COMMON /MATSI/adje(6,NELM),adjf(6,NELM),IZDPTH(NELM*6,2)
       integer adje,adjf

      common /crparam/ param(500),sarrow,hfish
      common /ciparam/ nparam,ncmd,nhispt
     $,irstt,irstv,irstim,isegdum
      common/crparam/rstim

C     MASK = 1 FOR NORMAL ELEMENT, = 0 FOR NO ELEMENT
C     IPER ARGUMENTS ARE IEL, ISIDE, III ,IFIELD. III=1 FOR IELP; =2 ISIDEP

      CHARACTER LINE*80,LINES(80),SESION*14,OCODE*2,CCURVE,REQD(100)
     $,CPARAM*40,CMD*10,OLDVER*4,RSTV*10,RSTT*10,X13*80,INBC*60
     $,PLANE*20,PLANE2*20,DRIVC*80,VCOND*80,INITC*80,INITP*80,VPROP*80
     $,QUANTY*20,COMPON*1,SCALPT*20,VECTPT*20,XYATTR*20,SUATTR*20
     $,PLFORM*20,CBC1,DARROW*20,PSNAME*5,DERIV*5,TEXTSW*40,TURBMOD*40
      CHARACTER CBC*3,CBC3*3,CRVTYP,ANS,HCODE,LETAPT,S*132,DATE*28
      CHARACTER*10 button,ARG(5),PLTMOD
      CHARACTER*26 ITEM(40),CHOICE
      EQUIVALENCE(LINE   ,LINES(1))
      EQUIVALENCE(LINE   ,ANS)
      EQUIVALENCE(EQTYPE,PARAM(16))
      EQUIVALENCE(AXIS  ,PARAM(17))
      EQUIVALENCE(GRID  ,PARAM(18))
      EQUIVALENCE(IFNAV ,IFADVC(1))
      EQUIVALENCE(TURBMOD,TEXTSW(1,1))
132   FORMAT(' ',2G14.6,A32)
C
C     Stuff that used to be in basicsp.inc
      COMMON /SPOST/ TDUMP(30)
      COMMON /POST/
     $ UMAX,VMAX,UVMAX,WMAX,TTMAX,TTMIN,WKMIN,WKMAX,WKDLT
     $  ,TIME,TIME0,DT,FINTIM
     $  ,XPLANE,YPLANE,ZPLANE
     $  ,XPOINT,YPOINT,ZPOINT
     $  ,XLINE(2),YLINE(2),ZLINE(2)
     $  ,VN(3),XV(3,3)

      COMMON /IPOST/
     $ ISTEP,NSTEPS,IDUMP
     $  ,IOSTEP,NHIS,NNHIS,NDUMPS,IINTEG,NINTEG
     $  ,IPLANE,JPLANE,KPLANE
     $  ,IPOINT,JPOINT,KPOINT,IELPNT,IHISPT
     $  ,IPLTYP,NXGRID,NXBAND,NXCONT,IROT

      COMMON /CPOST/ STRFLE,FLDFLE,FILENM,EXCODE(15)
      CHARACTER*2  EXCODE
      CHARACTER*17   STRFLE,FLDFLE,FILENM

      parameter (mobj=20,lobj=30000)
      common /rgeom/ zmouse,xobj(lobj),yobj(lobj),cobjs(6,mobj)
      common /igeom/ nobjs,iobjct,npts(mobj)
      common /lgeom/ ifobjs,ifobjg(mobj)
      common /cgeom/ ccobjs(mobj)
      logical ifobjs,ifobjg
      character*1 ccobjs

C     Window size
      COMMON /WIND/ WINDOWW,WINDOWH
      common /pixmp/ xpmn,ypmn,xpmx,ypmx
C
C     Arrays for face & edge manipulation
C
      COMMON /CFACES/ GROUP(6),SKPDAT(6,6)
     $               ,EFACE(6),EFACE1(6)
      INTEGER GROUP,SKPDAT,EFACE,EFACE1
C
C     New zipper and surface plot data
C
      PARAMETER (NSRFM=3*NXM*NELM)
      COMMON /QDPTS/ XYZQ(5,4,6,NELM)
      COMMON /QDPTI/ NLSTP,LIST(NSRFM),NMIRROR,LMIR(NSRFM)
     $                    ,LISW(NSRFM)
      COMMON /QDPTC/ NAMSEL,NAMELS(NSRFM),NAMESF(10)
      COMMON /SRFI/  NOPEN
      CHARACTER*16   NAMSEL,NAMELS,NAMESF
      COMMON /RCVOL/ VXPT1(3),VXPT2(3)
      COMMON /ICVOL/ NSLICE
      COMMON /CCVOL/ VPLANE
      COMMON /LCVOL/ IFCVOL
      CHARACTER*7 VPLANE
      LOGICAL IFCVOL

      common /forgotr/ cont_lev

      common /lat_r/ ulat1(3),ulat2(3),ulat3(3)
      common /lat_l/ if_lattice,if_lat_sph_cent
      logical        if_lattice,if_lat_sph_cent

      common /hexag/ long_edge(nelm)
     $             , neighb(6,nelm)

      parameter (lxybez=1000,lbez=100)
      common /cbezr/ xybez(2,lxybez)
      common /ibezr/ nbez,ibez(lbez)

c     These are hex27 arrays and also template arrays

      common /xyz27/ x27(27,nelm),y27(27,nelm),z27(27,nelm)
      common /ctemlate_r/ xtmpl_min, xtmpl_max
     $                  , ytmpl_min, ytmpl_max
     $                  , ztmpl_min, ztmpl_max
      common /ctemlate_i/ nel_templ,etemplate0,etemplate1
     $                  , nel_form ,eform0    ,eform1
      integer etemplate0,etemplate1,eform0,eform1

      integer icrn(4,6)
      save    icrn
      data    icrn / 1 , 2 , 6 , 5 
     $             , 2 , 3 , 7 , 6 
     $             , 3 , 4 , 8 , 7 
     $             , 4 , 1 , 5 , 8 
     $             , 1 , 2 , 3 , 4 
     $             , 5 , 6 , 7 , 8 /

      include 'devices.inc'
C     END OF BASICS.INC

