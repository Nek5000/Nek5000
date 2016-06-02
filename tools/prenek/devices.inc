C     Devices.inc
      COMMON /DEVS/ IF4107,IFLASE,IFLXY,IFGIGI,IFIRIS,IFMASS,IFPOST
     $             ,IFSUN,IFGKS,IFVMS,IFIBM,IF3D,IFNEK2,IFREAD,IFCEN
     $             ,IFNOSEG,IF12COL,IFBWGKS,IFPOSTS,IFXWIN
     $             ,IFGRID,IFGRDC,IFGRDP,IFGRDH
     $             ,IFHARD,IFHARD0,IFSAVE,IFFULL
     $             ,ifrevbk
      LOGICAL       IF4107,IFLASE,IFLXY,IFGIGI,IFIRIS,IFMASS,IFPOST
     $             ,IFSUN,IFGKS,IFVMS,IFIBM,IF3D,IFNEK2,IFREAD,IFCEN
     $             ,IFNOSEG,IF12COL,IFBWGKS,IFPOSTS,IFXWIN
     $             ,IFGRID,IFGRDC,IFGRDP,IFGRDH
     $             ,IFHARD,IFHARD0,IFSAVE,IFFULL
     $             ,ifrevbk
      COMMON/INITDT/ ESC,NULL
      CHARACTER*1    ESC,NULL
C
      COMMON/MODES/ IFDEMO,IFLEARN,ifgraf
      LOGICAL       IFDEMO,IFLEARN,ifgraf
      COMMON/GRIDS/ GRIDDX,GRIDDY,GRIDR,GRIDA,GRIDXP,GRIDYP
     $             ,XPHY0,YPHY0
     $             ,xhex(6),yhex(6)
      COMMON/IGRID/ NGRID
      COMMON/SCALE/XFAC ,YFAC ,XZERO ,YZERO
     $            ,XFACO,YFACO,XZEROO,YZEROO
      COMMON /CLIP/ XMOVEC,YMOVEC,WT,WB,WL,WR,WFR
C     Save 3D move and draws
      PARAMETER (NDRMAX=10000)
      COMMON /ISAVED/ IDRAWD,IDRCNT,LSTCLR
      COMMON /RSAVED/ XYZD(5,NDRMAX)
c
c
c
c     .plt file pointers
      common /fplotc/ fplotnam
      character*17    fplotnam
      common /fploti/ iplotnam,iplotnum

      real elec,elew,eleb,eledy ! Elevator coords for BC menu
      save elec,elew,eleb,eledy
      data elec,elew,eleb,eledy / 0.95 , 0.02 , 0.10 , 0.10 /


c     common /cfilli/ icfill,ifill,nfill
c     common /cfillr/ xfill(210),yfill(210)

c     End of Devices.inc***********
