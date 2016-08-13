      subroutine fcvfun_usr(ydot,jstart)
c                                               
c     user specific contributions to rhs     
c                                               
      include 'SIZE'                            
      include 'TOTAL'                          
                                               
      real ydot(1)
                                                
      return
      end
c----------------------------------------------------------------------
      subroutine cvunpack_usr(y,jstart)
c
c     copy the cvode solution (y) back to the internal nek array (t)
c
      include 'SIZE'
      include 'TOTAL'
      include 'CVODE'

      real y(1)

      return
      end
c----------------------------------------------------------------------
      subroutine cvpack_usr(y,jstart)
c
c     copy the internal nek array (t) to the cvode solution (y)
c
      include 'SIZE'
      include 'TOTAL'
      include 'CVODE'

      real y(1)

      return
      end
