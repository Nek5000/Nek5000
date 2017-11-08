      subroutine makeq_aux

      include 'SIZE'
      include 'TOTAL'

      logical  if_conv_std

      if_conv_std = .true.
      if (ifmhd.and.ifaxis) if_conv_std = .false. ! conv. treated in induct.f

      if (if_conv_std) call makeuq ! this will zero bq first

      if(filterType.eq.2) call make_hpf
    
      return
      end
