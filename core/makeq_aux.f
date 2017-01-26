      subroutine makeq_aux

      include 'SIZE'
      include 'TOTAL'

      logical  ifturb,if_conv_std

      if_conv_std = .true.
      if (ifmhd.and.ifaxis) if_conv_std = .false. ! conv. treated in induct.f

      call whatfld (ifturb)
      if (ifturb) call maketq ! zero bq
      if (.not.ifturb .and. if_conv_std)  call makeuq !zero bq

      call addhpf_q  ! Add high pass filtered fld for passive scalar   
   
      return
      end
