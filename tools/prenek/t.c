         call draw_dom_bdry
      call menu(xmouse,ymouse,button,'Edit Mesh')  ! Prompt for input
         call prs('No help yet.$')
         call get_domain
         call mv_dom_vtx
         call mv_dom_vtx
         call dom_delete
      call refresh
      call drgrid
         call drawel(e)
      call prs('Click points to define domain, menu area to close$.')
         call mouse(xmouse,ymouse,button)
         call prs(string)
            if (k.eq.1) call move(xmouse,ymouse)
            if (k.gt.1) call draw(xmouse,ymouse)
              call prsis(' You have < 3 points. Delete domain? (y/n)$')
              call res(ans,1)
           call prs('Click to define domain, menu area to close$.')
      call draw(xdom(1),ydom(1))
