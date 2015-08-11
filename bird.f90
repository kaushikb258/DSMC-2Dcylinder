!----------------------------------------------------------------
! DSMC SIMULATION OF SUPERSONIC FLOW PAST CYLINDER
! AUTHOR: Kaushik Balakrishnan, PhD
! Email: kaushikb258@gmail.com
!----------------------------------------------------------------


                       program dsmc2d

        implicit none

        ! Boltzmann constant
        real(kind=8), parameter :: kb = 1.38d-23

        integer, parameter :: maxnumpar = 1000000
        integer, parameter :: maxnumparcell = 10000

        integer, parameter :: timesteps = 10000

        real(kind=8), dimension(1:maxnumpar,1:2) :: x, xold 
        real(kind=8), dimension(1:maxnumpar,1:2) :: v, vold 

        ! Geometry
        real(kind=8), parameter :: xmin = 0.0d0 
        real(kind=8), parameter :: xmax = 1.0d0 
        real(kind=8), parameter :: ymin = 0.0d0 
        real(kind=8), parameter :: ymax = 0.5d0 
        real(kind=8), parameter :: h = 0.05d0 
        real(kind=8), parameter :: xcen = 0.5d0 
        real(kind=8), parameter :: ycen = 0.25d0 

        integer, parameter :: imax = 200
        integer, parameter :: jmax = 100

        real(kind=8), parameter :: dx = (xmax-xmin)/dble(imax) 
        real(kind=8), parameter :: dy = (ymax-ymin)/dble(jmax) 
        real(kind=8), parameter :: dz = dx
        real(kind=8), parameter :: vcell = dx*dy*dz

        ! time step   
        real(kind=8), parameter :: dt = 1.0d-7 

        ! molecule mass and dia
        real(kind=8), parameter :: mass = 6.6d-26  
        real(kind=8), parameter :: dia = 4.0d-10  

        real(kind=8), parameter :: pi = acos(-1.0d0)   
        real(kind=8), parameter :: crossarea = pi*dia*dia/4.0d0  

        ! orifice conditions (mks units) 
        real(kind=8), parameter :: rho_inject = 1.0d-8  
        real(kind=8), parameter :: v_inject = 2000.0d0  
        real(kind=8), parameter :: T_inject = 300.0d0  
        real(kind=8), parameter :: vrms = sqrt(3.0d0*kb*T_inject/mass)  

        ! number of real molecules per simulated molecule
        real(kind=8), parameter :: np = 1.0d9  

        integer, dimension(1:maxnumpar) :: irecycle
        integer :: nrecycle 

        integer :: numpar, n, ninject, i, j, ii, jj, ncoll, k, npcellmax
        real(kind=8) ::  mdot, area, r, vmax, vp, nc, pcoll, ke
        real(kind=8) :: vcm(2), dv(2), lam
         
        integer, dimension(1:imax,1:jmax) :: npcell        
        integer, dimension(1:imax,1:jmax,1:maxnumparcell) :: icell    

        real(kind=8), dimension(1:imax,1:jmax) :: density, uvel, vvel
        real(kind=8), dimension(1:imax,1:jmax) :: temperature, pressure 

        integer, parameter :: nprim = 5 
        real(kind=8), dimension(1:imax,1:jmax,1:nprim) :: prim
        real(kind=8), dimension(1:imax,1:jmax,1:3) :: xyz

        integer, parameter :: output_unit = 29
        character (len=100) title, filename
        integer :: iout
        real(kind=8) :: xnew, ynew, unew, vnew, rr



        iout = 0
        numpar = 0
        x = 0.0d0
        v = 0.0d0  


        area = (ymax-ymin)*dz 
        mdot = rho_inject*area*v_inject  
 
        print*, 'vrms: ', vrms
        call sleep(3) 

        if(dx.ne.dy) then 
         print*, 'dx, dy: ', dx, dy  
         stop
        endif


        do i = 1, imax
        do j = 1, jmax
         xyz(i,j,1) = (dble(i)-0.5d0)*dx 
         xyz(i,j,2) = (dble(j)-0.5d0)*dy
         xyz(i,j,3) = 0.0d0  
        enddo   
        enddo   



        do n = 1, timesteps 
 
! inject particles
        call random_number(r)
        ninject = int(mdot*dt/mass/np + r)          

        if(mod(n,10).eq.0) then
         print*, '---------------------------' 
         print*, 'iteration number: ', n
         print*, 'ninject: ', ninject 
        endif

        do i = 1, ninject
          numpar = numpar + 1
          x(numpar,1) = xmin + 1.0d-6
          call random_number(r)
          x(numpar,2) = r*(ymax-ymin)
          call bm(r)
          v(numpar,1) = v_inject + vrms*r 
          call bm(r)
          v(numpar,2) = vrms*r 
        enddo  



         if(numpar.gt.maxnumpar) then
          print*, 'no more particles available ', numpar, maxnumpar
          stop
         endif



! move particles
        do i = 1, numpar
         x(i,1) = x(i,1) + v(i,1)*dt  
         x(i,2) = x(i,2) + v(i,2)*dt  
        enddo


! apply BC
        nrecycle = 0
        irecycle = 0
        do i = 1, numpar

         ! left
         if(x(i,1).lt.xmin) then
          x(i,1) = xmin + 1.0d-6
          v(i,1) = abs(v(i,1))
         endif

         ! bottom
         if(x(i,2).lt.ymin) then
          ! periodic BC 
          x(i,2) = ymax - abs(x(i,2)-ymin)
         endif

         ! top
         if(x(i,2).gt.ymax) then
          ! periodic 
          x(i,2) = ymin + abs(x(i,2)-ymax)
         endif

         ! right
         if(x(i,1).gt.xmax) then
          x(i,1) = xmax-1.0d-6
          v(i,1) = -abs(v(i,1)) 

          ! make sure this particle is not already in the recycle list
           jj = 0
           do j = 1, nrecycle
            if(i.eq.irecycle(j)) then
             jj = 1
            endif
           enddo
           if(jj.eq.0) then
            nrecycle = nrecycle + 1
            irecycle(nrecycle) = i
           endif
         endif


!-------------
         ! cylinder

         rr = sqrt( (xcen-x(i,1))**2.0d0 + (ycen-x(i,2))**2.0d0 )
         if(rr.le.h) then
          ! bounce this particle from the cylinder surface
          call bounce_cylinder(x(i,1),x(i,2),v(i,1),v(i,2),dt, &
                  xcen,ycen,h,xnew,ynew,unew,vnew) 
          x(i,1) = xnew  
          x(i,2) = ynew  
          v(i,1) = unew  
          v(i,2) = vnew  
         endif  

!-------------

        enddo


!----------------------------------
       ! check to ensure no particles remain inside cylinder
 
       do i = 1, numpar  
        rr = sqrt( (xcen-x(i,1))**2.0d0 + (ycen-x(i,2))**2.0d0 )
         if(rr.le.h) then
           print*, 'particle inside cyl ', n, i, rr, h 
           stop
         endif
       enddo     

!----------------------------------

        if(mod(n,10).eq.0) then
         print*, '# particles in the system: ', numpar
         print*, '# particles recycled in this time step: ', nrecycle
         print*, 'irecycle ', irecycle(1:nrecycle)
        endif


 
! book-keeping
       npcell = 0
       icell = 0          
       do i = 1, numpar
        ii = int(x(i,1)/dx)+1  
        jj = int(x(i,2)/dy)+1  
        if(ii.gt.imax.or.jj.gt.jmax) then
          print*, 'error ', ii, jj, i, x(i,:) 
          stop
        endif 
        npcell(ii,jj) = npcell(ii,jj) + 1

        if(npcell(ii,jj).ge.maxnumparcell) then
         print*, 'too many particles in a cell ', ii, jj, npcell(ii,jj) 
         stop
        endif
  
        icell(ii,jj,npcell(ii,jj)) = i 
       enddo


       npcellmax = maxval(npcell)
       lam = vcell/(sqrt(2.0d0)*dble(npcellmax)*np*pi*dia*dia)  
       if(dx/lam.gt.0.3d0.or.dx/lam.lt.0.0d0) then
         print*, 'mean free path ', dx, lam, npcellmax 
         stop
       endif
       if(mod(n,10).eq.0) then
        print*, 'dx/lam: ', dx/lam
        print*, 'Knudsen number: ', lam/h
        print*, '--------------------------'   
       endif


! collide particles (hard sphere model)
       do i = 1, imax
        do j = 1, jmax

       ! find vmax for this cell
       vmax = 0.0d0  
       do ii = 1, npcell(i,j)
        vp = sqrt(v(icell(i,j,ii),1)**2.0d0 + v(icell(i,j,ii),2)**2.0d0)
        vmax = max(vmax,vp)
       enddo        

       if(vmax*dt/dx.gt.0.5d0) then
        print*, 'time step error ', dx, vmax, dt 
        stop
       endif 


       ! find number of collisions for this cell
       nc = dble(npcell(i,j)*(npcell(i,j)-1))*crossarea*dt/2.0d0/vcell  
       call random_number(r)
       ncoll = int(nc + r) 
       
       do k = 1, ncoll 
 
        ! randomly choose 2 particles
124     continue
        call random_number(r)
        ii = int(r*dble(npcell(i,j)))
        if(ii.eq.0) then
          goto 124
        endif
        ii = icell(i,j,ii) 
          
125     continue
        call random_number(r)
        jj = int(r*dble(npcell(i,j)))
        if(jj.eq.0) then
          goto 125
        endif
        jj = icell(i,j,jj) 
          

        ! reject self-collisions
        if(ii.eq.jj) then
         goto 124
        endif
 
        ! now try to collide ii and jj particles
        vp = (v(ii,1)-v(jj,1))**2.0d0 + (v(ii,2)-v(jj,2))**2.0d0  
        vp = sqrt(vp)

        ! probability of collision 
        pcoll = vp/vmax

        call random_number(r)
        if(r.le.pcoll) then
         ! we have a collision
         vcm(1) = (v(ii,1) + v(jj,1))/2.0d0  
         vcm(2) = (v(ii,2) + v(jj,2))/2.0d0  
         dv(1) = v(ii,1) - v(jj,1)
         dv(2) = v(ii,2) - v(jj,2)
         v(ii,1:2) = vcm(1:2) + 0.5d0*dv(1:2)
         v(jj,1:2) = vcm(1:2) - 0.5d0*dv(1:2)
        endif 

       enddo


        enddo
       enddo 




! discard the recycled particles
        if(nrecycle.gt.0) then
          xold = x
          vold = v
          x = 0.0d0
          v = 0.0d0

          if(nrecycle.ge.numpar) then
           print*, 'recycle error ', nrecycle, numpar
           stop
          endif

          ii = 0 
          do i = 1, numpar

          jj = 0 
          do j = 1, nrecycle         
            if(i.eq.irecycle(j)) then
             jj = 1
            endif
          enddo

          ! jj = 0 => not in recycle list
          ! jj = 1 => present in recycle list

            if(jj.eq.0) then
             ii = ii + 1
             x(ii,1:2) = xold(i,1:2)
             v(ii,1:2) = vold(i,1:2)
            endif

          enddo 

          if(ii.ne.numpar-nrecycle) then
           print*, 'bug in recycling ', ii, numpar, nrecycle, numpar-nrecycle
           print*, 'irecycle ', irecycle(1:nrecycle)
           stop
          endif

        numpar = ii  
        nrecycle = 0
        irecycle = 0

        endif 




! discard the recycled particles
!        if(nrecycle.gt.0) then
!          do i = 1, nrecycle
!            x(irecycle(i):numpar-1,1:2) = x(irecycle(i)+1:numpar,1:2) 
!            v(irecycle(i):numpar-1,1:2) = v(irecycle(i)+1:numpar,1:2) 
!            numpar = numpar - 1
!            do j = i+1, nrecycle
!             irecycle(j) = irecycle(j) - 1
!            enddo
!          enddo 
!        endif 



!--------------------------
            if(mod(n,1000).eq.0) then
        iout = iout + 1
! output data

       density = 0.0d0
       uvel = 0.0d0
       vvel = 0.0d0 
       temperature = 0.0d0 
       pressure = 0.0d0

       do i = 1, numpar
        ii = int(x(i,1)/dx)+1  
        jj = int(x(i,2)/dy)+1  
        density(ii,jj) = density(ii,jj) + np*mass
        uvel(ii,jj) = uvel(ii,jj) + v(i,1)
        vvel(ii,jj) = vvel(ii,jj) + v(i,2)
        
        ke = v(i,1)**2.0d0 + v(i,2)**2.0d0 
        temperature(ii,jj) = temperature(ii,jj) + ke 
       enddo

       do i = 1, imax
        do j = 1, jmax
          if(npcell(i,j).gt.0) then 
          density(i,j) = density(i,j)/vcell
          uvel(i,j) = uvel(i,j)/dble(npcell(i,j)) 
          vvel(i,j) = vvel(i,j)/dble(npcell(i,j))

          temperature(i,j) = temperature(i,j)/dble(npcell(i,j))
          ke = uvel(i,j)**2.0d0 + vvel(i,j)**2.0d0 
          temperature(i,j) = temperature(i,j) - ke
          temperature(i,j) = temperature(i,j)*mass/3.0d0/kb
 
          pressure(i,j) = dble(npcell(i,j))*np*kb*temperature(i,j)/vcell
          endif
        enddo  
       enddo

     
       prim(1:imax,1:jmax,1) = density(1:imax,1:jmax)  
       prim(1:imax,1:jmax,2) = uvel(1:imax,1:jmax)  
       prim(1:imax,1:jmax,3) = vvel(1:imax,1:jmax)  
       prim(1:imax,1:jmax,4) = temperature(1:imax,1:jmax)  
       prim(1:imax,1:jmax,5) = pressure(1:imax,1:jmax)

       title = 'dsmc2d'
       write(filename,'("output/output_",I3.3,".vtk")'),iout 

       print*, 'writing output file '

       call vtk_write(output_unit,filename,title,imax,jmax,xyz,nprim,prim)


            endif
!--------------------------

        enddo

                      end program

!------------------------------------------------
! Box Muller distribution

         subroutine bm(r)
    
       implicit none

       real(kind=8) :: r1, r2, r, pi

       pi = acos(-1.0d0)

       call random_number(r1) 
       call random_number(r2) 

       r = sqrt(-2.0d0*log(r1))*cos(2.0d0*pi*r2)

         return
         end subroutine 


!------------------------------------------------
