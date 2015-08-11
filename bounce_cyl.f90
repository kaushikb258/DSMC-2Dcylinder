                 subroutine bounce_cylinder(x,y,u,v,dt, & 
                  xcen,ycen,h,xnew,ynew,unew,vnew)

         implicit none

                real(kind=8), parameter :: pi = acos(-1.0d0)
                 
                real(kind=8) :: x,y,u,v,dt, & 
                  xcen,ycen,h,xnew,ynew,unew,vnew
                real(kind=8) :: xold, yold 
                real(kind=8) :: rr, tt, r1
                real(kind=8) :: xc, yc, theta, ur, ut
                real(kind=8) :: xin, yin, xout, yout, xmid, ymid
                integer :: iter

 
            ! one timestep ago
            xold = x - u*dt
            yold = y - v*dt
              


            ! find (xc, yc) -> coordinates on circumference
            ! bisection method
              
            xout = xold
            yout = yold
            xin = x
            yin = y

            iter = 0

23          continue             
 
            xmid = 0.5d0*(xout+xin)
            ymid = 0.5d0*(yout+yin)
 
            r1 = sqrt((xmid-xcen)**2.0d0 + (ymid-ycen)**2.0d0)

            if(r1.lt.h) then
             xin = xmid
             yin = ymid
            else 
             xout = xmid
             yout = ymid
            endif 

            iter = iter + 1
            if(iter.gt.1000) then
             print*, 'convergence failed ', iter 
             stop
            endif

            if(abs(r1-h)/h.gt.1.0d-6) then
             goto 23
            endif

            xc = xmid
            yc = ymid



            rr = sqrt((xc-xcen)**2.0d0 + (yc-ycen)**2.0d0)
            if(abs((rr-h)/h).ge.1.0d-4) then   
             print*, 'bug in cyl ', rr, h, x, y
             print*, 'cyl ', sqrt((x-xcen)**2.0d0 + (y-ycen)**2.0d0)
             stop
            endif


            ! find the excess time, tt
            r1 = sqrt((xc - x)**2.0d0 + (yc - y)**2.0d0)  
            tt = r1/sqrt(u*u + v*v)
            if(tt.gt.1.1d0*dt) then
             print*, 'bug in bounce cyl ', tt, dt, xold, yold
             print*,  xc, yc, x, y, u, v
             stop
            endif  


            theta = atan(abs(yc - ycen)/abs(xc - xcen))
            ! if first quadrant, theta is OK
            if(xc.lt.xcen.and.yc.gt.ycen) then
             ! second quadrant
             theta = pi - theta 
            else if(xc.lt.xcen.and.yc.lt.ycen) then    
             ! third quadrant
             theta = theta + pi
            else if(xc.gt.xcen.and.yc.lt.ycen) then
             ! fourth quadrant
             theta = 2.0d0*pi - theta 
            endif


            ur = cos(theta)*u + sin(theta)*v   
            ut = -sin(theta)*u + cos(theta)*v   

            ! reverse ur; keep ut the same
            ur = -ur
              

            ! compute new u and v
            unew = cos(theta)*ur - sin(theta)*ut
            vnew = sin(theta)*ur + cos(theta)*ut 


            ! compute new x and y
            xnew = xc + unew*tt 
            ynew = yc + vnew*tt 



              return
                 end subroutine
