***************************************************************
** The following code is for simulating the 2D flow past
*  a forward facing step using the high order accurate
*  WENO schemes, developed by Guang-Shan Jiang and Chi-Wang Shu
*  (see the paper on JCP, v126, pp. 202-228(1996)).
*
** The code is in Fortran 77 with single precision.
*  Compile with "f77 -r8" to achieve double precision.
* 
** Any comments on the performance of the code and the schemes
*  are very welcome.
*
** The authors are happy to provide this code to any potential
*  users. However, any redistribution of this code should 
*  carry this preamble.
*
** Please credit the authors for any work facilitated by this
*  code or the WENO schemes.
***************************************************************

      program main

      parameter(md=4,mn=4,nx=122,ny=39,nd=122)
      parameter(nxm=nx+md,nym=ny+md,ndm=nd+md)
      common/const/tend,cfl,epweno,gamma,gm1
      common/grid/dx,cdx,x(-md:nxm),dy,cdy,y(-md:nym),dt,em,tnum
      common/sol/ u(-md:nxm,-md:nym,mn),f(-md:nxm,-md:nym,mn)
     &           ,u0(-md:nxm,-md:nym,mn),uc(-md:nxm,-md:nym,mn)
      common/zerod/ul(-md:ndm,mn),ur(-md:ndm,mn)
     &            ,ut(-md:ndm,mn),ub(-md:ndm,mn),um(-md:ndm,mn)
      common/step/nxmid,nxl(-md:nym),nxr(-md:nym)
     &           ,nymid,nyl(-md:nxm),nyr(-md:nxm)

      call init

* Runge-Kutta time stepping:
      tnum=0.
      istop=0

        write(*,*) '4th order Runge-kutta in time'

      do 400 nt=1, 100000

      call bcx
      call enox
      call bcy
      call enoy
      dt=cfl/em/(cdx+cdy)
      if((tnum+dt).ge.tend) then
           istop=1
           dt=tend-tnum
      endif
      dt2=0.5*dt
      tnum=tnum+dt
      do 410 m=1,mn
      do 410 j=0,ny
      do 410 i=nxl(j),nxr(j)
        u0(i,j,m)=u(i,j,m)
         u(i,j,m)=u(i,j,m)+dt2*f(i,j,m)
410     continue

      call bcx
      call enox
      call bcy
      call enoy
      do 420 m=1,mn
      do 420 j=0,ny
      do 420 i=nxl(j),nxr(j)
        uc(i,j,m)=-u0(i,j,m)+u(i,j,m)
         u(i,j,m)= u0(i,j,m)+dt2*f(i,j,m)
420     continue

      call bcx
      call enox
      call bcy
      call enoy
      do 430 m=1,mn
      do 430 j=0,ny
      do 430 i=nxl(j),nxr(j)
        uc(i,j,m)=uc(i,j,m)+2.*u(i,j,m)
         u(i,j,m)=u0(i,j,m)+dt*f(i,j,m)
430     continue

      call bcx
      call enox
      call bcy
      call enoy
        do 440 m=1,mn
      do 440 j=0,ny
      do 440 i=nxl(j),nxr(j)
        u(i,j,m)=(1./3.)*(uc(i,j,m)+u(i,j,m)+dt2*f(i,j,m))
440      continue

      if(istop.eq.1) goto 999
400     continue
       
999    continue

       write(*,*) '  nt=',nt,'  t=', tnum

       call save

1      format(1x,'CPU time=',f12.4)

       stop
       end


      subroutine init

      parameter(md=4,mn=4,nx=122,ny=39,nd=122)
      parameter(nxm=nx+md,nym=ny+md,ndm=nd+md)
      common/const/tend,cfl,epweno,gamma,gm1
      common/grid/dx,cdx,x(-md:nxm),dy,cdy,y(-md:nym),dt,em,tnum
      common/sol/ u(-md:nxm,-md:nym,mn),f(-md:nxm,-md:nym,mn)
     &           ,u0(-md:nxm,-md:nym,mn),uc(-md:nxm,-md:nym,mn)
      common/zerod/ul(-md:ndm,mn),ur(-md:ndm,mn)
     &            ,ut(-md:ndm,mn),ub(-md:ndm,mn),um(-md:ndm,mn)
      common/step/nxmid,nxl(-md:nym),nxr(-md:nym)
     &           ,nymid,nyl(-md:nxm),nyr(-md:nxm)

* Setup parameters:

* Ratio of specific heats:
      gamma = 1.4
      gm1   = gamma - 1.0
* WENO constant
      epweno = 1.e-6
* Final time
      tend = 4.0
* CFL number
      cfl = 0.6

      xleft=0.0
      xright=3.0
      yleft=0.0
      yright=1.0
      xmid=0.6
      ymid=0.2

      write(*,*) ' '
      write(*,*) 'Flow past a forward facing step'
      write(*,*) ' '
      write(*,*) 'xleft=',xleft, '  xright=', xright
      write(*,*) 'yleft=',yleft, '  yright=', yright
      write(*,*) 'tend=', tend,  '  cfl=', cfl
      write(*,*) ' '
      write(*,*) 'constant initial condition'
      write(*,*) 'density x-velocity y-velocity pressure'
      write(*,*) '  1.4       3          0         1'
      write(*,*) '    '

      dx=(xright-xleft)/(nx+0.5)
      dy=(yright-yleft)/(ny+1)
      ii=int((xmid-xleft)/dx+1.e-4)      
      jj=int((ymid-yleft)/dy+1.e-4)-1      
      nxmid=ii
      nymid=jj
      cdx=1./dx
      cdy=1./dy
      
      do 2 i=-md,nxm 
        x(i)=xleft+i*dx
2     continue   
      do 3 i=-md,nym  
        y(i)=yleft+(i+0.5)*dy
3     continue   
      do 4 j=-md,nym
      if(j.le.jj) then
        nxl(j)=0
        nxr(j)=ii
        else
        nxl(j)=0
        nxr(j)=nx
        endif
4     continue
      do 5 i=-md,nxm
      if(i.le.ii) then
        nyl(i)=0
        nyr(i)=ny
      else
        nyl(i)=jj+1
        nyr(i)=ny
      endif
5     continue

c give initial conditions

      den=1.4
      vex=3.0
      vey=0.
      pre=1.0
      do 10 j=-md,nym                 
      do 10 i=-md,nxm                 
        u(i,j,1)=den                 
        u(i,j,2)=den*vex              
        u(i,j,3)=den*vey              
        u(i,j,4)=pre/gm1+0.5*den*(vex**2+vey**2)
10    continue

      do 20 j=0,ny
        ul(j,1)=den
        ul(j,2)=den*vex
        ul(j,3)=den*vey
        ul(j,4)=pre/gm1+0.5*den*(vex**2+vey**2)
        ur(j,1)=den
        ur(j,2)=den*vex
        ur(j,3)=den*vey
        ur(j,4)=pre/gm1+0.5*den*(vex**2+vey**2)
20     continue

      return  
      end   


      subroutine bcy

      parameter(md=4,mn=4,nx=122,ny=39,nd=122)
      parameter(nxm=nx+md,nym=ny+md,ndm=nd+md)
      common/const/tend,cfl,epweno,gamma,gm1
      common/grid/dx,cdx,x(-md:nxm),dy,cdy,y(-md:nym),dt,em,tnum
      common/sol/ u(-md:nxm,-md:nym,mn),f(-md:nxm,-md:nym,mn)
     &           ,u0(-md:nxm,-md:nym,mn),uc(-md:nxm,-md:nym,mn)
      common/zerod/ul(-md:ndm,mn),ur(-md:ndm,mn)
     &            ,ut(-md:ndm,mn),ub(-md:ndm,mn),um(-md:ndm,mn)
      common/step/nxmid,nxl(-md:nym),nxr(-md:nym)
     &           ,nymid,nyl(-md:nxm),nyr(-md:nxm)

c bottom boundary(reflective)

      do 3 j=1,md
      do 3 i=0,nxmid
        u(i,-j,1)= u(i,j-1,1)
        u(i,-j,2)= u(i,j-1,2)
        u(i,-j,3)=-u(i,j-1,3)
        u(i,-j,4)= u(i,j-1,4)
3     continue                                                          
      do 5 j=1,md
      do 5 i=nxmid+1,nx
        u(i,nymid+1-j,1)= u(i,nymid+j,1)
        u(i,nymid+1-j,2)= u(i,nymid+j,2)
        u(i,nymid+1-j,3)=-u(i,nymid+j,3)
        u(i,nymid+1-j,4)= u(i,nymid+j,4)
5     continue                                                          

c top boundary(reflective)

      do 4 j=1,md
      do 4 i=0,nx
        u(i,ny+j,1)= u(i,ny+1-j,1)
        u(i,ny+j,2)= u(i,ny+1-j,2)
        u(i,ny+j,3)=-u(i,ny+1-j,3)
        u(i,ny+j,4)= u(i,ny+1-j,4)
4     continue                                                          

      return
      end


      subroutine bcx

      parameter(md=4,mn=4,nx=122,ny=39,nd=122)
      parameter(nxm=nx+md,nym=ny+md,ndm=nd+md)
      common/const/tend,cfl,epweno,gamma,gm1
      common/grid/dx,cdx,x(-md:nxm),dy,cdy,y(-md:nym),dt,em,tnum
      common/sol/ u(-md:nxm,-md:nym,mn),f(-md:nxm,-md:nym,mn)
     &           ,u0(-md:nxm,-md:nym,mn),uc(-md:nxm,-md:nym,mn)
      common/zerod/ul(-md:ndm,mn),ur(-md:ndm,mn)
     &            ,ut(-md:ndm,mn),ub(-md:ndm,mn),um(-md:ndm,mn)
      common/step/nxmid,nxl(-md:nym),nxr(-md:nym)
     &           ,nymid,nyl(-md:nxm),nyr(-md:nxm)
      real utmp(-md:ndm,4)

c treat the singularity at the corner

      den=u(nxmid,nymid,1)
      vex=u(nxmid,nymid,2)/den
      vey=u(nxmid,nymid,3)/den
      q2=vex*vex+vey*vey
      pre=gm1*(u(nxmid,nymid,4)-0.5*den*q2)
      enth=(u(nxmid,nymid,4)+pre)/den
      ent=pre/abs(den)**gamma

      do 10 j=0,1
      do 10 i=nxmid+1,nxmid+4-2*j

        jj=j+nymid+1

        den=u(i,jj,1)
        vex=u(i,jj,2)/den
        vey=u(i,jj,3)/den
        q2=vex*vex+vey*vey
        pre=gm1*(u(i,jj,4)-0.5*den*q2)

        den=(abs(pre/ent))**(1./gamma)
        qq2=abs(2.*(enth-gamma*pre/den/gm1))
        t0=sqrt(qq2/q2) 
        vex=vex*t0
        vey=vey*t0
        u(i,jj,1)=den
        u(i,jj,2)=den*vex
        u(i,jj,3)=den*vey
        u(i,jj,4)=pre/gm1+0.5*den*qq2
10    continue

c left boundary(osher inflow b.c.)

      do 1 m=1,4                                                      
      do 1 j=0,ny
        utmp(j,m)=u(0,j,m)
1     continue                                                          

      call osher(0,ny,ul,utmp)

      do 7 m=1,4                                                      
      do 7 i=0,md
      do 7 j=0,ny
        u(-i,j,m)=um(j,m)
7     continue                                                          

c right boundary(part reflective, part osher outflow b.c.)

      do 6 i=1,md
      do 6 j=0,nymid
        u(nxmid+i,j,1)= u(nxmid+1-i,j,1)
        u(nxmid+i,j,2)=-u(nxmid+1-i,j,2)
        u(nxmid+i,j,3)= u(nxmid+1-i,j,3)
        u(nxmid+i,j,4)= u(nxmid+1-i,j,4)
6     continue                                                          

      do 2 m=1,4                                                      
      do 2 j=nymid+1,ny
        utmp(j,m)=u(nx,j,m)
2     continue                                                          

      call osher(nymid+1,ny,utmp,ur)

      do 8 m=1,4                                                      
      do 8 i=0,md
      do 8 j=nymid+1,ny
        u(nx+i,j,m)=um(j,m)
8     continue                                                          

      return
      end


      subroutine osher(loops,loope,vl,vr)

      parameter(md=4,mn=4,nx=122,ny=39,nd=122)
      parameter(nxm=nx+md,nym=ny+md,ndm=nd+md)
      common/const/tend,cfl,epweno,gamma,gm1
      common/grid/dx,cdx,x(-md:nxm),dy,cdy,y(-md:nym),dt,em,tnum
      common/sol/ u(-md:nxm,-md:nym,mn),f(-md:nxm,-md:nym,mn)
     &           ,u0(-md:nxm,-md:nym,mn),uc(-md:nxm,-md:nym,mn)
      common/zerod/ul(-md:ndm,mn),ur(-md:ndm,mn)
     &            ,ut(-md:ndm,mn),ub(-md:ndm,mn),um(-md:ndm,mn)
      common/step/nxmid,nxl(-md:nym),nxr(-md:nym)
     &           ,nymid,nyl(-md:nxm),nyr(-md:nxm)

      dimension vl(-md:ndm,4),vr(-md:ndm,4),um2(4),um3(4)

      do 10 j=loops,loope

      rho4=vl(j,1)
      u4=vl(j,2)/rho4   
      v4=vl(j,3)/rho4
      p4=gm1*(vl(j,4)-0.5*rho4*(u4*u4+v4*v4))
      c4=sqrt(gamma*p4/rho4)

      rho1=vr(j,1)
      u1=vr(j,2)/rho1
      v1=vr(j,3)/rho1
      p1=gm1*(vr(j,4)-0.5*rho1*(u1*u1+v1*v1))
      c1=sqrt(gamma*p1/rho1)

      c2=c1*(gm1/2.*(u4-u1)+(c4+c1))/(c1+c4*(p1/p4)**(gm1/2./gamma))
      p2=p1*(c2/c1)**(2.*gamma/gm1)
      p3=p2
      c3=c4*(p3/p4)**(gm1/2./gamma)
      u2=u1+2.*(c2-c1)/gm1
      u3=u2

      rho3=gamma*p3/c3**2
      rho2=gamma*p2/c2**2

      v3=v4
      v2=v1

      um3(1)=rho3
      um3(2)=rho3*u3
      um3(3)=rho3*v3
      um3(4)=p3/gm1+0.5*rho3*(u3*u3+v3*v3)

      um2(1)=rho2
      um2(2)=rho2*u2
      um2(3)=rho2*v2
      um2(4)=p2/gm1+0.5*rho2*(u2*u2+v2*v2)

      if(u4-c4.gt.0.) then
       um(j,1)=vl(j,1)
       um(j,2)=vl(j,2)
       um(j,3)=vl(j,3)
       um(j,4)=vl(j,4)
      else if(u2.gt.0.) then
       um(j,1)=um3(1)
       um(j,2)=um3(2)
       um(j,3)=um3(3)
       um(j,4)=um3(4)
      else if(u1+c1.gt.0.) then
       um(j,1)=um2(1)
       um(j,2)=um2(2)
       um(j,3)=um2(3)
       um(j,4)=um2(4)
      else
       um(j,1)=vr(j,1)
       um(j,2)=vr(j,2)
       um(j,3)=vr(j,3)
       um(j,4)=vr(j,4)
      end if

   10 continue
      
      return
      end

      subroutine save  

      parameter(md=4,mn=4,nx=122,ny=39,nd=122)
      parameter(nxm=nx+md,nym=ny+md,ndm=nd+md)
      common/const/tend,cfl,epweno,gamma,gm1
      common/grid/dx,cdx,x(-md:nxm),dy,cdy,y(-md:nym),dt,em,tnum
      common/sol/ u(-md:nxm,-md:nym,mn),f(-md:nxm,-md:nym,mn)
     &           ,u0(-md:nxm,-md:nym,mn),uc(-md:nxm,-md:nym,mn)
      common/zerod/ul(-md:ndm,mn),ur(-md:ndm,mn)
     &            ,ut(-md:ndm,mn),ub(-md:ndm,mn),um(-md:ndm,mn)
      common/step/nxmid,nxl(-md:nym),nxr(-md:nym)
     &           ,nymid,nyl(-md:nxm),nyr(-md:nxm)

c save the solution for graphing with Tecplot

      open(10,file='data')

      write(10,*)'   TITLE = "Flow Past a Forward Facing Step" '
      write(10,13)
      write(10,14) nxmid+1,nymid+2
13    format(3x,'VARIABLES = "X","Y","DENSITY","X-VELOCITY",'
     &                      ,'"Y-VELOCITY","PRESSURE","ENTROPY"')
14    format(3x,'ZONE T="Small Rect", I=',I3,' J=',I3)
      do 2 j=0,nymid+1
      do 2 i=0,nxmid
        den=u(i,j,1)
        xmt=u(i,j,2)
        ymt=u(i,j,3)
        vex=xmt/den
        vey=ymt/den
        pre=gm1*(u(i,j,4)-0.5*(xmt*vex+ymt*vey))
        ent=alog(abs(pre/abs(den)**gamma))
        write(10,12) x(i),y(j),den,vex,vey,pre,ent
2     continue    
      write(10,15) nx+1,(ny-nymid)
15    format(3x,'ZONE T="Large Rect", I=',I3,' J=',I3)
      do 3 j=nymid+1,ny
      do 3 i=0,nx
        den=u(i,j,1)
        xmt=u(i,j,2)
        ymt=u(i,j,3)
        vex=xmt/den
        vey=ymt/den
        pre=gm1*(u(i,j,4)-0.5*(xmt*vex+ymt*vey))
        ent=alog(abs(pre/abs(den)**gamma))
        write(10,12) x(i),y(j),den,vex,vey,pre,ent
3     continue    

      close(91)

12    format(7f11.6)    

      return
      end


      subroutine enox

      parameter(md=4,mn=4,nx=122,ny=39,nd=122)
      parameter(nxm=nx+md,nym=ny+md,ndm=nd+md)
      common/const/tend,cfl,epweno,gamma,gm1
      common/grid/dx,cdx,x(-md:nxm),dy,cdy,y(-md:nym),dt,em,tnum
      common/sol/ u(-md:nxm,-md:nym,mn),f(-md:nxm,-md:nym,mn)
     &           ,u0(-md:nxm,-md:nym,mn),uc(-md:nxm,-md:nym,mn)
      common/zerod/ul(-md:ndm,mn),ur(-md:ndm,mn)
     &            ,ut(-md:ndm,mn),ub(-md:ndm,mn),um(-md:ndm,mn)
      common/step/nxmid,nxl(-md:nym),nxr(-md:nym)
     &           ,nymid,nyl(-md:nxm),nyr(-md:nxm)
      real w(-md:nxm),h(-md:nxm),vx(-md:nxm),vy(-md:nxm)
     &    ,fp(-md:nxm,mn),uk(-md:nxm,mn),fh(-md:nxm,mn),am(mn)
     &    ,evl(-md:nxm,mn,mn),evr(-md:nxm,mn,mn)
     &    ,tf(-md:nxm,mn),tu(-md:nxm,mn)
     &    ,ff(-md:nxm,mn),gg(-md:nxm,mn,2),hh(-md:nxm,4,2)

* Use global Lax-Friedriches flux splitting:

*  Compute  f_x : 

              em=1.e-15

      do 100 j=0,ny

       am(1) = 1.e-15
       am(2) = 1.e-15
       am(4) = 1.e-15
       do 110 i=-md+nxl(j), nxr(j)+md
         den=u(i,j,1)
         xmt=u(i,j,2)
         ymt=u(i,j,3)
         eng=u(i,j,4)
         t0=1./den
         vex=xmt*t0
         vey=ymt*t0
         pre=gm1*(eng-0.5*(xmt*vex+ymt*vey))
         cvel=sqrt(abs(gamma*pre*t0))
         fp(i,1)=xmt
         fp(i,2)=pre+xmt*vex
         fp(i,3)=xmt*vey
         fp(i,4)=vex*(pre+eng)
         uk(i,1)=den
         uk(i,2)=xmt
         uk(i,3)=ymt
         uk(i,4)=eng
         w(i)=sqrt(abs(den))
         vx(i)=vex
         vy(i)=vey
         h(i)=(pre+eng)*t0
         am(1)=max(am(1),abs(vex-cvel))
         am(2)=max(am(2),abs(vex))
         am(4)=max(am(4),abs(vex+cvel))

110    continue
         am(1)=am(1)*1.1
         am(2)=am(2)*1.1
         am(3)=am(2)
         am(4)=am(4)*1.1
         em=max(em, max(am(1),am(4)))


       do 115 i=-1+nxl(j),nxr(j)  
* Compute e'vectors using Roe's average:
         ip=i+1
         t0=w(i)/(w(i)+w(ip)+1.e-15)
         t1=1.-t0
         vxm=t0*vx(i)+t1*vx(ip)
         vym=t0*vy(i)+t1*vy(ip)
         hm=t0*h(i)+t1*h(ip)
         qm=0.5*(vxm*vxm+vym*vym)
         cm=sqrt(abs(gm1*(hm-qm)))+1.e-15
         t0=vxm*cm
         evr(i,1,1)=1.0
         evr(i,1,2)=0.0
         evr(i,1,3)=1.0
         evr(i,1,4)=1.0
         evr(i,2,1)=vxm-cm
         evr(i,2,2)=0.0
         evr(i,2,3)=vxm
         evr(i,2,4)=vxm+cm
         evr(i,3,1)=vym
         evr(i,3,2)=1.0
         evr(i,3,3)=vym
         evr(i,3,4)=vym
         evr(i,4,1)=hm-t0
         evr(i,4,2)=vym
         evr(i,4,3)=qm
         evr(i,4,4)=hm+t0
         rcm=1./cm
         b1=gm1*rcm*rcm
         b2=qm*b1
         t0=vxm*rcm
         t1=b1*vxm
         t2=0.5*b1
         t3=b1*vym
         evl(i,1,1)=0.5*(b2+t0)
         evl(i,1,2)=-0.5*(t1+rcm)
         evl(i,1,3)=-0.5*t3
         evl(i,1,4)=t2
         evl(i,2,1)=-vym
         evl(i,2,2)=0.0
         evl(i,2,3)=1.0
         evl(i,2,4)=0.0
         evl(i,3,1)=1.-b2
         evl(i,3,2)=t1
         evl(i,3,3)=t3
         evl(i,3,4)=-b1
         evl(i,4,1)=0.5*(b2-t0)
         evl(i,4,2)=-0.5*(t1-rcm)
         evl(i,4,3)=-0.5*t3
         evl(i,4,4)=t2
115    continue

       do 120 m=1,mn
       do 120 i=-md+nxl(j),nxr(j)+md-1
        tf(i,m)=fp(i+1,m)-fp(i,m)
        tu(i,m)=uk(i+1,m)-uk(i,m)
120     continue

*  begin of the big loop

       do 125 m=1,mn
       do 125 i = -1, nx
125      ff(i,m) = 0.0

       do 160 m=1,mn

* Project the relevant first undivided differences into the 
* 'm'th characteristic field
       do 130 m1=1,mn
       do 130 i=-md+nxl(j),nxr(j)+md-1
       gg(i,m1,1)=0.5*(tf(i,m1)+am(m)*tu(i,m1))
       gg(i,m1,2)=gg(i,m1,1)-tf(i,m1)
130    continue

       do 140 m1=1,4
       k0=m1-3
       k1=3-m1
       do 140 i=-1+nxl(j),nxr(j)
         hh(i,m1,1)=evl(i,m,1)*gg(i+k0,1,1)+evl(i,m,2)*gg(i+k0,2,1)
     &             +evl(i,m,3)*gg(i+k0,3,1)+evl(i,m,4)*gg(i+k0,4,1)
       hh(i,m1,2)=evl(i,m,1)*gg(i+k1,1,2)+evl(i,m,2)*gg(i+k1,2,2)
     &             +evl(i,m,3)*gg(i+k1,3,2)+evl(i,m,4)*gg(i+k1,4,2)
140    continue

* Compute numerical flux in each characteristic field:

       do 150 m1=1,2
       do 150 i=-1+nxl(j),nxr(j)
           t1 = hh(i,1,m1) - hh(i,2,m1)
           t2 = hh(i,2,m1) - hh(i,3,m1)
           t3 = hh(i,3,m1) - hh(i,4,m1)
          tt1 = 13. * t1**2 + 3. * (   hh(i,1,m1) - 3*hh(i,2,m1) )**2
          tt2 = 13. * t2**2 + 3. * (   hh(i,2,m1) +   hh(i,3,m1) )**2
          tt3 = 13. * t3**2 + 3. * ( 3*hh(i,3,m1) -   hh(i,4,m1) )**2
           s1 =  1.0 / ( epweno + tt1 )**2
           s2 =  6.0 / ( epweno + tt2 )**2
           s3 =  3.0 / ( epweno + tt3 )**2
           t0 = 1. / ( s1 + s2 + s3 )
           s1 = s1 * t0
           s3 = s3 * t0
         ff(i,m) = ff(i,m) + ( s1*(t2-t1) + (0.5*s3-0.25)*(t3-t2) ) /3.
150    continue

160    continue

*  end of the big loop

* Project the numerical flux to the physical space:

        do 162 m=1,mn
        do 162 i=-1+nxl(j),nxr(j)
       fh(i,m)=(-fp(i-1,m)+7*(fp(i,m)+fp(i+1,m))-fp(i+2,m))/12.
     &          +evr(i,m,1)*ff(i,1)+evr(i,m,2)*ff(i,2)
     &          +evr(i,m,3)*ff(i,3)+evr(i,m,4)*ff(i,4)
162     continue

       do 166 m=1,mn
       do 166 i=nxl(j),nxr(j)
        f(i,j,m)=(fh(i-1,m)-fh(i,m))*cdx
166    continue
100    continue

      return
      end


      subroutine enoy

      parameter(md=4,mn=4,nx=122,ny=39,nd=122)
      parameter(nxm=nx+md,nym=ny+md,ndm=nd+md)
      common/const/tend,cfl,epweno,gamma,gm1
      common/grid/dx,cdx,x(-md:nxm),dy,cdy,y(-md:nym),dt,em,tnum
      common/sol/ u(-md:nxm,-md:nym,mn),f(-md:nxm,-md:nym,mn)
     &           ,u0(-md:nxm,-md:nym,mn),uc(-md:nxm,-md:nym,mn)
      common/zerod/ul(-md:ndm,mn),ur(-md:ndm,mn)
     &            ,ut(-md:ndm,mn),ub(-md:ndm,mn),um(-md:ndm,mn)
      common/step/nxmid,nxl(-md:nym),nxr(-md:nym)
     &           ,nymid,nyl(-md:nxm),nyr(-md:nxm)
      real w(-md:nym),h(-md:nym),vx(-md:nym),vy(-md:nym)
     &    ,fp(-md:nym,mn),uk(-md:nym,mn),fh(-md:nym,mn),am(mn)
     &    ,evl(-md:nym,mn,mn),evr(-md:nym,mn,mn)
     &    ,tf(-md:nym,mn),tu(-md:nym,mn)
     &    ,ff(-md:nym,mn),gg(-md:nym,mn,2),hh(-md:nym,4,2)

* Use global Lax-Friedriches flux splitting:

*  Compute  g_y :

      do 200 j=0,nx

       am(1) = 1.e-15
       am(2) = 1.e-15
       am(4) = 1.e-15
       do 210 i=-md+nyl(j), nyr(j)+md
          den=u(j,i,1)
          xmt=u(j,i,2)
          ymt=u(j,i,3)
          eng=u(j,i,4)
          t0=1./den
         vex=xmt*t0
         vey=ymt*t0
         pre=gm1*(eng-0.5*(xmt*vex+ymt*vey))
         cvel=sqrt(abs(gamma*pre*t0))
           fp(i,1)=ymt
           fp(i,2)=xmt*vey
           fp(i,3)=pre+ymt*vey
           fp(i,4)=vey*(pre+eng)
           uk(i,1)=den
           uk(i,2)=xmt
           uk(i,3)=ymt
           uk(i,4)=eng
         w(i)=sqrt(abs(den))
         vx(i)=vex
         vy(i)=vey
         h(i)=(pre+eng)*t0
           am(1)=max(am(1),abs(vey-cvel))
           am(2)=max(am(2),abs(vey))
           am(4)=max(am(4),abs(vey+cvel))
210     continue
           am(1)=am(1)*1.1
           am(2)=am(2)*1.1
           am(3)=am(2)
           am(4)=am(4)*1.1
           em=max(em,max(am(1),am(4)))

       do 215 i=-1+nyl(j),nyr(j)  
* Compute e'vectors using Roe's average:
          ip=i+1
          t0=w(i)/(w(i)+w(ip)+1.e-15)
          t1=1.-t0
         vxm=t0*vx(i)+t1*vx(ip)
         vym=t0*vy(i)+t1*vy(ip)
         hm=t0*h(i)+t1*h(ip)
         qm=0.5*(vxm*vxm+vym*vym)
         cm=sqrt(abs(gm1*(hm-qm)))+1.e-15
          t0=vym*cm
         evr(i,1,1)=1.0
         evr(i,1,2)=0.0
         evr(i,1,3)=1.0
         evr(i,1,4)=1.0
         evr(i,2,1)=vxm
         evr(i,2,2)=1.0
         evr(i,2,3)=vxm
         evr(i,2,4)=vxm
         evr(i,3,1)=vym-cm
         evr(i,3,2)=0.0
         evr(i,3,3)=vym
         evr(i,3,4)=vym+cm
         evr(i,4,1)=hm-t0
         evr(i,4,2)=vxm
         evr(i,4,3)=qm
         evr(i,4,4)=hm+t0
          rcm=1./cm
          b1=gm1*rcm*rcm
          b2=qm*b1
          t0=vym*rcm
          t1=b1*vym
          t2=0.5*b1
          t3=b1*vxm
         evl(i,1,1)=0.5*(b2+t0)
         evl(i,1,2)=-0.5*t3
         evl(i,1,3)=-0.5*(t1+rcm)
         evl(i,1,4)=t2
         evl(i,2,1)=-vxm
         evl(i,2,2)=1.0
         evl(i,2,3)=0.0
         evl(i,2,4)=0.0
         evl(i,3,1)=1.-b2
       evl(i,3,2)=t3
       evl(i,3,3)=t1
       evl(i,3,4)=-b1
       evl(i,4,1)=0.5*(b2-t0)
       evl(i,4,2)=-0.5*t3
       evl(i,4,3)=-0.5*(t1-rcm)
       evl(i,4,4)=t2
215    continue

       do 220 m=1,mn
       do 220 i=-md+nyl(j),nyr(j)+md-1
        tf(i,m)=fp(i+1,m)-fp(i,m)
        tu(i,m)=uk(i+1,m)-uk(i,m)
220     continue

*  begin of the big loop

       do 225 m=1,mn
       do 225 i = -1, ny
225      ff(i,m) = 0.0

       do 260 m=1,mn

* Project the relevant first undivided differences into the 
* 'm'th characteristic field
       do 230 m1=1,mn
       do 230 i=-md+nyl(j),nyr(j)+md-1
         gg(i,m1,1)=0.5*(tf(i,m1)+am(m)*tu(i,m1))
         gg(i,m1,2)=gg(i,m1,1)-tf(i,m1)
230    continue

       do 240 m1=1,4
         k0=m1-3
         k1=3-m1
       do 240 i=-1+nyl(j),nyr(j)
         hh(i,m1,1)=evl(i,m,1)*gg(i+k0,1,1)+evl(i,m,2)*gg(i+k0,2,1)
     &             +evl(i,m,3)*gg(i+k0,3,1)+evl(i,m,4)*gg(i+k0,4,1)
         hh(i,m1,2)=evl(i,m,1)*gg(i+k1,1,2)+evl(i,m,2)*gg(i+k1,2,2)
     &             +evl(i,m,3)*gg(i+k1,3,2)+evl(i,m,4)*gg(i+k1,4,2)
240    continue

* Compute numerical flux in each characteristic field:

       do 250 m1=1,2
       do 250 i=-1+nyl(j),nyr(j)
           t1 = hh(i,1,m1) - hh(i,2,m1)
           t2 = hh(i,2,m1) - hh(i,3,m1)
           t3 = hh(i,3,m1) - hh(i,4,m1)
          tt1 = 13. * t1**2 + 3. * (   hh(i,1,m1) - 3*hh(i,2,m1) )**2
          tt2 = 13. * t2**2 + 3. * (   hh(i,2,m1) +   hh(i,3,m1) )**2
          tt3 = 13. * t3**2 + 3. * ( 3*hh(i,3,m1) -   hh(i,4,m1) )**2
           s1 =  1.0 / ( epweno + tt1 )**2
           s2 =  6.0 / ( epweno + tt2 )**2
           s3 =  3.0 / ( epweno + tt3 )**2
            t0 = 1. / ( s1 + s2 + s3 )
            s1 = s1 * t0
            s3 = s3 * t0
         ff(i,m) = ff(i,m) + ( s1*(t2-t1) + (0.5*s3-0.25)*(t3-t2) ) /3.
250    continue

260    continue

*  end of the big loop

* Project the numerical flux to the physical space:

        do 262 m=1,mn
        do 262 i=-1+nyl(j),nyr(j)
       fh(i,m)=(-fp(i-1,m)+7*(fp(i,m)+fp(i+1,m))-fp(i+2,m))/12.
     &          +evr(i,m,1)*ff(i,1)+evr(i,m,2)*ff(i,2)
     &          +evr(i,m,3)*ff(i,3)+evr(i,m,4)*ff(i,4)
262     continue

       do 266 m=1,mn
       do 266 i=nyl(j),nyr(j)
        f(j,i,m)=f(j,i,m)+(fh(i-1,m)-fh(i,m))*cdy
266     continue
200   continue

      return
      end

