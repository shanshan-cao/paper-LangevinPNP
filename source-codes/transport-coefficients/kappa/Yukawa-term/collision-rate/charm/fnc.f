cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   This file contains the functions and subroutines for the calculation
c   of heavy quark scattering rate and qhat in a QGP medium
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine sigma_factor(argu,s,t,u,e2,e3,e4,sigFactor,
     &                        theta_cq_flag,theta_cg_flag)

      implicit none
      include 'coms.in'

      double precision argu(4),s,t,u,sigFactor
      integer theta_cq_flag,theta_cg_flag
      double precision t2,t4,p24,e2,c2,c4,s2s4,c24,p1,p2,e4,down,p4,
     &                 mud2,e3

      double precision px1,py1,pz1
      double precision px2,py2,pz2
      double precision px3,py3,pz3,p3
      double precision px4,py4,pz4
      double precision qperp2

      t2=argu(1)
      t4=argu(2)
      p24=argu(3)
      e2=argu(4)

      c2=cos(t2)
      c4=cos(t4)
      s2s4=sin(t2)*sin(t4)
      c24=c2*c4+s2s4*cos(p24)
     
      p1=sqrt(e1**2-m1**2)
      p2=sqrt(e2**2-m2**2)
      e4=(m1**2-m3**2)/2d0+e2*(e1-p1*c2)
      down=e1-p1*c4+e2*(1d0-c24)
      e4=e4/down
      p4=sqrt(e4**2-m4**2)

      e3=e1+e2-e4

      sigFactor=s2s4*e2*e4/down*1e6      

      s=m1**2+m2**2+2*e1*e2-2*p1*p2*c2
      t=m2**2+m4**2-2*e2*e4+2*p2*p4*c24
      u=m1**2+m4**2-2*e1*e4+2*p1*p4*c4

      mud2=6d0*PI*alpha*Te**2

      px1=0d0
      py1=0d0
      pz1=p1

      px2=p2*sin(t2)
      py2=0d0
      pz2=p2*cos(t2)

      px4=p4*sin(t4)*cos(p24)
      py4=p4*sin(t4)*sin(p24)
      pz4=p4*cos(t4)

      px3=px1+px2-px4
      py3=py1+py2-py4
      pz3=pz1+pz2-pz4
      p3=sqrt(px3**2+py3**2+pz3**2)

      qperp2=p3**2-((px1*px3+py1*py3+pz1*pz3)/p1)**2

c      if(s.ge.2d0*mud2.and.t.le.(-mud2).and.u.le.(-mud2))then
      if(s.gt.(m1+mdprime)**2d0
     #                .and.qperp2.le.(10d0*Te)**2) then
         theta_cq_flag=1
      else
         theta_cq_flag=0
      endif

      if(s.gt.(m1+mdprime)**2d0
     #                .and.qperp2.le.(10d0*Te)**2) then
         theta_cg_flag=1
      else
         theta_cg_flag=0
      endif

      end subroutine

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      function MM_cg2cg(s,t,u,e2,e3,e4)

      implicit none
      include 'coms.in'

      double precision s,t,u,MM,MM_cg2cg,m2m
      double precision e2,e3,e4,f2,f3,f4
      double precision md     

      m2m=m1*m1

      md=mdprime

      MM=64d0*((s-m2m)*(m2m-u)+2d0*m2m*(s+m2m))/(9d0*(s-m2m)**2)
      MM=MM+8d0*(5d0*m2m**2+3d0*m2m*t-10d0*m2m*u+4d0*t**2+5d0*t*u+5d0*u**2+(m2m-s)*(m2m-u))/(t-md**2)**2
      MM=MM+64d0*((s-m2m)*(m2m-u)+2d0*m2m*(u+m2m))/(9d0*(u-m2m)**2)
      MM=MM+16d0*(3d0*m2m**2-3d0*m2m*s-m2m*u+s**2)/((s-m2m)*(t-md**2))
      MM=MM+16d0/9d0*(m2m*(4d0*m2m-t))/((s-m2m)*(m2m-u))
      MM=MM+16d0*(3d0*m2m**2-m2m*s-3d0*m2m*u+u**2)/((t-md**2)*(u-m2m))


      f2=dofG/(exp(e2/Te)-1d0)
      f3=1d0/(exp(e3/Te)+1d0)
      f4=1d0/(exp(e4/Te)-1d0)

      MM_cg2cg=MM*f2*(1d0+f4)

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      function MM_cq2cq(s,t,u,e2,e3,e4)  ! DOF from flavor: 3*2=6

      implicit none
      include 'coms.in'

      double precision s,t,u,MM,MM_cq2cq,m2m
      double precision e2,e3,e4,f2,f3,f4
      double precision md

      m2m=m1*m1

      md=mdprime

      MM=((m2m-u)**2+(s-m2m)**2+2d0*m2m*t)/(t-md**2)**2

      f2=dofQ/(exp(e2/Te)+1d0)
      f3=1d0/(exp(e3/Te)+1d0)
      f4=1d0/(exp(e4/Te)+1d0)

      MM_cq2cq=MM*64d0/9d0*6d0*f2*(1d0-f4)

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      function fnc_cg(argu)

      implicit none

      include 'coms.in'

      double precision fnc_cg,argu(4)
      double precision MM,s,t,u,e2,e3,e4,sigFactor,MM_cg2cg
      integer theta_cq_flag,theta_cg_flag
      double precision aa
     

      call sigma_factor(argu,s,t,u,e2,e3,e4,sigFactor,theta_cq_flag,theta_cg_flag)

      if(theta_cg_flag.eq.0) then
         fnc_cg=0d0
         return
      endif

      aa=alpha**2/(256d0*PI*PI*e1)

      MM=MM_cg2cg(s,t,u,e2,e3,e4)

      fnc_cg=aa*MM*sigFactor

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      function fnc_cq(argu)

      implicit none

      include 'coms.in'

      double precision fnc_cq,argu(4)
      double precision MM,s,t,u,e2,e3,e4,sigFactor,MM_cq2cq
      integer theta_cq_flag,theta_cg_flag
      double precision aa
     

      call sigma_factor(argu,s,t,u,e2,e3,e4,sigFactor,theta_cq_flag,theta_cg_flag)

      if(theta_cq_flag.eq.0) then
         fnc_cq=0d0
         return
      endif

      aa=alpha**2/(256d0*PI*PI*e1)

      MM=MM_cq2cq(s,t,u,e2,e3,e4)

      fnc_cq=aa*MM*sigFactor

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      function fnc_cq_qhat(argu)

      implicit none

      include 'coms.in'

      double precision fnc_cq_qhat,argu(4)
      double precision MM,s,t,u,e2,e3,e4,sigFactor
      double precision MM_cg2cg,MM_cq2cq
      integer theta_cq_flag,theta_cg_flag
      double precision aa,f2g,f2q,f4g,f4q
      double precision px1,py1,pz1,px2,py2,pz2,px3,py3,pz3,px4,py4,pz4
      double precision p1,p2,p3,p4,theta2,theta4,phi4,qperp2

      call sigma_factor(argu,s,t,u,e2,e3,e4,sigFactor,theta_cq_flag,theta_cg_flag)

      if(theta_cq_flag.eq.0) then
         fnc_cq_qhat=0d0
         return
      endif

      theta2=argu(1)
      theta4=argu(2)
      phi4=argu(3)

      p1=sqrt(e1**2-m1**2)
      px1=0d0
      py1=0d0
      pz1=p1

      p2=sqrt(e2**2-m2**2)
      px2=p2*sin(theta2)
      py2=0d0
      pz2=p2*cos(theta2)

      p4=sqrt(e4**2-m4**2)
      px4=p4*sin(theta4)*cos(phi4)
      py4=p4*sin(theta4)*sin(phi4)
      pz4=p4*cos(theta4)

      px3=px1+px2-px4
      py3=py1+py2-py4
      pz3=pz1+pz2-pz4
      p3=sqrt(px3**2+py3**2+pz3**2)

      qperp2=p3**2-((px1*px3+py1*py3+pz1*pz3)/p1)**2
      
      aa=alpha**2/(256d0*PI*PI*e1)

      MM=MM_cq2cq(s,t,u,e2,e3,e4)

      fnc_cq_qhat=aa*sigFactor*qperp2*MM

      return
      end


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      function fnc_cg_qhat(argu)

      implicit none

      include 'coms.in'

      double precision fnc_cg_qhat,argu(4)
      double precision MM,s,t,u,e2,e3,e4,sigFactor
      double precision MM_cg2cg,MM_cq2cq
      integer theta_cq_flag,theta_cg_flag
      double precision aa,f2g,f2q,f4g,f4q
      double precision px1,py1,pz1,px2,py2,pz2,px3,py3,pz3,px4,py4,pz4
      double precision p1,p2,p3,p4,theta2,theta4,phi4,qperp2

      call sigma_factor(argu,s,t,u,e2,e3,e4,sigFactor,theta_cq_flag,theta_cg_flag)

      if(theta_cg_flag.eq.0) then
         fnc_cg_qhat=0d0
         return
      endif

      theta2=argu(1)
      theta4=argu(2)
      phi4=argu(3)

      p1=sqrt(e1**2-m1**2)
      px1=0d0
      py1=0d0
      pz1=p1

      p2=sqrt(e2**2-m2**2)
      px2=p2*sin(theta2)
      py2=0d0
      pz2=p2*cos(theta2)

      p4=sqrt(e4**2-m4**2)
      px4=p4*sin(theta4)*cos(phi4)
      py4=p4*sin(theta4)*sin(phi4)
      pz4=p4*cos(theta4)

      px3=px1+px2-px4
      py3=py1+py2-py4
      pz3=pz1+pz2-pz4
      p3=sqrt(px3**2+py3**2+pz3**2)

      qperp2=p3**2-((px1*px3+py1*py3+pz1*pz3)/p1)**2
      
      aa=alpha**2/(256d0*PI*PI*e1)

      MM=MM_cg2cg(s,t,u,e2,e3,e4)

      fnc_cg_qhat=aa*sigFactor*qperp2*MM

      return
      end


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
