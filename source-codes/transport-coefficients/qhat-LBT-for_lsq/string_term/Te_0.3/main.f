
      program main

      implicit none

      include 'coms.in'

      double precision fnc_cg,fnc_cq,fnc_gamma,fnc_cq_qhat,fnc_cg_qhat
      external fnc_cg,fnc_cq,fnc_gamma,fnc_cq_qhat,fnc_cg_qhat
      double precision p1,argu(4)
      integer ndim,itmx,ncall,nprn,init
      double precision region(8),sd,chi2a,tgral
      integer i,j,istart
      double precision gam_cq,gam_cg,gam_tot,qhat_over_T3
      double precision delta_Te
      
      double precision KPfactor,runKT,preKT,Kfactor

      double precision aaa_cq,aaa_cg
      
c      character*20 argFromCommand
c      integer numOfArg
c
c      numOfArg = IARGC()
c      if(numOfArg.ne.1) then
c         write(6,*) "Please specify i for T to start with."
c         stop
c      else
c         call GETARG(1,argFromCommand)
c         read(argFromCommand,*) istart
c      endif

      istart=0 ! start temperature at 0.1 GeV

      delta_Te=0.02d0

      ndim=4
      itmx=6
      ncall=1000000
      nprn=0
      init=0

      region(1)=0d0
      region(2)=0d0
      region(3)=0d0
      region(4)=0d0

      region(5)=PI
      region(6)=PI
      region(7)=2d0*PI
      region(8)=100d0

ccccccccccccccccccccccccccccccccccccccccccccccccc

ccccccccccccccccccccccccccccccccccccccccccccccccc
c calculate p1, T table for gamma and qhat

      open(20,file='qhat_vs_p-T_0.3.dat',status='replace')

         Te=0.3d0

         mdprime=a_md+b_md*Te
         msprime=dsqrt(as_ms+bs_ms*Te)

c         do j=0,np
         do j=1,50

            p1=0.5d0+j*2d0

            e1=sqrt(p1**2+m1**2)

            call vegas8(region,ndim,fnc_cq_qhat,init,ncall,
     &              itmx,nprn,tgral,sd,chi2a)
            aaa_cq=tgral/1e6/Te**3

            call vegas8(region,ndim,fnc_cg_qhat,init,ncall,
     &              itmx,nprn,tgral,sd,chi2a)
            aaa_cg=tgral/1e6/Te**3

            qhat_over_T3=aaa_cq+aaa_cg

            Kfactor=1d0
            qhat_over_T3=qhat_over_T3*Kfactor

            write(20,*) p1,qhat_over_T3

         enddo

      close(20)

 101  format(I4,2X,I4,4(2X,E12.6))
 102  format(E12.6,5(2X,E12.6))

      end
         
