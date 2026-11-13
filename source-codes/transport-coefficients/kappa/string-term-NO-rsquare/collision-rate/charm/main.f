
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
      double precision kappa
      double precision qhat_cq,qhat_cg
      
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

ccc////      call readdata

ccccccccccccccccccccccccccccccccccccccccccccccccc
c calculate p1, T table for gamma and qhat

      open(20,file='gamma-table.dat',status='replace')
c      open(21,file='gamma-check.dat',status='replace')
      open(21,file='kappa-string.dat',status='replace')

c      do i=istart,nt+istart-1
      do i=1, nt 

         Te=dble(i+5)*delta_Te
c         Te=0.3d0
         mdprime=a_md+b_md*Te;
         msprime=dsqrt(as_ms+bs_ms*Te);


c         do j=1,np 

          p1=10.0d0

c            if(j.eq.0) then
c                p1=0.01d0
c            else
cc                p1=exp(0.1d0*j)
cc                p1=exp(dble(j)-2d0)
c                p1=10.0d0*j
c            endif

            e1=sqrt(p1**2+m1**2)

            call vegas8(region,ndim,fnc_cq,init,ncall,
     &              itmx,nprn,tgral,sd,chi2a)
            gam_cq=tgral/1e6

            call vegas8(region,ndim,fnc_cg,init,ncall,
     &              itmx,nprn,tgral,sd,chi2a)
            gam_cg=tgral/1e6

            gam_tot=gam_cq+gam_cg

            call vegas8(region,ndim,fnc_cq_qhat,init,ncall,
     &              itmx,nprn,tgral,sd,chi2a)
c            qhat_cq=tgral/1e6/Te**3
            qhat_cq=tgral/1e6

            call vegas8(region,ndim,fnc_cg_qhat,init,ncall,
     &              itmx,nprn,tgral,sd,chi2a)
c            qhat_cg=tgral/1e6/Te**3
            qhat_cg=tgral/1e6

c            qhat_over_T3=qhat_cq+qhat_cg
            kappa=(qhat_cq+qhat_cg)/2d0

c            write(20,101) i,j,gam_tot,gam_cq,gam_cg,qhat_over_T3
c            write(21,102) Te,p1,gam_tot,gam_cq,gam_cg,qhat_over_T3
            write(21,102) Te,p1,e1,kappa
c         enddo
      enddo

      close(20)
      close(21)

 101  format(I4,2X,I4,4(2X,E12.6))
 102  format(E12.6,3(2X,E12.6))

      end
         
