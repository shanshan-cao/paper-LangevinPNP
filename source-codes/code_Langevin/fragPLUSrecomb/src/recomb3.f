cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       subroutine force_recomb3(hadrflag,lightFlavor1,lightFlavor2,
     &            orbital,massH,id_hadFO,mass_hadFO)
c do recombination to heavy baryon, only Lambda_c is available now
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

       include 'common.f' 
      
       integer hadrflag,i,loopcount,lightFlavor1,lightFlavor2,orbital
       double precision rlu,phi,utheta,temp,plight
       double precision kmin,kmin2,plmax,plmaxcm,gacell,gacm
       double precision becell,becm,Elmax,Elmaxcm,pHQ,dum_random
       double precision reducedM1,reducedM2
       double precision massH,mass_hadFO
       integer id_hadFO
c generate one light parton according to thermal distribution
       temp=Thydro(ipin)
       if(temp.lt.0.1d0) return

       if(lightFlavor1.eq.1) then
          mlight1=m_ud
       elseif(lightFlavor1.eq.3) then
          mlight1=m_s
       else
          write(6,*) "Flavor not available!"
          stop
       endif
       if(lightFlavor2.eq.1) then
          mlight2=m_ud
       elseif(lightFlavor2.eq.3) then
          mlight2=m_s
       else
          write(6,*) "Flavor not available!"
          stop
       endif

       ! estimation for 3 body, not sure the exact solution yet
       mlight=mlight1+mlight2
       plmax=15d0*temp
       Elmax=sqrt(plmax**2+mlight**2)
       becell=sqrt(c_vx(ipin)**2+c_vy(ipin)**2+c_vz(ipin)**2)
       gacell=1d0/sqrt(1d0-becell**2)
       plmaxcm=gacell*(plmax+becell*Elmax)
       Elmaxcm=gacell*(Elmax+becell*plmax)
       pHQ=sqrt(p_p0(ipin)**2-p_mass(ipin)**2)
       becm=(pHQ+plmaxcm)/(p_p0(ipin)+Elmaxcm)
       gacm=1d0/sqrt(1d0-becm**2)
       kmin=gacm*(pHQ-becm*p_p0(ipin))
       kmin2=kmin**2
       if(pHQ.lt.plmaxcm) kmin2=-1d0

c calculate parameters for Wigner function
       reducedM1=mlight1*mlight2/(mlight1+mlight2)
       reducedM2=(mlight1+mlight2)*p_mass(ipin)/
     &           (mlight1+mlight2+p_mass(ipin))
       sigma11=1d0/sqrt(reducedM1*omega_HB)
       sigma12=1d0/sqrt(reducedM2*omega_HB)

       reducedM1=mlight1*p_mass(ipin)/(mlight1+p_mass(ipin))
       reducedM2=(mlight1+p_mass(ipin))*mlight2/
     &           (mlight1+mlight2+p_mass(ipin))
       sigma21=1d0/sqrt(reducedM1*omega_HB)
       sigma22=1d0/sqrt(reducedM2*omega_HB)
       reducedM1=mlight2*p_mass(ipin)/(mlight2+p_mass(ipin))
       reducedM2=(mlight2+p_mass(ipin))*mlight1/
     &           (mlight1+mlight2+p_mass(ipin))
       sigma31=1d0/sqrt(reducedM1*omega_HB)
       sigma32=1d0/sqrt(reducedM2*omega_HB)

       loopcount=0
       do while(hadrflag.eq.1)
          call judgehadmech3(orbital,hadrflag,kmin2,massH)
          loopcount=loopcount+1
          if(loopcount.gt.1D8) then
             write(6,*) "too many loops, no recomb"
             return
          endif
       enddo

c do recombination if hadrflag=2
       op_num(ipout)=ipout

       dum_random=rlu(0)

       if(abs(p_id(ipin)).eq.4) then ! c quark
c          op_id(ipout)=sign(4122,p_id(ipin))
c          op_mass(ipout)=2.286d0
cccccccccccccccccccc out-charm-baryons  lisqccccccccccccccc
          op_id(ipout)=sign(id_hadFO,p_id(ipin))
          op_mass(ipout)=mass_hadFO
cccccccccccccccccccccccccccccccccccccccccccccccccccccc          
       elseif(abs(p_id(ipin)).eq.5) then ! b quark
          write(6,*) "Not ready for B baryon"
          stop
       else
          write(6,*) "Heavy meson not listed"
          stop
       endif

c in this version, recombine in 2-body c.m. first and then boost out
       op_px(ipout)=0d0
       op_py(ipout)=0d0
       op_pz(ipout)=0d0
       op_p0(ipout)=directM
       !op_p0(ipout)=massH

       call rotbos(0d0,0d0,vx_HM,vy_HM,vz_HM,
     &      op_px(ipout),op_py(ipout),op_pz(ipout),op_p0(ipout))

       op_r0(ipout)=p_r0(ipin)
       op_rx(ipout)=p_rx(ipin)
       op_ry(ipout)=p_ry(ipin)
       op_rz(ipout)=p_rz(ipin) 
       oThydro(ipout)=Thydro(ipin)
       osdensity(ipout)=sdensity(ipin)
       oc_vx(ipout)=c_vx(ipin)
       oc_vy(ipout)=c_vy(ipin)
       oc_vz(ipout)=c_vz(ipin)
       op_wt(ipout)=p_wt(ipin)
       op_ipx(ipout)=p_ipx(ipin)
       op_ipy(ipout)=p_ipy(ipin)
       op_ipz(ipout)=p_ipz(ipin)
       call rotbos(0d0,0d0,c_vx(ipin),c_vy(ipin),c_vz(ipin),
     &            op_px(ipout),op_py(ipout),op_pz(ipout),op_p0(ipout))

       call decay2body(op_px(ipout),op_py(ipout),op_pz(ipout),
     &      op_p0(ipout),op_mass(ipout))

       ipout=ipout+1
       count_recomb=count_recomb+1
       judge_HM=1
       frag=.true.

       return
       end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       subroutine judgehadmech3(orbital,hadrflag,kminsquare,massH)
c judge the hadronization mechanism, frag. or recomb.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

       include 'common.f'

       integer hadrflag,orbital
       double precision cmbetax,cmbetay,cmbetaz
       double precision plight1,plight2
       double precision q11x_pair,q11y_pair,q11z_pair,q11pair2
       double precision q12x_pair,q12y_pair,q12z_pair,q12pair2
       double precision q21x_pair,q21y_pair,q21z_pair,q21pair2
       double precision q22x_pair,q22y_pair,q22z_pair,q22pair2
       double precision q31x_pair,q31y_pair,q31z_pair,q31pair2
       double precision q32x_pair,q32y_pair,q32z_pair,q32pair2
       double precision probWigner1,probWigner2,probWigner3
       double precision probWigner,temp
       double precision rlu,utheta,phi
       double precision probmax,kminsquare,distribution,Egluon
       double precision massH,p0_HQ,px_HM,py_HM,pz_HM,p0_HM

c generate one light parton according to thermal distribution
       temp=Thydro(ipin)
c       temp = temp_common
c       temp=0.165d0

       plight1=rlu(0)*15d0*temp
       p0_light1=sqrt(plight1**2+mlight1**2)
       Egluon=sqrt(plight1**2+m_g**2)
       distribution=6d0*plight1**2/(exp(p0_light1/temp)+1d0)+
     &              plight1**2/(exp(Egluon/temp)-1d0)*16d0/6d0
       do while(rlu(0)*mxvalue1.gt.distribution)
          plight1=rlu(0)*15d0*temp
          p0_light1=sqrt(plight1**2+mlight1**2)
          Egluon=sqrt(plight1**2+m_g**2)
          distribution=6d0*plight1**2/(exp(p0_light1/temp)+1d0)+
     &              plight1**2/(exp(Egluon/temp)-1d0)*16d0/6d0
       enddo
       utheta=rlu(0)*2d0-1d0
       phi=rlu(0)*2d0*pi
       px_light1=plight1*sqrt(1-utheta*utheta)*cos(phi)
       py_light1=plight1*sqrt(1-utheta*utheta)*sin(phi)
       pz_light1=plight1*utheta

       plight2=rlu(0)*15d0*temp
       p0_light2=sqrt(plight2**2+mlight2**2)
       Egluon=sqrt(plight2**2+m_g**2)
       distribution=6d0*plight2**2/(exp(p0_light2/temp)+1d0)+
     &              plight2**2/(exp(Egluon/temp)-1d0)*16d0/6d0
       do while(rlu(0)*mxvalue2.gt.distribution)
          plight2=rlu(0)*15d0*temp
          p0_light2=sqrt(plight2**2+mlight2**2)
          Egluon=sqrt(plight2**2+m_g**2)
          distribution=6d0*plight2**2/(exp(p0_light2/temp)+1d0)+
     &              plight2**2/(exp(Egluon/temp)-1d0)*16d0/6d0
       enddo
       utheta=rlu(0)*2d0-1d0
       phi=rlu(0)*2d0*pi
       px_light2=plight2*sqrt(1-utheta*utheta)*cos(phi)
       py_light2=plight2*sqrt(1-utheta*utheta)*sin(phi)
       pz_light2=plight2*utheta

c       call rotbos(0d0,0d0,c_vx(ipin),c_vy(ipin),c_vz(ipin),
c     &            px_light,py_light,pz_light,p0_light)

c if use pure recombination mechanism, i.e. fW = 1
       if (hadr_flag.eq.2) then
          hadrflag = 2
          return
       endif

c boost heavy/light quark to the c.m. frame of the heavy-light system

       cmbetax=(p_px(ipin)+px_light1+px_light2)/
     &         (p_p0(ipin)+p0_light1+p0_light2)
       cmbetay=(p_py(ipin)+py_light1+py_light2)/
     &         (p_p0(ipin)+p0_light1+p0_light2)
       cmbetaz=(p_pz(ipin)+pz_light1+pz_light2)/
     &         (p_p0(ipin)+p0_light1+p0_light2)

       call rotbos(0d0,0d0,-cmbetax,-cmbetay,-cmbetaz,
     &            px_light1,py_light1,pz_light1,p0_light1)
       call rotbos(0d0,0d0,-cmbetax,-cmbetay,-cmbetaz,
     &            px_light2,py_light2,pz_light2,p0_light2)
       call rotbos(0d0,0d0,-cmbetax,-cmbetay,-cmbetaz,
     &            p_px(ipin),p_py(ipin),p_pz(ipin),p_p0(ipin))

c calculate the relative momentum

       ! config 1: combine 2 light quarks first, then combine HQ
       q11x_pair=(p0_light2*px_light1-p0_light1*px_light2)/
     &          (p0_light1+p0_light2)
       q11y_pair=(p0_light2*py_light1-p0_light1*py_light2)/
     &          (p0_light1+p0_light2)
       q11z_pair=(p0_light2*pz_light1-p0_light1*pz_light2)/
     &          (p0_light1+p0_light2)
       q11pair2=q11x_pair**2+q11y_pair**2+q11z_pair**2

       q12x_pair=(p_p0(ipin)*(px_light1+px_light2)-
     &          (p0_light1+p0_light2)*p_px(ipin))/
     &          (p0_light1+p0_light2+p_p0(ipin))
       q12y_pair=(p_p0(ipin)*(py_light1+py_light2)-
     &          (p0_light1+p0_light2)*p_py(ipin))/
     &          (p0_light1+p0_light2+p_p0(ipin))
       q12z_pair=(p_p0(ipin)*(pz_light1+pz_light2)-
     &          (p0_light1+p0_light2)*p_pz(ipin))/
     &          (p0_light1+p0_light2+p_p0(ipin))
       q12pair2=q12x_pair**2+q12y_pair**2+q12z_pair**2

       ! config 2: combine HQ with light 1 first, then combine 2 
       q21x_pair=(p_p0(ipin)*px_light1-p0_light1*p_px(ipin))/
     &          (p0_light1+p_p0(ipin))
       q21y_pair=(p_p0(ipin)*py_light1-p0_light1*p_py(ipin))/
     &          (p0_light1+p_p0(ipin))
       q21z_pair=(p_p0(ipin)*pz_light1-p0_light1*p_pz(ipin))/
     &          (p0_light1+p_p0(ipin))
       q21pair2=q21x_pair**2+q21y_pair**2+q21z_pair**2

       q22x_pair=(p0_light2*(px_light1+p_px(ipin))-
     &          (p0_light1+p_p0(ipin))*px_light2)/
     &          (p0_light1+p0_light2+p_p0(ipin))
       q22y_pair=(p0_light2*(py_light1+p_py(ipin))-
     &          (p0_light1+p_p0(ipin))*py_light2)/
     &          (p0_light1+p0_light2+p_p0(ipin))
       q22z_pair=(p0_light2*(pz_light1+p_pz(ipin))-
     &          (p0_light1+p_p0(ipin))*pz_light2)/
     &          (p0_light1+p0_light2+p_p0(ipin))
       q22pair2=q22x_pair**2+q22y_pair**2+q22z_pair**2

       ! config 3: combine HQ with light 2 first, then combine 1 
       q31x_pair=(p_p0(ipin)*px_light2-p0_light2*p_px(ipin))/
     &          (p0_light2+p_p0(ipin))
       q31y_pair=(p_p0(ipin)*py_light2-p0_light2*p_py(ipin))/
     &          (p0_light2+p_p0(ipin))
       q31z_pair=(p_p0(ipin)*pz_light2-p0_light2*p_pz(ipin))/
     &          (p0_light2+p_p0(ipin))
       q31pair2=q31x_pair**2+q31y_pair**2+q31z_pair**2

       q32x_pair=(p0_light1*(px_light2+p_px(ipin))-
     &          (p0_light2+p_p0(ipin))*px_light1)/
     &          (p0_light1+p0_light2+p_p0(ipin))
       q32y_pair=(p0_light1*(py_light2+p_py(ipin))-
     &          (p0_light2+p_p0(ipin))*py_light1)/
     &          (p0_light1+p0_light2+p_p0(ipin))
       q32z_pair=(p0_light1*(pz_light2+p_pz(ipin))-
     &          (p0_light2+p_p0(ipin))*pz_light1)/
     &          (p0_light1+p0_light2+p_p0(ipin))
       q32pair2=q32x_pair**2+q32y_pair**2+q32z_pair**2

       ! calculate Wigner function for recombination probability

       if(orbital.eq.0) then ! s wave
          probmax=max_3body_s
          probWigner1=(sigma11*sigma12)**3
     &               *exp(-q11pair2*sigma11**2-q12pair2*sigma12**2)
          probWigner2=(sigma21*sigma22)**3
     &               *exp(-q21pair2*sigma21**2-q22pair2*sigma22**2)
          probWigner3=(sigma31*sigma32)**3
     &               *exp(-q31pair2*sigma31**2-q32pair2*sigma32**2)
       else if(orbital.eq.1) then ! p wave
          ! symmetrize between l1=1 & l2=0 and l1=0 & l2=1
          probmax=max_3body_p
          probWigner1=(sigma11*sigma12)**3
     &               *exp(-q11pair2*sigma11**2-q12pair2*sigma12**2)
     &               *1d0/3d0*(q11pair2*sigma11**2+q12pair2*sigma12**2)
          probWigner2=(sigma21*sigma22)**3
     &               *exp(-q21pair2*sigma21**2-q22pair2*sigma22**2)
     &               *1d0/3d0*(q21pair2*sigma21**2+q22pair2*sigma22**2)
          probWigner3=(sigma31*sigma32)**3
     &               *exp(-q31pair2*sigma31**2-q32pair2*sigma32**2)
     &               *1d0/3d0*(q31pair2*sigma31**2+q32pair2*sigma32**2)
       else
          write(6,*) "only s and p waves are available"
          stop
       endif      

       ! symmetrize between 3 config's and add common the factor
       probWigner=(2d0*sqrt(PI))**6
     &           *(probWigner1+probWigner2+probWigner3)/3d0

       p0_HQ=p_p0(ipin)

c use Monte Carlo to decide whether to recombine
       if(rlu(0)*probmax.lt.probWigner) then
          hadrflag=2

          directM=p0_HQ+p0_light1+p0_light2
          vx_HM=cmbetax
          vy_HM=cmbetay
          vz_HM=cmbetaz

c          ! add/subtract a photon to conserve energy momentum
c          px_HM=0d0
c          py_HM=0d0
c          pz_HM=0d0
c          p0_HM=p0_HQ+p0_light1+p0_light2
c          call decay2body(px_HM,py_HM,pz_HM,p0_HM,massH)
c          call rotbos(0d0,0d0,cmbetax,cmbetay,cmbetaz,
c     &               px_HM,py_HM,pz_HM,p0_HM)
c          vx_HM=px_HM/p0_HM
c          vy_HM=py_HM/p0_HM
c          vz_HM=pz_HM/p0_HM
       endif 

       ! boost heavy/light quark back to the hydro frame (lrf)
       call rotbos(0d0,0d0,cmbetax,cmbetay,cmbetaz,
     &            px_light1,py_light1,pz_light1,p0_light1)
       call rotbos(0d0,0d0,cmbetax,cmbetay,cmbetaz,
     &            px_light2,py_light2,pz_light2,p0_light2)
       call rotbos(0d0,0d0,cmbetax,cmbetay,cmbetaz,
     &            p_px(ipin),p_py(ipin),p_pz(ipin),p_p0(ipin))

       return
       end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


