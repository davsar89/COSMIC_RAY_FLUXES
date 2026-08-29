! *******************************************************
function getSpec(ip,s,r,d,e,g)
! ip: particle ID
! s: W index
! r: cut-off rigidity in GV
! d: atmospheric depth in g/cm2
! e: energy in MeV/n
! g: local geometry effect
! *******************************************************
   implicit real*8 (a-h, o-z)

   getSpec=0.0
   if(ip.eq.0) then  ! neutron
      getSpec=getNeutSpec(s,r,d,g,e)
   elseif(ip.ge.1.and.ip.le.28) then ! proton to Ni
      getSpec=getIonSpec(ip,s,r,d,e)
   elseif(ip.eq.29.or.ip.eq.30) then ! Muon
      iptmp=ip-28
      getSpec=getMuonSpec(iptmp,s,r,d,e)
   elseif(ip.ge.31.and.ip.le.33) then
      iptmp=ip-28 ! 3:electron, 4:positron, 5:photon
      getSpec=getsecondary(iptmp,s,r,d,e)
   else
      ! Invalid particle identifiers are rejected by the caller.  Keep the
      ! library routine memory-safe for external callers by returning zero.
      getSpec=0.0d0
   endif

   return

end

! *******************************************************
function getFl(ip,s,r,d)  ! get Fl value, s:solar modulation potential, r:Cut off rigidity, d:depth
! *******************************************************
   implicit real*8 (a-h, o-z)
   parameter (npart=11) ! number of particle type, 0:neutron, 1:proton, 2:alpha, 3:electron, 4:positron, 5:photon, 6-11: Li-O
   parameter (nBdata=4) ! B(1) - B(4) : Fl= B(1)*(exp(-B(2)*d)-B(3)*exp(-B(4)*d))
   parameter (nAdata=10) ! A(1) - A(10) : Bmin = A(1)+A(2)*r+A(3)/(1+exp((r-A(4))/A(5))), Bmin = A(6)+A(7)*r+A(8)/(1+exp((r-A(9))/A(10)))
   parameter (nsor=2) ! solar minimum & maximum
   character chatmp1*1,chatmp5*5
   character, save:: pname(0:npart)*6
   real*8, save:: A(0:npart,nBdata,nAdata),B(nBdata)
   real*8, save:: spot(nsor),FL(nsor)
   integer*4, save:: ifirst
   data ifirst/0/
   data spot/0.0,150.0/
   data pname/'neutro','proton','alphaa','elemag','elemag','elemag', &
   & 'ions  ','ions  ','ions  ','ions  ','ions  ','ions  '/

   if(ifirst.eq.0) then
      ifirst=1
      do i=0,5
         if(i.le.2) then ! neutron, proton, alpha
            open(unit=15,file='input/'//pname(i)//'/Rigid-Dep.inp',status='old')
         elseif(i.eq.3) then
            open(unit=15,file='input/'//pname(i)//'/Rigid-Dep-EL.inp',status='old')
         elseif(i.eq.4) then
            open(unit=15,file='input/'//pname(i)//'/Rigid-Dep-PO.inp',status='old')
         elseif(i.eq.5) then
            open(unit=15,file='input/'//pname(i)//'/Rigid-Dep-PH.inp',status='old')
         endif
         read(15,'(a)') chatmp1
         do ib=1,nBdata
            read(15,1010) chatmp5,(A(i,ib,ia),ia=1,nAdata)
         enddo
         close(unit=15)
      enddo
      open(unit=15,file='input/ions/Rigid-Dep.inp',status='old')
      do i=npart-5,npart
         read(15,'(a)') chatmp1
         do ib=1,nBdata
            read(15,1010) chatmp5,(A(i,ib,ia),ia=1,nAdata)
         enddo
      enddo
      close(unit=15)
   endif
1010 format(a5,30es13.5)

   do is=1,nsor ! solar minimum and maximum
      if(ip.eq.1.or.ip.eq.2.or.ip.ge.npart-5) then ! need not correction
         r1=r
      else ! for neutron, electron, positron and photon, Rc for high-altitude should be corrected
         if(ip.eq.0) then
            ipidx=ip ! for neutron
         else
            ipidx=ip-2 ! 1:electron, 2:positron, 3:photon
         endif
         r1=r*getBestR(is,r,d,ipidx)
      endif
      do ib=1,nBdata
         if(is.eq.1) then ! solar minimum
            B(ib)=A(ip,ib,1)+A(ip,ib,2)*r1+A(ip,ib,3)/(1+exp((r1-A(ip,ib,4))/A(ip,ib,5)))
         else
            B(ib)=A(ip,ib,6)+A(ip,ib,7)*r1+A(ip,ib,8)/(1+exp((r1-A(ip,ib,9))/A(ip,ib,10)))
         endif
      enddo
      Fl(is)=B(1)*(exp(-B(2)*d)-B(3)*exp(-B(4)*d))
   enddo

   if(ip.le.5) then
      pow=getPow(ip,d,r) ! get Power index
   else
      pow=getPow(ip+2,d,r) ! in getPow, ip should be +2 for ions
   endif

   A2=(Fl(1)-Fl(2))/(getFFPfromW(spot(1))**pow-getFFPfromW(spot(2))**pow)
   A1=Fl(1)-A2*getFFPfromW(spot(1))**pow
   getFl=a1+a2*getFFPfromW(s)**pow

   return
end

! **********************************************************
function getFFPfromW(s) ! get FFP (MV) from W (sun spot number)
! **********************************************************
   implicit real*8 (a-h, o-z)
   if(s.ge.0) then
      getFFPfromW=370.0+3.0e-1*s**1.45 ! FFP in MV
   else
      getFFPfromW=370.0-3.0e-1*abs(s)**1.45 ! FFP in MV
   endif
   return
end

! **********************************************************
function getRfromE(iz,ia,Ek,Em) ! get Rigidity in MV from Kinetic Energy (MeV/n)
! **********************************************************
   implicit real*8 (a-h, o-z)
   getRfromE=sqrt((ia*Ek)**2+2*Ek*ia*Em)/iz
   return
end

! **********************************************************
function getEfromR(iz,Rm,COR) ! get Kinetic Energy (MeV) from Rigidity (MV)
! **********************************************************
   implicit real*8 (a-h, o-z)
   getEfromR=sqrt((iz*COR)**2+Rm**2)-Rm
   return
end

! *******************************************************
function getPow(ip,d,r)  ! get Power of solar modulation dependence, r:Cut off rigidity, d:depth
! *******************************************************
   implicit real*8 (a-h, o-z)
   parameter (npart=13) ! number of particle type, neutron, proton, alpha, photon, electron positron, mu+, mu-, Li-O
   parameter (nBdata=2) ! B(1) - B(4) : Pow = b1 + b2*d
   parameter (nAdata=5) ! A(1) - A(5) : B = A(1)+A(2)*r+A(3)/(1+exp((r-A(4))/A(5)))
   character chatmp1*1,chatmp5*5
   character, save:: pname(0:npart)*6
   real*8, save:: A(0:npart,nBdata,nAdata),B(nBdata)
   integer*4, save:: ifirst
   data ifirst/0/
   data pname/'neutro','proton','alphaa','elemag','elemag','elemag','muon--','muon--', &
   & 'ions  ','ions  ','ions  ','ions  ','ions  ','ions  '/

   if(ifirst.eq.0) then
      ifirst=1
      do i=0,7
         if(i.le.2) then ! neutron, proton, alpha
            open(unit=15,file='input/'//pname(i)//'/solar-dep.inp',status='old')
         elseif(i.eq.3) then
            open(unit=15,file='input/'//pname(i)//'/solar-dep-EL.inp',status='old')
         elseif(i.eq.4) then
            open(unit=15,file='input/'//pname(i)//'/solar-dep-PO.inp',status='old')
         elseif(i.eq.5) then
            open(unit=15,file='input/'//pname(i)//'/solar-dep-PH.inp',status='old')
         elseif(i.eq.6) then
            open(unit=15,file='input/'//pname(i)//'/solar-dep.plus',status='old')
         elseif(i.eq.7) then
            open(unit=15,file='input/'//pname(i)//'/solar-dep.mins',status='old')
         endif
         read(15,'(a)') chatmp1
         do ib=1,nBdata
            read(15,*) (A(i,ib,ia),ia=1,nAdata)
         enddo
         close(unit=15)
      enddo
      open(unit=15,file='input/ions/solar-dep.inp',status='old')
      read(15,'(a)') chatmp1
      do i=npart-5,npart
         do ib=1,nBdata
            read(15,*) (A(i,ib,ia),ia=1,nAdata)
         enddo
      enddo
      close(15)
   endif

   do ib=1,nBdata
      b(ib)=a(ip,ib,1)+a(ip,ib,2)*r+a(ip,ib,3)/(1.0+exp((r-a(ip,ib,4))/a(ip,ib,5)))
   enddo

   getpow=b(1)+b(2)*d
   return
end



! **********************************************************
function getBestR(is,r,d,ip) ! get best Rc data for high-altitude correction
! **********************************************************
   parameter(nBdata=6)
   parameter(ndep=26)
   parameter(npart=3) ! high-altitude correction is necessary only for electron, positron, photon
   implicit real*8 (a-h, o-z)
   real*8, save:: A(0:npart,nBdata,ndep),B(nBdata)
   real*8, save:: dep(ndep)
   character chatmp1*1,chatmp5*5
   character pname(npart)*2

   data ifirst/0/
   data pname/'EL','PO','PH'/

   if(ifirst.eq.0) then
      a(:,:,:)=0.0
      do i=0,npart
         if(i.eq.0) then ! neutron
            open(unit=15,file='input/neutro/bestR.inp',status='old')
         else ! electron, positron, photon
            open(unit=15,file='input/elemag/bestR-'//pname(i)//'.inp',status='old')
         endif
         read(15,'(A)') chatmp1
         do id=1,ndep
            read(15,*) dep(id),(A(i,ib,id),ib=1,nBdata)
         enddo
         close(unit=15)
      enddo
      ifirst=1
   endif

   do id=1,ndep
      if(d.lt.dep(id)) exit
   enddo
   if(id.eq.1) then
      ratio=0.0
      id=2
   elseif(id.eq.ndep+1) then
      ratio=1.0
      id=ndep
   else
      ratio=(d-dep(id-1))/(dep(id)-dep(id-1))
   endif

   do ib=1,nBdata
      B(ib)=A(ip,ib,id-1)+(A(ip,ib,id)-A(ip,ib,id-1))*ratio
   enddo

   if(is.eq.1) then ! solar minimum
      getBestR=10**(b(1)+b(2)*r+b(3)/r)
   else
      getBestR=10**(b(4)+b(5)*r+b(6)/r)
   endif

   return

end

! *******************************************************
function getNeutspec(s,r1,d1,g,e) ! get neutron flux
!     s:Wolf number
!     r:cut off rigidity (GV)
!     d:air depth (g/cm^2)
!     g:local geometry parameter, 0=< g =< 1: water weight fraction, 10:no-earth, 100:blackhole, -10< g < 0: pilot, g < -10: cabin
!     e:neutron energy (MeV)
! *******************************************************
   implicit real*8 (a-h, o-z)
   parameter(nA=12) ! number of basic spectrum parameter
   parameter(nG=6) ! number of geometry parameter
   real*8, save:: A(nA) ! basic spectrum parameter
   real*8, save:: geo(nG) ! geometry parameter

   data ifirst/0/
   data airbus/2.45/  ! weight of airbus340 (100t)
   data Eth/2.5e-8/  ! themal energy

   if(ifirst.eq.0) then ! first time, get universal parameter (i.e: independent of all parameters)
!      Read A parameter
      open(unit=3,file='input/neutro/fitting-lowspec.inp',status='old')
      read(3,'(a)') chatmp1
      read(3,1001) (A(ia),ia=1,nA) !A(4)&A(12) is s,r,d-dependence, so will be changed
      close(unit=3)
      ifirst=1
   endif
1001 format(30e13.5)

   r=max(1.0,r1) ! secondary particle fluxes are the same for Rc<1GV
   d=max(0.15,d1) ! secondary particle fluxes are the same for d < 0.15 g/cm2

!     get condition dependent parameters
   Fl=getFl(0,s,r,d)
   A(12)=getA12(r,d)
   A(4)=getA4(r,d)
   call getGpara(g,geo) ! obtain G parameters
!     calculate flux (/cm^2/s/lethargy)
   x=e  ! x is energy

   evap=a(1)*(x/a(2))**a(3)*exp(-x/a(2))
   gaus=a(4)*exp(-(log10(x)-log10(a(5)))**2/(2*log10(a(6))**2))
   conti=a(7)*log10(x/a(8))*(1+tanh(a(9)*log10(x/a(10))))*(1-tanh(a(11)*log10(x/a(12))))

   basic=conti+evap+gaus

   basic=basic*CorrNeut(s,r,d,e) ! correction for high altitude

   fG=geo(1)+geo(2)*log10(x/geo(3))*(1-tanh(geo(4)*log10(x/geo(5))))
   if(g.lt.0.0) then
      if(g.lt.-10.0) then ! cabin
         gtmp=g+10.0
      else ! pilot
         gtmp=g
      endif
      fG=fG*(abs(gtmp)-10*int(abs(gtmp)/10.0))/airbus  ! consider size of aircraft
   endif
   geofactor=10.0**fG
   ther=geo(6)*(x/Eth)**2*exp(-(x/Eth))

   getNeutspec=Fl*(basic*geofactor+ther)/e

   return
end

! *******************************************************
subroutine getGpara(g,geo) ! get surroudning environment parameters
! *******************************************************
   implicit real*8 (a-h, o-z)
   parameter(nG=6) ! number of geometry parameter
   dimension geo(nG) ! geometry parameter
   real*8, save:: P(24) ! Input parameters read from input file, P(1)-P(3) from Geo-Dep, P(4)-P(14) from Water-Dep, P(15)-P(24) from aircraft-dep
   logical, save:: initialized=.false.
   character chatmp1*1,chatmp4*4

1001 format(a4,1x,30e13.5)

   if(.not.initialized) then
      open(unit=3,file='input/neutro/Geo-Dep.inp',status='old')
      read(3,'(a)') chatmp1
      read(3,*) p(1),p(2),p(3)
      close(Unit=3)
      open(unit=3,file='input/neutro/Water-Dep.inp',status='old')
      read(3,'(a)') chatmp1
      read(3,1001) chatmp4,p(4),p(5),p(6)
      read(3,1001) chatmp4,p(7),p(8),p(9)
      read(3,1001) chatmp4,p(10),p(11),p(12),p(13),p(14)
      close(unit=3)
      open(unit=3,file='input/neutro/Aircraft-Dep.inp',status='old')
      read(3,'(a)') chatmp1
      read(3,*) (p(ip),ip=15,19)  ! for pilot location
      read(3,*) (p(ip),ip=20,24)  ! for passenger & small aircraft configuration
      close(Unit=3)
      initialized=.true.
   endif

   if(g.ge.10.0) then ! in semi-infite atmosphere
      do ig=1,nG
         geo(ig)=0.0
      enddo
      geo(3)=1.0  ! if geo(3)=0, the value should be in NaN
      geo(5)=1.0  ! if geo(5)=0, the value should be in NaN
   elseif(g.ge.0.0) then ! for normal ground case
      geo(1)=p(1)
      geo(2)=p(2)
      geo(3)=10.0**(p(4)+p(5)/(p(6)+g))
      geo(4)=p(3)
      geo(5)=p(7)+p(8)*g+p(9)*g**2
      geo(6)=(p(10)+p(11)*exp(-p(12)*g))/(1+p(13)*exp(-p(14)*g))
   else   ! pilot or cabin location
      is=14
      if(g.lt.-10.0) is=is+5 ! for passenger & small aircraft configuration, skip 5 more data
      do i=1,5
         geo(i)=p(is+i)
      enddo
      geo(6)=0.0  ! no thermal component
   endif
   return
end

! *******************************************************
function getA4(r,d)  ! get A4 value, s:solar modulation potential, r:Cut off rigidity, d:depth
! *******************************************************
   implicit real*8 (a-h, o-z)
   parameter (nBdata=4) ! B(1) - B(4) : Fl= B(1)*(exp(-B(2)*d)-B(3)*exp(-B(4)*d))
   parameter (nAdata=6) ! A(1) - A(6) : B = A(1)+A(2)/(1+exp((r-A(3))/A(4))), A(5):A(1) for APmax, A(6):A(3) for APmax
   character chatmp1*1,chatmp5*5
   real*8, save:: A(nAdata),B(nBdata)
   logical, save:: initialized=.false.

   if(.not.initialized) then ! first time
      open(unit=15,file='input/neutro/Depth-Dep-mid.out',status='old')  ! read B(2)-B(4) (independent of s,r,d)
      read(15,'(a)') chatmp1
      read(15,'(a)') chatmp1
      read(15,*) tmp1,tmp2,B(2),B(3),B(4)
      close(unit=15)
      open(unit=15,file='input/neutro/Rigid-Dep.inp',status='old')
      do i=1,5  ! skip 5 line, 1 header line + 4 Bdata
         read(15,'(a)') chatmp1
      enddo
      read(15,1010) chatmp5,(A(ia),ia=1,nAdata)
      close(unit=15)
      initialized=.true.
   endif
1010 format(a5,30e13.5)

   B(1)=A(1)+A(2)*r+A(3)/(1+exp((r-A(4))/A(5)))
   getA4=B(1)+B(2)*d/(1+B(3)*exp(B(4)*d))

   return
end

! *******************************************************
function getA12(r,d)  ! get Fl value, s:solar modulation potential, r:Cut off rigidity, d:depth
! *******************************************************
   implicit real*8 (a-h, o-z)
   parameter (nBdata=4) ! B(1) - B(4) : Fl= B(1)*(exp(-B(2)*d)-B(3)*exp(-B(4)*d))
   parameter (nAdata=6) ! A(1) - A(6) : B = A(1)+A(2)/(1+exp((r-A(3))/A(4))), A(5):A(1) for APmax, A(6):A(3) for APmax
   character chatmp1*1,chatmp5*5
   real*8, save:: A(nBdata,nAdata),B(nBdata)
   logical, save:: initialized=.false.

   if(.not.initialized) then ! first time
      open(unit=15,file='input/neutro/Depth-Dep-hig.out',status='old')  ! read B(2),B(4) (independent of s,r,d)
      read(15,'(a)') chatmp1
      read(15,'(a)') chatmp1
      read(15,*) tmp0,tmp1,tmp2,tmp3,B(4)
      close(unit=15)
      open(unit=15,file='input/neutro/Rigid-Dep.inp',status='old')
      do i=1,6  ! skip 5 line, 1 header line + 5 Bdata
         read(15,'(a)') chatmp1
      enddo
      do ib=1,3 ! for B1 to B3
         read(15,1010) chatmp5,(A(ib,ia),ia=1,nAdata)
      enddo
      close(unit=15)
      initialized=.true.
   endif
1010 format(a5,30es13.5)

   do ib=1,3
      B(ib)=A(ib,1)+A(ib,2)*r+A(ib,3)/(1+exp((r-A(ib,4))/A(ib,5)))
   enddo

   getA12=B(1)*(exp(-B(2)*d)+B(3)*exp(-B(4)*d))

   return
end

! *******************************************************
function CorrNeut(s,r,d,e)  ! get correction factor for high altitude data, only solar minimum and maximum
! *******************************************************
   parameter(nhensu=9) ! number of hensu
   parameter(mpara=5) ! number of parameters to COR dependence
   parameter(ndep=26) ! number of depth
   parameter(nsol=2) ! solar minimum and maximum
   implicit real*8 (a-h,o-z)
   dimension ainp(mpara,nhensu,ndep,nsol)
   dimension dep(ndep)
   dimension b(nhensu),c(nsol) ! temporary dimension
   dimension spot(nsol)
   character chatmp1*1

   save ainp,spot,ifirst,dep

   data ifirst/0/
   data spot/0.0,150.0/

   if(ifirst.eq.0) then ! first time call
!     Read parameters for depth-independent parameters
! open(15,file='correction-depth-rigid-final.inp',status='old') ! only high altidue correction mode
      open(15,file='input/neutro/correction-depth-rigid.inp',status='old') ! all altitude correction mode
      read(15,'(a)') chatmp1
      do ih=1,nhensu
         do id=1,ndep
            read(15,*) itmp,dep(id),((ainp(ip,ih,id,is),ip=1,mpara),is=1,nsol)
         enddo
      enddo
      close(15)
      ifirst=1
   endif

! find depth ID
   do id=1,ndep
      if(d.lt.dep(id)) exit
   enddo
   if(id.eq.1) then
      ratio=0.0
      id=2
   elseif(id.eq.ndep+1) then
      ratio=1.0
      id=ndep
   else
      ratio=(d-dep(id-1))/(dep(id)-dep(id-1))
   endif

   rc=max(1.0,r)
   do is=1,nsol
      do ih=1,nhensu ! determine 9 hensu used in the correction equation
         d1=ainp(1,ih,id-1,is)+Ainp(2,ih,id-1,is)*rc+Ainp(3,ih,id-1,is)/(1+exp((rc-Ainp(4,ih,id-1,is))/Ainp(5,ih,id-1,is)))
         d2=ainp(1,ih,id-0,is)+Ainp(2,ih,id-0,is)*rc+Ainp(3,ih,id-0,is)/(1+exp((rc-Ainp(4,ih,id-0,is))/Ainp(5,ih,id-0,is)))
         b(ih)=d1+(d2-d1)*ratio
      enddo
      C(is)=10**(b(1)+(b(2)*log10(e)+b(3))*(1-tanh(b(4)*log10(e/b(5))))+(b(6)*log10(e)+b(7))*(1+tanh(b(8)*log10(e/b(9)))))
   enddo

   ip=0 ! always neutron
   pow=getPow(ip,d,r) ! get Power index
   A2=(c(1)-c(2))/(getFFPfromW(spot(1))**pow-getFFPfromW(spot(2))**pow)
   A1=c(1)-A2*getFFPfromW(spot(1))**pow
   CorrNeut=a1+a2*getFFPfromW(s)**pow

   return
end


! ******************************************************
function getMuonSpec(ip,s,r1,d1,e) ! get muon spectrum, ip=30:mu+, =31:mu-
! ******************************************************
   implicit real*8 (a-h, o-z)
   parameter (nPart=2) ! Muon+ or Muon-
   parameter (nsol=2) ! solar minimum and maximum
   parameter (nAdata=7) ! number of A parameter

   dimension Acurr(nAdata)
   dimension Fl(nsol),spot(nsol)

   data restmass/105.6/
   data spot/0.0,150.0/
   data ethre/3.0e5/ ! threshold energy for high-energy muon correction

   if(e.lt.1.0e-2) then
      getMuonSpec=0.0
      return
   endif

   r=max(1.0,r1) ! muon fluxes are the same for Rc < 1 GV
   d=max(0.15,d1) ! secondary particle fluxes are the same for d < 0.15 g/cm2

   beta=sqrt(1.0-(restmass/(restmass+e))**2)
   tmp=max(2.0,log10(e)) ! below 100 MeV, this value should be constant

   do is=1,nsol
      iptmp=ip ! 1 for mu+, 2 for mu+
      call getAmuon(Acurr,iptmp,is,d,r)
      if(e.gt.ethre) then
         acurr(5)=acurr(5)+0.4 ! /(1+exp((5.8088849-tmp)/0.24867240)) ! high energy correction, see fit/muon/PowerCorrection
         acurr(1)=acurr(1)*ethre**0.4 ! /(1+exp((5.8088849-tmp)/0.24867240))
      endif
      Fl(is)=acurr(1)*(e+(acurr(2)+acurr(4)*tmp)/beta**acurr(3))**(-acurr(5))*(1+exp(-acurr(6)*(log(e)+acurr(7))))
   enddo

   iptmp=ip+5 ! 6 for mu+, 7 for mu-
   pow=getPow(iptmp,d,r) ! get Power index
   A2=(Fl(1)-Fl(2))/(getFFPfromW(spot(1))**pow-getFFPfromW(spot(2))**pow)
   A1=Fl(1)-A2*getFFPfromW(spot(1))**pow
   getMuonSpec=a1+a2*getFFPfromW(s)**pow

   return
end

! ******************************************************
subroutine getAmuon(Acurr,ip,is,d,r) ! get A parameter
! ******************************************************
   implicit real*8 (a-h, o-z)
   parameter (npart=2) ! mu+ and mu-
   parameter (nsol=2) ! solar minimum and maximum
   parameter (nAdata=7) ! number of A parameter A(1) to A(7)
   parameter (nBdata=10) ! A_min = B1+B2*r+B3/(1+exp((r-B4)/B5), A_max = B6+B7*r+B8/(1+exp((r-B9)/B10)
   parameter (ndep=26) ! number of depth
   dimension Acurr(nAdata)
   dimension Bdata(nAdata,nBdata,npart,ndep),dep(ndep)
   character chatmp1*1,chatmp4*4
   character charge(npart)*4
   character chanum(0:9)*1
   save Bdata,dep

   data chanum/'0','1','2','3','4','5','6','7','8','9'/
   data charge/'plus','mins'/
   data ifirst/0/

   if(ifirst.eq.0) then  ! first time called this routine
      ifirst=1
      do ip2=1,npart
         open(15,file='input/muon--/final135.'//charge(ip2),status='old') ! read A(1),A(3),A(5), they are solar independent
         read(15,*) chatmp1
         do id2=1,ndep
            read(15,*) dep(id2),Bdata(1,1,ip2,id2),Bdata(3,1,ip2,id2),Bdata(5,1,ip2,id2)
         enddo
         close(15)
         open(15,file='input/muon--/final2467.'//charge(ip2),status='old') ! read A(2),A(4),A(6),A(7)
         read(15,*) chatmp1
         do id2=1,ndep
            read(15,*) dep(id2),(Bdata(2,ib,ip2,id2),ib=1,nBdata)
         enddo
         read(15,*) chatmp1
         do id2=1,ndep
            read(15,*) dep(id2),(Bdata(4,ib,ip2,id2),ib=1,nBdata)
         enddo
         read(15,*) chatmp1
         do id2=1,ndep
            read(15,*) dep(id2),(Bdata(6,ib,ip2,id2),ib=1,nBdata)
         enddo
         read(15,*) chatmp1
         do id2=1,ndep
            read(15,*) dep(id2),(Bdata(7,ib,ip2,id2),ib=1,nBdata)
         enddo
         close(15)
      enddo
   endif

! find closest depth
   do id=1,ndep
      if(d.lt.dep(id)) exit
   enddo
   if(id.eq.1) then
      ratio=0.0
      id=2
   elseif(id.eq.ndep+1) then
      ratio=1.0
      id=ndep
   else
      ratio=(d-dep(id-1))/(dep(id)-dep(id-1))
   endif

   rc=max(1.0,r) ! minimum Rc = 1.0GV

   do ia=1,nAdata
      if(ia.eq.1) then
         Acurr(ia)=log(Bdata(ia,1,ip,id-1))+(log(Bdata(ia,1,ip,id))-log(Bdata(ia,1,ip,id-1)))*ratio ! log-interpolation
         Acurr(ia)=exp(Acurr(ia))
      elseif(ia.eq.3.or.ia.eq.5) then
         Acurr(ia)=Bdata(ia,1,ip,id-1)+(Bdata(ia,1,ip,id)-Bdata(ia,1,ip,id-1))*ratio
      else
         if(is.eq.1) then
            B1=Bdata(ia,1,ip,id-1)+Bdata(ia,2,ip,id-1)*rc+Bdata(ia,3,ip,id-1)/(1+exp((rc-Bdata(ia,4,ip,id-1))/Bdata(ia,5,ip,id-1)))
            B2=Bdata(ia,1,ip,id)+Bdata(ia,2,ip,id)*rc+Bdata(ia,3,ip,id)/(1+exp((rc-Bdata(ia,4,ip,id))/Bdata(ia,5,ip,id)))
         else
            B1=Bdata(ia,6,ip,id-1)+Bdata(ia,7,ip,id-1)*rc+Bdata(ia,8,ip,id-1)/(1+exp((rc-Bdata(ia,9,ip,id-1))/Bdata(ia,10,ip,id-1)))
            B2=Bdata(ia,6,ip,id)+Bdata(ia,7,ip,id)*rc+Bdata(ia,8,ip,id)/(1+exp((rc-Bdata(ia,9,ip,id))/Bdata(ia,10,ip,id)))
         endif
         Acurr(ia)=B1+(B2-B1)*ratio
      endif
   enddo

   return
end


! *******************************************************
function getIonSpec(iz,s,r,d,e) ! get Ion flux
!     s:Wolf number
!     r:cut off rigidity (GV)
!     d:air depth (g/cm^2)
!     e:ion energy (MeV/n)
! *******************************************************
   implicit real*8 (a-h, o-z)
   parameter(nAdata=6)
   parameter(npart=28)
   parameter(ngroup=6)

   character chatmp1*1

   real*8, save:: A(nAdata,ngroup) ! parameter used in combine.for
   real*8, save:: dEdx(npart)
   integer*4, save:: iAnum(npart)
   integer*4, save:: ifirst

   integer, save:: igidx(npart)    ! group index

   data igidx/1,2,3,3,3,4,4,4,4,5 &
   & ,5,5,5,5,5,5,5,5,5,6 &
   & ,6,6,6,6,6,6,6,6/

   data ianum/ 1, 4, 7, 9,11,12,14,16,19,20,23,24,27,28,31,32,35,40,39,40,45,48,51,52,55,56,59,59/

   data ifirst/0/
   data restmass/938.27d0/

   if(ifirst.eq.0) then  ! first time call
      ifirst=1
      open(unit=15,file='input/ions/Combine.inp',status='old')
      read(15,*) chatmp1
      do i=1,ngroup
         read(15,*) (A(ia,i),ia=1,nAdata)
      enddo
      close(unit=15)
   endif

   if(e.lt.1.0e-2) then  ! for lower energy, no output
      getIonSpec=0.0
      return
   endif

   x=e ! x is energy
   ig=igidx(iz)

   tmp=restmass*iAnum(iz)

   Ecut=getEfromR(iZ,tmp,r*1000.0)/iAnum(iz)-a(6,ig)*d  ! MeV/n
   EcPri=max(a(1,ig),Ecut*a(3,ig))
   EcSec=max(a(2,ig),Ecut*a(3,ig))

   if(iz.le.2) then
      ip=iz
   else
      ip=iz+3 ! in getsecondary, ip=iz+3 (electron, positron, photon) for ions
   endif

   getIonSpec=getPrimary(iz,iAnum(iz),s,d,x)*0.5*(tanh(a(4,ig)*(x/EcPri-1))+1.0) &
   & +getsecondary(ip,s,r,d,x)*0.5*(tanh(a(5,ig)*(1-x/EcSec))+1.0)

   return

end

! **********************************************************
function getsecondary(ip,s,r1,d1,e) ! get secondary particle flux (/cm^2/s/MeV)
! **********************************************************
   parameter(nBdata=8)
   parameter(ndep=26)
   parameter(npart=11) ! proton, alpha, electron, positron, photon, Li, Be, B, C, N, O
   implicit real*8 (a-h, o-z)
   real*8, save:: A(nBdata,npart,ndep),B(nBdata)
   real*8, save:: dep(ndep)
   integer, save:: ifirst
   character chatmp1*1,chatmp5*5

   data ifirst/0/

   r=max(1.0,r1) ! to get Fl, minimum Rc=1GV
   d=max(0.15,d1) ! secondary particle fluxes are the same for d < 0.15 g/cm2

   if(ifirst.eq.0) then
      do i=1,5 ! proton, alpha, electron, positron, photon
         if(i.eq.1) then
            open(unit=15,file='input/proton/fitting-lowspec.inp',status='old')
         elseif(i.eq.2) then
            open(unit=15,file='input/alphaa/fitting-lowspec.inp',status='old')
         elseif(i.eq.3) then
            open(unit=15,file='input/elemag/fitting-lowspec-EL.inp',status='old')
         elseif(i.eq.4) then
            open(unit=15,file='input/elemag/fitting-lowspec-PO.inp',status='old')
         elseif(i.eq.5) then
            open(unit=15,file='input/elemag/fitting-lowspec-PH.inp',status='old')
         endif
         read(15,'(A)') chatmp1
         do id=1,ndep
            read(15,*) dep(id),(A(ib,i,id),ib=1,nBdata)
         enddo
         close(unit=15)
      enddo
      open(15,file='input/ions/fitting-lowspec.inp',status='old')
      read(15,'(a)') chatmp1
      do i=npart-5,npart
         read(15,*) (A(ib,i,1),ib=1,nBdata)
         do ib=1,nBdata
            do id=1,ndep
               A(ib,i,id)=A(ib,i,1) ! for Li to O, depth independent
            enddo
         enddo
      enddo
      close(15)
      ifirst=1
   endif

   if(e.lt.1.0e-2.or.ip.gt.npart) then  ! for lower energy or heavier ions, no output
      getsecondary=0.0
      return
   endif

   do id=1,ndep
      if(d.lt.dep(id)) exit
   enddo
   if(id.eq.1) then
      ratio=0.0
      id=2
   elseif(id.eq.ndep+1) then
      ratio=1.0
      id=ndep
   else
      ratio=(d-dep(id-1))/(dep(id)-dep(id-1))
   endif

   do ib=1,nBdata
      B(ib)=A(ib,ip,id-1)+(A(ib,ip,id)-A(ib,ip,id-1))*ratio
   enddo

   if(ip.eq.1.or.ip.eq.2) then ! proton or alpha
      getsecondary=getFl(ip,s,r,d)*(b(1)*e**b(2))/(1+b(3)*e**b(4))/(1+b(5)*e**b(6))*(1+exp(-b(7)*(log(e)+b(8))))
   elseif(ip.eq.3.or.ip.eq.4) then ! electron or positron
      getsecondary=getFl(ip,s,r,d)*(b(1)*e**b(2))/(1+b(3)*e**b(4))/(1+b(5)*e**b(6))
   elseif(ip.eq.5) then ! photon
      getsecondary=getFl(ip,s,r,d)*(b(1)*e**b(2))*(1+b(3)*e**b(4))/(1+b(5)*e**b(6))/(1+exp(-b(7)*(log(e)+b(8))))
   else ! Li,Be,B,C,N,O
      getsecondary=getFl(ip,s,r,d)*(b(1)*e**b(2))/(1+b(3)*e**b(4))/(1+b(5)*e**b(6))/(1+exp(-b(7)*(log(e)+b(8))))
   endif

   return

end

! **********************************************************
function getTOAspec(iZ,iA,Ek,Spot) ! get TOA spectrum in (/(MeV/n)/s/m^2/sr)
!     Ek: Kinetic Energy in MeV/n
!     Spot: Wolf number estimated from count rate of neutron monitors
! **********************************************************
   parameter(npart=28)
   implicit real*8 (a-h,o-z)
   real*8, save:: Dpara(npart),alpha(npart),gamma(npart),bpara(npart) ! Data for Nymmik Model
   data Emp/938.27d0/ ! mass of proton, nucleus mass is simply assumed to be A*Emp

   data Dpara/1.85e4,3.69e3,19.50,17.70,49.20,103.00,36.70,87.40, &
   &           3.19,16.40,4.43,19.30,4.17,13.40,1.15,3.06,1.30, &
   &           2.33,1.87,2.17,0.74,2.63,1.23,2.12,1.14,9.32,0.10,0.48/ ! ISO-Model taken from Matthia-ASR2013

   data alpha/2.85,3.12,3.41,4.30,3.93,3.18,3.77,3.11,4.05,3.11,3.14, &
   &           3.65,3.46,3.00,4.04,3.30,4.40,4.33,4.49,2.93,3.78,3.79, &
   &           3.50,3.28,3.29,3.01,4.25,3.52/

   data gamma/2.74,2.77,2.82,3.05,2.96,2.76,2.89,2.70,2.82,2.76,2.84, &
   &           2.70,2.77,2.66,2.89,2.71,3.00,2.93,3.05,2.77,2.97,2.99, &
   &           2.94,2.89,2.74,2.63,2.63,2.63/

! ***** Determine LIS spectra based on DLR Model *****************
   R=getRfromE(iz,ia,Ek,Emp*iA)*0.001 ! rigidity in GV
   beta=sqrt(1-(Emp*iA/(Emp*iA+Ek*iA))**2)
   dR2dE=0.001/iZ/beta*iA ! convert GV to MV, MeV to MeV/n
   SpecLIS=Dpara(iz)*beta**alpha(iz)/R**gamma(iz)*dR2dE
! *******************************************************************

! consider solar modulation
   R0=getFFPfromW(Spot)*0.001 ! FFP in GV from W, taken from Matthia-ASR2013
   delta=0.02*Spot+4.7        ! taken from Matthia-ASR2013
   getTOAspec=SpecLIS*(R/(R+R0))**delta

   return
end


! *************************************************************
function getPrimary(iz,ia,s,d,e)  ! get Primary Flux
! *************************************************************
   parameter(nAdata=3)
   parameter(npart=28)
   parameter(nepoint=6)
   parameter(ngroup=6)
   parameter(ndep=26)
   implicit real*8 (a-h, o-z)
   real*8, save:: A(nAdata,ngroup,ndep)
   real*8, save:: dEdxTable(npart,nepoint),epoint(nepoint) ! dE/dx for each particle
   real*8, save:: dep(ndep) ! depth (g/cm2)
   real*8, save:: down
   integer, save:: igidx(npart)    ! group index
   dimension B(nAdata) ! temporary used dimension

   character gname(ngroup)*2
   character chatmp1*1

   data down/0.0/

   data gname/'H-','He','Be','N-','Si','Fe'/

   data igidx/1,2,3,3,3,4,4,4,4,5 &
   & ,5,5,5,5,5,5,5,5,5,6 &
   & ,6,6,6,6,6,6,6,6/

   if(down.eq.0.0) then ! first time call this routine
      down=1.720313d0
      do ig=1,ngroup
         open(unit=15,file='input/ions/primary-'//gname(ig)//'.inp',status='old')
         read(15,*) chatmp1
         do id=1,ndep
            read(15,*) dep(id),(a(i,ig,id),i=1,nAdata)
         enddo
         close(15)
      enddo
      open(15,file='input/ions/dEdx-table.inp',status='old')
      read(15,'(a1)') chatmp1
      read(15,'(a1)') chatmp1
      read(15,*) (epoint(ie),ie=1,nepoint)
      do ip=1,npart
         read(15,*) itmp,itmp,(dEdxTable(ip,ie),ie=1,nepoint)
      enddo
      close(15)
   endif

   ig=igidx(iz)

   do id=1,ndep
      if(d.lt.dep(id)) exit
   enddo
   if(id.eq.1) then
      ratio=0.0
      id=2
   elseif(id.eq.ndep+1) then
      ratio=1.0
      id=ndep
   else
      ratio=(d-dep(id-1))/(dep(id)-dep(id-1))
   endif

   do i=1,nAdata
      b(i)=a(i,ig,id-1)+(a(i,ig,id)-a(i,ig,id-1))*ratio
   enddo

! find dEdx
   do ie=2,nepoint-1
      if(e.le.epoint(ie)) exit
   enddo
   ratio=(log(e)-log(epoint(ie-1)))/(log(epoint(ie))-log(epoint(ie-1))) ! log-logg interpolation
   ratio=max(0.0,min(1.0,ratio))
   tmp=log(dEdxTable(iz,ie-1))+(log(dEdxTable(iz,ie))-log(dEdxTable(iz,ie-1)))*ratio ! log-log interpolation
   dEdx=exp(tmp)

   Eini=e+dEdx*d  ! Energy at the TOA
   getPrimary=getTOAspec(iz,ia,Eini,s)*(b(1)*exp(-b(2)*d)+(1.0-b(1))*exp(-b(3)*d))
   getPrimary=getPrimary*4.0*acos(-1.0)*1.0e-4/down  ! convert (/(MeV/n)/s/m^2/sr) to (/(MeV/n)/s/cm^2)

   return

end

! *************************************************************
function get511flux(s,r,d)  ! get 511 keV photon flux in (/cm2/s)
! *************************************************************
   parameter(ndep=26)
   implicit real*8 (a-h, o-z)
   real*8, save:: F511(ndep),dep(ndep)
   character chatmp1*1

   data ifirst/0/

   if(ifirst.eq.0) then
      ifirst=1
      open(15,file='input/elemag/flux511keV.inp')
      read(15,'(a1)') chatmp1
      do id=1,ndep
         read(15,*) dep(id),F511(id)
      enddo
      close(15)
   endif

   do id=1,ndep
      if(d.lt.dep(id)) exit
   enddo
   if(id.eq.1) then
      ratio=0.0
      id=2
   elseif(id.eq.ndep+1) then
      ratio=1.0
      id=ndep
   else
      ratio=(d-dep(id-1))/(dep(id)-dep(id-1))
   endif

   Fratio=F511(id-1)+(F511(id)-F511(id-1))*ratio

   iptmp=5 ! photon index
   ene=0.511 ! 511 keV
   get511flux=getsecondary(iptmp,s,r,d,ene)*Fratio

   return

end


! ******************************************************
function getd_model(alti,cido,iMSIS,istat) ! atmospheric depth in g/cm^2, altitude in km
! iMSIS=0: US Standard Atmosphere 1976; iMSIS=1: bundled NRLMSISE-00 table
! istat=0: success, 1: altitude too low, 2: altitude too high,
!       3: table I/O error, 4: invalid model
! ******************************************************
   parameter(maxUS=75) ! number of altitude bins for US Standard Atmosphere
   parameter(maxMSIS=129) ! number of altitude bins for NRLMSISE-00
   parameter(maxlat=36) ! number of latitude bins for NRLMSISE-00
   implicit real*8 (a-h, o-z)
   real*8, save:: altUS(maxUS)
   real*8, save:: altMSIS(maxMSIS)
   real*8, save:: depUS(maxUS)
   real*8, save:: depMSIS(maxMSIS,maxlat)
   real*8, save:: glat(maxlat)
   logical, save:: initialized=.false.
   character chatmp1*1
   integer ios

   getd_model=-1.0d0
   istat=0

   if(.not.initialized) then
      open(unit=15,file='input/AtomDepth.inp',status='old',action='read',iostat=ios)
      if(ios.ne.0) then
         istat=3
         return
      endif
      read(15,'(a)',iostat=ios) chatmp1
      if(ios.ne.0) goto 900
      do ia=1,maxUS
         read(15,*,iostat=ios) altUS(ia),depUS(ia)
         if(ios.ne.0) goto 900
      enddo
      read(15,'(a)',iostat=ios) chatmp1
      if(ios.ne.0) goto 900
      read(15,*,iostat=ios) (glat(ido),ido=1,maxlat)
      if(ios.ne.0) goto 900
      do ia=1,maxMSIS
         read(15,*,iostat=ios) altMSIS(ia),(depMSIS(ia,ido),ido=1,maxlat)
         if(ios.ne.0) goto 900
      enddo
      close(unit=15)
      initialized=.true.
   endif

   if(iMSIS.eq.0) then
      if(alti.lt.altUS(1)) then
         istat=1
         return
      elseif(alti.gt.altUS(maxUS)) then
         istat=2
         return
      endif
      do ia=2,maxUS
         if(alti.le.altUS(ia)) exit
      enddo
      ratio=(alti-altUS(ia-1))/(altUS(ia)-altUS(ia-1))
      getd_model=depUS(ia-1)+ratio*(depUS(ia)-depUS(ia-1))
   elseif(iMSIS.eq.1) then
      if(alti.lt.altMSIS(1)) then
         istat=1
         return
      elseif(alti.gt.altMSIS(maxMSIS)) then
         istat=2
         return
      endif
      do ia=2,maxMSIS
         if(alti.le.altMSIS(ia)) exit
      enddo
      ratio=(alti-altMSIS(ia-1))/(altMSIS(ia)-altMSIS(ia-1))

      if(cido.le.glat(1)) then
         ido=2
         ratio1=0.0d0
      elseif(cido.ge.glat(maxlat)) then
         ido=maxlat
         ratio1=1.0d0
      else
         do ido=2,maxlat
            if(cido.le.glat(ido)) exit
         enddo
         ratio1=(cido-glat(ido-1))/(glat(ido)-glat(ido-1))
      endif
      dep1=depMSIS(ia-1,ido-1)+ratio*(depMSIS(ia,ido-1)-depMSIS(ia-1,ido-1))
      dep2=depMSIS(ia-1,ido  )+ratio*(depMSIS(ia,ido  )-depMSIS(ia-1,ido  ))
      getd_model=dep1+ratio1*(dep2-dep1)
   else
      istat=4
      return
   endif
   return

900 continue
   close(unit=15,iostat=ios)
   istat=3
   getd_model=-1.0d0
   return
end

! ******************************************************
function getr(cido,ckei) ! cido and ckei are center of ido&keido of each grid
! ******************************************************
   implicit real*8 (a-h, o-z)
   real*8, save:: cordata(181,361) ! maximum 1 deg step, +1 mean
   real*8, save:: dpido(181),dpkei(361) ! data point ido & keido
   character chatmp1*1
   integer*4, save:: mkei,mido
   real*8, save:: skei,sido
   integer*4, save:: ifirst
   data ifirst/0/

   if(ifirst.eq.0) then ! come to this routine first, so read input data
      ifirst=1
      open(unit=15,file='input/CORdata.inp',status='old')
      read(15,*) mkei,mido ! read step size
      read(15,*) chatmp1
      skei=360.0/(mkei-1)  ! keido step
      sido=180.0/(mido-1)  ! ido step
      do id=mido,1,-1     ! read from 90 to -90 deg
         do ik=mkei,1,-1      ! read from 180 to -180 deg
            read(15,*) dpkei(ik),dpido(id),cordata(id,ik)
         enddo
      enddo
      close(unit=15)
   endif

! **** Determine Cut-off Rigidity ***********
   if(cido.lt.-90.0d0.or.cido.gt.90.0d0.or.ckei.lt.-180.0d0.or.ckei.gt.180.0d0) then
      getr=-1.0d0
      return
   endif
   id=max(1,min(mido-1,int((cido+90.0d0)/sido)+1))   ! lower latitude bin
   ik=max(1,min(mkei-1,int((ckei+180.0d0)/skei)+1))  ! lower longitude bin
   cor1=cordata(id,ik)*(dpido(id+1)-cido)/(dpido(id+1)-dpido(id))+cordata(id+1,ik)*(cido-dpido(id))/(dpido(id+1)-dpido(id))
   cor2=cordata(id,ik+1)*(dpido(id+1)-cido)/(dpido(id+1)-dpido(id))+cordata(id+1,ik+1)*(cido-dpido(id))/(dpido(id+1)-dpido(id))
   getr=cor1*(dpkei(ik+1)-ckei)/(dpkei(ik+1)-dpkei(ik))+cor2*(ckei-dpkei(ik))/(dpkei(ik+1)-dpkei(ik))

   return

end


function getHP(iy0,im0,id0,ic) ! get solar W index from bundled tables
!   ic=1: daily neutron-monitor data
!   ic=2: annual Usoskin/Wolf-number data
!   ic=3: suspected ground-level event (encoded offset removed)
!   ic=4: no data for the requested year/date
!   ic=5: invalid calendar date
!   ic=6: table I/O or format error
   implicit real*8 (a-h, o-z)
   integer, parameter:: nmonth=12, nday=31
   integer, parameter:: iymin=1500, iymax=2500
   real*8, parameter:: missing=-1.0d30
   integer iy0,im0,id0,ic
   integer iy,im,id,ios,itmp,itmp1,itmp2,maxday
   integer, save:: iystart=0,iyend=-1,iysUs=0,iyeUs=-1
   integer, dimension(12):: days_in_month
   real*8, save:: FFP(iymin:iymax,nmonth,nday)
   real*8, save:: FFPuso(iymin:iymax)
   logical, save:: initialized=.false.,tables_ok=.false.
   logical leap_year

   getHP=0.0d0

   if(.not.initialized) then
      FFP=missing
      FFPuso=missing
      tables_ok=.false.

      open(unit=15,file='input/FFPtable.day',status='old',action='read',iostat=ios)
      if(ios.ne.0) goto 900
      read(15,*,iostat=ios) iystart,iyend
      if(ios.ne.0.or.iystart.lt.iymin.or.iyend.gt.iymax.or.iystart.gt.iyend) goto 900
      do im=1,nmonth
         do id=1,nday
            read(15,*,iostat=ios) itmp1,itmp2,(FFP(iy,im,id),iy=iystart,iyend)
            if(ios.ne.0.or.itmp1.ne.im.or.itmp2.ne.id) goto 900
         enddo
      enddo
      close(unit=15,iostat=ios)

      open(unit=15,file='input/FFPtable.uso',status='old',action='read',iostat=ios)
      if(ios.ne.0) goto 900
      read(15,*,iostat=ios) iysUs,iyeUs
      if(ios.ne.0.or.iysUs.lt.iymin.or.iyeUs.gt.iymax.or.iysUs.gt.iyeUs) goto 900
      do iy=iysUs,iyeUs
         read(15,*,iostat=ios) itmp,FFPuso(iy)
         if(ios.ne.0.or.itmp.ne.iy) goto 900
      enddo
      close(unit=15,iostat=ios)

      tables_ok=.true.
900   continue
      if(.not.tables_ok) close(unit=15,iostat=ios)
      initialized=.true.
   endif

   if(.not.tables_ok) then
      ic=6
      return
   endif

! ****** Validate the calendar date **************
   if(iy0.lt.1.or.im0.lt.1.or.im0.gt.nmonth) then
      ic=5
      return
   endif
   leap_year=(mod(iy0,4).eq.0.and.mod(iy0,100).ne.0).or.mod(iy0,400).eq.0
   days_in_month=(/31,28,31,30,31,30,31,31,30,31,30,31/)
   if(leap_year) days_in_month(2)=29
   maxday=days_in_month(im0)
   if(id0.lt.1.or.id0.gt.maxday) then
      ic=5
      return
   endif

   if(iy0.lt.iymin.or.iy0.gt.iymax) then
      ic=4
      return
   endif

! ****** Prefer daily neutron-monitor data *************
   if(iy0.ge.iystart.and.iy0.le.iyend) then
      if(FFP(iy0,im0,id0).gt.-99.0d0) then
         getHP=FFP(iy0,im0,id0)
         if(getHP.gt.1000.0d0) then
            ic=3
            getHP=getHP-10000.0d0
         else
            ic=1
         endif
         return
      endif
   endif

! ****** Fall back to annual Usoskin data *************
   if(iy0.ge.iysUs.and.iy0.le.iyeUs) then
      if(FFPuso(iy0).gt.missing/2.0d0) then
         getHP=FFPuso(iy0)
         ic=2
         return
      endif
   endif

! ****** No data **************************
   ic=4
   getHP=0.0d0
   return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! module contains several utility functions

module utility

contains

   subroutine create_log_grid(min_val, max_val, grid_out, npoints)
      implicit none
      integer, intent(in) :: npoints
      real(kind=8), intent(in) :: min_val, max_val
      real(kind=8), intent(out), dimension(npoints) :: grid_out
      integer :: i
      real(kind=8) :: log_emin, log_emax, step

      if (npoints < 2) error stop 'create_log_grid requires at least two points'
      if (min_val <= 0.0d0 .or. max_val <= min_val) then
         error stop 'create_log_grid requires 0 < min_val < max_val'
      end if

      ! Compute the logarithms of the min and max energies.
      log_emin = log10(min_val)
      log_emax = log10(max_val)
      step = (log_emax - log_emin) / (real(npoints, kind=8) - 1.0d0)

      do i = 1, npoints
         grid_out(i) = 10.0d0**(log_emin + (real(i, kind=8) - 1.0d0) * step)
      end do
   end subroutine create_log_grid

   subroutine ID_to_string_mapping(n, label)
      implicit none
      integer, intent(in) :: n      ! Input integer ID
      character(len=*), intent(out) :: label   ! Output string label

      ! Associate integers with strings
      select case (n)
       case (0)
         label = 'neutron'
       case (1:28)
         label = 'H-Ni ions'
       case (29)
         label = 'muon+'
       case (30)
         label = 'muon-'
       case (31)
         label = 'electron'
       case (32)
         label = 'positron'
       case (33)
         label = 'photon'
       case default
         label = 'undefined'
      end select
   end subroutine ID_to_string_mapping

end module utility

! module contains several integration functions

module integration

contains

   function simpson_integration(energy_grid, flux_values, num_points) result(total_flux)
      implicit none
      integer, intent(in) :: num_points
      real(kind=8), dimension(num_points), intent(in) :: energy_grid, flux_values
      real(kind=8) :: total_flux, h0, h1, span
      integer :: i, last_simpson

      if (num_points < 2) error stop 'simpson_integration requires at least two points'
      total_flux = 0.0d0
      last_simpson = num_points
      if (mod(num_points, 2) == 0) last_simpson = num_points - 1

      ! Integrate each pair of unequal intervals using the exact integral of
      ! the quadratic Lagrange interpolant through the three points.
      do i = 1, last_simpson - 2, 2
         h0 = energy_grid(i+1) - energy_grid(i)
         h1 = energy_grid(i+2) - energy_grid(i+1)
         if (h0 <= 0.0d0 .or. h1 <= 0.0d0) error stop 'integration grid must increase'
         span = h0 + h1
         total_flux = total_flux + span / 6.0d0 * ( &
            (2.0d0 - h1/h0) * flux_values(i) + &
            (span*span/(h0*h1)) * flux_values(i+1) + &
            (2.0d0 - h0/h1) * flux_values(i+2))
      end do

      ! An even number of points leaves one final interval; do not omit it.
      if (last_simpson < num_points) then
         h0 = energy_grid(num_points) - energy_grid(num_points-1)
         if (h0 <= 0.0d0) error stop 'integration grid must increase'
         total_flux = total_flux + 0.5d0*h0 * &
            (flux_values(num_points-1) + flux_values(num_points))
      end if
   end function simpson_integration

   function log_simpson_integration(energy_grid, flux_values, num_points) result(total_flux)
      implicit none
      integer, intent(in) :: num_points
      real(kind=8), dimension(num_points), intent(in) :: energy_grid, flux_values
      real(kind=8) :: total_flux, log_step, weighted_sum, expected_log, tolerance
      integer :: i

      if (num_points < 3 .or. mod(num_points, 2) == 0) then
         error stop 'log_simpson_integration requires an odd number of at least three points'
      end if
      if (energy_grid(1) <= 0.0d0) error stop 'log integration requires positive energies'

      log_step = (log(energy_grid(num_points)) - log(energy_grid(1))) / &
         real(num_points - 1, kind=8)
      if (log_step <= 0.0d0) error stop 'integration grid must increase'
      tolerance = 1.0d3 * epsilon(1.0d0) * max(1.0d0, abs(log(energy_grid(num_points))))

      do i = 2, num_points
         if (energy_grid(i) <= energy_grid(i-1)) error stop 'integration grid must increase'
         expected_log = log(energy_grid(1)) + real(i - 1, kind=8) * log_step
         if (abs(log(energy_grid(i)) - expected_log) > tolerance) then
            error stop 'log_simpson_integration requires a logarithmically uniform grid'
         end if
      end do

      ! With x=ln(E), dE=E dx.  Simpson's rule is therefore applied to
      ! flux(E)*E on the uniformly spaced x grid.
      weighted_sum = flux_values(1)*energy_grid(1) + &
         flux_values(num_points)*energy_grid(num_points)
      do i = 2, num_points - 1
         if (mod(i, 2) == 0) then
            weighted_sum = weighted_sum + 4.0d0*flux_values(i)*energy_grid(i)
         else
            weighted_sum = weighted_sum + 2.0d0*flux_values(i)*energy_grid(i)
         end if
      end do
      total_flux = log_step * weighted_sum / 3.0d0
   end function log_simpson_integration

end module integration
