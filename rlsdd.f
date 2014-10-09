      real function rlsdd(enu,qsq,nubar,ma,force)
      implicit none
c---------------------------------------------------------------------------
c
c This function returns dsigma/dQ**2 from Llewellyn-Smith
c
c note: - Q**2 = -q**2 >0
c       - xsi = mu_p-mu_n = 1.793 - (-1.913) = 3.706 
c       - neglects psuedo-scalar form factor F_P
c
c 09/20/05: install kinematic constraint Q^2max=(4*enu**2)/(1+(2*enu)/mn)
c
c---------------------------------------------------------------------------
      real gm, ge, fv1, fv3, fv2, fa,fd, fp, fa3,gpi,gpi0,fp1,fpm,fpe
      real xsi, fa0, mv, mn0, mpi, mn, mmu, cost, gf, pi, hbc, me, ml
      real arb,mlmax,mlmin
      parameter (xsi=3.706, fa0=-1.267, mv=0.84,
c     +           mpi=0.1396, mn=0.9389, mmu=0.10566,
     +           mpi=0.1396, mn=0.9389, mmu=1.777,
c     +           mpi=0.1396, mn=0.9389, mmu=0.000511,
     +           cost=0.9738, gf=1.16637e-5, pi=3.14159, me=.000511,
c Uncomment for Normal Case
c     +           hbc=1.97327e-16, ml = .000511, mn0 = .9383) 
c     +           hbc=1.97327e-16, ml = .10566, mn0 = .9383) 
     +           hbc=1.97327e-16, ml = 1.777, mn0 = .9383) 
c Uncomment to vary min/max independently 
c     +           hbc=1.97327e-16, mlmin = .000511, 
c     +           mlmax = .10566, mn0 = .9383) 
c Declare Variables
      real a,b,c,su,tau 
      real enu,qsq,dsigma,qsqmax,qsqmin
      real s,enucm,elepcmin,plepcmin,rts,minmin,addmin
      real elepcmax,plepcmax,minmax,addmax
      real faboundary,ma
      logical nubar
      logical first
      logical force
      data first /.true./
      integer ipol

c Uncomment for Normal Case
      mlmin = ml
      mlmax = ml

c
c..initialization 
      rlsdd=0.
      if (first) then
        first=.false.
        print *,''
        print *,'RLSDD: xsi,ma,mv,fa0,',xsi,ma,mv,fa0
      endif
c
      tau=qsq/4./mn**2
      su = 4.*mn*enu-qsq-mmu**2
c
      s = 2.*enu*mn0+mn0*mn0
      rts = sqrt(s)/2.
      addmax = (1-((mlmax+mn0)*(mlmax+mn0))/s)
      addmin = (1-((mlmin+mn0)*(mlmin+mn0))/s)
      minmax = (1-((mlmax-mn0)*(mlmax-mn0))/s)
      minmin = (1-((mlmin-mn0)*(mlmin-mn0))/s)
      enucm = (s-mn0*mn0)/(2.*sqrt(s))
      elepcmax = (s+mlmax*mlmax-mn0*mn0)/(2.*sqrt(s))
      elepcmin = (s+mlmin*mlmin-mn0*mn0)/(2.*sqrt(s))
      plepcmax = rts*sqrt(addmax*minmax)
      plepcmin = rts*sqrt(addmin*minmin)
c      print *,plepcmax,elepcmax,sqrt(elepcmax**2-plepcmax**2)
c
      qsqmin = -mlmin*mlmin+2.*(enucm*elepcmin-enucm*plepcmin)
      qsqmax= (-mlmax*mlmax+2.*(enucm*elepcmax+enucm*plepcmax))
c      print *,'first',qsqmin,qsqmax,qsq,enu
c      print *,'second',rts,addm,minm,enucm,elepcm,plepcm
      if (enu.le.mlmax) then
        qsqmax = 0
      endif
c      print *,qsqmax,qsqmin
      ge = 1./(1.+qsq/mv**2)**2
      fpe = 1./(1.+qsq/mpi**2)**2
      gm = (1.+xsi)/(1.+qsq/mv**2)**2
      fpm = (1.+xsi)/(1.+qsq/mpi**2)**2
      fv1 = 1./(1.+tau)*(ge+tau*gm)
      fp1 = 1./(1.+tau)*(fpe+tau*fpm)
      fv2 = 1./(1.+tau)*(gm-ge)
      if (ma.le.0) then
         fa = 0.
      else
         fa = fa0/(1.+qsq/ma**2)**2
      endif
      fd = 1/(1.+qsq/(mv)**2)**2

c Model Dependent Variables
c GT unc   
c      fp = (2*(1.03)*fa*mn**2)/(mpi**2+qsq)
c Normal
      fp = (2*fa*mn**2)/((mpi)**2+qsq)
c mpi = .6mpi
c      fp = (2*fa*mn**2)/((.6*mpi)**2+qsq)
c mpi = 1.5mpi
c      fp = (2*fa*mn**2)/((1.5*mpi)**2+qsq)
c Full
c      gpi0 = .97*fa0 
c      if (qsq.le..6) then
c        gpi = (((.03/.6)*qsq)+.97)*fa  
c      else
c        gpi = fa 
c      endif
c      fp = (2*fa0*mn**2/(qsq))*((fa/fa0)
c     +    -(gpi/(gpi0*(1+(qsq/(mpi*mpi))))))
c Dipole
c      fp = (-(1.414*13.4*.092)/(mpi**2))/(1.+qsq/(mpi)**2)**2
c Dipole f1v
c       fp = fp1*(-1.414*13.4*.092/(mpi**2)) 
c Axial Second Class
c      fa3 = (.1*fa*(xsi))/(fa0*2)
c      fa3 = (.32*fa)/(-2)
c      fa3 = 0
c      fp = 0
c
c Vector Second Class
c      fv3 = (mn/me)*(.0024*fd)

      a=(mmu**2+qsq)/mn**2*(
     +  (1.+tau)*abs(fa)**2- (1.-tau)*abs(fv1)**2 
     +  + tau*(1.-tau)*abs(fv2)**2 + 4.*tau*fv1*fv2
     +  - 4.*tau*(1.+tau)*abs(fa3)**2
     +  - (mmu/mn/2.)**2*(abs(fv1+fv2)**2 + abs(fa+2*fp)**2
     +  -4.*(1.+tau)*(abs(fv3**2)+abs(fp**2))))
      b=(qsq/mn**2*fa*(fv1+fv2)+((mmu/mn)**2)*((fa-2.*tau*fp)*fa3 
     + - (fv1 -tau*fv2)*fv3))
      c=(abs(fa)**2+abs(fv1)**2+tau*abs(fv2)**2 + 4.*tau*abs(fa3)**2)/4.
c      
c      print *,qsq,a,b,c
      ipol=1     
      if (nubar) ipol=-1
c   
      if ((qsq.le.qsqmax).or.force) then
        if ((qsq.gt.qsqmin).or.force) then
        arb = abs(fa+2*fp)**2-4.*(1.+tau)*abs(fp**2)
c        print *,enu,tau,abs(fa+2*fp)**2,4.*(1.+tau)*abs(fp**2),arb
        rlsdd = 10000*mn**2*(gf*hbc)**2*cost**2/8./pi/enu**2 * 
     +          (a-ipol*b*su/mn**2+c*su**2/mn**4)
        endif
      endif

c      if (enu.eq.4.5) then
        print *,enu,a,b,c,qsq,qsqmin,qsqmax,rlsdd
c      endif
c
      return
      end
