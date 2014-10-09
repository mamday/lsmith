      real function rlsdd(enu,qsq,nubar,ma)
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
      real gm, ge, fv1, fv2, fa, fp
      real xsi, fa0, mv, mpi, mn, mmu, cost, gf, pi, hbc
      parameter (xsi=3.706, fa0=-1.267, mv=0.84,
     +           mpi=0.1396, mn=0.9389, mmu=0.10566,
     +           cost=0.9738, gf=1.16637e-5, pi=3.14159,
     +           hbc=1.97327e-16) 
      real a,b,c,su,tau 
      real enu,qsq,dsigma,qsqmax
      real faboundary,ma
      logical nubar
      logical first
      data first /.true./
      integer ipol
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
      qsqmax=(4.*enu**2)/(1.+(2.*enu)/mn)
c
      ge = 1./(1.+qsq/mv**2)**2
      gm = (1.+xsi)/(1.+qsq/mv**2)**2
      fv1 = 1./(1.+tau)*(ge+tau*gm)
      fv2 = 1./(1.+tau)*(gm-ge)
      if (ma.le.0) then
         fa = 0.
      else
         fa = fa0/(1.+qsq/ma**2)**2
      endif
c
      a=(mmu**2+qsq)/mn**2*(
     +  (1.+tau)*abs(fa)**2- (1.-tau)*abs(fv1)**2 
     +  + tau*(1.-tau)*abs(fv2)**2 + 4.*tau*fv1*fv2
     +  - (mmu/mn/2.)**2*(abs(fv1+fv2)**2 + abs(fa)**2))
      b=qsq/mn**2*fa*(fv1+fv2)
      c=(abs(fa)**2+abs(fv1)**2+tau*abs(fv2)**2)/4.
c      
      ipol=1      
      if (nubar) ipol=-1
c   
      if (qsq.le.qsqmax) then
        rlsdd = mn**2*(gf*hbc)**2*cost**2/8./pi/enu**2 * 
     +          (a-ipol*b*su/mn**2+c*su**2/mn**4)
      endif
c
      return
      end
