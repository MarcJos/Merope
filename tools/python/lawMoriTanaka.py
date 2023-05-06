#!/usr/bin/python
# 2021
# Homogenization of two-phases isotropic composites, general relations
# + implementation of the Mori-tanaka estimate for isotropic two-phases solids
# + next step: nonlinear problem
# 
#
# Copyright : see License.txt


from math import *

# Usefull functions for elastic moduli (isotropic case)
# elas[0] = k, elas[1] = mu
def PoissonYoung(elas):
   res=[]
   Poisson = (3*elas[1]-2.0*elas[0])/(2.0*(3.0*elas[1]+elas[0]))
   res.append( Poisson )
   res.append( 3.0*elas[1]*(1-2.0*Poisson) )
   return res
#
def ShearBulk(elas):
   Poisson = elas[0]
   Young = elas[1]
   res = [Young/(2.0*(1.0+Poisson)),Young/(3.0*(1.0-2.0*Poisson))]
   return res

# Usefull tensorial operations (second-order symétrical tensors)
# XX, YY, ZZ, YZ, XZ, XY
def zerot2():
   zero = [0.0,0.0,0.0,0.0,0.0,0.0]
   return zero  

def idt2(moy=1):
   ii = [moy,moy,moy,0.0,0.0,0.0]
   return ii

# Spherical part
def epsm(eps):
   trace=0.0
   for i in range(3):
      trace+=eps[i]
   return (trace/3.0)

# Deviatoric part
def deviator(eps):
   epsm_ = epsm(eps)
   dev=zerot2()
   for i in range(0,6):
      if i<=2:
         dev[i] = eps[i]-epsm_
      else:
         dev[i]=eps[i]
   return dev

# Double contracted product
def contract2(eps,sig):
   res = 0.0
   coeff=[1.0,1.0,1.0,2.0,2.0,2.0]
   for j in range(0,6):
      res+=coeff[j]*eps[j]*sig[j] 
   return res

# von Mises
def epseq(eps):
   deps = deviator(eps)
   res=contract2(deps,deps)
   return sqrt(2.0*res/3.0)

def t4contract2(mod,eps):
   res=zerot2()
   eps_dev = deviator(eps)
   epsm_ = epsm(eps)
   for j in range(0,6):
      res[j] = mod[0]*eps_dev[j]
      if j<=2:
         res[j]+=mod[1]*epsm_
   return res

def souplesse(mod):
   souplesse=[1.0/(2.0*mod[0]),1.0/(3.0*mod[1])]
   return souplesse

def testt2():
   eps=[0.0,0.0,0.0,0.0,0.0,sqrt(3)]
   deps = deviator(eps)
   print("testt2, tensorial operations (invariants) = ")
   print("testt2, strain = ",eps)
   print("testt2, deviatoric part = ",deps)
   print("testt2, invariants, epsm = ",epsm(eps),", epseq=",epseq(eps))
   sig=[0.0,0.0,0.0,0.0,0.0,1.0/sqrt(3)]
   print("testt2, sig = ",sig)
   print("testt2, eps:sig = ",contract2(eps,sig))

# Mori-Tanaka model
# (1) => matrix phase
# (2) => inclusions
# Effective moduli
def elaseff(c2,elas1,elas2):
   c1 = 1-c2
   mu1 = elas1[0]
   k1 = elas1[1]
   mu2=elas2[0]
   k2 = elas2[1]
   mustar = (mu1/6)*(9*k1+8*mu1)/(k1+2*mu1)
   kstar = 4*mu1/3.0
   eff=[]
   eff.append( mu1 + c2*( (mu2-mu1)/(1+c1*(mu2-mu1)/(mustar+mu1)) ) ) 
   eff.append( k1 + c2*( (k2-k1)/(1+c1*(k2-k1)/(kstar+k1)) ) )
   return eff

# General expressions true for isotropic two-phases composites
# Whatever the chosen estimate of the effective behaviour
# Effective polarization
def taueff(c2,elas1,elas2,eff,tau1,tau2):
   c1 = 1.0-c2
   eff = elaseff(c2,elas1,elas2)
   locmat = loc_mat(c2,elas1,elas2,eff)
   # Strain localization in inclusions
   locinc=[0.0,0.0]
   for i in range(0,2):
      locinc[i] = (1.0-c1*locmat[i])/c2
   # Spherical part
   tautilde_tr = c1*locmat[1]*epsm(tau1) + c2*locinc[1]*epsm(tau2)
   # All components
   dtau1 = deviator(tau1)
   dtau2 = deviator(tau2)
   res=zerot2()
   for j in range(0,6):
      res[j] = c1*locmat[0]*dtau1[j] + c2*locinc[0]*dtau2[j]
      if j<=2:
         res[j] += tautilde_tr
   return res

# Average per phase of the strain localization tensor, purely mechanical loadings
# loc_mat[0] = deviatoric component
# loc_mat[1] = bulk component
def loc_mat(c2,elas1,elas2,eff):
   c1=1.0-c2
   moy=[0.0,0.0]
   aa=[0.0,0.0]
   # i = phase index
   for i in range(0,2):
      moy[i] = c1*elas1[i] + c2*elas2[i]
   for i in range(0,2):
      aa[i] = 1.0 + (1.0/c1)*(eff[i]-moy[i])/(elas1[i]-elas2[i])
   return aa
 
def loc_inc(c2,elas1,elas2,eff):
   c1=1.0-c2
   aa2=[0.0,0.0]
   locmat = loc_mat(c2,elas1,elas2,eff)
   for i in range(0,2):
      aa2[i] = (1.0-c1*locmat[i])/c2
   return aa2
 
# Average per phase of the strain localization tensor, overal strain = 0 (purely polarization loading) 
def loc_mat_pol(c2,elas1,elas2,eff,tau1,tau2):
   c1=1.0-c2
   moy=[0.0,0.0]
   aa1=zerot2()
   # i = phase index
   for i in range(0,2):
      moy[i] = c1*elas1[i] + c2*elas2[i]
   # Bulk component
   aa1_tr = (1.0/(3.0*c1))*(eff[1]-moy[1])*(epsm(tau1)-epsm(tau2))/((elas1[1]-elas2[1])**2)
   # All components
   dtau1 = deviator(tau1)
   dtau2 = deviator(tau2)
   for j in range(0,6):
      aa1[j] = (1.0/(2.0*c1))*(eff[0]-moy[0])*(dtau1[j]-dtau2[j])/((elas1[0]-elas2[0])**2)
      if j<=2:
         aa1[j]+=aa1_tr
   return aa1

def loc_inc_pol(c2,elas1,elas2,eff,tau1,tau2):
   c1=1.0-c2
   moy=[0.0,0.0]
   aa2=zerot2()
   # i = phase index
   for i in range(0,2):
      moy[i] = c1*elas1[i] + c2*elas2[i]
   # Bulk component
   aa2_tr = (1.0/(3.0*c2))*(eff[1]-moy[1])*(epsm(tau2)-epsm(tau1))/((elas2[1]-elas1[1])**2)
   # All components
   dtau1 = deviator(tau1)
   dtau2 = deviator(tau2)
   for j in range(0,6):
      aa2[j] = (1.0/(2.0*c2))*(eff[0]-moy[0])*(dtau2[j]-dtau1[j])/((elas2[0]-elas1[0])**2)
      if j<=2:
         aa2[j]+=aa2_tr
   return aa2

def weff(c2,elas1,elas2,epsmacro,tau1,tau2):
   c1 = 1-c2
   eff = elaseff(c2,elas1,elas2)
   # purely mechanical part
   res = weff_meca(eff[0],eff[1],epseq(epsmacro),epsm(epsmacro))
   tautilde_ = taueff(c2,elas1,elas2,eff,tau1,tau2)
   # additional terms related to the prescribed polarizations
   aamat = loc_mat_pol(c2,elas1,elas2,eff,tau1,tau2)
   aainc = loc_inc_pol(c2,elas1,elas2,eff,tau1,tau2)
   res = res + contract2(epsmacro,tautilde_)
   res = res + 0.5*(c1*contract2(aamat,tau1) + c2*contract2(aainc,tau2))
   souplesse1 = souplesse(elas1)
   souplesse2 = souplesse(elas2)
   res = res + 0.5*(c1*contract2(tau1,t4contract2(souplesse1,tau1))+c2*contract2(tau2,t4contract2(souplesse2,tau2)))
   return res

