# These are the Kosugi (200x) equations
# The input is matric potential, psi and the hydraulic parameters.
# psi must be sent in as a numpy array.
# The pars variable is like a MATLAB structure.

# Need to also pack these functions into pars
#	erf: from scipy.special import erf
#   log: from pylab import log

def thetaFun(psi,pars):
  Se=0.5+0.5*pars.erf(-(pars.log(psi/pars.psi0)/pars.s-pars.s)/2**0.5)
  Se[psi>=0]=1.
  return pars.thetaR+(pars.thetaS-pars.thetaR)*Se
  
def CFun(psi,pars):
  Se=0.5+0.5*pars.erf(-(pars.log(psi/pars.psi0)/pars.s-pars.s)/2**0.5)
  Se[psi>=0]=1.
  dSedh=1/(2*pars.pi)**0.5/pars.s/(-psi)*pars.exp(-[
  # dSedh(psi>=0)=0;
  return Se*pars.Ss+(pars.thetaS-pars.thetaR)*dSedh
  
def KFun(psi,pars):
  Se=(1+abs(psi*pars.alpha)**pars.n)**(-pars.m)
  Se[psi>=0]=1.
  return pars.Ks*Se**pars.neta*(1-(1-Se**(1/pars.m))**pars.m)**2
  
def setpars():
  from scipy.special import erf
  from pylab import log
  from pylab import pi
  from pylab import exp
  class Args: pass
  pars=Args()
  pars.erf=erf
  pars.log=log
  pars.exp=exp
  pars.pi=pi
  pars.thetaR=float(raw_input("thetaR = "))
  pars.thetaS=float(raw_input("thetaS = "))
  alpha=float(raw_input("alpha = "))
  n=float(raw_input("n = "))
  m=1-1/n
  pars.s=((1-m)*log((2**(1/m)-1)/m))**0.5
  pars.psi0=-m**(1-m)/alpha
  pars.Ks=float(raw_input("Ks = "))
  pars.neta=float(raw_input("neta = "))
  pars.Ss=float(raw_input("Ss = "))
  
  return pars
  
def PlotProps(pars):
  import numpy as np
  import pylab as pl
  import kosugi as ks
  psi=np.linspace(-10,2,200)
  pl.figure
  pl.subplot(3,1,1)
  pl.plot(psi,vg.thetaFun(psi,pars))
  pl.ylabel(r'$\theta(\psi) [-]$')
  pl.subplot(3,1,2)
  pl.plot(psi,vg.CFun(psi,pars))
  pl.ylabel(r'$C(\psi) [1/m]$')
  pl.subplot(3,1,3)
  pl.plot(psi,vg.KFun(psi,pars))
  pl.xlabel(r'$\psi [m]$')
  pl.ylabel(r'$K(\psi) [m/d]$')
  #pl.show()
  
def HygieneSandstone():
  from scipy.special import erf
  from pylab import log
  from pylab import pi
  from pylab import exp
  class Args: pass
  pars=Args()
  pars.erf=erf
  pars.log=log
  pars.exp=exp
  pars.pi=pi
  pars.thetaR=0.153
  pars.thetaS=0.25
  alpha=0.79
  n=10.4
  m=1-1/n
  pars.s=((1-m)*log((2**(1/m)-1)/m))**0.5
  pars.psi0=-m**(1-m)/alpha
  pars.Ks=1.08
  pars.neta=0.5
  pars.Ss=0.000001
  return pars
  
def TouchetSiltLoam():
  from scipy.special import erf
  from pylab import log
  from pylab import pi
  from pylab import exp
  class Args: pass
  pars=Args()
  pars.erf=erf
  pars.log=log
  pars.exp=exp
  pars.pi=pi
  pars.thetaR=0.19
  pars.thetaS=0.469
  alpha=0.5
  n=7.09
  m=1-1/n
  pars.s=((1-m)*log((2**(1/m)-1)/m))**0.5
  pars.psi0=-m**(1-m)/alpha
  pars.Ks=3.03
  pars.neta=0.5
  pars.Ss=0.000001
  return pars

def SiltLoamGE3():
  from scipy.special import erf
  from pylab import log
  from pylab import pi
  from pylab import exp
  class Args: pass
  pars=Args()
  pars.erf=erf
  pars.log=log
  pars.exp=exp
  pars.pi=pi
  pars.thetaR=0.131
  pars.thetaS=0.396
  alpha=0.423
  n=2.06
  m=1-1/n
  pars.s=((1-m)*log((2**(1/m)-1)/m))**0.5
  pars.psi0=-m**(1-m)/alpha
  pars.Ks=0.0496
  pars.neta=0.5
  pars.Ss=0.000001
  return pars
  
def GuelphLoamDrying():
  from scipy.special import erf
  from pylab import log
  from pylab import pi
  from pylab import exp
  class Args: pass
  pars=Args()
  pars.erf=erf
  pars.log=log
  pars.exp=exp
  pars.pi=pi
  pars.thetaR=0.218
  pars.thetaS=0.520
  alpha=1.15
  n=2.03
  m=1-1/n
  pars.s=((1-m)*log((2**(1/m)-1)/m))**0.5
  pars.psi0=-m**(1-m)/alpha
  pars.Ks=0.316
  pars.neta=0.5
  pars.Ss=0.000001
  return pars
  
def GuelphLoamWetting():
  from scipy.special import erf
  from pylab import log
  from pylab import pi
  from pylab import exp
  class Args: pass
  pars=Args()
  pars.erf=erf
  pars.log=log
  pars.exp=exp
  pars.pi=pi
  pars.thetaR=0.218
  pars.thetaS=0.434
  alpha=2.0
  n=2.76
  m=1-1/n
  pars.s=((1-m)*log((2**(1/m)-1)/m))**0.5
  pars.psi0=-m**(1-m)/alpha
  pars.Ks=0.316
  pars.neta=0.5
  pars.Ss=0.000001
  return pars
  
def BeitNetofaClay():
  from scipy.special import erf
  from pylab import log
  from pylab import pi
  from pylab import exp
  class Args: pass
  pars=Args()
  pars.erf=erf
  pars.log=log
  pars.exp=exp
  pars.pi=pi
  pars.thetaR=0.
  pars.thetaS=0.446
  alpha=0.152
  n=1.17
  m=1-1/n
  pars.s=((1-m)*log((2**(1/m)-1)/m))**0.5
  pars.psi0=-m**(1-m)/alpha
  pars.Ks=0.00082
  pars.neta=0.5
  pars.Ss=0.000001
  return pars
  
