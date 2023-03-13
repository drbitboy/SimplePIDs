import os
import sys

def float_arg(arg,default):
  try: return float(arg)
  except: return float(default)

class PROCESS_MODEL:
  def __init__(self):
    self.pv = self.cv = 0.
    self.tau = 10.
  def setparms(self,tau=None,**kwargs):
    self.tau = float_arg(tau,self.tau)
    return self
  def update(self,cv):
    self.cv = cv
    self.pv += (self.cv - self.pv) / self.tau

class PID:
  def __init__(self,MODEL_CLASS,**kwargs):
    self.model = MODEL_CLASS().setparms(**kwargs)
    self.direct = False
    self.last_pv = self.sp = self.model.pv
    self.cv = self.model.cv
    self.last_error,self.Kp,self.Ti,self.Td = [0.]*4
    self.tune(**kwargs)
  def tune(self,Kp=None,Ti=None,Td=None,direct=None,**kwargs):
    self.Kp = float_arg(Kp,self.Kp)
    self.Ti = float_arg(Ti,self.Ti)
    self.Td = float_arg(Td,self.Td)
    if not (direct is None):
      self.direct = (True is direct) and True or False
    return self
  def setsp(self,sp=None,**kwargs):
    self.sp = float_arg(sp,self.sp)
    return self
  def update(self):
    error = self.model.pv - self.sp
    proportional = (     error - self.last_error) * self.Kp
    integral     = (                       error) / self.Ti
    derivative   = (self.model.pv - self.last_pv) * self.Td
    dcv = proportional + integral + derivative
    self.cv += (dcv * (self.direct and -1. or +1.))
    self.last_error,self.last_pv = error,self.model.pv
    self.model.update(self.cv)
    return self.last_error

def getargs(args,dikt=dict()):
  for arg in args:
    if arg.startswith('--'):
      toks = arg[2:].split('=')
      dikt[toks[0]] = (1==len(toks)) and True or '='.join(toks[1:])
  return dikt

if "__main__" == __name__:
  import matplotlib.pyplot as plt
  ga = getargs(sys.argv[1:],dict(Kp=1,Ti=3,sp=10,direct=True))
  for Ti in [16,8,2,4,1,.5,.25]:
    ga['Ti'] = Ti
    pid = PID(PROCESS_MODEL).tune(**ga).setsp(**ga)
    errors = [1.]
    while max(map(abs,errors[-10:])) > 1e-3: errors.append(pid.update())
    plt.plot(errors[1:],label=f"Ti={Ti}")
  plt.xlabel('Update (~ time)')
  plt.ylabel('Error')
  plt.legend()
  plt.show()
