import os
import sys
import math

def float_arg(arg,default):
  try: return float(arg)
  except: return float(default)
def bool_arg(arg,default):
  c1 = str(arg).strip()[:1].lower()
  if c1 in 'ty1': return True
  if c1 in 'fn0' and (not (None is arg)): return False
  return default and True or False

class SHELL_BALL_BEAM_MODEL:
  def __init__(self):
    self.velocity = self.position = self.theta = 0.
    self.g = -9.8
    self.channel_angle_deg = 90.0
    self.rate_limit = 0.0
    self.use_sine = True

  def setparms(self
              ,angle_deg=None
              ,rate_limit=None
              ,use_sine=None
              ,position=None
              ,velocity=None
              ,theta=None
              ,**kwargs
              ):
    self.use_sine = bool_arg(use_sine,self.use_sine)
    self.rate_limit = float_arg(rate_limit,self.rate_limit)
    self.channel_angle_deg = float_arg(angle_deg,self.channel_angle_deg)
    self.position = float_arg(position,self.position)
    self.velocity = float_arg(velocity,self.velocity)
    self.theta = float_arg(theta,self.theta)
    assert self.channel_angle_deg >= 1.0
    assert self.channel_angle_deg <= 180.0
    sine_half = math.sin(self.channel_angle_deg * math.pi / 360.0)
    k = (2.0 / ( 3.0 * sine_half * sine_half))
    self.accel_coeff = self.g / (k + 1)
    return self

  def CV2theta(self,CV):
    """Convert from PID output CV to target theta, radians"""
    ### No scaling for now:  assume CV is theta in radians
    return CV

  def getCV(self): return self.theta
  def getPV(self): return self.position

  def update(self,CV,ts):
    """Integrate over timestep (ts) assuming constant dtheta/dt"""

    ### Save theta at beginning of, calculate target theta at end of, ts
    old_theta = self.theta
    target_theta = self.CV2theta(CV)

    ### Impose rate limit, if present, on change in theta

    if self.rate_limit <= 0.0:
      ### No rate limit, assume instantaneous response
      self.theta = target_theta
      tknee = 0.0

    else:
      ### Check change in theta over timestep against rate limit:

      ### - First calculate targeted change
      dtheta = target_theta - old_theta

      if abs(dtheta) > 1e-9:
        ### - If change is not small, then calculate how long that
        ###   change would take (tknee)
        tknee = abs(dtheta / self.rate_limit)

        if tknee > ts:
          ### - If that takes longer than the timestep (ts), then scale
          ###   back to how much change is possible over the timestep
          dtheta *= (ts / tknee)
          self.theta = old_theta + dtheta
          tknee = ts

        else:
          ### - Target will be reached by time tknee
          self.theta = target_theta

        ### Calculate rate of change of theta, positive or negative
        dthetadt = dtheta < 0.0 and -self.rate_limit or self.rate_limit

        cosT0,cosT1 = map(math.cos,(old_theta,self.theta,))
        sinT0,sinT1 = map(math.sin,(old_theta,self.theta,))
        gok = self.accel_coeff / dthetadt

        ### Integrate velocity from 0.0 to tknee (tknee <= ts)
        old_velocity = self.velocity
        self.velocity -= gok * (cosT1-cosT0)

        ### Integrate position from 0.0 to tknee
        gok2 = gok / dthetadt
        self.position += tknee * ( old_velocity + (gok * cosT0))
        self.position -= gok2 * (sinT1-sinT0)

      else:
        ### - Small change:  assume constant theta over entire timestep
        tknee = 0.0
        self.theta = target_theta

    ### Integrate velocity & position @ constant theta, from tknee to ts
    if tknee < ts:
      sinT1 = math.sin(self.theta)
      dt = ts - tknee
      delta_velocity = self.accel_coeff * sinT1 * dt
      self.position += (self.velocity + (delta_velocity / 2)) * dt
      self.velocity += delta_velocity

class FIRST_ORDER_MODEL:
  def __init__(self):
    self.state = self.CV = 0.
    self.tau = 16.
  def setparms(self,tau=None,**kwargs):
    self.tau = float_arg(tau,self.tau)
    return self
  def getCV(self): return self.CV
  def getPV(self): return self.state
  def update(self,CV,ts):
    self.CV = CV
    self.state += (self.CV - self.state) / self.tau

class INDEPENDENTPID:
  def __init__(self,MODEL_CLASS,**kwargs):
    self.model = MODEL_CLASS().setparms(**kwargs)
    self.direct = False
    self.last_pv = self.sp = self.model.getPV()
    self.CV = self.model.getCV()
    self.last_error,self.Kp,self.Ki,self.Kd,self.last_derivterm = [0.]*5
  def tune(self,Kp=None,Ki=None,Kd=None,direct=None,**kwargs):
    self.Kp = float_arg(Kp,self.Kp)
    self.Ki = float_arg(Ki,self.Ki)
    self.Kd = float_arg(Kd,self.Kd)
    if not (direct is None):
      self.direct = (True is direct) and True or False
    return self
  def setsp(self,sp=None,**kwargs):
    self.sp = float_arg(sp,self.sp)
    return self
  def update(self,ts):
    model_pv = self.model.getPV()
    error = model_pv - self.sp
    proportional = (error - self.last_error) * self.Kp
    integral     = (                  error) * self.Ki * ts
    derivative   = (model_pv - self.last_pv) * self.Kd / ts
    dCV = proportional + integral + derivative - self.last_derivterm
    self.CV += (dCV * (self.direct and -1. or +1.))
    self.last_error = error
    self.last_pv = model_pv
    self.last_derivterm = derivative
    self.model.update(self.CV,ts)
    return self.last_error

def getargs(args,dikt=dict()):
  for arg in args:
    if arg.startswith('--'):
      toks = arg[2:].split('=')
      toks0 = toks.pop(0)
      dikt[toks0] = val = (not len(toks)) and True or '='.join(toks)
      dikt[toks0.replace('-','_')] = dikt[toks0.replace('_','-')] = val
  return dikt

if "__main__" == __name__:
  import matplotlib.pyplot as plt

  if getargs(sys.argv[1:]).get('first-order',False):
    ga = getargs(sys.argv[1:],dict(Kp=1,Ki=10,sp=10,direct=True))
    for Ki in [32, 16,4,1,.25,.0625]:
      ga['Ki'] = 1.0 / Ki
      model = FIRST_ORDER_MODEL
      pid = INDEPENDENTPID(model,**ga).tune(**ga).setsp(**ga)
      errors = [1.]
      while max(map(abs,errors[-10:])) > 1e-3:
        errors.append(pid.update(1.0))
      plt.plot(errors[1:],label=f"Ki={Ki}")
    plt.xlabel('Update (~ time)')
    plt.ylabel('Error')

  if getargs(sys.argv[1:]).get('shell-ball-beam',False):
    ga = getargs(sys.argv[1:],dict(Kp=3,Kd=2,sp=0,direct=False))
    model = SHELL_BALL_BEAM_MODEL
    pid = INDEPENDENTPID(model,**ga).tune(**ga).setsp(**ga)
    p0,v0, = pid.model.getPV(),pid.model.velocity
    errors = [1.]
    pvs,cvs,times, = list(),list(),list()
    timestep = ga.get('timestep',0.01)
    while max(map(abs,errors[-10:])) > 1e-4:
      errors.append(pid.update(timestep))
      cvs.append(pid.CV)
      pvs.append(pid.last_pv)
      times.append(times and (times[-1]+timestep) or 0.0)
      if len(errors) > 10000: break

    """
    plt.ylabel('PV (Position, m)\nCV (Tilt angle, rad)')
    plt.axhline(0.0,linestyle='dotted',color='black')
    ##plt.plot(times,errors[1:],label='Error')
    plt.plot(times,pvs,label='Position')
    plt.plot(times,cvs,label='Tilt angle')
    """

    fig,ax1=plt.subplots()

    c1 = 'tab:red'
    c2 = 'tab:blue'

    ax1.set_xlabel('Time, s')
    ax1.set_ylabel('CV (Tilt angle, rad)',color=c1)
    ax1.axhline(0,linestyle='dotted',color=c1)
    ax1.tick_params(axis='y',labelcolor=c1)
    ax1.set_ylim([-1.0, 0.5])

    ax2=ax1.twinx()

    ax2.set_ylabel('PV (Position, m)',color=c2)
    ax2.axhline(0,linestyle='dotted',color=c2)
    ax2.tick_params(axis='y',labelcolor=c2)
    ax2.set_ylim([-.151, 0.01])

    ang, = ax1.plot(times,cvs,color=c1)
    pos, = ax2.plot(times,pvs,color=c2)

    plt.title(f"""
Model[Accel={pid.model.accel_coeff}(m/s$^{{2}}$)/rad, Pos$_{{Initial}}$={p0}m, Vel$_{{Initial}}$={v0}m/s]
PID[K$_{{P}}$={pid.Kp}, K$_{{I}}$={pid.Ki}, K$_{{D}}$={pid.Kd}]
""".strip())

  try:
    assert len(errors) > 1
    plt.legend(handles=[pos,ang],loc='lower right')
    plt.show()
  except:
    pass
