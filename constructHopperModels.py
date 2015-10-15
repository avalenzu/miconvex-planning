from __future__ import division
import numpy as np
from uuid import uuid4
from pyomo.environ import *
from pyomo.opt import SolverFactory

from math import sqrt
from hopper import Hopper
from hopperUtil import *

desiredPrecision = 2
N = 25
tf = 2*1.6
legLength = 0.16
r0 = [0, legLength/2]
rf = [1.0, legLength]
v0 = [0, 0]
w0 = 0
hipOffset = {'front': {'x': 0.5, 'z': -0.25}, 'hind': {'x': -0.5, 'z': -0.25}}

matlab_hopper = eng.Hopper(legLength, hipOffset)
hop = Hopper(N, eng, matlab_hopper)
# hop.mdt_precision = int(ceil(-np.log2(desiredPrecision)))
hop.dtBounds = tuple((1/sqrt(legLength/9.81))*np.array([0.05, 0.2]))
hop.dtNom = 0.05*(1/sqrt(legLength/9.81))
hop.rotationMax = np.pi/8
hop.nOrientationSectors = 1 #int(floor(np.pi/8/desiredPrecision))
print 'hop.nOrientationSectors = %d' % hop.nOrientationSectors
hop.velocityMax = 3.
hop.positionMax = 1.5*rf[0]/legLength
hop.forceMax = 3.
addThreePlatfomWorld(hop, legLength, 0.25*legLength)
#addFlatWorld(hop, legLength)
hop.constructVisualizer()
m_nlp = hop.constructPyomoModel()
def normL2(m, var):
    index = var.index_set()
    return sum(var[i]**2 for i in index)

def normL1(m, var):
    index = var.index_set()
    slackName = '%sSlacks' % var.cname()
    lbName = '%sLB' % slackName
    ubName = '%sUB' % slackName
    slackMax = max([max(np.abs(var[i].bounds)) for i in index])
    setattr(m, slackName, Var(index, bounds=(0.0, slackMax)))

    def _lbRule(m, *args):
        return var[args] <= getattr(m, slackName)[args]
    setattr(m, lbName, Constraint(index, rule=_lbRule))

    def _ubRule(m, *args):
        return var[args] >= -getattr(m, slackName)[args]
    setattr(m, ubName, Constraint(index, rule=_ubRule))

    return summation(getattr(m, slackName))

def normLInfinity(m, var):
    index = var.index_set()
    slackName = '%sSlack' % var.cname()
    lbName = '%sLB' % slackName
    ubName = '%sUB' % slackName
    slackMax = max([max(np.abs(var[i].bounds)) for i in index])
    setattr(m, slackName, Var(bounds=(0.0, slackMax)))

    def _lbRule(m, *args):
        return var[args] <= getattr(m, slackName)
    setattr(m, lbName, Constraint(index, rule=_lbRule))

    def _ubRule(m, *args):
        return var[args] >= -getattr(m, slackName)
    setattr(m, ubName, Constraint(index, rule=_ubRule))

    return getattr(m, slackName)

def exprNormLInfinity(m, expr, slackMax):
    slackName = 'slack_%s' % str(uuid4()).replace('-','')
    lbName = '%sLB' % slackName
    ubName = '%sUB' % slackName
    setattr(m, slackName, Var(bounds=(0.0, slackMax)))

    def _lbRule(m):
        return expr <= getattr(m, slackName)
    setattr(m, lbName, Constraint(rule=_lbRule))

    def _ubRule(m):
        return expr >= -getattr(m, slackName)
    setattr(m, ubName, Constraint(rule=_ubRule))

    return getattr(m, slackName)


#norm = normL1;
norm = normL2;
#norm = normLInfinity;

def objRule(m):
    #     return sum(m.beta[foot, bv, ti]**2 for foot in m.feet for bv in m.BV_INDEX for ti in m.t)
    #     + sum(m.pdd[foot, i, j]**2 for foot in m.feet for i in m.R2_INDEX for j in m.t)
    #return sum(m.f[foot, i, j]**2 + m.pd[foot, i, j]**2 + m.pdd[foot, i, j]**2 for foot in m.feet for i in m.R2_INDEX for j in m.t) + sum(m.T[ti]**2 for ti in m.t)
    #return sum(m.f[foot, i, j]**2 for foot in m.feet for i in m.R2_INDEX for j in m.t) + sum(m.hipTorque[foot, ti]**2 for foot in m.feet for ti in m.t)
    #return sum(m.f[foot, i, j]**2 + m.pd[foot, i, j]**2  + m.pdd[foot, i, j]**2 for foot in m.feet for i in m.R2_INDEX for j in m.t) + sum(m.hipTorque[foot, ti]**2 for foot in m.feet for ti in m.t) + summation(m.dt)

    #return sum(m.f[foot, i, j]**2 + m.pdd[foot, i, j]**2 for foot in m.feet for i in m.R2_INDEX for j in m.t) + sum(m.hipTorque[foot, ti]**2 for foot in m.feet for ti in m.t) + summation(m.dt)
    footRegionChanges = 0.0
    for t in m.t:
        if t != m.t[-1]:
            for region in m.REGION_INDEX:
                if hop.regions[region]['mu'] != 0.:
                    for foot in m.feet:
                        current_indicator = getattr(m, 'footRegionConstraints[%d,%s,%d]indicator_var' % (region, foot, t));
                        next_indicator = getattr(m, 'footRegionConstraints[%d,%s,%d]indicator_var' % (region, foot, t+1));
                        footRegionChanges += exprNormLInfinity(m, next_indicator - current_indicator, 1.0)
    #return footRegionChanges + norm(m, m.pd) + norm(m, m.f) + norm(m, m.hipTorque)
    return 1e1*footRegionChanges + norm(m, m.pdd) + norm(m, m.beta) + norm(m, m.hipTorque)

m_nlp.Obj = Objective(rule=objRule, sense=minimize)

m_nlp.rx0 = Constraint(expr=m_nlp.r['x',m_nlp.t[1]] == r0[0]/legLength)
#m_nlp.rz0 = Constraint(expr=m_nlp.r['z',m_nlp.t[1]] <= 1)

m_nlp.th0 = Constraint(expr=m_nlp.th[m_nlp.t[1]] == 0)

m_nlp.vx0 = Constraint(expr=m_nlp.v['x',m_nlp.t[1]] == 0)
m_nlp.vz0 = Constraint(expr=m_nlp.v['z',m_nlp.t[1]] == 0)

m_nlp.w0 = Constraint(expr=m_nlp.w[m_nlp.t[1]] == 0)

m_nlp.Fx0 = Constraint(expr=m_nlp.F['x', m_nlp.t[1]] == 0)
m_nlp.Fz0 = Constraint(expr=m_nlp.F['z', m_nlp.t[1]] == 0)
m_nlp.T0 = Constraint(expr=m_nlp.T[m_nlp.t[1]] == 0)

m_nlp.rxf = Constraint(expr=m_nlp.r['x',m_nlp.t[-1]] >= rf[0]/legLength)
#m_nlp.rzf = Constraint(expr=m_nlp.r['z',m_nlp.t[-1]] == m_nlp.r['z', m_nlp.t[1]])

m_nlp.thf = Constraint(expr=m_nlp.th[m_nlp.t[-1]] == 0)

m_nlp.vxf = Constraint(expr=m_nlp.v['x',m_nlp.t[-1]] == m_nlp.v['x',m_nlp.t[1]])
m_nlp.vzf = Constraint(expr=m_nlp.v['z',m_nlp.t[-1]] == 0)

m_nlp.wf = Constraint(expr=m_nlp.w[m_nlp.t[-1]] == 0)

m_nlp.Fxf = Constraint(expr=m_nlp.F['x', m_nlp.t[-1]] == 0)
m_nlp.Fzf = Constraint(expr=m_nlp.F['z', m_nlp.t[-1]] == 0)
m_nlp.Tf = Constraint(expr=m_nlp.T[m_nlp.t[-1]] == 0)

def _maxVerticalVelocityRule(m, t):
    return m.v['z', t] <= 0.5

m_nlp.maxVerticalVelocityConstraint = Constraint(m_nlp.t, rule=_maxVerticalVelocityRule)

def _periodicFootPosition(m, foot, xz):
    return m.p[foot, xz, m.t[1]] == m.p[foot, xz, m.t[-1]]

#m_nlp.periodicFootPosition = Constraint(m_nlp.feet, m_nlp.R2_INDEX, rule=_periodicFootPosition)

#m = constructMDTModel(m_nlp, desiredPrecision)
m = constructRelaxedModel(m_nlp)
#for z_data in m.z.values():
    #z_data._component().branchPriority = 1
m_nlp_orig = m_nlp.clone()
#m.dt.fix()

#def _momentRule(m, t):
    #return m.T[t] == sum(m.footRelativeToCOM[foot,'x',t]*m.f[foot,'z',t] - m.footRelativeToCOM[foot,'z',t]*m.f[foot, 'x',t] for foot in m.feet)

#m_nlp.momentAbountCOM = Constraint(m_nlp.t, rule=_momentRule)

def _hipTorqueRule(m, foot, t):
    return m.hipTorque[foot, t] == m.p[foot,'x',t]*m.f[foot,'z',t] - m.p[foot,'z',t]*m.f[foot, 'x',t]

#m_nlp.hipTorqueConstraint = Constraint(m_nlp.feet, m_nlp.t, rule=_hipTorqueRule)


m_nlp.pwSin.deactivate()
m_nlp.pwCos.deactivate()

def _cos(m, t):
    return m.cth[t] == cos(m.th[t])
m_nlp.Cos = Constraint(m_nlp.t, rule=_cos)

def _sin(m, t):
    return m.sth[t] == sin(m.th[t])
m_nlp.Sin = Constraint(m_nlp.t, rule=_sin)

opt_nlp = SolverFactory('ipopt')
opt_minlp = constructCouenneSolver()

#opt = constructGurobiSolver(mipgap=0.8, MIPFocus=1, TimeLimit=90., Threads=11)
opt = constructGurobiSolver(TimeLimit=480., Threads=11)
#opt = constructGurobiSolver(TimeLimit=50., Threads=11)

hop.constructVisualizer()
