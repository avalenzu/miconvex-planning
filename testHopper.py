from __future__ import division
import numpy as np
from pyomo.environ import *
from pyomo.opt import SolverFactory

from math import sqrt
from hopper import Hopper
from hopperUtil import *

desiredPrecision = 4
N = 20
tf = 2*1.6
legLength = 0.3
r0 = [0, legLength/2]
rf = [2, legLength]
v0 = [0, 0]
w0 = 0
hipOffset = {'front': {'x': 0.5, 'z': -0.25}, 'hind': {'x': -0.5, 'z': -0.25}}
step_height = 0.1
platform1_start = -0.5
platform1_end = 0.3
platform1_height = 0*step_height
platform2_start = 0.7
platform2_end = 1.3
platform2_height = step_height
platform3_start = 1.7
platform3_end = 9.5
platform3_height = 2*step_height

matlab_hopper = eng.Hopper(legLength, hipOffset)
hop = Hopper(N, eng, matlab_hopper)
# hop.mdt_precision = int(ceil(-np.log2(desiredPrecision)))
hop.dt = tf/N/sqrt(legLength/9.81)
hop.rotationMax = np.pi/8
hop.nOrientationSectors = 1 #int(floor(np.pi/8/desiredPrecision))
print 'hop.nOrientationSectors = %d' % hop.nOrientationSectors
hop.velocityMax = 3
hop.positionMax = 7
hop.forceMax = 3
hop.addPlatform(platform1_start/legLength, platform1_end/legLength, platform1_height/legLength, 1)
hop.addPlatform(platform2_start/legLength, platform2_end/legLength, platform2_height/legLength, 1)
hop.addPlatform(platform3_start/legLength, platform3_end/legLength, platform3_height/legLength, 1)
hop.addFreeBlock(bottom=platform1_height/legLength, right=platform2_start/legLength)
hop.addFreeBlock(bottom=platform2_height/legLength, right=platform3_start/legLength)
hop.addFreeBlock(bottom=platform3_height/legLength)
hop.constructVisualizer()
hop.constructPyomoModel()

m_nlp = hop.model

def objRule(m):
    #     return sum(m.beta[foot, bv, ti]**2 for foot in m.feet for bv in m.BV_INDEX for ti in m.t)
    #     + sum(m.pdd[foot, i, j]**2 for foot in m.feet for i in m.R2_INDEX for j in m.t)
        return sum(m.f[foot, i, j]**2 + m.pdd[foot, i, j]**2 for foot in m.feet for i in m.R2_INDEX for j in m.t) + sum(m.T[ti]**2 for ti in m.t)

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

# m_nlp.maxVerticalVelocityConstraint = Constraint(m_nlp.t, rule=_maxVerticalVelocityRule)

def _periodicFootPosition(m, foot, xz):
    return m.p[foot, xz, m.t[1]] == m.p[foot, xz, m.t[-1]]

#m_nlp.periodicFootPosition = Constraint(m_nlp.feet, m_nlp.R2_INDEX, rule=_periodicFootPosition)

m = constructMDTModel(m_nlp, desiredPrecision)
for z_data in m.z.values():
    z_data._component().branchPriority = 1

opt_nlp = constructMinotaurSolver()

opt = constructGurobiSolver(MIPFocus=1, TimeLimit=40., Threads=11)

hop.constructVisualizer()

#results = opt_nlp.solve(m_nlp, tee=True)
#m_nlp.solutions.load_from(results)
#hop.loadResults(m_nlp)

#results = opt.solve(m, tee=True)
#m.solutions.store_to(results)
#m_nlp.solutions.load_from(results, ignore_invalid_labels=True)
#fixIntegerVariables(m_nlp)

#results_nlp = opt_nlp.solve(m_nlp, tee=True)

#hop.loadResults(m_nlp)
