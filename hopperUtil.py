from __future__ import division
import numpy as np
from pyomo.opt import SolverFactory
from pyomo.core.plugins.transform.radix_linearization import *

def constructMDTModel(m_nlp, desiredPrecision):
    RadixLinearization.radix = 2
    mdt = RadixLinearization()
    #discretizationVar = m_nlp.footRelativeToCOM
    discretizationVar = m_nlp.f
    fBounds = m_nlp.f.values()[0].bounds
    print fBounds
    footBounds = m_nlp.footRelativeToCOM.values()[0].bounds
    print footBounds
    maxVal = (fBounds[1] - fBounds[0])*(footBounds[1] - footBounds[0])
    precision = int(np.ceil(-np.log2(desiredPrecision/maxVal)))
    print 'MDT precision: %d' % precision
    m = mdt.create_using(m_nlp, verbose=True, precision=precision, discretize=[discretizationVar])

    for z_data in m.z.values():
        z_data._component().branchPriority = 1
    return m

def constructCouenneSolver(**kwargs):
    opt = SolverFactory('couenne')
    return opt

def constructMinotaurSolver():
    opt = SolverFactory('qpd')
    opt.set_options('--AMPL=1')
    opt.set_options('--nlp_engine=IPOPT')
    opt.set_options('--bnb_time_lmit=3600.')
    # opt.set_options('--linfpump=1')
    return opt

def constructGurobiSolver(**kwargs):
    opt = SolverFactory('_gurobi_direct')
    for key, value in kwargs.iteritems():
        opt.set_options('%s=%f' % (key, value))
    return opt
    # opt.set_options('mipgap=0.05')
    #if timeout > 0:
        #opt.set_options('TimeLimit=%f' % timeout)
    #opt.set_options('Threads=%f' % threads)
    # opt.set_options('Seed=0')
    #opt.set_options('Presolve=2')

def extractPostition(m):
    return np.vstack([np.array([m.r[xz, ti].value for ti in m.t]) for xz in m.R2_INDEX])

def extractOrientation(m):
    return np.atleast_2d(np.array([m.th[ti].value for ti in m.t]))

def extractHipPosition(m):
    return np.dstack([np.vstack([np.array([m.hip[foot, xz, ti].value for ti in m.t]) for xz in m.R2_INDEX]) for foot in m.feet])

def extractRelativeFootPosition(m):
    return np.dstack([np.vstack([np.array([m.p[foot, xz, ti].value for ti in m.t]) for xz in m.R2_INDEX]) for foot in m.feet])

def extractFootForce(m):
    return np.dstack([np.vstack([np.array([m.f[foot, xz, ti].value for ti in m.t]) for xz in m.R2_INDEX]) for foot in m.feet])

def extractTotalTorque(m):
    return np.atleast_2d(np.array([m.T[ti].value for ti in m.t]))

def fixIntegerVariables(m):
    for var in m.component_data_objects(Var):
        if not var.is_continuous():
            print 'Fixing %s' % var.cname()
            var.fixed = True
