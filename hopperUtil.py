from __future__ import division
import numpy as np
from pyomo.opt import SolverFactory
from pyomo.core.plugins.transform.radix_linearization import *
from mccormick_envelope import *
from pyomo.core.base.component import register_component, Component, ComponentUID


def constructRelaxedModel(m_nlp, dt=None):
    m = m_nlp.clone()
    if dt is not None:
        m.dt.fix(dt)
    else:
        m.dt.fix()
    mccormick = McCormickEnvelope()
    m = mccormick.create_using(m, verbose=True)
    return m

def constructMDTModel(m_nlp, desiredPrecision, dt=None):
    m = m_nlp.clone()
    if dt is not None:
        m.dt.fix(dt)
    else:
        m.dt.fix()
    RadixLinearization.radix = 2
    mdt = RadixLinearization()
    #discretizationVar = m.footRelativeToCOM
    discretizationVar = m.f
    fBounds = m.f.values()[0].bounds
    print fBounds
    footBounds = m.footRelativeToCOM.values()[0].bounds
    print footBounds
    maxVal = (fBounds[1] - fBounds[0])*(footBounds[1] - footBounds[0])
    precision = 1+0*int(np.ceil(-np.log2(desiredPrecision/maxVal)))
    print 'MDT precision: %d' % precision
    m = mdt.create_using(m, verbose=True, precision=desiredPrecision, discretize=[discretizationVar])

    for z_data in m.z.values():
        z_data._component().branchPriority = 1
    return m

def constructCouenneSolver(**kwargs):
    opt = SolverFactory('couenne')
    return opt

def constructMinotaurSolver():
    opt = SolverFactory('qpd')
    opt.set_options('--ampl=1')
    opt.set_options('--nlp_engine=IPOPT')
    opt.set_options('--bnb_time_limit=3600.')
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

def fixIntegerVariables(m):
    for var in m.component_data_objects(Var):
        if not var.is_continuous():
            #print 'Fixing %s to %s' % (ComponentUID(var), var.value)
            var.fixed = True

def unfixIntegerVariables(m):
    for var in m.component_data_objects(Var):
        if not var.is_continuous():
            #print 'Fixing %s to %s' % (ComponentUID(var), var.value)
            var.fixed = False

def addThreePlatfomWorld(hop, legLength, step_height):
    step_length = 2.0*legLength
    gap_length = 0.65*step_length
    platform1_start = -1*legLength
    platform1_end = platform1_start + step_length
    platform1_height = 0*step_height
    platform2_start = platform1_end + gap_length
    platform2_end = platform2_start + step_length
    platform2_height = 2*step_height
    platform3_start = platform2_end + gap_length
    platform3_end = platform3_start + step_length
    platform3_height = -1*step_height
    hop.addPlatform(platform1_start/legLength, platform1_end/legLength, platform1_height/legLength, 1, 0.5*4.78*step_length, -0.5*4.78*step_length)
    #hop.addPlatform(platform2_start/legLength, platform2_end/legLength, platform2_height/legLength, 1, 0.5*4.78*step_length, 0.5*4.78*step_length - 3.25*step_length)
    hop.addPlatform(platform2_start/legLength, platform2_end/legLength, platform2_height/legLength, 1, 0.5*4.78*step_length, 0.5*4.78*step_length - 2.25*step_length)
    hop.addPlatform(platform3_start/legLength, platform3_end/legLength, platform3_height/legLength, 1, 0.5*4.78*step_length, -0.5*4.78*step_length)
    hop.addFreeBlock(bottom=platform1_height/legLength, right=platform2_start/legLength)
    hop.addFreeBlock(bottom=platform2_height/legLength, left=platform1_end/legLength, right=platform3_start/legLength)
    hop.addFreeBlock(bottom=platform3_height/legLength, left=platform2_end/legLength)

def addFlatWorld(hop, legLength):
    hop.addPlatform(0., 1.8/legLength, 0., 1, 0.6/legLength, 0.)
    hop.addFreeBlock(bottom=0., left=0)


