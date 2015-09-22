#  _________________________________________________________________________
#
# Derived from:
#  Pyomo: Python Optimization Modeling Objects
#  Copyright (c) 2014 Sandia Corporation.
#  Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
#  the U.S. Government retains certain rights in this software.
#  This software is distributed under the BSD License.
#  _________________________________________________________________________

from pyomo.util.plugin import alias
from pyomo.core import Binary, value, as_numeric
from pyomo.core.base import Transformation, Var, Constraint, ConstraintList, Block, RangeSet
from pyomo.core.base.expr import _ProductExpression, _PowExpression
from pyomo.core.base.var import _VarData

from scipy.constants import golden

from six import iteritems

import logging
logger = logging.getLogger(__name__)


class McCormickEnvelope(Transformation):
    """
    This plugin generates linear relaxations of bilinear problems using McCormick envelopes.
    """

    alias("core.mccormick_envelope",
           doc="Linearize bilinear and quadratic terms through "
           "McCormick envelopes" )
    radix = 2
    def _create_using(self, model, **kwds):
        verbose = kwds.pop('verbose',False)

        M = model.clone()

        # Iterate over all Constraints and identify the bilinear and
        # quadratic terms
        bilinear_terms = []
        quadratic_terms = []
        for constraint in M.component_map(Constraint, active=True).itervalues():
            for cname, c in constraint._data.iteritems():
                if c.body.polynomial_degree() != 2:
                    continue
                self._collect_bilinear(c.body, bilinear_terms, quadratic_terms)


        #
        # Relax things
        #

        # Define a block (namespace) for holding the disaggregated
        # variables and new constraints
        if False: # Set to true when the LP writer is fixed
            M._radix_linearization = Block()
            _block = M._radix_linearization
        else:
            _block = M

        _known_bilinear = {}
        # For each quadratic term, if it hasn't been discretized /
        # generated, do so, and remember the resulting W term for later
        # use...
        #for _expr, _x1 in quadratic_terms:
        #    self._discretize_term( _expr, _x1, _x1,
        #                           _block, _discretize, _known_bilinear )
        # For each bilinear term, if it hasn't been discretized /
        # generated, do so, and remember the resulting W term for later
        # use...
        for _expr, _x1, _x2 in bilinear_terms:
            self._relax_term(_expr, _x1, _x2, _block, _known_bilinear)
        for _expr, _x1 in quadratic_terms:
            self._relax_term(_expr, _x1, _x1, _block, _known_bilinear)

        # Return the relaxed instance!
        return M

    def _relax_term(self, _expr, _u, _v, _block, _known_bilinear):
        _id = (id(_v), id(_u))
        if _id not in _known_bilinear:
            _known_bilinear[_id] = self._relax_bilinear(_block, _u, _v)
        # _expr should be a "simple" product expression; substitute
        # in the bilinear "W" term for the raw bilinear terms
        _expr._numerator = [ _known_bilinear[_id] ]


    def _relax_bilinear(self, b, u, v):
        u_lb, u_ub = u.bounds
        v_lb, v_ub = v.bounds
        if u_lb is None or u_ub is None:
             raise RuntimeError("Couldn't relax variable %s: missing "
                               "finite lower/upper bounds." % (u.cname(True)))
        if v_lb is None or v_ub is None:
             raise RuntimeError("Couldn't relax variable %s: missing "
                               "finite lower/upper bounds." % (v.cname(True)))
        w = Var(bounds=(min(u_lb*v_lb, u_lb*v_ub, u_ub*v_lb, u_ub*v_ub),
                         max(u_lb*v_lb, u_lb*v_ub, u_ub*v_lb, u_ub*v_ub)))
        b.add_component("w_%s_%s" % (u.cname(), v.cname()), w)

        _c = ConstraintList(noruleinit=True)
        b.add_component( "c_mccormick_%s_%s" % (u.cname(), v.cname()), _c )

        _c.add(expr=w >= u * v_lb + u_lb * v - u_lb*v_lb)
        _c.add(expr=w >= u * v_ub + u_ub * v - u_ub*v_ub)
        _c.add(expr=w <= u * v_lb + u_ub * v - u_ub*v_lb)
        _c.add(expr=w <= u * v_ub + u_lb * v - u_lb*v_ub)

        return w

    def _collect_bilinear(self, expr, bilin, quad):
        if not expr.is_expression():
            return
        if type(expr) is _ProductExpression:
            if len(expr._numerator) != 2:
                for e in expr._numerator:
                    self._collect_bilinear(e, bilin, quad)
                # No need to check denominator, as this is poly_degree==2
                return
            if not isinstance(expr._numerator[0], _VarData) or \
                    not isinstance(expr._numerator[1], _VarData):
                raise RuntimeError("Cannot yet handle complex subexpressions")
            if expr._numerator[0] is expr._numerator[1]:
                quad.append( (expr, expr._numerator[0]) )
            else:
                bilin.append( (expr, expr._numerator[0], expr._numerator[1]) )
            return
        if type(expr) is _PowExpression and value(expr._args[1]) == 2:
            # Note: directly testing the value of the exponent above is
            # safe: we have already verified that this expression is
            # polynominal, so the exponent must be constant.
            tmp = _ProductExpression()
            tmp._numerator = [ expr._args[0], expr._args[0] ]
            tmp._denominator = []
            expr._args = (tmp, as_numeric(1))
            #quad.append( (tmp, tmp._args[0]) )
            self._collect_bilinear(tmp, bilin, quad)
            return
        # All other expression types
        for e in expr._args:
            self._collect_bilinear(e, bilin, quad)

