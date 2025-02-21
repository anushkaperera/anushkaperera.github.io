from pyomo.environ import *
from pyomo.dae import *


mu = 1.0

m = ConcreteModel()
m.t = ContinuousSet(bounds=(0.0, 1.0))
m.x = Var(m.t)
m.y = Var(m.t)
m.x_dot = DerivativeVar(m.x)
m.y_dot = DerivativeVar(m.y)

@m.Constraint(m.t)
def ode_x(m, t):
    return m.x_dot[t] == m.y[t]

@m.Constraint(m.t)
def ode_y(m, t):
    return m.y_dot[t] == mu * (1.0 - m.x[t] * m.x[t]) * m.y[t] - m.x[t]

m.pc = ConstraintList()
m.pc.add(m.x[0]==0)
m.pc.add(m.y[0]==0)

TransformationFactory('dae.finite_difference').apply_to(m, wrt=m.t, nfe=30)
SolverFactory('ipopt').solve(m).write()
