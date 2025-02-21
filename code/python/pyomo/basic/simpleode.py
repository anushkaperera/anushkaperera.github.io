from pyomo.environ import *
from pyomo.dae import *
from pyomo.dae.simulator import Simulator


tau = 1.0
k = 1.0

m = ConcreteModel()
m.t = ContinuousSet(bounds=(0.0, 10.0))
m.x = Var(m.t)
m.y = Var(m.t)
m.dx = DerivativeVar(m.x)
m.dy = DerivativeVar(m.y)

m.x[0.0].fix(0.0)
m.y[0.0].fix(0.0)

def u(t):
    return 1.0

@m.Constraint(m.t)
def ode1(m, t):
    return tau * m.dx[t] + m.x[t] == k * u(t)

@m.Constraint(m.t)
def ode2(m, t):
    return m.y[t] == m.x[t]

tsim, profiles = Simulator(m, package='casadi').simulate(numpoints=100)
