model pidsiso
  // Model parameters: Time constant and gain
  parameter Real tau = 1;
  parameter Real kp = 1;
  parameter Real taud = 1;
  
  // Steady state
  parameter Real u_ss = 1;
  parameter Real x_ss = kp * u_ss;
  
  // Input and state variables
  Real u;
  Real x(start=x_ss);

  // PI controller parameters: Use some tuning rules here, for example, SIMC.
  parameter Real tauc = taud;
  parameter Real k = (1 / kp) * tau / (tauc + taud);
  parameter Real Ti = min(tau, 4 * (tauc + taud));
  
  // PI controller state and setpoint
  Real z(start=(Ti / k) * u_ss);
  Real sp;
equation
  // Controller
  sp = if time < 10 then x_ss else 1.5 * x_ss;
  der(z) = sp - x;
  u = k * (sp - x) + (k / Ti) * z;
  
  // Model
  tau * der(x) + x = kp * delay(u, taud);
end pidsiso;
