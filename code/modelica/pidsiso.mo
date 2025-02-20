model pidsiso
  // Model parameters: Time constant and gain
  parameter Real tau = 1;
  parameter Real kp = 1;

  // Steady state
  parameter Real u_ss = 1;
  parameter Real x_ss = kp * u_ss;
  
  // Input and state variables
  Real u;
  Real x(start=x_ss);

  // PI controller parameters
  parameter Real k = 1;
  parameter Real Ti = 1;
  Modelica.Blocks.Continuous.LimPID PID(controllerType = Modelica.Blocks.Types.SimpleController.PI, k = k, Ti = Ti, initType = Modelica.Blocks.Types.Init.InitialOutput, y_start = x_ss);
equation
  // Controller
  PID.u_s = if time < 10 then 1 else 1.5;
  PID.u_m = x;
  PID.y = u;
  // Model
  tau * der(x) + x = kp * u;
end pidsiso;
