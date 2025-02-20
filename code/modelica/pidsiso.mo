model pidsiso
  // Model parameters: Time constant and gain
  parameter Real tau = 1;
  parameter Real kp = 1;

  // Input and state variables
  Real u;
  Real x(start=-10);

  // PI controller parameters
  parameter Real k = 1;
  parameter Real Ti = 1;
  Modelica.Blocks.Continuous.LimPID PID(controllerType = Modelica.Blocks.Types.SimpleController.PI, k = k, Ti = Ti);
equation
  PID.u_s = 1;
  PID.u_m = x; //u = if time < 10 then 0 else 1; // for open loop test
  u = PID.y;  
  tau * der(x) + x = kp * u;
end pidsiso;
