model pidsiso
  parameter Real tau = 1;
  parameter Real kp = 1;
  Real u;
  Real x(start=-10);
  parameter Real k = 1;
  parameter Real Ti = 1;
  parameter Boolean auto = false;
  Modelica.Blocks.Continuous.LimPID PID(controllerType = Modelica.Blocks.Types.SimpleController.PI, k = k, Ti = Ti);
equation
  PID.u_s = 1;
  PID.u_m = x;
  //u = if time < 10 then 0 else 1; // for open loop test
  u = PID.y;  
  tau * der(x) + x = kp * u;
annotation(
    uses(Modelica(version = "4.0.0")));
end pidsiso;
