function dxdt = epileptor(t,x,x0,y0,Irest1,Irest2,tau0,tau2,gamma,vals)
%% Epileptor model

x1 = x(1); y1 = x(2); z = x(3); x2 = x(4); y2 = x(5); g = x(6);

dx1_dt = y1 - f1(x1,x2,z) - z + Irest1;
dy1_dt = y0 - 5 * x1^2 - y1;
dz_dt  = 1/tau0 * (4 * ( x1 - x0 ) - z);
dx2_dt = - y2 + x2 - x2^3 + Irest2 + 2 * g - 0.3 * (z - 3.5);
dy2_dt = (1/tau2) * (- y2 + f2(x2));
dg_dt  = - gamma * (g - 0.1 * x1);

dxdt = [dx1_dt;dy1_dt;dz_dt;dx2_dt;dy2_dt;dg_dt];

end

function val = f1(x1,x2,z)
if x1 < 0
    val = x1^3 - 3 * x1^2;
else
    val = (x2 - 0.6 * (z - 4)^2) * x1;
end

end

function val = f2(x2)
if x2 < -0.25
    val = 0;
else
    val = 6 * (x2 + 0.25);
end

end
