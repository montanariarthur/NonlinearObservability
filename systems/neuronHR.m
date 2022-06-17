function dx = neuronHR(t,x,a,b,c,d,I,r,s,xR)

dx =[x(2) - a*x(1)^3 + b*x(1)^2 + I - x(3);
     c - d*x(1)^2 - x(2);
     r*(s*(x(1)-xR) - x(3))];