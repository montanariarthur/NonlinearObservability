function dx = cord(t,x,a,b,F,G,order)

dx =[-x(2).^order - x(3).^order - a*x(1) + a*F;
    x(1)*x(2) - b*x(1)*x(3) - x(2) + G;
    b*x(1)*x(2) + x(1)*x(3) - x(3)];
