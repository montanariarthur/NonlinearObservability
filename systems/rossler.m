function dx = rossler(t,x,a,b,c)

dx =[-x(2) - x(3);
    x(1) + a*x(2); 
    b + x(3)*(x(1)-c)];

