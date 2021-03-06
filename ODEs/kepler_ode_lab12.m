function dYdt = kepler_ode_lab12(i, Y)
% assume 'Y' is a vector of form y = [ x, y, z ]

global a b c;
a = 0.125053830;
b = a;


dx = -Y(2) - Y(3);
dy = Y(1) + a*Y(2);
dz = b + Y(3)*(Y(1) - c(i));

dYdt = [dx; dy; dz];