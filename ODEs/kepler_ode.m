function dYdt = kepler_ode(t,Y)
% assume 'Y' is a vector of form y = [ x ;  y ; vx ; vy ]

r = [Y(1);Y(2); Y(5); Y(6); Y(9); Y(10)];
v = [Y(3);Y(4); Y(7); Y(8); Y(11); Y(12)];
rS = [r(1), r(2)];
vS = [v(1), v(2)];
rP = [r(3), r(4)];
vP = [v(3), v(4)];
rA = [r(5), r(6)];
vA = [v(5), v(6)];

% use the following three variables from the main worksheet
% must declare them 'global' here and there
global mSun mPlanet G  mAsteroid

% assumes the sun is located at r=(0 0)
FSP = - mSun * mPlanet * G / norm(rP - rS)^3 .* (rP - rS);
a1 = [FSP./mSun , FSP ./ mPlanet];
FAP = - mPlanet*mAsteroid*G/ norm(rP - rA)^3 .* (rA - rP);
a2 = FAP./mAsteroid;

dYdt = [ vS(1) ; vS(2) ; a1(1) ; a1(2) ; vP(1); v(2); a1(3); a1(4); vA(1); vA(2); a2(1); a2(2)];  % must column vector 


        
