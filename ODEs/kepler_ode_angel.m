function dYdt = kepler_ode_angel(t,Y)
% assume 'Y' is a vector of form y = [ x ;  y ; vx ; vy ]
% input:
% Y = [xS; yS; vSx; vSy; xP; yP; vPx; vPy; xA; yA; vAx; vAy]; 
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
%F = - mSun * mPlanet * G / norm(rP - rS)^3 .* (rP - rS);
%a = [F./mSun , F ./ mPlanet, F./mAsteroid];

%from stephanie -- you need to calculate SIX forces and then get the
%accelerations!

F_SP = - mSun * mPlanet * G / (norm(rS-rP))^3 .* (rS-rP);
F_PS = -F_SP;

F_SA = - mSun * mAsteroid * G / (norm(rS-rA))^3 .* (rS-rA);
F_AS = -F_SA;

F_PA = - mPlanet * mAsteroid * G / (norm(rP-rA))^3 .* (rP-rA);
F_AP = -F_PA;

aS = (F_SP+F_SA)./mSun;
aP = (F_PS+F_PA)./mPlanet;
aA = (F_AS+F_AP)./mAsteroid;
dYdt = [ vS(1) ; vS(2) ; aS(1) ; aS(2) ; vP(1); vP(2); aP(1); aP(2); vA(1); vA(2); aA(1); aA(2)];  % must column vector 


        
