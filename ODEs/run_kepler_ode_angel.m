clear
clf
global mSun mPlanet G mAsteroid
mSun    = 1e+7;
mPlanet = 1e+4; % this value does not matter as long as mPlanet<<mSun
G       = 1;
mAsteroid = 1;

mu = (mSun*mPlanet)/(mSun + mPlanet);
R = 20 ; %abs(rSun - rPlanet)
%w = sqrt(G *mSun*mPlanet/R^3* (mu));
% from stephanie -- your w is very wrong order of magnitude
w = sqrt( (G*mSun*mPlanet) / (mu*R^3) );

%rCM = [(mPlanet/(mSun + mPlanet))*sqrt(20), 0]; % center of mass
% from stephanie -- your rCM was wrong.
rSun = [0 0]; rPlanet = [20 0];
rCM = (rPlanet*mPlanet + rSun*mSun) * (1/ (mPlanet + mSun));
% initial planet position
%xP = sqrt(20) - rCM(1);
%yP = 0 - rCM(2);
% from stephanie - no, the planet is initially [20 0] in sun-centered
% system or [19.98 0] in center-of-mass-centered system!
% why take square root of 20?
xP = 20-rCM(1); yP = 0;

% initial sun position
%xS = -rCM(1);
%yS = -rCM(2);

% from stephanie
xS = -rCM(1); yS = 0; %look at the picture.

% initial asteroid position
% from stephanie --- USE cosd and sind!!
%xA = R*cos(45);
%yA = R*sin(45);
xA = (R-rCM(1))*cosd(45);
yA = (R-rCM(1))*sind(45);

% initial planet velocities
vPx = 0;
vPy = w*xP;

% initial asteroid velocity
%vAx = vPx;
%vAy = vPy;
% from steph -- asteroid velocities are cross product of omega and r
% also check picture -- asteroid velocity should be + in y and - in x
vAx = -w*yA;
vAy = w*xA;

%Sun Velocity 
vSx = 0;
vSy = w*xS; %from steph -- make sure vSy is negative (see picture)

Y = [xS; yS; vSx; vSy; xP; yP; vPx; vPy; xA; yA; vAx; vAy]; % by convention, should be a column vector
% looked at your IC's "Y" -- they are all really big, like a million. this
% may be the problem.
%pause
T = sqrt((4*pi^2)/(G*(mSun + mPlanet))*R^3);
t    =  0.0;
tmax = T/4;
dt   =  0.00001;
n    = floor(tmax/dt); % number of steps

tt = zeros(n+1,1);
Yt = zeros(n+1,size(Y,1));
tt(1)   = t;
Yt(1,:) = Y';


for i=1:n

    f1 = kepler_ode_angel(t     , Y                );
    f2 = kepler_ode_angel(t+dt/2, Y+f1 .*  (dt/2)  );
    f3 = kepler_ode_angel(t+dt/2, Y+f2 .*  (dt/2)  );
    f4 = kepler_ode_angel(t+dt  , Y+f3 .*   dt     );
   
    Y = Y + (f1 + f2.*2 + f3.*2 + f4) .* (dt/6);
    
    t = t+dt;
    tt(i+1)   = t;  % this records the time at each step
    Yt(i+1,:) = Y'; % this records the vector Y at each step --> Yt(:,1) has all x coordinates
end

%plot(0,0,'X',Yt(:,1),Yt(:,2),'-or');
%title('Kepler XY Plot ','fontsize',20);
% from steph - make sure you plot all the things!!
plot(0,0,'X',Yt(:,1),Yt(:,2),'-or',Yt(:,5),Yt(:,6),'k-',Yt(:,9),Yt(:,10),'b-');
title('Kepler XY Plot', 'fontsize',20);
axis equal

