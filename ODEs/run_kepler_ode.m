clear
clf
global mSun mPlanet G mAsteroid
mSun    = 1e+7;
mPlanet = 1e+4; % this value does not matter as long as mPlanet<<mSun
G       = 1;
mAsteroid = 1;

mu = (mSun*mPlanet)/(mSun + mPlanet);
R = 20 ; %abs(rSun - rPlanet)
w = sqrt(G *mSun*mPlanet/R^3* (mu));
rCM = [(mPlanet/(mSun + mPlanet))*sqrt(20), 0]; % center of mass
% initial planet position
xP = sqrt(20) - rCM(1);
yP = 0 - rCM(2);


% initial sun position
xS = -rCM(1);
yS = -rCM(2);


% initial asteroid position
xA = R*cos(45);
yA = R*sin(45);

% initial planet velocities
vPx = 0;
vPy = w*xP;

% initial asteroid velocity
vAx = vPx;
vAy = vPy;


%Sun Velocity 
vSx = 0;
vSy = w*xS;

Y = [xS; yS; vSx; vSy; xP; yP; vPx; vPy; xA; yA; vAx; vAy]; % by convention, should be a column vector


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

    f1 = kepler_ode(t     , Y                );
    f2 = kepler_ode(t+dt/2, Y+f1 .*  (dt/2)  );
    f3 = kepler_ode(t+dt/2, Y+f2 .*  (dt/2)  );
    f4 = kepler_ode(t+dt  , Y+f3 .*   dt     );
   
    Y = Y + (f1 + f2.*2 + f3.*2 + f4) .* (dt/6);
    
    t = t+dt;
    tt(i+1)   = t;  % this records the time at each step
    Yt(i+1,:) = Y'; % this records the vector Y at each step --> Yt(:,1) has all x coordinates
end

plot(0,0,'X',Yt(:,1),Yt(:,2),'-or');
title('Kepler XY Plot ','fontsize',20);
axis equal


