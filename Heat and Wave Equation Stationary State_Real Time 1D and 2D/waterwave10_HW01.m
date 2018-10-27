function waterwave
% WATER WAVE
% 2D Shallow Water Model
%
% Lax-Wendroff finite difference method.
% Reflective boundary conditions.
% Random water drops initiate gravity waves.
% Surface plot displays height colored by momentum.
% Plot title shows t = simulated time and tv = a measure of total variation.
% An exact solution to the conservation law would have constant tv.
% Lax-Wendroff produces nonphysical oscillations and increasing tv.
%
% See:
%    http://en.wikipedia.org/wiki/Shallow_water_equations
%    http://www.amath.washington.edu/~rjl/research/tsunamis
%    http://www.amath.washington.edu/~dgeorge/tsunamimodeling.html
%    http://www.amath.washington.edu/~claw/applications/shallow/www

% Parameters

clf
clear
f = 0;                   % gif creation
n = 128*1.5;                  % grid size
g = 9.8;                 % gravitational constant
dt = 0.01;               % hardwired timestep
dx = 1.0;
dy = 1.0;
nplotstep = 8;           % plot interval
ndrops = 1;              % maximum number of drops
dropstep = 500;          % drop interval

G = zeros(n+2,n+2);
% G1 = zeros(n+2, n+2);
% G2 = zeros(n +2, n + 2);
for i=1:4
   G = initial_conditions_add_hill(G,n,n,rand()*n,rand()*n,rand()*n,rand()-.5);
end
for i=1:4
   G = initial_conditions_add_hill(G,n,n,rand()*n,rand()*n,rand()*n/5,rand()-.5);
end


% % for i=1:(n+2) 
% %   for j=1:(n+2)
% %       G1(i,j) = -1.5 + (j -1)*(2/(n+1)); % Enter you ocean floor elevation her
% %       G2(i,j) = -1.5 + (n+1 - j)*(2/(n+1));
% % Ridge 
% %       G(i, j) = min([G1(i, j), G2(i, j)])
% % Trench
% %       G(i,j) = max([G1(i,j), G2(i,j)]); 
% %       G(i, j) = G(i, j) - 1; 
% %       
% %   end
% % end
x=0:1/(n +1): 1;
      surf(x, x, G);
      pause  

% G(:,n/2-5:n/2+5) = -0.75;
% G(1:n+2,:) = -10;

% Initialize graphics

[surfplot,top] = initgraphics(n);

H = ones(n+2,n+2);   U = zeros(n+2,n+2);  V = zeros(n+2,n+2);
Hx = zeros(n+1,n+1); Ux = zeros(n+1,n+1); Vx = zeros(n+1,n+1);
Hy = zeros(n+1,n+1); Uy = zeros(n+1,n+1); Vy = zeros(n+1,n+1);
ndrop = ceil(rand*ndrops);
nstep = 0;

% Set inital conditions
for i=2:n+1 
    for j=2:n+1
       H(i,j) = H(i,j)+5*exp(-0.2*i^2);
    end
end

% loop over time steps.
while (nstep<3000)   

   nstep = nstep + 1

   % Reflective boundary conditions
   H(:,1) = H(:,2);      U(:,1) = U(:,2);       V(:,1) = -V(:,2);
   H(:,n+2) = H(:,n+1);  U(:,n+2) = U(:,n+1);   V(:,n+2) = -V(:,n+1);
   H(1,:) = H(2,:);      U(1,:) = -U(2,:);      V(1,:) = V(2,:);
   H(n+2,:) = H(n+1,:);  U(n+2,:) = -U(n+1,:);  V(n+2,:) = V(n+1,:);

   % First half step

   % x direction
   i = 1:n+1;
   j = 1:n;

   Horg = H;
   H = H - G;
   % height
   Hx(i,j) = (H(i+1,j+1)+H(i,j+1))/2 - dt/(2*dx)*(U(i+1,j+1)-U(i,j+1));
   Gx(i,j) = (G(i+1,j+1)+G(i,j+1))/2;

   % x momentum
   Ux(i,j) = (U(i+1,j+1)+U(i,j+1))/2 -  ...
             dt/(2*dx)*((U(i+1,j+1).^2./H(i+1,j+1) + g.*(H(i,j+1)+H(i+1,j+1))./2.*Horg(i+1,j+1)) - ...
                        (U(i  ,j+1).^2./H(i  ,j+1) + g.*(H(i,j+1)+H(i+1,j+1))./2.*Horg(i  ,j+1)));

   % y momentum
   Vx(i,j) = (V(i+1,j+1)+V(i,j+1))/2 - ...
             dt/(2*dx)*((U(i+1,j+1).*V(i+1,j+1)./H(i+1,j+1)) - ...
                        (U(i,j+1).*V(i,j+1)./H(i,j+1)));

   % y direction
   i = 1:n;
   j = 1:n+1;

   % height
   Hy(i,j) = (H(i+1,j+1)+H(i+1,j))/2 - dt/(2*dy)*(V(i+1,j+1)-V(i+1,j));
   Gy(i,j) = (G(i+1,j+1)+G(i+1,j))/2;

   % x momentum
   Uy(i,j) = (U(i+1,j+1)+U(i+1,j))/2 - ...
             dt/(2*dy)*((V(i+1,j+1).*U(i+1,j+1)./H(i+1,j+1)) - ...
                        (V(i+1,j).*U(i+1,j)./H(i+1,j)));
   % y momentum
   Vy(i,j) = (V(i+1,j+1)+V(i+1,j))/2 - ...
             dt/(2*dy)*((V(i+1,j+1).^2./H(i+1,j+1) + g.*(H(i+1,j+1)+H(i+1,j))./2.*Horg(i+1,j+1)) - ...
                        (V(i+1,j  ).^2./H(i+1,j  ) + g.*(H(i+1,j+1)+H(i+1,j))./2.*Horg(i+1,j  )));

   % Second half step
   i = 2:n+1;
   j = 2:n+1;

   % height
   H(i,j) = H(i,j) - (dt/dx)*(Ux(i,j-1)-Ux(i-1,j-1)) - ...
                     (dt/dy)*(Vy(i-1,j)-Vy(i-1,j-1));
   % x momentum
   U(i,j) = U(i,j) - (dt/dx)*((Ux(i  ,j-1).^2./Hx(i  ,j-1) + g.*(Hx(i-1,j-1)+Hx(i,j-1))./2.*(Hx(i  ,j-1)+Gx(i  ,j-1)) ) - ...
                              (Ux(i-1,j-1).^2./Hx(i-1,j-1) + g.*(Hx(i-1,j-1)+Hx(i,j-1))./2.*(Hx(i-1,j-1)+Gx(i-1,j-1)))) ...
                   - (dt/dy)*((Vy(i-1,j  ).*Uy(i-1,j  )./Hy(i-1,j  )) - ...
                              (Vy(i-1,j-1).*Uy(i-1,j-1)./Hy(i-1,j-1)));
   % y momentum
   V(i,j) = V(i,j) - (dt/dx)*((Ux(i  ,j-1).*Vx(i  ,j-1)./Hx(i  ,j-1)) - ...
                              (Ux(i-1,j-1).*Vx(i-1,j-1)./Hx(i-1,j-1))) ...
                   - (dt/dy)*((Vy(i-1,j  ).^2./Hy(i-1,j  ) + g.*(Hy(i-1,j)+Hy(i-1,j-1))./2.*(Hy(i-1,j  )+Gy(i-1,j  )) ) - ...
                              (Vy(i-1,j-1).^2./Hy(i-1,j-1) + g.*(Hy(i-1,j)+Hy(i-1,j-1))./2.*(Hy(i-1,j-1)+Gy(i-1,j-1))));

   H = H + G;

   % Update plot
   if mod(nstep,nplotstep) == 0
      C = abs(U(i,j)) + abs(V(i,j));  % Color shows momemtum
      t = nstep*dt;
      tv = norm(C,'fro');
      set(surfplot,'zdata',H(i,j),'cdata',C);
      set(top,'string',sprintf('t = %6.2f,  tv = %6.2f',t,tv))
      drawnow
%GIF Creation

%       nplotstep = nplotstep + 20;
%       saveas(gcf, ['heat_', sprintf('%03d', nplotstep), '.png'])
%       pause(0.01)
   end

   if all(all(isnan(H))), break, end  % Unstable, quit
end

%close(gcf)

function [surfplot,top] = initgraphics(n);
% INITGRAPHICS  Initialize graphics for waterwave.
% [surfplot,top,start,stop] = initgraphics(n)
% returns handles to a surface plot, its title, and two uicontrol toggles.

   clf
   shg
   set(gcf,'numbertitle','off','name','Shallow_water')
   x = (0:n-1)/(n-1);
   surfplot = surf(x,x,ones(n,n),zeros(n,n));
   grid off
   axis([0 1 0 1 -1 3])
   caxis([-1 1])
   shading faceted
   c = (1:64)'/64;
   cyan = [0*c c c];
   colormap(cyan)
   top = title('Click restart');
   