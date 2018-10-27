N = 100; %Number of grid points (J in textbook)
x = 0:1/(N-1):1; %values of x at grid points
dx = 1/(N-1); %Spatial step
s = .99 ; % < 1 for stability
dt = sqrt(s)*dx; %Time step for stability

T = 30; %Final time
M = floor(T/dt); %Number of iterations

u0 = x.*(1-x).*exp(-10*(x-0.5).^2); %Initial wave profile
u1 = x.^3 + sin(x); %Initial Velocity

ymin = min(u0);
ymax = max(u0);

un = [u0(2:end),pi]; %u(i+1)
up = [pi,u0(1:end-1)]; %u(i-1)
u_prev = u0;
u = s*(un + up)/2 + (1-s)*u0 + 5*u1*dt ; %last term is non-zero initial velocity

plot(x,u);
pause;
key_pressed=0;
set(gcf, 'KeyPressFcn', 'key_pressed=mod(key_pressed+1,2);')

%Press any key to stop/start simulation
for i=1:M
    
    un = [u(2:end),pi]; %u(i+1)
    up = [pi,u(1:end-1)]; %u(i-1)
    
    uxx = (un - 2*u + up)/dx^2;
    v = 2*u - u_prev + dt^2*uxx;
    u_prev = u;  %Update previous values of u
    u = v;  %Update current value of u
    
    %Draw u
    plot([-dx,x,1+dx],[0,u,0]);  %Draw the string plus the fixed u=0 endpoints
    ymin = min(min(u),ymin);
    ymax = max(max(u),ymax);
    ylim([ymin ymax]);
    drawnow;
    
    %Pause when key is pressed
    while key_pressed == 1
        pause(0.1);
    end
end
