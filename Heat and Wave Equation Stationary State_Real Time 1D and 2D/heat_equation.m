N = 100; %Number of grid points (J in textbook)
x = 0:1/(N-1):1; %values of x at grid points
dx = 1/(N-1); %Spatial step
s=.49;  %<0.5 for stability
dt = s*dx^2; %Time step for stability

T = 0.1; %Final time
M = floor(T/dt); %Number of iterations

u0 = x.^2;
%u0 = 1./(x - .5  ); %Initial heat distribution


ymin = min(u0);
ymax = max(u0);

u=u0;


plot(x,u);
pause;
key_pressed=0;
set(gcf, 'KeyPressFcn', 'key_pressed=mod(key_pressed+1,2);') 

%Press any key to stop/start simulation
for i=1:M
    
    un = [u(2:end),0]; %u(i+1)
    up = [0,u(1:end-1)]; %u(i-1)
    
    uxx = (un - 2*u + up)/dx^2;
    u = u + dt*uxx;
    
    %Draw u
    plot(x, u) %Nuemann BC Plot
    plot([-dx,x,1+dx],[0,u,0]);
    ylim([ymin ymax]);
    drawnow; 
    
    %Pause when key is pressed
    while key_pressed == 1
        pause(0.1);
    end
end
