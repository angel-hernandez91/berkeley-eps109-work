N = 30; %Number of grid points (in x and y directions)
[x, y] = meshgrid(0:1/(N-1):1); %values of (x,y) at grid points
dx = 1/(N-1); %Spatial step
dy = 1/(N-1);
err = dx^2; %Error tolerance at which iterations stop
            %Since finite differences are O(dx^2) accurate
            %there is no need to solve scheme at an accuracy beyond dx^2

f = (x -5).^2 ; %RHS f
u0 = exp(10*x) + exp(10*y); %Boundary values

zmin = min(u0(:));
zmax = max(u0(:));

u=u0;

surf(x,y,u);
pause;
key_pressed=0;
set(gcf, 'KeyPressFcn', 'key_pressed=mod(key_pressed+1,2);')

e=1;
i=1;

%Press any key to stop/start simulation
while e > err  %Loop until scheme error e is less than err=dx^2
    
    i
    v = (u(3:end,2:end-1) + u(1:end-2,2:end-1) + u(2:end-1,3:end) + u(2:end-1,1:end-2) + dx^2*f(2:end-1,2:end-1))/4;
    
    %Compute error in solving scheme
    e = max(max(abs(u(2:end-1,2:end-1) -v)))/dx^2
    u(2:end-1,2:end-1) = v;
    
    %Draw u
    surf(x,y,u);
    zmin=min(min(u(:)),zmin);
    zmax=max(max(u(:)),zmax);
    zlim([zmin zmax]);
    drawnow;
    
    %Pause when key is pressed
    while key_pressed == 1
        pause(0.1);
    end
    i = i+1;
end
