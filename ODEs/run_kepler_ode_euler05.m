clear
clf
global c n;
% initial conditions
x = 0.1;
y = 0.2;
z = 0.3;

Y = [x; y; z]; % by convention, should be a column vector


t    =  0.0;
tmax = 10.0;
dt   =  1/100;
n    = tmax/dt; % number of steps

tt = zeros(n+1,1);
Yt = zeros(n+1,size(Y,1));
tt(1)   = t;
Yt(1,:) = Y';
c = zeros(1, n);
c(1) = 2;
cstep = 0.01;
for i=1:n
    c(i+1) = c(i) + cstep;
    for k = 1:100
        
        f1 = kepler_ode_lab12(i, Y                );
        f2 = kepler_ode_lab12(i, Y+f1 .*  (dt/2)  );
        f3 = kepler_ode_lab12(i, Y+f2 .*  (dt/2)  );
        f4 = kepler_ode_lab12(i, Y+f3 .*   dt     );

        Y = Y + (f1 + f2.*2 + f3.*2 + f4) .* (dt/6);

        t = t+dt;
        tt(i+1)   = t;  % this records the time at each step
        Yt(i+1,:) = Y'; % this records the vector Y at each step --> Yt(:,1) has all x coordinates
        
        if sign(Yt(i, 2)) ~= sign(Yt(i+1, 2)) && tt(i + 1) > 0
            plot(c(i), Yt(i, 1), 's', 'markersize', 3,'MarkerFaceColor','b','MarkerEdgeColor','b')
            hold on
            
        end
    end
end


% subplot(2, 2, 1);
% plot(tt(30000:end), Yt(30000:end,1));
% title('c = 20')
% subplot(2, 2, 2)
% plot(Yt(30000:end, 1), Yt(30000:end, 2));
% title('c = 20')
% subplot(2, 2, 3);
% plot3(Yt(30000:end, 1), Yt(30000:end, 2), Yt(30000:end,3));
% title('c = 20')

% plot(0,0,'X',Yt(:,1),Yt(:,2),'-or');



