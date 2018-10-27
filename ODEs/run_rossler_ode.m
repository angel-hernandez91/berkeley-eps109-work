clear
clf
global a b c;
a = 0.125053830;
b = a;
c = 13;
% cstep = 0.1;
% for c = 2:cstep:10
    
    % initial conditions
    x = 0.1;
    y = 0.2;
    z = 0.3;
    
    Y = [x; y; z]; % by convention, should be a column vector
    
    
    t    =  0.0;
    tmax = 200.0;
    dt   =  1/100;
    n    = tmax/dt; % number of steps
    
    tt = zeros(n+1,1);
    Yt = zeros(n+1,size(Y,1));
    tt(1)   = t;
    Yt(1,:) = Y';
    
    for i=1:n
        f1 = rossler_ode(i,Y                );
        f2 = rossler_ode(i,Y+f1 .*  (dt/2)  );
        f3 = rossler_ode(i,Y+f2 .*  (dt/2)  );
        f4 = rossler_ode(i,Y+f3 .*   dt     );
        
        Y = Y + (f1 + f2.*2 + f3.*2 + f4) .* (dt/6);
        
        t = t+dt;
        tt(i+1)   = t;  % this records the time at each step
        Yt(i+1,:) = Y'; % this records the vector Y at each step --> Yt(:,1) has all x coordinates
        
%         if sign(Yt(i, 2)) ~= sign(Yt(i+1, 2)) && tt(i + 1) > 100
%             plot(c, Yt(i, 1), 's', 'markersize', 3,'MarkerFaceColor','b','MarkerEdgeColor','b')
%             hold on
%         end
        
    end
% end


subplot(2, 2, 1);
plot(tt(15000:end), Yt(15000:end,1));
title('x(t) versus t with c = 13')
subplot(2, 2, 2)
plot(Yt(15000:end, 1), Yt(15000:end, 2));
title('x(t) versus y(t) with c = 13')
subplot(2, 2, 3);
plot3(Yt(15000:end, 1), Yt(15000:end, 2), Yt(15000:end,3));
title('x(t), y(t), and z(t) with c = 13')

% plot(0,0,'X',Yt(:,1),Yt(:,2),'-or');



