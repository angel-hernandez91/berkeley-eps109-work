N = 50;
L = 10;
x = 0: L/N: L;
f = 0;
%creation of boundary conditions
T = zeros(N+1, N + 1);
T(N/2, 2:N) = 100;
T(2:N, 2) = 100;
T(2:N, N) = 100;


      

%time step stuff
eta = 0.24;
dx = L/N;
k = 1;
dt = eta/(25*k);

for k = 1: 5000
    T_new = T;
    
    for i = 2: N 
        for j = 2: N
            T_new(i, 1) = T_new(i, 2);
            T_new(N +1, i) = T_new(i, N);
            T_new(1, j) = T_new(2, j);
            T_new(N +1, j) = T_new(N, j);
            T_new(i, j) = (.24)*(T(i -1, j) + T(i +1, j) + T(i, j -1) + T(i, j + 1) -4*T(i, j)) + T(i, j);
            
        end
    end
    T = T_new;
%     if mod(k + 9, 10) == 0
%         f = f + 1;
        xx = 0:L/N:L;
        yy = 0:L/N:L;
        [X, Y] = meshgrid(yy, xx);
        surf(X, Y, T)
        axis([0 L 0 L 0 50])
%         saveas(gcf, ['heat_', sprintf('%03d', f), '.png'])
          pause(0.1)
%     end
end
    