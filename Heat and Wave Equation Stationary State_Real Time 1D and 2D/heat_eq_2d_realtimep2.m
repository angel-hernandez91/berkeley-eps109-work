N = 200;
L = 10;
x = 0: L/N: L;
f = 0;

%initialize vectors for boundary conditions
T = MidpointCircle(0, 25, N - 10, N - 10, 100);


%time step stuff
eta = 0.24;
dx = L/N;
k = 1;
dt = eta/(25*k);

for k = 1: 3000
    T_new = T;
    
    for i = 2: N 
        for j = 2: N
            T_new(i, j) = (.24)*(T(i -1, j) + T(i +1, j) + T(i, j -1) + T(i, j + 1) -4*T(i, j)) + T(i, j);
        end
    end
    T = T_new;
    %if mod(k + 9, 10) == 0
        xx = 0:L/(size(T, 1) -1):L;
        yy = 0:L/(size(T, 1) -1):L;
        [X, Y] = meshgrid(yy, xx);
        f = f + 1;
        surf(X, Y, T)
        axis([0 L 0 L 0 100])
        saveas(gcf, ['heat_', sprintf('%03d', f), '.png'])
        pause(0.01)
    %end
end
    