N = 50;
L = 10;

T = zeros(101, N+1);
T(2:100, 1) = 8;
T(1, 2:floor(N/2)) = 3;
T(1, (N+ 2)/2 + 1: N) = 5;
T(2: 100, N) = 0;
T(101, 2:floor(N/2)) = 2;
T(101, (N + 2)/2 + 1: N) = 5;
T(1, 1) = 5.5;
T(101, 1) = 5;
T(1, N+1) = 2.5;
T(101, N+1) = 2.5;
T(101, (N + 2)/2) = 3.5;
T(1, (N + 2)/2) = 4;

for k = 1: 5000
    T_new = T;
    for i = 2: N
        for j = 2: N
            T_new(i, j) = (1/4)*(T(i -1, j) + T(i +1, j) + T(i, j -1) + T(i, j + 1));
        end
    end
    T = T_new;
    xx = 0:100;
    yy = 0:L/N:L;
    [X, Y] = meshgrid(yy, xx);
    surf(X, Y, T);
    view([+34.5 14]);
    title('s_i = 25053830')
    pause(0.1);
end
pcolor(X, Y, T)
shading interp
hold on
contour(X, Y, T, 30, 'k');
    