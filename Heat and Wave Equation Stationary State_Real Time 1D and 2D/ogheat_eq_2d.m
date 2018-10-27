N = 51;
L = 10;

T = zeros(N+1, N+1);
T(:, 1) = 100;

for k = 1: 5000
    T_new = T;
    for i = 2: N
        for j = 2: N
            T_new(i, j) = (1/4)*(T(i -1, j) + T(i +1, j) + T(i, j -1) + T(i, j + 1));
        end
    end
    T = T_new;
    xx = 0:L/N:L;
    yy = 0:L/N:L;
    [X, Y] = meshgrid(yy, xx);
    surf(X, Y, T);
    view([+34.5 14]);
    pause(0.1);
end
pcolor(X, Y, T)
shading interp
hold on
contour(X, Y, T, 30, 'k');
    