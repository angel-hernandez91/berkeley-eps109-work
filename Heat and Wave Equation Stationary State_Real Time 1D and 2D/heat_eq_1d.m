N = 25;
L = 10;

x = [0: L/N: L];
T = zeros(1, N + 1);
T(1) = 2;
T(N + 1) = 5;


for k = 1:500
    T_new = T;
    for i = 2:N
        if mod(N, 2) == 0
            T_new(N/2) = 0;
        else
            T_new(N/2 + 1/2) = 0;
        end
        T_new(i) = (1/2)*(T(i -1) + T(i + 1));
    end
    T = T_new;
    plot(x, T)
    title('s_1 = 2, s_2 = 5, s_3 = 0')
    pause(0.1)
end
