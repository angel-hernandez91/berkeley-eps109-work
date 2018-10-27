clear
%host rock/dike spacing variables, default var = 2
var = 2;
d = 4;

%intervals
N = 50;
M = 5*N;
L = 2*(d* var) + d;
x = 0: L/M: L;
T = zeros(size(x));

%step size
dx = L/N;

%heat coefficient
k = 8.8183e-05;

%time step
eta = 0.49;
dt = (eta*(dx^2))/k;

%initial Conditions
T(2*N + 1: 3*N + 1) = 1175;

iterations = 3600*24/dt;

for m = 1:floor(iterations)
    T_new = T;
    for i = 2:M
        T_new(i) = eta*(T(i + 1) - 2*T(i) + T(i -1)) + T(i);
    end
    T = T_new;
%     plot(x, T)
%     axis([0 L 0 1200]);
%     drawnow;
%     pause(0.1)
end



%integration
%endpts
endpts = T(2*N + 1) + T(3*N + 1);
I = dx*(endpts +2*sum(T(2*N + 2: 3*N)))/2;
disp(I)
disp(T(M/2 + 1))
disp(I/T(M/2 + 1))

