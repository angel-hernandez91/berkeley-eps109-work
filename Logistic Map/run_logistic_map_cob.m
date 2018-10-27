clear
clf;


 % waits until you press RETURN

x = 0.5;  % starting point for iteration
r = 3.83;
for i = 1:100
    xx(i) = i;
    yy(i) = x;
    x = logistic_map(x, r); % assign the next value in the iteration to 'x'
end
x = 0: 0.01: 1;
y = x;
fx = r*(x - x.^2)';
plot(x,y, 'k' ) % plot 2
hold on
plot(x, fx, 'k' )
hold on


u = zeros(1, 100);
v = zeros(1, 100);
u(1) = 0.5;
u(2) = 0.5;

for i = 1: 99
    
    v(i+1) = r*(u(i) - u(i)^2);
    u(i+2) = v(i + 1);
    
    hold on
    str = sprintf('r = %f', r);
    title(str);
    plot([u(i), u(i + 1)],[v(i), v(i + 1)], 'r-', 'Markersize', 20)

end

