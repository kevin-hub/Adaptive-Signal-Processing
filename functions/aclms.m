function [p,e,h,g] = aclms(x,d,u,y)

[order, N] = size(x);
p = zeros(1,N);
e = zeros(1,N);
h = zeros(order,N+1);
g = zeros(order,N+1);

for i = 1:N
    p(i) = h(:,i)' * x(:,i) + g(:,i)' * conj(x(:,i));
    e(i) = d(i) - p(i);
    h(:,i+1) = (1-u*y)*h(:,i) + u .* conj(e(i)) * x(:,i);
    g(:,i+1) = (1-u*y)*g(:,i) + u .* conj(e(i)) * conj(x(:,i));
end

h = h(:,2:end);
g = g(:,2:end);

end