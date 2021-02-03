function [y,e,w] = clms(x,d,u,g)

[order, N] = size(x);
y = zeros(1,N);
e = zeros(1,N);
w = zeros(order,N+1);

for i = 1:N
    y(i) = w(:,i)' * x(:,i);
    e(i) = d(i) - y(i);
    w(:,i+1) = (1-u*g)*w(:,i) + u .* conj(e(i)) * x(:,i);
end

w = w(:,2:end);

end

