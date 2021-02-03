function [y,e,w] = tanhptron(x,d,u,g,a,w)

[order, N] = size(x);
w(order,N+1) = 0;
y = zeros(1,N);
e = zeros(1,N);

for i = 1:N
    y(i) = a*tanh(w(:,i)' * x(:,i));
    e(i) = d(i) - y(i);
    w(:,i+1) = (1-u*g)*w(:,i) + u .* e(i) * x(:,i);
end

w = w(:,2:end);

end