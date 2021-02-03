function [y,e,w] = gass(x,d,u,g,p,method)

[order, N] = size(x);
w = zeros(order,N+1);
y = zeros(1,N);
e = zeros(1,N);
u = [u zeros(1,N)];
phi = zeros(order,N+1);

for i = 1:N
    y(i) = w(:,i)' * x(:,i);
    e(i) = d(i) - y(i);
    
    w(:,i+1) = (1-u(i)*g)*w(:,i) + u(i) .* e(i) * x(:,i);
    u(i+1) = u(i) + p * e(i) * x(:,i)' * phi(:,i);
    
    switch method.name
        case 'B'
            phi(:,i+1) = (eye(order) - u(i)*x(:,i)*x(:,i)')*phi(:,i)+e(i)*x(:,i);
        case 'AF'
            phi(:,i+1) = method.param * phi(:,i)+e(i)*x(:,i);
        case 'MX'
            phi(:,i+1) = e(i) * x(:,i);
    end
    
end

w = w(:,2:end);

end

