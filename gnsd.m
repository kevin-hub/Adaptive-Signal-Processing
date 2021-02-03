function [y,e,w] = gnsd(x,d,u,g,p)

[order, N] = size(x);
w = zeros(order,N+1);
y = zeros(1,N);
e = zeros(1,N);
r = ones(1,N+1)/u;

for i = 1:N
    y(i) = w(:,i)' * x(:,i);
    e(i) = d(i) - y(i);
    
    w(:,i+1) = (1-g/r(i))*w(:,i) + 1 / (r(i) + x(:,i)' * x(:,i)) * e(i) * x(:,i);
    
    if i>1 %as reg requires past samples (i.e. 2 previous samples to calc current sample)
        r(i+1) = r(i) - p * u * e(i) * e(i-1) * x(:,i)' * x(:,i-1) / (r(i-1) + x(:,i-1)' * x(:,i-1)) ^ 2;
    end
end

