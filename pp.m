function [output] = pp(input,order,delay)

nSamples = length(input);
output = zeros(order,nSamples);

for i = 1:order
    % grouped samples to approximate the value at certain instant
    output(i,:) = [zeros(1, i+delay-1), input(1: nSamples-(i+delay-1))];
end
end