function p = circularity(x)
p = abs(mean(x.^2)/mean(abs(x.^2)));
end

