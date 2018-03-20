function g= sigmoidb(B, y)
% The sigmoid/logistic function
% g(y) = 1/(1 + e^{-by})
% Returns the sigmoid of y, where y can be a single value, vector or a matrix
% If b and y are single values we directly multiply them. If b is a scalar and y is a vector, we multiply the scalar with the vector.
% B is the same dimension as y!

g = zeros(size(y));
y([1,2],[1:3])
(-(B.*y))([1,2],1:3)
exp(-(B.*y))([1,2],1:3)

g = 1 ./ (1 + exp(- (B .* y)));

end
