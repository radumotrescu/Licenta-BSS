function [grad, delmyW, delb] = wgradientbeta(eta, kappa, b, Y, W)
%Gradient function, it gives the value by which W must be updated to get new W
% \del W = eta (I + (1-2Z)Y')W where
%	eta - learning rate
%	Z = sigmoid(Y)
%	W is the predicted matrix
%	b is a column vector with one value for each source

delW = zeros(size(W));

B = repmat(b, 1, size(Y,2));
Z = sigmoidb(B,Y);
Z([1,2],1:3);
Id = eye(size(Y,1));
#Id
grad = eta * (Id + (B .* (1-2*Z))*Y') * W;
#eta
#(B.*(1-2*Z))*Y'
#size(W)
#size(grad)

size(B .* (1-2*Z));
delmyW = eta * (((B .* (1-2*Z))*Y') +pinv(W)') ;
delb = kappa * (sum((Y .* (1-2*Z) *Y')'))' ;

end
