load('sounds.mat')
load('icaTest.mat')
A(1:2,1:2);
m1=audioread('init-music-1.wav');
m2=audioread('init-music-2.wav');

e1=audioread('expected-music-1.wav');
e2=audioread('expected-music-2.wav');

M=[m1,m2];
B=[1.1 , 0.7 ; 1.5 , 0.5];
 
M=[e1,e2];
%A=rand(numSrc, numSrc);
%M*A;

%X=B*M.';
X=B*M.';
%X=M.';
subplot (4, 1, 1)
audiowrite('miii1.wav',X(1,:),16000)
plot(X(1,:))
subplot (4, 1, 2)
audiowrite('miii2.wav',X(2,:),16000)
plot(X(2,:))

eta = 0.01;
eta0 = eta;
T=1000;

%Make some random guess of mix-matrix inverse
W = rand(size(B))./50;
W
size(X)
numSrc = size(X,1);
numSrc

b = ones(numSrc,1);
b
kappa=0.01;	%0.0001-no change, 0.01-good
num_iter=1;

size(W)
size(X)
for i=0:num_iter,
	Y = W*X;			% predict source matrix based on guessed mix matrix
  size(Y);
	[delW, delmyW, delb] = wgradientbeta(eta, kappa, b, Y, W);	% gradient descent - shift by delta
	W = W + delW;			% update W

	%W = W + (delmyW * 0.001);			% update W
	eta = eta0 / (1 + (i/T));	% annealing - learning rate
  eta;
	if(mod(i,100)==0),
		b;
    delb;
		b = b + delb;
    W
    delW
		%fprintf('runs %d \n',i);
		W;
	%	corrMat = correlations(srcMat,Y)
	%	%fflush(stdout);
	end;
end;


W
Y = W*X;				% predict source matrix based on guessed mix matrix

Y = (Y - min(min(Y))) ./ (max(max(Y)) - min(min(Y)));
Y2=Y;
%Y2 = Y .* 2.0;

subplot (4, 1, 3)
plot(Y2(1,:))
subplot (4, 1, 4)
plot(Y2(2,:))


audiowrite('m1.wav',Y2(1,:),16000)
audiowrite('m2.wav',Y2(2,:),16000)



%for i=0:num_iter,
%	Y = W*X;			% predict source matrix based on guessed mix matrix
%	[delW, mygrad] = gradient(eta, Y, W);	% gradient descent - shift by delta
%	W = W + delW;			% update W
%	%W = W + (mygrad * 0.001);			% update W
%	eta = eta0 / (1 + (i/T));	% annealing - learning rate
%	%if(mod(i,1000)==0),
%	%	fprintf('runs %d \n',i);
%	%	W
%	%	corrMat = correlations(srcMat,Y)
%	%	%fflush(stdout);
%	%end;
%end;
