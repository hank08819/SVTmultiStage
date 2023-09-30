%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SVT function
% Authors: Xuemei Yang, Chun Li, Henry Han
% Last update: July 07, 2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ X,iterations ] = SVT(M,P,T,delta,itermax,tol)

% Single value thresholding algorithm SVT
% function solves the following optimization problem
%                  min  ||X||*
%               s.t. P(X-M) = 0
% X: recovered matrix
% M: observed matrix
% T: single value threshold
% delta: step size
% output£ºX,iterations

% initialization
Y = zeros(size(M));
iterations = 0;

if nargin < 3
    T =  sqrt(n1*n2);
end
if nargin < 4
    delta = 1;
end
if nargin < 5
    itermax = 1000 ;
end
if nargin < 6
    tol = 1e-5;
end

for ii = 1:itermax
    % singular value threshold operation
    [U, S, V] = svd(Y, 'econ') ;
    S = sign(S) .* max(abs(S) - T, 0) ;
    X = U*S*V' ;
    % update auxiliary matrix Y
    Y = Y + delta* P.* (M-X);
    Y = P.*Y ;
    % computer error
    error= norm( P.* (M-X),'fro' )/norm( P.* M,'fro' );
    if error<tol
        break;
    end
    % update iterations
    iterations = ii ;
    if mod(ii,10000)==0
        fprintf('  *');
    end
end
error;
end
