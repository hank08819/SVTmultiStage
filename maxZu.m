function [ X2max ] = maxZu(YcellZu, X1, P,delta,cishu,tol)
% 对最大块/组 二次填充
 
M2k = X1(YcellZu,:) ; P2k = P(YcellZu,:); [a1, a2] = size(M2k); T = sqrt(a1*a2); 
[ X2k,iterations ] = SVT(M2k,P2k,T, delta,cishu,tol) ;  % 
X2max = X2k;