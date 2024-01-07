% Train & Test using tow-stage matrix filling by SVT algorithm,data:CCLE_AA
%input:original data M with missing values
%output:the best filled matrix X2 using tow-stage matrix filling
%output:indicators for measuring filling effect,RMSE2,P2,S2,R2
clc; clear;close all;
fdata = 'rM-294C-23D_noNA.csv';  %  original data
M = csvread(fdata,1,0);   % original data
M = M + 1;  
[dm, dn] = size(M);   
numAll = dm*dn;
 
% construct index matrix P £¨0:missing value£©
P = ones(numAll,1) ;  
pos =xlsread('pos.xls'); 
i=1;p0=floor(numAll/10);%i=1:10, 10-fold cross validation test
P(pos(p0*(i-1)+1:p0*i)) = 0;
P=reshape(P,dm,dn) ;
 
T = sqrt(numAll);% set threshold 
delta = 1.85 ;  % % step size (0, 2)
ddcishu = 10^7;
tol = 1e-6;
 
[ X,iterations ] = SVT(M,P,T,delta,ddcishu,tol) ;%X:predicted value
%Measurement indicators of the first filling,RMSE1,pearson(spearman)
%correlation coefficient,coefficient of determination
rmse=zeros(1,dn);
for j=1:dn
    Vx = X(:,j);
    Vy = M(:,j);
    rmse(j) = norm( (Vx-Vy),'fro' )/ sqrt(dm) ;  
end
RMSE1=mean(rmse);
pp=zeros(1,dn);%pearson correlation coefficient
for j=1: dn
    Vx = X(:,j);
    Vy = M(:,j);
    [Co_r,Co_p] = corr(Vx,Vy,'type','Pearson');
    pp(j) = Co_r;    
end
P1= mean(pp);%pearson correlation coefficient
s=zeros(1,dn);%spearman correlation coefficient
for j=1: dn
    Vx = X(:,j);
    Vy = M(:,j);
    [Co_r,Co_p] = corr(Vx,Vy,'type','spearman');
    s(j) = Co_r;    
end
S1= mean(s);%spearman correlation coefficient
r=zeros(1,dn);%Coefficient of determination 
for j=1:dn
    Vx = X(:,j);
    Vy = M(:,j);Vy_m=mean(Vy');
    r(j) =1-(norm((Vy-Vx),'fro')/ norm((Vy-Vy_m),'fro'))^2;  
end
R1=mean(r);%Coefficient of determination

% the second stage filling 
X1 = X; 
X2 = 0*X1 ;
 
% cluster and block 
Xjilu = cell(21,1); 
qcshu = 0 ;
ZUnum = 8;    % number of groups
TC2flag = 1 ; % marker of the second filling
P22 = zeros(21,dn);s22 = zeros(21,dn);r22=zeros(21,dn); rmse22=zeros(21,dn);
while TC2flag > 0
    for q = 0: 0.05: 1
        qcshu = qcshu +1 ;
        Correlation_M = eye(dm,dm);  
        for i =  1:(dm-1)
            for j = (i+1):dm
                Vx = X1(i,:);  Vy = X1(j,:);
                wq = ( P(i,:)+P(j,:) )/2 ; %If a certain element of wq is equal to 1, it indicates that the corresponding element of vx and vy at the position is not a missing value  
                wq(wq<1) = q ;
                Vx = Vx .* wq;  Vy = Vy .* wq;
                [Correlation_r,Correlation_p] = corr(Vx',Vy','type','Pearson');
                if Correlation_p < 0.05 
                    Correlation_M(i,j) = Correlation_r;
                end
            end
        end
        for i = 2: dm
            for j = 1: (i-1)
                Correlation_M(i,j) = Correlation_M(j,i) ;
            end
        end
 
        C_M = Correlation_M ; 
        [ Yzu ] = fenzu(C_M, ZUnum) ; 
        YcellC = cell(ZUnum,1); 
        LenZu = zeros(1,ZUnum); 
        for k = 1: ZUnum
            YcellC{k,1}=find(Yzu==k);%Marking of rows
            LenZu(k) = length(YcellC{k,1});% the numbers of the row of each group
        end
 
        X2 = X1 ;
        [maxZhi,maxP] = max(LenZu);     
        [ X2max ] = maxZu(YcellC{maxP,1}, X1, P,delta,ddcishu,tol) ; % the second filling only for the maximum block          
        X2(YcellC{maxP,1},:) = X2max ; % Embed the maximum block of secondary filling into X2    
        %Measurement indicators of the second filling,RMSE1,pearson(spearman)
        %correlation coefficient,coefficient of determination
       for j=1:dn
            Vx = X2(:,j);
            Vy = M(:,j);
            rmse22(qcshu,j) = norm( (Vx-Vy),'fro' )/ sqrt(dm) ;  
        end
        for j=1: dn
            Vx = X2(:,j);
            Vy = M(:,j);
            [Co_r,Co_p] = corr(Vx,Vy,'type','Pearson');
            P22(qcshu,j) = Co_r;    
        end   
        for j=1: dn
            Vx = X2(:,j);
            Vy = M(:,j);
            [Co_r,Co_p] = corr(Vx,Vy,'type','spearman');
            s22(qcshu,j) = Co_r;    
        end
        for j=1:dn
           Vx = X2(:,j);
           Vy = M(:,j);Vy_m=mean(Vy');
           r22(qcshu,j) =1-(norm((Vy-Vx),'fro')/ norm((Vy-Vy_m),'fro'))^2;  
        end
        Xjilu{qcshu,1} = X2 ;  
    end
    TC2flag = 0 ;     
end
S2=mean(s22,2);P2=mean(P22,2);R2=mean(r22,2);RMSE2=mean(rmse22,2);
