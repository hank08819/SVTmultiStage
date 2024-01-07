function [ Yzu ] = fenzu(simM, ZUnum)
% 根据相似性 分组/聚类
% simM 是一个已知相似性矩阵，下面先要将其转化为距离矩阵，
%并形成 linkage 函数指定的格式
% ZUnum 指定组/类别 数
% Y=pdist(X1);
    Y = squareform( 1 - simM ) ;
    YZ=linkage(Y,'average');
%     dendrogram(YZ);
    YT=cluster(YZ,'maxclust',ZUnum);
    Yzu = YT';