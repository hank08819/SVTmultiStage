function [ Yzu ] = fenzu(simM, ZUnum)
% ���������� ����/����
% simM ��һ����֪�����Ծ���������Ҫ����ת��Ϊ�������
%���γ� linkage ����ָ���ĸ�ʽ
% ZUnum ָ����/��� ��
% Y=pdist(X1);
    Y = squareform( 1 - simM ) ;
    YZ=linkage(Y,'average');
%     dendrogram(YZ);
    YT=cluster(YZ,'maxclust',ZUnum);
    Yzu = YT';