function [tags,SD1,center] = kmean(data,K,iter)
datalength = length(data);

%%���ѡ���ʼK������ֵ
center = zeros(K);
randnumber = randperm(datalength,K);
for i=1:K
    center(i) = data(randnumber(1,i));
end
%%����
for i=1:iter
    distance = zeros(K,datalength);
    tags = zeros(datalength,K);
    for j = 1:datalength
        for m = 1:K
            distance(m,j) = abs(data(j)-center(m));%�������
        end
        [minDistance,index]=min(distance(:,j));% Ѱ�Ҿ�����С���������
        tags(j,index)=1;% �����С����������λ�ã����
    end
    center = (data.'*tags).';%���¸�ֵ����ֵ
    for j = 1:K
        center(j) = center(j)/sum(tags(:,j));
    end
end
%%������
DISTANCE = zeros(K,datalength);
for j = 1:datalength
    for m = 1:K
        DISTANCE(m,j) = power(distance(m,j),2);%�������
    end
end
SD = DISTANCE*tags;
for m = 1:K
    SD1(m) = sqrt(SD(m,m)/(sum(tags(:,m))-1));%���㷽��
end
end