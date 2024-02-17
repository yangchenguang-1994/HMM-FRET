function [tags,SD1,center] = kmean(data,K,iter)
datalength = length(data);

%%随机选择初始K个中心值
center = zeros(K);
randnumber = randperm(datalength,K);
for i=1:K
    center(i) = data(randnumber(1,i));
end
%%迭代
for i=1:iter
    distance = zeros(K,datalength);
    tags = zeros(datalength,K);
    for j = 1:datalength
        for m = 1:K
            distance(m,j) = abs(data(j)-center(m));%计算距离
        end
        [minDistance,index]=min(distance(:,j));% 寻找距离最小的类别索引
        tags(j,index)=1;% 标记最小距离所处的位置（类别）
    end
    center = (data.'*tags).';%重新赋值中心值
    for j = 1:K
        center(j) = center(j)/sum(tags(:,j));
    end
end
%%输出结果
DISTANCE = zeros(K,datalength);
for j = 1:datalength
    for m = 1:K
        DISTANCE(m,j) = power(distance(m,j),2);%计算距离
    end
end
SD = DISTANCE*tags;
for m = 1:K
    SD1(m) = sqrt(SD(m,m)/(sum(tags(:,m))-1));%计算方差
end
end