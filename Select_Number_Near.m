function Num=Select_Number_Near(j,Link)         %找出节点j的邻居节点
p=0;
for i=1:length(Link)  %length：行数与列数的较大值 size：获取数组的行数与列数 numel：元素总数
    if Link(i,j)==1
        p=p+1;
        %Num(p)=i;
    else
        continue;
    end
end
Num = zeros(1,p);
k=0;
for m=1:length(Link)  %length：行数与列数的较大值 size：获取数组的行数与列数 numel：元素总数
    if Link(m,j)==1
        k=k+1;
        Num(k)=m;
    else
        continue;
    end
end