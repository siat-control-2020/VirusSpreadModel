function Num=Select_Number_Near(j,Link)         %�ҳ��ڵ�j���ھӽڵ�
p=0;
for i=1:length(Link)  %length�������������Ľϴ�ֵ size����ȡ��������������� numel��Ԫ������
    if Link(i,j)==1
        p=p+1;
        %Num(p)=i;
    else
        continue;
    end
end
Num = zeros(1,p);
k=0;
for m=1:length(Link)  %length�������������Ľϴ�ֵ size����ȡ��������������� numel��Ԫ������
    if Link(m,j)==1
        k=k+1;
        Num(k)=m;
    else
        continue;
    end
end