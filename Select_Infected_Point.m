function  [Number,Sta]=Select_Infected_Point(M,Sta,De)
%M：维数；Sta：用户状态；De：度
%0:易感状态S（Susceptible）  P_0_1; （P_0_3:预免疫系数）
%1:潜伏状态E（Exposed）      P_1_0；  P_1_2；P_1_3
%2:染病状态I（Infected）     P_2_0；  P_2_3
%3:免疫状态R（Recovered）    P_3_0
%取度最大的节点为病源节点，如果已免疫，则选次大的，一次下去
flag=1;
while flag
    Number=max(find(max(De)==De));
    if Sta(Number)==3
        De(Number)=0;
        continue;
    else
        Sta(Number)=2;
        flag=0;
    end
end
        
    