function  [Number,Sta]=Select_Infected_Point(M,Sta,De)
%M��ά����Sta���û�״̬��De����
%0:�׸�״̬S��Susceptible��  P_0_1; ��P_0_3:Ԥ����ϵ����
%1:Ǳ��״̬E��Exposed��      P_1_0��  P_1_2��P_1_3
%2:Ⱦ��״̬I��Infected��     P_2_0��  P_2_3
%3:����״̬R��Recovered��    P_3_0
%ȡ�����Ľڵ�Ϊ��Դ�ڵ㣬��������ߣ���ѡ�δ�ģ�һ����ȥ
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
        
    