function Japan4

tic %保存当前时间
%load 'C:\Users\lqs\Desktop\新冠肺炎\cellular-automata\cellular-automata\Simulation of short message network virus propagation based on Cellular Automata\Data\Link.txt';    %读入邻接矩阵
%-------------------------------------------------------------------------%
%状态分布及状态转移概率SEIR
%0:易感状态未隔离S（Susceptible）  P_0_1; （P_0_3:预免疫系数）
%1:潜伏状态未隔离E（Exposed）      P_1_0;  P_1_2；P_1_3
%2:染病状态未隔离I（Infected）     P_2_0；  P_2_3
%3:免疫状态R（Recovered）    P_3_0
%4:易感被隔离SQ
%5:潜伏期被隔离EQ
%6:患病被隔离H
%7:死亡D
%-------------------------------------------------------------------------%
M=10000;
Link=zeros(M,M);
%计算各用户节点的度
%De=sum(Link);                                                              %用户节点的度，对矩阵按列求和
%------------——————----参数设置与说明--------------------------------%

%《《感染率设置为相同
I_E=0.05;                                                                   %潜伏期E用户的传染强度
I_I=0.05;                                                                   %发病期I用户的传染强度

%《《接触人数由邻接矩阵确定
%lamda=sum(De)/M;                                                           %平均接触人数
%State：手机用户所处状态State=zeros(1,M); 0:表示易感状态（Susceptible）
%---------------------------------1---------------------------------------%

%parameters%

q1 = 0;         %易感人群接触患者，被隔离概率
q2 = 0;        %患病被隔离概率
contactNum=4;   %平常接触人数
contactNum2=60; %聚会接触
Init = 100;  %设置初始感染人数
death_rate = 0.1; %死亡率
input = 3; %输入型病例
file_name = '3.csv';

%----------%
TimeStep=100;
%P_m1=0.1;         %抵抗病毒概率
P_m2 = 0.005;     %潜伏期自愈的概率
P_m3 = 0.002;     %治愈后再患病的概率

E_threshold = 14;   % 潜伏期最大时限
touch = 4;

TimeLong_F=zeros(1,M); %用户处于染病期的时间长短
TimeLong_E=zeros(1,M); %用户处于潜伏期的时间长短
qTime=zeros(1,M); %隔离时间
Sta=zeros(1,M);        %状态
%不进行预免疫设定

%状态转换
%初始随机选择一个节点为病源点
%取度最大的节点为病源节点
%[Number,Sta]=Select_Infected_Point(M,Sta,De);
%Number：病源节点
%State ：确定病源节点以后的节点状态矩阵
State=zeros(TimeStep,M); %模拟时间内的状态矩阵，每一行代表当前时间片内每个人的状态
Init_E=randperm(M,Init*0.7); % 初始化潜伏者
Init_I=randperm(M,Init*0.3); % 初始化患病者
%初始化感染人群状态
Sta(Init_E)=1;
Sta(Init_I)=2;
%根据初始病源点初始Link邻接矩阵
flag = 0;
for x=1:M
    if Sta(x)==1
        temprow = randperm(M,contactNum);  %假设潜伏者每天随机接触50人,考虑间接接触
        Link(x,temprow)=1;
        Link(temprow,x)=1;
    end
    if Sta(x)==2
        temprow = randperm(M,contactNum);  %假设感染者每天随机接触50人
        Link(x,temprow)=1;
        Link(temprow,x)=1;
    end
end


counter = zeros(M,M);  % 计数器，表示人之间上一次接触后的天数
index = find(Link==1);
% 初始化计数器
%counter(index) = counter(index)+1;

Number_State=zeros(8,TimeStep);  %用户处于4个状态的统计数量

for t=1:TimeStep
    if t==1
        State(t,:)=Sta; % 复制初始状态
    else
        % 重新计算度
        De=sum(Link);
        lamda=sum(De)/M;
        
        % 遍历所有人
        for j=1:M
            %判断用户节点处于什么状态，然后根据其状态确定其转变情况
            if State(t-1,j)==0                          %此时处于易感状态0，可能向潜伏期转移
                Num=Select_Number_Near(j,Link);         %找出节点j的邻居节点,Num是一个数组
                P=zeros(1,length(Num));                 %每个邻居节点感染该节点的概率
                for k=1:length(Num)
                    if State(t-1,Num(k))==1             %若邻居节点处于潜伏期E（1）
                        %if TimeLong_E(Num(k)) >= counter(Num(k),j)  % 若潜伏期期间接触了
                        %P(k)=I_E/De(Num(k))*sum((lamda.^(1:De(Num(k))).*exp(-lamda))./(factorial(1:De(Num(k))-1)));
                        P(k)=I_E;
                        if P(k) >=rand                %若感染概率大于抵抗概率
                            State(t,j)=1;      % 转变为潜伏状态
                            TimeLong_E(j)=1;
                        else
                            State(t,j)=State(t-1,j);        %维持原来状态
                        end
                        %end
                    else
                        if State(t-1,Num(k))==2          %节点处于染病期I（2）
                            %P(k)=I_I/De(Num(k))*sum((lamda.^(1:De(Num(k))).*exp(-lamda))./(factorial(1:De(Num(k))-1)));
                            P(k)=I_I;
                            if P(k) >=rand                %若感染概率大于抵抗概率
                                State(t,j)=1;      % 转变为潜伏状态
                                TimeLong_E(j)=1;
                            else
                                State(t,j)=State(t-1,j);        %维持原来状态
                            end
                        else
                            continue;
                        end
                    end
                end
                
            else
                % 这部分的逻辑修改过了
                Num=Select_Number_Near(j,Link);
                if State(t-1,j)==1         %此时处于潜伏状态E，可能向染病I和免疫R转移
                    if rand<=exp((TimeLong_E(j)-E_threshold)/TimeLong_E(j))     %向染病状态I转移
                        if rand<q2
                            State(t,j)=6;
                            qTime(j)=qTime(j)+1;
                            TimeLong_E(j)=0;
                            %向染病状态转移，曾经接触的部分人群会被隔离
                            for m=1:length(Num)
                                if State(t-1,Num(m))==0             %若邻居节点处于易感（1）
                                    if counter(Num(m),j)<14 && rand < q1 % 若接触时间少于14天
                                        State(t-1,Num(m))=4; % 状态转移为易感被隔离
                                        qTime(Num(m))=qTime(Num(m))+1; %隔离时间加1天
                                    end
                                else
                                    if State(t-1,Num(m))==1          %节点处于潜伏未隔离状态（2）
                                        if counter(Num(m),j)<14 && rand < q1  % 若接触时间少于14天
                                            State(t-1,Num(m))=5;  % 状态转移为潜伏被隔离
                                            qTime(Num(m))=qTime(Num(m))+1; %隔离时间加1天
                                        end
                                    else
                                        if State(t-1,Num(m))==2  %节点处于患病但未隔离
                                            if counter(Num(m),j)<14 && rand < q2 % 若接触时间少于14天
                                                State(t-1,Num(m))=6; % 状态转移为患病被隔离
                                                qTime(Num(m))=qTime(Num(m))+1; %隔离时间加1天
                                            end
                                        end
                                    end
                                end
                            end
                        else
                            State(t,j)=2;
                            TimeLong_F(j)=TimeLong_F(j)+1;         %用户j处于染病状态的时间长短
                            TimeLong_E(j)=0;
                        end
                    else
                        if rand<=P_m2             %向免疫状态R转移
                            State(t,j)=3;
                        else
                            State(t,j)=State(t-1,j); %保持原状态
                            TimeLong_E(j)=TimeLong_E(j)+1;
                        end
                    end
                else
                    if State(t-1,j)==2        %此时处于确诊状态I，可以向免疫R转移
                        if rand > 1/(1+exp((10-TimeLong_F(j))/2))         % (1+exp((10-TimeLong_F(j))/2))表示治愈率，随着患病天数会越来越大
                            if rand < q2
                                State(t,j)=5;
                                qTime(j)=qTime(j)+1;
                                TimeLong_F(j)=TimeLong_F(j)+1;
                                %患病状态被隔离，曾经接触的部分人群会被隔离
                                for m=1:length(Num)
                                    if State(t-1,Num(m))==0             %若邻居节点处于易感（1）
                                        if counter(Num(m),j)<14 && rand < q1 % 若接触时间少于14天
                                            State(t-1,Num(m))=4; % 状态转移为易感被隔离
                                            qTime(Num(m))=qTime(Num(m))+1; %隔离时间加1天
                                        end
                                    else
                                        if State(t-1,Num(m))==1          %节点处于潜伏未隔离状态（2）
                                            if counter(Num(m),j)<14 && rand < q1  % 若接触时间少于14天
                                                State(t-1,Num(m))=5;  % 状态转移为潜伏被隔离
                                                qTime(Num(m))=qTime(Num(m))+1; %隔离时间加1天
                                            end
                                        else
                                            if State(t-1,Num(m))==2  %节点处于患病但未隔离
                                                if counter(Num(m),j)<14 && rand < q2 % 若接触时间少于14天
                                                    State(t-1,Num(m))=6; % 状态转移为潜伏被隔离
                                                    qTime(Num(m))=qTime(Num(m))+1; %隔离时间加1天
                                                end
                                            end
                                        end
                                    end
                                end
                            else
                                State(t,j)=State(t-1,j);           %维持患病状态
                                TimeLong_F(j)=TimeLong_F(j)+1;
                            end
                        else
                            if rand > death_rate
                                State(t,j)=3;
                                TimeLong_F(j)=0; %处于感染期（中毒状态）的时间长度
                            else
                                State(t,j)=7;
                                TimeLong_F(j)=0;
                            end
                        end
                    else
                        %处于免疫期,很小的概率会再次染病
                        if rand <= P_m3
                            State(t,j)=2;
                            TimeLong_F(j)=TimeLong_F(j)+1;
                        else
                            State(t,j)=State(t-1,j);
                        end
                    end
                end
            end
            %易感已被隔离的状态转移
            if State(t-1,j)==4
                if qTime(j)>14
                    State(t-1,j)=0;
                else
                    qTime(j)=qTime(j)+1;
                end
            end
            %潜伏已被隔离的状态转移
            if State(t-1,j)==5
                if rand<=exp((TimeLong_E(j)-E_threshold)/TimeLong_E(j)) %向染病状态转移，必被隔离，隔离时间+1，潜伏时间置0，染病时间+1
                    State(t-1,j)=6;
                    qTime(j) = qTime(j)+1;
                    TimeLong_F(j)=TimeLong_F(j)+1;
                    TimeLong_E(j)=0;
                else    %保持原状态
                    qTime(j) = qTime(j)+1;
                    TimeLong_E(j)=TimeLong_E(j)+1;
                end
            end
            %确诊已被隔离的状态转移
            if State(t-1,j)==6
                if rand > 1/(1+exp((10-TimeLong_F(j))/2))         % (1+exp((10-TimeLong_F(j))/2))表示治愈率，随着患病天数会越来越大
                    State(t,j)=State(t-1,j);           %维持患病状态
                    TimeLong_F(j)=TimeLong_F(j)+1;
                    qTime(j)=qTime(j)+1;
                else
                    if rand > death_rate    % 死亡率
                        State(t,j)=3;
                        TimeLong_F(j)=0; %处于感染期的时间长度
                    else
                        State(t,j)=7;
                        TimeLong_F(j)=0;
                        qTime(j)=0;
                    end
                end
            end
        end
        
        tempIndex = find(Link==1);
        %counter(tempIndex)=0;
        % 更新计数器和邻接矩阵
        index2 = find(counter>0);
        counter(tempIndex) = 1;
        counter(index2) = counter(index2)+1;
        Link(tempIndex)=0;
        %
        %touch = (Number_State(2,t-1)+Number_State(1,t-1))*10; %接触次数，与潜伏期和患病但未隔离的数量有关,假设每个人接触10人
        
       % 每一周进行一次宗教集会
        if mod(t,7)==0
            newTouchNum = round(contactNum2*sum(State(t,:)==0)/M);
        else
            newTouchNum = contactNum;
        end
        %更新邻接矩阵
        for x=1:M
            if State(t-1,x)==1
                temprow = randperm(M,newTouchNum);  %假设潜伏者每天随机接触10人
                Link(x,temprow)=1;
                Link(temprow,x)=1;
            end
            if State(t-1,x)==2
                temprow = randperm(M,newTouchNum);  %假设感染者每天随机接触10人
                Link(x,temprow)=1;
                Link(temprow,x)=1;
            end
        end
        for x=1:M
            if qTime(x)>0  %假设处于被隔离状态则无接触
                Link(x,:)=0;
                Link(:,x)=0;
            end
        end
        % 新增输入病例逻辑
        % 状态为潜伏状态未隔离
        % 实现思路：不增加原矩阵维度，选择将易感状态转变过来
%         tempInputNum = round(rand*2); %输入型病例数量
%         potentialIndex = find(State(t,:)==0); %易感人群
%         inputIndex = potentialIndex(randperm(numel(potentialIndex),tempInputNum)); %作为输入型的易感人群
%         State(t,inputIndex)= 1;
%         TimeLong_E(inputIndex)=3;  %潜伏期置3
    end
    
    %统计各状态的节点数量
    Number_State(1,t)=sum(State(t,:)==0);%处于易感状态S的总节点数量
    Number_State(2,t)=sum(State(t,:)==1);%处于潜伏未隔离E的总节点数量
    Number_State(3,t)=sum(State(t,:)==2);%处于患病未隔离状态I的总节点数量
    Number_State(4,t)=sum(State(t,:)==3);%处于免疫状态R的总节点数量
    Number_State(5,t)=sum(State(t,:)==4);%处于易感被隔离状态SQ的总节点数量
    Number_State(6,t)=sum(State(t,:)==5);%处于潜伏被隔离状态EQ的总节点数量
    Number_State(7,t)=sum(State(t,:)==6);%处于患病被隔离状态H的总节点数量
    Number_State(8,t)=sum(State(t,:)==7);%处于死亡状态D的总节点数量
    Number_State(9,t)=sum(State(t,:)==2)+sum(State(t,:)==6);%处于患病的总节点数量
    
    % 疫情稳定后，出现输入型
    if t>30
        if Number_State(6,t) <60
            flag = flag+1;
            tempInputNum = round(rand*input); %输入型病例数量
            potentialIndex = find(State(t,:)==0); %易感人群
            inputIndex = potentialIndex(randperm(numel(potentialIndex),tempInputNum)); %作为输入型的易感人群
            State(t,inputIndex)= 1;
            TimeLong_E(inputIndex)=3;  %潜伏期置3
        end
        if flag>20
            q1=0.0;
            q2=0.0;
        end
    end
    %set(gca,'YTick',0:10:100);
    figure(1)
    if rem(t,3)==0
        plot([t-1 t],[Number_State(1,t-1) Number_State(1,t)],'md-'),hold on
        plot([t-1 t],[Number_State(2,t-1) Number_State(2,t)],'gh:'),hold on
        plot([t-1 t],[Number_State(3,t-1) Number_State(3,t)],'bs-.'),hold on
        plot([t-1 t],[Number_State(4,t-1) Number_State(4,t)],'k.-'),hold on
        plot([t-1 t],[Number_State(5,t-1) Number_State(5,t)],'rx-'),hold on
        plot([t-1 t],[Number_State(6,t-1) Number_State(6,t)],'y+-'),hold on
        plot([t-1 t],[Number_State(7,t-1) Number_State(7,t)],'cv-'),hold on
        plot([t-1 t],[Number_State(8,t-1) Number_State(8,t)],'kp-'),hold on
        plot([t-1 t],[Number_State(9,t-1) Number_State(9,t)],'k-'),hold on
    else
        continue;
    end
    legend('易感未隔离S','潜伏状态E','染病状态I','免疫状态R','易感被隔离SQ','潜伏期被隔离EQ','确诊被隔离H','死亡D','患病')
    xlabel('模拟时间')
    ylabel('各状态的数量')
end

csvwrite(file_name,Number_State);