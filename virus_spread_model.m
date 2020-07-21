function Japan4

tic %���浱ǰʱ��
%load 'C:\Users\lqs\Desktop\�¹ڷ���\cellular-automata\cellular-automata\Simulation of short message network virus propagation based on Cellular Automata\Data\Link.txt';    %�����ڽӾ���
%-------------------------------------------------------------------------%
%״̬�ֲ���״̬ת�Ƹ���SEIR
%0:�׸�״̬δ����S��Susceptible��  P_0_1; ��P_0_3:Ԥ����ϵ����
%1:Ǳ��״̬δ����E��Exposed��      P_1_0;  P_1_2��P_1_3
%2:Ⱦ��״̬δ����I��Infected��     P_2_0��  P_2_3
%3:����״̬R��Recovered��    P_3_0
%4:�׸б�����SQ
%5:Ǳ���ڱ�����EQ
%6:����������H
%7:����D
%-------------------------------------------------------------------------%
M=10000;
Link=zeros(M,M);
%������û��ڵ�Ķ�
%De=sum(Link);                                                              %�û��ڵ�Ķȣ��Ծ��������
%------------������������----����������˵��--------------------------------%

%������Ⱦ������Ϊ��ͬ
I_E=0.05;                                                                   %Ǳ����E�û��Ĵ�Ⱦǿ��
I_I=0.05;                                                                   %������I�û��Ĵ�Ⱦǿ��

%�����Ӵ��������ڽӾ���ȷ��
%lamda=sum(De)/M;                                                           %ƽ���Ӵ�����
%State���ֻ��û�����״̬State=zeros(1,M); 0:��ʾ�׸�״̬��Susceptible��
%---------------------------------1---------------------------------------%

%parameters%

q1 = 0;         %�׸���Ⱥ�Ӵ����ߣ����������
q2 = 0;        %�������������
contactNum=4;   %ƽ���Ӵ�����
contactNum2=60; %�ۻ�Ӵ�
Init = 100;  %���ó�ʼ��Ⱦ����
death_rate = 0.1; %������
input = 3; %�����Ͳ���
file_name = '3.csv';

%----------%
TimeStep=100;
%P_m1=0.1;         %�ֿ���������
P_m2 = 0.005;     %Ǳ���������ĸ���
P_m3 = 0.002;     %�������ٻ����ĸ���

E_threshold = 14;   % Ǳ�������ʱ��
touch = 4;

TimeLong_F=zeros(1,M); %�û�����Ⱦ���ڵ�ʱ�䳤��
TimeLong_E=zeros(1,M); %�û�����Ǳ���ڵ�ʱ�䳤��
qTime=zeros(1,M); %����ʱ��
Sta=zeros(1,M);        %״̬
%������Ԥ�����趨

%״̬ת��
%��ʼ���ѡ��һ���ڵ�Ϊ��Դ��
%ȡ�����Ľڵ�Ϊ��Դ�ڵ�
%[Number,Sta]=Select_Infected_Point(M,Sta,De);
%Number����Դ�ڵ�
%State ��ȷ����Դ�ڵ��Ժ�Ľڵ�״̬����
State=zeros(TimeStep,M); %ģ��ʱ���ڵ�״̬����ÿһ�д���ǰʱ��Ƭ��ÿ���˵�״̬
Init_E=randperm(M,Init*0.7); % ��ʼ��Ǳ����
Init_I=randperm(M,Init*0.3); % ��ʼ��������
%��ʼ����Ⱦ��Ⱥ״̬
Sta(Init_E)=1;
Sta(Init_I)=2;
%���ݳ�ʼ��Դ���ʼLink�ڽӾ���
flag = 0;
for x=1:M
    if Sta(x)==1
        temprow = randperm(M,contactNum);  %����Ǳ����ÿ������Ӵ�50��,���Ǽ�ӽӴ�
        Link(x,temprow)=1;
        Link(temprow,x)=1;
    end
    if Sta(x)==2
        temprow = randperm(M,contactNum);  %�����Ⱦ��ÿ������Ӵ�50��
        Link(x,temprow)=1;
        Link(temprow,x)=1;
    end
end


counter = zeros(M,M);  % ����������ʾ��֮����һ�νӴ��������
index = find(Link==1);
% ��ʼ��������
%counter(index) = counter(index)+1;

Number_State=zeros(8,TimeStep);  %�û�����4��״̬��ͳ������

for t=1:TimeStep
    if t==1
        State(t,:)=Sta; % ���Ƴ�ʼ״̬
    else
        % ���¼����
        De=sum(Link);
        lamda=sum(De)/M;
        
        % ����������
        for j=1:M
            %�ж��û��ڵ㴦��ʲô״̬��Ȼ�������״̬ȷ����ת�����
            if State(t-1,j)==0                          %��ʱ�����׸�״̬0��������Ǳ����ת��
                Num=Select_Number_Near(j,Link);         %�ҳ��ڵ�j���ھӽڵ�,Num��һ������
                P=zeros(1,length(Num));                 %ÿ���ھӽڵ��Ⱦ�ýڵ�ĸ���
                for k=1:length(Num)
                    if State(t-1,Num(k))==1             %���ھӽڵ㴦��Ǳ����E��1��
                        %if TimeLong_E(Num(k)) >= counter(Num(k),j)  % ��Ǳ�����ڼ�Ӵ���
                        %P(k)=I_E/De(Num(k))*sum((lamda.^(1:De(Num(k))).*exp(-lamda))./(factorial(1:De(Num(k))-1)));
                        P(k)=I_E;
                        if P(k) >=rand                %����Ⱦ���ʴ��ڵֿ�����
                            State(t,j)=1;      % ת��ΪǱ��״̬
                            TimeLong_E(j)=1;
                        else
                            State(t,j)=State(t-1,j);        %ά��ԭ��״̬
                        end
                        %end
                    else
                        if State(t-1,Num(k))==2          %�ڵ㴦��Ⱦ����I��2��
                            %P(k)=I_I/De(Num(k))*sum((lamda.^(1:De(Num(k))).*exp(-lamda))./(factorial(1:De(Num(k))-1)));
                            P(k)=I_I;
                            if P(k) >=rand                %����Ⱦ���ʴ��ڵֿ�����
                                State(t,j)=1;      % ת��ΪǱ��״̬
                                TimeLong_E(j)=1;
                            else
                                State(t,j)=State(t-1,j);        %ά��ԭ��״̬
                            end
                        else
                            continue;
                        end
                    end
                end
                
            else
                % �ⲿ�ֵ��߼��޸Ĺ���
                Num=Select_Number_Near(j,Link);
                if State(t-1,j)==1         %��ʱ����Ǳ��״̬E��������Ⱦ��I������Rת��
                    if rand<=exp((TimeLong_E(j)-E_threshold)/TimeLong_E(j))     %��Ⱦ��״̬Iת��
                        if rand<q2
                            State(t,j)=6;
                            qTime(j)=qTime(j)+1;
                            TimeLong_E(j)=0;
                            %��Ⱦ��״̬ת�ƣ������Ӵ��Ĳ�����Ⱥ�ᱻ����
                            for m=1:length(Num)
                                if State(t-1,Num(m))==0             %���ھӽڵ㴦���׸У�1��
                                    if counter(Num(m),j)<14 && rand < q1 % ���Ӵ�ʱ������14��
                                        State(t-1,Num(m))=4; % ״̬ת��Ϊ�׸б�����
                                        qTime(Num(m))=qTime(Num(m))+1; %����ʱ���1��
                                    end
                                else
                                    if State(t-1,Num(m))==1          %�ڵ㴦��Ǳ��δ����״̬��2��
                                        if counter(Num(m),j)<14 && rand < q1  % ���Ӵ�ʱ������14��
                                            State(t-1,Num(m))=5;  % ״̬ת��ΪǱ��������
                                            qTime(Num(m))=qTime(Num(m))+1; %����ʱ���1��
                                        end
                                    else
                                        if State(t-1,Num(m))==2  %�ڵ㴦�ڻ�����δ����
                                            if counter(Num(m),j)<14 && rand < q2 % ���Ӵ�ʱ������14��
                                                State(t-1,Num(m))=6; % ״̬ת��Ϊ����������
                                                qTime(Num(m))=qTime(Num(m))+1; %����ʱ���1��
                                            end
                                        end
                                    end
                                end
                            end
                        else
                            State(t,j)=2;
                            TimeLong_F(j)=TimeLong_F(j)+1;         %�û�j����Ⱦ��״̬��ʱ�䳤��
                            TimeLong_E(j)=0;
                        end
                    else
                        if rand<=P_m2             %������״̬Rת��
                            State(t,j)=3;
                        else
                            State(t,j)=State(t-1,j); %����ԭ״̬
                            TimeLong_E(j)=TimeLong_E(j)+1;
                        end
                    end
                else
                    if State(t-1,j)==2        %��ʱ����ȷ��״̬I������������Rת��
                        if rand > 1/(1+exp((10-TimeLong_F(j))/2))         % (1+exp((10-TimeLong_F(j))/2))��ʾ�����ʣ����Ż���������Խ��Խ��
                            if rand < q2
                                State(t,j)=5;
                                qTime(j)=qTime(j)+1;
                                TimeLong_F(j)=TimeLong_F(j)+1;
                                %����״̬�����룬�����Ӵ��Ĳ�����Ⱥ�ᱻ����
                                for m=1:length(Num)
                                    if State(t-1,Num(m))==0             %���ھӽڵ㴦���׸У�1��
                                        if counter(Num(m),j)<14 && rand < q1 % ���Ӵ�ʱ������14��
                                            State(t-1,Num(m))=4; % ״̬ת��Ϊ�׸б�����
                                            qTime(Num(m))=qTime(Num(m))+1; %����ʱ���1��
                                        end
                                    else
                                        if State(t-1,Num(m))==1          %�ڵ㴦��Ǳ��δ����״̬��2��
                                            if counter(Num(m),j)<14 && rand < q1  % ���Ӵ�ʱ������14��
                                                State(t-1,Num(m))=5;  % ״̬ת��ΪǱ��������
                                                qTime(Num(m))=qTime(Num(m))+1; %����ʱ���1��
                                            end
                                        else
                                            if State(t-1,Num(m))==2  %�ڵ㴦�ڻ�����δ����
                                                if counter(Num(m),j)<14 && rand < q2 % ���Ӵ�ʱ������14��
                                                    State(t-1,Num(m))=6; % ״̬ת��ΪǱ��������
                                                    qTime(Num(m))=qTime(Num(m))+1; %����ʱ���1��
                                                end
                                            end
                                        end
                                    end
                                end
                            else
                                State(t,j)=State(t-1,j);           %ά�ֻ���״̬
                                TimeLong_F(j)=TimeLong_F(j)+1;
                            end
                        else
                            if rand > death_rate
                                State(t,j)=3;
                                TimeLong_F(j)=0; %���ڸ�Ⱦ�ڣ��ж�״̬����ʱ�䳤��
                            else
                                State(t,j)=7;
                                TimeLong_F(j)=0;
                            end
                        end
                    else
                        %����������,��С�ĸ��ʻ��ٴ�Ⱦ��
                        if rand <= P_m3
                            State(t,j)=2;
                            TimeLong_F(j)=TimeLong_F(j)+1;
                        else
                            State(t,j)=State(t-1,j);
                        end
                    end
                end
            end
            %�׸��ѱ������״̬ת��
            if State(t-1,j)==4
                if qTime(j)>14
                    State(t-1,j)=0;
                else
                    qTime(j)=qTime(j)+1;
                end
            end
            %Ǳ���ѱ������״̬ת��
            if State(t-1,j)==5
                if rand<=exp((TimeLong_E(j)-E_threshold)/TimeLong_E(j)) %��Ⱦ��״̬ת�ƣ��ر����룬����ʱ��+1��Ǳ��ʱ����0��Ⱦ��ʱ��+1
                    State(t-1,j)=6;
                    qTime(j) = qTime(j)+1;
                    TimeLong_F(j)=TimeLong_F(j)+1;
                    TimeLong_E(j)=0;
                else    %����ԭ״̬
                    qTime(j) = qTime(j)+1;
                    TimeLong_E(j)=TimeLong_E(j)+1;
                end
            end
            %ȷ���ѱ������״̬ת��
            if State(t-1,j)==6
                if rand > 1/(1+exp((10-TimeLong_F(j))/2))         % (1+exp((10-TimeLong_F(j))/2))��ʾ�����ʣ����Ż���������Խ��Խ��
                    State(t,j)=State(t-1,j);           %ά�ֻ���״̬
                    TimeLong_F(j)=TimeLong_F(j)+1;
                    qTime(j)=qTime(j)+1;
                else
                    if rand > death_rate    % ������
                        State(t,j)=3;
                        TimeLong_F(j)=0; %���ڸ�Ⱦ�ڵ�ʱ�䳤��
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
        % ���¼��������ڽӾ���
        index2 = find(counter>0);
        counter(tempIndex) = 1;
        counter(index2) = counter(index2)+1;
        Link(tempIndex)=0;
        %
        %touch = (Number_State(2,t-1)+Number_State(1,t-1))*10; %�Ӵ���������Ǳ���ںͻ�����δ����������й�,����ÿ���˽Ӵ�10��
        
       % ÿһ�ܽ���һ���ڽ̼���
        if mod(t,7)==0
            newTouchNum = round(contactNum2*sum(State(t,:)==0)/M);
        else
            newTouchNum = contactNum;
        end
        %�����ڽӾ���
        for x=1:M
            if State(t-1,x)==1
                temprow = randperm(M,newTouchNum);  %����Ǳ����ÿ������Ӵ�10��
                Link(x,temprow)=1;
                Link(temprow,x)=1;
            end
            if State(t-1,x)==2
                temprow = randperm(M,newTouchNum);  %�����Ⱦ��ÿ������Ӵ�10��
                Link(x,temprow)=1;
                Link(temprow,x)=1;
            end
        end
        for x=1:M
            if qTime(x)>0  %���账�ڱ�����״̬���޽Ӵ�
                Link(x,:)=0;
                Link(:,x)=0;
            end
        end
        % �������벡���߼�
        % ״̬ΪǱ��״̬δ����
        % ʵ��˼·��������ԭ����ά�ȣ�ѡ���׸�״̬ת�����
%         tempInputNum = round(rand*2); %�����Ͳ�������
%         potentialIndex = find(State(t,:)==0); %�׸���Ⱥ
%         inputIndex = potentialIndex(randperm(numel(potentialIndex),tempInputNum)); %��Ϊ�����͵��׸���Ⱥ
%         State(t,inputIndex)= 1;
%         TimeLong_E(inputIndex)=3;  %Ǳ������3
    end
    
    %ͳ�Ƹ�״̬�Ľڵ�����
    Number_State(1,t)=sum(State(t,:)==0);%�����׸�״̬S���ܽڵ�����
    Number_State(2,t)=sum(State(t,:)==1);%����Ǳ��δ����E���ܽڵ�����
    Number_State(3,t)=sum(State(t,:)==2);%���ڻ���δ����״̬I���ܽڵ�����
    Number_State(4,t)=sum(State(t,:)==3);%��������״̬R���ܽڵ�����
    Number_State(5,t)=sum(State(t,:)==4);%�����׸б�����״̬SQ���ܽڵ�����
    Number_State(6,t)=sum(State(t,:)==5);%����Ǳ��������״̬EQ���ܽڵ�����
    Number_State(7,t)=sum(State(t,:)==6);%���ڻ���������״̬H���ܽڵ�����
    Number_State(8,t)=sum(State(t,:)==7);%��������״̬D���ܽڵ�����
    Number_State(9,t)=sum(State(t,:)==2)+sum(State(t,:)==6);%���ڻ������ܽڵ�����
    
    % �����ȶ��󣬳���������
    if t>30
        if Number_State(6,t) <60
            flag = flag+1;
            tempInputNum = round(rand*input); %�����Ͳ�������
            potentialIndex = find(State(t,:)==0); %�׸���Ⱥ
            inputIndex = potentialIndex(randperm(numel(potentialIndex),tempInputNum)); %��Ϊ�����͵��׸���Ⱥ
            State(t,inputIndex)= 1;
            TimeLong_E(inputIndex)=3;  %Ǳ������3
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
    legend('�׸�δ����S','Ǳ��״̬E','Ⱦ��״̬I','����״̬R','�׸б�����SQ','Ǳ���ڱ�����EQ','ȷ�ﱻ����H','����D','����')
    xlabel('ģ��ʱ��')
    ylabel('��״̬������')
end

csvwrite(file_name,Number_State);