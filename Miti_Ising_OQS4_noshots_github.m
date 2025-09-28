%Error mitigation of open quantum systems for 5-qubit transverse-field Ising model with ​​periodic boundary conditions​​, where only the first qubit is subject to ​​Amplitude Damping (AD) noise​​.
%The corresponding data file: miti_Ising_OQS4_sp5n500gam1p5rho1O1_NMnoshots

%​​​​Note 1:​​ To rapidly test the program on a personal computer, it is recommended to reduce the key variable 'num_NMt', which controls the total number of iterations, to a value below 5*10^2 from the original setting of 5*10^6 as used in our paper. 
%Note 2: The program was last modified on June 3, 2025. This document was prepared on September 28, 2025.

%Error mitigation of Ising circle for OQS

clear all;close all;clc%
global I X Y Z
I=eye(2); X=[0 1;1 0]; Y=[0 -i;i 0]; Z=[1 0;0 -1];
o=[1;0];l=[0;1];add0=(o+l)/sqrt(2);

n_dim_Ising=1;%定义Ising model维数
num_spin=5; %

for n=1:num_spin %n表示第1行n列处格点
    Z_Pauli{n}=Pauli_product(Z,num_spin,n); %生成广义Z; Pauli_product子程序产生I...A...I张量积
    X_Pauli{n}=Pauli_product(X,num_spin,n); %生成广义X;
end

n_circ=circshift(1:num_spin,1); %目的是为后面设置周期性边界条件而提前建立索引！
for n=1:num_spin
    ZZ{n}=Z_Pauli{n}*Z_Pauli{n_circ(n)};
end

clear n

hx=2;gam=1.5;            %参数初始化
rho0=zeros(2^(num_spin-1));
rho0(1,1)=1;
b1=[0; 1];
b2=[0;1];
rho11=kron(b2*b2',rho0);%布居数观测量
rho0=kron(b1*b1',rho0);%态初始化:|1>|0>^(num_spin-1)
% rho0=kron(l*l',rho0);%态初始化:|1>|0>^(num_spin-1)
rho0_vec=reshape(rho0.',2^(2*num_spin),1); %注意：要注意rho的正确矢量化规则！

Ham=zeros(2^(num_spin)); %Hamiltonian初始化
for n=1:num_spin
    Ham=-0.1*( ZZ{n} )+0.1*hx*X_Pauli{n}+Ham;
    cZZ(n)=-1; %哈密顿量Pauli分解后对应的系数
    cX(n)=hx;

    cHa0(n)=0.1*cZZ(n); %cHα
    cHa0(num_spin+n)=0.1*cX(n);

    PH{n}=ZZ{n};
    PH{num_spin+n}=X_Pauli{n};
end
test_H=0;
for n=1:length(cHa0)
    test_H=test_H+cHa0(n)*PH{n};
end
if all(all(test_H==Ham))==0 %判断系数有无错设。
   warning('哈密顿量分解或初始化错误','on'); %警告提醒
   beep;%警报声音
   errordlg('哈密顿量分解或初始化错误！','输入错误');%弹出警告窗口
   return % %以列为单位（检测矩阵A是否全是逻辑1，是则返1，否则返0）（（any检测矩阵A是否存在部分是逻辑1，是则返1，否则返0））
end

dt=0.02;
L=sqrt(gam)*kron(o*l',eye(2^(num_spin-1)));num_L=1;
cla0=sqrt(gam)/2*[1 1i]; %耗散算子pauli分解系数
Jl_H=kron(o*l',L')+kron(l*o',L);
Jl=expm(-1i*Jl_H*sqrt(dt));

Lin_vec=kron(Ham,eye(2^num_spin))-kron(eye(2^num_spin),Ham.')+...
    1i*(    kron(L,conj(L))-1/2*kron( L'*L,eye(2^num_spin) )-1/2*kron( eye(2^num_spin),conj(L'*L) )   );


trans = @(s) base2dec(strrep(strrep(strrep(strrep(s, 'I', '0'), 'X', '1'), 'Y', '2'), 'Z', '3'), 4)+1;
%该句柄函数可将Pauli串如XXY转换为4进制数112再转为10进制数22+1，之所以要加1，是为了避免matlab中诸如a(0)这样的错误。注意：诸如s=’XYZ’必须加’’号！

K_truc=7; %表A的整体截断（trucation）到的dt^K的项
%--------------------------------------------------------------------
cj=[1 -1/2 -1/2];
C{1}=cj/1;
cl{1}=[1 0 0]; cl{2}=[0 1 0]; cl{3}=[0 0 1];

for k=2:K_truc
    C{k}=kron(C{k-1},cj)/k;
end

B1{1}=(-1)^1/(factorial(2*0)*factorial(2*1-2*0) )*cl{3}+(-1)^1/(factorial(2*1)*factorial(2*1-2*1) )*cl{2};
B2{1}=(-1)^(1-1)/(factorial(2*0+1)*factorial(2*1-2*0-1) )*cl{1};
B{1}=B1{1}+B2{1};
for k=2:K_truc
    for j=0:k

        L2=1;
        m=0;
        while(m<j)
            m=m+1;
            L2=kron(L2,cl{2});%共进行kron运算j次,%[]^j
        end

        L3=1;
        m=0;
        while(m<k-j)
            m=m+1;
            L3=kron(L3,cl{3});%[]^(k-j),共进行kron运算k-j次,%[]^j
        end
        BL23=kron(L2,L3);%临时变量，以免数据编号影响
        B1p{j+1}=BL23/(factorial(2*j)*factorial(2*k-2*j))*(-1)^k; %共k+1项,%B1p类角标表B的左边乘积项

        if j<=k-1
            L3_2=1;
            m=0;
            while(m<k-j-1)
                m=m+1;
                L3_2=kron(L3_2,cl{3});%[]^(k-j-1),共进行kron运算k-j次,%[]^j
            end
            BL123=kron( cl{1},kron(L2,L3_2) );
            B2p{j+1}=BL123/(factorial(2*j+1)*factorial(2*k-2*j-1))*(-1)^(k-1);%共k项,%B2p类角标表B的右边乘积项
        end

    end

    B1s=0;
    for j1=1:k+1
        B1s=B1s+B1p{j1}; %B1s类角标表B的左边求和项
    end

    B2s=0;
    for j1=1:k
        B2s=B2s+B2p{j1};%B2s类角标表B的右边求和项
    end
    B1{k}=B1s;
    B2{k}=B2s;
    B{k}=B1{k}+B2{k};
end

A{1}=C{1}-B{1};
A_abs(1)=sum(abs(A{1}));
B_abs(1)=sum(abs(B{1}));
C_abs(1)=sum(abs(C{1}));
for k=2:K_truc
    AB=0;
    for i0=1:k-1
        AB=AB+kron(A{k-i0},B{i0});
    end
    A{k}=C{k}-B{k}-AB;
    A_abs(k)=sum(abs(A{k}));
    B_abs(k)=sum(abs(B{k}));
    C_abs(k)=sum(abs(C{k}));
end
disp(['|C|^(n)的各阶值是:',num2str(C_abs)])
disp(['|B|^(n)的各阶值是:',num2str(B_abs)])
disp(['|A|^(n)的各阶值是:',num2str(A_abs)])

cla1(1)=sqrt(gam)/2;cla1(2)=1i*sqrt(gam)/2;     %L()L'项
num_cla1=length(cla1);
PL{1}=Pauli_product(X,num_spin,1); %耗散项分解的Pauli算子
PL{2}=Pauli_product(Y,num_spin,1);
p_cla1=abs(cla1)/sum(abs(cla1));

xl_ab{1}=kron(cla1,conj(cla1)); %L()L'项
xl_ab{2}=conj(xl_ab{1});   %L'L()项
xl_ab{3}=conj(xl_ab{1});   %()L'L项

trans2=@(x) dec2base(x,2)-'0'; %将10进制数转成2进制字符串，然后转成2进制字符串数组。-'0'是为了把三进制数组如'101'转成数组[1 0 1];
transA=@(x) dec2base(x,3)-'0'; %将10进制数转成3进制字符串，然后转成三进制字符串数组。-'0'是为了把三进制数组如'101'转成数组[1 0 1];

for n1_3=1:3 %生成xl_ab列矢量，便于编号
    vec_xl_ab{n1_3}=reshape(xl_ab{n1_3}.',num_cla1^2,1);
    xl_ab_abssum(n1_3)=sum(abs(vec_xl_ab{n1_3}));
    p_xl_ab{n1_3}=abs(vec_xl_ab{n1_3})/xl_ab_abssum(n1_3);
end
xl_abssum=sum(xl_ab_abssum);

vec_A_xl_prodt=1;
for k=1:K_truc %张量思想做初始化
    A_xl_pro{k}=zeros((num_cla1^2)^k,length(A{k})); %length(A{k})实际上等于3^k
    posA_xl_pro=trans1tombit_add1(transA(0:length(A{k})-1),k,3); %转为3进制后统一加1的过程已在，相当于转变为了{1,2,3}组成的自定义3进制

    for nn=1:length(A{k})
        vec_x1_pro_tem=1;
        for nnk=1:k
            vec_x1_pro_tem=kron(vec_x1_pro_tem,vec_xl_ab{posA_xl_pro(nn,nnk)});
        end
        A_xl_pro{k}(:,nn)=A_xl_pro{k}(:,nn)+A{k}(nn)*vec_x1_pro_tem; %A_xl_pro{k}是(num_cla1^2)^k×3^k矩阵，即4^k×3^k=12^k个矩阵元！
    end
    A_xl_prodt{k}=A_xl_pro{k}*dt^k;
    vec_A_xl_prodt_tem{k}=reshape(A_xl_prodt{k},1,(num_cla1^2)^k*3^k); %A_xl_prodt列按序堆积在一起，最后转成行矩阵。%12^k个矩阵元
    vec_A_xl_prodt=[vec_A_xl_prodt vec_A_xl_prodt_tem{k}];

    mun(k)=sum(sum(abs(A_xl_pro{k})));
    mundt(k)=sum(sum(abs(A_xl_prodt{k}))); %mundt将会被重新更新为mundt1
end
mun0=1;
mundt0=dt^0;
mundt1=[mundt0 mundt]; %重新更新mundt
disp(['μN的dt的0-K次方阶是: ',num2str(mundt1)])

xl(1)=2*sum(abs(cla1))^2; %xl_ab有3种构型，但Dl=Dl1-Dl2/2-Dl3/3,后两个前应乘以-1/2
xl0=2*sum(abs(cHa0));%H超算子有2种构型
xl=[xl0 xl(1:end)];%xl0表Hamiliton项相关
for q=1:K_truc
    mum(q)=(2*sum(xl))^q/factorial(q);
    mumdt(q)=mum(q)*dt^q;
end
mumdt0=1;
mumdt1=[mumdt0 mumdt];
mumdt1(2)=0; %之所以强制置0，是因为mum得dt^1项实际可以抵消而不为0.
disp(['μMdt的dt的0-K阶次为:',num2str(mumdt1)]);

for k=0:K_truc
    sumeqk{k+1}=generate_samples(3,k);
    p_ss_tem=[];
    for i = 1:size(sumeqk{k+1}, 1)
        ss_tem = sumeqk{k+1}(i, :);
        factorial_term = factorial(k) / (2^k * factorial(ss_tem(1))*prod(factorial(ss_tem(2:end))));
        sum_term = sum(xl(1:end)); % 因为 x_l 都是 1，所以这个和就是 m+1
        xl_term = prod(xl(1:end).^ss_tem(2:end)) / sum_term^(k - ss_tem(1));  % 我们只需要保留 (1/sum_term)^(kM-s) 的部分，但因为 sum_term 是常数（m+1），所以可以提前算出并归一化
        p_ss_tem = [p_ss_tem factorial_term*xl_term]; % x_term 在这里为 1，因为所有 x_l 都是 1
    end
    p_ss{k+1}=p_ss_tem;%正常情况下是已经自动归一化了的！
end

t_vec=0:dt:ceil(1.5/dt)*dt;

tic

num_NMt=5*10^6;  %num_NMt=5*10^6 in our paper.

count_qNMt=0;
HPmap_value_NMt=ones(num_NMt,length(t_vec));%预分配内存
overlap_rhoNM=zeros(num_NMt,length(t_vec));%预分配内存

for nn=1:length(t_vec)
    if nn==1
        rho0_Hdt_tem=eye(2^num_spin);
        rho0_J0_tem=rho0;
    else
        rho0_Hdt_tem=expm(-1i*Ham*dt)*rho0_J0_tem*expm(1i*Ham*dt);
        rho0_J0_tem=partialtrace1(Jl*kron(o*o',rho0_Hdt_tem)*Jl');
    end
    overlap_rho0J0(nn)=trace(rho11*rho0_J0_tem);%rho11是布局数观测量，刚好也是shots的概率
end
toc0=toc

toc1=toc

parfor qNMt=1:num_NMt
    count_qNMt=count_qNMt+1;
    count_tt=0;
    overlap_rhoNM_tem=zeros(1,length(t_vec)); %预分配内存

    for tt=t_vec  %tt必须从第一个dt开始，因为后面作用的Jl里面包含了一个dt
        count_tt=count_tt+1;
% %         [count_qNMt count_tt];
        if count_tt==1
            rho0_Hdt=rho0;
            rho_NM=rho0;
            overlap_rhoNM_tem(count_tt)=trace(rho11*rho0);
            continue;

        else
            rho0_Hdt=expm(-1i*Ham*dt)*rho_NM*expm(1i*Ham*dt);
        end

        rho0_J0=partialtrace1(Jl*kron(o*o',rho0_Hdt)*Jl'); %注意：这里的rho0_J0其实是t时刻时经过随机mitigation之后的rho_NM，此处rho_NM未必是合法量子态，未必厄米。

        p_kmundt1=mundt1/sum(mundt1);
        k_mundt1=randsrc(1,1,[1:length(mundt1);p_kmundt1]); %生成k阶次
        Pa=eye(2^num_spin); Pb=eye(2^num_spin);

        kA=k_mundt1-1;
        if kA==0
            xlab_prod=1; %这段有可能由于历史上用的是xl_pro=1而有误
            Pa;Pb;
        else
            %---------------------------%↓生成A{kA}的抽样及相关编号------------------------------
            pA=abs(A{kA})/sum(abs(A{kA}));
            length_A = 1:length(A{kA});
            pos_A=randsrc(1,1,[length_A;pA]);
            samples_A=A{kA}(pos_A);

            pos_A3=transA(pos_A-1); %转换成3进制数字编号.与trans1tombit_add1函数共用时，千万别忘记有可能需要-1处理！
            pos_A3add1=trans1tombit_add1(pos_A3,kA,3);%转成{1,2,3}构成的3进制形式
            %---------------------------%↑生成A{kA}的抽样及相关编号------------------------------

            %----------%↓生成chi_ab^{1}chi_ab^{2}...chi_ab^{k}连乘以及P1P2 ρ P2P3...的连乘--------
            xl_ab_tem=1;
% %             rho_tem=rho0_J0; %重置变量
            for k=kA:-1:1%k表示采样到A^{k}阶  %注意：当for循环降序循环时，中间的值-1必不可少！！！
                class_ADi=pos_A3add1(k); %class_ADi表抽样到的ji编号，ji∈{1,2,3}
                a_xlab=randsrc(1,1,[1:length(cla1);p_cla1]);b_xlab=randsrc(1,1,[1:length(cla1);p_cla1]);
                cla1_randa=cla1(a_xlab); cla1_randb=cla1(b_xlab);
                xl_ab_tem=xl_ab_tem*cla1_randa*cla1_randb';
                % %                             xl_ab_tem=exp( 1i*angle(cla1_randa*cla1_randb') );
                if class_ADi==1  %P1P2 ρ P2P3...的连乘的构型判断
                    xl_ab_tem;
                    Pa=PL{a_xlab}*Pa;Pb=PL{b_xlab}*Pb;
                elseif class_ADi==2
                    xl_ab_tem;%注：此处不是xl_ab_tem=sign(-1/2)*xl_ab_tem，因为正负号信息已经包含在后续对应的A_j1j2...jkA里了！
                    Pa=PL{b_xlab}'*PL{a_xlab}*Pa;Pb=eye(2^num_spin)*Pb;
                else             %==3
                    xl_ab_tem;
                    Pa=eye(2^num_spin)*Pa;Pb=(PL{b_xlab}'*PL{a_xlab})'*Pb;
                end
            end
            xlab_prod=xlab_prod*A{kA}( transkAto10(pos_A3add1) )*xl_ab_tem;         %chi_ab^{1}chi_ab^{2}...chi_ab^{k}的连乘结果
            %----------%↑生成chi_ab^{1}chi_ab^{2}...chi_ab^{k}连乘以及P1P2 ρ P2P3...的连乘--------
            %----------------------------------------------------------------------------

        end


        p_kmumdt1=mumdt1/sum(mumdt1);
        pos_mumdt1=randsrc(1,1,[1:length(mumdt1);p_kmumdt1]);

        kM=pos_mumdt1-1; %一共C_(kM+2)^2个。kM表dt^kM项

        if kM==0 %[s s0 s1...]采样程序，储存结果为samples_ss。
            ss=sumeqk{kM+1}; %ss=[0 0 0]
            pos_ss=1;
            samples_ss=ss;
        else
            ss=sumeqk{kM+1};    %生成[s s0 s1...sm]序列,要求sum(s)=k。此例m=1.%sumeqk事先已生成好，此处仅查询以提高效率。%一共C_(kM+2)^2个
            pos_ss=randsrc(1,1,[1:size(ss, 1);p_ss{kM+1}]);
            samples_ss=ss(pos_ss,:);
        end
        samples_ss;
        %----------s1,s2,...sm段信息↓-----------------
        if samples_ss(3)==0
            xlab_prod=1*xlab_prod;
        else %samples_ss(3)>=1
            xlab_prod=1*xlab_prod;
            for qsm=samples_ss(3):-1:1 %qsm表sm段的编号
                %产生xl_ab编号及PρP...操作
                class_xlab=randsrc(1,1,[1:3;1/2 1/4 1/4]); %从xl_ab的3种构型中选择1种
                a_xlab=randsrc(1,1,[1:length(cla1);p_cla1]);b_xlab=randsrc(1,1,[1:length(cla1);p_cla1]);
                cla1_randa=cla1(a_xlab); cla1_randb=cla1(b_xlab);
                xlab_prod= xlab_prod*cla1_randa*cla1_randb'; %从具体的10进制编号中获取xlab信息,连乘
                if class_xlab==1
                    xlab_prod;
                    Pa=PL{a_xlab}*Pa;Pb=PL{b_xlab}*Pb;
                elseif class_xlab==2
                    xlab_prod=sign(-1/2)*xlab_prod;
                    Pa=PL{b_xlab}'*PL{a_xlab}*Pa;Pb=eye(2^num_spin)*Pb;
                else %class_xlab==3
                    xlab_prod=sign(-1/2)*xlab_prod;
                    Pa=eye(2^num_spin)*Pa;Pb=(PL{b_xlab}'*PL{a_xlab})'*Pb;
                end
            end
        end
        %----------s1,s2,...sm段信息↑-----------------

        %----------s0段信息↓-----------------
        if samples_ss(2)==0 %各变量保持不变
            xlab_prod;
            Pa;Pb;
% %             rho_Mtem;
        else %samples_ss(2)>=1
            for qH=samples_ss(2):-1:1 %qH表H段的编号;%产生xl_ab编号及PρP...操作
                p_H=abs(cHa0)/sum(abs(cHa0));
                pos_H=randsrc(1,1,[1:length(cHa0);p_H]);
                class_H=randsrc(1,1,[1:2;1/2 1/2]);%超算子H共2中构型，分别是-iHρ与(-iHρ)',两种构型出现概率相同！
                if class_H==1
                    xlab_prod=-1i*cHa0(pos_H)*xlab_prod;
                    Pa=PH{pos_H}*Pa;Pb=eye(2^num_spin)*Pb;
                else %class_H==2
                    xlab_prod=( -1i*cHa0(pos_H) )'*xlab_prod;
                    Pa=eye(2^num_spin)*Pa;Pb=PH{pos_H}*Pb;
                end
            end
        end
        %----------s0段信息↑-----------------

        %----------s段信息↓-----------------
        if samples_ss(1)==0
            xlab_prod;
            Pa;Pb;
        else %samples_ss(1)>=1
            for qHsm=samples_ss(1):-1:1 %qH表H段的编号;%产生xl_ab编号及PρP...操作
                Prji_s=xl/sum(xl);
                class_Hsm=randsrc(1,1,[1:length(xl);Prji_s]); %产生s位ji信息，ji=0,1,...m。
                if class_Hsm==1     %ji_s表s段的编号对应的信息值
                    p_H=abs(cHa0)/sum(abs(cHa0));
                    pos_H=randsrc(1,1,[1:length(cHa0);p_H]);
                    class_H=randsrc(1,1,[1:2;1/2 1/2]);%超算子H共2中构型，分别是-iHρ与(-iHρ)',两种构型出现概率相同！
                    if class_H==1 %此编号装填的是xl0部分信息，即H超算子部分信息
                        xlab_prod=-1i*cHa0(pos_H)*xlab_prod;
                        Pa=PH{pos_H}*Pa;Pb=eye(2^num_spin)*Pb;
                    else %class_H==2
                        xlab_prod=(-1i*cHa0(pos_H))'*xlab_prod;
                        Pa=eye(2^num_spin)*Pa;Pb=PH{pos_H}*Pb;
                    end
                else %class_Hsm>=2。%此编号装填的是xl(2:end)部分信息，即L1...Lm耗散超算子部分信息
                    %产生xl_ab编号及PρP...操作
                    class_xlab=randsrc(1,1,[1:3;1/2 1/4 1/4]); %从xl_ab的3种构型中选择1种
                    a_xlab=randsrc(1,1,[1:length(cla1);p_cla1]);b_xlab=randsrc(1,1,[1:length(cla1);p_cla1]);
                    cla1_randa=cla1(a_xlab); cla1_randb=cla1(b_xlab);
                    xlab_prod= xlab_prod*cla1_randa*cla1_randb'; %从具体的10进制编号中获取xlab信息,连乘
                    if class_xlab==1
                        xlab_prod;
                        Pa=PL{a_xlab}*Pa;Pb=PL{b_xlab}*Pb;
                    elseif class_xlab==2
                        xlab_prod=sign(-1/2)*xlab_prod;
                        Pa=PL{b_xlab}'*PL{a_xlab}*Pa;Pb=eye(2^num_spin)*Pb;
                    else %class_xlab==3
                        xlab_prod=sign(-1/2)*xlab_prod;
                        Pa=eye(2^num_spin)*Pa;Pb=(PL{b_xlab}'*PL{a_xlab})'*Pb;
                    end
                end
            end
        end
        %----------s段信息↑-----------------
        xlab_prod=xlab_prod*(-1)^( kM-samples_ss(1) );%之所以有该行，是因为M前有系数(-1)^(k-s)。kM表dt^kM项
        phi=angle(xlab_prod);
        rho_M_tem=exp(1i*phi)*Pa*rho0_J0*Pb'; %做shots时不能乘以系数sum(mumdt1)，因为有可能使overlap>1，使shots失败
        rho_M=(rho_M_tem+rho_M_tem')/2; 

        rho_NM=rho_M; 

        overlap_rhoNM_tem(count_tt)=trace(rho11*rho_NM); %加abs是为了防止出现特别接近于0的负数出现！
    end
    overlap_rhoNM(qNMt,:)=overlap_rhoNM_tem;%overlap_rhoNM实际上是观测量O的均值！只是此时恰好O=|1><1|，均值刚好是overlap.
end

%------------------shots模块，按需打开-------------------------------------------
shots_rhoNM=zeros(num_NMt,length(t_vec));
for count_tt=1:length(t_vec)
    parfor qNMt=1:num_NMt
        overlap_rhoNM2(qNMt,count_tt)=(sum(mundt1)*sum(mumdt1))^(count_tt-1)*overlap_rhoNM(qNMt,count_tt);
    end
end

toc2=toc
sampletime=toc2-toc1

m=0;
for tt=t_vec
    m=m+1;
    rho_exact=reshape(expm(-1i*Lin_vec*tt)*rho0_vec,2^num_spin,2^num_spin).' ;
    overlap_rhoexact(m)=trace(rho11*rho_exact);
end

%-------------------------------------------------------------
overlap_rhoexact_mean=mean(overlap_rhoexact,1);

mean_overlap_rho0J0=mean(overlap_rho0J0,1);
mean_overlap_rhoNM=mean(overlap_rhoNM2,1);%overlap_rhoNM2才是乘以系数prod(HPmap_value_NMt(qNMt,1:count_tt))*(sum(mundt1)*sum(mumdt1))^(nn-1)的结果！
deltaJ0_noshots=abs(mean_overlap_rho0J0-overlap_rhoexact_mean);
deltaNM_noshots=abs(mean_overlap_rhoNM-overlap_rhoexact_mean);

%-------------------------------------------------------------------------

a2=A_abs(2);lambda=1.7; Dl=sum(abs(cla1)); %注：文章里给出的lambda>=2.5,数值拟合的结果是lambda=1.6765,取1.7.
c=3.5;T=1;
syms K tau positive;
spin_total=30;
for nn=3:spin_total
    nspin=nn;
    H_norm1=1.5*nspin;
    Lin_norm=H_norm1+Dl^2;
    %----------------------miti方法的系统误差bound↓----------------------------------
    epsH{nn-2}=2*( (2*exp(1)*H_norm1*tau)/(K+1) )^(K+1)+ ( (2*exp(1)*H_norm1*tau)/(K+1) )^(2*K+2);
    epsN{nn-2}=2*a2*lambda^(K-1)*Dl^(2*K+2)*tau^(K+1); %12.12之前版本是lambda^K，这里已改正
    epsM{nn-2}=( 4*exp(1)*Lin_norm*tau/(K+1) )^(K+1);
    epsHNM{nn-2}=vpa(simplify(epsH{nn-2}+epsN{nn-2}+epsM{nn-2}),5);
    %-----------------------miti方法的系统误差bound↑-----------------------------------------------
    %----------------------miti方法的统计误差bound↓----------------------------------
    muH=exp( 2*c*(2*H_norm1*tau)^4 ); %12.12之前版本是exp里有lambda，这里已去掉
    muN=exp(2*a2*Dl^4*tau^2);
    muM=exp(16*Lin_norm^2*tau^2);
    munsteps=simplify( (muH*muN*muM)^(T/tau) );
    munstepslog=simplify( log(munsteps) );
    %----------------------miti方法的统计误差bound（设定值为exp(3)）↑----------------------------------
    tau_miti(nn-2)=vpa(solve(munstepslog-3==0),5); %把miti法的统计方差bound为exp(2)=7.39<10，因此这方面的误差可以看成足够小，因此可以靠增加采样数而归零。

    trotter_err{nn-2}=simplify(( 0.04*(nn-1)+2*a2*gam^2+ ( (4*exp(1)*(0.3*nn+gam))/2 )^2 )*tau*T);
    trotter_err{nn-2}=vpa(trotter_err{nn-2},5); %代入tau1_miti会发现，相应的trotter error普遍在0.132-0.143(0.608-0.189，随qubit数递减)之间.
                                        % 即，能保持miti的统计error
                                        % bound在exp(2)的tau，能使得trotter法的error能达到0.132-0.143(0.608-0.189)之间.
                                        % 而miti法的系统误差又可通过提高截断K的阶次来实现，因此可以看成两种误差源都可以为0.
    muti_map2_trotter_err(nn-2)=vpa(subs(trotter_err{nn-2},tau_miti(nn-2)),5);
end

epslog10=-2:-0.5:-8;
for mm=1:length(epslog10)  %注意调用2024.12.23的miti_Ising_OQS4_sp5n500gam1p5rho1O1_NMnoshots时该段应重新生成，因为以前循环是nn=1:spin_total-2，现已改正。
    mm
    parfor nn=3:spin_total
        tau_trotter_solve(nn-2,mm)=solve(trotter_err{nn-2}==10^epslog10(mm)); %nn表示真实qubit数
        subs1=vpa(subs( epsHNM{nn-2},tau,tau_miti(nn-2) ),5); %得到关于K的函数
        K_solve(nn-2,mm)=solve(subs1==10^epslog10(mm));
        %     tau_miti_solve(nn-2,mm)=solve(epsHNM{nn-2}==10^epslog10(mm));
        K_solve2(nn-2,mm)=ceil(K_solve(nn-2,mm));

        miti_cnot(nn-2,mm)=( 2*nn+8+(4+2+2)*K_solve2(nn-2,mm)*2 )*ceil((T/tau_miti(nn-2))); %K_solve2×2是因为每一个补偿有2项Pauli项。
        miti_cnot_noceilK(nn-2,mm)=( 2*nn+8+(4+2+2)*K_solve(nn-2,mm)*2 )*ceil((T/tau_miti(nn-2)));
        miti_Rz(nn-2,mm)=( 66*(2*nn+11) )*ceil((T/tau_miti(nn-2)));

        trotter_cnot(nn-2,mm)=( 2*nn+4 )*ceil((T/tau_trotter_solve(nn-2,mm)));
        trotter_Rz(nn-2,mm)=( 66*(2*nn+2) )*ceil((T/tau_trotter_solve(nn-2,mm)));
    end
end
miti_cnot18=miti_cnot(18,:);% 红色变成#fdbb84=[253,187,132],蓝色变成#3182bd=[49,130,189].matlab都要除以255。
%-------------------------------------出图区-----------------------------------------------------------
%-------------------------------------出图区-----------------------------------------------------------
% 创建一个图形
figure;
% 设置图形的整体大小
figWidth = 12;  % 图形宽度（可以调整）
figHeight = 6;  % 图形高度（可以调整）
set(gcf, 'Position', [0, 0, figWidth*100, figHeight*100]);  % 设置图形窗口的大小
% 子图1宽度是子图2的两倍
subplotWidth1 = 2;  % 第一个子图宽度的比例
subplotWidth2 = 1;  % 第二个子图宽度的比例
gap = 0.1;  % 子图间距，设为总图宽度的5%

% 设置坐标轴距离窗口边缘的距离
paddingA1 = [0.075, 0.1];  % 内边距，可以调整，避免坐标轴与窗口边缘重合.分别表相对于左上角，[左边距, 上边距] 
paddingB1 = [0.47, 1-0.50]; %[宽度,高度]
paddingA2 = [0.66, 0.1];  
paddingB2 = [0.31, 1-0.50];
% 创建第一个子图 (宽度是2/3)
scalebox_fig1=[paddingA1,paddingB1];
scalebox_fig2=[paddingA2,paddingB2];

ax1 = axes('Position', scalebox_fig1); %[left, bottom, width, height]
loglog(10.^epslog10,trotter_Rz(18,:),'d','Color',[35,130,205]/255,'MarkerFaceColor',[35,130,205]/255,'LineWidth',2.5,'MarkerSize', 10)
hold on
loglog(10.^epslog10,trotter_cnot(18,:),'^','Color',[35,130,205]/255,'MarkerFaceColor',[35,130,205]/255,'LineWidth',2.5,'MarkerSize', 10)  %蓝色
hold on
loglog(10.^epslog10,miti_Rz(18,:),'d','Color',[233,92,74]/255,'MarkerFaceColor',[233,92,74]/255,'LineWidth',2.5,'MarkerSize', 10);%'MarkerFaceColor','b'设置符号内部填充色（实心色）
hold on
loglog(10.^epslog10, miti_cnot(18,:),'^','Color',[233,92,74]/255,'MarkerFaceColor',[233,92,74]/255,'LineWidth',2.5,'MarkerSize', 10); %红色
hold on

set(gca,'FontName','Times New Roman','FontSize',20,'LineWidth',2.5); %设置坐标轴字体，字号，及线宽%正常是字体18，线型2.5
hold on

% 使用线性拟合（一阶多项式）
pmiti_cnot = polyfit(epslog10, double(miti_cnot_noceilK(18,:)) , 1); %K不取整时的情形，注意，因为理论上miti方法的cnot门数量miti_cnot与log(ε)成正比关系，而非像其他的那样是log(miti_cnot)。
pmiti_Rz = polyfit(epslog10, double( log10(miti_Rz(18,:))) , 1); 
ptrotter_cnot = polyfit(epslog10, double( log10(trotter_cnot(18,:))) , 1); 
ptrotter_Rz = polyfit(epslog10, double( log10(trotter_Rz(18,:))) , 1); 

% 生成拟合线的数据点（可选，仅用于可视化）
n_fit = linspace(min(epslog10), max(epslog10), 100); % 生成100个点用于更平滑的拟合线
miti_cnot_fit = polyval(pmiti_cnot, n_fit); % 计算这些点的拟合值
miti_Rz_fit = polyval(pmiti_Rz, n_fit); % 计算这些点的拟合值
trotter_cnot_fit = polyval(ptrotter_cnot, n_fit); % 计算这些点的拟合值
trotter_Rz_fit = polyval(ptrotter_Rz, n_fit); % 计算这些点的拟合值

% 绘制拟合线
loglog(10.^n_fit, miti_cnot_fit, 'k--', 'LineWidth', 2.0, 'DisplayName', '拟合线'); 
hold on
loglog(10.^n_fit, 10.^miti_Rz_fit, 'k--', 'LineWidth', 2.0, 'DisplayName', '拟合线'); 
hold on
loglog(10.^n_fit, 10.^trotter_cnot_fit, 'k--', 'LineWidth', 2.0, 'DisplayName', '拟合线'); 
hold on
loglog(10.^n_fit, 10.^trotter_Rz_fit, 'k--', 'LineWidth', 2.0, 'DisplayName', '拟合线'); 
xlim(10.^[epslog10(end)-0.1,epslog10(1)+0.1]) %该行必须放在hold off之前，否则图像将会被覆盖！
hold off
legend('$\mathrm{R_z}$ (Trotter-like)','CNOT (Trotter-like)','$\mathrm{R_z}$ (Ours)','CNOT (Ours)','Fitting','interpreter','latex',...
    'FontName','Times New Roman','FontSize', 18, 'Location','northwest','NumColumns',1)

str_miti_cnot=strcat('$\mathcal{O}\!\left(\log_{10}\varepsilon\right)$'); 
str_miti_Rz=strcat( '$\mathcal{O}(\varepsilon^{0.0})$' );
str_trotter_cnot=strcat( '$\mathcal{O}(\varepsilon^{-1.0})$' )
str_trotter_Rz=strcat( '$\mathcal{O}(\varepsilon^{-1.0})$' );

text(10^(epslog10(end)-0.11), miti_cnot_fit(1), str_miti_cnot,'interpreter','latex','FontName','Times New Roman','FontSize', 20,'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle');
text(10^(epslog10(end)-0.11), 10^(miti_Rz_fit(1)), str_miti_Rz,'interpreter','latex','FontName','Times New Roman','FontSize', 20, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');%'left''bottom'表坐标点对齐的text box位置为左下点。
text(10^(epslog10(end)-0.11), 10^(trotter_cnot_fit(1)), str_trotter_cnot,'interpreter','latex','FontName','Times New Roman','FontSize', 20, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');
text(10^(epslog10(end)-0.11), 10^(trotter_Rz_fit(1)), str_trotter_Rz,'interpreter','latex','FontName','Times New Roman','FontSize', 20, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');

% 设置坐标轴从小到大显示
set(gca, 'XDir', 'rev');
set(gcf, 'Color', 'white'); % gcf代表当前图形窗口，'Color'属性用于设置背景色
                                                            
xlabel(['Precision, ','$\varepsilon$'],'interpreter','latex','FontName','Times New Roman','FontSize',20,'LineWidth',2.5)
ylabel('Number of gates per sample','interpreter','latex','FontName','Times New Roman','FontSize',20,'LineWidth',2.5)
set(gca,'FontName','Times New Roman','FontSize',20,'LineWidth',2.5); %设置坐标轴字体，字号，及线宽

%-----------------------------------------------------------------------------
% 创建第二个子图 (宽度是1/3)
ax2 = axes('Position', scalebox_fig2);

plot(t_vec([1:2:end end]),deltaJ0_noshots([1:2:end end]),'o','Color',[35,130,205]/255,'MarkerFaceColor',[35,130,205]/255,'LineWidth',2.5,'MarkerSize', 10);%红色
hold on
figure3=plot(t_vec([1:2:end end]),deltaNM_noshots([1:2:end end]),'o','Color',[233,92,74]/255,'MarkerFaceColor',[233,92,74]/255,'LineWidth',2.5,'MarkerSize', 10);%蓝色
xlabel(['Simulation time, ','$T$'],'interpreter','latex','FontName','Times New Roman','FontSize',20,'LineWidth',2.5) %与Miti_Ising_OQS.m配套命令
ylabel(['Precision, ','$\varepsilon$'],'interpreter','latex','FontName','Times New Roman','FontSize',20,'LineWidth',2.5)
xlim([0 1.55])
legend('Trotter-like','Ours','interpreter','latex','FontName','Times New Roman','FontSize', 18, 'Location','northwest','NumColumns',1)
set(gca,'FontName','Times New Roman','FontSize',20,'LineWidth',2.5); %设置坐标轴字体，字号，及线宽
set(0,'defaultfigurecolor','w') %设置默认图片背景为白色


%  export_fig Ising_model_gates.pdf %专用输出函数，去白边


beep

% %    save miti_Ising_OQS4_sp5n500gam1p5rho1O1_NMnoshots.mat

%-------------------------------------------------------------------------------------------------------------
%------------------------------------------------------------------------------------------------------------------
function P_pro= Pauli_product(A,n_total,n)%注：该程序产生I...A...I张量积
global I X Y Z
for j=1:n_total
    if j==1%定Pa连乘的初值
        if n==1
            P_pro=A;
        else
            P_pro=I;
        end
    elseif j>=2
        if j~=n %若j≠n,B定为I，若j=n，则B定为A
            B=I;
        else
            B=A;
        end
        P_pro=kron(P_pro,B);
    end
end
end

function result = multiKronecker(matrices)
% 递归地计算多个矩阵的外积
if nargin < 1 || isempty(matrices)
    error('需要提供至少一个矩阵');
end

% 如果只有一个矩阵，直接返回该矩阵
if numel(matrices) == 1
    result = matrices{1};
else
    % 递归计算剩余矩阵的外积，并与第一个矩阵做外积
    result = kron(matrices{1}, multiKronecker(matrices(2:end)));
end
end

function result = nkron(AA, BB, n)
% 输入：A 和 BB 是矩阵，n 是非负自然数
% 输出：A 与 n 个 BB 的连续 Kronecker 直积

% 初始化结果为 AA
result = AA;

% 迭代计算 n 次 Kronecker 积
for i = 1:n
    result = kron(result, BB);
end
end

function result2=transkAto10(Mat) %功能是将自定义三进制数矩阵如Mat=[3 2 1]转成十进制数
result2=0;
leng_Mat=length(Mat);
for m1=1:leng_Mat
    result2=result2+( Mat(m1)-1 )*3^(leng_Mat-m1);
end
result2=result2+1;
end

function result3=trans1tombit_add1(mat,mbit,nary) %使用该函数前提是已经将mat转化为nary进制形式了%作用是将mat中不到mbit列的扩展到mbit列，并用自定义{1,2,...}的nary进制。如只有1列的扩展为2列：[1 mat]，防止后续pos_xl_ab2add1{class_ADi}(n,2)的索引溢出。
[row_mat, col_mat]=size(mat);
result3=mat+1;
if col_mat<=mbit
    result3=padarray(result3,[0 mbit-col_mat],1,'pre'); %将矩阵mat扩展0行，1列，扩展方向为前方（即左上，mposta为右下，tboth0为双向），扩展位置填充矩阵元为1.
else %mat也可能超出mbit，如3进制形式为[1 0 0 0],实际是27
    for mm=1:row_mat
        if result3(mm,1)==2
            result3(mm,:)=padarray(nary,[0 col_mat-1],nary,'pre'); %要注意mat可能是n行，mbit+1列矩阵。作用如：将[1 0 0 0]转为nary进制的[nary nary nary]
        end
    end
    result3(:,1)=[];%删列
end
end

function result4=partialtrace1(mat) %对qubit 1求偏迹
[row_mat, col_mat]=size(mat);
assert(row_mat==col_mat,'输入错误！'); %
kron0=kron([1;0],eye(row_mat/2));
kron1=kron([0;1],eye(row_mat/2));
result4=kron0'*mat*kron0+kron1'*mat*kron1;
end

function samples = generate_samples(m, k) %生成长为m的s=(s,s0,...sm-2)的所有集合，要求sum(s)=k
assert(m>0&k>=0&k<=9,'M样本生成程序输入错误！');%程序最多允许k+1进制
samples1=zeros((k+1)^m,m);
for j=0:(k+1)^m-1
    if k==0
        sam=0;
    else
        sam= dec2base(j,k+1)-'0'; %注意dec2base(10,16)结果是'A'，但'A’-‘0’却是17，而非10!
    end
    if length(sam)<=m
        samples1(j+1,:)=padarray(sam,[0 m-length(sam)],0,'pre');
    end
end
samples1_sum=sum(samples1,2);
row_samples1=samples1_sum==k;
samples=samples1(row_samples1,:);
end

function results=S2k_generate(k,x,num_spin,cHa0,PH)
funuk=@(y) 1/( 4-4^(1/(2*y-1)) );
if k==1
    uk=1/2;
    D2k=eye(2^num_spin);
    D2kdag=eye(2^num_spin);
    for m=1:length(cHa0)
        D2k= expm(-1i*cHa0(m)*PH{m}*uk*x)*D2k;%左乘.D2=S1(x/2)
        D2kdag= D2kdag*expm(-1i*cHa0(m)*PH{m}*uk*x);%右乘
    end
    results=D2kdag*D2k;

else %k>=2
    uk=funuk(k);
    results=S2k_generate(k-1,uk*x,num_spin,cHa0,PH)^2*S2k_generate(k-1,(1-4*uk)*x,num_spin,cHa0,PH)*S2k_generate(k-1,uk*x,num_spin,cHa0,PH)^2;
end
end

