clc,clear;

% rho—控制变量，J—目标泛函值

% 参数设置
mu=0.3; beta=0.7; delta=0.8; lamda=0.9; %a=0.2445; b=1.076; %p=0.8;没用到

Du = [0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 26 27 28 29 30 31 32 33 34 38 39 59 134]; % 度
Pk = [0.114638 0.232804 0.130511 0.105820 0.059965 0.072310 0.049383 0.024691 0.044092 0.021164 0.026455 0.017637 0.008818 0.005291 0.008818 0.007055 0.005291 0.003527 0.001764 ...
      0.007055 0.003527 0.001764 0.007055 0.001764 0.007055 0.001764 0.005291 0.001764 0.001764 0.003527 0.001764 0.001764 0.003527 0.001764 0.003527 0.001764 0.001764 0.001764];

k2 = 0.0;
for i=1:length(Du)
   k2 = k2+Du(i)*Pk(i); 
end

time = 5;            %最大24，超过则无穷大
M = 1000; % 时间区间等分数
h = time / M; % 时间步长
h2 = h / 2;
T = 0:h:time; % 时间网格剖分

% 求解一阶最优系统循环控制量
weight = 0.0001;

S_cell = {};
E_cell = {};
I_cell = {};
J_cell = {};
rho_cell = {};
simulate_rho = [0.0 0.2 0.4 0.6 0.8];
% 求解最优系统
for simulate=1:6
    if simulate~=6
        rho_min = simulate_rho(simulate);
        rho_max = simulate_rho(simulate);
    else
        rho_min = 0.0; % 控制变量rho的下界
        rho_max = 1.0; % 控制变量rho的上界
    end
    
    %初始化
    test = -1;
% 控制变量的初始猜测值
rho = zeros(M+1,length(Du));        %创建M+1行Du_Num列的0矩阵
% 为状态变量分配内存
S = zeros(M+1,length(Du));         %创建M+1行Du_Num列的0矩阵
E = zeros(M+1,length(Du));
I = zeros(M+1,length(Du));
Lamda1 = zeros(M+1,length(Du));
Lamda2 = zeros(M+1,length(Du));
Lamda3 = zeros(M+1,length(Du));

B = 0.5*ones(length(Du),1);
y_rho_opt = zeros(M+1,1);

% 为状态变量赋初值
S(1,:) = 467/567*ones(length(Du),1);     %S矩阵的第一行赋值
E(1,:) = 0*ones(length(Du),1);
I(1,:) = 100/567*ones(length(Du),1);
theta = zeros(M+1,1);

num = 0; % 统计数值解达到精度要求时的迭代次数
while (test < 0 && num < 100)
    old_rho = rho;
    
    old_S = S;
    old_E = E;
    old_I = I;
    
    old_Lamda1 = Lamda1;
    old_Lamda2 = Lamda2;
    old_Lamda3 = Lamda3;
    
    %求解状态方程
    %前向更新
    for t=1:M        %每一步
        
       theta(t) = 0.0;
       for i=1:length(Du)
           theta(t) = theta(t)+Du(i)*Pk(i)*I(t,i);
       end
       theta(t) = theta(t)/k2;
       
       for k=1:length(Du)     %求解每个度的状态
            KS_1 = (1-rho(t,k))*(-Du(k)*lamda*S(t,k)*theta(t)+beta*E(t,k))+rho(t,k)*delta*E(t,k)+rho(t,k)*delta*I(t,k);
            KE_1 = (1-rho(t,k))*(Du(k)*lamda*S(t,k)*theta(t)-beta*E(t,k)-mu*E(t,k))-rho(t,k)*delta*E(t,k);
            KI_1 = (1-rho(t,k))*mu*E(t,k)-rho(t,k)*delta*I(t,k);
            
            KS_2 = (1-0.5*(rho(t,k)+rho(t+1,k)))*(-Du(k)*lamda*(S(t,k)+h2*KS_1)*theta(t)+beta*(E(t,k)+h2*KE_1))+0.5*(rho(t,k)+rho(t+1,k))*delta*(E(t,k)+h2*KE_1)+0.5*(rho(t,k)+rho(t+1,k))*delta*(I(t,k)+h2*KI_1);
            KE_2 = (1-0.5*(rho(t,k)+rho(t+1,k)))*(Du(k)*lamda*(S(t,k)+h2*KS_1)*theta(t)-beta*(E(t,k)+h2*KE_1)-mu*(E(t,k))+h2*KS_1)-0.5*(rho(t,k)+rho(t+1,k))*delta*(E(t,k)+h2*KS_1);
            KI_2 = (1-0.5*(rho(t,k)+rho(t+1,k)))*mu*(E(t,k)+h2*KE_1)-0.5*(rho(t,k)+rho(t+1,k))*delta*(I(t,k)+h2*KI_1);
            
            KS_3 = (1-0.5*(rho(t,k)+rho(t+1,k)))*(-Du(k)*lamda*(S(t,k)+h2*KS_2)*theta(t)+beta*(E(t,k)+h2*KE_2))+0.5*(rho(t,k)+rho(t+1,k))*delta*(E(t,k)+h2*KE_2)+0.5*(rho(t,k)+rho(t+1,k))*delta*(I(t,k)+h2*KI_2);
            KE_3 = (1-0.5*(rho(t,k)+rho(t+1,k)))*(Du(k)*lamda*(S(t,k)+h2*KS_2)*theta(t)-beta*(E(t,k)+h2*KE_2)-mu*(E(t,k)+h2*KE_2))-0.5*(rho(t,k)+rho(t+1,k))*delta*(E(t,k)+h2*KE_2);
            KI_3 = (1-0.5*(rho(t,k)+rho(t+1,k)))*mu*(E(t,k)+h2*KE_2)-0.5*(rho(t,k)+rho(t+1,k))*delta*(I(t,k)+h2*KI_2);

            KS_4 = (1-0.5*(rho(t,k)+rho(t+1,k)))*(-Du(k)*lamda*(S(t,k)+h*KS_3)*theta(t)+beta*(E(t,k)+h*KE_3))+0.5*(rho(t,k)+rho(t+1,k))*delta*(E(t,k)+h*KE_3)+0.5*(rho(t,k)+rho(t+1,k))*delta*(I(t,k)+h*KI_3) ;
            KE_4 = (1-0.5*(rho(t,k)+rho(t+1,k)))*(Du(k)*lamda*(S(t,k)+h*KS_3)*theta(t)-beta*(E(t,k)+h*KE_3)-mu*(E(t,k)+h*KE_3))-0.5*(rho(t,k)+rho(t+1,k))*delta*(E(t,k)+h*KE_3);
            KI_4 = (1-0.5*(rho(t,k)+rho(t+1,k)))*mu*(E(t,k)+h*KE_3)-0.5*(rho(t,k)+rho(t+1,k))*delta*(I(t,k)+h*KI_3);

            S(t+1,k) = S(t,k)+h/6*(KS_1+2*KS_2+2*KS_3+KS_4);
            E(t+1,k) = E(t,k)+h/6*(KE_1+2*KE_2+2*KE_3+KE_4);
            I(t+1,k) = I(t,k)+h/6*(KI_1+2*KI_2+2*KI_3+KI_4);
            
            y_rho_opt(t) = 0;
            for i=1:length(Du)
               y_rho_opt(t) = y_rho_opt(t)+Pk(i)*I(t,i);
            end 
       end
    end
    
    %前后更新
    for t=M+1:-1:2         %每一步
       for k=1:length(Du)     %求解每个度的状态
           theta(t) = 0.0;
           for i=1:length(Du)
               theta(t) = theta(t)+Du(i)*Pk(i)*I(t,i);
           end
           theta(t) = theta(t)/k2;
           
            KLamda1_1 = Du(k)*lamda*(1-rho(t,k))*theta(t)*Lamda1(t,k)-Du(k)*lamda*(1-rho(t,k))*theta(t)*Lamda2(t,k);
            KLamda2_1 = -((1-rho(t,k))*beta+rho(t,k)*delta)*Lamda1(t,k)+((1-rho(t,k))*(beta+mu)+rho(t,k)*delta)*Lamda2(t,k)-(1-rho(t,k))*mu*Lamda3(t,k);
            sum_temp = 0.0;
            for i=1:length(Du)
                sum_temp = sum_temp+(Lamda1(t,i)-Lamda2(t,i))*Du(i)*S(t,i);
            end
            KLamda3_1 = ((1-rho(t,k))*lamda*Du(k)*Pk(k)*sum_temp-k2)/k2-rho(t,k)*delta*Lamda1(t,k)+rho(t,k)*delta*Lamda3(t,k);
            
            KLamda1_2 = Du(k)*lamda*(1-0.5*(rho(t,k)+rho(t-1,k)))*theta(t)*(Lamda1(t,k)-h2*KLamda1_1)-Du(k)*lamda*(1-0.5*(rho(t,k)+rho(t-1,k)))*theta(t)*(Lamda2(t,k)-h2*KLamda2_1);
            KLamda2_2 = -((1-0.5*(rho(t,k)+rho(t-1,k)))*beta+0.5*(rho(t,k)+rho(t-1,k))*delta)*(Lamda1(t,k)-h2*KLamda1_1)+((1-0.5*(rho(t,k)+rho(t-1,k)))*(beta+mu)+0.5*(rho(t,k)+rho(t-1,k))*delta)*(Lamda2(t,k)-h2*KLamda2_1)-(1-0.5*(rho(t,k)+rho(t-1,k)))*mu*(Lamda3(t,k)-h2*KLamda3_1);
            sum_temp = 0.0;
            for i=1:length(Du)
                sum_temp = sum_temp+((Lamda1(t,i)-h2*KLamda1_1)-(Lamda2(t,i)-h2*KLamda2_1))*Du(i)*S(t,i);
            end
            KLamda3_2 = ((1-0.5*(rho(t,k)+rho(t-1,k)))*lamda*Du(k)*Pk(k)*sum_temp-k2)/k2-0.5*(rho(t,k)+rho(t-1,k))*delta*(Lamda1(t,k)-h2*KLamda1_1)+0.5*(rho(t,k)+rho(t-1,k))*delta*(Lamda3(t,k)-h2*KLamda3_1);
            
            KLamda1_3 = Du(k)*lamda*(1-0.5*(rho(t,k)+rho(t-1,k)))*theta(t)*(Lamda1(t,k)-h2*KLamda1_2)-Du(k)*lamda*(1-0.5*(rho(t,k)+rho(t-1,k)))*theta(t)*(Lamda2(t,k)-h2*KLamda2_2);
            KLamda2_3 = -((1-0.5*(rho(t,k)+rho(t-1,k)))*beta+0.5*(rho(t,k)+rho(t-1,k))*delta)*(Lamda1(t,k)-h2*KLamda1_2)+((1-0.5*(rho(t,k)+rho(t-1,k)))*(beta+mu)+0.5*(rho(t,k)+rho(t-1,k))*delta)*(Lamda2(t,k)-h2*KLamda2_2)-(1-0.5*(rho(t,k)+rho(t-1,k)))*mu*(Lamda3(t,k)-h2*KLamda3_2);
            sum_temp = 0.0;
            for i=1:length(Du)
                sum_temp = sum_temp+((Lamda1(t,i)-h2*KLamda1_2)-(Lamda2(t,i)-h2*KLamda2_2))*Du(i)*S(t,i);
            end
            KLamda3_3 = ((1-0.5*(rho(t,k)+rho(t-1,k)))*lamda*Du(k)*Pk(k)*sum_temp-k2)/k2-0.5*(rho(t,k)+rho(t-1,k))*delta*(Lamda1(t,k)-h2*KLamda1_2)+0.5*(rho(t,k)+rho(t-1,k))*delta*(Lamda3(t,k)-h2*KLamda3_2);
            
            KLamda1_4 = Du(k)*lamda*(1-0.5*(rho(t,k)+rho(t-1,k)))*theta(t)*(Lamda1(t,k)-h2*KLamda1_3)-Du(k)*lamda*(1-0.5*(rho(t,k)+rho(t-1,k)))*theta(t)*(Lamda2(t,k)-h2*KLamda2_3);
            KLamda2_4 = -((1-0.5*(rho(t,k)+rho(t-1,k)))*beta+0.5*(rho(t,k)+rho(t-1,k))*delta)*(Lamda1(t,k)-h2*KLamda1_3)+((1-0.5*(rho(t,k)+rho(t-1,k)))*(beta+mu)+0.5*(rho(t,k)+rho(t-1,k))*delta)*(Lamda2(t,k)-h2*KLamda2_3)-(1-0.5*(rho(t,k)+rho(t-1,k)))*mu*(Lamda3(t,k)-h2*KLamda3_3); 
            sum_temp = 0.0;
            for i=1:length(Du)
                sum_temp = sum_temp+((Lamda1(t,i)-h2*KLamda1_3)-(Lamda2(t,i)-h2*KLamda2_3))*Du(i)*S(t,i);
            end
            KLamda3_4 = ((1-0.5*(rho(t,k)+rho(t-1,k)))*lamda*Du(k)*Pk(k)*sum_temp-k2)/k2-0.5*(rho(t,k)+rho(t-1,k))*delta*(Lamda1(t,k)-h2*KLamda1_3)+0.5*(rho(t,k)+rho(t-1,k))*delta*(Lamda3(t,k)-h2*KLamda3_3);
            
            Lamda1(t-1,k) = Lamda1(t,k) - h/6*(KLamda1_1+2*KLamda1_2+2*KLamda1_3+KLamda1_4);
            Lamda2(t-1,k) = Lamda2(t,k) - h/6*(KLamda2_1+2*KLamda2_2+2*KLamda2_3+KLamda2_4);
            Lamda3(t-1,k) = Lamda3(t,k) - h/6*(KLamda3_1+2*KLamda3_2+2*KLamda3_3+KLamda3_4);
       end
        
    end 
    
    for t=1:M+1 
        for k=1:length(Du)   
            rho(t,k) = min(rho_max,max(rho_min,1/B(k)*((Du(k)*lamda*S(t,k)*theta(t)-beta*E(t,k)+mu*E(t,k)+delta*E(t,k))*Lamda2(t,k)+(mu*E(t,k)+delta*I(t,k))*Lamda3(t,k)-...
                       (Du(k)*lamda*S(t,k)*theta(t)-beta*E(t,k)*delta*E(t,k)+delta*I(t,k))*Lamda1(t,k))));
            rho(t,k) = 0.5*(old_rho(t,k) + rho(t,k));
        end
    end

    % 确定是否循环继续求解最优系统
    temp(1) = weight*sum(abs(S(:))) - sum(abs(old_S(:) - S(:)));
    temp(2) = weight*sum(abs(E(:))) - sum(abs(old_E(:) - E(:)));
    temp(3) = weight*sum(abs(I(:))) - sum(abs(old_I(:) - I(:)));
    temp(4) = weight*sum(abs(Lamda1(:))) - sum(abs(old_Lamda1(:) - Lamda1(:)));
    temp(5) = weight*sum(abs(Lamda2(:))) - sum(abs(old_Lamda2(:) - Lamda2(:)));
    temp(6) = weight*sum(abs(Lamda3(:))) - sum(abs(old_Lamda3(:) - Lamda3(:)));
    
    test = min(temp);
    num = num + 1; 
    
end

J = 0.0;
for k=1:length(Du)
    temp_J = trapz(T, I(:,k)+1/2*B(k)*rho(:,k).^2);
    J = J+temp_J;
end
integralAverage = trapz(T, sum(rho,2)) / time;
if simulate~=6
    disp('==============================================================================')
    fprintf('控制rho=%.2f, 迭代次数 = %d, 目标泛函 = %f, 控制变量积分平均 = %f\n', ... 
        simulate_rho(simulate), num, J, integralAverage);
    disp('==============================================================================')
else
    disp('==============================================================================')
    fprintf('控制rho=opt, 迭代次数 = %d, 目标泛函 = %f, 控制变量积分平均 = %f\n', ... 
        num, J, integralAverage);
    disp('==============================================================================')
end

S_cell = [S_cell S];
E_cell = [E_cell E];
I_cell = [I_cell I];
J_cell = [J_cell J];
rho_cell = [rho_cell rho];
end

save data