clc,clear,close all;

load data;

%不同rho下的J
figure('name','图1_J')
    x = [0.0 0.2 0.4 0.6 0.8 1.0];
    y = [];
    for i=1:length(J_cell)
       y = [y J_cell{i}]; 
    end
    bar(y);
    set(gca,'xticklabel',{'0.0' '0.2' '0.4' '0.6' '0.8' 'opt'})
    xlabel('rho')
    ylabel('J')
  
%不同rho下S的变化趋势
figure('name','图2_S')
    for i=1:length(S_cell)
        plot(sum(S_cell{i},2))
        hold on;
    end
    legend('rho=0.0', 'rho=0.2', 'rho=0.4', 'rho=0.6', 'rho=0.8', 'rho=opt')
    xlabel('time')
    ylabel('S')
    
%不同rho下S的变化趋势
figure('name','图3_E')
    for i=1:length(E_cell)
        plot(sum(E_cell{i},2))
        hold on;
    end
    legend('rho=0.0', 'rho=0.2', 'rho=0.4', 'rho=0.6', 'rho=0.8', 'rho=opt')
    xlabel('time')
    ylabel('E')

%不同rho下I的变化趋势
figure('name','图4_I')
    y = [];
    for i=1:length(I_cell)
        plot(sum(I_cell{i},2))
        hold on;
    end
    legend('rho=0.0', 'rho=0.2', 'rho=0.4', 'rho=0.6', 'rho=0.8', 'rho=opt')
    xlabel('time')
    ylabel('I')

%最优控制下不同度的rho的变化趋势
figure('name','图5_rho')
    temp_Du = [5 10 20 38 134];
    index = [];
    for i=1:length(Du)
        for j=1:length(temp_Du)
            if Du(i)==temp_Du(j)
                plot(rho_cell{6}(:,i))
                hold on;
            end
        end
    end
    text_legend = {};
    for i=1:length(temp_Du)
        text_legend = [text_legend ['Du=' num2str(temp_Du(i))]];
    end
    legend(text_legend)
    xlabel('time')
    ylabel('rho_Du')
    axis([0 1100 0 1.1]) 
    
%不同rho下P(k)*I(k)和的变化趋势
figure('name','图6_P(k)*I(k)和')
    for i=1:length(I_cell)
        y = [];
        for j=1:600
            y = [y I_cell{i}(j,:)*Pk'];
        end
        plot(y)
        hold on;
    end
    legend('rho=0.0', 'rho=0.2', 'rho=0.4', 'rho=0.6', 'rho=0.8', 'rho=opt')
    xlabel('time')
    ylabel('sum(P(k)*I(k))')  
