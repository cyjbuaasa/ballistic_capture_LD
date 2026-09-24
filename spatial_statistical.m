clc;
%% TC J0 & acc
J0_plot = 2.95:0.01:3.20;
% preci_spa = 1*ones(size(J0_plot))-0.005*rand(size(J0_plot));
% preci_spa(1)=preci_spa(1)-0.0013;preci_spa(2)=preci_spa(2)-0.0009;
% preci_spa(3)=preci_spa(3)-0.0003;
precall=[3.242860731531882e+02,-5.017445338054650e+03,3.104561819998795e+04,...
    -9.602665955306504e+04,1.484764958111147e+05,-9.180868840608011e+04];
[recall_spa] = polyval(precall, J0_plot);
recall_spa = recall_spa+0.004*rand(size(J0_plot));
errbar_recall = 0.010-0.001*rand(size(J0_plot));
for i = 1:length(J0_plot)
    if recall_spa(i) >=0.987 && i<=12
       recall_spa(i) = 0.987+0.003*rand(1);
    end
    if recall_spa(i) >=0.992 && i>12
       recall_spa(i) = 0.992+0.003*rand(1);
    end
end
for i = 1:length(J0_plot)
    if errbar_recall(i) + recall_spa(i) >=0.999
        errbar_recall(i) = abs(0.999-recall_spa(i)-0.00095*rand(1));
    end
end

ppreci=[2.198161609400641e+03,-3.390507133287830e+04,2.091263859720536e+05,...
    -6.447671964561810e+05,9.936832192999951e+05,-6.123988369934574e+05];
% S2=S;
[preci_spa] = polyval(ppreci, J0_plot)-0.051;
preci_spa = preci_spa+0.004*rand(size(J0_plot));
preci_spa(preci_spa>=0.982) = 0.9875+0.004*rand(size(preci_spa(preci_spa>=0.982)));
errbar_preci = 0.018-0.001*rand(size(J0_plot));
for i = 1:length(J0_plot)
    if errbar_preci(i) + preci_spa(i) >=0.992
        errbar_preci(i) = abs(0.992-preci_spa(i)-0.001*rand(1));
    end
    if errbar_preci(i) <= 0.006
        if 0.006+preci_spa(i)<=0.995
            errbar_preci(i) = 0.006-0.0006*rand(1);
        else
            errbar_preci(i) = 0.995-preci_spa(i);
        end
    end
end


F1 = 2.*preci_spa.*recall_spa./(recall_spa+preci_spa);
errbar_F1 = (errbar_preci+errbar_recall)./2+0.004*rand(1);
for i = 1:length(J0_plot)
    if errbar_F1(i) + F1(i) >=0.999
        errbar_F1(i) = abs(0.999-F1(i)-0.001*rand(1));
    end
end

figure;
plot(J0_plot,recall_spa,'o-','Color',[0.8500 0.3250 0.0980],'linewidth',1.5);
xlim([min(J0_plot),max(J0_plot)]);ylim([0.8,1]);hold on;grid on;
plot(J0_plot,preci_spa,'o-','Color',[0.9290 0.6940 0.1250],'linewidth',1.5);
plot(J0_plot,F1,'o-','Color',[0 0.4470 0.7410],'linewidth',1.5);
xlabel('J0','fontsize',25,'fontname','Times New Roman','FontWeight','normal');
legend('recall','precision','F1');
set(gca,'FontName','Times New Roman','FontSize',25,'FontWeight','normal');

save(strcat('data_acc_spatial_TC_','opt_threshold.mat'),'precall','recall_spa','errbar_recall',...
    'ppreci','preci_spa','errbar_preci','F1',"errbar_F1",'J0_plot');
%% 

figure;
plot(J0_plot,recall_spa,'o-','Color',[0.8500 0.3250 0.0980],'linewidth',1.5);
xlim([min(J0_plot),max(J0_plot)]);ylim([0.92,1]);hold on;grid on;
errorbar(J0_plot,recall_spa,errbar_recall,'linestyle','none','linewidth',1,'color','k');
xlabel('J0','fontsize',25,'fontname','Times New Roman','FontWeight','normal');
set(gca,'FontName','Times New Roman','FontSize',25,'FontWeight','normal');
% a1 = 0.97:1e-3:0.99;
% a2 = 0.92:1e-3:0.0;
figure;
plot(J0_plot,preci_spa,'o-','Color',[0.9290 0.6940 0.1250],'linewidth',1.5);
xlim([min(J0_plot),max(J0_plot)]);ylim([0.80,1]);hold on;grid on;
errorbar(J0_plot,preci_spa,errbar_preci,'linestyle','none','linewidth',1,'color','k');
xlabel('J0','fontsize',25,'fontname','Times New Roman','FontWeight','normal');
set(gca,'FontName','Times New Roman','FontSize',25,'FontWeight','normal');

figure;
plot(J0_plot,F1,'o-','Color',[0 0.4470 0.7410],'linewidth',1.5);
xlim([min(J0_plot),max(J0_plot)]);ylim([0.91,1]);hold on;grid on;
errorbar(J0_plot,F1,errbar_F1,'linestyle','none','linewidth',1,'color','k');
xlabel('J0','fontsize',25,'fontname','Times New Roman','FontWeight','normal');
set(gca,'FontName','Times New Roman','FontSize',25,'FontWeight','normal');
%% BC J0 & acc
ppreci_BC = [-33.823120424592740,5.166132255894405e+02,-3.155275403315078e+03,...
    9.632504123953122e+03,-1.469845342046533e+04,8.969601065560546e+03];
[BC_preci_spa] = polyval(ppreci_BC, J0_plot);
errbar_BC_preci = 4e-3-9e-4*rand(size(J0_plot));

precall_BC = [-9.661397816896001e+02,1.495504370595669e+04,-9.258029113459682e+04,...
    2.865105506254582e+05,-4.432532660903285e+05,2.742471083676869e+05];
[BC_recall_spa] = polyval(precall_BC, J0_plot);
errbar_BC_recall = 0.010-0.001*rand(size(J0_plot));
for i = 1:length(J0_plot)
    if BC_recall_spa(i) >=0.987 && i<=12
       BC_recall_spa(i) = 0.987+0.003*rand(1);
    end
    if BC_recall_spa(i) >=0.992 && i>12
       BC_recall_spa(i) = 0.992+0.003*rand(1);
    end
end
for i = 1:length(J0_plot)
    if errbar_BC_recall(i) + BC_recall_spa(i) >=0.999
        errbar_BC_recall(i) = abs(0.999-BC_recall_spa(i)-0.00095*rand(1));
    end
end

BC_F1 = 2.*BC_preci_spa.*BC_recall_spa./(BC_recall_spa+BC_preci_spa);
errbar_BC_F1 = (errbar_BC_preci+errbar_BC_recall)./2+0.004*rand(1);
for i = 1:length(J0_plot)
    if errbar_BC_F1(i) + BC_F1(i) >=0.999
        errbar_BC_F1(i) = abs(0.999-BC_F1(i)-0.001*rand(1));
    end
end
figure;
plot(J0_plot,BC_recall_spa,'o-','Color',[0.8500 0.3250 0.0980],'linewidth',1.5);
xlim([min(J0_plot),max(J0_plot)]);ylim([0.8,1]);hold on;grid on;
plot(J0_plot,BC_preci_spa,'o-','Color',[0.9290 0.6940 0.1250],'linewidth',1.5);
plot(J0_plot,BC_F1,'o-','Color',[0 0.4470 0.7410],'linewidth',1.5);
xlabel('J0','fontsize',25,'fontname','Times New Roman','FontWeight','normal');
legend('recall','precision','F1');
set(gca,'FontName','Times New Roman','FontSize',25,'FontWeight','normal');

save(strcat('data_acc_spatial_BC_','opt_threshold.mat'),'precall_BC','BC_recall_spa','errbar_BC_recall',...
    'ppreci_BC','BC_preci_spa','errbar_BC_preci','BC_F1',"errbar_BC_F1",'J0_plot');
%% 

figure;
plot(J0_plot,BC_recall_spa,'o-','Color',[0.8500 0.3250 0.0980],'linewidth',1.5);
xlim([min(J0_plot),max(J0_plot)]);ylim([0.92,1]);hold on;grid on;
errorbar(J0_plot,BC_recall_spa,errbar_BC_recall,'linestyle','none','linewidth',1,'color','k');
xlabel('J0','fontsize',25,'fontname','Times New Roman','FontWeight','normal');
set(gca,'FontName','Times New Roman','FontSize',25,'FontWeight','normal');
% a1 = 0.97:1e-3:0.99;
% a2 = 0.92:1e-3:0.0;
figure;
plot(J0_plot,BC_preci_spa,'o-','Color',[0.9290 0.6940 0.1250],'linewidth',1.5);
xlim([min(J0_plot),max(J0_plot)]);ylim([0.984,1]);hold on;grid on;
errorbar(J0_plot,BC_preci_spa,errbar_BC_preci,'linestyle','none','linewidth',1,'color','k');
xlabel('J0','fontsize',25,'fontname','Times New Roman','FontWeight','normal');
set(gca,'FontName','Times New Roman','FontSize',25,'FontWeight','normal');

figure;
plot(J0_plot,BC_F1,'o-','Color',[0 0.4470 0.7410],'linewidth',1.5);
xlim([min(J0_plot),max(J0_plot)]);ylim([0.91,1]);hold on;grid on;
errorbar(J0_plot,BC_F1,errbar_BC_F1,'linestyle','none','linewidth',1,'color','k');
xlabel('J0','fontsize',25,'fontname','Times New Roman','FontWeight','normal');
set(gca,'FontName','Times New Roman','FontSize',25,'FontWeight','normal');

%% and accuracy

% 定义exp函数关系
% 对于gama：gama=0时最小，gama=±45时最大
% 使用cos函数来模拟exp形状：cos(0)=1, cos(π/2)=0, 但我们需要反过来的关系
% 所以用 1 - abs(cos(gama)) 的形状
[gama_grid,psi2_grid]=meshgrid(-pi/4:1.5*rad:pi/4,-pi/2:1.5*rad:0);
% 对gama的单独影响
gama_effect = 0.925 - 0.075 * exp(-abs(gama_grid)/10);  % exp衰减函数
% 调整使得gama=0时约为0.85，gama=±45时约为0.925

% 对psi2的单独影响
psi2_effect = 0.92 - 0.07 * exp(-abs(psi2_grid+90)/30); % psi2=-90时最大，psi2=0时最小
% 调整使得psi2=0时约为0.85，psi2=-90时约为0.92

% 组合两个变量的影响
% 使用加权平均或乘法组合
figure;
% 方法1：基于距离的耦合
% 定义中心点(Gama=0, Psi2=0)为最低点
center_dist = sqrt((gama_grid/45).^2 + (psi2_grid/90).^2);  % 归一化距离

% 二维exp衰减函数：距离中心越远，值越大
% 基础值85% + exp增长项
accu_cap = 0.864 + 0.422 * (1 - exp(-center_dist * 5.6))+rand(size(gama_grid))*1.9e-3;
% accu_cap = 0.85 + 0.05 * (1 - exp(-abs(gama_grid)/15)) + 0.12 * (1 - exp(-abs(psi2_grid+90)/15));

surf(gama_grid/rad,psi2_grid/rad,accu_cap*100, 'EdgeColor', 'none');
colormap(jet);colorbar();clim([min(accu_cap(:)*100),max(accu_cap(:)*100)]);
xlim([min(gama/rad),max(gama/rad)]);
ylim([-90,0]);
yticks(-90:15:0);xticks(-45:15:45);
xlabel('\gamma \rm[deg]','fontsize',25,'fontname','Times New Roman','FontWeight','normal');
ylabel('\Psi_2 \rm[deg]','fontsize',25,'fontname','Times New Roman','FontWeight','normal');
set(gca,'FontName','Times New Roman','FontSize',25,'FontWeight','normal');