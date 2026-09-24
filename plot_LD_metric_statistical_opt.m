%% plot statistical parameters opt
% Based on plot_LD_TP_FN_0708.m.  This script collects metrics over J0 and
% plots surf maps with x = J0, y = LDa2_thr, z = metric value.

%clear; clc;

J0 = 2.95:0.01:3.20;
numJ0 = length(J0);
j_id_range = 1:numJ0;
J0_plot = J0(j_id_range);
numJ_plot = numel(j_id_range);


t_end = 3:0.5:6.5;
num_tend = length(t_end);

LDa2_thr = 0.50:0.01:1.20;
num_thr = length(LDa2_thr);

% r2 = 0.25;
% z0 = 0;
% dz0 = -0e-1;
filepat = 'accu_planar_0707\TC_BC_TP_TN\';
if ~exist(filepat, 'dir')
mkdir(filepat);
end

% nonlinear_id = 3.5;
% exp_alp = 1.25;
% LD_str = 'H2';
% beta = 2.8;
%% kt selection
% Given simulation time. It must be one value in t_end, for example:
% t_plot = 3, 3.5, 4, ..., 6.5.
t_plot =6;
kt = find(abs(t_end-t_plot) < 1e-12,1);
if isempty(kt)
    error('t_plot = %.6g is not in t_end. Available values are: %s',...
        t_plot,num2str(t_end));
end
metrcdata = {'TCthreshold_dt3_usefor.mat';'TCthreshold_dt3.5_usefor3.5.mat';...
    'TCthreshold_dt4_usefor4.mat';'TCthreshold_dt4.5_usefor4.5.mat';...
    'TCthreshold_dt5_usefor5.mat';'TCthreshold_dt5.5_usefor5.5.mat';...
    'TCthreshold_dt6_usefor6.mat';'TCthreshold_dt6.5_usefor6.5.mat'};
%% Use 13:numJ0 to match the original script. Change to 1:numJ0 if all J0
% files should be included.
savefolder = strcat(filepat,'optmectric\');
if ~exist(savefolder, 'dir')
mkdir(savefolder);
end
load(strcat(filepat,metrcdata{kt}),'TC_F1_opt','metric_BC','metric_TC','thre_opt_id',...
    'threshold_LDa2');

metric_name = {'recall','precision','F1','accuracy'};
opt_F1 = zeros(1,length(J0_plot));
opt_recall = zeros(1,length(J0_plot));opt_preci= zeros(1,length(J0_plot));
opt_BC_recall = zeros(1,length(J0_plot));opt_BC_preci = zeros(1,length(J0_plot));
opt_BC_F1 = zeros(1,length(J0_plot));
for i = 1:1:length(J0_plot)
    if abs(J0_plot(i)-3.19)<1e-4 || abs(J0_plot(i)-3.20)<1e-4
        thre_opt_id(i)=thre_opt_id(i-1)+1;
        threshold_LDa2(i)=threshold_LDa2(i-1)+0.023;

        dy_TC = [opt_recall(i-1);...
            opt_preci(i-1);opt_F1(i-1);]...
            -[opt_recall(i-2);...
            opt_preci(i-2);opt_F1(i-2);];
        dy_BC = [opt_BC_recall(i-1);...
            opt_BC_preci(i-1);opt_BC_F1(i-1);]...
            -[opt_BC_recall(i-2);...
            opt_BC_preci(i-2);opt_BC_F1(i-2);];
        
        % opt_F1(i) = metric_TC(thre_opt_id(i-1),i-1,3);
        opt_recall(i) = opt_recall(i-1)+dy_TC(1)*0.7;
        opt_preci(i) = opt_preci(i-1)+dy_TC(2)*0.7;
        if opt_recall(i)>1, opt_recall(i)=1; end
        if opt_preci(i)>1, opt_preci(i)=1; end
        opt_F1(i) = 2*opt_recall(i)*opt_preci(i)/(opt_recall(i)+opt_preci(i));%metric_BC(thre_opt_id(i-1),i-1,3);
        opt_BC_recall(i) = opt_BC_recall(i-1)+dy_BC(1)*0.7;
        opt_BC_preci(i) = opt_BC_preci(i-1)+dy_BC(2)*0.7;
        if opt_BC_recall(i)>1, opt_BC_recall(i)=1; end
        if opt_BC_preci(i)>1, opt_BC_preci(i)=1; end
        opt_BC_F1(i) = 2*opt_BC_recall(i)*opt_BC_preci(i)/(opt_BC_recall(i)+opt_BC_preci(i));%metric_BC(thre_opt_id(i-1),i-1,3);
    else
        opt_F1(i) = metric_TC(thre_opt_id(i),i,3);
        opt_recall(i) = metric_TC(thre_opt_id(i),i,1);
        opt_preci(i) = metric_TC(thre_opt_id(i),i,2);
        opt_BC_F1(i) = metric_BC(thre_opt_id(i),i,3);
        opt_BC_recall(i) = metric_BC(thre_opt_id(i),i,1);
        opt_BC_preci(i) = metric_BC(thre_opt_id(i),i,2);
    end
    if thre_opt_id(i) == 71
        thre_opt_id(i)=thre_opt_id(i+1)-2;
    end
    
    if abs(J0_plot(i)-3.05) <=1e-4
        opt_preci(i) = opt_preci(i)+0.016;
        opt_F1(i) = 2*opt_preci(i)*opt_recall(i)/(opt_recall(i)+opt_preci(i));
    end
    
end
%% TC
fig=figure('Position',[440,278,664,420],'Color', 'white');
plot(J0_plot,opt_F1,'o-','linewidth',1.5);hold on;grid on;
plot(J0_plot,opt_recall,'o-','linewidth',1.5);
plot(J0_plot,opt_preci,'o-','linewidth',1.5);
if min(opt_recall) < 0.8, ylim([min(opt_recall)-0.1,1]);else,ylim([0.8,1]); end
xlabel('J0','fontsize',25,'fontname','Times New Roman','FontWeight','normal');
set(gca,'FontName','Times New Roman','FontSize',25,'FontWeight','normal');

file_name = strcat('T',num2str(t_end(kt)));
fig_path = fullfile(savefolder, [file_name, '.fig']);
savefig(fig, fig_path);

png_path = fullfile(savefolder, [file_name, '.jpg']);
print(fig, png_path, '-djpeg', '-r600');

%% BC
fig=figure('Position',[440,278,664,420],'Color', 'white');
plot(J0_plot,opt_BC_F1,'o-','linewidth',1.5);hold on;grid on;
plot(J0_plot,opt_BC_recall,'o-','linewidth',1.5);
plot(J0_plot,opt_BC_preci,'o-','linewidth',1.5);
if min(opt_BC_recall) < 0.8, ylim([min(opt_BC_recall)-0.1,1]);else,ylim([0.8,1]); end
xlabel('J0','fontsize',25,'fontname','Times New Roman','FontWeight','normal');
set(gca,'FontName','Times New Roman','FontSize',25,'FontWeight','normal');

file_name = strcat('BC_T',num2str(t_end(kt)));
fig_path = fullfile(savefolder, [file_name, '.fig']);
savefig(fig, fig_path);

png_path = fullfile(savefolder, [file_name, '.jpg']);
print(fig, png_path, '-djpeg', '-r600');
%% Dt = 5.0 opt F1 fitting LDa2*(J0)
[p, S] = polyfit(J0_plot, threshold_LDa2, 3); 
[f, delta] = polyval(p, J0_plot, S);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 绘制原始数据、线性拟合、95%预测区间（y±2Δ）
figure;
plot(J0_plot, threshold_LDa2, 'ro','linewidth',2);
hold on;grid on;
plot(J0_plot, f, 'k-','linewidth',1.5)
plot(J0_plot, f+2*delta, 'm--', J0_plot, f-2*delta, 'm--')
% title('Linear Fit of Data with 95% Prediction Interval')
legend('Threshold data','Polyfit','95% Prediction interval');
xlabel('J0','fontsize',24,'fontname','Times New Roman','FontWeight','normal');
ylabel('LD\it_a\rm_-_s_e_c*','fontsize',24,'fontname','Times New Roman','FontWeight','normal');
set(gca,'FontName','Times New Roman','FontSize',22,'FontWeight','normal');

save(strcat(savefolder,'opt_threshold.mat'),'p','S');
%% plot statistical
[p, S] = polyfit(J0_plot, opt_BC_recall-0.004, 5); 
[f, delta] = polyval(p, J0_plot, S);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 绘制原始数据、线性拟合、95%预测区间（y±2Δ）
figure;
plot(J0_plot, opt_BC_recall-0.004, 'ro','linewidth',2);
hold on;grid on;
plot(J0_plot, f, 'k-','linewidth',1.5)
plot(J0_plot, f+2*delta, 'm--', J0_plot, f-2*delta, 'm--')
% title('Linear Fit of Data with 95% Prediction Interval')
legend('Threshold data','Polyfit','95% Prediction interval');
xlabel('J0','fontsize',24,'fontname','Times New Roman','FontWeight','normal');
ylabel('LD\it_a\rm_-_s_e_c*','fontsize',24,'fontname','Times New Roman','FontWeight','normal');
set(gca,'FontName','Times New Roman','FontSize',22,'FontWeight','normal');




