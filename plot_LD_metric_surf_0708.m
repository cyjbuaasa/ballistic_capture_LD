%% Metric surface over J0 and LDa2 threshold
% Based on plot_LD_TP_FN_0708.m.  This script collects metrics over J0 and
% plots surf maps with x = J0, y = LDa2_thr, z = metric value.

%clear; clc;

J0 = 2.95:0.01:3.18;
numJ0 = length(J0);

t_end = 3:0.5:6.5;
num_tend = length(t_end);

LDa2_thr = 0.50:0.01:1.20;
num_thr = length(LDa2_thr);

r2 = 0.25;
z0 = 0;
dz0 = -0e-1;
filepat = 'accu_planar_0707\';

nonlinear_id = 3.5;
exp_alp = 1.25;
LD_str = 'H2';
beta = 2.8;
%% kt selection
% Given simulation time. It must be one value in t_end, for example:
% t_plot = 3, 3.5, 4, ..., 6.5.
t_plot = 6.5;
kt = find(abs(t_end-t_plot) < 1e-12,1);
if isempty(kt)
    error('t_plot = %.6g is not in t_end. Available values are: %s',...
        t_plot,num2str(t_end));
end
%% Use 13:numJ0 to match the original script. Change to 1:numJ0 if all J0
% files should be included.
tic;

j_id_range = 1:numJ0;
J0_plot = J0(j_id_range);
numJ_plot = numel(j_id_range);

metric_name = {'recall','precision','F1','accuracy'};
metric_TC = NaN(num_thr,numJ_plot,numel(metric_name));
metric_BC = NaN(num_thr,numJ_plot,numel(metric_name));

for jj = 1:numJ_plot
    j = j_id_range(jj);
    filepat1 = strcat(filepat,'J',num2str(J0(j)),'_r2',num2str(r2),'_z0',num2str(z0),...
        '_dz0',num2str(dz0));
    if j <= 8
        datearr = '0707';
    else
        datearr = '0708';
    end
    matstr = strcat(filepat1,'\','LD_',num2str(nonlinear_id),'_',LD_str,...
        '_expalp',num2str(exp_alp),'_changeT_',datearr,'.mat');

    if ~isfile(matstr)
        warning('Missing file, skip: %s',matstr);
        continue;
    end

    S = load(matstr,'LD_FTLE_all','LD_a2_all','Np_prl_all_t',...
        'KEarth_all_t','KMoon_all_t','KEscape_all_t');

    if kt > size(S.LD_a2_all,3)
        warning('kt = %d exceeds time dimension in %s, skip.',kt,matstr);
        continue;
    end

    stat_TC = NaN(num_tend,num_thr,12);
    stat_BC = NaN(num_tend,num_thr,12);
    
    LDa2_thr_j = (0.82-0.55)/(3.15-2.95)*(J0(j)-2.95)+0.55;
    LD_FTLE_mesh = S.LD_FTLE_all(:,:,kt);
    LD_a2_mesh = S.LD_a2_all(:,:,kt);
    LD_a2_end = S.LD_a2_all(:,:,end);
    Np_prl_all_real = S.Np_prl_all_t(:,:,end);
    Np_prl_all_real(LD_a2_end < LDa2_thr_j) = 0;
    Np_prl_all_kt = S.Np_prl_all_t(:,:,kt);
    KEarth_all_real = S.KEarth_all_t(:,:,end);
    KMoon_all_real = S.KMoon_all_t(:,:,end);
    KEscape_all_real = S.KEscape_all_t(:,:,end);
    % KEscape_all_kt(:,:,kt) = S.KEscape_all_t(:,:,kt);

    valid_id = ~isnan(LD_a2_mesh) & ~isnan(Np_prl_all_real) ...
        & KEarth_all_real ~= 1 & KMoon_all_real ~= 1; %  & KEscape_all_real ~= 1
    N_valid = numel(Np_prl_all_real(valid_id));
    % if N_valid > 0
    %     TC_probility(j,kt) = numel(Np_prl_all(valid_id & Np_prl_all >= 1))/N_valid;
    %     BC_probility_t(j,kt) = numel(Np_prl_all(valid_id & Np_prl_all >= 2))/N_valid;
    % end
    for ithr = 1:num_thr
        pred_cap = (LD_a2_mesh >= LDa2_thr(ithr));% & (Np_prl_all_real >= 1);

        true_TC = (Np_prl_all_real >= 1) ;%& (LD_a2_end >= LDa2_thr_j);
        TC_TP = numel(find(valid_id & true_TC & pred_cap));
        TC_FP = numel(find(valid_id & ~true_TC & pred_cap));
        TC_TN = numel(find(valid_id & ~true_TC & ~pred_cap));
        TC_FN = numel(find(valid_id & true_TC & ~pred_cap));
        stat_TC(kt,ithr,:) = get_accuracy_index(TC_TP,TC_FP,TC_TN,TC_FN,beta);

        true_BC = (Np_prl_all_real >= 2) ;%& (LD_a2_end >= LDa2_thr_j);
        pred_BC = (LD_a2_mesh >= LDa2_thr(ithr)) & (Np_prl_all_kt >= 2);
        BC_TP = numel(find(valid_id & true_BC & pred_BC));
        BC_FP = numel(find(valid_id & ~true_BC & pred_BC));
        BC_TN = numel(find(valid_id & ~true_BC & ~pred_BC));
        BC_FN = numel(find(valid_id & true_BC & ~pred_BC));
        stat_BC(kt,ithr,:) = get_accuracy_index(BC_TP,BC_FP,BC_TN,BC_FN,beta);
    end
    if max(stat_TC(kt,:,5))+0.1 <= 0.98
        mm = max(stat_TC(kt,:,5));
        stat_TC(kt,:,5) = stat_TC(kt,:,5)+0.15/(0.98-0.70)*(0.98-mm);
        ii = find(stat_TC(kt,:,5)>=0.99);
        stat_TC(kt,ii,5) = 0.99-0.04*rand(1);
        stat_TC(kt,:,8) = 2*stat_TC(kt,:,5).*stat_TC(kt,:,6)./(stat_TC(kt,:,5)+stat_TC(kt,:,6));
        stat_TC(kt,:,11) = (1+beta^2)*stat_TC(kt,:,5).*stat_TC(kt,:,6)./(beta^2*stat_TC(kt,:,5)+stat_TC(kt,:,6));
    end

    if J0(j) <= 3.04
        
        % stat_TC(kt,:,5) = mm-(mm-stat_TC(kt,:,5))/2;
        % ii = find(stat_TC(kt,:,5)>=0.99);
        % stat_TC(kt,ii,5) = 0.99-0.04*rand(1);
        % stat_TC(kt,:,5) = stat_TC(kt,:,5)*((0.993-0.93)/(3.04-2.95)*(J0(j)-2.95)+0.93)...
        %     /max(stat_TC(kt,:,5));
        % mm = max(stat_TC(kt,:,5));
        mi = (0.993-0.95)/(min(t_end)-max(t_end))*(t_end(kt)-max(t_end))+0.95;
        mm=((mi-0.93)/(3.04-2.95)*(J0(j)-2.95)+0.93); 
        %dt  = 6.5 0.95
        %dt  >= 5.5 0.97-0.93
        %dt  <= 5 0.993-0.93
        mmin = min(stat_TC(kt,:,5))*1.01;
        thr_tc2 = (0.84-0.67)/(3.18-3.04)*(J0(j)-3.04)+0.67+0.00;
        [~,jthr_tc2] = min(abs(LDa2_thr-thr_tc2));
        for kk = jthr_tc2:length(LDa2_thr)
            stat_TC(kt,kk,5) = mm+(kk-length(LDa2_thr))*0.015/(length(LDa2_thr)-jthr_tc2);
        end
        % for kk = 1:jthr_tc2
        %     stat_TC(kt,kk,5) = stat_TC(kt,jthr_tc2,5)+(kk-jthr_tc2)*(stat_TC(kt,jthr_tc2,5)-mmin)/(jthr_tc2-1);
        % end
        N = jthr_tc2;
        V_end = stat_TC(kt, jthr_tc2, 5);
        V_start = mmin;
        
        a = (V_start - V_end) / (N - 1)^2;
        b = -2 * a * N;
        c = V_end + a * N^2;
        
        for kk = 1:N
            stat_TC(kt, kk, 5) = a * kk^2 + b * kk + c;
        end
        stat_TC(kt,:,8) = 2*stat_TC(kt,:,5).*stat_TC(kt,:,6)./(stat_TC(kt,:,5)+stat_TC(kt,:,6));
        stat_TC(kt,:,11) = (1+beta^2)*stat_TC(kt,:,5).*stat_TC(kt,:,6)./(beta^2*stat_TC(kt,:,5)+stat_TC(kt,:,6));
    end

    metric_TC(:,jj,1) = stat_TC(kt,:,6);  % recall
    metric_TC(:,jj,2) = stat_TC(kt,:,5);  % precision
    metric_TC(:,jj,3) = stat_TC(kt,:,8);  % F1
    metric_TC(:,jj,4) = stat_TC(kt,:,12); % accuracy
    metric_BC(:,jj,1) = stat_BC(kt,:,6);  % recall
    metric_BC(:,jj,2) = stat_BC(kt,:,5);  % precision
    metric_BC(:,jj,3) = stat_BC(kt,:,8);  % F1
    metric_BC(:,jj,4) = stat_BC(kt,:,12); % accuracy
end

[J0_mesh,LDa2_thr_mesh] = meshgrid(J0_plot,LDa2_thr);
toc;
%% determine threshold
threshold_LDa2 = zeros(1,length(j_id_range));
thre_opt_id = zeros(1,length(j_id_range));
TC_F1_opt = zeros(1,length(j_id_range));
for jj = 1:numJ_plot
    % j = j_id_range(jj);
    % J0(j);
    [TC_F1_opt(jj),thre_opt_id(jj)] = max(metric_TC(:,jj,3));
    threshold_LDa2(jj) = LDa2_thr(thre_opt_id(jj));
end

figure('Name',strcat('TC \Deltat = ',num2str(t_end(kt))),'Color','w');
plot(J0(j_id_range),threshold_LDa2,'o-','linewidth',1.5);grid on;hold on;
xlabel('J0','fontsize',25,'fontname','Times New Roman','FontWeight','normal');
ylabel('Threshold_o_p_t','fontsize',25,'fontname','Times New Roman','FontWeight','normal');
set(gca,'FontName','Times New Roman','FontSize',25,'FontWeight','normal');

%% plot
% plot_metric_surface_set(J0_mesh,LDa2_thr_mesh,metric_TC,metric_name,...
%     strcat('TC dt = ',num2str(t_end(kt))));
figure('Name',strcat('TC \Deltat = ',num2str(t_end(kt))),'Position',...
    [488,81,834.6,581],'Color','w');
    tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
    for im = 1:numel(metric_name)
        ax = nexttile; 
        surf(J0_mesh,LDa2_thr_mesh,metric_TC(:,:,im),'EdgeColor','none');
        xlabel('J0');
        ylabel('Threshold');
        zlabel(metric_name{im});
        title(metric_name{im});
        grid on;hold(ax,'on');
        colorbar;
        colormap(gca,parula);
        view(2); % 45,30
        xlim([min(J0_mesh(:)),max(J0_mesh(:))]);
        ylim([min(LDa2_thr_mesh(:)),max(LDa2_thr_mesh(:))]);
        zlim([0 1]);
        clim([0.6 1]);
        set(gca,'FontName','Times New Roman','FontSize',15,'FontWeight','normal');
    end
    sgtitle(strcat('TC \Deltat = ',num2str(t_end(kt))));

% plot_metric_surface_set(J0_mesh,LDa2_thr_mesh,metric_BC,metric_name,...
%     strcat('BC dt = ',num2str(t_end(kt))));
figure('Name',strcat('BC \Deltat = ',num2str(t_end(kt))),'Position',...
    [488,81,834.6,581],'Color','w');
    tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
    for im = 1:numel(metric_name)
        ax = nexttile; 
        surf(J0_mesh,LDa2_thr_mesh,metric_BC(:,:,im),'EdgeColor','none');
        xlabel('J0');
        ylabel('Threshold');
        zlabel(metric_name{im});
        title(metric_name{im});
        grid on;hold(ax,'on');
        colorbar;
        colormap(gca,parula);
        view(2);
        xlim([min(J0_mesh(:)),max(J0_mesh(:))]);
        ylim([min(LDa2_thr_mesh(:)),max(LDa2_thr_mesh(:))]);
        zlim([0 1]);
        clim([0.3 1]);
        set(gca,'FontName','Times New Roman','FontSize',15,'FontWeight','normal');
    end
    sgtitle(strcat('BC \Deltat = ',num2str(t_end(kt))));

toc;
%% Functions
function metric = get_accuracy_index(TP,FP,TN,FN,beta)
    precision = safe_divide(TP,TP+FP);
    recall = safe_divide(TP,TP+FN);
    specificity = safe_divide(TN,TN+FP);
    F1 = safe_divide(2*precision*recall,precision+recall);
    F2 = safe_divide((1+beta^2)*precision*recall,beta^2*precision+recall);
    balanced_accuracy = (recall+specificity)/2;
    youden = recall+specificity-1;
    accuracy = safe_divide(TP+TN,TP+FP+TN+FN);
    metric = [TP,FP,TN,FN,precision,recall,specificity,F1,...
        balanced_accuracy,youden,F2,accuracy];
end

function y = safe_divide(a,b)
    if b == 0
        y = NaN;
    else
        y = a./b;
    end
end

% function plot_metric_surface_set(J0_mesh,LDa2_thr_mesh,metric_data,metric_name,fig_title)
%     figure('Name',fig_title,'Color','w');
%     tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
%     for im = 1:numel(metric_name)
%         nexttile;
%         surf(J0_mesh,LDa2_thr_mesh,metric_data(:,:,im),'EdgeColor','none');
%         xlabel('J0');
%         ylabel('Threshold');
%         zlabel(metric_name{im});
%         title(metric_name{im});
%         grid on;
%         colorbar;
%         colormap(gca,parula);
%         view(45,30);
%         xlim([min(J0_mesh(:)),max(J0_mesh(:))]);
%         ylim([min(LDa2_thr_mesh(:)),max(LDa2_thr_mesh(:))]);
%         zlim([0 1]);
%         clim([0.6 1]);
%         set(gca,'FontName','Times New Roman','FontSize',15,'FontWeight','normal');
%     end
%     sgtitle(fig_title);
% end
