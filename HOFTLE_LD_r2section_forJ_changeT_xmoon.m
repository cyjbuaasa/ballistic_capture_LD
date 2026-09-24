%% multibody orbital dynamics
clc;clear 
%% parameter
addpath('D:\matlab_program\cyj_program\CRTBP\TO_evolution');
basic_CRTBP_etal;
rmpath('D:\matlab_program\cyj_program\CRTBP\TO_evolution');
% 
mu_EMS = parameter.u.EM_Sun;
day = 86400;rad = pi/180;
mu = parameter.u.Moon;
TU = parameter.t.TU;
DU = parameter.d.DU;
DUS = parameter.d.as*DU;
TUS = parameter.t.TUS;
VU = parameter.v.VU;
VUS = parameter.v.VUS;
RE = parameter.d.REarth;
RM = parameter.d.RMoon;
m_s = parameter.m.Sun;       %%
rho_sun = parameter.d.as;   %%
rEM_hill = parameter.d.rEM_hill;
rm_hill = parameter.d.rm_hill;
ome_sun = parameter.t.ome_sun;
%% Unit dimensionless and Lagrangian points
start_time = cputime;
T=2*pi*TU; % 
LL=get_LL(mu,mu_EMS,TU);
option_LL = optimoptions('fsolve','StepTolerance',1e-12,'FunctionTolerance',1e-12,'OptimalityTolerance',1e-12,'Display','none');

%% add STT & tensor path
addpath('D:\matlab_program\cyj_program\CRTBP\STT');
% addpath(genpath('F:\cyj_CRTBP\tensor_toolbox-v3.8'));
%% de2/dt equation
tic;
% H0 = -847.500; % Initial Hamiltonian in BCRFBP
J0 = 2.95:0.01:3.18;
numJ0 = length(J0);
eps=1e-3;
theta_s = pi/2;
t0 = 0;t_end = 3:0.5:6.5;
num_tend = length(t_end);
LDa2_thr = 0.50:0.01:1.20;
num_thr = length(LDa2_thr);
min_recall_TC = 0.90;
min_recall_BC = 0.90;

% z component setting
r2 = 0.25;
z0 = 0;dz0 = -0e-1;
r_sec = 1.50*66180/DU;

%% ini parameters
tic
filepat = 'accu_planar_0709_xmoon\';
if ~exist(filepat, 'dir')
   mkdir(filepat);
end
error_log_dir = fullfile(filepat, ['parfor_error_log_', datestr(now,'yyyymmdd_HHMMSS')]);
if ~exist(error_log_dir, 'dir')
   mkdir(error_log_dir);
end

theta = asin(z0/r2);
alp = [0:0.5:360]*rad;
% ballistic angle
Psi = [0:0.5:360]*rad;

num_all = length(alp)*length(Psi);
[alp_mesh,Psi_mesh] = meshgrid(alp,Psi);

options_LD = odeset('RelTol',1e-9,'AbsTol',1e-9,'Events',@only_hit_LD);
% v-planar
v_plan = zeros(size(alp_mesh));
v_vert = dz0;

dx_eps = 1e-4;
nonlinear_id = 3.5;
% 1--FTLE; 2--HOFTLE(.1-STT; .2-DSTT; .3-TDSTT); 
% 3--LD(.1-E2; .2-H2; .3-v2; .4-a2; .5-H2|E2|r2|v2|onlyv(v-norm) & a2)
exp_alp = 1.25;
keta = 0;
m_STT = 2; % m order STT
dyn_str = 'EM'; % 'EM' 'diff_sbcm'
LD_str = 'H2';
beta = 2;

BC_probility = zeros(numJ0,1);
TC_probility = zeros(numJ0,num_tend);
BC_probility_t = zeros(numJ0,num_tend);
best_thr_TC = NaN(numJ0,num_tend);
best_thr_BC = NaN(numJ0,num_tend);
best_metric_TC = NaN(numJ0,num_tend,11);
best_metric_BC = NaN(numJ0,num_tend,11);
metric_name = {'TP','FP','TN','FN','precision','recall',...
    'specificity','F1','balanced_accuracy','youden','F2'};
save_fig_each_tend = true;
%% accuracy
% True positive: It is a temporary capture, predicted to be a temporary capture;
% False positive: It is not a temporary capture, but is predicted to be temporarily captured;
% True negative: It is not a temporary capture, and it is predicted that it is not a temporary capture;
% False negative: It is a temporary capture, but it is predicted that it is not a temporary capture;
% temporary capture(TC): Np_prl >=1
% ballistic capture(BC): Np_prl >=2
% predicted-LD-based TC: LD >= LD_a_sec* (Np_prl >=1 && E_sec<=0 at prl epoch)
% predicted-LD-based BC: LD >= LD_a_sec* && Np_prl >=2 (&& E_sec<=0 at prl epoch)

% precision = TP/(TP+FP); 5
% recall = TP/(TP+FN); 6
% specificity = TN/(TN+FP); 7
% F1 = 2*precision*recall/(precision+recall); 8
% F2 = (1+beta^2)*precision*recall/(beta^2*precision+recall); 11
% balanced_accuracy = (recall+specificity)/2; 9
% youden = recall+specificity-1; 10
%% for J0
tic; 
for j = 1:numJ0
    LD_FTLE_flat = NaN(num_all,num_tend);
    LD_a2_flat = NaN(num_all,num_tend);
    Np_prl_flat = NaN(num_all,num_tend);
    KEarth_flat = NaN(num_all,num_tend);
    KMoon_flat = NaN(num_all,num_tend);
    KEscape_flat = NaN(num_all,num_tend);
    KXMoon_flat = NaN(num_all,num_tend);
    eve_struct_all = cell(length(alp),length(Psi));
    fprintf("J0 = %f\n",J0(j));
    parfor i = 1:num_all
        LD_FTLE_i = NaN(1,num_tend);
        LD_a2_i = NaN(1,num_tend);
        Np_prl_i = NaN(1,num_tend);
        KEarth_i = NaN(1,num_tend);
        KMoon_i = NaN(1,num_tend);
        KEscape_i = NaN(1,num_tend);
        KXMoon_i = NaN(1,num_tend);
        eve_struct = struct('tend',0,'xend',0,'eve',0);
        eve = 0;tend = 0;xend = 0;

        y0 = NaN(6,1); vv = NaN; sig = NaN;
        try
            r2_vec = r2*[cos(alp_mesh(i))*cos(theta);sin(alp_mesh(i))*cos(theta);sin(theta)];
            r1_vec = [1;0;0]+r2_vec;
            r1 = norm(r1_vec);

            vv = fsolve(@(v)(1*(v^2+v_vert^2)-((1-mu+r2_vec(1))^2+r2_vec(2)^2) ...
                -2*(1-mu)./r1-2*mu./r2-mu*(1-mu)+J0(j)),1e-2,option_LL);
            if vv >= 1e-3
                sig = alp_mesh(i) + Psi_mesh(i);
                y0 = [[r2_vec+[1-mu;0;0]];vv*cos(sig);vv*sin(sig);v_vert];
                if nonlinear_id == 3.5
                    [LD_FTLE_i,LD_a2_i,Np_prl_i,KEarth_i,KMoon_i,KEscape_i,KXMoon_i,tend,xend,eve] = ...
                        get_Dy_LD_cons_tend_xmoon_mex(t0,t_end,num_tend,y0,...
                            dyn_str,LD_str,eps,exp_alp,keta,r_sec);
                elseif nonlinear_id >1 && nonlinear_id < 3
                    [orbDEX,Np_prl,eve_struct,t,X,tend,xend,eve] = get_Dyindex(t0,t_end(end),y0,nonlinear_id,m_STT,...
                        dyn_str,eps,exp_alp,keta,dx_eps,LD_str,options_LD,mu);
                    LD_FTLE_i(end) = orbDEX;Np_prl_i(end) = Np_prl;
                else
                    [orbDEX,Np_prl,eve_struct,tend,xend,eve] = get_Dyindex_1_3_mex(t0,t_end(end),y0,nonlinear_id,m_STT,...
                        dyn_str,eps,exp_alp,keta,dx_eps,LD_str,mu);
                    LD_FTLE_i(end) = orbDEX(1);Np_prl_i(end) = Np_prl;
                end
                eve_struct = struct('tend',tend,'xend',xend,'eve',eve);
            end
        catch ME
            task = getCurrentTask();
            if isempty(task)
                worker_id = 0;
            else
                worker_id = task.ID;
            end
            [psi_idx, alp_idx] = ind2sub(size(alp_mesh), i);
            err_file = fullfile(error_log_dir, sprintf('fail_J%0.4f_i%06d_w%03d.txt', J0(j), i, worker_id));
            fid = fopen(err_file, 'w');
            if fid ~= -1
                fprintf(fid, 'Error in HOFTLE_LD_r2section_forJ_changeT_xmoon parfor sample\n');
                fprintf(fid, 'message: %s\nidentifier: %s\n', ME.message, ME.identifier);
                fprintf(fid, 'j=%d, J0=%.16g, i=%d, alp_idx=%d, psi_idx=%d, worker=%d\n', j, J0(j), i, alp_idx, psi_idx, worker_id);
                fprintf(fid, 'alp=%.16g rad, Psi=%.16g rad, sig=%.16g rad, vv=%.16g\n', alp_mesh(i), Psi_mesh(i), sig, vv);
                fprintf(fid, 'y0=[%.16g %.16g %.16g %.16g %.16g %.16g]\n', y0);
                for kk = 1:numel(ME.stack)
                    fprintf(fid, 'stack(%d): %s, line %d\n', kk, ME.stack(kk).name, ME.stack(kk).line);
                end
                fclose(fid);
            end
            fprintf(2, 'parfor failed: J0=%.6f, i=%d, alp=%.6f, Psi=%.6f, worker=%d. Log: %s\n', ...
                J0(j), i, alp_mesh(i), Psi_mesh(i), worker_id, err_file);
        end
        LD_FTLE_flat(i,:) = LD_FTLE_i;
        LD_a2_flat(i,:) = LD_a2_i;
        Np_prl_flat(i,:) = Np_prl_i;
        KEarth_flat(i,:) = KEarth_i;
        KMoon_flat(i,:) = KMoon_i;
        KEscape_flat(i,:) = KEscape_i;
        KXMoon_flat(i,:) = KXMoon_i;
        eve_struct_all{i} = eve_struct;
        if mod(i,1e4) == 0, fprintf("%d done\n",i); end
    end
    toc;

    filepat1 = strcat(filepat,'J',num2str(J0(j)),'_r2',num2str(r2),'_z0',num2str(z0),...
        '_dz0',num2str(dz0));
    
    if ~exist(filepat1, 'dir')
       mkdir(filepat1);
    end
    LD_FTLE_all = NaN(length(alp),length(Psi),num_tend);
    LD_a2_all = NaN(length(alp),length(Psi),num_tend);
    Np_prl_all_t = NaN(length(alp),length(Psi),num_tend);
    KEarth_all_t = NaN(length(alp),length(Psi),num_tend);
    KMoon_all_t = NaN(length(alp),length(Psi),num_tend);
    KEscape_all_t = NaN(length(alp),length(Psi),num_tend);
    KXMoon_all_t = NaN(length(alp),length(Psi),num_tend);
    stat_TC = NaN(num_tend,num_thr,11);
    stat_BC = NaN(num_tend,num_thr,11);
%%
    for kt = 1:num_tend
        LD_FTLE_mesh = reshape(LD_FTLE_flat(:,kt),size(alp_mesh));
        LD_a2_mesh = reshape(LD_a2_flat(:,kt),size(alp_mesh));
        Np_prl_all = reshape(Np_prl_flat(:,kt),size(alp_mesh));
        Np_prl_all_real = reshape(Np_prl_flat(:,end),size(alp_mesh));
        KEarth_all = reshape(KEarth_flat(:,kt),size(alp_mesh));
        KMoon_all = reshape(KMoon_flat(:,kt),size(alp_mesh));
        KEscape_all = reshape(KEscape_flat(:,kt),size(alp_mesh));
        % KEscape_all_real = reshape(KEscape_flat(:,end),size(alp_mesh));
        KXMoon_all = reshape(KXMoon_flat(:,kt),size(alp_mesh));
        % KXMoon_all_real = reshape(KXMoon_flat(:,end),size(alp_mesh));

        LD_FTLE_mesh(LD_FTLE_mesh == 0) = nan;
        LD_a2_mesh(LD_a2_mesh == 0) = nan;
        Np_prl_all(isnan(LD_a2_mesh)) = nan;
        % Np_prl_all(KEarth_all == 1) = nan;Np_prl_all(KMoon_all == 1) = nan;

        LD_FTLE_all(:,:,kt) = LD_FTLE_mesh;
        LD_a2_all(:,:,kt) = LD_a2_mesh;
        Np_prl_all_t(:,:,kt) = Np_prl_all;
        KEarth_all_t(:,:,kt) = KEarth_all;
        KMoon_all_t(:,:,kt) = KMoon_all;
        KEscape_all_t(:,:,kt) = KEscape_all;
        KXMoon_all_t(:,:,kt) = KXMoon_all;

        valid_id = ~isnan(LD_a2_mesh) & ~isnan(Np_prl_all_real) ...
            & KEarth_all ~= 1 & KMoon_all ~= 1;% ...
            % & KEscape_all ~= 1 & KEscape_all_real ~= 1;
        N_valid = numel(find(valid_id));
        if N_valid > 0
            % every tend calculating TC&BC probability
            TC_probility(j,kt) = numel(Np_prl_all(valid_id & Np_prl_all >= 1))/N_valid;
            BC_probility_t(j,kt) = numel(Np_prl_all(valid_id & Np_prl_all >= 2))/N_valid;
        end
        % for ithr = 1:num_thr
        %     pred_cap = (LD_a2_mesh >= LDa2_thr(ithr));% & (Np_prl_all_real >= 1);
        % 
        %     true_TC = Np_prl_all_real >= 1;
        %     TC_TP = numel(find(valid_id & true_TC & pred_cap));
        %     TC_FP = numel(find(valid_id & ~true_TC & pred_cap));
        %     TC_TN = numel(find(valid_id & ~true_TC & ~pred_cap));
        %     TC_FN = numel(find(valid_id & true_TC & ~pred_cap));
        %     stat_TC(kt,ithr,:) = get_accuracy_index(TC_TP,TC_FP,TC_TN,TC_FN,beta);
        % 
        %     true_BC = Np_prl_all_real >= 2;
        %     pred_BC = (LD_a2_mesh >= LDa2_thr(ithr)) & (Np_prl_all >= 2);
        %     BC_TP = numel(find(valid_id & true_BC & pred_BC));
        %     BC_FP = numel(find(valid_id & ~true_BC & pred_BC));
        %     BC_TN = numel(find(valid_id & ~true_BC & ~pred_BC));
        %     BC_FN = numel(find(valid_id & true_BC & ~pred_BC));
        %     stat_BC(kt,ithr,:) = get_accuracy_index(BC_TP,BC_FP,BC_TN,BC_FN,beta);
        % end
        % 
        % TC_recall = stat_TC(kt,:,6);TC_youden = stat_TC(kt,:,10);
        % TC_F1 = stat_TC(kt,:,8);TC_F2 = stat_TC(kt,:,11);
        % id_good = find(TC_recall >= min_recall_TC & ~isnan(TC_youden));
        % if ~isempty(id_good)
        %     [~,id_loc] = max(TC_youden(id_good));id_best = id_good(id_loc);
        % else
        %     [~,id_best] = max(TC_F1);
        % end
        % % [~,id_best] = max(TC_F2);
        % best_thr_TC(j,kt) = LDa2_thr(id_best);
        % best_metric_TC(j,kt,:) = stat_TC(kt,id_best,:);
        % 
        % BC_recall = stat_BC(kt,:,6);BC_youden = stat_BC(kt,:,10);
        % BC_F1 = stat_BC(kt,:,8); BC_F2 = stat_BC(kt,:,11);
        % id_good = find(BC_recall >= min_recall_BC & ~isnan(BC_youden));
        % if ~isempty(id_good)
        %     [~,id_loc] = max(BC_youden(id_good));id_best = id_good(id_loc);
        % else
        %     [~,id_best] = max(BC_F1);
        % end
        % % [~,id_best] = max(BC_F2);
        % best_thr_BC(j,kt) = LDa2_thr(id_best);
        % best_metric_BC(j,kt,:) = stat_BC(kt,id_best,:);

        if  kt == num_tend || save_fig_each_tend % && false
            filepat_t = fullfile(filepat1,strcat('T',num2str(t_end(kt))));
            if ~exist(filepat_t, 'dir')
                mkdir(filepat_t);
            end
            fig_LD = figure('Position', [46.60,192.2,949.60,560], 'Color', 'white');
            hAx = axes;
            LD_field=surf(hAx,alp_mesh/rad,Psi_mesh/rad,LD_FTLE_mesh, 'EdgeColor', 'none');
            colormap(hAx,jet);colorbar(hAx);
            xlim([min(alp/rad),max(alp/rad)]); ylim([min(Psi/rad),max(Psi/rad)]);
            xlabel('\alpha \rm[deg]','fontsize',25,'fontname','Times New Roman','FontWeight','normal');
            ylabel('\Psi \rm[deg]','fontsize',25,'fontname','Times New Roman','FontWeight','normal');
            set(gca,'FontName','Times New Roman','FontSize',25,'FontWeight','normal');
            set(LD_field,'ButtonDownFcn',{@mycall_map0709,hAx,...
                t0,t0+t_end(kt),alp_mesh,Psi_mesh,LD_FTLE_mesh,r2,theta,v_vert,...
                J0(j),nonlinear_id,m_STT,dyn_str,eps,exp_alp,keta,dx_eps,LD_str});
            view(2);
            figname_LD = fullfile(filepat_t,strcat('LD_',num2str(nonlinear_id),'_',LD_str,...
                '_T',num2str(t_end(kt)),'_expalp',num2str(exp_alp),'_',datestr(datetime('today'),'mmdd'),'.fig'));

            figa2 = figure('Position', [46.60,192.2,949.60,560], 'Color', 'white');
            surf(alp_mesh/rad,Psi_mesh/rad,LD_a2_mesh, 'EdgeColor', 'none');
            colormap(jet);colorbar;
            xlim([min(alp/rad),max(alp/rad)]); ylim([min(Psi/rad),max(Psi/rad)]);
            xlabel('\alpha \rm[deg]','fontsize',25,'fontname','Times New Roman','FontWeight','normal');
            ylabel('\Psi \rm[deg]','fontsize',25,'fontname','Times New Roman','FontWeight','normal');
            set(gca,'FontName','Times New Roman','FontSize',25,'FontWeight','normal');
            view(2);
            figname_a2 = fullfile(filepat_t,strcat('LDa2_T',num2str(t_end(kt)),'_',datestr(datetime('today'),'mmdd'),'.fig'));

            fig_Np = figure('Position', [46.60,192.2,949.60,560], 'Color', 'white');
            surf(alp_mesh/rad,Psi_mesh/rad,Np_prl_all, 'EdgeColor', 'none');
            colormap(jet);colorbar();
            max_Np = max(Np_prl_all(:),[],'omitnan');
            if isempty(max_Np) || isnan(max_Np) || max_Np <= 0, max_Np = 1; end
            clim([0,max_Np]);
            xlim([min(alp/rad),max(alp/rad)]); ylim([min(Psi/rad),max(Psi/rad)]);
            xlabel('\alpha \rm[deg]','fontsize',25,'fontname','Times New Roman','FontWeight','normal');
            ylabel('\Psi \rm[deg]','fontsize',25,'fontname','Times New Roman','FontWeight','normal');
            set(gca,'FontName','Times New Roman','FontSize',25,'FontWeight','normal');
            view(2);
            figname_Np = fullfile(filepat_t,strcat('Np_',num2str(nonlinear_id),'_',LD_str,...
                '_T',num2str(t_end(kt)),'_expalp',num2str(exp_alp),'_',datestr(datetime('today'),'mmdd'),'.fig'));
            savefig(fig_LD, figname_LD);savefig(figa2, figname_a2);savefig(fig_Np, figname_Np);
            close(fig_LD);close(figa2);close(fig_Np);
        end
    end
    %%
    BC_probility(j) = BC_probility_t(j,end);

    % save(strcat(filepat1,'\','LD_',num2str(nonlinear_id),'_',LD_str,'_expalp',num2str(exp_alp),...
    %         '_changeT_',datestr(datetime('today'),'mmdd'),'.mat'),'LD_FTLE_all',...
    %         'LD_a2_all','Psi_mesh','alp_mesh','theta','BC_probility','TC_probility',...
    %         'BC_probility_t','alp','Psi','t0','t_end','nonlinear_id','m_STT','J0',...
    %         'Np_prl_all_t','r2','r_sec','v_vert','dyn_str','eps','exp_alp','dx_eps','LD_str','z0','dz0',...
    %         'KMoon_all_t','KEarth_all_t','KEscape_all_t','KXMoon_all_t','LDa2_thr','stat_TC','stat_BC','metric_name',...
    %         'best_thr_TC','best_thr_BC','best_metric_TC','best_metric_BC','-v7.3');
    save(strcat(filepat1,'\','LD_',num2str(nonlinear_id),'_',LD_str,'_expalp',num2str(exp_alp),...
            '_changeT_',datestr(datetime('today'),'mmdd'),'.mat'),'LD_FTLE_all',...
            'LD_a2_all','Psi_mesh','alp_mesh','theta','BC_probility','TC_probility',...
            'BC_probility_t','alp','Psi','t0','t_end','nonlinear_id','m_STT','J0',...
            'Np_prl_all_t','r2','r_sec','v_vert','dyn_str','eps','exp_alp','dx_eps','LD_str','z0','dz0',...
            'KMoon_all_t','KEarth_all_t','KEscape_all_t','KXMoon_all_t','LDa2_thr','-v7.3');
    % no eve_struct_all Too large!!!
    toc;
end

%% LD_threshold_changeT_summary
tic;
fit_order = 2;
fit_thr_TC = NaN(num_tend,fit_order+1);
fit_thr_BC = NaN(num_tend,fit_order+1);
fit_thr_TC_value = NaN(numJ0,num_tend);
fit_thr_BC_value = NaN(numJ0,num_tend);
J0_fine = 2.95:0.001:3.18;
fit_thr_TC_fine = NaN(length(J0_fine),num_tend);
fit_thr_BC_fine = NaN(length(J0_fine),num_tend);
for kt = 1:num_tend
    id_fit = ~isnan(best_thr_TC(:,kt));
    if numel(find(id_fit)) >= fit_order+1
        fit_thr_TC(kt,:) = polyfit(J0(id_fit),best_thr_TC(id_fit,kt)',fit_order);
        fit_thr_TC_value(:,kt) = polyval(fit_thr_TC(kt,:),J0);
        fit_thr_TC_fine(:,kt) = polyval(fit_thr_TC(kt,:),J0_fine);
    end
    id_fit = ~isnan(best_thr_BC(:,kt));
    if numel(find(id_fit)) >= fit_order+1
        fit_thr_BC(kt,:) = polyfit(J0(id_fit),best_thr_BC(id_fit,kt)',fit_order);
        fit_thr_BC_value(:,kt) = polyval(fit_thr_BC(kt,:),J0);
        fit_thr_BC_fine(:,kt) = polyval(fit_thr_BC(kt,:),J0_fine);
    end
end
save(strcat(filepat,'LD_threshold_changeT_xmoon_summary_',datestr(datetime('today'),'mmdd'),'.mat'),...
    'J0','t_end','LDa2_thr','metric_name','TC_probility','BC_probility_t',...
    'best_thr_TC','best_thr_BC','best_metric_TC','best_metric_BC',...
    'fit_order','J0_fine','fit_thr_TC','fit_thr_BC','fit_thr_TC_value','fit_thr_BC_value',...
    'fit_thr_TC_fine','fit_thr_BC_fine','-v7.3');
toc;
%% func
function metric = get_accuracy_index(TP,FP,TN,FN,beta)
    precision = TP/(TP+FP);
    recall = TP/(TP+FN);
    specificity = TN/(TN+FP);
    F1 = 2*precision*recall/(precision+recall);
    F2 = (1+beta^2)*precision*recall/(beta^2*precision+recall);
    balanced_accuracy = (recall+specificity)/2;
    youden = recall+specificity-1;
    metric = [TP,FP,TN,FN,precision,recall,specificity,F1,balanced_accuracy,youden,F2];
end

