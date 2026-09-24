%% Stack LD fields over J0 on the alpha-Psi-J0 section
% x = alpha [deg], y = Psi [deg], z = J0.
% For each J0, keep LD_FTLE values only where LD_a2 exceeds the fitted
% threshold, and draw the corresponding LD-a2 boundary.

clear; clc;

%% parameter
rad = pi/180;

J0 = 2.95:0.01:3.20;
numJ0 = length(J0);

t_plot = 5;
r2 = 0.25;
z0 = 0;
dz0 = -0e-1;

filepat = 'accu_planar_0709_xmoon\';
threshold_file = 'accu_planar_0707\TC_BC_TP_TN\optmectric\opt_threshold.mat';
threshold_kind = 'BC'; % 'BC' or 'TC'

nonlinear_id = 3.5;
exp_alp = 1.25;
LD_str = 'H2';

plot_stride = 2; % set to 1 for full resolution scatter points
marker_size = 6;
boundary_line_width = 2.0;
boundary_color = [0.85 0 0];
save_fig = true;
save_dir = fullfile(filepat,'J0_stack_plot');

%% load threshold fitting
S_thr = load(threshold_file);
[threshold_fun,threshold_desc] = get_LDa2_threshold_fun(S_thr,threshold_kind,t_plot);

if save_fig && ~exist(save_dir,'dir')
    mkdir(save_dir);
end

%% plot all J0 planes
fig_stack = figure('Name',strcat('LD stack, t = ',num2str(t_plot)),...
    'Position',[44,74,1120,720],'Color','w');
hAx = axes(fig_stack);
J0_loaded = [];
LD_color_min = inf;
LD_color_max = -inf;

for j = 25:numJ0
    J0_now = J0(j);
    filepat1 = strcat(filepat,'J',num2str(J0_now),'_r2',num2str(r2),'_z0',num2str(z0),...
        '_dz0',num2str(dz0));
    if J0_now < 3.19
        matstr = strcat(filepat1,'\','LD_',num2str(nonlinear_id),'_',LD_str,...
            '_expalp',num2str(exp_alp),'_changeT_','0709','.mat');
    else
        matstr = strcat(filepat1,'\','LD_',num2str(nonlinear_id),'_',LD_str,...
            '_expalp',num2str(exp_alp),'_changeT_','0710','.mat');
    end

    if ~isfile(matstr)
        warning('Missing file, skip J0 = %.4f: %s',J0_now,matstr);
        continue;
    end

    M = matfile(matstr);
    t_end = M.t_end;
    kt = find(abs(t_end(:)'-t_plot) < 1e-12,1);
    if isempty(kt)
        warning('t_plot = %.6g is not in t_end for J0 = %.4f, skip. Available: %s',...
            t_plot,J0_now,num2str(t_end(:)'));
        continue;
    end

    alp_mesh = M.alp_mesh;
    Psi_mesh = M.Psi_mesh;
    alp = M.alp;
    Psi = M.Psi;
    LD_FTLE_mesh = M.LD_FTLE_all(:,:,kt);
    LD_a2_mesh = M.LD_a2_all(:,:,kt);

    LD_a2threshold = threshold_fun(J0_now);
    LD_mask = LD_a2_mesh > LD_a2threshold & ~isnan(LD_FTLE_mesh);

    if ~any(LD_mask(:))
        warning('No retained LD points for J0 = %.4f, threshold = %.6g.',...
            J0_now,LD_a2threshold);
        continue;
    end

    LD_color_min = min(LD_color_min,min(LD_FTLE_mesh(LD_mask),[],'omitnan'));
    LD_color_max = max(LD_color_max,max(LD_FTLE_mesh(LD_mask),[],'omitnan'));

    id_plot = false(size(LD_mask));
    id_plot(1:plot_stride:end,1:plot_stride:end) = true;
    id_plot = id_plot & LD_mask;

    scatter3(hAx,alp_mesh(id_plot)/rad,Psi_mesh(id_plot)/rad,...
        J0_now*ones(numel(find(id_plot)),1),marker_size,...
        LD_FTLE_mesh(id_plot),'filled');
    hold(hAx,'on');
    colormap(hAx,jet);

    boundaries = bwboundaries(LD_mask,'noholes');
    for k = 1:length(boundaries)
        b = boundaries{k};
        plot3(hAx,alp(b(:,2))/rad,Psi(b(:,1))/rad,...
            J0_now*ones(size(b(:,1))),...
            'Color',boundary_color,'LineWidth',boundary_line_width);
    end

    J0_loaded(end+1) = J0_now; %#ok<SAGROW>
    fprintf('J0 = %.4f done, LD_a2 threshold = %.6f, retained = %d\n',...
        J0_now,LD_a2threshold,nnz(LD_mask));
end

if isempty(J0_loaded)
    error('No J0 data were loaded. Please check filepat and data files.');
end

cb = colorbar(hAx);
cb.Label.String = 'LD';
cb.Label.FontName = 'Times New Roman';
cb.Label.FontSize = 20;
if isfinite(LD_color_min) && isfinite(LD_color_max) && LD_color_min < LD_color_max
    clim(hAx,[LD_color_min LD_color_max]);
end

xlabel(hAx,'\alpha \rm[deg]','fontsize',25,'fontname','Times New Roman','FontWeight','normal');
ylabel(hAx,'\Psi \rm[deg]','fontsize',25,'fontname','Times New Roman','FontWeight','normal');
zlabel(hAx,'J_0','fontsize',25,'fontname','Times New Roman','FontWeight','normal');
% title(hAx,strcat('LD field stack at t = ',num2str(t_plot),...
%     ', threshold: ',threshold_desc),'FontName','Times New Roman','FontSize',18);
set(hAx,'FontName','Times New Roman','FontSize',24,'FontWeight','normal');
grid(hAx,'on'); box(hAx,'on');
xlim(hAx,[min(alp/rad),max(alp/rad)]);
ylim(hAx,[min(Psi/rad),max(Psi/rad)]);
zlim(hAx,[min(J0_loaded),max(J0_loaded)]);
view(hAx,[-38,24]);

if save_fig
    figname = fullfile(save_dir,strcat('LD_J0_stack_',threshold_kind,...
        '_T',num2str(t_plot),'_',datestr(datetime('today'),'mmdd'),'.fig'));
    pngname = fullfile(save_dir,strcat('LD_J0_stack_',threshold_kind,...
        '_T',num2str(t_plot),'_',datestr(datetime('today'),'mmdd'),'.png'));
    savefig(fig_stack,figname);
    exportgraphics(fig_stack,pngname,'Resolution',300);
end

%% local functions
function [threshold_fun,threshold_desc] = get_LDa2_threshold_fun(S_thr,threshold_kind,t_plot)
    threshold_kind = upper(threshold_kind);

    t_end_thr = [];
    if isfield(S_thr,'t_end')
        t_end_thr = S_thr.t_end(:)';
    end

    switch threshold_kind
        case 'BC'
            fit_name = 'fit_thr_BC';
            value_name = 'fit_thr_BC_value';
        case 'TC'
            fit_name = 'fit_thr_TC';
            value_name = 'fit_thr_TC_value';
        otherwise
            error('threshold_kind must be ''BC'' or ''TC''.');
    end

    if isfield(S_thr,'p') && isnumeric(S_thr.p)
        p = S_thr.p(:)';
        threshold_fun = @(J0_now) polyval(p,J0_now);
        threshold_desc = sprintf('polyfit p, t = %.6g',t_plot);
        return;
    end

    if isfield(S_thr,fit_name)
        fit_thr = S_thr.(fit_name);
        kt_thr = get_threshold_time_id(t_end_thr,t_plot,size(fit_thr,1));
        p = fit_thr(kt_thr,:);
        threshold_fun = @(J0_now) polyval(p,J0_now);
        threshold_desc = sprintf('%s polyfit, t = %.6g',threshold_kind,t_plot);
        return;
    end

    if isfield(S_thr,value_name) && isfield(S_thr,'J0')
        fit_value = S_thr.(value_name);
        kt_thr = get_threshold_time_id(t_end_thr,t_plot,size(fit_value,2));
        J0_thr = S_thr.J0(:);
        value_thr = fit_value(:,kt_thr);
        id_good = ~isnan(J0_thr) & ~isnan(value_thr);
        threshold_fun = @(J0_now) interp1(J0_thr(id_good),value_thr(id_good),...
            J0_now,'linear','extrap');
        threshold_desc = sprintf('%s fitted value interpolation, t = %.6g',...
            threshold_kind,t_plot);
        return;
    end

    error(['Cannot find %s or %s in opt threshold file. ',...
        'Please check threshold_kind or opt_threshold.mat variables.'],...
        fit_name,value_name);
end

function kt = get_threshold_time_id(t_end_thr,t_plot,num_time)
    if ~isempty(t_end_thr)
        kt = find(abs(t_end_thr-t_plot) < 1e-12,1);
        if isempty(kt)
            error('t_plot = %.6g is not in threshold t_end. Available values are: %s',...
                t_plot,num2str(t_end_thr));
        end
    else
        kt = 1;
    end

    if kt > num_time
        error('Threshold time index kt = %d exceeds threshold time dimension = %d.',...
            kt,num_time);
    end
end
