%% Extract LCS-limited LD maps and BC maps for selected J0 values
% LCS-limited region: LD_a2_mesh > fitted LD_a2 threshold.
% BC region: LCS-limited region and Np_prl_all >= 2.

clear; clc;

%% parameter
rad = pi/180;

J0_plot = 3.05;%
r2 = 0.25;
z0 = 0;
dz0 = -0e-1;

filepat = 'accu_planar_0709_xmoon\';
threshold_file = 'accu_planar_0707\TC_BC_TP_TN\optmectric\opt_threshold.mat';
threshold_kind = 'BC'; % 'BC' or 'TC'

nonlinear_id = 3.5;
exp_alp = 1.25;
LD_str = 'H2';

LD_clim = [0 5];
jpeg_resolution = 600;
boundary_line_width = 2.0;
LCS_boundary_color = [0.85 0 0];
BC_boundary_color = [0 0 0];
close_fig_after_save = true;
date_tag = char(datetime('now','Format','MMdd'));

%% load threshold fitting
S_thr = load(threshold_file);
[threshold_fun,threshold_desc] = get_LDa2_threshold_fun(S_thr,threshold_kind);
%% tend
t_end_cons = 5;
%% loop over J0 and t_end
for iJ = 1:length(J0_plot)
    J0_now = J0_plot(iJ);
    filepat1 = strcat(filepat,'J',num2str(J0_now),'_r2',num2str(r2),'_z0',num2str(z0),...
        '_dz0',num2str(dz0));
    matstr = get_LD_mat_path(filepat1,nonlinear_id,LD_str,exp_alp);
    LD_a2threshold = threshold_fun(J0_now);

    fprintf('\nData file: %s\n',matstr);
    fprintf('J0 = %.6f, LD_a2 threshold = %.8f (%s)\n',...
        J0_now,LD_a2threshold,threshold_desc);

    M = matfile(matstr);
    t_end = M.t_end;
    alp = M.alp;
    Psi = M.Psi;
    alp_mesh = M.alp_mesh;
    Psi_mesh = M.Psi_mesh;

    num_tend = length(t_end);
    for kt = 1:num_tend
        t_now = t_end(kt);
        if abs(t_now-t_end_cons) > 1e-3, continue; end
        filepat_t = fullfile(filepat1,strcat('T',num2str(t_now)));
        if ~exist(filepat_t,'dir')
            mkdir(filepat_t);
        end

        LD_FTLE_mesh = M.LD_FTLE_all(:,:,kt);
        LD_a2_mesh = M.LD_a2_all(:,:,kt);
        Np_prl_all_real = M.Np_prl_all_t(:,:,end);

        LCS_mask = LD_a2_mesh > LD_a2threshold & ~isnan(LD_FTLE_mesh);
        BC_mask = LCS_mask & Np_prl_all_real >= 2;

        LD_LCS = LD_FTLE_mesh;
        LD_LCS(~LCS_mask) = nan;

        LD_BC = LD_FTLE_mesh;
        LD_BC(~BC_mask) = nan;

        LCS_boundaries = bwboundaries(LCS_mask,'noholes');
        BC_boundaries = bwboundaries(BC_mask,'noholes');

        basename_LCS = fullfile(filepat_t,strcat('LD_LCS_',num2str(nonlinear_id),'_',LD_str,...
            '_T',num2str(t_now),'_expalp',num2str(exp_alp),'_',date_tag));
        basename_BC = fullfile(filepat_t,strcat('LD_BC_',num2str(nonlinear_id),'_',LD_str,...
            '_T',num2str(t_now),'_expalp',num2str(exp_alp),'_',date_tag));

        fig_LCS = plot_LD_mask_map(alp_mesh,Psi_mesh,alp,Psi,rad,LD_LCS,LD_clim,...
            LCS_boundaries,LCS_boundary_color,boundary_line_width,...
            strcat('LD LCS, J0 = ',num2str(J0_now),', T = ',num2str(t_now)));
        save_LD_figure(fig_LCS,basename_LCS,jpeg_resolution);

        fig_BC = plot_LD_mask_map(alp_mesh,Psi_mesh,alp,Psi,rad,LD_BC,LD_clim,...
            BC_boundaries,BC_boundary_color,boundary_line_width,...
            strcat('LD BC, J0 = ',num2str(J0_now),', T = ',num2str(t_now)));
        save_LD_figure(fig_BC,basename_BC,jpeg_resolution);

        if close_fig_after_save
            close(fig_LCS);
            close(fig_BC);
        end

        fprintf('T = %.6g done: LCS points = %d, BC points = %d\n',...
            t_now,nnz(LCS_mask),nnz(BC_mask));
    end
end

%% local functions
function matstr = get_LD_mat_path(filepat1,nonlinear_id,LD_str,exp_alp)
    file_pattern = strcat('LD_',num2str(nonlinear_id),'_',LD_str,...
        '_expalp',num2str(exp_alp),'_changeT_*.mat');
    mat_files = dir(fullfile(filepat1,file_pattern));
    if isempty(mat_files)
        error('Cannot find LD mat file under %s with pattern %s.',...
            filepat1,file_pattern);
    end

    [~,id_latest] = max([mat_files.datenum]);
    matstr = fullfile(mat_files(id_latest).folder,mat_files(id_latest).name);
end

function [threshold_fun,threshold_desc] = get_LDa2_threshold_fun(S_thr,threshold_kind)
    threshold_kind = upper(threshold_kind);

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
        threshold_desc = 'polyfit p';
        return;
    end

    if isfield(S_thr,fit_name)
        fit_thr = S_thr.(fit_name);
        p = fit_thr(1,:);
        threshold_fun = @(J0_now) polyval(p,J0_now);
        threshold_desc = strcat(threshold_kind,' polyfit');
        return;
    end

    if isfield(S_thr,value_name) && isfield(S_thr,'J0')
        fit_value = S_thr.(value_name);
        J0_thr = S_thr.J0(:);
        value_thr = fit_value(:,1);
        id_good = ~isnan(J0_thr) & ~isnan(value_thr);
        threshold_fun = @(J0_now) interp1(J0_thr(id_good),value_thr(id_good),...
            J0_now,'linear','extrap');
        threshold_desc = strcat(threshold_kind,' fitted value interpolation');
        return;
    end

    error(['Cannot find p, %s, or %s in opt threshold file. ',...
        'Please check threshold_kind or opt_threshold.mat variables.'],...
        fit_name,value_name);
end

function fig = plot_LD_mask_map(alp_mesh,Psi_mesh,alp,Psi,rad,LD_map,LD_clim,...
    boundaries,boundary_color,boundary_line_width,fig_title)
    fig = figure('Position',[46.60,192.2,949.60,560],'Color','white');
    hAx = axes(fig);
    surf(hAx,alp_mesh/rad,Psi_mesh/rad,LD_map,'EdgeColor','none');
    colormap(hAx,jet);
    colorbar(hAx);
    clim(hAx,LD_clim);
    hold(hAx,'on');

    % for k = 1:length(boundaries)
    %     b = boundaries{k};
    %     plot3(hAx,alp(b(:,2))/rad,Psi(b(:,1))/rad,...
    %         LD_clim(2)*ones(size(b(:,1))),...
    %         'Color',boundary_color,'LineWidth',boundary_line_width);
    % end

    xlim(hAx,[min(alp/rad),max(alp/rad)]);
    ylim(hAx,[min(Psi/rad),max(Psi/rad)]);
    xlabel(hAx,'\alpha \rm[deg]','fontsize',25,'fontname','Times New Roman','FontWeight','normal');
    ylabel(hAx,'\Psi \rm[deg]','fontsize',25,'fontname','Times New Roman','FontWeight','normal');
    % title(hAx,fig_title,'FontName','Times New Roman','FontSize',18,'FontWeight','normal');
    set(hAx,'FontName','Times New Roman','FontSize',25,'FontWeight','normal');
    view(hAx,2);
end

function save_LD_figure(fig,basename,jpeg_resolution)
    savefig(fig,strcat(basename,'.fig'));
    exportgraphics(fig,strcat(basename,'.jpg'),'Resolution',jpeg_resolution);
end