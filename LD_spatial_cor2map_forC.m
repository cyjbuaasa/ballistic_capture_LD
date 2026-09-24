%% multibody orbital dynamics
clc;clear 
%% parameter
addpath('D:\matlab_program\cyj_program\CRTBP\TO_evolution');
basic_parameter;
rmpath('D:\matlab_program\cyj_program\CRTBP\TO_evolution');
% 无量纲参数
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
m_s = parameter.m.Sun;       %%太阳质量比
rho_sun = parameter.d.as;   %%太阳与地月质心距离
rEM_hill = parameter.d.rEM_hill;
rm_hill = parameter.d.rm_hill;
ome_sun = parameter.t.ome_sun;
%% Unit dimensionless and Lagrangian points
start_time = cputime;
T=2*pi*TU; % 
LL=get_LL(mu,mu_EMS,TU);
option_LL = optimoptions('fsolve','StepTolerance',1e-12,'FunctionTolerance',1e-12,...
    'OptimalityTolerance',1e-12,'Display','none');
%% add STT & tensor path
addpath('D:\matlab_program\cyj_program\CRTBP\STT');
addpath(genpath('D:\MATLAB2023a\R2023a\tensor_toolbox-v3.8'));
%% de2/dt equation
tic;
% H0 = -847.500; % Initial Hamiltonian in BCRFBP
J0 = [2.95:0.005:3.045,3.055:0.005:3.15];
numJ0 = length(J0);
eps=1e-3;
t0 = 0;t_end = 4.5;

% sphere setting
r2 = 0.25;
psi2 = [-88*rad,[-85*rad:5*rad:0]]; gama = -75*rad:5*rad:75*rad;
% z0 = 0;dz0 = -0e-1;
%% ini parameters
tic
filepat = 'spatial_new0526\';

alp = [0:1:360]*rad;
% ballistic angle
psi1 = [0:1:360]*rad;

num_all = length(alp)*length(psi1);
[alp_mesh,psi1_mesh] = meshgrid(alp,psi1);

options_LD = odeset('RelTol',1e-9,'AbsTol',1e-9,'Events',@only_hit_LD);
% v-planar
v_rot = zeros(size(alp_mesh));

dx_eps = 1e-4;
nonlinear_id = 3.5;
% 1--FTLE; 2--HOFTLE(.1-STT; .2-DSTT; .3-TDSTT); 
% 3--LD(.1-E2; .2-H2; .3-v2; .4-a2; .5-H2|E2|r2|v2|onlyv(v-norm) & a2)
exp_alp = 1.25;
m_STT = 2; % m order STT
dyn_str = 'EM'; % 'EM' 'diff_sbcm'
LD_str = 'H2';

BC_spa_probility = zeros(length(gama),length(psi2));
%% for J0
tic;
for j = 2:numJ0
    tic;
    fprintf("J0 = %f\n",J0(j));

    BC_mapnew = containers.Map('KeyType', 'char', 'ValueType', 'any');
     filepat1 = strcat(filepat,'J',num2str(J0(j)),'_r2',num2str(r2));
        mkdir(filepat1);  %%不存在时候，创建文件夹
    for k = 1:length(psi2)
    for n = 1:length(gama)
     fprintf('psi2 = %f, gama = %f \n',psi2(k)/rad,gama(n)/rad);
   
    % filepat2 = strcat('spatial\','J',num2str(J0(j)),'_r2',num2str(r2),...
    %     '_gama',num2str(gama(n)/rad),'_psi2',num2str(psi2(k)/rad));

    key = strcat('gama',num2str(gama(n)/rad),'_psi2',num2str(psi2(k)/rad));

    LD_FTLE_mesh = zeros(size(alp_mesh));
    LD_a2_mesh = zeros(size(alp_mesh));
    Np_prl_all = zeros(length(alp),length(psi1));
    % eve_struct_all = cell(length(alp),length(psi1));
    parfor i = 1:num_all
        r2_vec = r2*[cos(alp_mesh(i))*cos(gama(n));sin(alp_mesh(i))*cos(gama(n));sin(gama(n))];
        r1_vec = [1;0;0]+r2_vec;
        r1 = norm(r1_vec);
        
        vv = fsolve(@(v)(1*(v^2)-((1-mu+r2_vec(1))^2+r2_vec(2)^2) ...
            -2*(1-mu)./r1-2*mu./r2-mu*(1-mu)+J0(j)),1e-2,option_LL);
        Np_prl=0;% eve_struct=struct('tend',0,'xend',0,'eve',0);
        if vv >= 1e-3
            v_rot(i) = vv;
            % 转换 rot-v local
            T = [cos(gama(n)),0,sin(gama(n));0,1,0;-sin(gama(n)),0,cos(gama(n))]*...
            [cos(alp_mesh(i)),sin(alp_mesh(i)),0;-sin(alp_mesh(i)),cos(alp_mesh(i)),0;0,0,1];
            v_vlocal = v_rot(i)*[sin(psi2(k));cos(psi2(k))*cos(psi1_mesh(i));...
                cos(psi2(k))*sin(psi1_mesh(i))];
            v0 = T\v_vlocal;
            
            y0 = [[r2_vec+[1-mu;0;0]];v0];
            if nonlinear_id >1 && nonlinear_id < 3
                [orbDEX,Np_prl,aa,t,X,tend,xend,eve] = get_Dyindex(t0,t_end,y0,nonlinear_id,m_STT,...
                dyn_str,eps,exp_alp,dx_eps,LD_str,options_LD,mu);
            else
                [orbDEX,Np_prl,aa,tend,xend,eve] = get_Dyindex_1_3_mex(t0,t_end,y0,nonlinear_id,m_STT,...
                dyn_str,eps,exp_alp,dx_eps,LD_str,mu);
            end
            if nonlinear_id == 3.5
                LD_FTLE_mesh(i) = orbDEX(1);
                LD_a2_mesh(i) = orbDEX(2);
            else
                LD_FTLE_mesh(i) = orbDEX;
            end
        end
        Np_prl_all(i) = Np_prl;% eve_struct_all{i} = eve_struct;
        if mod(i,1e4) == 0, fprintf("%d done\n",i); end
    end
    toc;
    % 数据处理
    LD_FTLE_mesh(LD_FTLE_mesh == 0) = nan;% 有效 FTLE map
    LD_a2_mesh(LD_a2_mesh == 0) = nan; % 有效a2 LD map
    Np_prl_all(isnan(LD_a2_mesh)) = nan;
    Np_prl_all((LD_a2_mesh <= 1.0)) = nan;% 有效近月点map
     % 捕获可能性
    BC_spa_probility(n,k) = (numel(Np_prl_all(Np_prl_all>=2)))...
        /(numel(Np_prl_all)-numel(Np_prl_all(isnan(Np_prl_all))));
    
    LD_BC = LD_FTLE_mesh;LD_BC((LD_a2_mesh <= 1.0)) = nan;
    LD_BC(Np_prl_all<2)=nan;
    
    [row_TC, col_TC] = find(LD_a2_mesh > 1.0);
    ind_TC = sub2ind(size(LD_BC), row_TC, col_TC);
    TC_ang = [alp(col_TC);psi1(row_TC);]';
    LD_TC = LD_FTLE_mesh(ind_TC);

    [row_BC, col_BC] = find(Np_prl_all>=2);
    BC_ang = [alp(col_BC);psi1(row_BC);]';
    ind_BC = sub2ind(size(LD_BC), row_BC, col_BC);
    LD_BC2 = LD_BC(ind_BC);

    BC_mapnew(key) = struct('LD_TC',LD_TC,'TC_ang',TC_ang,'LD_BC',LD_BC2,'BC_ang',BC_ang,...
        'Np',Np_prl_all);
    %% LD
    % fig_LD = figure('Position', [46.60,192.2,949.60,560], 'Color', 'white');
    % hAx = axes;
    % LD_field=surf(hAx,alp_mesh/rad,psi1_mesh/rad...
    %     ,LD_FTLE_mesh, 'EdgeColor', 'none');
    % colormap(hAx,jet);c=colorbar(hAx);
    % xlim([min(alp/rad),max(alp/rad)]); ylim([min(psi1/rad),max(psi1/rad)]);
    % xlabel('\alpha \rm[deg]','fontsize',25,'fontname','Times New Roman','FontWeight','normal');
    % ylabel('\Psi_1 \rm[deg]','fontsize',25,'fontname','Times New Roman','FontWeight','normal');
    % set(gca,'FontName','Times New Roman','FontSize',25,'FontWeight','normal');
    % set(LD_field,'ButtonDownFcn',{@mycall_map_spatial,hAx,...
    %         t0,t0+t_end,alp_mesh,psi1_mesh,LD_FTLE_mesh,r2,gama(n),psi2(k),...
    %         J0(j),nonlinear_id,m_STT,dyn_str,eps,exp_alp,dx_eps,LD_str});
    % figname_LD = fullfile(filepat1,strcat('LD_',num2str(nonlinear_id),'_',LD_str,...
    %     '_expalp',num2str(exp_alp),'_',datestr(datetime('today'),'mmdd'),...
    %     '.fig'));
    %% if nonlinear_id == 3.5   
    % LD_FTLE_fixfroma2 = LD_FTLE_mesh;
    % LD_FTLE_fixfroma2(LD_a2_mesh <= 0.917) = nan;
    % figLDcap = figure('Position', [46.60,192.2,949.60,560], 'Color', 'white');
    % surf(alp_mesh/rad,psi1_mesh/rad,LD_FTLE_fixfroma2, 'EdgeColor', 'none');
    % colormap(jet);colorbar();
    % xlim([min(alp/rad),max(alp/rad)]); ylim([min(psi1/rad),max(psi1/rad)]);
    % xlabel('\alpha \rm[deg]','fontsize',25,'fontname','Times New Roman','FontWeight','normal');
    % ylabel('\Psi_1 \rm[deg]','fontsize',25,'fontname','Times New Roman','FontWeight','normal');
    % set(gca,'FontName','Times New Roman','FontSize',25,'FontWeight','normal');
    % figname_LD_cap = fullfile(filepat1,strcat('LDcap_',num2str(nonlinear_id),'_',LD_str,...
    %     '_expalp',num2str(exp_alp),'_',datestr(datetime('today'),'mmdd'),...
    %     '.fig'));
    % 
    % LD_FTLE_fixfroma3 = LD_FTLE_fixfroma2;
    % LD_FTLE_fixfroma3(LD_FTLE_fixfroma2 <= 4.05) = nan;
    % figLDBC = figure('Position', [46.60,192.2,949.60,560], 'Color', 'white');
    % surf(alp_mesh/rad,psi1_mesh/rad,LD_FTLE_fixfroma3, 'EdgeColor', 'none');
    % colormap(jet);colorbar();
    % if ~all(isnan(LD_FTLE_mesh(:)))
    %     clim([min(LD_FTLE_mesh(:)),max(LD_FTLE_mesh(:))]);
    % end
    % xlim([min(alp/rad),max(alp/rad)]); ylim([min(psi1/rad),max(psi1/rad)]);
    % xlabel('\alpha \rm[deg]','fontsize',25,'fontname','Times New Roman','FontWeight','normal');
    % ylabel('\Psi_1 \rm[deg]','fontsize',25,'fontname','Times New Roman','FontWeight','normal');
    % set(gca,'FontName','Times New Roman','FontSize',25,'FontWeight','normal');
    % figname_LD_BC = fullfile(filepat1,strcat('BChighLD_',num2str(nonlinear_id),'_',LD_str,...
    %         '_expalp',num2str(exp_alp),'_',datestr(datetime('today'),'mmdd'),...
    %         '.fig'));
%% Np map
    % fig_Np = figure('Position', [46.60,192.2,949.60,560], 'Color', 'white');
    % % Np_prl_all(LD_a2_mesh <= 1.0) = nan;
    % surf(alp_mesh/rad,psi1_mesh/rad,Np_prl_all, 'EdgeColor', 'none');
    % colormap(jet);colorbar();
    % if ~all(isnan(Np_prl_all(:)))
    %     clim([0,max(Np_prl_all(:))]);
    % end
    % xlim([min(alp/rad),max(alp/rad)]); ylim([min(psi1/rad),max(psi1/rad)]);
    % xlabel('\alpha \rm[deg]','fontsize',25,'fontname','Times New Roman','FontWeight','normal');
    % ylabel('\Psi_1 \rm[deg]','fontsize',25,'fontname','Times New Roman','FontWeight','normal');
    % set(gca,'FontName','Times New Roman','FontSize',25,'FontWeight','normal');
    % figname_Np = fullfile(filepat1,strcat('Np_',num2str(nonlinear_id),'_',LD_str,...
    %     '_expalp',num2str(exp_alp),'_',datestr(datetime('today'),'mmdd'),...
    %     '.fig'));

    % save(strcat(filepat1,'\','LD_',num2str(nonlinear_id),'_',LD_str,'_expalp',num2str(exp_alp),...
    %         '_',datestr(datetime('today'),'mmdd'),'.mat'),'LD_FTLE_mesh','LD_field',...
    %         'LD_a2_mesh','hAx','psi1_mesh','alp_mesh','gama','psi2','BC_probility',...
    %         'alp','psi1','t0','t_end','nonlinear_id','m_STT','J0','Np_prl_all','eve_struct_all',...
    %          'r2','dyn_str','eps','exp_alp','dx_eps','LD_str');
    
    % savefig(fig_LD, figname_LD);
    % savefig(figLDcap, figname_LD_cap); savefig(figLDBC, figname_LD_BC);
    % savefig(fig_Np, figname_Np); 
    close all;
    end
    end
save_name_base = strcat('LD_BCmapnew_',num2str(nonlinear_id),'_',LD_str,'_expalp',num2str(exp_alp),...
    '_',datestr(datetime('today'),'mmdd'));
save_name = strcat(save_name_base,'.mat');
save_path = fullfile(filepat1,save_name);
if exist(save_path,'file') == 2
    save_name = strcat(save_name_base,'.mat');
    save_path = fullfile(filepat1,save_name);
end
save(save_path,'BC_mapnew','BC_spa_probility',...
    't0','t_end','nonlinear_id','J0','alp_mesh','psi1_mesh','-v7.3');
toc;

end
toc;
%% plot 正确做法不再需要太多数据
n=18;k=6;
key = strcat('gama',num2str(gama(n)/rad),'_psi2',num2str(psi2(k)/rad));
tmp=BC_mapnew(key);

% BC LD
fig_LD = figure('Position', [46.60,192.2,949.60,560], 'Color', 'white');
 hAx = axes;
 LD_field=scatter3(hAx,tmp.BC_ang(:,1)/rad,tmp.BC_ang(:,2)/rad,tmp.LD_BC(:),36,tmp.LD_BC(:),'filled');
 colormap(hAx,jet);c=colorbar(hAx);
 xlim([min(alp/rad),max(alp/rad)]); ylim([min(psi1/rad),max(psi1/rad)]);
 xlabel('\alpha \rm[deg]','fontsize',25,'fontname','Times New Roman','FontWeight','normal');
 ylabel('\Psi_1 \rm[deg]','fontsize',25,'fontname','Times New Roman','FontWeight','normal');
 set(gca,'FontName','Times New Roman','FontSize',25,'FontWeight','normal');
 set(LD_field,'ButtonDownFcn',{@mycall_map_spatiallessdata,hAx,...
         t0,t0+t_end,tmp.BC_ang(:,1),tmp.BC_ang(:,2),tmp.LD_BC,r2,gama(n),psi2(k),...
         J0(j),nonlinear_id,m_STT,dyn_str,eps,exp_alp,dx_eps,LD_str});

 % TC LD
fig_LDtc = figure('Position', [46.60,192.2,949.60,560], 'Color', 'white');
 hAx = axes;
 LD_field=scatter3(hAx,tmp.TC_ang(:,1)/rad,tmp.TC_ang(:,2)/rad,tmp.LD_TC(:),36,tmp.LD_TC(:),'filled');
 colormap(hAx,jet);c=colorbar(hAx);
 xlim([min(alp/rad),max(alp/rad)]); ylim([min(psi1/rad),max(psi1/rad)]);
 xlabel('\alpha \rm[deg]','fontsize',25,'fontname','Times New Roman','FontWeight','normal');
 ylabel('\Psi_1 \rm[deg]','fontsize',25,'fontname','Times New Roman','FontWeight','normal');
 set(gca,'FontName','Times New Roman','FontSize',25,'FontWeight','normal');
 set(LD_field,'ButtonDownFcn',{@mycall_map_spatiallessdata,hAx,...
         t0,t0+t_end,tmp.BC_ang(:,1),tmp.BC_ang(:,2),tmp.LD_TC,r2,gama(n),psi2(k),...
         J0(j),nonlinear_id,m_STT,dyn_str,eps,exp_alp,dx_eps,LD_str});

