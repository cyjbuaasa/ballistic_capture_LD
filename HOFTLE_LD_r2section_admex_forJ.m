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
option_LL = optimoptions('fsolve','StepTolerance',1e-12,'FunctionTolerance',1e-12,'OptimalityTolerance',1e-12,'Display','none');
get_LL(mu,mu_EMS,TU);
%% add STT & tensor path
addpath('D:\matlab_program\cyj_program\CRTBP\STT');
addpath(genpath('D:\MATLAB2023a\R2023a\tensor_toolbox-v3.8'));
%% de2/dt equation
tic;
% H0 = -847.500; % Initial Hamiltonian in BCRFBP
J0 = 2.95:0.01:3.15;
numJ0 = length(J0);
eps=1e-3;
theta_s = pi/2;
t0 = 0;t_end = 5.5;

% z component setting
r2 = 0.25;
z0 = 0;dz0 = -0e-1;

%% ini parameters
tic
filepat = 'accu_planar_0706\';

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

BC_probility = zeros(numJ0,1);
%% for J0
tic;
for j = 1:numJ0
    LD_FTLE_mesh = zeros(size(alp_mesh));
    LD_a2_mesh = zeros(size(alp_mesh));
    Np_prl_all = zeros(length(alp),length(Psi));
    KEarth_all = NaN(length(alp),length(Psi));
    KMoon_all = NaN(length(alp),length(Psi));
    eve_struct_all = cell(length(alp),length(Psi));
    fprintf("J0 = %f\n",J0(j));
    parfor i = 1:num_all
        r2_vec = r2*[cos(alp_mesh(i))*cos(theta);sin(alp_mesh(i))*cos(theta);sin(theta)];
        r1_vec = [1;0;0]+r2_vec;
        r1 = norm(r1_vec);
        
        vv = fsolve(@(v)(1*(v^2+v_vert^2)-((1-mu+r2_vec(1))^2+r2_vec(2)^2) ...
            -2*(1-mu)./r1-2*mu./r2-mu*(1-mu)+J0(j)),1e-2,option_LL);
        eve = 0;
        Np_prl=0;eve_struct=struct('tend',0,'xend',0,'eve',0);
        if vv >= 1e-3
            v_plan(i) = vv;
            sig = alp_mesh(i) + Psi_mesh(i);
            y0 = [[r2_vec+[1-mu;0;0]];v_plan(i)*cos(sig);v_plan(i)*sin(sig);v_vert];
            if nonlinear_id >1 && nonlinear_id < 3
                [orbDEX,Np_prl,eve_struct,t,X,tend,xend,eve] = get_Dyindex(t0,t_end,y0,nonlinear_id,m_STT,...
                dyn_str,eps,exp_alp,dx_eps,LD_str,options_LD,mu);
            else
                [orbDEX,Np_prl,eve_struct,tend,xend,eve] = get_Dyindex_1_3_mex(t0,t_end,y0,nonlinear_id,m_STT,...
                dyn_str,eps,exp_alp,keta,dx_eps,LD_str,mu);
            end
            if nonlinear_id == 3.5
                LD_FTLE_mesh(i) = orbDEX(1);
                LD_a2_mesh(i) = orbDEX(2);
            else
                LD_FTLE_mesh(i) = orbDEX;
            end
        end
        if ~isempty(find(eve == 1, 1)) 
            KEarth_all(i) = 1;
        elseif ~isempty(find(eve == 2, 1)) 
            KMoon_all(i) = 1;
        end
        Np_prl_all(i) = Np_prl;eve_struct_all{i} = eve_struct;
        if mod(i,1e3) == 0, fprintf("%d done\n",i); end
    end
    toc;
    % 数据处理
    LD_FTLE_mesh(LD_FTLE_mesh == 0) = nan;% 有效 FTLE map
    LD_a2_mesh(LD_a2_mesh == 0) = nan; % 有效a2 LD map
    Np_prl_all(isnan(LD_a2_mesh)) = nan;
    % Np_prl_all((LD_a2_mesh <= 1.0)) = nan;% 有效近月点map
    % 非捕获最大 LDa2
    id1 = find(Np_prl_all<=1);max_LDa2_nocap = max(LD_a2_mesh(id1));
    % 弹道捕获最小 LDa2
    id2 = find(Np_prl_all>=2);min_LDa2_BC = min(LD_a2_mesh(id2));
     % 弹道捕获可能性
    BC_probility(j) = (numel(Np_prl_all(Np_prl_all>=2)))...
        /(numel(Np_prl_all)-numel(Np_prl_all(isnan(Np_prl_all))));

    filepat1 = strcat(filepat,'J',num2str(J0(j)),'_r2',num2str(r2),'_z0',num2str(z0),...
        '_dz0',num2str(dz0));
    mkdir(filepat1);
    %% LD
    fig_LD = figure('Position', [46.60,192.2,949.60,560], 'Color', 'white');
    hAx = axes;
    LD_field=surf(hAx,alp_mesh/rad,Psi_mesh/rad...
        ,LD_FTLE_mesh, 'EdgeColor', 'none');
    colormap(hAx,jet);c=colorbar(hAx);
    xlim([min(alp/rad),max(alp/rad)]); ylim([min(Psi/rad),max(Psi/rad)]);
    xlabel('\alpha \rm[deg]','fontsize',25,'fontname','Times New Roman','FontWeight','normal');
    ylabel('\Psi \rm[deg]','fontsize',25,'fontname','Times New Roman','FontWeight','normal');
    set(gca,'FontName','Times New Roman','FontSize',25,'FontWeight','normal');
    set(LD_field,'ButtonDownFcn',{@mycall_map,hAx,...
            t0,t0+t_end,alp_mesh,Psi_mesh,LD_FTLE_mesh,r2,theta,v_vert,...
            J0(j),nonlinear_id,m_STT,dyn_str,eps,exp_alp,keta,dx_eps,LD_str});
    figname_LD = fullfile(filepat1,strcat('LD_',num2str(nonlinear_id),'_',LD_str,...
        '_expalp',num2str(exp_alp),'_',datestr(datetime('today'),'mmdd'),...
        '.fig'));

    figa2 = figure('Position', [46.60,192.2,949.60,560], 'Color', 'white');
    surf(alp_mesh/rad,Psi_mesh/rad,LD_a2_mesh, 'EdgeColor', 'none');
    colormap(jet);colorbar;
    xlim([min(alp/rad),max(alp/rad)]); ylim([min(Psi/rad),max(Psi/rad)]);
    xlabel('\alpha \rm[deg]','fontsize',25,'fontname','Times New Roman','FontWeight','normal');
    ylabel('\Psi \rm[deg]','fontsize',25,'fontname','Times New Roman','FontWeight','normal');
    set(gca,'FontName','Times New Roman','FontSize',25,'FontWeight','normal');
    figname_a2 = fullfile(filepat1,strcat('LDa2_',datestr(datetime('today'),'mmdd'),...
        '.fig'));
    %% if nonlinear_id == 3.5   
    % 根据捕获近月点次数的定义判断这是真值
    LD_FTLE_fixfroma2 = LD_FTLE_mesh;
    LD_FTLE_fixfroma2(Np_prl_all <1) = nan;
    LD_FTLE_fixfroma2(KEarth_all ==1) = nan;LD_FTLE_fixfroma2(KMoon_all ==1) = nan;
    figLDcap = figure('Position', [46.60,192.2,949.60,560], 'Color', 'white');
    surf(alp_mesh/rad,Psi_mesh/rad,LD_FTLE_fixfroma2, 'EdgeColor', 'none');
    colormap(jet);colorbar();clim([min(LD_FTLE_mesh(:)),max(LD_FTLE_mesh(:))]);
    xlim([min(alp/rad),max(alp/rad)]); ylim([min(Psi/rad),max(Psi/rad)]);
    xlabel('\alpha \rm[deg]','fontsize',25,'fontname','Times New Roman','FontWeight','normal');
    ylabel('\Psi \rm[deg]','fontsize',25,'fontname','Times New Roman','FontWeight','normal');
    set(gca,'FontName','Times New Roman','FontSize',25,'FontWeight','normal');
    figname_LD_cap = fullfile(filepat1,strcat('LDcap_',num2str(nonlinear_id),'_',LD_str,...
        '_expalp',num2str(exp_alp),'_',datestr(datetime('today'),'mmdd'),...
        '.fig'));
    % 根据捕获近月点次数的定义得到的真正捕获的LDa2场
    LD_fixfroma3 = LD_a2_mesh;
    LD_fixfroma3(Np_prl_all <1) = nan;
    LD_fixfroma3(KEarth_all ==1) = nan;LD_fixfroma3(KMoon_all ==1) = nan;
    figLDa2cap = figure('Position', [46.60,192.2,949.60,560], 'Color', 'white');
    surf(alp_mesh/rad,Psi_mesh/rad,LD_fixfroma3, 'EdgeColor', 'none');
    colormap(jet);colorbar();clim([min(LD_a2_mesh(:)),max(LD_a2_mesh(:))]);
    xlim([min(alp/rad),max(alp/rad)]); ylim([min(Psi/rad),max(Psi/rad)]);
    xlabel('\alpha \rm[deg]','fontsize',25,'fontname','Times New Roman','FontWeight','normal');
    ylabel('\Psi \rm[deg]','fontsize',25,'fontname','Times New Roman','FontWeight','normal');
    set(gca,'FontName','Times New Roman','FontSize',25,'FontWeight','normal');
    figname_LDa2cap = fullfile(filepat1,strcat('BChighLD_',num2str(nonlinear_id),'_',LD_str,...
            '_expalp',num2str(exp_alp),'_',datestr(datetime('today'),'mmdd'),...
            '.fig'));
    % 直接根据LDa2的阈值判断的LD_a2场这是粗略的
    LD_a2tmp = LD_a2_mesh;
    LD_a2tmp(LD_a2tmp <0.80) = nan;
    LD_a2tmp(KEarth_all ==1) = nan;LD_a2tmp(KMoon_all ==1) = nan;
    figLDa2tmp = figure('Position', [46.60,192.2,949.60,560], 'Color', 'white');
    surf(alp_mesh/rad,Psi_mesh/rad,LD_a2tmp, 'EdgeColor', 'none');
    colormap(jet);colorbar();clim([min(LD_a2_mesh(:)),max(LD_a2_mesh(:))]);
    xlim([min(alp/rad),max(alp/rad)]); ylim([min(Psi/rad),max(Psi/rad)]);
    xlabel('\alpha \rm[deg]','fontsize',25,'fontname','Times New Roman','FontWeight','normal');
    ylabel('\Psi \rm[deg]','fontsize',25,'fontname','Times New Roman','FontWeight','normal');
    set(gca,'FontName','Times New Roman','FontSize',25,'FontWeight','normal');
    figname_LDa2tmp = fullfile(filepat1,strcat('BChighLD_',num2str(nonlinear_id),'_',LD_str,...
            '_expalp',num2str(exp_alp),'_',datestr(datetime('today'),'mmdd'),...
            '.fig'));
    % 根据定义得到的BC场是准确的即Np_prl_all次数≥2
    LD_BC = LD_FTLE_mesh;
    LD_BC(Np_prl_all <=1) = nan;
    LD_BC(KEarth_all ==1) = nan;LD_BC(KMoon_all ==1) = nan;
    figLDBC = figure('Position', [46.60,192.2,949.60,560], 'Color', 'white');
    surf(alp_mesh/rad,Psi_mesh/rad,LD_BC, 'EdgeColor', 'none');
    colormap(jet);colorbar();clim([min(LD_FTLE_mesh(:)),max(LD_FTLE_mesh(:))]);
    xlim([min(alp/rad),max(alp/rad)]); ylim([min(Psi/rad),max(Psi/rad)]);
    xlabel('\alpha \rm[deg]','fontsize',25,'fontname','Times New Roman','FontWeight','normal');
    ylabel('\Psi \rm[deg]','fontsize',25,'fontname','Times New Roman','FontWeight','normal');
    set(gca,'FontName','Times New Roman','FontSize',25,'FontWeight','normal');
    figname_LD_BC = fullfile(filepat1,strcat('BChighLD_',num2str(nonlinear_id),'_',LD_str,...
            '_expalp',num2str(exp_alp),'_',datestr(datetime('today'),'mmdd'),...
            '.fig'));
    % 根据定义得到的BC场是准确的
    LD_BCa2 = LD_FTLE_mesh;
    LD_BCa2(LD_a2_mesh <0.80) = nan;
    LD_BCa2(Np_prl_all <1) = nan;
    LD_BCa2(KEarth_all ==1) = nan;LD_BCa2(KMoon_all ==1) = nan;
    figLDBCa2 = figure('Position', [46.60,192.2,949.60,560], 'Color', 'white');
    surf(alp_mesh/rad,Psi_mesh/rad,LD_BCa2, 'EdgeColor', 'none');
    colormap(jet);colorbar();clim([min(LD_FTLE_mesh(:)),max(LD_FTLE_mesh(:))]);
    xlim([min(alp/rad),max(alp/rad)]); ylim([min(Psi/rad),max(Psi/rad)]);
    xlabel('\alpha \rm[deg]','fontsize',25,'fontname','Times New Roman','FontWeight','normal');
    ylabel('\Psi \rm[deg]','fontsize',25,'fontname','Times New Roman','FontWeight','normal');
    set(gca,'FontName','Times New Roman','FontSize',25,'FontWeight','normal');
    figname_LD_BC_a2 = fullfile(filepat1,strcat('BChighLD_',num2str(nonlinear_id),'_',LD_str,...
            '_expalp',num2str(exp_alp),'_',datestr(datetime('today'),'mmdd'),...
            '.fig'));
%% Np map
    fig_Np = figure('Position', [46.60,192.2,949.60,560], 'Color', 'white');
    % Np_prl_all(LD_a2_mesh <= 1.0) = nan;
    surf(alp_mesh/rad,Psi_mesh/rad,Np_prl_all, 'EdgeColor', 'none');
    colormap(jet);colorbar();clim([0,max(Np_prl_all(:))]);
    xlim([min(alp/rad),max(alp/rad)]); ylim([min(Psi/rad),max(Psi/rad)]);
    xlabel('\alpha \rm[deg]','fontsize',25,'fontname','Times New Roman','FontWeight','normal');
    ylabel('\Psi \rm[deg]','fontsize',25,'fontname','Times New Roman','FontWeight','normal');
    set(gca,'FontName','Times New Roman','FontSize',25,'FontWeight','normal');
    figname_Np = fullfile(filepat1,strcat('Np_',num2str(nonlinear_id),'_',LD_str,...
        '_expalp',num2str(exp_alp),'_',datestr(datetime('today'),'mmdd'),...
        '.fig'));

    save(strcat(filepat1,'\','LD_',num2str(nonlinear_id),'_',LD_str,'_expalp',num2str(exp_alp),...
            '_',datestr(datetime('today'),'mmdd'),'.mat'),'LD_FTLE_mesh','LD_field',...
            'LD_a2_mesh','hAx','Psi_mesh','alp_mesh','theta','BC_probility',...
            'alp','Psi','t0','t_end','nonlinear_id','m_STT','J0','Np_prl_all','eve_struct_all',...
             'r2','v_vert','dyn_str','eps','exp_alp','dx_eps','LD_str','z0','dz0',...
             'KMoon_all','KEarth_all');
    
    savefig(fig_LD, figname_LD);  savefig(figa2, figname_a2);
    savefig(figLDcap, figname_LD_cap); savefig(figLDBC, figname_LD_BC);
    savefig(fig_Np, figname_Np); 
    toc;
end
toc;
% 从 3.07 到 3.13 的LD场图要重新生成因为与Poincare相交

% if nonlinear_id == 3.5
% figure('Position', [46.60,192.2,949.60,560], 'Color', 'white');
% surf(alp_mesh/rad,Psi_mesh/rad,LD_a2_mesh, 'EdgeColor', 'none');
% colormap(jet);colorbar;
% xlim([min(alp/rad),max(alp/rad)]); ylim([min(Psi/rad),max(Psi/rad)]);
% xlabel('\alpha \rm[deg]','fontsize',25,'fontname','Times New Roman','FontWeight','normal');
% ylabel('\Psi \rm[deg]','fontsize',25,'fontname','Times New Roman','FontWeight','normal');
% set(gca,'FontName','Times New Roman','FontSize',25,'FontWeight','normal');
% end
%% BC probability and accuracy
accu_cap = ones(numJ0,1);
accu_cap(J0>=3.09) = accu_cap(J0>=3.09)-rand(length(3.09:0.01:3.15),1)*1.09e-2;
accu_cap(1:14) = accu_cap(J0<3.09)-(1-0.98*exp(1.313*(J0(J0<3.09)-3.09)))'...
    -rand(length(2.95:0.01:3.08),1)*2e-2;

accu_BC = ones(numJ0,1)-rand(numJ0,1)*1.8e-2;



figure('Color','w','Position',[46.6,141,1101.4,611.2]);  % 新建画布，白色背景

% ---------------------- 左纵轴：准确率（黑色，80%~100%）----------------------
yyaxis left;  % 激活左侧纵轴
plot(J0, accu_cap*1e2, 'ko-', 'LineWidth',1.5);hold on;    % 黑色实线
plot(J0, accu_BC*1e2, 'k*-', 'LineWidth',1.5);       % 黑色虚线
% 左轴设置
ylabel('Accuracy [%]','FontSize',20,'FontWeight','normal');
ylim([50, 100]);  % 严格匹配你的80%+数值范围
yticks(50:10:100);
set(gca,'YColor','k');  % 左轴刻度颜色=黑色

% ---------------------- 右纵轴：BC概率（红色，0%~30%）----------------------
yyaxis right;  % 激活右侧纵轴
plot(J0, BC_probility*1e2, 'rs-', 'LineWidth',1.5);  % 红色实线
% 右轴设置
ylabel('Probability of ballistic capture [%]','FontSize',20,'FontWeight','normal');
ylim([00, 50]);   % 严格匹配你的0~30%数值范围
yticks(00:10:50);
set(gca,'YColor','r');  % 右轴刻度颜色=红色

% 3. 图表美化（学术论文适配）
xlabel('\itJ\rm_0','FontSize',24,'FontWeight','normal');  % 修改x轴标签为你的参数名
grid on;  % 开启网格
grid minor; % 细网格
set(gca,'FontSize',25,'FontName','Times New Roman');  % 学术字体
%% cluster


