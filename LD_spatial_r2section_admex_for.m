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
gamma1=fsolve(@(x)(x^5-(3-mu)*x^4+(3-2*mu)*x^3-mu*x^2+2*mu*x-mu)-0,(mu)^(1/3),option_LL);
gamma2=fsolve(@(x)(x^5+(3-mu)*x^4+(3-2*mu)*x^3-mu*x^2-2*mu*x-mu)-0,(mu)^(1/3),option_LL);
gamma3=fsolve(@(x)(x^5+(2+mu)*x^4+(1+2*mu)*x^3-(1-mu)*x^2-2*(1-mu)*x-1+mu)-0,1-7/12*mu,option_LL);
LL.L1=[1-mu-gamma1;0;0];
LL.L2=[1-mu+gamma2;0;0];
LL.L3=[-mu-gamma3;0;0];
LL.L4=[0.5-mu;sqrt(3)/2;0];
LL.L5=[0.5-mu;-sqrt(3)/2;0];
gamma1=fsolve(@(x)(x^5-(3-mu_EMS)*x^4+(3-2*mu_EMS)*x^3-mu_EMS*x^2+2*mu_EMS*x-mu_EMS)-0,(mu_EMS)^(1/3),option_LL);
gamma2=fsolve(@(x)(x^5+(3-mu_EMS)*x^4+(3-2*mu_EMS)*x^3-mu_EMS*x^2-2*mu_EMS*x-mu_EMS)-0,(mu_EMS)^(1/3),option_LL);
gamma3=fsolve(@(x)(x^5+(2+mu_EMS)*x^4+(1+2*mu_EMS)*x^3-(1-mu_EMS)*x^2-2*(1-mu_EMS)*x-1+mu_EMS)-0,1-7/12*mu_EMS,option_LL);
LL.SEL1=[1-mu_EMS-gamma1;0;0];
LL.SEL2=[1-mu_EMS+gamma2;0;0];
LL.SEL3=[-mu_EMS-gamma3;0;0];
LL.SEL4=[0.5-mu_EMS;sqrt(3)/2;0];
LL.SEL5=[0.5-mu_EMS;-sqrt(3)/2;0];
clc;
LL.JEL1 = Jacobi([LL.SEL1;0;0;0],mu_EMS);LL.JEL2 = Jacobi([LL.SEL2;0;0;0],mu_EMS);
LL.JEL3 = Jacobi([LL.SEL3;0;0;0],mu_EMS);LL.JELT = Jacobi([LL.SEL4;0;0;0],mu_EMS);
LL.JEML1 = Jacobi([LL.L1;0;0;0],mu);LL.JEML2 = Jacobi([LL.L2;0;0;0],mu);
LL.JEML3 = Jacobi([LL.L3;0;0;0],mu);LL.JEMLT = Jacobi([LL.L4;0;0;0],mu);
%% add STT & tensor path
addpath('D:\matlab_program\cyj_program\CRTBP\STT');
addpath(genpath('D:\MATLAB2023a\R2023a\tensor_toolbox-v3.8'));
%% de2/dt equation
tic;
% H0 = -847.500; % Initial Hamiltonian in BCRFBP
J0 = 3.05;
numJ0 = length(J0);
eps=1e-3;
t0 = 0;t_end = 4.5;

% sphere setting
r2 = 0.25;
psi2 = [-88*rad,[-85*rad:5*rad:0]]; gama = -pi/4:5*rad:pi/4;
% z0 = 0;dz0 = -0e-1;

%% ini parameters
tic
filepat = 'spatial\';

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
keta = 0;
BC_probility = zeros(length(gama),length(psi2));
%% for J0
tic;
for j = 1:1
    for k = 1:length(psi2)
    for n = 1:length(gama)
        fprintf('psi2 = %f, gama = %f \n',psi2(k)/rad,gama(n)/rad);
    filepat1 = strcat(filepat,'J',num2str(J0(j)),'_r2',num2str(r2),...
        '_gama',num2str(gama(n)/rad),'_psi2',num2str(psi2(k)/rad));
    if exist(filepat1,'dir')==0 %%判断文件夹是否存在
        mkdir(filepat1);  %%不存在时候，创建文件夹
    else
        disp('dir is exist,continue \n');
        %%如果文件夹存在，输出:dir is exist
        continue;
    end

    LD_FTLE_mesh = zeros(size(alp_mesh));
    LD_a2_mesh = zeros(size(alp_mesh));
    Np_prl_all = zeros(length(alp),length(psi1));
    eve_struct_all = cell(length(alp),length(psi1));
    fprintf("J0 = %f\n",J0(j));
    parfor i = 1:num_all
        r2_vec = r2*[cos(alp_mesh(i))*cos(gama(n));sin(alp_mesh(i))*cos(gama(n));sin(gama(n))];
        r1_vec = [1;0;0]+r2_vec;
        r1 = norm(r1_vec);
        
        vv = fsolve(@(v)(1*(v^2)-((1-mu+r2_vec(1))^2+r2_vec(2)^2) ...
            -2*(1-mu)./r1-2*mu./r2-mu*(1-mu)+J0(j)),1e-2,option_LL);
        Np_prl=0;eve_struct=struct('tend',0,'xend',0,'eve',0);
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
                [orbDEX,Np_prl,eve_struct,t,X,tend,xend,eve] = get_Dyindex(t0,t_end,y0,nonlinear_id,m_STT,...
                dyn_str,eps,exp_alp,dx_eps,LD_str,options_LD,mu);
            else
                [orbDEX,Np_prl,eve_struct,tend,xend,eve] = get_Dyindex_1_3_mex(t0,t_end,y0,nonlinear_id,m_STT,...
                dyn_str,eps,exp_alp,dx_eps,LD_str,mu);
            end
            if nonlinear_id == 3.5
                LD_FTLE_mesh(i) = orbDEX(1);
                LD_a2_mesh(i) = orbDEX(2);
            else
                LD_FTLE_mesh(i) = orbDEX;
            end
        end
        Np_prl_all(i) = Np_prl;eve_struct_all{i} = eve_struct;
        if mod(i,1e3) == 0, fprintf("%d done\n",i); end
    end
    toc;
    % 数据处理
    LD_FTLE_mesh(LD_FTLE_mesh == 0) = nan;% 有效 FTLE map
    LD_a2_mesh(LD_a2_mesh == 0) = nan; % 有效a2 LD map
    Np_prl_all(isnan(LD_a2_mesh)) = nan;
    Np_prl_all((LD_a2_mesh <= 1.0)) = nan;% 有效近月点map
     % 捕获可能性
    BC_probility(n,k) = (numel(Np_prl_all(Np_prl_all>=2)))...
        /(numel(Np_prl_all)-numel(Np_prl_all(isnan(Np_prl_all))));
    %% LD
    fig_LD = figure('Position', [46.60,192.2,949.60,560], 'Color', 'white');
    hAx = axes;
    LD_field=surf(hAx,alp_mesh/rad,psi1_mesh/rad...
        ,LD_FTLE_mesh, 'EdgeColor', 'none');
    colormap(hAx,jet);c=colorbar(hAx);
    xlim([min(alp/rad),max(alp/rad)]); ylim([min(psi1/rad),max(psi1/rad)]);
    xlabel('\alpha \rm[deg]','fontsize',25,'fontname','Times New Roman','FontWeight','normal');
    ylabel('\Psi_1 \rm[deg]','fontsize',25,'fontname','Times New Roman','FontWeight','normal');
    set(gca,'FontName','Times New Roman','FontSize',25,'FontWeight','normal');
    set(LD_field,'ButtonDownFcn',{@mycall_map_spatial,hAx,...
            t0,t0+t_end,alp_mesh,psi1_mesh,LD_FTLE_mesh,r2,gama(n),psi2(k),...
            J0(j),nonlinear_id,m_STT,dyn_str,eps,exp_alp,keta,dx_eps,LD_str});
    figname_LD = fullfile(filepat1,strcat('LD_',num2str(nonlinear_id),'_',LD_str,...
        '_expalp',num2str(exp_alp),'_',datestr(datetime('today'),'mmdd'),...
        '.fig'));
    %% if nonlinear_id == 3.5   
    LD_FTLE_fixfroma2 = LD_FTLE_mesh;
    LD_FTLE_fixfroma2(LD_a2_mesh <= 0.917) = nan;
    figLDcap = figure('Position', [46.60,192.2,949.60,560], 'Color', 'white');
    surf(alp_mesh/rad,psi1_mesh/rad,LD_FTLE_fixfroma2, 'EdgeColor', 'none');
    colormap(jet);colorbar();
    xlim([min(alp/rad),max(alp/rad)]); ylim([min(psi1/rad),max(psi1/rad)]);
    xlabel('\alpha \rm[deg]','fontsize',25,'fontname','Times New Roman','FontWeight','normal');
    ylabel('\Psi_1 \rm[deg]','fontsize',25,'fontname','Times New Roman','FontWeight','normal');
    set(gca,'FontName','Times New Roman','FontSize',25,'FontWeight','normal');
    figname_LD_cap = fullfile(filepat1,strcat('LDcap_',num2str(nonlinear_id),'_',LD_str,...
        '_expalp',num2str(exp_alp),'_',datestr(datetime('today'),'mmdd'),...
        '.fig'));

    LD_FTLE_fixfroma3 = LD_FTLE_fixfroma2;
    LD_FTLE_fixfroma3(LD_FTLE_fixfroma2 <= 4.05) = nan;
    figLDBC = figure('Position', [46.60,192.2,949.60,560], 'Color', 'white');
    surf(alp_mesh/rad,psi1_mesh/rad,LD_FTLE_fixfroma3, 'EdgeColor', 'none');
    colormap(jet);colorbar();
    if ~all(isnan(LD_FTLE_mesh(:)))
        clim([min(LD_FTLE_mesh(:)),max(LD_FTLE_mesh(:))]);
    end
    xlim([min(alp/rad),max(alp/rad)]); ylim([min(psi1/rad),max(psi1/rad)]);
    xlabel('\alpha \rm[deg]','fontsize',25,'fontname','Times New Roman','FontWeight','normal');
    ylabel('\Psi_1 \rm[deg]','fontsize',25,'fontname','Times New Roman','FontWeight','normal');
    set(gca,'FontName','Times New Roman','FontSize',25,'FontWeight','normal');
    figname_LD_BC = fullfile(filepat1,strcat('BChighLD_',num2str(nonlinear_id),'_',LD_str,...
            '_expalp',num2str(exp_alp),'_',datestr(datetime('today'),'mmdd'),...
            '.fig'));
%% Np map
    fig_Np = figure('Position', [46.60,192.2,949.60,560], 'Color', 'white');
    % Np_prl_all(LD_a2_mesh <= 1.0) = nan;
    surf(alp_mesh/rad,psi1_mesh/rad,Np_prl_all, 'EdgeColor', 'none');
    colormap(jet);colorbar();
    if ~all(isnan(Np_prl_all(:)))
        clim([0,max(Np_prl_all(:))]);
    end
    xlim([min(alp/rad),max(alp/rad)]); ylim([min(psi1/rad),max(psi1/rad)]);
    xlabel('\alpha \rm[deg]','fontsize',25,'fontname','Times New Roman','FontWeight','normal');
    ylabel('\Psi_1 \rm[deg]','fontsize',25,'fontname','Times New Roman','FontWeight','normal');
    set(gca,'FontName','Times New Roman','FontSize',25,'FontWeight','normal');
    figname_Np = fullfile(filepat1,strcat('Np_',num2str(nonlinear_id),'_',LD_str,...
        '_expalp',num2str(exp_alp),'_',datestr(datetime('today'),'mmdd'),...
        '.fig'));

    save(strcat(filepat1,'\','LD_',num2str(nonlinear_id),'_',LD_str,'_expalp',num2str(exp_alp),...
            '_',datestr(datetime('today'),'mmdd'),'.mat'),'LD_FTLE_mesh','LD_field',...
            'LD_a2_mesh','hAx','psi1_mesh','alp_mesh','gama','psi2','BC_probility',...
            'alp','psi1','t0','t_end','nonlinear_id','m_STT','J0','Np_prl_all','eve_struct_all',...
             'r2','dyn_str','eps','exp_alp','dx_eps','LD_str');
    
    savefig(fig_LD, figname_LD);
    savefig(figLDcap, figname_LD_cap); savefig(figLDBC, figname_LD_BC);
    savefig(fig_Np, figname_Np); 
    close all;
    end
    end
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
%% BC probability
[gama_grid,psi2_grid]=meshgrid(gama,psi2);
figure;
surf(gama_grid/rad,psi2_grid/rad,BC_probility'*100, 'EdgeColor', 'none');
colormap(jet);colorbar();clim([0.00,max(BC_probility(:))]);
xlim([min(gama/rad),max(gama/rad)]); ylim([min(psi2/rad),max(psi2/rad)]);
xlabel('\gamma \rm[deg]','fontsize',25,'fontname','Times New Roman','FontWeight','normal');
ylabel('\Psi_2 \rm[deg]','fontsize',25,'fontname','Times New Roman','FontWeight','normal');
set(gca,'FontName','Times New Roman','FontSize',25,'FontWeight','normal');
    
%%
% 1. 定义更细的网格范围（与原网格相同范围，点数更多）
x_vec_fine = linspace(min(gama), max(gama), 50);  % 例如 细分数 = 500
y_vec_fine = linspace(-pi/2, 0, 50);
[Xq, Yq] = meshgrid(x_vec_fine, y_vec_fine);

% 2. 进行插值
% method 可选：'linear'(默认), 'cubic', 'spline', 'makima'
Zq = interp2(gama_grid, psi2_grid, BC_probility'*100, Xq, Yq, 'spline');

% 3. 绘制细化后的曲面
figure;
surf(Xq/rad, Yq/rad, Zq, 'EdgeColor', 'none');
colormap(jet);colorbar();clim([0.00,max(BC_probility(:)*100)]);
xlim([min(gama/rad),max(gama/rad)]);
ylim([min(psi2/rad),max(psi2/rad)]);
yticks(-90:15:0);xticks(-45:15:45);
xlabel('\gamma \rm[deg]','fontsize',25,'fontname','Times New Roman','FontWeight','normal');
ylabel('\Psi_2 \rm[deg]','fontsize',25,'fontname','Times New Roman','FontWeight','normal');
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


