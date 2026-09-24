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
LL=get_LL(mu,mu_EMS,TU);
%% add STT & tensor path
addpath('D:\matlab_program\cyj_program\CRTBP\STT');
addpath(genpath('D:\MATLAB2023a\R2023a\tensor_toolbox-v3.8'));
%% de2/dt equation
tic;
% H0 = -847.500; % Initial Hamiltonian in BCRFBP
J0 = 3.125;
eps=1e-3;
theta_s = pi/2;
t0 = 0;t_end = 4.5;

% z component setting
r2 = 0.35;
z0 = 0;dz0 = -0e-1;
% r2*sin(pi/6)
%% hill region
[x,y]=meshgrid(linspace(-0.4,0.4,1000));
x = x+1.0;
rr1=((x+mu).^2+y.^2+z0^2).^0.5;
rr2=((x+mu-1).^2+y.^2+z0^2).^0.5;
rS = ((x-rho_sun*cos(theta_s)).^2+(y-rho_sun*sin(theta_s)).^2).^0.5;
H = -(dz0^2)+(x.^2+y.^2)+2*(1-mu)./rr1+2*mu./rr2+1*mu*(1-mu);

figure();
v=J0*[1, 1];
contour(x,y,H,v,'k-','linewidth',1);hold on;grid on;
% contourf(x,y,H,v)
% colormap(gray)

axis([-2,2,-1.5,1.5])
plot(1-mu,0,'r.','MarkerSize',16);
plot(-mu,0,'b.','MarkerSize',16); % earth
scatter(LL.L1(1),LL.L1(2),50,'k.');
scatter(LL.L2(1),LL.L2(2),50,'k.');
set(gca,'FontName','Times New Roman','FontSize',10);
xlabel('\itx \rm[DU]','fontsize',24,'fontname','Times New Roman','FontWeight','normal');
ylabel('\ity \rm[DU]','fontsize',24,'fontname','Times New Roman','FontWeight','normal');
set(gca,'FontName','Times New Roman','FontSize',24,'FontWeight','normal');
%% 
tic
theta = asin(z0/r2);
alp = [0:0.5:360]*rad;
% ballistic angle
Psi = [0:0.5:360]*rad;

num_all = length(alp)*length(Psi);
[alp_mesh,Psi_mesh] = meshgrid(alp,Psi);
LD_FTLE_mesh = zeros(size(alp_mesh));
LD_a2_mesh = zeros(size(alp_mesh));
options_LD = odeset('RelTol',1e-9,'AbsTol',1e-9,'Events',@only_hit_LD);
% v-planar
v_plan = zeros(size(alp_mesh));
v_vert = dz0;

dx_eps = 1e-4;
nonlinear_id = 3.5;
% 1--FTLE; 2--HOFTLE(.1-STT; .2-DSTT; .3-TDSTT); 
% 3--LD(.1-E2; .2-H2; .3-v2; .4-a2; .5-H2|E2|r2|v2|onlyv(v-norm) & a2)
exp_alp = 1.25;
m_STT = 2; % m order STT
dyn_str = 'EM'; % 'EM' 'diff_sbcm'
LD_str = 'H2';

Np_prl_all = zeros(length(alp),length(Psi));
eve_struct_all = cell(length(alp),length(Psi));
parfor i = 1:num_all
    r2_vec = r2*[cos(alp_mesh(i))*cos(theta);sin(alp_mesh(i))*cos(theta);sin(theta)];
    r1_vec = [1;0;0]+r2_vec;
    r1 = norm(r1_vec);
    
    vv = fsolve(@(v)(1*(v^2+v_vert^2)-((1-mu+r2_vec(1))^2+r2_vec(2)^2) ...
        -2*(1-mu)./r1-2*mu./r2-mu*(1-mu)+J0),1e-2,option_LL);
    Np_prl=0;eve_struct=struct('tend',0,'xend',0,'eve',0);
    if vv >= 1e-3
        v_plan(i) = vv;
        sig = alp_mesh(i) + Psi_mesh(i);
        y0 = [[r2_vec+[1-mu;0;0]];v_plan(i)*cos(sig);v_plan(i)*sin(sig);v_vert];
        if nonlinear_id >1 && nonlinear_id < 3
            [orbDEX,Np_prl,eve_struct,t,X,tend,xend,eve] = get_Dyindex(t0,t_end,y0,nonlinear_id,m_STT,...
            dyn_str,eps,exp_alp,0,dx_eps,LD_str,options_LD,mu);
        else
            [orbDEX,Np_prl,eve_struct,tend,xend,eve] = get_Dyindex_1_3_mex(t0,t_end,y0,nonlinear_id,m_STT,...
            dyn_str,eps,exp_alp,0,dx_eps,LD_str,mu);
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
toc

LD_FTLE_mesh(LD_FTLE_mesh == 0) = nan;% 有效 FTLE map
LD_a2_mesh(LD_a2_mesh == 0) = nan; % 有效a2 LD map
Np_prl_all(isnan(LD_a2_mesh)) = nan;
Np_prl_all((LD_a2_mesh <= 1.0)) = nan;% 有效近月点map
%% plot LD field
filepat1 = strcat('J',num2str(J0),'_r2',num2str(r2),'_z0',num2str(z0),...
        '_dz0',num2str(dz0));
if ~exist(filepat1, 'dir')
mkdir(filepat1);
end

fig = figure('Position', [46.60,192.2,949.60,560], 'Color', 'white');
hAx = axes;
LD_field=surf(hAx,alp_mesh/rad,Psi_mesh/rad...
    ,LD_FTLE_mesh, 'EdgeColor', 'none');
colormap(hAx,jet);c=colorbar(hAx);
xlim([min(alp/rad),max(alp/rad)]); ylim([min(Psi/rad),max(Psi/rad)]);
view(2);
xlabel('\alpha \rm[deg]','fontsize',25,'fontname','Times New Roman','FontWeight','normal');
ylabel('\Psi \rm[deg]','fontsize',25,'fontname','Times New Roman','FontWeight','normal');
set(gca,'FontName','Times New Roman','FontSize',25,'FontWeight','normal');
set(LD_field,'ButtonDownFcn',{@mycall_map,hAx,...
        t0,t0+t_end,alp_mesh,Psi_mesh,LD_FTLE_mesh,r2,theta,v_vert,...
        J0,nonlinear_id,m_STT,dyn_str,eps,exp_alp,dx_eps,LD_str});
figname_LD = fullfile(filepat1,strcat('LD_',num2str(nonlinear_id),'_',LD_str,...
                '_T',num2str(t_end),'_expalp',num2str(exp_alp),'_',datestr(datetime('today'),'mmdd'),'.fig'));

%% save 
save(strcat(filepat1,'\','LD_',num2str(nonlinear_id),'_',LD_str,'_expalp',num2str(exp_alp),...
    '_changeT_',datestr(datetime('today'),'mmdd'),'.mat'),'LD_FTLE_mesh','LD_field',...
        'LD_a2_mesh','hAx','Psi_mesh','alp_mesh','theta',...
        'alp','Psi','t0','t_end','nonlinear_id','m_STT','J0','Np_prl_all','eve_struct_all',...
         'r2','v_vert','dyn_str','eps','exp_alp','dx_eps','LD_str','z0','dz0');

if nonlinear_id == 3.5
figa2=figure('Position', [46.60,192.2,949.60,560], 'Color', 'white');
surf(alp_mesh/rad,Psi_mesh/rad,LD_a2_mesh, 'EdgeColor', 'none');
colormap(jet);colorbar;
xlim([min(alp/rad),max(alp/rad)]); ylim([min(Psi/rad),max(Psi/rad)]);
xlabel('\alpha \rm[deg]','fontsize',25,'fontname','Times New Roman','FontWeight','normal');
ylabel('\Psi \rm[deg]','fontsize',25,'fontname','Times New Roman','FontWeight','normal');
set(gca,'FontName','Times New Roman','FontSize',25,'FontWeight','normal');
view(2);
figname_a2 = fullfile(filepat1,strcat('LDa2_T',num2str(t_end),'_expalp',...
    num2str(exp_alp),'_',datestr(datetime('today'),'mmdd'),'.fig'));
savefig(figa2, figname_a2);
savefig(fig, figname_LD);
end
%% if nonlinear_id == 3.5   
LD_FTLE_fixfroma2 = LD_FTLE_mesh;
LD_FTLE_fixfroma2(LD_a2_mesh <= 1.0) = nan;
figure('Position', [46.60,192.2,949.60,560], 'Color', 'white');
surf(alp_mesh/rad,Psi_mesh/rad,LD_FTLE_fixfroma2, 'EdgeColor', 'none');
colormap(jet);colorbar();
xlim([min(alp/rad),max(alp/rad)]); ylim([min(Psi/rad),max(Psi/rad)]);
xlabel('\alpha \rm[deg]','fontsize',25,'fontname','Times New Roman','FontWeight','normal');
ylabel('\Psi \rm[deg]','fontsize',25,'fontname','Times New Roman','FontWeight','normal');
set(gca,'FontName','Times New Roman','FontSize',25,'FontWeight','normal');

LD_FTLE_fixfroma3 = LD_FTLE_fixfroma2;
LD_FTLE_fixfroma3(LD_FTLE_fixfroma2 <= 3) = nan;
figure('Position', [46.60,192.2,949.60,560], 'Color', 'white');
surf(alp_mesh/rad,Psi_mesh/rad,LD_FTLE_fixfroma3, 'EdgeColor', 'none');
colormap(jet);colorbar();clim([min(LD_FTLE_mesh(:)),max(LD_FTLE_mesh(:))]);
xlim([min(alp/rad),max(alp/rad)]); ylim([min(Psi/rad),max(Psi/rad)]);
xlabel('\alpha \rm[deg]','fontsize',25,'fontname','Times New Roman','FontWeight','normal');
ylabel('\Psi \rm[deg]','fontsize',25,'fontname','Times New Roman','FontWeight','normal');
set(gca,'FontName','Times New Roman','FontSize',25,'FontWeight','normal');
%% polar fig
LD_FTLE_fixfroma2 = LD_FTLE_mesh;
% LD_FTLE_fixfroma2(LD_a2_mesh <= 1.0) = nan;
% 1. 将网格数据 (矩阵) 转换为向量
theta_vec = alp_mesh(:);   % 极角，转为弧度
rho_vec = Psi_mesh(:) / rad;     % 极径
c_vec = LD_FTLE_fixfroma2(:);    % 颜色值
LD_a2_vec = LD_a2_mesh(:);

theta_vec(LD_a2_vec<=1.0)=[];
rho_vec(LD_a2_vec<=1.0)=[];c_vec(LD_a2_vec<=1.0)=[];
% 2. 创建极坐标散点图
figure('Position',[511.5,155,671,549.5],'Color', 'white');
pax = polaraxes;
hold(pax, 'on');
% 关键: 将颜色数据 c_vec 作为第4个参数传入，即可实现基于数值的颜色映射[reference:0][reference:1]
% 'filled' 参数让标记填充颜色[reference:2]
polarscatter(pax, theta_vec, rho_vec, ...
             3, c_vec, ...
             'filled', ...
             'MarkerFaceAlpha', 0.85, ...
             'MarkerEdgeColor', 'none');

% 3. 设置颜色和颜色条
colormap(jet);          % 设置颜色图为 jet
colorbar;               % 显示颜色条
pax.RTick = [];
pax.ThetaTick = 0:30:330;
pax.ThetaTickLabel = {'0^\circ','30^\circ','60^\circ','90^\circ', ...
                      '120^\circ','150^\circ','180^\circ','210^\circ', ...
                      '240^\circ','270^\circ','300^\circ','330^\circ'};

pax.ThetaZeroLocation = 'right';
pax.ThetaDir = 'counterclockwise';

set(gca,'FontName','Times New Roman','FontSize',25,'FontWeight','normal');
%% binary plot
figure;
LD_FTLE_mesh2 = LD_FTLE_mesh;
LD_FTLE_mesh2(LD_a2_mesh<=1.0) = 0;
LD_FTLE_mesh2(LD_a2_mesh>1.0) = 3;
scatter3(alp_mesh(:)/rad, Psi_mesh(:)/rad, LD_FTLE_mesh2(:), 6, LD_FTLE_mesh(:), 'filled');
colormap(jet);colorbar;
ylim([0,360]);yticks(0:90:360);
xlim([0,360]);xticks(0:90:360);
clim([min(LD_FTLE_mesh(:)), max(LD_FTLE_mesh(:))]);
grid on;
xlabel('\alpha \rm[deg]','fontsize',25,'fontname','Times New Roman','FontWeight','normal');
ylabel('\Psi \rm[deg]','fontsize',25,'fontname','Times New Roman','FontWeight','normal');
% zlabel('LD value','fontsize',25,'fontname','Times New Roman','FontWeight','normal');
set(gca,'FontName','Times New Roman','FontSize',25,'FontWeight','normal');


%% Np map
figure('Position', [46.60,192.2,949.60,560], 'Color', 'white');
% Np_prl_all(LD_a2_mesh <= 1.0) = nan;

surf(alp_mesh/rad,Psi_mesh/rad,Np_prl_all, 'EdgeColor', 'none');
colormap(jet);colorbar();clim([0,max(Np_prl_all(:))]);
xlim([min(alp/rad),max(alp/rad)]); ylim([min(Psi/rad),max(Psi/rad)]);
xlabel('\alpha \rm[deg]','fontsize',25,'fontname','Times New Roman','FontWeight','normal');
ylabel('\Psi \rm[deg]','fontsize',25,'fontname','Times New Roman','FontWeight','normal');
set(gca,'FontName','Times New Roman','FontSize',25,'FontWeight','normal');

%% only compare HOFTLE & FTLE
HO_FTLE_mesh = LD_FTLE_mesh;
load('J3.125_r20.25_z2e-3_dz3e-3\FTLE_1_E2_expalp1.25.mat','LD_FTLE_mesh');

figure('Position', [46.60,192.2,949.60,560], 'Color', 'white');
LD_field=surf(alp_mesh/rad,Psi_mesh/rad...
    ,HO_FTLE_mesh-LD_FTLE_mesh, 'EdgeColor', 'none');
colormap(jet);colorbar();
xlim([min(alp/rad),max(alp/rad)]); ylim([min(Psi/rad),max(Psi/rad)]);
xlabel('\alpha \rm[deg]','fontsize',25,'fontname','Times New Roman','FontWeight','normal');
ylabel('\Psi \rm[deg]','fontsize',25,'fontname','Times New Roman','FontWeight','normal');
set(gca,'FontName','Times New Roman','FontSize',25,'FontWeight','normal');


%% cluster
% 目前的聚类分析还不完善
[Z_continuous, idx, bestK, posit, points_LD] = ...
    class_kmeans_adaptive(LD_FTLE_fixfroma2,alp_mesh,Psi_mesh,4,10);

figure('Position', [46.60,192.2,949.60,560]); 
scatter(points_LD(:,1)/rad, points_LD(:,2)/rad, 10, idx, 'filled');
colormap(parula(bestK)); % 用不同颜色区分类别
c2=colorbar;
%c2.Label.String='Cluster';c2.Label.Rotation=0;
%c2.Label.Position=[0.383068757166193,bestK+1/bestK,0];% [0.012698387153564,13.228855570588351,0]; % 0.383068757166193,27.72885608673097,0
clim([1,bestK]);grid on;
xlim([0,360]);ylim([0,360]);
xticks(0:45:360);yticks(0:45:360);
xlabel('\alpha \rm[deg]','fontsize',25,'fontname','Times New Roman','FontWeight','normal');
ylabel('\Psi \rm[deg]','fontsize',25,'fontname','Times New Roman','FontWeight','normal');
set(gca,'FontName','Times New Roman','FontSize',25,'FontWeight','normal');


%% save
filepat = strcat('J',num2str(J0),'_r20.25_z0',num2str(z0),...
    '_dz0',num2str(dz0),'\');
save(strcat(filepat,'LD_',num2str(nonlinear_id),'_',LD_str,'_expalp',num2str(exp_alp),...
        '_',datestr(datetime('today'),'mmdd'),'.mat'),'LD_FTLE_mesh','LD_field',...
        'LD_a2_mesh','hAx','Psi_mesh','alp_mesh','theta',...
        'alp','Psi','t0','t_end','nonlinear_id','m_STT','J0','Np_prl_all','eve_struct_all',...
         'r2','v_vert','dyn_str','eps','exp_alp','dx_eps','LD_str','z0','dz0');



