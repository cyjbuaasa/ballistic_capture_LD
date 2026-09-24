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
DatL1 = load('D:\matlab_program\cyj_program\CRTBP\periodic_orb_2506\Lyap_bifurcation\Lyapunov_planar\Lyapunov_L1_honrizontal_initial_data.mat');
DatL2 = load('D:\matlab_program\cyj_program\CRTBP\periodic_orb_2506\Lyap_bifurcation\Lyapunov_planar\Lyapunov_L2_honrizontal_initial_data.mat');
addpath('D:\matlab_program\cyj_program\CRTBP\periodic_orb_2506\Lyap_bifurcation\Lyapunov_planar');
%% de2/dt equation
tic;
% H0 = -847.500; % Initial Hamiltonian in BCRFBP
J0 = 3.125;
eps=1e-3;
theta_s = pi/2;
t0 = 0;t_end = 4.5;
% z component setting
r2 = 0.25;
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
% axis([-2,2,-1.5,1.5])
axis equal;% xlim([-0.4,0.4])
plot(1-mu,0,'r.','MarkerSize',16);
plot(-mu,0,'b.','MarkerSize',16); % earth
scatter(LL.L1(1),LL.L1(2),50,'k.');
scatter(LL.L2(1),LL.L2(2),50,'k.');
set(gca,'FontName','Times New Roman','FontSize',10);
xlabel('\itx \rm[DU]','fontsize',24,'fontname','Times New Roman','FontWeight','normal');
ylabel('\ity \rm[DU]','fontsize',24,'fontname','Times New Roman','FontWeight','normal');
set(gca,'FontName','Times New Roman','FontSize',22,'FontWeight','normal');
%% 
tic
theta = asin(z0/r2);
alp = [0:0.1:360]*rad;
flag = nan(size(alp));
for i = 1:length(alp)
    r2_vec = r2*[cos(alp(i));sin(alp(i));0];
    r1_vec = [1;0;0]+r2_vec;
    r1 = norm(r1_vec);
    
    vv = fsolve(@(v)(1*(v^2+0^2)-((1-mu+r2_vec(1))^2+r2_vec(2)^2) ...
        -2*(1-mu)./r1-2*mu./r2-mu*(1-mu)+J0),1e-2,option_LL);
    if vv >= 1e-3
        flag(i) = 1;
    end
end

% r2=0.35
% plot(1-mu+r2*cos(alp(flag==1)),r2*sin(alp(flag==1)),'r.','linewidth',2.0);
% r2=0.25
plot(1-mu+r2*cos(alp(flag==1)),r2*sin(alp(flag==1)),'m.','linewidth',2.0);
%% L1 || L2
typ = 'L2';
%% Lyapnov orb continuation
options = odeset('RelTol',1e-12,'AbsTol',1e-12);
if strcmp(typ,'L1')
    % L1
    [Jic,jmin] = min(abs(DatL1.C_family-J0));
    x0 = DatL1.x_start(1:6,jmin);
    T0 = DatL1.x_start(7,jmin);
else
    % L2
    [Jic,jmin] = min(abs(DatL2.C_family-J0));
    x0 = DatL2.x_start(1:6,jmin);
    T0 = DatL2.x_start(7,jmin);
end
X=singleshoot_half_period_C([x0(1),x0(5),J0,T0], mu);
tspan = [0 2*X(end)];
x0_shot = [X(1);0;0;0;X(2);0];
[t1,xx] = ode45(@(t,x)orbitdym(t, x, mu),tspan,x0_shot,options);
% figure;
% plot3(xx(:,1),xx(:,2),xx(:,3),'g-','LineWidth',2); hold on;
% axis equal;grid on;
%% manifold 
%截取若干点
[la,lb]=size(xx);
n=50;x_E = zeros(n,6);
for i=1:n
    j=round(la/n*i);
    x_E(i,:)=xx(j,:);
end
r_tol = 1e-6;
% figure;
for i=1:n
    I=eye(6);
    stm0=reshape(I,[],1);
    x1=[x_E(i,:)';stm0];
    [t1_C,xx_C] = ode45(@(t,x)dym(t, x, mu),tspan,x1,options);
    %单值矩阵
    stm=[xx_C(end,(7:12))' xx_C(end,(13:18))' xx_C(end,(19:24))' ...
        xx_C(end,(25:30))' xx_C(end,(31:36))' xx_C(end,(37:42))'];
    %特征值
    ee=eig(stm);
    %特征值对应的特征向量
    [aa,~]=eig(stm);
    
    idx_s=find(ee(:,1)==min(ee));
    E_s=aa(:,idx_s);
    
    V_s=E_s(1:6)/norm(E_s(1:3),2);
    if V_s(1)<0
    V_s=-V_s;
    end
    
    %绘制稳定流形
    if strcmp(typ,'L1')
        x1_s=x_E(i,:)'-r_tol*V_s;
    else
        x1_s=x_E(i,:)'+r_tol*V_s;
    end
    I=eye(6);
    stm0=reshape(I,[],1);
    % x1_s=[x1_s;stm0];
    Tspan=-50*tspan;
    options=odeset('RelTol',1e-12,'AbsTol',1e-12,'Events',@onsection);
    [t1,xx_s,tend,xend,eve]= ode45(@orbitdym,Tspan,x1_s,options,mu); % dym
    
    plot(xx_s(:,1),xx_s(:,2),'b');grid on;hold on;
end
%% event func
function [value,isterminal,direction] = onsection(t,X,~,~,~,~,~)
%%% onsection detection stop

    u=0.0121506683;
    DU=3.84405*10^5;
    r2 = 0.35;
    Re=[-u;0;0];                                           
    Rm=[1-u;0;0];
    Rs=[X(1);X(2);X(3)];
    rEM_hill = 1.501292597582450e+06*2/DU;

    R2=Rs-Rm;
    R1=Rs-Re;
    % 计算地心距和径向速度
    v1 = X(4:6)-[X(2);-X(1)-u;0]; % 惯性系
    v2 = X(4:6)-[X(2);-X(1)+1-u;0];
    value = [norm(R2)-r2;];% norm(R2)-(6*66180/DU)
        isterminal = [1;]; 
        direction = [1];

    % value = [norm(R2)-((1738)/DU);norm(R1)-((6378)/DU);radial_velocity_m];
    % isterminal = [1;1;0]; 
    % direction = [0;0;1];   % all direction
    % % 3:perilune
    % 1 -- 正向穿越; -1 -- 负向穿越; 

end



