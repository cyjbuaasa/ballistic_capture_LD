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
J0 = 3.201;
eps=1e-1;
theta_s = pi/2;x0 = 0.7;
t0 = 0;t_end = 5;

% z component setting
z0 = 0;% 2e-3
dz0 = 0;% 3e-3
r2 = 0.23;
%% hill region
[x,y]=meshgrid(linspace(-0.4,0.4,1000));
x = x+1.0-mu;
rr1=((x+mu).^2+y.^2+z0^2).^0.5;
rr2=((x+mu-1).^2+y.^2+z0^2).^0.5;
% rS = ((x-rho_sun*cos(theta_s)).^2+(y-rho_sun*sin(theta_s)).^2).^0.5;
H = -(dz0^2)+(x.^2+y.^2)+2*(1-mu)./rr1+2*mu./rr2+1*mu*(1-mu);
%%
theta = asin(z0/r2);
alp = [0:0.5:360]*rad;
sigma = 0:pi/100:2*pi;
% de2_dt = zeros(length(sigma),1);
v_vert = dz0;
% v-planar
v_plan = zeros(length(alp),1);
% ballistic angle
PPsi = zeros(length(alp),2);
for i = 1:length(alp)
    r2_vec = r2*[cos(alp(i))*cos(theta);sin(alp(i))*cos(theta);sin(theta)];
    r1_vec = [1;0;0]+r2_vec;
    r1 = norm(r1_vec);
    
    vv = fsolve(@(v)(1*(v^2+v_vert^2)-((1-mu+r2_vec(1))^2+r2_vec(2)^2) ...
        -2*(1-mu)./r1-2*mu./r2-mu*(1-mu)+J0),...
        1e-2,option_LL);
    if vv >= 1e-3
        v_plan(i) = vv;
    end
end

% 以 dE2/dt < 0 为准则并不全，得到的条件太少
figure();
v=J0*[1, 1];
contour(x,y,H,v,'k-','linewidth',1);hold on;grid on;axis equal;
plot(1-mu,0,'r.','MarkerSize',16);
% plot(-mu,0,'b.','MarkerSize',16); % earth
scatter(LL.L1(1),LL.L1(2),50,'k.');
scatter(LL.L2(1),LL.L2(2),50,'k.');
plot(1-mu+r2*cos(alp(v_plan>0))*cos(theta),r2*sin(alp(v_plan>0))*cos(theta),'m.','LineWidth',1.5);
set(gca,'FontName','Times New Roman','FontSize',10);
xlabel('\itx \rm[DU]','fontsize',24,'fontname','Times New Roman','FontWeight','normal');
ylabel('\ity \rm[DU]','fontsize',24,'fontname','Times New Roman','FontWeight','normal');
set(gca,'FontName','Times New Roman','FontSize',24,'FontWeight','normal');
%% validation
j = 100;
sig = alp(j)+pi/2+0.6;
r2_vec = r2*[cos(alp(j))*cos(theta);sin(alp(j))*cos(theta);sin(theta)];
Y0 = [[r2_vec+[1-mu;0;0]];v_plan(j)*cos(sig);v_plan(j)*sin(sig);v_vert];
Jacobi(Y0,mu);

[t,X,tend,xend,eve] = ode45(@orbitdym,[0,t_end]+t0,Y0,options_LD,mu);

E2 = zeros(length(t),1);H2 = zeros(length(t),1);
a2 = zeros(length(t),1);
for j = 1:length(t)
    E2(j)=E_secondary(X(j,1:6)',mu);
    H2(j) = norm(cross(X(j,1:3)'-[1-mu;0;0],X(j,4:6)'-[X(j,2);-X(j,1)+1-mu;0]));
    r=[((X(j,1)+mu)^2+X(j,2)^2+X(j,3)^2)^0.5;((X(j,1)-1+mu)^2+X(j,2)^2+X(j,3)^2)^0.5];
    a2(j) = norm([-mu*(X(j,1)-1+mu)/(r(2)^3);-mu*X(j,2)/(r(2)^3);-mu*X(j,3)/(r(2)^3)]);
end
figure;
subplot(2,1,1);plot(t,E2,'LineWidth',1.5);grid on;
xlabel('\itt \rm[TU]','fontsize',23,'fontname','Times New Roman','FontWeight','normal');
ylabel('\itE_s_e_c','fontsize',23,'fontname','Times New Roman','FontWeight','normal');
set(gca,'FontName','Times New Roman','FontSize',23,'FontWeight','normal');
subplot(2,1,2);plot(t,H2,'LineWidth',1.5);grid on;
xlabel('\itt \rm[TU]','fontsize',23,'fontname','Times New Roman','FontWeight','normal');
ylabel('\itH_s_e_c','fontsize',23,'fontname','Times New Roman','FontWeight','normal');
set(gca,'FontName','Times New Roman','FontSize',23,'FontWeight','normal');

figure;plot(t,a2,'LineWidth',1.5);grid on;
xlabel('\itt \rm[TU]','fontsize',23,'fontname','Times New Roman','FontWeight','normal');
ylabel('\ita_s_e_c','fontsize',23,'fontname','Times New Roman','FontWeight','normal');
set(gca,'FontName','Times New Roman','FontSize',23,'FontWeight','normal');

figure;plot3(X(:,1),X(:,2),X(:,3),'LineWidth',1.5);hold on;grid on;axis equal;
plot3(1-mu,0,0,'r.','MarkerSize',16);
plot3(-mu,0,0,'b.','MarkerSize',16); % earth
xlabel('\itx \rm[DU]','fontsize',23,'fontname','Times New Roman','FontWeight','normal');
ylabel('\ity \rm[DU]','fontsize',23,'fontname','Times New Roman','FontWeight','normal');
zlabel('\itz \rm[DU]','fontsize',23,'fontname','Times New Roman','FontWeight','normal');
set(gca,'FontName','Times New Roman','FontSize',23,'FontWeight','normal');









