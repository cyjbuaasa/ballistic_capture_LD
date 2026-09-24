function mycall_map(hObject,~,hAx,t0,t_end,alp_mesh,Psi_mesh,orbDEX,r2,theta,v_vert,J0,nonlinear_id,m_STT,...
    dyn_str,eps,alp,dx_eps,LD_str)
    % load('orbit_data_initial.mat');
    rad = pi/180;day = 86400;
    pt = get(hAx,'CurrentPoint');
    x = pt(1,1)*rad;
    y = pt(1,2)*rad;
    z = pt(1,3);
    DU=3.84405*10^5; % km
    TU=4.34811305; % day
    DUS = 1.494609474249150e+08;
    option_LL = optimoptions('fsolve','StepTolerance',1e-12,'FunctionTolerance',1e-12,'OptimalityTolerance',1e-12,'Display','none');
    keta=0;
    % 寻找与点击对象最近的目标值
    % d = abs(theta/rad - x)
    % d2 = abs(phi/rad - y);
    % [v,i] = min(d);[v2,i2] = min(d2);
    [v, idx] = min((alp_mesh(:)-x).^2 + (Psi_mesh(:)-y).^2);
    [row, col] = ind2sub(size(alp_mesh), idx);
    fprintf(' alp: %f\n psi: %f\n', alp_mesh(row, col)/rad, Psi_mesh(row, col)/rad);
    fprintf(' orbDex = %f\n ',orbDEX(row, col));
    % 小于3像素点再绘图防止误点击
    mu = 0.0121506683;
    ome_sun = -0.925195985;
    mu_EMS=3.003480569974130e-06;
    optionsode=odeset('RelTol',1e-10,'AbsTol',1e-10,'Events',@only_hit_simple);

    if v < 3
        r2_vec = r2*[cos(alp_mesh(row, col))*cos(theta);...
            sin(alp_mesh(row, col))*cos(theta);sin(theta)];
        r1_vec = [1;0;0]+r2_vec;
        r1 = norm(r1_vec);
        
        vv = fsolve(@(v)(1*(v^2+v_vert^2)-((1-mu+r2_vec(1))^2+r2_vec(2)^2) ...
            -2*(1-mu)./r1-2*mu./r2-mu*(1-mu)+J0),...
            1e-2,option_LL);
        sig = alp_mesh(row, col) + Psi_mesh(row, col);
        y0 = [[r2_vec+[1-mu;0;0]];vv*cos(sig);vv*sin(sig);v_vert];
        % %% DRO轨道
        % x0=initial.dro_2_1.RV;
        % T0=initial.dro_2_1.T;
        % tspan0=[0 T0];
        % [~,xx0] =ode45(@dym,tspan0,x0,optionsode);

        %% 转移轨道
        options_LD = odeset('RelTol',1e-9,'AbsTol',1e-9,'Events',@only_hit_new);
        [orb_dynam_index,Np_prl,eve_struct,t,X,tend,xend,eve] = get_Dyindex(t0,t_end*1.8,y0,nonlinear_id,m_STT,...
    dyn_str,eps,alp,keta,dx_eps,LD_str,options_LD,mu);
        [tend,xend,eve] = postprocess_xmoon_events(tend,xend,eve);
    %     [orb_dynam_index,Np_prl,eve_struct,t,X,tend,xend,eve] = get_Dyindex_1_3(t0,t_end*2.0,y0,nonlinear_id,m_STT,...
    % dyn_str,eps,alp,dx_eps,LD_str,mu);
        % alpha=0:0.01:2*pi;
        % xx1(:,1)=h_e/DU*cos(alpha)-mu;
        % xx1(:,2)=h_e/DU*sin(alpha);
        % xx3(:,1)=h_m/DU*cos(alpha)+1-mu;
        % xx3(:,2)=h_m/DU*sin(alpha);
        %% 绘制转移轨道
        figure;
        plot3(X(:,1),X(:,2),X(:,3), 'LineWidth', 1.5);
        hold on;
        plot3(1-mu+r2*cos(0:1e-2:2*pi),r2*sin(0:1e-2:2*pi),0*cos(0:1e-2:2*pi),'k--','LineWidth', 1.5);
        plot3(1-mu,0,0,'r.','MarkerSize',14); % moon;
        % plot3(-mu,0,0,'b.','MarkerSize',16); % earth
        % plot3(xx1(:,1),xx1(:,2),zeros(length(alpha),1),'k-', 'LineWidth', 1.5);
        % plot3(xx3(:,1),xx3(:,2),zeros(length(alpha),1),'k-', 'LineWidth', 1.5);
        % infor=['Earth-Moon rotating, ','{\theta}=',num2str(roundn(theta(i)/rad,-3)),'[deg], {\phi}=',num2str(roundn(phi(i)/rad,-3)),'[deg]'];
        % title(infor);hold off;
        xlabel('\itx \rm[DU]','fontsize',18,'fontname','Times New Roman','FontWeight','normal');
        ylabel('\ity \rm[DU]','fontsize',18,'fontname','Times New Roman','FontWeight','normal');
        zlabel('\itz \rm[DU]','fontsize',18,'fontname','Times New Roman','FontWeight','normal');
        set(gca,'FontName','Times New Roman','FontSize',24,'FontWeight','normal');
        box on;grid on;axis equal;
        
        E2 = zeros(length(t),1);H2 = zeros(length(t),1);
        a2 = zeros(length(t),1);X_secE = zeros(6,length(t));
        for j = 1:length(t)
            E2(j)=E_secondary(X(j,1:6)',mu);
            H2(j)=norm(cross([X(j,1)-1+mu;X(j,2);X(j,3)],...
                [X(j,4)-X(j,2);X(j,5)+X(j,1)+mu-1;X(j,6)]));
            r=[((X(j,1)+mu)^2+X(j,2)^2+X(j,3)^2)^0.5;((X(j,1)-1+mu)^2+X(j,2)^2+X(j,3)^2)^0.5];
            a2(j) = norm([-mu*(X(j,1)-1+mu)/(r(2)^3);-mu*X(j,2)/(r(2)^3);-mu*X(j,3)/(r(2)^3)]);
            X_secE(:,j)=trans_REM_secE(X(j,1:6)',t(j));
        end
        figure;plot(t,E2,'LineWidth',1.5);grid on;
        xlabel('\itt \rm[TU]','fontsize',23,'fontname','Times New Roman','FontWeight','normal');
        ylabel('\itE_s_e_c','fontsize',23,'fontname','Times New Roman','FontWeight','normal');
        set(gca,'FontName','Times New Roman','FontSize',23,'FontWeight','normal');
        figure;plot(t,H2,'LineWidth',1.5);grid on;
        xlabel('\itt \rm[TU]','fontsize',23,'fontname','Times New Roman','FontWeight','normal');
        ylabel('\itH_s_e_c','fontsize',23,'fontname','Times New Roman','FontWeight','normal');
        set(gca,'FontName','Times New Roman','FontSize',23,'FontWeight','normal');
        
        figure;plot(t,a2,'LineWidth',1.5);grid on;
        xlabel('\itt \rm[TU]','fontsize',23,'fontname','Times New Roman','FontWeight','normal');
        ylabel('\ita_s_e_c','fontsize',23,'fontname','Times New Roman','FontWeight','normal');
        set(gca,'FontName','Times New Roman','FontSize',23,'FontWeight','normal');
        figure;plot(t,X(:,end),'LineWidth',1.5);grid on;
        xlabel('\itt \rm[TU]','fontsize',23,'fontname','Times New Roman','FontWeight','normal');
        % ylabel('\ita_s_e_c','fontsize',23,'fontname','Times New Roman','FontWeight','normal');
        set(gca,'FontName','Times New Roman','FontSize',23,'FontWeight','normal');
        
        figure;
        plot3(X_secE(1,:),X_secE(2,:),X_secE(3,:), 'LineWidth', 1.5);
        hold on;
        plot3(0,0,0,'r.','MarkerSize',14); % moon;
        xlabel('\itx \rm[DU]','fontsize',18,'fontname','Times New Roman','FontWeight','normal');
        ylabel('\ity \rm[DU]','fontsize',18,'fontname','Times New Roman','FontWeight','normal');
        zlabel('\itz \rm[DU]','fontsize',18,'fontname','Times New Roman','FontWeight','normal');
        set(gca,'FontName','Times New Roman','FontSize',24,'FontWeight','normal');
        box on;grid on;axis equal;
        xlim([-0.25,0.25]); ylim([-0.25,0.25]);
        % n=length(t); x_IE = zeros(6,n);x_RSE=zeros(6,n);
        % for i=1:n
        %         x_IE(:,i)=trans_REM_IE(xx2(i,1:6)',t(i)-t(1));
        %         x_RSE(:,i)=trans_REM_RSE(xx2(i,1:6)',t(i));
        % end
        % axes(left2);
        % plot3(x_IE(1,:),x_IE(2,:),x_IE(3,:), 'LineWidth', 1.5);hold on; grid on;
        % plot3(0,0,0,'b.','MarkerSize',8);
        % plot3(cos(alpha),sin(alpha),zeros(length(alpha),1),'k--', 'LineWidth', 1.5);axis equal;
        % % plot3(xx1(:,1)+mu,xx1(:,2),zeros(length(alpha),1),'k-', 'LineWidth', 1.25);
        % title('Earth inertial');        hold off;
        % xlabel('\itx \rm[DU]','fontsize',18,'fontname','Times New Roman','FontWeight','normal');
        % ylabel('\ity \rm[DU]','fontsize',18,'fontname','Times New Roman','FontWeight','normal');
        % set(gca,'FontName','Times New Roman','FontSize',18,'FontWeight','normal');
        
    end
end

function [value,isterminal,direction] = only_hit_simple(t,X)
%%% 防撞击检测
%%% 一但航天器与地球月球距离小于地球半径或月球半径，就停止积分
% 分别检测：1 撞月、2 撞地、3 到达逃逸边界、4 近地点、5 远地点、6 近月点/7 离轨时刻

    u=0.0121506683;
    DU=3.84405*10^5;
    
    Re=[-u;0;0];                                           
    Rm=[1-u;0;0];
    Rs=[X(1);X(2);X(3)];
    rEM_hill = 1.501292597582450e+06*2/DU;

    R2=Rs-Rm;
    R1=Rs-Re;
    % 计算地心距和径向速度
    v1 = X(4:6)-[X(2);-X(1)-u;0]; % 惯性系
    v2 = X(4:6)-[X(2);-X(1)+1-u;0];
    radial_velocity_e = dot(R1,v1)/norm(R1); % 这就是 dr1/dt
    radial_velocity_m = dot(R2,v2)/norm(R2); % 这就是 dr2/dt
    % H_t = norm(cross(R2,v2));
    % [~,ind] = min(abs(mod(t-t_nomi(1),T)-t_nomi));
    % H_t_nomi = H_nomi(ind);
    % departure_index = abs(H_t-H_t_nomi)/abs(H_t_nomi);
    value = [norm(R2)-((1738/2)/DU);norm(R1)-((6378*2/3)/DU);norm(Rs)-rEM_hill;...
        radial_velocity_e;radial_velocity_e;radial_velocity_m;];
    isterminal = [1;1;0;0;0;0]; 
    direction = [0;0;1;1;-1;1];   % all direction
    % 1 -- 正向穿越; -1 -- 负向穿越; 

end
