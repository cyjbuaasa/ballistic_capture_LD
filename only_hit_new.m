function [value,isterminal,direction] = only_hit_new(t,X,~,~,~,~,~,~)
%%% 防撞击检测
%%% 一但航天器与地球月球距离小于地球半径或月球半径，就停止积分
% 分别检测：1 撞地、2 撞月、3 近月点、4 到达逃逸边界、5 transit

    u=0.0121506683;
    DU=3.84405*10^5;
    
    Re=[-u;0;0];                                           
    Rm=[1-u;0;0];
    Rs=[X(1);X(2);X(3)];
    rEM_hill = 1.501292597582450e+06*2/DU;

    R2=Rs-Rm;
    R1=Rs-Re;
    % 
    v1 = X(4:6)-[X(2);-X(1)-u;0]; %
    v2 = X(4:6)-[X(2);-X(1)+1-u;0];
    radial_velocity_e = dot(R1,v1)/norm(R1); % ?dr1/dt
    radial_velocity_m = dot(R2,v2)/norm(R2); % ?dr2/dt
    % % H_t = norm(cross(R2,v2));
    % % [~,ind] = min(abs(mod(t-t_nomi(1),T)-t_nomi));
    % % H_t_nomi = H_nomi(ind);
    % % departure_index = abs(H_t-H_t_nomi)/abs(H_t_nomi);
    % M_t = 1/2*(norm(X(1:3))^2-norm(X0(1:3))^2);% X(8);
    % [~,ind] = min(abs(mod(t-t_nomi(1),T)-t_nomi));
    % M_t_nomi = Mom_nomi(ind);departure_index = abs(M_t-M_t_nomi)/(max(Mom_nomi)-min(Mom_nomi))*2;
    % Esec=E_secondary(X, u);
    % if norm(X(1:3)-Rm) <= 1.95*66180/DU && Esec <= 0
        value = [norm(R2)-((1738)/DU);norm(R1)-((6378)/DU);radial_velocity_m;...
            norm(R2)-(6*66180/DU);Rs(1)-1+u];% norm(R2)-(6*66180/DU)
        isterminal = [1;1;0;0;0]; 
        direction = [0;0;1;1;0];   % 
        % 3:perilune
    % else
    %     % if radial_velocity_m<0, rvmtmp=-1; 
    %     % else, rvmtmp=1; end
    %     % if Rs(1)-1+u<0, rtranstmp=-1; 
    %     % else, rtranstmp=1; end
    %     value = [norm(R2)-((1738)/DU);norm(R1)-((6378)/DU);nan;...
    %         norm(R2)-(6*66180/DU);nan]; % disabled perilune/xmoon events stay finite for odezero
    %     isterminal = [1;1;0;0;0]; 
    %     direction = [0;0;1;1;0];   % 
    %     % 3:perilune
    % end
    % value = [norm(R2)-((1738)/DU);norm(R1)-((6378)/DU);radial_velocity_m];
    % isterminal = [1;1;0]; 
    % direction = [0;0;1];   % all direction
    % % 3:perilune
    % 1 -- posi; -1 -- nega; 

end
