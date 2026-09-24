function [orb_dynam_index,Np_prl,eve_struct,tend,xend,eve] = get_Dyindex_1_3(t0,t_end,y0,nonlinear_id,m_STT,...
    dyn_str,eps,alp,keta,dx_eps,LD_str,mu)
% 应该加上计算环绕次星圈数的角度量    

    % nonlinear_id = 3.1;
    % 1--FTLE; 2--HOFTLE(.1-STT; .2-DSTT; .3-TDSTT);
    % 3--LD(.1-E2; .2-H2; .3-v2; .4-a2; .5-H2|E2|r2|v2|onlyv(v-norm) & a2)
    % m_STT = 2; % m order STT
    % dyn_str = 'EM'; % 'EM' 'diff_sbcm'
    options_LD = odeset('RelTol',1e-9,'AbsTol',1e-9,'Events',@only_hit_LD);
    DU = 384405;
    rm_hill = 66180/DU;
    init_cond = [];
    if nonlinear_id == 1
        % FTLE cost-expensive
        [t,X,tend,xend,eve] = ode45(@dym,[0,t_end]+t0,[y0;reshape(eye(6),36,1)],options_LD,mu);
        Phi =reshape(X(end,7:42),6,6); CGT_end = Phi'*Phi;
        [eig_vec,lamda]=eig(CGT_end);
        LD_FTLE_mesh = 1/t_end*log(sqrt(max(lamda(:))));
        % A = tenzeros([6 6]);
        % for ii = 1:6
        %     for jj = 1:6
        %         A(perms([ii jj])) = CGT_end(ii,jj);
        %     end
        % end
        % rng('default');
        % [lambda_max,eig_vec_max, flag, it] = eig_sshopm(A, 'MaxIts', 100, 'Shift', 1,'Display',0);
        % LD_FTLE_mesh(i) = 1/t_end*log(sqrt(lambda_max));
    elseif nonlinear_id == 2.3
        [t,X,tend,xend,eve] = ode45(@dym,[0,t_end],[y0;reshape(eye(6),36,1)],options_LD,mu);
        lamda_TV = zeros(length(t),1); 
        for ii = 1:length(t)
            Phi =reshape(X(ii,7:42),6,6); CGTi = Phi'*Phi;
                [eig_vec,lamda]=eig(CGTi);
            lamda_TV(ii) = max(abs(diag(lamda)));
        end
        % FTLE_TV = 1./(t).*log(sqrt(lamda_TV));
        % LD_FTLE_mesh(i) = max(FTLE_TV);
        [max_lamda_TV,j]=max(lamda_TV);
        LD_FTLE_mesh = 1/(t(j))*log(sqrt(max_lamda_TV));
    elseif nonlinear_id == 3.1
        % LD cost-cheaper
        % velocity-based
        % [t,X,tend,xend,eve] = ode45(@LD_v_rot,[0,t_end]+theta_s/ome_sun,[y0;0],options_LD,eps,alp);
        % energy-based
        [t,Xi,tend,xend,eve] = ode45(@LD_E2_rot,[0,t_end]+t0,[y0;0],options_LD,dyn_str,eps,alp);
        LD_FTLE_mesh = Xi(end,7);
        % id = find(eve == 3);
        % if isempty(id)
        %     LD_FTLE_mesh = LD_FTLE_mesh/2;
        % end
    elseif nonlinear_id == 3.2
        % orb-angular-momentum-based
        [t,Xi,tend,xend,eve] = ode45(@LD_H2_rot,[0,t_end]+t0,[y0;0],options_LD,dyn_str,eps,alp);
        LD_FTLE_mesh = Xi(end,7);
        % fix: for perilune within Lunar Hill's Sphere
        % id = find(eve == 3);
    elseif nonlinear_id == 3.3
        % velocity-based
        % this method is not fixed
        [t,Xi,tend,xend,eve] = ode45(@LD_v2_rot,[0,t_end]+t0,[y0;0],options_LD,dyn_str,eps,alp);
        LD_FTLE_mesh = Xi(end,7);
    elseif nonlinear_id == 3.4
        % sec-accelaration-based
        % this method is not fixed
        [t,Xi,tend,xend,eve] = ode45(@LD_a2_rot,[0,t_end]+t0,[y0;0],options_LD,dyn_str,keta);
        LD_FTLE_mesh = Xi(end,7);
    elseif nonlinear_id == 3.5
        % H2|E2 & sec-accelaration-based
        % this method is not fixed
        [t,Xi,tend,xend,eve] = ode45(@LD_syn_rot,[0,t_end]+t0,[y0;0;0],options_LD,dyn_str,LD_str,eps,alp,keta);
        LD_FTLE_mesh = Xi(end,7:8);
    else
        t=0;X = zeros(6,1);
        LD_FTLE_mesh=0;
        tend=0;xend=zeros(6,1);eve = 0;
    end
    orb_dynam_index = LD_FTLE_mesh;
    eve_struct = struct('tend',tend,'xend',xend,'eve',eve);
    Np_prl = length(find(eve == 3));
end




