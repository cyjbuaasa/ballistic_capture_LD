function [orb_dynam_index,Np_prl,eve_struct,t,X,tend,xend,eve] = get_Dyindex(t0,t_end,y0,nonlinear_id,m_STT,...
    dyn_str,eps,alp,keta,dx_eps,LD_str,options_LD,mu)
    
    % nonlinear_id = 3.1;
    % % 1--FTLE; 2--HOFTLE(.1-STT; .2-DSTT; .3-TDSTT); 3--LD(.1-E2; .2-H2; .3-v2; .4-a2)
    % m_STT = 2; % m order STT
    % dyn_str = 'EM'; % 'EM' 'diff_sbcm'
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
    elseif nonlinear_id == 2.1
        % HOFTLE -- time increase, accuracy increase
        %   t_span: 积分时间（与输入一致）
        %   phi1: 1阶STT（STM，6×6）
        %   phi2: 2阶STT（6×6×6）
        %   phi3: 3阶STT（6×6×6×6）
        %   phi4: 4阶STT（6×6×6×6×6）
        % initialization: STM + highorser STT(m>=2)
        n = 6; % dimension-state
        phi1_0 = eye(n);
        phi2_0 = zeros(n, n, n);
        phi3_0 = zeros(n, n, n, n);
        phi4_0 = zeros(n, n, n, n, n);
        if m_STT == 1
            init_cond = [y0;reshape(phi1_0, [], 1)];
        elseif m_STT ==2
            init_cond = [y0;reshape(phi1_0, [], 1);reshape(phi2_0, [], 1)];
        elseif m_STT == 3
            init_cond = [y0;reshape(phi1_0, [], 1);reshape(phi2_0, [], 1);reshape(phi3_0, [], 1)];
        elseif m_STT == 4
        init_cond = [y0;reshape(phi1_0, [], 1);reshape(phi2_0, [], 1); 
                     reshape(phi3_0, [], 1);reshape(phi4_0, [], 1)];
        end

        [t,X,tend,xend,eve] = ode45(@crtbp_variational_equations,[0,t_end]+t0,...
            init_cond,options_LD,mu,m_STT);
        [xf, HOSTT_cell] = crt_hostt_perms(y0,X,mu,m_STT);
        HOCGT_cell = cell(4,1);
        [C2, C3, C4, C5] = HOCGT_perms(y0,HOSTT_cell,mu,m_STT);
        HOCGT_cell{1} = C2;HOCGT_cell{2} =C3;
        HOCGT_cell{3} = C4;HOCGT_cell{4} =C5;

        rng('default');
        lambda_max = zeros(m_STT,1);eig_vec_max=zeros(6,m_STT);
        sum_lamda_max = 0;
        for m_id = 1:m_STT
            [lambda_max(m_id),eig_vec_max(:,m_id), flag, it] = eig_sshopm(HOCGT_cell{m_id},...
                'MaxIts', 100, 'Shift', 1,'Display',0);
            sum_lamda_max = sum_lamda_max+abs(lambda_max(m_id))*dx_eps^(m_id+1-2);
        end
        sum_lamda_max = sqrt(sum_lamda_max);
        LD_FTLE_mesh = 1/t_end*log(sum_lamda_max);
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
        [t,X,tend,xend,eve] = ode45(@LD_E2_rot,[0,t_end]+t0,[y0;0],options_LD,dyn_str,eps,alp);
        LD_FTLE_mesh = X(end,end);
        % id = find(eve == 3);
        % if isempty(id)
        %     LD_FTLE_mesh = LD_FTLE_mesh/2;
        % end
        % if ~isempty(id)
        %     len = length(id); r2n = zeros(len,1);
        %     for k=1:len
        %         % attention the vector!
        %         r2n(k) = norm(xend(k,1:3)-[1-mu,0,0]);
        %     end
        %     LD_FTLE_mesh = LD_FTLE_mesh*rm_hill/min(r2n);
        % end
    elseif nonlinear_id == 3.2
        % orb-angular-momentum-based
        [t,X,tend,xend,eve] = ode45(@LD_H2_rot,[0,t_end]+t0,[y0;0],options_LD,dyn_str,eps,alp);
        LD_FTLE_mesh = X(end,end);
        % fix: for perilune within Lunar Hill's Sphere
        % id = find(eve == 3);
        % if ~isempty(id)
        %     len = length(id); r2n = zeros(len,1);
        %     for k=1:len
        %         r2n(k) = norm(xend(k,1:3)-[1-mu,0,0]);
        %     end
        %     LD_FTLE_mesh = LD_FTLE_mesh*rm_hill/min(r2n);
        % end
    elseif nonlinear_id == 3.3
        % velocity-based
        % this method is not fixed
        [t,X,tend,xend,eve] = ode45(@LD_v2_rot,[0,t_end]+t0,[y0;0],options_LD,dyn_str,eps,alp);
        LD_FTLE_mesh = X(end,end);
    elseif nonlinear_id == 3.4
        % sec-accelaration-based
        % this method is not fixed
        [t,X,tend,xend,eve] = ode45(@LD_a2_rot,[0,t_end]+t0,[y0;0],options_LD,dyn_str,keta);
        LD_FTLE_mesh = X(end,7);
    elseif nonlinear_id == 3.5
        % H2|E2 & sec-accelaration-based
        % this method is not fixed
        [t,X,tend,xend,eve] = ode45(@LD_syn_rot,[0,t_end]+t0,[y0;0;0],options_LD,dyn_str,LD_str,eps,alp,keta);
        LD_FTLE_mesh = X(end,7:8);
    else
        t=0;X = zeros(6,1);
        LD_FTLE_mesh=0;
        tend=0;xend=zeros(6,1);eve = 0;
    end
    orb_dynam_index = LD_FTLE_mesh;
    eve_struct = struct('tend',tend,'xend',xend,'eve',eve);
    Np_prl = length(find(eve == 3));
end




