function [LD_FTLE_i,LD_a2_i,Np_prl_i,KEarth_i,KMoon_i,KEscape_i,KXMoon_i,tend,xend,eve] = get_Dy_LD_cons_tend_xmoon(t0,t_end,num_tend,y0,...
    dyn_str,LD_str,eps,exp_alp,keta,r_sec)

% only_hit_new records raw events; event 3/5 are filtered after ode45.
% The filtered tend/xend/eve are used below.
options_LD = odeset('RelTol',1e-9,'AbsTol',1e-9,'Events',@only_hit_new);
LD_FTLE_i = NaN(1,num_tend);
LD_a2_i = NaN(1,num_tend);
Np_prl_i = NaN(1,num_tend);
KEarth_i = NaN(1,num_tend);
KMoon_i = NaN(1,num_tend);
KEscape_i = NaN(1,num_tend);
KXMoon_i = NaN(1,num_tend);

[t,Xi,tend,xend,eve] = ode45(@LD_syn_rot,[t0,t0+t_end],[y0;0;0],...
    options_LD,dyn_str,LD_str,eps,exp_alp,keta,r_sec);
[tend,xend,eve] = postprocess_xmoon_events(tend,xend,eve);

escape_first = inf;
escape_before_prl = false;
if ~isempty(tend)
    t_escape = tend(eve == 4);
    t_prl = tend(eve == 3);
    if ~isempty(t_escape)
        escape_first = min(t_escape);
        if isempty(t_prl) || isempty(find(t_prl < escape_first,1))
            escape_before_prl = true;
        end
    end
end

for kt = 1:num_tend
    t_now = t0+t_end(kt);
    id_t = find(abs(t-(t_now)) < 1e-8,1,'last');
    if ~isempty(id_t)
        LD_FTLE_i(kt) = Xi(id_t,7);
        LD_a2_i(kt) = Xi(id_t,8);
    else
        LD_FTLE_i(kt) = Xi(end,7);
        LD_a2_i(kt) = Xi(end,8);
    end

    if ~isempty(tend)
        if ~isempty(find(eve == 1 & tend <= t_now, 1))
            KEarth_i(kt) = 1;
        end
        if ~isempty(find(eve == 2 & tend <= t_now, 1))
            KMoon_i(kt) = 1;
        end
        if ~isempty(find(eve == 5 & tend <= t_now, 1))
            KXMoon_i(kt) = 1;
        end
        if escape_before_prl && escape_first <= t_now
            KEscape_i(kt) = 1;
        end

        if KXMoon_i(kt) == 1
            if isfinite(escape_first)
                Np_prl_i(kt) = length(find(eve == 3 & tend <= t_now & tend < escape_first));
            else
                Np_prl_i(kt) = length(find(eve == 3 & tend <= t_now));
            end
        else
            Np_prl_i(kt) = 0;
        end
    else
        Np_prl_i(kt) = 0;
    end
end

end
