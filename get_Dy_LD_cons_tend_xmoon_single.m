function [LD_FTLE,LD_a2,Np_prl,KEarth,KMoon,KEscape,KXMoon,tend,xend,eve] = ...
    get_Dy_LD_cons_tend_xmoon_single(t0,t_end,y0,...
    dyn_str,LD_str,eps,exp_alp,keta,r_sec)
% Compute the LD and event classifications at one scalar terminal time.

% t_end is a scalar duration measured from t0. Constructing tspan through a
% fixed-size array keeps its size at 1-by-2 during MATLAB Coder generation.
assert(isscalar(t_end) && isfinite(t_end) && t_end > 0, ...
    't_end must be a positive finite scalar.');

options_LD = odeset('RelTol',1e-9,'AbsTol',1e-9,'Events',@only_hit_new);

tspan = zeros(1,2);
tspan(1) = t0;
tspan(2) = t0 + t_end;
[t,Xi,tend,xend,eve] = ode45(@LD_syn_rot,tspan,[y0;0;0],...
    options_LD,dyn_str,LD_str,eps,exp_alp,keta,r_sec);
[tend,xend,eve] = postprocess_xmoon_events(tend,xend,eve);

% If integration ends at an event before tspan(2), use the final available
% state, matching the fallback behavior of the multi-time implementation.
t_final = tspan(2);
id_t = find(abs(t-t_final) < 1e-8,1,'last');
if isempty(id_t)
    id_t = size(Xi,1);
end
LD_FTLE = Xi(id_t,7);
LD_a2 = Xi(id_t,8);

KEarth = 0;
KMoon = 0;
KEscape = 0;
KXMoon = 0;
Np_prl = 0;

if isempty(tend)
    return
end

if ~isempty(find(eve == 1 & tend <= t_final,1))
    KEarth = 1;
end
if ~isempty(find(eve == 2 & tend <= t_final,1))
    KMoon = 1;
end
if ~isempty(find(eve == 5 & tend <= t_final,1))
    KXMoon = 1;
end

t_escape = tend(eve == 4);
t_prl = tend(eve == 3);
escape_first = inf;
escape_before_prl = false;
if ~isempty(t_escape)
    escape_first = min(t_escape);
    if isempty(t_prl) || isempty(find(t_prl < escape_first,1))
        escape_before_prl = true;
    end
end

if escape_before_prl && escape_first <= t_final
    KEscape = 1;
end

if KXMoon == 1
    if isfinite(escape_first)
        Np_prl = length(find(eve == 3 & tend <= t_final & ...
            tend < escape_first));
    else
        Np_prl = length(find(eve == 3 & tend <= t_final));
    end
end

end
