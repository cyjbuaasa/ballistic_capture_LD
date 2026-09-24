function [tend_f,xend_f,eve_f,keep_event] = postprocess_xmoon_events(tend,xend,eve)
% Keep event 3/5 only when the xmoon capture gate is satisfied.

    u = 0.0121506683;
    DU = 3.84405*10^5;
    Rm = [1-u;0;0];

    r_event3 = 1.95*66180/DU;
    r_event5 = 1.25*66180/DU;

    keep_event = true(size(eve));
    for ii = 1:numel(eve)
        if eve(ii) == 3 || eve(ii) == 5
            Xi = xend(ii,:).';
            Rs = Xi(1:3);
            Esec = E_secondary(Xi, u);
            if eve(ii) == 3
                keep_event(ii) = norm(Rs-Rm) <= r_event3 && Esec <= 0;
            else
                keep_event(ii) = norm(Rs-Rm) <= r_event5 ;
            end
        end
    end

    tend_f = tend(keep_event,:);
    xend_f = xend(keep_event,:);
    eve_f = eve(keep_event,:);

end
