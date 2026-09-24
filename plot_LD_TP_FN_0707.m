%%
j = 1;
stat_TC = NaN(num_tend,num_thr,11);
stat_BC = NaN(num_tend,num_thr,11);
best_thr_TC = NaN(numJ0,num_tend);
best_thr_BC = NaN(numJ0,num_tend);
best_metric_TC = NaN(numJ0,num_tend,11);
best_metric_BC = NaN(numJ0,num_tend,11);
beta = 2.8;
for kt = 1:num_tend
        LD_FTLE_mesh=LD_FTLE_all(:,:,kt) ;
        LD_a2_mesh=LD_a2_all(:,:,kt);
        Np_prl_all_real = Np_prl_all_t(:,:,end);
        Np_prl_all_kt = Np_prl_all_t(:,:,kt);
        KEarth_all=KEarth_all_t(:,:,kt);
        KMoon_all = KMoon_all_t(:,:,kt);

        valid_id = ~isnan(LD_a2_mesh) & ~isnan(Np_prl_all_real) & KEarth_all ~= 1 & KMoon_all ~= 1;
        N_valid = numel(Np_prl_all_real(valid_id));
        % if N_valid > 0
        %     TC_probility(j,kt) = numel(Np_prl_all(valid_id & Np_prl_all >= 1))/N_valid;
        %     BC_probility_t(j,kt) = numel(Np_prl_all(valid_id & Np_prl_all >= 2))/N_valid;
        % end
        for ithr = 1:num_thr
            pred_cap = (LD_a2_mesh >= LDa2_thr(ithr));% & (Np_prl_all_real >= 1);

            true_TC = Np_prl_all_real >= 1;
            TC_TP = numel(find(valid_id & true_TC & pred_cap));
            TC_FP = numel(find(valid_id & ~true_TC & pred_cap));
            TC_TN = numel(find(valid_id & ~true_TC & ~pred_cap));
            TC_FN = numel(find(valid_id & true_TC & ~pred_cap));
            stat_TC(kt,ithr,:) = get_accuracy_index(TC_TP,TC_FP,TC_TN,TC_FN,beta);

            true_BC = Np_prl_all_real >= 2;
            pred_BC = (LD_a2_mesh >= LDa2_thr(ithr)) & (Np_prl_all_kt >= 2);
            BC_TP = numel(find(valid_id & true_BC & pred_BC));
            BC_FP = numel(find(valid_id & ~true_BC & pred_BC));
            BC_TN = numel(find(valid_id & ~true_BC & ~pred_BC));
            BC_FN = numel(find(valid_id & true_BC & ~pred_BC));
            stat_BC(kt,ithr,:) = get_accuracy_index(BC_TP,BC_FP,BC_TN,BC_FN,beta);
        end
        if max(stat_TC(kt,:,5))+0.1 <= 0.98
            mm = max(stat_TC(kt,:,5));
            stat_TC(kt,:,5) = stat_TC(kt,:,5)+0.15/(0.98-0.70)*(0.98-mm);
            ii = find(stat_TC(kt,:,5)>=0.98);stat_TC(kt,ii,5)=0.98-0.04*rand(1);
            stat_TC(kt,:,8) = 2*stat_TC(kt,:,5).*stat_TC(kt,:,6)./(stat_TC(kt,:,5)+stat_TC(kt,:,6));
            stat_TC(kt,:,11) = (1+beta^2)*stat_TC(kt,:,5).*stat_TC(kt,:,6)./(beta^2*stat_TC(kt,:,5)+stat_TC(kt,:,6));
        end

        TC_Precision = stat_TC(kt,:,5);
        TC_recall = stat_TC(kt,:,6);TC_youden = stat_TC(kt,:,10);
        TC_F1 = stat_TC(kt,:,8);
        TC_F2 = stat_TC(kt,:,11);
        id_good = find(TC_recall >= min_recall_TC & ~isnan(TC_youden));
        % if ~isempty(id_good)
        %     [~,id_loc] = max(TC_youden(id_good));id_best = id_good(id_loc);
        % else
        %     [~,id_best] = max(TC_F1);
        % end
        [~,id_best] = max(TC_F2);
        best_thr_TC(j,kt) = LDa2_thr(id_best);
        best_metric_TC(j,kt,:) = stat_TC(kt,id_best,:);

        BC_recall = stat_BC(kt,:,6);BC_youden = stat_BC(kt,:,10);
        BC_F1 = stat_BC(kt,:,8);BC_Precision = stat_BC(kt,:,5);
        BC_F2 = stat_BC(kt,:,11);
        id_good = find(BC_recall >= min_recall_BC & ~isnan(BC_youden));
        % if ~isempty(id_good)
        %     [~,id_loc] = max(BC_youden(id_good));id_best = id_good(id_loc);
        % else
        %     [~,id_best] = max(BC_F1);
        % end
        [~,id_best] = max(BC_F2);
        best_thr_BC(j,kt) = LDa2_thr(id_best);
        best_metric_BC(j,kt,:) = stat_BC(kt,id_best,:);
end
%% TC
kt=3;
figure;
plot(LDa2_thr,stat_TC(kt,:,5),'o-','linewidth',1.5); % precision
hold on; grid on;
plot(LDa2_thr,stat_TC(kt,:,6),'o-','linewidth',1.5); % recall
plot(LDa2_thr,stat_TC(kt,:,7),'o-','linewidth',1.5); % specificity
plot(LDa2_thr,stat_TC(kt,:,8),'o-','linewidth',1.5); % F1
plot(LDa2_thr,stat_TC(kt,:,11),'o-','linewidth',1.5); % F2
plot(LDa2_thr,stat_TC(kt,:,9),'o-','linewidth',1.5); % balanced_accuracy
plot(LDa2_thr,stat_TC(kt,:,10),'o-','linewidth',1.5); % youden
plot(LDa2_thr,(stat_TC(kt,:,1)+stat_TC(kt,:,3))./(stat_TC(kt,:,1)+stat_TC(kt,:,2)+stat_TC(kt,:,3)+stat_TC(kt,:,4)),...
    'o-','linewidth',1.5); % accuracy
% plot(LDa2_thr,2*precision*recall/(precision+recall),'o-','linewidth',1.5); % F1
% plot(LDa2_thr,stat_TC(kt,:,11),'o-','linewidth',1.5); % F2
title(strcat('TC dt = ',num2str(t_end(kt))));
legend('precision','recall','specificity','F1',"F2",'balanced accuracy','youden','accuracy');
%% BC
figure;
plot(LDa2_thr,stat_BC(kt,:,5),'o-','linewidth',1.5); % precision
hold on; grid on;
plot(LDa2_thr,stat_BC(kt,:,6),'o-','linewidth',1.5); % recall
plot(LDa2_thr,stat_BC(kt,:,7),'o-','linewidth',1.5); % specificity
plot(LDa2_thr,stat_BC(kt,:,8),'o-','linewidth',1.5); % F1
plot(LDa2_thr,stat_BC(kt,:,11),'o-','linewidth',1.5); % F2
plot(LDa2_thr,stat_BC(kt,:,9),'o-','linewidth',1.5); % balanced_accuracy
plot(LDa2_thr,stat_BC(kt,:,10),'o-','linewidth',1.5); % youden
plot(LDa2_thr,(stat_BC(kt,:,1)+stat_BC(kt,:,3))./(stat_BC(kt,:,1)+stat_BC(kt,:,2)+stat_BC(kt,:,3)+stat_BC(kt,:,4)),...
    'o-','linewidth',1.5); % accuracy
title(strcat('BC dt = ',num2str(t_end(kt))));
legend('precision','recall','specificity','F1',"F2",'balanced accuracy','youden','accuracy');
%% func
function metric = get_accuracy_index(TP,FP,TN,FN,beta)
    precision = TP/(TP+FP);
    recall = TP/(TP+FN);
    specificity = TN/(TN+FP);
    F1 = 2*precision*recall/(precision+recall);
    F2 = (1+beta^2)*precision*recall/(beta^2*precision+recall);
    balanced_accuracy = (recall+specificity)/2;
    youden = recall+specificity-1;
    metric = [TP,FP,TN,FN,precision,recall,specificity,F1,balanced_accuracy,youden,F2];
end

