%% pooled estimate and bootstrap CI demo
% This demo shows how to compute pooled recall, precision and F1 over
% spatial slices, e.g. different (gamma, psi2) conditions.
%
% For each fixed J0:
%   1. Each (gamma, psi2) slice has its own (alpha, psi1) map.
%   2. TP/FP/TN/FN are computed inside each slice.
%   3. The slice confusion matrices are summed to obtain the pooled estimate.
%   4. Bootstrap resampling is performed over slices to obtain 95% CIs.

clear; clc;

%% Example settings
J0_all = 2.95:0.01:3.20;
numJ = length(J0_all);
B = 1000;              % bootstrap replicates
alpha_CI = 0.05;       % 95% CI

recall_pool = NaN(numJ,1);
precision_pool = NaN(numJ,1);
F1_pool = NaN(numJ,1);
recall_CI = NaN(numJ,2);
precision_CI = NaN(numJ,2);
F1_CI = NaN(numJ,2);

for j = 1:numJ
    J0 = J0_all(j);

    % -------------------------------------------------------------
    % Replace this block by your own data loading loop.
    %
    % true_cell{s}: logical alpha x psi1 map for true TC/BC
    % pred_cell{s}: logical alpha x psi1 map for LD-predicted TC/BC
    % valid_cell{s}: logical alpha x psi1 map for valid samples
    %
    % Here s is the index of a spatial slice, i.e. one (gamma, psi2).
    % Different slices may have different alpha x psi1 sizes.
    % -------------------------------------------------------------
    [true_cell,pred_cell,valid_cell] = load_one_J0_demo(J0);

    [metric,CI,conf_pool,conf_slice] = pooled_metric_bootstrap( ...
        true_cell,pred_cell,valid_cell,B,alpha_CI);

    recall_pool(j) = metric.recall;
    precision_pool(j) = metric.precision;
    F1_pool(j) = metric.F1;
    recall_CI(j,:) = CI.recall;
    precision_CI(j,:) = CI.precision;
    F1_CI(j,:) = CI.F1;

    fprintf('J0 = %.3f, TP=%d, FP=%d, TN=%d, FN=%d, recall=%.4f, precision=%.4f, F1=%.4f\n', ...
        J0,conf_pool.TP,conf_pool.FP,conf_pool.TN,conf_pool.FN, ...
        metric.recall,metric.precision,metric.F1);
end

%% Plot with asymmetric error bars
figure('Color','w'); hold on; box on; grid on;
errorbar(J0_all,recall_pool,recall_pool-recall_CI(:,1),recall_CI(:,2)-recall_pool,'o-','LineWidth',1.5);
errorbar(J0_all,precision_pool,precision_pool-precision_CI(:,1),precision_CI(:,2)-precision_pool,'o-','LineWidth',1.5);
errorbar(J0_all,F1_pool,F1_pool-F1_CI(:,1),F1_CI(:,2)-F1_pool,'o-','LineWidth',1.5);
xlabel('J0');
ylabel('metric');
legend('recall','precision','F1','Location','best');

%% Local functions
function [metric,CI,conf_pool,conf_slice] = pooled_metric_bootstrap(true_cell,pred_cell,valid_cell,B,alpha)
    num_slice = numel(true_cell);
    conf_slice.TP = zeros(num_slice,1);
    conf_slice.FP = zeros(num_slice,1);
    conf_slice.TN = zeros(num_slice,1);
    conf_slice.FN = zeros(num_slice,1);

    for s = 1:num_slice
        true_id = logical(true_cell{s});
        pred_id = logical(pred_cell{s});
        valid_id = logical(valid_cell{s});

        conf_slice.TP(s) = sum(valid_id &  true_id &  pred_id,'all');
        conf_slice.FP(s) = sum(valid_id & ~true_id &  pred_id,'all');
        conf_slice.TN(s) = sum(valid_id & ~true_id & ~pred_id,'all');
        conf_slice.FN(s) = sum(valid_id &  true_id & ~pred_id,'all');
    end

    conf_pool.TP = sum(conf_slice.TP);
    conf_pool.FP = sum(conf_slice.FP);
    conf_pool.TN = sum(conf_slice.TN);
    conf_pool.FN = sum(conf_slice.FN);
    metric = metric_from_conf(conf_pool.TP,conf_pool.FP,conf_pool.TN,conf_pool.FN);

    boot_recall = NaN(B,1);
    boot_precision = NaN(B,1);
    boot_F1 = NaN(B,1);

    for b = 1:B
        id = randi(num_slice,num_slice,1);
        TP = sum(conf_slice.TP(id));
        FP = sum(conf_slice.FP(id));
        TN = sum(conf_slice.TN(id));
        FN = sum(conf_slice.FN(id));
        m = metric_from_conf(TP,FP,TN,FN);
        boot_recall(b) = m.recall;
        boot_precision(b) = m.precision;
        boot_F1(b) = m.F1;
    end

    CI.recall = quantile(boot_recall,[alpha/2,1-alpha/2]);
    CI.precision = quantile(boot_precision,[alpha/2,1-alpha/2]);
    CI.F1 = quantile(boot_F1,[alpha/2,1-alpha/2]);
end

function metric = metric_from_conf(TP,FP,TN,FN)
    metric.TP = TP;
    metric.FP = FP;
    metric.TN = TN;
    metric.FN = FN;
    metric.recall = safe_div(TP,TP+FN);
    metric.precision = safe_div(TP,TP+FP);
    metric.F1 = safe_div(2*metric.precision*metric.recall,metric.precision+metric.recall);
end

function y = safe_div(a,b)
    if b == 0
        y = NaN;
    else
        y = a/b;
    end
end

function [true_cell,pred_cell,valid_cell] = load_one_J0_demo(~)
    % Demo data only. Replace this function with loading your real
    % (gamma, psi2) slices for one fixed J0.
    num_gamma = 5;
    num_psi2 = 6;
    num_slice = num_gamma*num_psi2;
    true_cell = cell(num_slice,1);
    pred_cell = cell(num_slice,1);
    valid_cell = cell(num_slice,1);

    for s = 1:num_slice
        n_alpha = randi([250,360]);
        n_psi1 = randi([240,360]);
        valid = rand(n_alpha,n_psi1) > 0.03;
        true_id = rand(n_alpha,n_psi1) > 0.80;
        miss = rand(n_alpha,n_psi1) < 0.05;
        false_alarm = rand(n_alpha,n_psi1) < 0.02;
        pred_id = (true_id & ~miss) | (~true_id & false_alarm);

        true_cell{s} = true_id;
        pred_cell{s} = pred_id;
        valid_cell{s} = valid;
    end
end
