function [score, precision, recall] = avg_f1(Thetas, invcov_series)
[n, ~, T] = size(Thetas);

total_tp = 0;
total_fp = 0;
total_fn = 0;
for i=1:T
    preds = Thetas(:,:,i);
    actuals = invcov_series(:,:,i);
    pred_true = (preds ~= 0);
    actual_true = (actuals ~= 0);
    
    true_positive = sum(sum(pred_true .* actual_true));
    false_positive = sum(sum(pred_true .* ~actual_true));
    false_negative = sum(sum(~pred_true .* actual_true));
    
    total_tp = total_tp + true_positive;
    total_fp = total_fp + false_positive;
    total_fn = total_fn + false_negative;
end

precision = total_tp / (total_tp + total_fp);
recall = total_tp / (total_tp + total_fn);

score = 2 * precision * recall / (precision + recall);
end