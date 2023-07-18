function Rsquare = calculate_Rsquare(fitted, target)
    ss_total = sum((target - mean(target)).^2);
    ss_res = sum((target - fitted).^2);
    Rsquare = 1 - ss_res / ss_total;
end