function rmspe = calculate_rmspe(x, y)
    rmspe = sqrt(mean(((x - y)./ y).^2)) * 100;
end