clc;

%% plot epsilon
min_wl = 500;
max_wl = 500;
min_used_wl = 0;
max_used_wl = 2000;
display_name = {'HbO2', 'Hb', 'Bilirubin', 'Biliverdin', 'Methemoglobin'};
plot_arr = [];
for chromophore_index = 1:length(epsilon_raw)
    plot_arr(chromophore_index) = plot(epsilon_raw{chromophore_index}(:, 1), epsilon_raw{chromophore_index}(:, 2), ...
                                        'LineWidth', 2);
    hold on
    
    % judge x axis range
    if min_wl > min(epsilon_raw{chromophore_index}(:, 1))
        min_wl = min(epsilon_raw{chromophore_index}(:, 1));
    end
    if max_wl < max(epsilon_raw{chromophore_index}(:, 1))
        max_wl = max(epsilon_raw{chromophore_index}(:, 1));
    end
    
    % judge used wl range (overlapped wl)
    if min_used_wl < min(epsilon_raw{chromophore_index}(:, 1))
        min_used_wl = min(epsilon_raw{chromophore_index}(:, 1));
    end
    if max_used_wl > max(epsilon_raw{chromophore_index}(:, 1))
        max_used_wl = max(epsilon_raw{chromophore_index}(:, 1));
    end    
end
title('Chromophores\_usedRange');
xlabel('wavelength (nm)');
ylabel('molar extinction (cm-1/M)');
% xlim([min_wl, max_wl]);
xlim([min_used_wl, max_used_wl]);
xline(min_used_wl, '--k', 'LineWidth', 2);
xline(max_used_wl, '--k', 'LineWidth', 2);
legend(plot_arr, display_name)
% saveas(gcf,'Methemoglobin_epsilon.png')

%% plot liveData
fig_title = ['DRS spectrum - ', subject{sheet}];
t = tiledlayout(2, week_num);
for index_week = 1:week_num
    for index_category=1:length(liveData_categories)
        nexttile(week_num + index_week - (index_category-1)*week_num)
        plot(liveData_seg.(subject{sheet}).(['week_', num2str(index_week)]).(liveData_categories{index_category})(:, 1), ...
             liveData_seg.(subject{sheet}).(['week_', num2str(index_week)]).(liveData_categories{index_category})(:, 2) ...
         );
        title(['week\_', num2str(index_week), ' ', liveData_categories{index_category}])
    end
end
ylabel(t, 'intensity')
xlabel(t, 'wavelength (nm)')
title(t,fig_title)
% saveas(gcf,fullfile(fitting_folder, [fig_title, '.png']));

%% plot used liveData
fig_title = ['DRS spectrum - ', subject{sheet}, ' (used wl points, ', num2str(wl_res), 'nm)'];
t = tiledlayout(2, week_num);
for index_week = 1:week_num
    for index_category=1:length(liveData_categories)
        nexttile(week_num + index_week - (index_category-1)*week_num)
        % initial graph
        plot(liveData_seg.(subject{sheet}).(['week_', num2str(index_week)]).(liveData_categories{index_category})(:, 1), ...
             liveData_seg.(subject{sheet}).(['week_', num2str(index_week)]).(liveData_categories{index_category})(:, 2) ...
         );
        hold on
        % mark sampling point
        plot(liveData_used.(subject{sheet}).(['week_', num2str(index_week)]).(liveData_categories{index_category})(:, 1), ...
             liveData_used.(subject{sheet}).(['week_', num2str(index_week)]).(liveData_categories{index_category})(:, 2), ...
             'ro');
        title(['week\_', num2str(index_week), ' ', liveData_categories{index_category}])
    end
end
ylabel(t, 'intensity')
xlabel(t, 'wavelength (nm)')
title(t, fig_title)
% saveas(gcf,[fig_title, '.png'])

%% plot each week delta_c ---- no used
xtick_names = {'bilirubin', 'biliverdin', 'methemoglobin'};
color_set = {'#FFD2D2', '#FD8A8A', '#FD0000', '#747474', '#000000'};

for index = 1:week_num
    plot(liveData_used.(subject{sheet}).(['week_', num2str(index)]).fitted_delta_c, '-o', 'Color', color_set{index}, 'LineWidth', 2);
    hold on
end
legend('week_1', 'week_2', 'week_3', 'week_4', 'week_5');
set(gca, 'xtick', (1:length(epsilon_files)), 'xticklabel', xtick_names);
xlabel('Chromophore species')
ylabel('mol/L (M)')
title('Variation of chromophore concentration (each week)')
grid on

%% plot each week delta_c ----------- new
sheet = 3; week_num = 4;
for index = 1:week_num
    liveData_used.(subject{sheet}).HbO2(index) = liveData_used.(subject{sheet}).(['week_', num2str(index)]).fitted_delta_c(1);
    liveData_used.(subject{sheet}).Hb(index) = liveData_used.(subject{sheet}).(['week_', num2str(index)]).fitted_delta_c(2);
    liveData_used.(subject{sheet}).Bilirubin(index) = liveData_used.(subject{sheet}).(['week_', num2str(index)]).fitted_delta_c(3);
    liveData_used.(subject{sheet}).Biliverdin(index) = liveData_used.(subject{sheet}).(['week_', num2str(index)]).fitted_delta_c(4);
    liveData_used.(subject{sheet}).Methemoglobin(index) = liveData_used.(subject{sheet}).(['week_', num2str(index)]).fitted_delta_c(5);
end

display_name = {'HbO2', 'Hb', 'Bilirubin', 'Biliverdin', 'Methemoglobin'};
for idx = 1:length(display_name)
    plot(1:week_num, liveData_used.(subject{sheet}).(display_name{idx}), '-o', 'LineWidth', 2)
    hold on
end

legend(display_name);
set(gca, 'xtick', (1:week_num), 'xticklabel', {'week\_1', 'week\_2', 'week\_3', 'week\_4', 'week\_5'});
xlabel('Time (week)')
ylabel('mol/L (M)')
xlim([0.5, week_num+0.5]);
title([subject{sheet}, ' - Variation of chromophore concentration (each week) - ', num2str(wl_res), 'nm'])
grid on

%% plot delta_OD v.s. fitted_delta_OD
sheet = 3; week_num = 4;
for week_index = 1:week_num
    % plot delta_OD
    plot(liveData_used.(subject{sheet}).(['week_', num2str(week_index)]).delta_OD(:, 1), ...
        liveData_used.(subject{sheet}).(['week_', num2str(week_index)]).delta_OD(:, 2), ...
        '-o', 'LineWidth', 1)
    
    hold on
    
    % plot fitted delta_OD
    plot(liveData_used.(subject{sheet}).(['week_', num2str(week_index)]).fitted_delta_OD(:, 1), ...
        liveData_used.(subject{sheet}).(['week_', num2str(week_index)]).fitted_delta_OD(:, 2), ...
        '-o', 'LineWidth', 1)
    
    % calculate R-square and rmspe (root mean square percentage error)
    % R-square
    R_square = calculate_Rsquare(liveData_used.(subject{sheet}).(['week_', num2str(week_index)]).fitted_delta_OD(:, 2), ...
                                 liveData_used.(subject{sheet}).(['week_', num2str(week_index)]).delta_OD(:, 2));
    liveData_used.(subject{sheet}).(['week_', num2str(week_index)]).R_square = R_square;
    
    % rmspe
    bound_percentile = 35;
    delta_OD_lowerBound = prctile(abs(liveData_used.(subject{sheet}).(['week_', num2str(week_index)]).delta_OD(:, 2)), bound_percentile);
    cal_points = find(abs(liveData_used.(subject{sheet}).(['week_', num2str(week_index)]).delta_OD(:, 2)) > delta_OD_lowerBound);
    
    rmspe = calculate_rmspe(liveData_used.(subject{sheet}).(['week_', num2str(week_index)]).fitted_delta_OD(cal_points, 2),...
                            liveData_used.(subject{sheet}).(['week_', num2str(week_index)]).delta_OD(cal_points, 2));
    liveData_used.(subject{sheet}).(['week_', num2str(week_index)]).rmspe = rmspe;
    
    % adjust appearance
    digit_num = 3;
    yline(delta_OD_lowerBound, '--k', 'LineWidth', 2);
    % ylim([-1, 15]);
    xlabel('wavelength (nm)')
    ylabel('{\Delta}OD (-)')
    legend('measured {\Delta}OD', 'fitted {\Delta}OD', ['rmspe calculation Bound' 10 '(measured {\Delta}OD = +/-' num2str(round(delta_OD_lowerBound, digit_num)) ')']);
    title([subject{sheet}, ' - {\Delta}OD fitting result, ', num2str(wl_res),'nm (week\_', num2str(week_index), ')'])
    fig_title = [subject{sheet}, '_sameScale_add_Hb_', num2str(wl_res), 'nm_deltaOD_fitting_week_', num2str(week_index), '_addrmspeLowerBound_shift_', num2str(wl_shift)];
    txt = {['R\_square = ', num2str(round(R_square, digit_num))],['rmspe = ', num2str(round(rmspe, digit_num)), '%']};
    text(0.68, 0.65, txt, 'Units', 'normalized');
    saveas(gcf,[fig_title, '.png']);
    hold off
end

%% plot delta_mua v.s. fitted_delta_mua
sheet = 3; week_num = 4;
for week_index = 1:week_num
    % plot delta_mua
    plot(liveData_used.(subject{sheet}).(['week_', num2str(week_index)]).delta_mua(:, 1), ...
         liveData_used.(subject{sheet}).(['week_', num2str(week_index)]).delta_mua(:, 2), ...
        '-o', 'LineWidth', 1)
    
    hold on
    
    % plot fitted delta_mua
    plot(liveData_used.(subject{sheet}).(['week_', num2str(week_index)]).fitted_delta_mua(:, 1), ...
        liveData_used.(subject{sheet}).(['week_', num2str(week_index)]).fitted_delta_mua(:, 2), ...
        '-o', 'LineWidth', 1)
    
    % calculate R-square and rmspe (root mean square percentage error)
    % R-square
    R_square = calculate_Rsquare(liveData_used.(subject{sheet}).(['week_', num2str(week_index)]).fitted_delta_mua(:, 2), ...
                                 liveData_used.(subject{sheet}).(['week_', num2str(week_index)]).delta_mua(:, 2));
    liveData_used.(subject{sheet}).(['week_', num2str(week_index)]).R_square = R_square;
    
    % rmspe
    bound_percentile = 35;
    delta_mua_lowerBound = prctile(abs(liveData_used.(subject{sheet}).(['week_', num2str(week_index)]).delta_mua(:, 2)), bound_percentile);
    cal_points = find(abs(liveData_used.(subject{sheet}).(['week_', num2str(week_index)]).delta_mua(:, 2)) > delta_mua_lowerBound);

    rmspe = calculate_rmspe(liveData_used.(subject{sheet}).(['week_', num2str(week_index)]).fitted_delta_mua(cal_points, 2),...
                            liveData_used.(subject{sheet}).(['week_', num2str(week_index)]).delta_mua(cal_points, 2));
    liveData_used.(subject{sheet}).(['week_', num2str(week_index)]).rmspe = rmspe;
    
    % adjust appearance
    digit_num = 3;
    yline(delta_mua_lowerBound, '--k', 'LineWidth', 2);
    % ylim([-1, 15]);
    xlabel('wavelength (nm)')
    ylabel('{\Delta}mua (1/cm)')
    legend('measured {\Delta}mua', 'fitted {\Delta}mua', ['rmspe calculation Bound' 10 '(measured {\Delta}mua = +/-' num2str(round(delta_mua_lowerBound, digit_num)) ')'], 'Location','northeast');
    title([subject{sheet}, ' - {\Delta}mua fitting result, ', num2str(wl_res),'nm (week\_', num2str(week_index), ')'])    
    fig_title = [subject{sheet}, '_sameScale_add_Hb_', num2str(wl_res), 'nm_deltaMua_fitting_week_', num2str(week_index), '_addrmspeLowerBound_shift_', num2str(wl_shift)];
    txt = {['R\_square = ', num2str(round(R_square, digit_num))],['rmspe = ', num2str(round(rmspe, digit_num)), '%']};
    text(0.65, 0.6, txt, 'Units', 'normalized');
    saveas(gcf,[fig_title, '.png']);
    hold off
end

%% plot measured delta_mua - fitted delta_mua
sheet = 3; week_num = 4;
for week_index = 1:week_num
    % plot measured - fitted
    plot(liveData_used.(subject{sheet}).(['week_', num2str(week_index)]).delta_mua(:, 1), ...
         liveData_used.(subject{sheet}).(['week_', num2str(week_index)]).delta_mua(:, 2) - ...
         liveData_used.(subject{sheet}).(['week_', num2str(week_index)]).fitted_delta_mua(:, 2), ...
         '-o', 'LineWidth', 1)
    
    % adjust appearance
    xlabel('wavelength (nm)')
    ylabel('measured {\Delta}mua - fitted {\Delta}mua (1/cm)')
    legend('difference')
    fig_title = [subject{sheet}, ' - Evaluation of residuals (measured - fitted), ', num2str(wl_res),'nm (week\_', num2str(week_index), ')'];
    title(fig_title) 
    saveas(gcf, [[subject{sheet}, ' - Evaluation of residuals (measured - fitted), ', num2str(wl_res),'nm (week_', num2str(week_index), ')_shift_', num2str(wl_shift)], '.png']);
end

%% plot each week fitting error ---- no used
xtick_names = {'week\_1', 'week\_2', 'week\_3', 'week\_4', 'week\_5'};
color_set = {'#FFD2D2', '#FD8A8A', '#FD0000', '#747474', '#000000'};

for index = 1:week_num
    liveData_used.(subject{sheet}).error(index) = liveData_used.(subject{sheet}).(['week_', num2str(index)]).('error');
end

plot(liveData_used.(subject{sheet}).error, '-o')
set(gca, 'xtick', (1:week_num), 'xticklabel', xtick_names);
xlabel('Week')
ylabel('Least square error')
title('Each week error')

%% compare chromophore variation between different subject
display_name = {'HbO2', 'Hb', 'Bilirubin', 'Biliverdin', 'Methemoglobin'};
t = tiledlayout(2, 3);
for chro_idx=1:length(display_name)
    nexttile(chro_idx)
    for sheet_idx=1:2:3        
        if sheet_idx == 1
            week_num = 5;
            plot(1:week_num, liveData_used.(subject{sheet_idx}).(display_name{chro_idx}), '-o', 'LineWidth', 2)
            hold on
        end
        if sheet_idx == 3
            week_num = 4;
            plot(1:week_num, liveData_used.(subject{sheet_idx}).(display_name{chro_idx}), '-o', 'LineWidth', 2)
            hold on
        end        
    end
    yline(0, '--k', 'LineWidth', 2);
    xlim([0.5, 5+0.5]);
    xlabel('Time (week)')
    ylabel('{\Delta}C (mol/L)')    
    set(gca, 'xtick', (1:5), 'xticklabel', {'week\_1', 'week\_2', 'week\_3', 'week\_4', 'week\_5'});
    legend('FY', 'ZO', '{\Delta}C = 0');
    title(display_name{chro_idx})
    grid on
    hold off
end
title(t, ['Compare chromophore {\Delta}C between subject FY & ZO (', num2str(wl_res), 'nm)']);






