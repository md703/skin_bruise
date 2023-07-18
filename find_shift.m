%% find the max value and corresponding id of HbO2 spectrum
% set wl pivot
HbO2_wl_1 = 500;
HbO2_wl_2 = 550;
HbO2_wl_3 = 600;

% find id of wl pivot
HbO2_wl_1_id = find(epsilon_interp{1}(:, 1) == HbO2_wl_1);
HbO2_wl_2_id = find(epsilon_interp{1}(:, 1) == HbO2_wl_2);
HbO2_wl_3_id = find(epsilon_interp{1}(:, 1) == HbO2_wl_3);

% get peak of HbO2, and corresponding wavelength
[HbO2_peak_1_value, HbO2_peak_1_id] = max(epsilon_interp{1}(HbO2_wl_1_id:HbO2_wl_2_id, 2));
[HbO2_peak_2_value, HbO2_peak_2_id] = max(epsilon_interp{1}(HbO2_wl_2_id:HbO2_wl_3_id, 2));
HbO2_wl_Peak_1_id = HbO2_wl_1_id + HbO2_peak_1_id - 1;
HbO2_wl_Peak_2_id = HbO2_wl_2_id + HbO2_peak_2_id - 1;

HbO2_wl_peak = [epsilon_interp{1}(HbO2_wl_Peak_1_id, 1), ...
                epsilon_interp{1}(HbO2_wl_Peak_2_id, 1)]; % first wl is first peak, second wl is second peak       
HbO2_value_peak = [HbO2_peak_1_value, HbO2_peak_2_value]; % first wl is first peak, second wl is second peak

%% set target
Target = 'mua'; % ('OD' or 'mua')

%% find the max value and corresponding id of FY measured delta_mua
FY_mea_wl = []; FY_mea_value = [];
% set pivot
sheet = 1; week_num = 5;
FY_mea_wl_1 = 525;
FY_mea_wl_2 = 560;
FY_mea_wl_3 = 600;

% find id of wl pivot
FY_mea_wl_1_id = find(target_wl == FY_mea_wl_1);
FY_mea_wl_2_id = find(target_wl == FY_mea_wl_2);
FY_mea_wl_3_id = find(target_wl == FY_mea_wl_3);

for week_index = 2:week_num
    % get peak of FY_delta_mua, and corresponding wavelength
    [FY_mea_peak_1_value, FY_mea_peak_1_id] = max(liveData_used.(subject{sheet}).(['week_', num2str(week_index)]).(['delta_', Target])(FY_mea_wl_1_id:FY_mea_wl_2_id, 2));
    [FY_mea_peak_2_value, FY_mea_peak_2_id] = max(liveData_used.(subject{sheet}).(['week_', num2str(week_index)]).(['delta_', Target])(FY_mea_wl_2_id:FY_mea_wl_3_id, 2));
    FY_mea_wl_Peak_1_id = FY_mea_wl_1_id + FY_mea_peak_1_id - 1;
    FY_mea_wl_Peak_2_id = FY_mea_wl_2_id + FY_mea_peak_2_id - 1;

    FY_mea_wl_peak = [liveData_used.(subject{sheet}).(['week_', num2str(week_index)]).(['delta_', Target])(FY_mea_wl_Peak_1_id, 1), ...
                      liveData_used.(subject{sheet}).(['week_', num2str(week_index)]).(['delta_', Target])(FY_mea_wl_Peak_2_id, 1)]; % first wl is first peak, second wl is second peak       
    FY_mea_value_peak = [FY_mea_peak_1_value, FY_mea_peak_2_value]; % first wl is first peak, second wl is second peak

    FY_mea_wl = [FY_mea_wl; FY_mea_wl_peak];
    FY_mea_value = [FY_mea_value; FY_mea_value_peak];
end

%% find the max value and corresponding id of FY fitted delta_mua
FY_fit_wl = []; FY_fit_value = [];
% set pivot
sheet = 1; week_num = 5;
FY_fit_wl_1 = 525;
FY_fit_wl_2 = 560;
FY_fit_wl_3 = 600;

% find id of wl pivot
FY_fit_wl_1_id = find(target_wl == FY_fit_wl_1);
FY_fit_wl_2_id = find(target_wl == FY_fit_wl_2);
FY_fit_wl_3_id = find(target_wl == FY_fit_wl_3);

for week_index = 2:week_num
    % get peak of FY_delta_mua, and corresponding wavelength
    [FY_fit_peak_1_value, FY_fit_peak_1_id] = max(liveData_used.(subject{sheet}).(['week_', num2str(week_index)]).(['fitted_delta_', Target])(FY_fit_wl_1_id:FY_fit_wl_2_id, 2));
    [FY_fit_peak_2_value, FY_fit_peak_2_id] = max(liveData_used.(subject{sheet}).(['week_', num2str(week_index)]).(['fitted_delta_', Target])(FY_fit_wl_2_id:FY_fit_wl_3_id, 2));
    FY_fit_wl_Peak_1_id = FY_fit_wl_1_id + FY_fit_peak_1_id - 1;
    FY_fit_wl_Peak_2_id = FY_fit_wl_2_id + FY_fit_peak_2_id - 1;

    FY_fit_wl_peak = [liveData_used.(subject{sheet}).(['week_', num2str(week_index)]).(['fitted_delta_', Target])(FY_fit_wl_Peak_1_id, 1), ...
                      liveData_used.(subject{sheet}).(['week_', num2str(week_index)]).(['fitted_delta_', Target])(FY_fit_wl_Peak_2_id, 1)]; % first wl is first peak, second wl is second peak       
    FY_fit_value_peak = [FY_fit_peak_1_value, FY_fit_peak_2_value]; % first wl is first peak, second wl is second peak

    FY_fit_wl = [FY_fit_wl; FY_fit_wl_peak];
    FY_fit_value = [FY_fit_value; FY_fit_value_peak];
end

%% find the max value and corresponding id of ZO measured delta_mua
ZO_mea_wl = []; ZO_mea_value = [];
% set pivot
sheet = 3; week_num = 4;
ZO_mea_wl_1 = 525;
ZO_mea_wl_2 = 560;
ZO_mea_wl_3 = 600;

% find id of wl pivot
ZO_mea_wl_1_id = find(target_wl == ZO_mea_wl_1);
ZO_mea_wl_2_id = find(target_wl == ZO_mea_wl_2);
ZO_mea_wl_3_id = find(target_wl == ZO_mea_wl_3);

for week_index = 1:week_num-1
    % get peak of ZO_delta_mua, and corresponding wavelength
    [ZO_mea_peak_1_value, ZO_mea_peak_1_id] = max(liveData_used.(subject{sheet}).(['week_', num2str(week_index)]).(['delta_', Target])(ZO_mea_wl_1_id:ZO_mea_wl_2_id, 2));
    [ZO_mea_peak_2_value, ZO_mea_peak_2_id] = max(liveData_used.(subject{sheet}).(['week_', num2str(week_index)]).(['delta_', Target])(ZO_mea_wl_2_id:ZO_mea_wl_3_id, 2));
    ZO_mea_wl_Peak_1_id = ZO_mea_wl_1_id + ZO_mea_peak_1_id - 1;
    ZO_mea_wl_Peak_2_id = ZO_mea_wl_2_id + ZO_mea_peak_2_id - 1;

    ZO_mea_wl_peak = [liveData_used.(subject{sheet}).(['week_', num2str(week_index)]).(['delta_', Target])(ZO_mea_wl_Peak_1_id, 1), ...
                      liveData_used.(subject{sheet}).(['week_', num2str(week_index)]).(['delta_', Target])(ZO_mea_wl_Peak_2_id, 1)]; % first wl is first peak, second wl is second peak       
    ZO_mea_value_peak = [ZO_mea_peak_1_value, ZO_mea_peak_2_value]; % first wl is first peak, second wl is second peak

    ZO_mea_wl = [ZO_mea_wl; ZO_mea_wl_peak];
    ZO_mea_value = [ZO_mea_value; ZO_mea_value_peak];
end

%% find the max value and corresponding id of ZO fitted delta_mua
ZO_fit_wl = []; ZO_fit_value = [];
% set pivot
sheet = 3; week_num = 4;
ZO_fit_wl_1 = 530;
ZO_fit_wl_2 = 560;
ZO_fit_wl_3 = 600;

% find id of wl pivot
ZO_fit_wl_1_id = find(target_wl == ZO_fit_wl_1);
ZO_fit_wl_2_id = find(target_wl == ZO_fit_wl_2);
ZO_fit_wl_3_id = find(target_wl == ZO_fit_wl_3);

for week_index = 1:week_num-1
    % get peak of ZO_delta_mua, and corresponding wavelength
    [ZO_fit_peak_1_value, ZO_fit_peak_1_id] = max(liveData_used.(subject{sheet}).(['week_', num2str(week_index)]).(['fitted_delta_', Target])(ZO_fit_wl_1_id:ZO_fit_wl_2_id, 2));
    [ZO_fit_peak_2_value, ZO_fit_peak_2_id] = max(liveData_used.(subject{sheet}).(['week_', num2str(week_index)]).(['fitted_delta_', Target])(ZO_fit_wl_2_id:ZO_fit_wl_3_id, 2));
    ZO_fit_wl_Peak_1_id = ZO_fit_wl_1_id + ZO_fit_peak_1_id - 1;
    ZO_fit_wl_Peak_2_id = ZO_fit_wl_2_id + ZO_fit_peak_2_id - 1;

    ZO_fit_wl_peak = [liveData_used.(subject{sheet}).(['week_', num2str(week_index)]).(['fitted_delta_', Target])(ZO_fit_wl_Peak_1_id, 1), ...
                      liveData_used.(subject{sheet}).(['week_', num2str(week_index)]).(['fitted_delta_', Target])(ZO_fit_wl_Peak_2_id, 1)]; % first wl is first peak, second wl is second peak       
    ZO_fit_value_peak = [ZO_fit_peak_1_value, ZO_fit_peak_2_value]; % first wl is first peak, second wl is second peak

    ZO_fit_wl = [ZO_fit_wl; ZO_fit_wl_peak];
    ZO_fit_value = [ZO_fit_value; ZO_fit_value_peak];
end

%% calculate shift amount (delta_mua)
FY_shift = FY_mea_wl - HbO2_wl_peak;
ZO_shift = ZO_mea_wl - HbO2_wl_peak;


