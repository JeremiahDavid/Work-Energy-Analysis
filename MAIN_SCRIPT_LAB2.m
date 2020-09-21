%% start Program
clear
clc
close all
fig = 1;
%% Part 2 ( torque vs crank shaft angle)
load('IC_Engines_air_injection.mat');
% air from 55-75 plots
timer = IC_Engines_air_injection.time.time;
% displacement data
disp_lowr = IC_Engines_air_injection.displacement.low;
disp_medr = IC_Engines_air_injection.displacement.med;
disp_hir = IC_Engines_air_injection.displacement.high;
disp_fwr = IC_Engines_air_injection.displacement.flywheel;
disp_vecr = [disp_lowr, disp_medr, disp_hir, disp_fwr];
% angle data
angle_lowr = IC_Engines_air_injection.angle.low;
angle_medr = IC_Engines_air_injection.angle.med;
angle_hir = IC_Engines_air_injection.angle.high;
angle_fwr = IC_Engines_air_injection.angle.flywheel;
angle_vecr = [angle_lowr, angle_medr, angle_hir, angle_fwr];
% torque data
torque_lowr = IC_Engines_air_injection.torque.low;
torque_medr = IC_Engines_air_injection.torque.med;
torque_hir = IC_Engines_air_injection.torque.high;
torque_fwr = IC_Engines_air_injection.torque.flywheel;
torque_vecr = [torque_lowr, torque_medr, torque_hir, torque_fwr];

% Plot torque Data
figure(fig)
fig = fig + 1;
hold on
scatter(angle_lowr,torque_lowr,'r.','DisplayName', 'Low RPM')
scatter(angle_medr,torque_medr,'b.','DisplayName', 'Medium RPM')
scatter(angle_hir,torque_hir,'g.','DisplayName', 'High RPM')
legend
xlabel('Angle [deg]')
ylabel('Torque [N-m]')

[torque_low, angle_lowt] = data_fixer_torque(angle_lowr,torque_lowr);
[torque_med, angle_medt] = data_fixer_torque(angle_medr,torque_medr);
[torque_hi, angle_hit] = data_fixer_torque(angle_hir,torque_hir);
figure(fig)
fig = fig + 1;
hold on
scatter(angle_lowt,torque_low,'r.','DisplayName', 'Low RPM')
scatter(angle_medt,torque_med,'b.','DisplayName', 'Medium RPM')
scatter(angle_hit,torque_hi,'g.','DisplayName', 'High RPM')
legend
xlabel('Angle [deg]')
ylabel('Torque [N-m]')

% Plot the angular position data
figure(fig)
fig = fig + 1;
hold on
scatter(timer,angle_lowr,'r.','DisplayName', 'Low RPM')
scatter(timer,angle_medr,'b.','DisplayName', 'Medium RPM')
scatter(timer,angle_hir,'g.','DisplayName', 'High RPM')
legend
xlabel('Time [s]')
ylabel('Angle [deg]')

% Plot the linear and adjusted angular position data
[time_low, angle_lowa] = data_fixer(timer,angle_lowr);
[time_med, angle_meda] = data_fixer(timer,angle_medr);
[time_hi, angle_hia] = data_fixer(timer,angle_hir);

figure(fig)
fig = fig + 1;
hold on
scatter(time_low,angle_lowa,'r.','DisplayName', 'Low RPM')
scatter(time_med,angle_meda,'b.','DisplayName', 'Medium RPM')
scatter(time_hi,angle_hia,'g.','DisplayName', 'High RPM')
legend
xlabel('Time [s]')
ylabel('Angle [deg]')

% Plot the displacement data
figure(fig)
fig = fig + 1;
hold on
scatter(timer,disp_lowr,'r.','DisplayName', 'Low RPM')
scatter(timer,disp_medr,'b.','DisplayName', 'Medium RPM')
scatter(timer,disp_hir,'g.','DisplayName', 'High RPM')
legend
xlabel('Time [s]')
ylabel('Displacememt [mm]')

figure(fig)
fig = fig + 1;
hold on
% scatter(angle_lowr,disp_lowr,'r.','DisplayName', 'Low RPM')
scatter(angle_medr,disp_medr,'b.','DisplayName', 'Medium RPM')
% scatter(angle_hir,disp_hir,'g.','DisplayName', 'High RPM')
legend
xlabel('Angle [deg]')
ylabel('Displacememt [mm]')

%% Work from air
    % function to calculate the displacement
tdc_m = .00065; 
bdc_m = .03837;
piston_m = @(x) .5*(bdc_m - tdc_m).*sind(x - 90) + .5*(bdc_m - tdc_m);
    % constants
bore = 0.05198;
stroke = 0.03772;
gamma = 1.4;
P_atm = 101325; % pascals
P_line = 90; % psi
P_line = P_line*6894.76; % pascals

    % compression stroke
Pa_c = P_atm; % pascals
vh1_vh2 = 1e-5;
vh3 = 0.56*pi*0.0095^2/4;
V_tdc = vh1_vh2 + vh3;
V_bdc = V_tdc + stroke*(pi*bore^2)/4;
Pb_c = (Pa_c*V_bdc^gamma)/V_tdc^gamma - P_atm;
W_air_comp = (Pa_c*(V_bdc)^gamma/(1-gamma))*(V_tdc^(1-gamma) ...
    - V_bdc^(1-gamma)) - P_atm*(V_tdc - V_bdc);

    % power stroke
disp_ab = piston_m(180) - piston_m(190); % distance the piston travels from tdc to 10 deg atdc
disp_bc = piston_m(190) - piston_m(260); % distance the piston travels from 10 deg atdc to 80 deg atdc

V_a = V_tdc;
V_b = V_tdc + disp_ab*(pi*bore^2)/4;
V_c = V_tdc + disp_bc*(pi*bore^2)/4;
V_d = V_bdc;
Pa_p = Pb_c;
Pb_p = P_line-P_atm;
Pc_p = P_line-P_atm;

W_air_power_ab = (Pa_p*(V_a)^gamma/(1-gamma))*(V_b^(1-gamma) ...
    - V_a^(1-gamma)) - P_atm*(V_b - V_a);
W_air_power_bc = (Pb_p)*(V_c - V_b);
W_air_power_cd = (Pc_p*(V_c)^gamma/(1-gamma))*(V_d^(1-gamma) ...
    - V_c^(1-gamma)) - P_atm*(V_d - V_c);
W_air_power = W_air_power_ab + W_air_power_bc + W_air_power_cd;

W_air_intake = 0;
W_air_exhaust = 0;
% create a struct for work of air
W_air.intake = W_air_intake;
W_air.comp = W_air_comp;
W_air.power = W_air_power;
W_air.exhaust = W_air_exhaust;


%% Find the angular velocity for the begining and end of the 4 strokes
% angular velocities are roughly constant between 400 and 600 deg.
% estimate the rpm for each trial to be used in the model 
positions = {};
positions_base = [720-350, 180+720-350, 360+720-350, 540+720-350, 720+720-350];
positions_fw = [680-350, 180+680-350, 360+680-350, 540+680-350, 720+680-350];
positions{1} = positions_base;
positions{2} = positions_base;
positions{3} = positions_base;
positions{4} = positions_fw;

[time_lowd, angle_lowd] = data_fixer(timer,angle_lowr);
angle_lowd_rad = angle_lowd(1:end-1);
omega_low = (pi/180)*diff(angle_lowd)./diff(time_lowd);

[time_medd, angle_medd] = data_fixer(timer,angle_medr);
angle_medd_rad = angle_medd(1:end-1);
omega_med = (pi/180)*diff(angle_medd)./diff(time_medd);

[time_hid, angle_hid] = data_fixer(timer,angle_hir);
angle_hid_rad = angle_hid(1:end-1);
omega_hi = (pi/180)*diff(angle_hid)./diff(time_hid);

[time_fwd, angle_fwd] = data_fixer(timer,angle_fwr);
angle_fwd_rad = angle_fwd(1:end-1);
omega_fw = (pi/180)*diff(angle_fwd)./diff(time_fwd);

% Plot the angular velocity data
figure(fig)
fig = fig + 1;
hold on
scatter(angle_lowd_rad,omega_low,'r.','DisplayName', 'Low RPM')
scatter(angle_medd_rad,omega_med,'b.','DisplayName', 'Medium RPM')
scatter(angle_hid_rad,omega_hi,'g.','DisplayName', 'High RPM')
legend
xlabel('Angle [deg]')
ylabel('Angular Velocity [rad/s]')

% for power stroke: tdc -> bdc
omega_data_comp = {};
[omega_low_comp_bdc, omega_low_comp_bdc_ind] = max(omega_low(30:60));
[omega_low_comp_tdc, omega_low_comp_tdc_ind] = min(omega_low(30:60));
omega_data_comp = [omega_data_comp, [omega_low_comp_bdc, omega_low_comp_tdc]];
[omega_med_comp_bdc, omega_med_comp_bdc_ind]  = max(omega_med(30:60));
[omega_med_comp_tdc, omega_med_comp_tdc_ind] = min(omega_med(30:60));
omega_data_comp = [omega_data_comp, [omega_med_comp_bdc, omega_med_comp_tdc]];
[omega_hi_comp_bdc, omega_hi_comp_bdc_ind] = max(omega_hi(15:60));
[omega_hi_comp_tdc, omega_hi_comp_tdc_ind] = min(omega_hi(15:60));
omega_data_comp = [omega_data_comp, [omega_hi_comp_bdc, omega_hi_comp_tdc]];
[omega_fw_comp_bdc, omega_fw_comp_bdc_ind] = max(omega_fw(15:60));
[omega_fw_comp_tdc, omega_fw_comp_tdc_ind] = min(omega_fw(15:60));
omega_data_comp = [omega_data_comp, [omega_fw_comp_bdc, omega_fw_comp_tdc]];

% omega_low_comp_bdc = angle_speed_low(find(angle_test_low>=positions_low(2),1));
% omega_low_comp_tdc = angle_speed_low(find(angle_test_low>=positions_low(3),1));
% omega_med_comp_bdc = angle_speed_med(find(angle_test_med>=positions_med(2),1));
% omega_med_comp_tdc = angle_speed_med(find(angle_test_med>=positions_med(3),1));
% omega_hi_comp_bdc = angle_speed_hi(find(angle_test_hi>=positions_hi(2),1));
% omega_hi_comp_tdc = angle_speed_hi(find(angle_test_hi>=positions_hi(3),1));

% for power stroke: bdc -> tdc
omega_data_power = {};
omega_low_power_tdc = min(omega_low(omega_low_comp_tdc_ind:omega_low_comp_tdc_ind+60));
omega_low_power_bdc = max(omega_low(omega_low_comp_tdc_ind:omega_low_comp_tdc_ind+60));
omega_data_power = [omega_data_power, [omega_low_power_tdc, omega_low_power_bdc]];
omega_med_power_tdc = min(omega_med(omega_med_comp_tdc_ind:omega_med_comp_tdc_ind+60));
omega_med_power_bdc = max(omega_med(omega_med_comp_tdc_ind:omega_med_comp_tdc_ind+60));
omega_data_power = [omega_data_power, [omega_med_power_tdc, omega_med_power_bdc]];
omega_hi_power_tdc = min(omega_hi(omega_hi_comp_tdc_ind:omega_hi_comp_tdc_ind+60));
omega_hi_power_bdc = max(omega_hi(omega_hi_comp_tdc_ind:omega_hi_comp_tdc_ind+60));
omega_data_power = [omega_data_power, [omega_hi_power_tdc, omega_hi_power_bdc]];
omega_fw_power_tdc = min(omega_fw(omega_fw_comp_tdc_ind:omega_fw_comp_tdc_ind+60));
omega_fw_power_bdc = max(omega_fw(omega_fw_comp_tdc_ind:omega_fw_comp_tdc_ind+60));
omega_data_power = [omega_data_power, [omega_fw_power_tdc, omega_fw_power_bdc]];

% omega_low_power_tdc = angle_speed_low(find(angle_test_low>=positions_low(3),1));
% omega_low_power_bdc = angle_speed_low(find(angle_test_low>=positions_low(4),1));
% omega_med_power_tdc = angle_speed_med(find(angle_test_med>=positions_med(3),1));
% omega_med_power_bdc = angle_speed_med(find(angle_test_med>=positions_med(4),1));
% omega_hi_power_tdc = angle_speed_hi(find(angle_test_hi>=positions_hi(3),1));
% omega_hi_power_bdc = angle_speed_hi(find(angle_test_hi>=positions_hi(4),1));

%% Calculate the change in energy for the 4 strokes
I_fly = 0.0114; %from part 1 calculations
I_fly_1fw = I_fly/4;
delta_E_intake = 0;
delta_E_exhaust = 0;
delta_E_comp_vec = [];
for i = 1:length(omega_data_comp)-1
    E_a_comp = 0.5*I_fly*omega_data_comp{i}(1)^2; 
    E_b_comp = 0.5*I_fly*omega_data_comp{i}(2)^2;
    delta_E_comp = E_b_comp - E_a_comp;
    delta_E_comp_vec = [delta_E_comp_vec, delta_E_comp];
end
delta_E_power_vec = [];
for i = 1:length(omega_data_comp)-1
    E_a_power = 0.5*I_fly*omega_data_power{i}(1)^2; 
    E_b_power = 0.5*I_fly*omega_data_power{i}(2)^2;
    delta_E_power = E_b_power - E_a_power;
    delta_E_power_vec = [delta_E_power_vec, delta_E_power];
end
E_a_comp = 0.5*I_fly_1fw*omega_data_comp{end}(1)^2;
E_b_comp = 0.5*I_fly_1fw*omega_data_comp{end}(2)^2;
delta_E_comp_1fw = E_b_comp - E_a_comp;

E_a_power = 0.5*I_fly_1fw*omega_data_power{end}(1)^2;
E_b_power = 0.5*I_fly_1fw*omega_data_power{end}(2)^2;
delta_E_power_1fw = E_b_power - E_a_power;

% combine all the change in energy data into a struct
delta_E.comp.low = delta_E_comp_vec(1);
delta_E.comp.med = delta_E_comp_vec(2);
delta_E.comp.hi = delta_E_comp_vec(3);
delta_E.comp.fw = delta_E_comp_1fw;

delta_E.power.low = delta_E_power_vec(1);
delta_E.power.med = delta_E_power_vec(2);
delta_E.power.hi = delta_E_power_vec(3);
delta_E.power.fw = delta_E_power_1fw;

%% Find the Frictional Torque
% use part 1 frictional torque values and rpm's to create a linear
% regession modeling torque as a function of rpm
frictional_torque_rpm = [200;300;500];
frictoinal_torque = [-0.7262;-0.8148;-0.9341];
[torq_fric_model, stats] = fit(frictional_torque_rpm,frictoinal_torque,'poly1');
fric_torque = @(rpm) torq_fric_model.p1*rpm + torq_fric_model.p2;

rpm_plot = linspace(frictional_torque_rpm(1),frictional_torque_rpm(end),1001);
fric_torque_plot = fric_torque(rpm_plot);
figure(fig)
fig = fig + 1;
hold on
scatter(frictional_torque_rpm,frictoinal_torque,'DisplayName','Data from Week 1')
plot(rpm_plot,fric_torque_plot,'r--','DisplayName','Linear Regression')

start_low = find(angle_lowa>positions_base(1),1);
finish_low = find(angle_lowa>positions_base(2),1);
low_rpm = (1/6)*(angle_lowa(finish_low) - angle_lowa(start_low))/ ...
    (time_low(finish_low) - time_low(start_low));
start_med = find(angle_meda>positions_base(1),1);
finish_med = find(angle_meda>positions_base(2),1);
med_rpm = (1/6)*(angle_meda(finish_med) - angle_meda(start_med))/ ...
    (time_med(finish_med) - time_med(start_med));
start_hi = find(angle_hia>positions_base(1),1);
finish_hi = find(angle_hia>positions_base(2),1);
hi_rpm = (1/6)*(angle_hia(finish_hi) - angle_hia(start_hi))/ ...
    (time_hi(finish_hi) - time_hi(start_hi));

% use the model to estimate the frictional torque
frictional_torque_low = fric_torque(low_rpm);
frictional_torque_med = fric_torque(med_rpm);
frictional_torque_hi = fric_torque(hi_rpm);
frictional_torques = [frictional_torque_low,frictional_torque_med, ...
    frictional_torque_hi,frictional_torque_hi];

scatter([low_rpm;med_rpm;hi_rpm], ...
    [frictional_torque_low;frictional_torque_med;frictional_torque_hi], ...
    'k*','DisplayName','Interpolated values')
legend
xlabel('RPM')
ylabel('Frictional Torque')

%% calculate the work of the motor for each of the 4 strokes
% set indicators for stroke
intake = 1;
comp = 2;
power = 3;
exhaust = 4;
frictional_torque_intake = 1;
W_motor_intake_vec = {};
W_motor_comp_vec = {};
W_motor_power_vec = {};
W_motor_exhaust_vec = {};

[~, iteration_length] = size(torque_vecr);
for i = 1:iteration_length
    [W_motor_intake, W_fric_intake, fric_torque, angle_W_motor_intake] = ...
        data_fixer_work_stroke(i,intake,positions{i},frictional_torques(i), ...
        angle_vecr(:,i),torque_vecr(:,i));
    W_motor_intake = W_motor_intake(end);
    W_motor_intake_vec = [W_motor_intake_vec, [W_motor_intake,W_fric_intake]];
    [W_motor_comp, W_fric_comp, ~, angle_W_motor_comp] = ...
        data_fixer_work_stroke(i,comp,positions{i},frictional_torques(i), ...
        angle_vecr(:,i),torque_vecr(:,i));
    W_motor_comp = W_motor_comp(end);
    W_motor_comp_vec = [W_motor_comp_vec, [W_motor_comp,W_fric_comp]];
    [W_motor_power, W_fric_power, ~, angle_W_motor_power] = ...
        data_fixer_work_stroke(i,power,positions{i},frictional_torques(i), ...
        angle_vecr(:,i),torque_vecr(:,i));
    W_motor_power = W_motor_power(end);
    W_motor_power_vec = [W_motor_power_vec, [W_motor_power,W_fric_power]];
    [W_motor_exhaust, W_fric_exhaust, ~, angle_W_motor_exhaust] = ...
        data_fixer_work_stroke(i,exhaust,positions{i},frictional_torques(i), ...
       angle_vecr(:,i),torque_vecr(:,i));
    W_motor_exhaust = W_motor_exhaust(end);
    W_motor_exhaust_vec = [W_motor_exhaust_vec, [W_motor_exhaust,W_fric_exhaust]];
end


% place engine work data in a struct
    % intake
W_motor.intake.low = W_motor_intake_vec{1}(1);
W_motor.intake.med = W_motor_intake_vec{2}(1);
W_motor.intake.hi = W_motor_intake_vec{3}(1);
W_motor.intake.fw = W_motor_intake_vec{4}(1);
    % compression
W_motor.comp.low = W_motor_comp_vec{1}(1);
W_motor.comp.med = W_motor_comp_vec{2}(1);
W_motor.comp.hi = W_motor_comp_vec{3}(1);
W_motor.comp.fw = W_motor_comp_vec{4}(1);  
    % power
W_motor.power.low = W_motor_power_vec{1}(1);
W_motor.power.med = W_motor_power_vec{2}(1);
W_motor.power.hi = W_motor_power_vec{3}(1);
W_motor.power.fw = W_motor_power_vec{4}(1);       
    % exhaust
W_motor.exhaust.low = W_motor_exhaust_vec{1}(1);
W_motor.exhaust.med = W_motor_exhaust_vec{2}(1);
W_motor.exhaust.hi = W_motor_exhaust_vec{3}(1);
W_motor.exhaust.fw = W_motor_exhaust_vec{4}(1); 


%% Calculate work due to friction for stroke
% place friction work data in a struct
    % intake
W_fric.intake.low = W_motor_intake_vec{1}(2);
W_fric.intake.med = W_motor_intake_vec{2}(2);
W_fric.intake.hi = W_motor_intake_vec{3}(2);
W_fric.intake.fw = W_motor_intake_vec{4}(2);
    % compression
W_fric.comp.low = W_motor_comp_vec{1}(2);
W_fric.comp.med = W_motor_comp_vec{2}(2);
W_fric.comp.hi = W_motor_comp_vec{3}(2);
W_fric.comp.fw = W_motor_comp_vec{4}(2); 
    % power
W_fric.power.low = W_motor_power_vec{1}(2);
W_fric.power.med = W_motor_power_vec{2}(2);
W_fric.power.hi = W_motor_power_vec{3}(2);
W_fric.power.fw = W_motor_power_vec{4}(2);       
    % exhaust
W_fric.exhaust.low = W_motor_exhaust_vec{1}(2);
W_fric.exhaust.med = W_motor_exhaust_vec{2}(2);
W_fric.exhaust.hi = W_motor_exhaust_vec{3}(2);
W_fric.exhaust.fw = W_motor_exhaust_vec{4}(2); 


%% solve for the Work of error
%     % intake error
% W_error_intake.low = delta_E_intake - W_air_intake - W_motor.intake.low - W_fric.intake.low;
% W_error_intake.med = delta_E_intake - W_air_intake - W_motor.intake.med - W_fric.intake.med;
% W_error_intake.hi = delta_E_intake - W_air_intake - W_motor.intake.hi - W_fric.intake.hi;
% W_error_intake.fw = delta_E_intake - W_air_intake - W_motor.intake.fw - W_fric.intake.fw;
% 
%     % compression error 
% W_error_comp.low = delta_E.comp.low - W_air.comp - W_motor.comp.low - W_fric.comp.low;
% W_error_comp.med = delta_E.comp.med - W_air.comp - W_motor.comp.med - W_fric.comp.med;
% W_error_comp.hi = delta_E.comp.hi - W_air.comp - W_motor.comp.hi - W_fric.comp.hi;
% W_error_comp.fw = delta_E.comp.fw - W_air.comp - W_motor.comp.fw - W_fric.comp.fw;
% 
%     % power error 
% W_error_power.low = delta_E.power.low - W_air.power - W_motor.power.low - W_fric.power.low;
% W_error_power.med = delta_E.power.med - W_air.power - W_motor.power.med - W_fric.power.med;
% W_error_power.hi = delta_E.power.hi - W_air.power - W_motor.power.hi - W_fric.power.hi;
% W_error_power.fw = delta_E.power.fw - W_air.power - W_motor.power.fw - W_fric.power.fw;
% 
%     % exhaust error
% W_error_exhaust.low = delta_E_exhaust - W_air_exhaust - W_motor.exhaust.low - W_fric.exhaust.low;
% W_error_exhaust.med = delta_E_exhaust - W_air_exhaust - W_motor.exhaust.med - W_fric.exhaust.med;
% W_error_exhaust.hi = delta_E_exhaust - W_air_exhaust - W_motor.exhaust.hi - W_fric.exhaust.hi;
% W_error_exhaust.fw = delta_E_exhaust - W_air_exhaust - W_motor.exhaust.fw - W_fric.exhaust.fw;
% 
% W_error.intake = W_error_intake;
% W_error.comp = W_error_comp;
% W_error.power = W_error_power;
% W_error.exhaust = W_error_exhaust;

%% Calculate the net work for a full cycle
    % intake 
W_total_low.intake = W_air.intake + W_motor.intake.low + W_fric.intake.low;% + W_error.intake.low;
W_total_med.intake = W_air.intake + W_motor.intake.med + W_fric.intake.med;% + W_error.intake.med;
W_total_hi.intake = W_air.intake + W_motor.intake.hi + W_fric.intake.hi;% + W_error.intake.hi;
W_total_fw.intake = W_air.intake + W_motor.intake.fw + W_fric.intake.fw;% + W_error.intake.fw;

    % compression
W_total_low.comp = W_air.comp + W_motor.comp.low + W_fric.comp.low;% + W_error.comp.low;
W_total_med.comp = W_air.comp + W_motor.comp.med + W_fric.comp.med;% + W_error.comp.med;
W_total_hi.comp = W_air.comp + W_motor.comp.hi + W_fric.comp.hi;% + W_error.comp.hi;
W_total_fw.comp = W_air.comp + W_motor.comp.fw + W_fric.comp.fw;% + W_error.comp.fw;

    % power
W_total_low.power = W_air.power + W_motor.power.low + W_fric.power.low;% + W_error.power.low;
W_total_med.power = W_air.power + W_motor.power.med + W_fric.power.med;% + W_error.power.med;
W_total_hi.power = W_air.power + W_motor.power.hi + W_fric.power.hi;% + W_error.power.hi;
W_total_fw.power = W_air.power + W_motor.power.fw + W_fric.power.fw;% + W_error.power.fw;

    % exhaust
W_total_low.exhaust = W_air.exhaust + W_motor.exhaust.low + W_fric.exhaust.low;% + W_error.exhaust.low;
W_total_med.exhaust = W_air.exhaust + W_motor.exhaust.med + W_fric.exhaust.med;% + W_error.exhaust.med;
W_total_hi.exhaust = W_air.exhaust + W_motor.exhaust.hi + W_fric.exhaust.hi;% + W_error.exhaust.hi;
W_total_fw.exhaust = W_air.exhaust + W_motor.exhaust.fw + W_fric.exhaust.fw;% + W_error.exhaust.fw;

    % full cycle
W_total_cycle_low = W_total_low.intake + W_total_low.comp + W_total_low.power + W_total_low.exhaust  
W_total_cycle_med = W_total_med.intake + W_total_med.comp + W_total_med.power + W_total_med.exhaust 
W_total_cycle_hi = W_total_hi.intake + W_total_hi.comp + W_total_hi.power + W_total_hi.exhaust 
W_total_cycle_fw = W_total_fw.intake + W_total_fw.comp + W_total_fw.power + W_total_fw.exhaust 

%% Net work from Energy
W_energy_cycle_low = delta_E_intake + delta_E.comp.low + delta_E.power.low + delta_E_exhaust
W_energy_cycle_med = delta_E_intake + delta_E.comp.med + delta_E.power.med + delta_E_exhaust
W_energy_cycle_hi = delta_E_intake + delta_E.comp.hi + delta_E.power.hi + delta_E_exhaust
W_energy_cycle_fw = delta_E_intake + delta_E.comp.fw + delta_E.power.fw + delta_E_exhaust

%% total net work
W_error_low = W_energy_cycle_low - W_total_cycle_low
W_error_med = W_energy_cycle_med - W_total_cycle_med 
W_error_hi = W_energy_cycle_hi - W_total_cycle_hi
W_error_fw = W_energy_cycle_fw - W_total_cycle_fw

load('IC_Engines_Results')
W_net_low = W_total_cycle_low-W_net_no_air.low;
W_net_med = W_total_cycle_med-W_net_no_air.med;





