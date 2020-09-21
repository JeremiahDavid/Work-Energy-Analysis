%% start Program
clear
clc
close all
%% Obtain Torque Equation
load('IC_Engines_Part1')
voltage_offset = 1.0204425;
mass_numbers = [0;1;2;3;4;5;6];
torque_Nm_cal = 0.2264*9.81*moment_m.*(mass_numbers);
% voltage = voltage - voltage_offset;
[cal, stats] = fit(voltage,torque_Nm_cal,'Poly1'); % calibration fit
% cal.p1
% torque = @(v) 1.64.*v -1.678;
torque = @(v) cal.p1.*(v - voltage_offset);
torque = @(v) cal.p1.*v + cal.p2;
%% Plot Calibration Data
volt_plot = linspace(1,2.6,1001);
torque_plot = torque(volt_plot);
up_bound = 2.*stats.rmse + torque_plot;
low_bound = torque_plot - 2.*stats.rmse;
figure(8)
hold on
p = zeros(1,3);
p(1) = scatter(voltage,torque_Nm_cal,'b*','DisplayName','Raw Calibration Data');
p(2) = plot(volt_plot,torque_plot,'r-','DisplayName','Calibration Regression');
p(3) = plot(volt_plot,up_bound,'g--','DisplayName','95.4% Confidence Bounds');
p(4) = plot(volt_plot,low_bound,'g--');
legend(p(1:3))
xlabel('Voltage [V]')
ylabel('Torque [N*m]')
%% Load Saved Data for Angle 
load('IC_Engines_Part2_1')
[time_ss_200, angle_ss_200] = data_fixer(data1.time_ss_200,data1.tooth_count_200.*10);
[time_ss_300, angle_ss_300] = data_fixer(data1.time_ss_300,data1.tooth_count_300.*10);
[time_ss_500, angle_ss_500] = data_fixer(data1.time_ss_500,data1.tooth_count_500.*10);
% Plot Steady State Data Angle
RPM_Engine_200 = sprintf('Steady State RPM: %.2f',200);
RPM_Engine_300 = sprintf('Steady State RPM: %.2f',300);
RPM_Engine_500 = sprintf('Steady State RPM: %.2f',500);
figure(1)
hold on
scatter(time_ss_200,angle_ss_200,'r.','DisplayName', RPM_Engine_200)
scatter(time_ss_300,angle_ss_300,'b.','DisplayName', RPM_Engine_300)
scatter(time_ss_500,angle_ss_500,'g.','DisplayName', RPM_Engine_500) 
legend
xlabel('Time [s]')
ylabel('Angle [deg]')
title('Engine - Angular Position | Steady State')

% get data from function
[time_tr_flywheels, angle_tr_flwheels] = data_fixer(data1.time_tr_flywheels,data1.tooth_count_tr_flwheels.*10);
[time_tr_no_flywheels, angle_tr_no_flwheels] = data_fixer(data1.time_tr_no_flywheels,data1.tooth_count_tr_no_flwheels.*10);
% find the average angular velocity for the transient data (this is for
% the future part to calculate the moment of inertia)
        % pull the set of values from 3-4.5 seconds
i1 = find(time_tr_flywheels>3,1);
i2 = find(time_tr_flywheels>4.5,1);
ss_angle_tr_flywheels = angle_tr_flwheels(i1:i2);
ss_time_tr_flywheels = time_tr_flywheels(i1:i2);
ss_angle_tr_no_flywheels = angle_tr_no_flwheels(i1:i2);
ss_time_tr_no_flywheels = time_tr_no_flywheels(i1:i2);
omega_b_fw_vec = [];
omega_b_nfw_vec = [];
        % loop through the values calcluating each angular velocity and
        % then take the average
for i = 1:length(ss_angle_tr_flywheels)-1
    omega_b_fw = (ss_angle_tr_flywheels(i+1) - ss_angle_tr_flywheels(i)) ...
        /(ss_time_tr_flywheels(i+1) - ss_time_tr_flywheels(i));
    omega_b_fw_rad_s = (pi/180)*omega_b_fw;
    omega_b_fw_vec = [omega_b_fw_vec; omega_b_fw_rad_s];
    omega_b_nfw = (ss_angle_tr_no_flywheels(i+1) - ss_angle_tr_no_flywheels(i)) ...
        /(ss_time_tr_no_flywheels(i+1) - ss_time_tr_no_flywheels(i));
    omega_b_nfw_rad_s = (pi/180)*omega_b_nfw;
    omega_b_nfw_vec = [omega_b_nfw_vec; omega_b_nfw_rad_s];
end
omega_b_fw_rad_s = mean(omega_b_fw_vec)
omega_b_nfw_rad_s = mean(omega_b_nfw_vec)

% plot transient data - Angle
figure (2)
hold on
scatter(time_tr_flywheels,angle_tr_flwheels,'r.','DisplayName', 'Transient Flywheels')
scatter(time_tr_no_flywheels,angle_tr_no_flwheels,'b.','DisplayName', 'Transient No Flywheels')
legend
xlabel('Time [s]')
ylabel('Angle [deg]')
title('Engine - Angular Position | Transient')
%% Load Saved Data for Torque
load('IC_Engines_Part2_2')
[torque_200, angle_200] = data_fixer_torque(data1.tooth_count_200.*10,torque(data.voltage_200));
[torque_300, angle_300] = data_fixer_torque(data1.tooth_count_300.*10,torque(data.voltage_300));
[torque_500, angle_500] = data_fixer_torque(data1.tooth_count_500.*10,torque(data.voltage_500));
torque_200 = -3.*torque_200;
torque_300 = -3.*torque_300;
torque_500 = -3.*torque_500;
torque_friction_200 = -mean(torque_200(1:120))
torque_friction_300 = -mean(torque_300(1:175))
torque_friction_500 = -mean(torque_500(1:250))
angle_200 = (pi/180)*angle_200;
angle_300 = (pi/180)*angle_300;
angle_500 = (pi/180)*angle_500;
% plot Steady State torque
figure (3)
hold on
scatter(angle_200,torque_200,'r.','DisplayName', RPM_Engine_200)
scatter(angle_300,torque_300,'b.','DisplayName', RPM_Engine_300)
scatter(angle_500,torque_500,'g.','DisplayName', RPM_Engine_500) 
legend
xlabel('Angle [rad]')
ylabel('Torque [N-m]')
title('Engine - Torque | Steady State')

% get data from function
[torque_tr_flywheels, angle_tr_flwheels] = data_fixer_torque(data1.tooth_count_tr_flwheels.*10,torque(data.voltage_tr_flwheels));
[torque_tr_no_flywheels, angle_tr_no_flwheels] = data_fixer_torque(data1.tooth_count_tr_no_flwheels.*10,torque(data.voltage_tr_no_flwheels));
torque_tr_flywheels = -3.*torque_tr_flywheels;
torque_tr_no_flywheels = -3.*torque_tr_no_flywheels;
angle_tr_flwheels = (pi/180)*angle_tr_flwheels;
angle_tr_no_flwheels = (pi/180)*angle_tr_no_flwheels;
% plot transient data - Angle
figure(4)
hold on
scatter(angle_tr_flwheels,torque_tr_flywheels,'r.','DisplayName', 'Transient Flywheels')
scatter(angle_tr_no_flwheels,torque_tr_no_flywheels,'b.','DisplayName', 'Transient No Flywheels')
legend
xlabel('Angle [rad]')
ylabel('Torque [N-m]')
title('Engine - Torque | Transient')
%% plot torque vs time
figure(5)
hold on
scatter(time_tr_flywheels,torque_tr_flywheels,'r.','DisplayName', 'Transient Flywheels')
scatter(time_tr_no_flywheels,torque_tr_no_flywheels,'b.','DisplayName', 'Transient No Flywheels')
legend
xlabel('Time [s]')
ylabel('Torque [N-m]')
title('Engine - Torque | Transient')
%% Find work
[work_tr_flywheels, angle_tr_flywheels_work] = data_fixer_work(torque_friction_500,data1.tooth_count_tr_flwheels.*10,3.*torque(data.voltage_tr_flwheels));
[work_tr_no_flywheels, angle_tr_no_flywheels_work] = data_fixer_work(torque_friction_500,data1.tooth_count_tr_no_flwheels.*10,3.*torque(data.voltage_tr_no_flwheels));
W_friction_fw = -0.8661*(angle_tr_flywheels_work(end)-angle_tr_flywheels_work(1))
W_friction_no_fw = -0.7856*(angle_tr_no_flywheels_work(end)-angle_tr_no_flywheels_work(1))
figure(6)
hold on
scatter(angle_tr_flywheels_work,work_tr_flywheels,'r.','DisplayName', 'Transient Flywheels')
scatter(angle_tr_no_flywheels_work,work_tr_no_flywheels,'b.','DisplayName', 'Transient No Flywheels')
legend
xlabel('Angle [rad]')
ylabel('Work [J]')
title('Cumulative Work')
work_no_flywheel = work_tr_no_flywheels(end)
work_flywheels = work_tr_flywheels(end)
omega_fly_b = 500; % RPM
omega_fly_b = omega_fly_b*pi/30;
I_fly = 2*work_flywheels/omega_b_fw_rad_s^2
% I_fly = 2*work_flywheels/omega_fly_b^2
I_fly_no = 2*work_no_flywheel/omega_b_nfw_rad_s^2
% I_fly_no = 2*work_no_flywheel/omega_fly_b^2
I_fly_net = I_fly - I_fly_no
%% theoretical Flywheel
I_fly_theoretical = 2*1.0141*0.07^2
%% Calculate piston displacement
tdc_m = .00065; 
bdc_m = .03837;
piston_m = @(x) .5*(bdc_m - tdc_m).*sind(x - 90) + .5*(bdc_m - tdc_m);
piston_distance_200rpm = piston_m(angle_ss_200);
piston_distance_300rpm = piston_m(angle_ss_300);
piston_distance_500rpm = piston_m(angle_ss_500);
figure(7)
hold on
scatter(time_ss_200,piston_distance_200rpm,'r.','DisplayName',RPM_Engine_200)
scatter(time_ss_300,piston_distance_300rpm,'b.','DisplayName',RPM_Engine_300)
scatter(time_ss_500,piston_distance_500rpm,'g.','DisplayName',RPM_Engine_500)
legend
xlabel('Time [s]')
ylabel('Piston Displacement [m]')
title('Piston Displacement')

% %% Plot data 
% figure(1)
% hold on
% scatter(data1.time_ss_200,data1.tooth_count_200.*10,'b.','DisplayName','Steady State 200 RPM')
% scatter(data1.time_ss_300,data1.tooth_count_300.*10,'r.','DisplayName','Steady State 300 RPM')
% scatter(data1.time_ss_500,data1.tooth_count_500.*10,'g.','DisplayName','Steady State 500 RPM')
% legend
% figure(2)
% hold on
% scatter(data1.time_tr_flywheels,data1.tooth_count_tr_flwheels.*10,'b.','DisplayName','Transient All Flywheels')
% scatter(data1.time_tr_no_flywheels,data1.tooth_count_tr_no_flwheels.*10,'r.','DisplayName','Transient No Flywheels')
% xlabel('Time [s]')
% ylabel('Angle [deg]')
% legend