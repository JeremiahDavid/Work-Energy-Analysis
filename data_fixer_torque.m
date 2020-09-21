function [torque_uniq, angle_uniq_torque] = data_fixer_torque(angle,torque)
%% Adjust Data for Continuous Plot
angle_cont = angle;
identifier = diff(angle_cont); % the difference between i and i+1 in angle_count
flags = find(identifier<0) + 1; % find the indexes where the angle resets to 10
flags = [flags; length(angle_cont)]; % add the last index to the list of flags
for i = 1:length(flags)-1 % loop through the indexes that were flagged
    previous_max = angle_cont(flags(i)-1); % the index right before it was flagged
    if flags(i) ~= flags(end-1) % test if the current flag index = the second to last flag index
        angle_cont(flags(i):flags(i+1)-1) = angle_cont(flags(i):flags(i+1)-1)+previous_max;
        % ^ add the previous max to the section of angle values between the
        % current flag index and the next flag index (includes the current
        % flag index and excludes the next flag index)
    else % the last section of angle values in the data set
        angle_cont(flags(i):flags(i+1)) = angle_cont(flags(i):flags(i+1))+previous_max;
        % ^ add the previous max to the section of angle values between the
        % current flag index and the next flag index (includes both the
        % current and next flag indexes)
    end
end

%% Average Repeating Angles - Torque
[torque_uniq,angle_uniq_torque] = Averager(torque,angle_cont);
% angle_uniq_torque_rad = angle_uniq_torque.*(pi/180); % convert the angle to radians
% %% Numerical Integration
% integration_values = [0]; % initialize the vector to hold the integratoin values
% angle_uniq_torque_rad = [0; angle_uniq_torque_rad]; % add 0 to the beginning
% torque_uniq = [0; torque_uniq]; % add 0 to the beginning
% for i = 1:length(angle_uniq_torque_rad)-1
%     x1 = angle_uniq_torque_rad(i); % beginning of integration for segment
%     x2 = angle_uniq_torque_rad(i+1); % end of integration for segment
%     x_dif = x2-x1; % integration spacing
%     y1 = -torque_uniq(i); % beginning of y
%     y2 = -torque_uniq(i+1); % end of y
%     integration = trapz(x_dif,[y1 y2]); % perform the integration
%     integration_values = [integration_values; integration_values(end)+integration] % add the current integration the vector
% end
% % angle_uniq_torque_rad(1) = []; % remove 0 from the beginning


    function [normal_data_uniq,stepped_data_uniq] = Averager(normal_data,stepped_data)
        % this function will take a data set that has repeating values in the data
        % over an increasing time, and represent the repeating values with the last
        % value that occurs. It will also record the time this occured, and all
        % this new data is stored in 2 new vectors.
        stepped_data_uniq = []; % initialize the new vector of unique data values
        normal_data_uniq = []; % initialize the new vector of time values for the new unique data values
        average_normal = 0;
        start = 1;
        for i = 1:length(stepped_data)-1
            if stepped_data(i) == stepped_data(end) % test for the end case
                stepped_data_uniq = [stepped_data_uniq; stepped_data(end)]; % grab the last value in continuous data set
                normal_data_uniq = [normal_data_uniq; average_normal]; % time value at this value
            elseif stepped_data(i) == stepped_data(i+1) % test if the last value was the same
                average_normal = mean(normal_data(start:i));
            else % if the last value was not the same and it is not the last value in the data set
                start = i+1; % assign a new start value for the average calculater
                stepped_data_uniq = [stepped_data_uniq; stepped_data(i)]; % add the unique value to the vector
                normal_data_uniq = [normal_data_uniq; average_normal]; % time value at this value
            end
        end
    end
end
% %% Average Repeating Angles - Torque
% [torque_uniq,angle_uniq_torque] = Averager(torque,angle_cont);
% torque_uniq = torque_uniq*-1; % correct the direction of positive torque
% angle_uniq_torque_rad = angle_uniq_torque.*(pi/180); % convert the angle to radians
% %% Plot Data with Repeating Torques Averaged
% figure(4)
% scatter(angle_uniq_torque_rad,torque_uniq,'r.')
% xlabel('Angular Position [rad]')
% ylabel('Torque [N-m]')
% title('All Data - Averaged Reapeating Angle - Torque')
% %% Numerical Integration
% integration_values = [0]; % initialize the vector to hold the integratoin values
% angle_uniq_torque_rad = [0; angle_uniq_torque_rad]; % add 0 to the beginning
% torque_uniq = [0; torque_uniq]; % add 0 to the beginning
% for i = 1:length(angle_uniq_torque_rad)-1
%     x1 = angle_uniq_torque_rad(i); % beginning of integration for segment
%     x2 = angle_uniq_torque_rad(i+1); % end of integration for segment
%     x_dif = x2-x1; % integration spacing
%     y1 = -torque_uniq(i); % beginning of y
%     y2 = -torque_uniq(i+1); % end of y
%     integration = trapz(x_dif,[y1 y2]); % perform the integration
%     integration_values = [integration_values; integration_values(end)+integration] % add the current integration the vector
% end
% % angle_uniq_torque_rad(1) = []; % remove 0 from the beginning
% %% Plot Data with Repeating Torques Averaged
% figure(5)
% scatter(angle_uniq_torque_rad,integration_values,'r.')
% xlabel('Angular Position [rad]')
% ylabel('Work [Joules]')
% title('All Data - Averaged Reapeating Angle - Work')



