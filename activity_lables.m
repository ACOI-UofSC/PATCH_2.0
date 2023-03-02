% First read in the activity log file
logfile_name = strcat(participant_num, "_time log.xlsx");
opts = detectImportOptions(logfile_name); % Automatically detects the format of file
% Read in each column as a separate VAR. This is ideal because of the mixed
% data types
[labels, acti_start, acti_end, flag_start,flag_end, flag, ppg_begin, ppg_end] = readvars(logfile_name);
% t_PPG begins at the beginning of the trial this corresponds to the time of the first activity
patch_start = acti_start(1); 
acti_times = cat(2, acti_start, acti_end);
flag_times = cat(2, flag_start, flag_end);
% Subtract the patch start from acti_start
acti_times = acti_times - patch_start;
% Intialize the number corresponding to 1 second
base_second = datenum(duration(0,0,1));
% divide acti_start's number by base_sec
% base_second to determine how many seconds
% have passed since the start of the PPG sensor.
acti_times =round(acti_times / base_second);

% Do the same process for the flag_times if it isn't empty
if iscell(flag_times) == 0
    flag_times = flag_times - patch_start;
    flag_times = round(flag_times / base_second);
end
