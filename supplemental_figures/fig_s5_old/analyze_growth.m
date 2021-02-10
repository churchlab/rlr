function analyze_growthMax(filename, opt_blank_wells, opt_blank_value, ...
        opt_interval, opt_mid_log_interval, opt_hide_plots, ...
        opt_output_override, opt_name)
% ANALYZE_GROWTH Compute doubling time given kinetic read time series.
%
%     Args:
%         filename: Full path to kinetic read data. Tab-delimited. First row is
%             well names. Each row is the the value of reads at each time point.
%         opt_blank_wells: Optional array of wells left blank for calibration.
%             Provide as integers. E.g. [48, 96] for wells D12 and H12.
%         opt_blank_value: Optional value of blank read. If provided,
%             opt_blank_wells is ignored.
%         opt_interval: Optional. Interval, in minutes, between reads.
%             Default 5 min.
%         opt_mid_log_interval: Optional. Size of window, in minutes, for which
%             we measure linear growth. Default 40 minutes.
%         opt_hide_plots: Optional. Boolean. If true, hide plots.
%         opt_output_override: Optional. Override the name of the output file
%           opt_name (MAX) name used to title plots
%
%     The output is written to a new file in the same location as the input
%     a text file, with extension '.analyzed_growth.csv'. This can be imported
%     into Excel as a tab-delimited file.
%
%     Example usage:
%
%         analyze_growth('/home/glebk/Data/2015_06_10_growth_test.txt')
%
%     Blank wells may be provided to be sourced as the average blank read.
%     For example, if well H12 (96th well) is blank:
%
%         analyze_growth('/home/glebk/Data/2015_06_10_growth_test.txt', [96])
%
%     Alternatively, one can provide a global blank value to adjust all reads by:
%
%         analyze_growth('/home/glebk/Data/2015_06_10_growth_test.txt', [], 0.09)
%
%     Note for using optional arguments: User must provide values for all
%     arguments up to the optional argument you would like to use.
%     For example, to set opt_mid_log_interval to 30 min, the command is:
%
%         analyze_growth('/home/glebk/Data/2015_06_10_growth_test.txt', [], 0, 5, 30)
%
%     Authors:
%         Jaron Mercer - ported to Matlab
%         Harris Wang - original draft
%         Gleb Kuznetsov (kuznetsov@g.harvard.edu) - primary maintainer
%         Tim Wannier


% Close any open windows.
close all;

% Default interval, in minutes, between measurements in the input data.
DEFAULT_INTERVAL = 5;

% Default window, in minutes, when we measure log-linear growth.
DEFAULT_MID_LOG_INTERVAL = 40;


% Initial data import.
input_data = importdata(filename, '\t', 1);

% Matrix where rows are consecutive time measurements and each column
% corresponds to a well.
data = input_data.data;

% Well names.
headers = input_data.colheaders;

% Size of data.
num_points = size(data, 1);
num_wells = size(data, 2);


% Determine blank measurement. Default 0.
%
% Warn user if they have not specified blank wells, as this will lead to
% incorrect absolute doubling times.
%
% Blank reads are removed from data.

% Default blank value.
blank_avg = 0.0;

if exist('opt_blank_value', 'var') & opt_blank_value > 0
    blank_avg = opt_blank_value;
elseif (exist('opt_blank_wells', 'var') & length(opt_blank_wells) > 0 & ...
        opt_blank_wells(1) > 0)
    % The number of blanks.
    num_blanks = size(opt_blank_wells, 2);

    % Sort blanks so that data can be removed from matrix.
    opt_blank_wells = sort(opt_blank_wells, 'descend');

    % Initiate blanks matrix.
    blank_reads = zeros(num_points,num_blanks);

    % Copy blanks to blanks matrix and remove from data and headers.
    for well = 1:num_blanks
        blank_reads(:, well) = data(:, opt_blank_wells(well));
        data(:, opt_blank_wells(well)) = [];
        headers(:, opt_blank_wells(well)) = [];
    end

    % Columns of data with no blanks.
    num_wells = size(data, 2);

    % Average blanks matrix.
    blank_avg = mean(mean(blank_reads));
else
    prompt = strcat(...
        '\n   -------------\n   -- WARNING --\n   -------------\n\n' ...
        , 'No blank wells have been specified.\n' ...
        , 'Running a growth rate analysis on unblanked data will return erroneous results.\n');
    fprintf(prompt);
end

% Subtract blank (might be 0, in which case this is NO-OP).
data = data - blank_avg;

% Make sure values are positive following blanking.
data = max(data, 0.001);


% Set interval between reads.
interval = DEFAULT_INTERVAL;
if exist('opt_interval', 'var')
    interval = opt_interval;
end


% Set delta (number of timepoints) based on mid-log window and interval
% between measurements.
mid_log_interval = DEFAULT_MID_LOG_INTERVAL;
if exist('opt_mid_log_interval', 'var')
    mid_log_interval = opt_mid_log_interval;
end
delta = round(mid_log_interval / interval);


% Determine whether to show plots or not.
show_plots = true;
if exist('opt_hide_plots', 'var') & opt_hide_plots
    show_plots = false;
end


% Show data dimensions.
whos data;

% Plot data.
if show_plots
    plot(data);
    title(strcat(headers{1,1}, ' - ', headers{1,size(headers,2)}));
    xlabel('Time (min/5)');
    ylabel('OD 600');
    title(opt_name);
end


% We fit the natural log of the data.
ln_data = log(data);
ln_data(isinf(ln_data)) = NaN;
ln_data = real(ln_data);

% Plot log of data.
if show_plots
    figure(2);
    plot(ln_data);
    xlabel('Time (min/5)');
    ylabel('ln(OD 600)');
    title(strcat(opt_name,'ln(', headers{1,1}, ' - ', headers{1,size(headers,2)}, ')'));
end


%%% Find the linear portion of ln(OD 600). Then do linear regression.
%
% We dynamically search for the limits of the linear portion of the log-OD
% plot. We do this by performing a linear regression on a sliding interval
% window, and keeping track of the window with the greatest slope.
%
% This is the main "algorithmic" part of the script and the primary change
% Jaron Mercer made to Harris Wang's version of this script.

mus = zeros(size(ln_data,2),1);
doubleTs = zeros(size(ln_data,2),1);
rSqrs = zeros(size(ln_data,2),1);
maxODs = zeros(size(ln_data,2),1);
deltas = zeros(size(ln_data,2),1);
window_starts = zeros(size(ln_data,2),1);
maxSlopeODs = zeros(size(ln_data,2),1);
warnings = zeros(size(ln_data,2),1);

for well = 1:num_wells
    % Initialize state variables.
    maxSlope = 0;
    maxRsquared = 0;
    maxDelta = 0;
    maxStart = 0;
    maxSlopeOD = 0;

    % Find the greatest slope.
    intervals_since_greatest = 0;

    for window_start = 1:(size(ln_data, 1) - delta)
        % We expect the greatest interval early on so no need to go more
        % than 2 hours past max.
        % Only break if the correlation is good, else search exhausively.
        %  -- Condition added by Tim W.
        if intervals_since_greatest > 24 & maxRsquared > .995
            break
        end

        % Only bother fitting if the interval OD is between 0.2 and 0.7.
        if data(window_start, well) > 0.2 & data(window_start, well) < 0.7
            x = (linspace(window_start, window_start + delta - 1, delta))';
            y = ln_data(window_start: window_start + delta - 1, well);

            % returns 1x2 matrix: [slope, y-intercept]
            line = polyfit(x,y,1);

            corrcoefMatrix = corrcoef(x,y);
            % Use either of the off-axis values to calculate R-squared.
            rSquared = corrcoefMatrix(1, 2) ^ 2;
            if line(1, 1) > maxSlope & rSquared > 0.95
                maxSlope = line(1, 1);
                maxRsquared = rSquared;
                maxDelta = delta;
                maxStart = window_start;
                maxSlopeOD = data(window_start + delta / 2, well);
                intervals_since_greatest = 0;
            else
                intervals_since_greatest = intervals_since_greatest + 1;
            end
        end
    end

    % Save output data.
    % To prevent erroring on blank wells that weren't specified or wells with
    % insufficient maxODs, we will zero out these wells. -- Tim W.
    if maxSlope == 0
        mus(well, 1) = 0;
        doubleTs(well, 1) = 0;
        rSqrs(well, 1) = 0;
        maxODs(well, 1) = 0;
        window_starts(well, 1) = 0;
        deltas(well, 1) = 0;
        maxSlopeODs(well, 1) = 0;
        % display warning
        disp(strcat('Warning: no signal on well --', headers(1, well)));
        warnings(well, 1) = 1;
    else
        % Added µ to the list of outputs -- Tim W.
        mus(well, 1) = (maxSlope / interval) * 60;
        doubleTs(well, 1) = (log(2) / maxSlope) * interval;
        rSqrs(well, 1) = maxRsquared; % save r-squared
        maxODs(well, 1) = max(data(:, well));
        window_starts(well, 1) = maxStart * interval;
        deltas(well, 1) = maxDelta;
        maxSlopeODs(well, 1) = maxSlopeOD;
    end

    % Output a warning for low r-squared values.
    if rSqrs(well, 1) < .95
        disp(strcat('Warning: low confidence on well --', headers(1, well)));
        if show_plots
            figure(3); hold  on; plot(ln_data(:, well));
        end
        warnings(well, 1) = 1;
    end
end


%%% Save data to csv.

if exist('opt_output_override', 'var')
    output_filename = opt_output_override;
else
    output_filename = strcat(filename(1:size(filename, 2) - 4), '.analyzed_growth.csv');
end

fid = fopen(output_filename, 'w');

% Write the header row.
fprintf(fid, 'id,well,mu_hourly,doubling_time_min,r_sqrd,max_OD,window_start_time_min,delta,max_slope_od,warnings\n');

% Iterate through the well data, writing one row at a time.
for well = 1:num_wells
    fprintf(fid, '%d,%s,%f,%f,%f,%f,%f,%f,%f,%f\n', ...
            well, headers{well}, mus(well), doubleTs(well), rSqrs(well), ...
            maxODs(well), window_starts(well), deltas(well), ...
            maxSlopeODs(well), warnings(well));
end

fclose(fid);
