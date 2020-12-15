source = ["fMRI" "FRED" "weather"];
lambda_list = [1.5 0.5 0.75]; % 0.04 for hourly, 0.75 for daily
beta_list = [1e3 1e3 1e3];
rho_list = [1 1 1];
norm_type = 6;

refit = true;
data_idx = 1;
fMRI_idx = 1;
fMRI_num_nodes = 15;
weather_daily = true;

disp(['Data from: ', source(data_idx)]);

if data_idx == 1
    load('NV_COS_data.mat');
    data_cell = NVdata.timeseries(fMRI_idx);
    data = data_cell{1};
    [n, ~] = size(data);
    rand_idx = randperm(n);
    data = data(rand_idx(1:fMRI_num_nodes), :);
elseif data_idx == 2
    load Data_USEconModel;
    [row, col] = find(isnan(Data));
    cutoff = max(row);
    Data = Data(cutoff+1:end,:);
    data = Data.';
else
    if weather_daily
        file = importdata('weather_daily.csv'); % Temperature	Humidity	Pressure	Wind Speed	Maximum Temperature	Minimum Temperature	Peak Wind Speed	Precipitation	Snow Depth	Snowfall	Sustained Wind Speed
    else
        file = importdata('weather_hourly.csv'); % Temperature	Precipitation	Humidity	Pressure	Visibility	Wind Gust Speed	Wind Speed
    end
    text_data = string(file.textdata);
    headers = text_data(1,:);
    data = file.data;
    data = data';
    dates = text_data(2:end,1);
    dates = datetime(dates,'InputFormat','yyyy-MM-dd''T''HH:mm:ss','TimeZone','local');
    dates = fillmissing(dates,'linear');
    [~, unique_ind] = unique(dates);
    duplicate_ind = setdiff(1:length(dates), unique_ind);
    dates(duplicate_ind) = [];
    data(:, duplicate_ind) = [];
    data = fillmissing(data,'linear',2,'SamplePoints',dates,'EndValues','nearest');
end

data = data-mean(data,2);
data = shiftdim(data,-1);
[~, n, T] = size(data);
lambda = lambda_list(data_idx);
beta = beta_list(data_idx);
rho = rho_list(data_idx);

if refit
    start = tic;
    [Thetas, ~] = tvgl_self(data, lambda, beta, rho, norm_type);
    run_time = toc(start);
end

time_idx = round(1:T/9:T);
figure;
for i = 1:9
    time = time_idx(i);
    G = graph(Thetas(:,:,time), 'upper', 'OmitSelfLoops');
    LWidths = 5*G.Edges.Weight / max(G.Edges.Weight);
    subplot(3, 3, i);
    if data_idx == 1
        NLabels = rand_idx(1:fMRI_num_nodes);
    elseif data_idx == 2
        NLabels = DataTable.Properties.VariableNames;
    else
        NLabels = cellstr(headers(2:end));
    end
    if isempty(LWidths)
        plot(G, 'NodeLabel', NLabels, 'MarkerSize', 7, 'Layout', 'circle');
    else
        plot(G, 'LineWidth', LWidths, 'NodeLabel', NLabels, 'MarkerSize', 7, 'Layout', 'circle');
    end
    if data_idx == 1
        title(['T = ', num2str(time)]);
    elseif data_idx == 2
        formatOut = 'QQ-yy';
        quarters = DataTable.Time(cutoff+1:end);
        title(['T = ', datestr(quarters(time), formatOut)]);
    else
        if weather_daily
            formatOut = 'mm/dd';
        else
            formatOut = 'mm/dd HH:MM AM';
        end
        title(['T = ', datestr(dates(time), formatOut)]);
    end
    hold on;
end
hold off;

td = zeros(T-1,1);
for i = 2:T
    td(i-1) = norm(Thetas(:,:,i)-Thetas(:,:,i-1), 'fro');
end
figure;
if data_idx == 1
    semilogy(2:T, td);
elseif data_idx == 2
    semilogy(DataTable.Time(cutoff+2:end), td);
else
    semilogy(dates(2:end), td);
end
xlabel('Time');
ylabel('Temporal Deviation');
hold off;

figure;
[max_td, max_td_idx] = max(td(2:end));
max_change_time = max_td_idx-1:max_td_idx;
max_change_time = max_change_time + 1;
for i = 1:2
    t = max_change_time(i);
    G = graph(Thetas(:,:,t+1), 'upper', 'OmitSelfLoops');
    LWidths = 5*G.Edges.Weight / max(G.Edges.Weight);
    subplot(1, 2, i);
    if data_idx == 1
        NLabels = rand_idx(1:fMRI_num_nodes);
    elseif data_idx == 2
        NLabels = DataTable.Properties.VariableNames;
    else
        NLabels = cellstr(headers(2:end));
    end
    if isempty(LWidths)
        plot(G, 'NodeLabel', NLabels, 'MarkerSize', 7, 'Layout', 'circle');
    else
        plot(G, 'LineWidth', LWidths, 'NodeLabel', NLabels, 'MarkerSize', 7, 'Layout', 'circle');
    end
    if data_idx == 1
        title(['T = ', num2str(t+1)]);
    elseif data_idx == 2
        formatOut = 'QQ-yy';
        quarters = DataTable.Time(cutoff+1:end);
        title(['T = ', datestr(quarters(t+1), formatOut)]);
    else
        if weather_daily
            formatOut = 'mm/dd';
        else
            formatOut = 'mm/dd HH:MM AM';
        end
        title(['T = ', datestr(dates(t+1), formatOut)]);
    end
    hold on;
end
hold off;