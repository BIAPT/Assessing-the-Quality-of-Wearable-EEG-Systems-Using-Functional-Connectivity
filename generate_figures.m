%{
    Yacine Mahdid March 19 2020 (COVID-19 Outbreak home time)
    In this script we are processing the data in the /result folder
    in order to generate the figure we need for the paper.
    This make use of NeuroAlgo

    Figure we need:
    for all headset
    1) Spectrogram with temporal smoothing of order 10
    2) Topographic Maps
    3) wPLI reordered FTCPO
    4) dPLI reordered FTCPO
    
    for the downsampled headset
    5) wPLI reordered FTCPO
    6) dPLI reordered FTCPO
%}

%% General Parameters
DATA_DIR = "/home/yacine/Documents/BIAPT/Assessing-the-Quality-of-Wearable-EEG-Systems-Using-Functional-Connectivity/result/";
OUT_DIR = "/home/yacine/Documents/BIAPT/Assessing-the-Quality-of-Wearable-EEG-Systems-Using-Functional-Connectivity/figure/";

STATES = {'closed','open'};
HEADSETS = {'egi','quick_30','dsi','open_bci','emotiv'};
FREQUENCIES = {'alpha', 'beta', 'theta', 'delta'}; 

% These will be transformed into loops when ready
for f = 1:length(FREQUENCIES)
    frequency = FREQUENCIES{f};
    for h = 1:length(HEADSETS)
        headset = HEADSETS{h};
        
        % We are doing both the open and closed state at the same time to
        % be able to normalize correctly across the matrices

        filename = strcat(DATA_DIR,frequency, filesep, headset, "_", "open", ".mat");
        result_open = load(filename);

        filename = strcat(DATA_DIR,frequency, filesep, headset, "_", "closed", ".mat");
        result_close = load(filename);
        
         %% Weighted Phase Lag Index
        handle_wpli = make_figure_pli(result_open, result_close, headset, frequency, "wpli");
        filename = strcat(OUT_DIR, frequency, filesep, headset, "_wpli.png");
        saveas(handle_wpli,filename);
        handle_dpli = make_figure_pli(result_open, result_close, headset, frequency, "dpli");
        filename = strcat(OUT_DIR, frequency, filesep, headset, "_dpli.png");
        saveas(handle_dpli,filename);
        close all;
        
            
    end

end

function [handle] = make_figure_pli(result_open, result_close, headset, frequency, type)
        switch type
            case "wpli"
                pli_open = result_open.wpli;
                pli_close = result_close.wpli;
                
                if ~strcmp(headset, "egi")
                    pli_down_open = result_open.down_sampled_wpli;
                    pli_down_close = result_close.down_sampled_wpli;
                end
            case "dpli"
                pli_open = result_open.dpli;
                pli_close = result_close.dpli;
                if ~strcmp(headset, "egi")
                    pli_down_open = result_open.down_sampled_dpli;
                    pli_down_close = result_close.down_sampled_dpli;
                end
        end

        [p_pli_open, ~, p_regions_open] = process_pli(pli_open, headset, type);
        [p_pli_close, ~, p_regions_close] = process_pli(pli_close, headset, type);


        if strcmp(headset, "egi")
            p_pli_all = [p_pli_open(:); p_pli_close(:)];
            handle = figure;
            subplot(2,1,1)
            plot_pli(p_pli_close,p_regions_close,p_pli_all)
            title(strcat(headset, " for close at ", frequency, " for ", type));
            colorbar;
            subplot(2,1,2)
            plot_pli(p_pli_open,p_regions_open,p_pli_all)
            title(strcat(headset, " for open at ", frequency, " for ", type));
            colorbar;
            set(handle, 'Position', [976, 1, 509, 785]); % this is to have normalized size across figures
        
        else
            [p_pli_down_open, ~, p_regions_down_open] = process_pli(pli_down_open, "egi", type);
            [p_pli_down_close, ~, p_regions_down_close] = process_pli(pli_down_close, "egi", type);            
            
            p_pli_all = [p_pli_open(:); p_pli_close(:); p_pli_down_open(:); p_pli_down_close(:)];
            
            handle = figure;
            subplot(2,2,1)
            plot_pli(p_pli_down_close, p_regions_down_close, p_pli_all)
            title(strcat("egi-", headset," downsampled for close at ", frequency, " for ", type));
            subplot(2,2,2)
            plot_pli(p_pli_close, p_regions_close, p_pli_all)
            title(strcat(headset, " for close at ", frequency, " for ", type));
            colorbar;
            subplot(2,2,3)
            plot_pli(p_pli_down_open, p_regions_down_open, p_pli_all)
            title(strcat("egi-", headset," downsampled for open at ", frequency, " for ", type));
            subplot(2,2,4)
            plot_pli(p_pli_open, p_regions_open, p_pli_all)
            title(strcat(headset, " for open at ", frequency, " for ", type));
            colorbar;
            set(handle, 'Position', [587, 55, 976, 686]);
        end
end


function plot_pli(pli,regions,pli_all)
    imagesc(pli);
    xtickangle(90)
    xticklabels(regions);
    yticklabels(regions);  
    xticks(1:length(regions));
    yticks(1:length(regions));
    min_color = min(pli_all);
    max_color = max(pli_all);
    caxis([min_color max_color])
    colormap("jet");    
end

function [p_pli, p_location, p_regions] = process_pli(pli_struct, headset, type)
        % Select the casset that will be used for the reordering /
        % filtering of non-scalp electrodes
        switch headset
            case "egi"
                cassette = "biapt_egi129.csv";
            case "quick_30"
                cassette = "biapt_quick_30.csv";
            case "emotiv"
                cassette = "biapt_quick_30.csv";
            case "dsi"
                cassette = "biapt_dsi24.csv";
            case "open_bci"
                cassette = "biapt_dsi24.csv";
            otherwise
                disp("Error wrong headset name");
                return
        end
        
        if strcmp(type, "wpli")
            avg_pli = pli_struct.data.avg_wpli;            
        else
            avg_pli = pli_struct.data.avg_dpli;
        end

        location = pli_struct.metadata.channels_location;

        % Filtering step where we take only the left handside
        [avg_pli, location] = get_left_hemisphere(avg_pli, location);

        [p_pli, p_location, p_regions] = reorder_channels(avg_pli, location, cassette);
end


function [f_pli, f_location] = get_left_hemisphere(pli, location)

    % Make a mask for the left side
    is_left = ([location.is_left] == 1);
    f_pli = pli(is_left, is_left);

    % Iterate over the location and get the left side of the brain
    f_location = [];
    for i = 1:length(location)
        if location(i).is_left
            f_location = [f_location, location(i)]; 
        end
    end
end










function power_map_figure = topographic_map(power,location)
    power_map_figure = figure;
    topoplot(power, location,'maplimits','absmax', 'electrodes', 'off');
    min_color = min(power);
    max_color = max(power);
    caxis([min_color max_color])
    colorbar;
end

% Not used since we do not need to reproduce this part
function make_spectrogram(sp)
    spectrum = sp.data.spectrums;
    timestamps = sp.data.timestamps;
    frequencies = sp.data.frequencies;

    % From EEGapp
    smooth_spectrum = medfilt1(spectrum, 10);  % This performs temporal smoothing with a median filter of order tso7
    %Create the figure for the spectrogram
    figure;
    plot_matrix(smooth_spectrum, timestamps, frequencies);
    title(' Spectrogram(Average of all electrodes)');
    ylabel('Frequency','fontsize',12);
    xlabel('Time','fontsize',12);
    colormap(jet);
end