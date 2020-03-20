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
        for s = 1:length(STATES)
            state = STATES{s};

            filename = strcat(DATA_DIR,frequency, filesep, headset, "_", state, ".mat");
            result = load(filename);

            %% Spectrogram
            sp = result.sp;

            %% Power Topographic Map
            %{
            td = result.td;
            power = td.data.power;
            location = td.metadata.channels_location;

            if strcmp(headset, "egi")
                [power, location] = filter_channels(power, location);    
            end
            power_map_figure = topographic_map(power,location);
            title(strcat(headset, " for ", state, " at ", frequency));
            filename = strcat(OUT_DIR,frequency, filesep, headset, "_", state, "_topo.fig");
            saveas(power_map_figure,filename,'fig')
            close all; % to clear the figures that pops up
            %}
            
            %% Weighted Phase Lag Index
            wpli = result.wpli;
            avg_wpli = wpli.data.avg_wpli;
            location = wpli.metadata.channels_location;
            
            if strcmp(headset, "egi")
                [r_wpli, r_location, r_regions] = reorder_channels(avg_wpli, location);
                figure;
                imagesc(r_wpli);
                %xtickangle(90)
                %xticklabels(r_regions);
                %yticklabels(r_regions);  
                %xticks(1:length(r_regions));
                %yticks(1:length(r_regions));
                colormap("jet");
                colorbar;
            end
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