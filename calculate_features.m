%{
    Yacine Mahdid March 18 2020 (COVID-19 Outbreak home time)
    In this script we are processing the data in the /data folder
    to yield the result we need for the experiment.

    This make use of NeuroAlgo
%}

%% General Parameters
DATA_DIR = "/home/yacine/Documents/Assessing-the-Quality-of-Wearable-EEG-Systems-Using-Functional-Connectivity/data/";
OUT_DIR = "/home/yacine/Documents/Assessing-the-Quality-of-Wearable-EEG-Systems-Using-Functional-Connectivity/result/";

STATES = {'closed','open'};
HEADSETS = {'egi','quick_30','dsi','open_bci','emotiv'};
FREQUENCIES = {'alpha', 'beta', 'theta', 'delta'}; 

%% Experimental Parameters
% Spectral Power
SP = struct();
SP.TBP = 2;
SP.NUM_TAPER = 3;
SP.SPR_WINDOW_SIZE = 2; % in seconds
SP.STEP_SIZE = 0.1; % in seconds
SP.BP = [1 30];

% Topographic Distribution
TD = struct();
TD.BP = [8 13]; % in Hz

% wPLI and dPLI
PLI = struct();
PLI.BP = [8 13]; % This is in Hz
PLI.WINDOW_SIZE = 10; % This is in seconds and will be how we chunk the whole dataset
PLI.NUMBER_SURROGATE = 20; % Number of surrogate wPLI to create
PLI.P_VALUE = 0.05; % the p value to make our test on
PLI.STEP_SIZE = PLI.WINDOW_SIZE;


for f = 1:length(FREQUENCIES)
    frequency = FREQUENCIES{f};
    bandpass = [];
    switch frequency
        case "alpha"
            bandpass = [8 13];
        case "beta"
            bandpass = [14 30];
        case "theta"
            bandpass = [4 7];
        case "delta"
            bandpass = [1 4];
        otherwise 
            disp("Wrong frequency!")
            return
    end
    for h = 1:length(HEADSETS)
        headset = HEADSETS{h};
        for s = 1:length(STATES)
            state = STATES{s};

            out_filename = strcat(OUT_DIR, frequency,filesep, headset,"_",state,".mat");
            result = struct();

            disp(strcat("Analyzing ", headset , " At state: ", state));
            filename = strcat(DATA_DIR, headset, filesep, headset, "_", state, ".mat");
            recording = load_mat(filename);
        
            % Spectral Power

            SP.window_size = floor(recording.length_recording/recording.sampling_rate); %dynamic parameter

            result.sp = na_spectral_power(recording, SP.window_size, SP.TBP, ... 
                                SP.NUM_TAPER, SP.SPR_WINDOW_SIZE, SP.BP, SP.STEP_SIZE);

            % Power Topographic Map
            TD.window_size = SP.window_size; % dynamic parameter
            TD.step_size = SP.window_size;
            TD.BP = bandpass;
            result.td = na_topographic_distribution(recording, TD.window_size, ...
                                                    TD.step_size, TD.BP);
                                                
            % Functional Connectivity
            PLI.BP = bandpass;
            result.wpli = na_wpli(recording, PLI.BP, PLI.WINDOW_SIZE, PLI.STEP_SIZE, ...
                                  PLI.NUMBER_SURROGATE, PLI.P_VALUE);

            result.dpli = na_dpli(recording, PLI.BP, PLI.WINDOW_SIZE, PLI.STEP_SIZE, ...
                                  PLI.NUMBER_SURROGATE, PLI.P_VALUE);


            % If we are not working with the egi headset we should also load
            % the downsampled egi headset
            if strcmp(headset, "egi") == 0
               disp(strcat("Calculating functional connectivity on the downsampled ", headset));
               filename = strcat(DATA_DIR, headset, filesep, "egix",headset, "_", state, ".mat");
               recording = load_mat(filename);

               result.down_sampled_wpli = na_wpli(recording, PLI.BP, PLI.WINDOW_SIZE, PLI.STEP_SIZE, ...
                                  PLI.NUMBER_SURROGATE, PLI.P_VALUE);
               result.down_sampled_dpli = na_dpli(recording, PLI.BP, PLI.WINDOW_SIZE, PLI.STEP_SIZE, ...
                                  PLI.NUMBER_SURROGATE, PLI.P_VALUE);
            end


            % Save the data
            save(out_filename, '-struct', 'result');
        
        end
    end
end

% Note this pipeline is not perfect with respect to step size
% for spectral power, will need to average out some of the residual window