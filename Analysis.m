%% Analysis

% Script to analyze all animals
% - The code relies on the following functions : -> Open_Files.m
%                                                -> Extracting_LFPs_and_events_from_all.m
%                                                -> ... other desired functions
% by Flavio Mourao.
% email: mourao.fg@illinois.edu
% Maren Lab - Beckman Institute for Advanced Science and Technology
% University of Illinois Urbana-Champaign

% Started in:  12/2023
% Last update: 10/2024

%%

tic
set(groot,'DefaultFigureColormap',jet)

[files] = Open_files();

for ms = 1:length(files.FilesLoaded{1, 1}) % Loop over files path

    % Settings
    id          = files.id(ms).name;
    FilesLoaded = {files.FilesLoaded{1,1}(ms).name};
    Path        = files.FilesLoaded{1,1}(ms).folder;

    %% Open and extract raw files

    % Load function
    % [data,parameters] = Extracting_LFPs_and_events(id,FilesLoaded,Path);

    %% Scripts pre-processing

    %    Pre_processing_;        % Pre processing
    %    Behavior                % Behavior Analyse
    %    Behavior_exposure       % Behavior Analyse forexposure sessions

    %% Load. Choose data from data.lfp according pre_processing.m and main plots script
    name = strcat(Path,'/',id);
    fprintf('\n loading data... \n');

    %     load(name,'data','parameters');
    load(name);

    %% Context Exposure

    % 1) Baseline timeepochs -> first minute. All animals checked. No noise
    % B_clean{1} = dsearchn(data.timev_decimated',[10 170]');      % MT1
    % B_clean{2} = dsearchn(data.timev_decimated',[10 170]');      % MT3
    % B_clean{3} = dsearchn(data.timev_decimated',[10 170]');      % MT4
    % B_clean{4} = dsearchn(data.timev_decimated',[10 170]');      % MT5
    % B_clean{5} = dsearchn(data.timev_decimated',[10 170]');      % MT6
    % B_clean{6} = dsearchn(data.timev_decimated',[10 170]');      % MT7

    % B_clean{1} = dsearchn(data.timev_decimated',[1010 1190]');   % MT1
    % B_clean{2} = dsearchn(data.timev_decimated',[1010 1190]');   % MT3
    % B_clean{3} = dsearchn(data.timev_decimated',[1010 1190]');   % MT4
    % B_clean{4} = dsearchn(data.timev_decimated',[1010 1190]');   % MT5
    % B_clean{5} = dsearchn(data.timev_decimated',[1010 1190]');   % MT6
    % B_clean{6} = dsearchn(data.timev_decimated',[1010 1190]');   % MT7


    %% Extinction Session 1

    % Baseline
    % This follows the order of animals according to the main loop over ms values (see Analysis.m). Each session needs to be double-checked.

    % - IMPORTANT -
    % The baseline is being analyzed without noise or freezing periods.
    % All CS and ITI are being saved entirely. Only the plots are being saved without noise trials.


    % 2) Baseline timeepochs -> first minute. All animals checked. No noise
    % B_clean{1} = dsearchn(data.timev_decimated',[10 170]');    % MT1
    % B_clean{2} = dsearchn(data.timev_decimated',[10 170]');    % MT3
    % B_clean{3} = dsearchn(data.timev_decimated',[10 170]');    % MT4
    % B_clean{4} = dsearchn(data.timev_decimated',[10 170]');    % MT5
    % B_clean{5} = dsearchn(data.timev_decimated',[10 170]');    % MT6
    % B_clean{6} = dsearchn(data.timev_decimated',[10 170]');    % MT7


    %% Retrieval Session

    % Baseline
    % This follows the order of animals according to the main loop over ms values (see Analysis.m). Each session needs to be double-checked.

    % - IMPORTANT -
    % The baseline is being analyzed without noise or freezing periods.
    % All CS and ITI are being saved entirely. Only the plots are being saved without noise trials.

    % 1) Baseline timeepochs -> first minute. All animals checked. No noise
    %     B_clean{1} = dsearchn(data.timev_decimated',[10 170]');   % MT1
    %     B_clean{2} = dsearchn(data.timev_decimated',[10 170]');   % MT3
    %     B_clean{3} = dsearchn(data.timev_decimated',[10 170]');   % MT4
    %     B_clean{4} = dsearchn(data.timev_decimated',[10 170]');   % MT5
    %     B_clean{5} = dsearchn(data.timev_decimated',[10 170]');   % MT6
    %     B_clean{6} = dsearchn(data.timev_decimated',[10 170]');   % MT7


    %% Scripts
    p_welch_full_trials
    sFFT_spectrogram

    Coherence_full_trials
    PLV_phase_Full_Trials
    PLV_phase_permutat_

    Granger_mvgc_full_trials
    Granger_mvgc_full_trials_permutation
    Granger_mvgc_overtime

    %% Clear
    if ms < length(files.FilesLoaded{1, 1})
        %          clear('FilesLoaded','Path','data','parameters','newStr1','path', 'name' )
        clear('FilesLoaded','Path','newStr1','path', 'name' )

    else
        clear('id','FilesLoaded','Path','ms')
    end

end
toc

fprintf('\n Done. \n');

%% last update 12/10/2024
%  listening:
