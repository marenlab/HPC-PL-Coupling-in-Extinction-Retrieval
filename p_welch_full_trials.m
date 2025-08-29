%%  Welch power spectral density estimate 

% by Flavio Mourao.
% email: mourao.fg@tamu.edu
% Maren Lab - Department of Psychological and Brain Sciences
% Texas A&M University

% Started in:  12/2023
% Last update: 06/2024


% First extract the data with: Extracting_LFPs_and_events.m
% ... then organize the data with the script: Pre_processing.m

%%
fprintf('\n Welch power spectral density estimate ... \n');

%% pwelch
%  [pxx,f] = pwelch(x,window,noverlap,f,fs)

% - cell          
% - 1st column    - > baseline
% - 2nd column    - > cs events

% in each cell
% - rows          - > Hz
% - columns       - > channels
% 3th dimentions  - > CS-Trials

pw = [];

% Time window
pw.full_trial.baseline_timewin    = 2048; % in ms
pw.full_trial.CS_Trials_timewin   = 2048; % in ms
pw.full_trial.ITI_timewin         = 2048; % in ms

% Convert time window to points
pw.full_trial.baseline_timewinpnts   = hamming(round(pw.full_trial.baseline_timewin/(1000/parameters.decimated_srate)));
pw.full_trial.CS_Trials_timewinpnts  = hamming(round(pw.full_trial.CS_Trials_timewin/(1000/parameters.decimated_srate)));
pw.full_trial.ITI_timewinpnts        = hamming(round(pw.full_trial.ITI_timewin/(1000/parameters.decimated_srate)));

% nFFT
pw.full_trial.nFFT = 2^15; %4096; %2^nextpow2(pw.full_trial.baseline_timewinpnts));

% Number of overlap samples
pw.full_trial.overlap = 90;
pw.full_trial.baseline_noverlap = floor(pw.full_trial.overlap*0.01 * pw.full_trial.baseline_timewin);
pw.full_trial.CS_Trials_noverlap = floor(pw.full_trial.overlap*0.01 * pw.full_trial.CS_Trials_timewin);
pw.full_trial.ITI_noverlap = floor(pw.full_trial.overlap*0.01 * pw.full_trial.ITI_timewin);


% Baseline
not1 = 6;

for ii = 1:size(data.lfp{not1, 1},1)

    if ii == 1
        [pw.full_trial.Pxx{1,1}(ii,:),pw.full_trial.freq_baseline] = pwelch(data.lfp{not1, 1}(ii,B_clean{ms}(1):B_clean{ms}(2)),pw.full_trial.baseline_timewinpnts,pw.full_trial.baseline_noverlap,pw.full_trial.nFFT,parameters.decimated_srate);

    else
        pw.full_trial.Pxx{1,1}(ii,:) = pwelch(data.lfp{not1, 1}(ii,B_clean{ms}(1):B_clean{ms}(2)),pw.full_trial.baseline_timewinpnts,pw.full_trial.baseline_noverlap,pw.full_trial.nFFT,parameters.decimated_srate);

    end

end

clear ('ii')


CS-Trials or freezing epochs
Choose data from data.lfp according pre_processing.m define
not2 = 7;

for jj = 1:size(data.lfp{not2, 1},3)
    for ii = 1:size(data.lfp{not2, 1},1)
        if ii == 1
            [pw.full_trial.Pxx{2,1}(ii,:,jj),pw.full_trial.freq_CS_trials] = pwelch(data.lfp{not2, 1}(ii,:,jj),pw.full_trial.CS_Trials_timewinpnts,pw.full_trial.CS_Trials_noverlap,pw.full_trial.nFFT,parameters.decimated_srate);

        else
            pw.full_trial.Pxx{2,1}(ii,:,jj) = pwelch(data.lfp{not2, 1}(ii,:,jj),pw.full_trial.CS_Trials_timewinpnts,pw.full_trial.CS_Trials_noverlap,pw.full_trial.nFFT,parameters.decimated_srate);

        end

    end
end

clear ('jj','ii')

% % ITI or non-freezing epochs
% Choose data from data.lfp according pre_processing.m define
not3 = 8;

for jj = 1:size(data.lfp{not3, 1},3)
    for ii = 1:size(data.lfp{not3, 1},1)
        if ii == 1
            [pw.full_trial.Pxx{3,1}(ii,:,jj),pw.full_trial.freq_ITI] = pwelch(data.lfp{not3, 1}(ii,:,jj),pw.full_trial.ITI_timewinpnts,pw.full_trial.ITI_noverlap,pw.full_trial.nFFT,parameters.decimated_srate);

        else
            pw.full_trial.Pxx{3,1}(ii,:,jj) = pwelch(data.lfp{not3, 1}(ii,:,jj),pw.full_trial.ITI_timewinpnts,pw.full_trial.ITI_noverlap,pw.full_trial.nFFT,parameters.decimated_srate);

        end

    end
end

clear ('jj','ii')

% figure
% plot(pw.full_trial.freq_baseline,pw.full_trial.Pxx{1,1}(1,:))
% hold on
% plot(pw.full_trial.freq_CS_trials,mean(pw.full_trial.Pxx{2,1}(1,:,1:3),3))
% plot(pw.full_trial.freq_CS_trials,mean(pw.full_trial.Pxx{3,1}(1,:,1:3),3))
% 
% xlim([3 12])

clear ('not1','not2','not3')

%% Set range

steps                                    = diff(pw.full_trial.freq_baseline); % according to the fft time window

pw.full_trial.parameters.frex_2_12Hz     = 2:steps(1):12;
pw.full_trial.parameters.frex_idx_2_12Hz = dsearchn(pw.full_trial.freq_baseline,pw.full_trial.parameters.frex_2_12Hz');

%% Normalization and trials average

pw.full_trial.Pxx_Baseline_norm         = [];
pw.full_trial.Pxx_TotalPower_norm       = [];
pw.full_trial.Pxx_Baseline_norm_mean_first_trials    = [];
pw.full_trial.Pxx_TotalPower_norm_mean_first_trials  = [];
pw.full_trial.Pxx_Baseline_norm_mean_last_trials     = [];
pw.full_trial.Pxx_TotalPower_norm_mean_last_trials   = [];

% Plexon rescale to uV
% The original files *.nex5 did not have the correct scale.
scale = 1/0.153;

for jj = 1:size(pw.full_trial.Pxx,1)

    %Normalization
    %pw.full_trial.Pxx_Baseline_norm{jj,1}   = scale.*(pw.full_trial.Pxx{jj,1}(:,pw.full_trial.parameters.frex_idx_2_12Hz,:)./pw.full_trial.Pxx{1, 1}(:,pw.full_trial.parameters.frex_idx_2_12Hz)); % Baseline  normalization
    pw.full_trial.Pxx_TotalPower_norm{jj,1} = scale.*(pw.full_trial.Pxx{jj,1}(:,pw.full_trial.parameters.frex_idx_2_12Hz,:)./sum(pw.full_trial.Pxx{jj,1}(:,pw.full_trial.parameters.frex_idx_2_12Hz,:),2)); % total power normalization


    % Averaging trials

    if jj == 1 % baseline condition
        pw.full_trial.Pxx_Baseline_norm_mean{jj,1}   = pw.full_trial.Pxx_Baseline_norm{jj,1};   % baseline cell will in 1,1 of course....but just to keep the organization
        pw.full_trial.Pxx_TotalPower_norm_mean{jj,1} = pw.full_trial.Pxx_TotalPower_norm{jj,1}; % baseline cell will in 1,1 of course....but just to keep the organization

    else

        pw.full_trial.Pxx_Baseline_norm_mean{jj,1}    = mean(pw.full_trial.Pxx_Baseline_norm{jj,1}(:,:,CSIT{ms}),3);   % averaging first 10 trials
        pw.full_trial.Pxx_TotalPower_norm_mean{jj,1}  = mean(pw.full_trial.Pxx_TotalPower_norm{jj,1}(:,:,CSIT{ms}),3); % averaging first 10 trials
 
    end

end

clear('jj')

%% Plot to check - for now I am considring just CS-Tone

%data_2_plot = pw.full_trial.Pxx_stats_baseline;
%data_2_plot      = pw.full_trial.Pxx_TotalPower_norm;

data_2_plot{1,1}     = pw.full_trial.Pxx{1,1}(:,pw.full_trial.parameters.frex_idx_2_12Hz);

data_2_plot_mean = pw.full_trial.Pxx_TotalPower_norm_mean;

% choose channel
ch1 = 1;
ch2 = 2;
ch3 = 3;


f = figure;
set(gcf, 'Position', get(0, 'Screensize'));
set(gcf,'color','w');
axis('square')
sgtitle('\fontsize{18} \bf Extinction')

% Baseline
b(1) = subplot(3,7,1);
plot(pw.full_trial.freq_baseline(pw.full_trial.parameters.frex_idx_2_12Hz),data_2_plot{1, 1}(ch1,:),'linewidth',3,'color',[.6 .6 .6] )
xlim([2 12])
ylabel({'\fontsize{16} \bf mPFC (PL)';'Normalized Amplitude (a.u.)';[]})
xlabel('Hz')
%yticks([0 0.05 0.1 0.15 0.2])
box off
title({'\fontsize{16} \bf Baseline';[]})

b(2) = subplot(3,7,8);
plot(pw.full_trial.freq_baseline(pw.full_trial.parameters.frex_idx_2_12Hz),data_2_plot{1, 1}(ch2,:),'linewidth',3,'color',[.6 .6 .6] )
xlim([2 12])
ylabel({'\fontsize{16} \bf mPFC (IL)';'Normalized Amplitude (a.u.)';[]})
box off
xlabel('Hz')
%yticks([0 0.05 0.1 0.15 0.2])

b(3) = subplot(3,7,15);
plot(pw.full_trial.freq_baseline(pw.full_trial.parameters.frex_idx_2_12Hz),data_2_plot{1, 1}(ch3,:),'linewidth',3,'color',[.6 .6 .6] )
xlim([2 12])
box off
ylabel({'\fontsize{16} \bf dHPC';'Normalized Amplitude (a.u.)';[]})
xlabel('Hz')
%yticks([0 0.05 0.1 0.15 0.2])


% CS-Tones

for ii = 1:size(CSIT{ms},2)

    % First 10 trials
    t(ii+1) = subplot(3,7,ii+1);
    plot(pw.full_trial.freq_baseline(pw.full_trial.parameters.frex_idx_2_12Hz),data_2_plot{2, 1}(ch1,:,CSIT{ms}(1,ii)),'linewidth',3,'Color',[0 0 .6])

    xlim([2 12])
    %ylabel('Amplitude')
    xlabel('Hz')
    axis('square')
    box off
    title({['CS-Tone ' num2str(CSIT{ms}(1,ii))];[]})

    t(ii+13) = subplot(3,7,ii+8);
    plot(pw.full_trial.freq_baseline(pw.full_trial.parameters.frex_idx_2_12Hz),data_2_plot{2, 1}(ch2,:,CSIT{ms}(1,ii)),'linewidth',3,'Color',[0 0 .6])

    xlim([2 12])
    %ylabel('Amplitude')
    axis('square')
    box off
    xlabel('Hz')
    title({['CS-Tone ' num2str(CSIT{ms}(1,ii))];[]})

    t(ii+25) = subplot(3,7,ii+15);
    plot(pw.full_trial.freq_baseline(pw.full_trial.parameters.frex_idx_2_12Hz),data_2_plot{2, 1}(ch3,:,CSIT{ms}(1,ii)),'linewidth',3,'Color',[0 0 .6])

    xlim([2 12])
    %ylabel('Amplitude')
    axis('square')
    box off
    xlabel('Hz')
    title({['CS-Tone ' num2str(CSIT{ms}(1,ii))];[]})

end


% Trial averaged

m(1) = subplot(3,7,7);
plot(pw.full_trial.freq_baseline(pw.full_trial.parameters.frex_idx_2_12Hz),data_2_plot_mean{2, 1}(ch1,:),'linewidth',3,'color',[0 0 .6] )

xlim([2 12])
xlabel('Hz')
yticks([0 0.05 0.1 0.15 0.2])
box off
title({'\fontsize{16} \bf Averaged Trials';[]})


m(2) = subplot(3,7,14);
plot(pw.full_trial.freq_baseline(pw.full_trial.parameters.frex_idx_2_12Hz),data_2_plot_mean{2, 1}(ch2,:),'linewidth',3,'color',[0 0 .6] )

xlim([2 12])
box off
xlabel('Hz')
yticks([0 0.05 0.1 0.15 0.2])
%title({'Averaged Trials'})


m(3) = subplot(3,7,21);
plot(pw.full_trial.freq_baseline(pw.full_trial.parameters.frex_idx_2_12Hz),data_2_plot_mean{2, 1}(ch3,:),'linewidth',3,'color',[0 0 .6] )

xlim([2 12])
box off
xlabel('Hz')
yticks([0 0.05 0.1 0.15 0.2])
%title({'Averaged Trials'})
legend('\fontsize{14} First 10 Trials','\fontsize{14} Last 10 Trials')
lg = legend;
set(lg,...
    'Position',[0.8859375 0.290602184305253 0.0497395833333333 0.0241788321167883]);

linkaxes([b t m],'xy');
%b(1).YLim  = [0 15];
b(1).YLim  = [0 0.20];


clear ('ii','ch1','ch2','ch3','f','lg','m','t','b','data_2_plot','data_2_plot_mean1','data_2_plot_mean2')

%% Save figure

newStr = id(1:end-4);

path = '/Users/flavio/Desktop';
%path = files.FilesLoaded{1, 1}.folder;

%name = strcat(path,'/',newStr,'_pw_Total_Power_best_5_min_blocks_baseline_norm');
name = strcat(path,'/',newStr,'_pw_full_trials_Total_Power');

saveas(gcf,name,'png')

set(gcf,'renderer','Painters')
exportgraphics(gcf,strcat(name,'.eps'),'Resolution', 300)

close all

clear('name','newStr1','path')


%% Extract data from Total power normalization

% Set freuency range
pw.full_trial.parameters.frex_3_6Hz     = 3:steps(1):6;
pw.full_trial.parameters.frex_idx_3_6Hz = dsearchn(pw.full_trial.parameters.frex_2_12Hz',pw.full_trial.parameters.frex_3_6Hz');

pw.full_trial.parameters.frex_6_8Hz     = 6:steps(1):8;
pw.full_trial.parameters.frex_idx_6_8Hz = dsearchn(pw.full_trial.parameters.frex_2_12Hz',pw.full_trial.parameters.frex_6_8Hz');

pw.full_trial.parameters.frex_7_10Hz     = 7:steps(1):10;
pw.full_trial.parameters.frex_idx_7_10Hz = dsearchn(pw.full_trial.parameters.frex_2_12Hz',pw.full_trial.parameters.frex_7_10Hz');


% Integrated band frequencies
% Integrated band frequency - All Trials
pw.full_trial.Pxx_TotalPower_norm{1,2} = squeeze(mean(pw.full_trial.Pxx_TotalPower_norm{1,1}(:,pw.full_trial.parameters.frex_idx_3_6Hz,:),2)); % Baseline
pw.full_trial.Pxx_TotalPower_norm{2,2} = squeeze(mean(pw.full_trial.Pxx_TotalPower_norm{2,1}(:,pw.full_trial.parameters.frex_idx_3_6Hz,:),2)); % CS-Trials
pw.full_trial.Pxx_TotalPower_norm{3,2} = squeeze(mean(pw.full_trial.Pxx_TotalPower_norm{3,1}(:,pw.full_trial.parameters.frex_idx_3_6Hz,:),2)); % ITI

pw.full_trial.Pxx_TotalPower_norm{1,3} = squeeze(mean(pw.full_trial.Pxx_TotalPower_norm{1,1}(:,pw.full_trial.parameters.frex_idx_6_8Hz,:),2)); % Baseline
pw.full_trial.Pxx_TotalPower_norm{2,3} = squeeze(mean(pw.full_trial.Pxx_TotalPower_norm{2,1}(:,pw.full_trial.parameters.frex_idx_6_8Hz,:),2)); % CS-Trials
pw.full_trial.Pxx_TotalPower_norm{3,3} = squeeze(mean(pw.full_trial.Pxx_TotalPower_norm{3,1}(:,pw.full_trial.parameters.frex_idx_6_8Hz,:),2)); % ITI

pw.full_trial.Pxx_TotalPower_norm{1,4} = squeeze(mean(pw.full_trial.Pxx_TotalPower_norm{1,1}(:,pw.full_trial.parameters.frex_idx_7_10Hz,:),2)); % Baseline
pw.full_trial.Pxx_TotalPower_norm{2,4} = squeeze(mean(pw.full_trial.Pxx_TotalPower_norm{2,1}(:,pw.full_trial.parameters.frex_idx_7_10Hz,:),2)); % CS-Trials
pw.full_trial.Pxx_TotalPower_norm{3,4} = squeeze(mean(pw.full_trial.Pxx_TotalPower_norm{3,1}(:,pw.full_trial.parameters.frex_idx_7_10Hz,:),2)); % ITI


% Integrated band frequency averaging
pw.full_trial.Pxx_TotalPower_norm_mean{1,2} = squeeze(mean(pw.full_trial.Pxx_TotalPower_norm_mean{1,1}(:,pw.full_trial.parameters.frex_idx_3_6Hz,:),2)); % Baseline
pw.full_trial.Pxx_TotalPower_norm_mean{2,2} = squeeze(mean(pw.full_trial.Pxx_TotalPower_norm_mean{2,1}(:,pw.full_trial.parameters.frex_idx_3_6Hz,:),2)); % CS-Trials
pw.full_trial.Pxx_TotalPower_norm_mean{3,2} = squeeze(mean(pw.full_trial.Pxx_TotalPower_norm_mean{3,1}(:,pw.full_trial.parameters.frex_idx_3_6Hz,:),2)); % ITI

pw.full_trial.Pxx_TotalPower_norm_mean{1,3} = squeeze(mean(pw.full_trial.Pxx_TotalPower_norm_mean{1,1}(:,pw.full_trial.parameters.frex_idx_6_8Hz,:),2)); % Baseline
pw.full_trial.Pxx_TotalPower_norm_mean{2,3} = squeeze(mean(pw.full_trial.Pxx_TotalPower_norm_mean{2,1}(:,pw.full_trial.parameters.frex_idx_6_8Hz,:),2)); % CS-Trials
pw.full_trial.Pxx_TotalPower_norm_mean{3,3} = squeeze(mean(pw.full_trial.Pxx_TotalPower_norm_mean{3,1}(:,pw.full_trial.parameters.frex_idx_6_8Hz,:),2)); % ITI

pw.full_trial.Pxx_TotalPower_norm_mean{1,4} = squeeze(mean(pw.full_trial.Pxx_TotalPower_norm_mean{1,1}(:,pw.full_trial.parameters.frex_idx_7_10Hz,:),2)); % Baseline
pw.full_trial.Pxx_TotalPower_norm_mean{2,4} = squeeze(mean(pw.full_trial.Pxx_TotalPower_norm_mean{2,1}(:,pw.full_trial.parameters.frex_idx_7_10Hz,:),2)); % CS-Trials
pw.full_trial.Pxx_TotalPower_norm_mean{3,4} = squeeze(mean(pw.full_trial.Pxx_TotalPower_norm_mean{3,1}(:,pw.full_trial.parameters.frex_idx_7_10Hz,:),2)); % ITI



% Peak Frequency
% Peak frequency range - All Trials
pw.full_trial.Pxx_TotalPower_norm_peak{1,2} = squeeze(max(pw.full_trial.Pxx_TotalPower_norm{1,1}(:,pw.full_trial.parameters.frex_idx_3_6Hz,:),[],2)); % Baseline
pw.full_trial.Pxx_TotalPower_norm_peak{2,2} = squeeze(max(pw.full_trial.Pxx_TotalPower_norm{2,1}(:,pw.full_trial.parameters.frex_idx_3_6Hz,:),[],2)); % CS-Trials
pw.full_trial.Pxx_TotalPower_norm_peak{3,2} = squeeze(max(pw.full_trial.Pxx_TotalPower_norm{3,1}(:,pw.full_trial.parameters.frex_idx_3_6Hz,:),[],2)); % ITI

pw.full_trial.Pxx_TotalPower_norm_peak{1,3} = squeeze(max(pw.full_trial.Pxx_TotalPower_norm{1,1}(:,pw.full_trial.parameters.frex_idx_6_8Hz,:),[],2)); % Baseline
pw.full_trial.Pxx_TotalPower_norm_peak{2,3} = squeeze(max(pw.full_trial.Pxx_TotalPower_norm{2,1}(:,pw.full_trial.parameters.frex_idx_6_8Hz,:),[],2)); % CS-Trials
pw.full_trial.Pxx_TotalPower_norm_peak{3,3} = squeeze(max(pw.full_trial.Pxx_TotalPower_norm{3,1}(:,pw.full_trial.parameters.frex_idx_6_8Hz,:),[],2)); % ITI

pw.full_trial.Pxx_TotalPower_norm_peak{1,4} = squeeze(max(pw.full_trial.Pxx_TotalPower_norm{1,1}(:,pw.full_trial.parameters.frex_idx_7_10Hz,:),[],2)); % Baseline
pw.full_trial.Pxx_TotalPower_norm_peak{2,4} = squeeze(max(pw.full_trial.Pxx_TotalPower_norm{2,1}(:,pw.full_trial.parameters.frex_idx_7_10Hz,:),[],2)); % CS-Trials
pw.full_trial.Pxx_TotalPower_norm_peak{3,4} = squeeze(max(pw.full_trial.Pxx_TotalPower_norm{3,1}(:,pw.full_trial.parameters.frex_idx_7_10Hz,:),[],2)); % ITI

clear('steps')

%% Save data

%newStr = regexprep(files.id.name,'.mat','_');
%newStr = files.id(ms).name(1:end-8);
newStr = id(1:end-4);

path = '/Users/flavio/Desktop';
%path = files.FilesLoaded{1, 1}.folder;

name = strcat(path,'/',newStr,'_pw_full_trials_Total_Power.mat');

% save data
save(name,'pw','-v7.3')

clear('name','newStr','path')

%% last update 06/06/2024
%  listening:
