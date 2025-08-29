
%% Phase Coherence analyses based on Hilbert Transform (Instantaneous phase)
%  - Performs analysis considering the trial periods

% The code relies on the following package:
% --> circular-statistics-toolbox
%     https://github.com/circstat/circstat-matlab


% by Flavio Mourao.
% email: mourao.fg@tamu.edu
% Maren Lab - Department of Psychological and Brain Sciences
% Texas A&M University

% Started in:  04/2024
% Last update: 05/2024

%%
% First extract the data with: Extracting_LFPs_and_events.m
% ... then organize the data with the script: Pre_processing.m

%% Run each session sequentially

%%
fprintf('\n Phase analyses based on Hilbert Transform... \n');

%% Hilbert Transform

% Loop over channels and make hilbert transform
% Here, periods without noise are considered for both the baseline and the trials.
% Defined by the variables: Baseline: B_clean / CS-Tone + 8Hz opto stimulation CSIT /
% ITI:CSIT_1 --> not defined yet


hilb = [];

hilb.parameters.filters = [3 6;6 8;7 10];

% Baseline
for ii = 1:size(data.lfp{6,1},1)
    for jj = 1:size(hilb.parameters.filters,1)

        temp = eegfilt(data.lfp{6,1}(ii,B_clean{1,ms}(1,1):B_clean{1,ms}(2,1)),parameters.decimated_srate,hilb.parameters.filters(jj,1),hilb.parameters.filters(jj,2)); % just filtering
        hilb.amplitude{1,jj}(ii,:) = abs(hilbert(temp)); % getting the amplitude envelope
        hilb.phase{1,jj}(ii,:) = angle(hilbert(temp)); % getting the angles

    end

end


%  CS-Tone
for ii = 1:size(data.lfp{7,1},1)
    for tt = 1:size(CSIT{ms},2)
        for jj = 1:size(hilb.parameters.filters,1)

            temp = eegfilt(data.lfp{7,1}(ii,:,CSIT{ms}(1,tt)),parameters.decimated_srate,hilb.parameters.filters(jj,1),hilb.parameters.filters(jj,2)); % just filtering
            hilb.amplitude{2,jj}(ii,:,tt) = abs(hilbert(temp)); % getting the amplitude envelope
            hilb.phase{2,jj}(ii,:,tt) = angle(hilbert(temp)); % getting the angles

        end
    end

end



% ITI
for ii = 1:size(data.lfp{8,1},1)
    for tt = 1:size(CSIT{ms},2)
        for jj = 1:size(hilb.parameters.filters,1)

            temp = eegfilt(data.lfp{8,1}(ii,:,CSIT{ms}(1,tt)),parameters.decimated_srate,hilb.parameters.filters(jj,1),hilb.parameters.filters(jj,2)); % just filtering
            hilb.amplitude{3,jj}(ii,:,tt) = abs(hilbert(temp)); % getting the amplitude envelope
            hilb.phase{3,jj}(ii,:,tt) = angle(hilbert(temp)); % getting the angles

        end
    end

end


clear('ii','jj','tt','temp')

%% Linear timing normalization. To compare phase clusters baseline with CS-tone and ITI period

data_iterp_amplitude = cell(size(hilb.amplitude));

for ii = 1:size(hilb.amplitude,2)
    for cc = 1:size(hilb.amplitude{2,ii},1)

        for tt = 1:size(hilb.amplitude{2,ii},3)
            % CS-Tone + 8Hz opto stimulation from beginning to end ( 5s opto --> 10s opto + CS-Tone --> 5s opto)
            xq1 = linspace(1,length(hilb.amplitude{2,ii}),length(hilb.amplitude{1,ii}));
            data_iterp_amplitude{2,ii}(cc,:,tt) = interp1(1:length(hilb.amplitude{2,ii}(cc,:,tt)),hilb.amplitude{2,ii}(cc,:,tt),xq1,'spline');
        end

        for tt = 1:size(hilb.amplitude{3,ii},3)
            % CS-Tone + 8Hz opto stimulation (10s opto + CS-Tone)
            xq2 = linspace(1,length(hilb.amplitude{3,ii}),length(hilb.amplitude{1,ii}));
            data_iterp_amplitude{3,ii}(cc,:,tt) = interp1(1:length(hilb.amplitude{3,ii}(cc,:,tt)),hilb.amplitude{3,ii}(cc,:,tt),xq2,'spline');
        end

    end
end

% Substitute previous values for interpoleted ones
hilb.amplitude(2:3,:) = data_iterp_amplitude(2:3,:);


data_iterp_phase = cell(size(hilb.phase));

for ii = 1:size(hilb.phase,2)
    for cc = 1:size(hilb.phase{2,ii},1)

        for tt = 1:size(hilb.phase{2,ii},3)
            % CS-Tone + 8Hz opto stimulation from beginning to end ( 5s opto --> 10s opto + CS-Tone --> 5s opto)
            xq1 = linspace(1,length(hilb.phase{2,ii}),length(hilb.phase{1,ii}));
            data_iterp_phase{2,ii}(cc,:,tt) = interp1(1:length(hilb.phase{2,ii}(cc,:,tt)),hilb.phase{2,ii}(cc,:,tt),xq1);
        end

        for tt = 1:size(hilb.phase{3,ii},3)
            % CS-Tone + 8Hz opto stimulation (10s opto + CS-Tone)
            xq2 = linspace(1,length(hilb.phase{3,ii}),length(hilb.phase{1,ii}));
            data_iterp_phase{3,ii}(cc,:,tt) = interp1(1:length(hilb.phase{3,ii}(cc,:,tt)),hilb.phase{3,ii}(cc,:,tt),xq2);

        end
    end
end

% Substitute previous values for interpoleted ones
hilb.phase(2:3,:) = data_iterp_phase(2:3,:);


clear('ii','cc','tt','xq1','xq2','data_iterp_amplitude','data_iterp_phase')

%%
hilb.parameters.time_pre = 5;
hilb.parameters.time_pre_idx = hilb.parameters.time_pre * parameters.decimated_srate;

% time_pos = 5;
% time_pos_idx = time_pos * parameters.decimated_srate;

% CS-Tone + 8Hz opto stimulation from beginning to end ( 5s opto --> 10s opto + CS-Tone --> 5s opto)
for ii = 1:size(data.lfp{7,1},1)
    for tt = 1:size(CSIT{ms},2)
        for jj = 1:size(hilb.parameters.filters,1)

            temp = eegfilt(data.lfp{5,1}(ii, data.events{2, 1}(CSIT{ms}(1,tt),1) - hilb.parameters.time_pre_idx : data.events{2, 1}(CSIT{ms}(1,tt),2) + hilb.parameters.time_pre_idx),...
                parameters.decimated_srate,hilb.parameters.filters(jj,1),hilb.parameters.filters(jj,2)); % just filtering

            hilb.amplitude{4,jj}(ii,:,tt) = abs(hilbert(temp)); % getting the amplitude envelope
            hilb.phase{4,jj}(ii,:,tt) = angle(hilbert(temp)); % getting the angles
        end
    end

end


parameters.timev = linspace(-hilb.parameters.time_pre, size(data.lfp{6,1},2)./parameters.decimated_srate + hilb.parameters.time_pre, size(data.lfp{6,1},2)+ hilb.parameters.time_pre_idx);

clear('ii','jj','tt','temp')

%% - CHANNELS MAP -
% *just in case to check

% mPFC
% .Row 1 -> pre limbic
% .Row 2 -> infra limbic

% Hippocampus
% .Row 3 -> CA1

%% Delta phase from Euler representation of angles. All samples. - Full trials


% Difference between channels

% all possible combinations between all channels:
hilb.parameters.combinations  = nchoosek(1:size(hilb.phase{1,1},1),2);

% Initilize variable

% Delta Phase
hilb.delta_phase = cell(size(hilb.phase));

% Euler differences
hilb.euler_differences = cell(size(hilb.phase));

% All time series
hilb.PLV_trial = cell(size(hilb.phase));
hilb.delta_phase_trial = cell(size(hilb.phase));

% Average to get trials
hilb.PLV_mean_trial = cell(size(hilb.phase));
hilb.delta_phase_mean_trial = cell(size(hilb.phase));

% Average to get session
hilb.PLV_mean_session = cell(size(hilb.phase));
hilb.delta_phase_mean_session = cell(size(hilb.phase));


% Extracts relative phase and circular variance length (PLV)
for hh = 1:size(hilb.phase,1)
    for ii = 1:size(hilb.phase,2)

        for jj = 1:size(hilb.parameters.combinations,1)


            if hh == 1 % condition for baseline
                hilb.delta_phase{hh,ii}(jj,:)       = squeeze(hilb.phase{hh,ii}(hilb.parameters.combinations(jj,1),:) - hilb.phase{hh,ii}(hilb.parameters.combinations(jj,2),:))';
                hilb.euler_differences{hh,ii}(jj,:) = exp(1i*(hilb.delta_phase{hh,ii}(jj,:)));

                hilb.PLV_trial{hh,ii}(jj,:)                 = abs(hilb.euler_differences{hh,ii}(jj,:)); % Be CAREFULL. Not work considering isolated samples. Over time is necessary a sliding windown. See section below
                hilb.delta_phase_trial{hh,ii}(jj,:)         = angle(hilb.euler_differences{hh,ii}(jj,:));

                hilb.PLV_mean_trial{hh,ii}(jj,:)            = abs(mean(hilb.euler_differences{hh,ii}(jj,:),2));
                hilb.delta_phase_mean_trial{hh,ii}(jj,:)    = angle(mean(hilb.euler_differences{hh,ii}(jj,:),2));

                hilb.PLV_mean_session{hh,ii}(jj,:)            = abs(mean(mean(hilb.euler_differences{hh,ii}(jj,:),2),2));
                hilb.delta_phase_mean_session{hh,ii}(jj,:)    = angle(mean(mean(hilb.euler_differences{hh,ii}(jj,:),2),2));

            else % condition for CS and ITI
                hilb.delta_phase{hh,ii}(jj,:,:)        = squeeze(hilb.phase{hh,ii}(hilb.parameters.combinations(jj,1),:,:) - hilb.phase{hh,ii}(hilb.parameters.combinations(jj,2),:,:));
                hilb.euler_differences{hh,ii}(jj,:,:)  = exp(1i*(hilb.delta_phase{hh,ii}(jj,:,:)));

                hilb.PLV_trial{hh,ii}(jj,:,:)                 = abs(hilb.euler_differences{hh,ii}(jj,:,:)); % Be CAREFULL. Not work considering isolated samples. Over time is necessary a sliding windown See section below
                hilb.delta_phase_trial{hh,ii}(jj,:,:)         = angle(hilb.euler_differences{hh,ii}(jj,:,:));

                hilb.PLV_mean_trial{hh,ii}(jj,:,:)           = abs(mean(hilb.euler_differences{hh,ii}(jj,:,:),2));
                hilb.delta_phase_mean_trial{hh,ii}(jj,:,:)   = angle(mean(hilb.euler_differences{hh,ii}(jj,:,:),2));

                hilb.PLV_mean_session{hh,ii}(jj,:)            = abs(mean(mean(hilb.euler_differences{hh,ii}(jj,:,:),2),3));
                hilb.delta_phase_mean_session{hh,ii}(jj,:)    = angle(mean(mean(hilb.euler_differences{hh,ii}(jj,:,:),2),3));

            end
        end
    end
end

clear('hh','jj','ii')

%% Delta phase from Euler representation of angles. Sliding window. Time bins and overlap.
% Extracts relative phase and length of circular variance (PLV)

% Define time bins
time_bins     = 3; % seconds
time_bins_idx = time_bins * parameters.decimated_srate;

% Overlap
timeoverlap    = .95; % percentage
overlap = round((time_bins_idx)-(timeoverlap*time_bins_idx));

% Euler differences
hilb.euler_differences_bins = cell(size(hilb.phase));

% All time series
hilb.PLV_trial_bins         = cell(size(hilb.phase));
hilb.delta_phase_trial_bins = cell(size(hilb.phase));

% Average to get trials
hilb.PLV_mean_trial_bins         = cell(size(hilb.phase));
hilb.delta_phase_mean_trial_bins = cell(size(hilb.phase));

% Average to get session
hilb.PLV_mean_session_bins         = cell(size(hilb.phase));
hilb.delta_phase_mean_session_bins = cell(size(hilb.phase));


% Extracts relative phase and circular variance length (PLV)
for ii = 1:size(hilb.phase,2)
    for hh = 1:size(hilb.phase,1)
        for jj = 1:size(hilb.phase{hh,ii},1)

            bins = (2:overlap:size(hilb.phase{hh,ii},2)-time_bins_idx);

            % loop and overlap over time bins
            for bb = 2:size(bins,2)

                if hh == 1 % condition for baseline
                    hilb.euler_differences_bins{hh,ii}(jj,bb-1)        = mean(hilb.euler_differences{hh,ii}(jj,bins(bb-1):bins(bb-1) + time_bins_idx -1),2); % euler_differences already extracted from section above

                    hilb.PLV_trial_bins{hh,ii}(jj,bb-1)                = abs(hilb.euler_differences_bins{hh,ii}(jj,bb-1));
                    hilb.delta_phase_trial_bins{hh,ii}(jj,bb-1)        = angle(hilb.euler_differences_bins{hh,ii}(jj,bb-1));


                else % condition for CS and ITI
                    hilb.euler_differences_bins{hh,ii}(jj,bb-1,:) = mean(hilb.euler_differences{hh,ii}(jj,bins(bb-1):bins(bb-1) + time_bins_idx -1,:),2);

                    hilb.PLV_trial_bins{hh,ii}(jj,bb-1,:)              = abs(hilb.euler_differences_bins{hh,ii}(jj,bb-1,:));
                    hilb.delta_phase_trial_bins{hh,ii}(jj,bb-1,:)      = angle(hilb.euler_differences_bins{hh,ii}(jj,bb-1,:));


                end
            end

            if hh == 1 % condition for baseline

                hilb.PLV_mean_trial_bins{hh,ii}(jj,:)           = abs(mean(hilb.euler_differences_bins{hh,ii}(jj,:),2));
                hilb.delta_phase_mean_trial_bins{hh,ii}(jj,:)   = angle(mean(hilb.euler_differences_bins{hh,ii}(jj,:),2));

                hilb.PLV_mean_session_bins{hh,ii}(jj,:)             = abs(mean(mean(hilb.euler_differences_bins{hh,ii}(jj,:),2),3));
                hilb.delta_phase_mean_session_bins{hh,ii}(jj,:)     = angle(mean(mean(hilb.euler_differences_bins{hh,ii}(jj,:),2),3));

            else % condition for CS and ITI

                hilb.PLV_mean_trial_bins{hh,ii}(jj,:)           = abs(mean(hilb.euler_differences_bins{hh,ii}(jj,:,:),2));
                hilb.delta_phase_mean_trial_bins{hh,ii}(jj,:)   = angle(mean(hilb.euler_differences_bins{hh,ii}(jj,:),2));

                hilb.PLV_mean_session_bins{hh,ii}(jj,:)           = abs(mean(mean(hilb.euler_differences_bins{hh,ii}(jj,:,:),2),3));
                hilb.delta_phase_mean_session_bins{hh,ii}(jj,:)   = angle(mean(mean(hilb.euler_differences_bins{hh,ii}(jj,:,:),2),3));
            end

        end
    end
end

clear('time_bins','time_bins_idx','hh','jj','ii','bb','timeoverlap','overlap','bins')


%% 5-Trials average

hilb.delta_phase_5_trial = [];
hilb.delta_phase_mean_5_trial = [];
hilb.PLV_mean_5_trial = [];

for ii = 1:size(hilb.delta_phase,1)
    for jj = 1:size(hilb.delta_phase,2)
        
        if ii == 1
            hilb.delta_phase_5_trial{ii,jj}      = hilb.delta_phase{ii,jj};
            hilb.delta_phase_mean_5_trial{ii,jj} = hilb.delta_phase_mean_trial{ii,jj};
            hilb.PLV_mean_5_trial{ii,jj}         = hilb.PLV_mean_trial{ii,jj};
        else
              
                hilb.delta_phase_5_trial{ii,jj}        = circ_mean(hilb.delta_phase{ii,jj},[],3);
                hilb.delta_phase_mean_5_trial{ii,jj}   = circ_mean(hilb.delta_phase_mean_trial{ii,jj},[],2);
                hilb.PLV_mean_5_trial{ii,jj}           = mean(hilb.PLV_mean_trial{ii,jj},2);


        end
    end
end


clear('ii','jj','mm','c')

%% Polar PLots. Baseline and Mean Trials. 4 Hz and 8 Hz

% Define data to plot
data_phase = hilb.delta_phase_5_trial;

data_phase_mean = hilb.delta_phase_mean_5_trial  ;
data_PLV_mean = hilb.PLV_mean_5_trial ;

sgtitle_ = {'3-6 Hz Bandpass filter';'6-8 Hz Bandpass filter';'7-10 Hz Bandpass filter'};

% Baseline
for ii = 1:size(data_phase,2)

    figure(ii);
    set(gcf, 'Position', get(0, 'Screensize'));
    set(gcf,'color','w');

    sgtitle({sgtitle_{ii};''},'FontWeight','bold','FontSize',16);

    % mPFC PL <--> mPFC IL
    subplot(3,2,1)
    polarplot([zeros(1,size(data_phase{1,ii}(1,1:240:end),2)) data_phase{1,ii}(1,1:240:end)]',repmat([0 1],1,size(data_phase{1,ii}(1,1:240:end),2))','color',[.6 .6 .6]);
    hold on
    I=polarplot([0 data_phase_mean{1,ii}(1,1)],[0 data_PLV_mean{1,ii}(1,1)]);
    set(I,'linewidth',6,'color',[0 0 0])
    thetaticks([0 90 180 270])
    title({['\fontsize{16} Baseline'];[];['\fontsize{14} PLV = ', num2str(data_PLV_mean{1,ii}(1,1))];[]})

    % mPFC PL <--> dHPC
    subplot(3,2,3)
    polarplot([zeros(1,size(data_phase{1,ii}(2,1:240:end),2)) data_phase{1,ii}(2,1:240:end)]',repmat([0 1],1,size(data_phase{1,ii}(2,1:240:end),2))','color',[.6 .6 .6]);
    hold on
    I=polarplot([0 data_phase_mean{1,ii}(2,1)],[0 data_PLV_mean{1,ii}(2,1)]);
    set(I,'linewidth',6,'color',[0 0 0])
    thetaticks([0 90 180 270])
    title({['\fontsize{14} PLV = ', num2str(data_PLV_mean{1,ii}(2,1))];[]})

    % mPFC IL <--> dHPC
    subplot(3,2,5)
    polarplot([zeros(1,size(data_phase{1,ii}(3,1:240:end),2)) data_phase{1,ii}(3,1:240:end)]',repmat([0 1],1,size(data_phase{1,ii}(3,1:240:end),2))','color',[.6 .6 .6]);
    hold on
    I=polarplot([0 data_phase_mean{1,ii}(3,1)],[0 data_PLV_mean{1,ii}(3,1)]);
    set(I,'linewidth',6,'color',[0 0 0])
    thetaticks([0 90 180 270])
    title({['\fontsize{14} PLV = ', num2str(data_PLV_mean{1,ii}(3,1))];[]})


    annotation(figure(ii),'textbox',...
        [0.0172619047619056 0.807610993657511 0.101190476190476 0.0243699788583516],...
        'String','mPFC PL <--> mPFC IL',...
        'FontWeight','bold',...
        'FontSize',14,...
        'FontName','Helvetica Neue',...
        'FitBoxToText','off',...
        'EdgeColor',[1 1 1]);


    annotation(figure(ii),'textbox',...
        [0.022619047619049 0.524312896405924 0.101190476190477 0.0243699788583518],...
        'String','mPFC PL <--> dHPC',...
        'FontWeight','bold',...
        'FontSize',14,...
        'FontName','Helvetica Neue',...
        'FitBoxToText','off',...
        'EdgeColor',[1 1 1]);

    annotation(figure(ii),'textbox',...
        [0.0267857142857152 0.250528541226218 0.101190476190476 0.0243699788583517],...
        'String','mPFC IL <--> dHPC',...
        'FontWeight','bold',...
        'FontSize',14,...
        'FontName','Helvetica Neue',...
        'FitBoxToText','off',...
        'EdgeColor',[1 1 1]);



    annotation(figure(ii),'textbox',...
    [0.0148809523809532 0.957716701902754 0.101190476190476 0.0243699788583516],...
    'String',id(1:end-11),...
    'FontWeight','bold',...
    'FontSize',18,...
    'FontName','Helvetica Neue',...
    'FitBoxToText','off',...
    'EdgeColor',[1 1 1]);

end

% Mean
for ii = 1:size(data_phase,2)
    figure(ii)
    for tt = 1:size(data_phase{3,ii},3)

        if tt == 1
            color = [0 0 .6];
            titl_ = '\fontsize{16} Averaged Trials';
        else
            color = [0 0 .6];
            titl_ = '\fontsize{16} Averaged Trials';
        end

        % mPFC PL <--> mPFC IL
        subplot(3,2,tt+1)
        polarplot([zeros(1,size(data_phase{2,ii}(1,1:240:end,tt),2)) data_phase{2,ii}(1,1:240:end,tt)]',repmat([0 1],1,size(data_phase{2,ii}(1,1:240:end,tt),2))','color',[.6 .6 .6]);
        hold on
        I=polarplot([0 data_phase_mean{2,ii}(1,tt)],[0 data_PLV_mean{2,ii}(1,tt)]);
        set(I,'linewidth',6,'color',color)
        thetaticks([0 90 180 270])
        title({[titl_];[];['\fontsize{14} PLV = ' num2str(data_PLV_mean{2,ii}(1,tt))];[]},'color',color)

        % mPFC PL <--> dHPC
        subplot(3,2,tt+3)
        polarplot([zeros(1,size(data_phase{2,ii}(2,1:240:end,tt),2)) data_phase{2,ii}(2,1:240:end,tt)]',repmat([0 1],1,size(data_phase{2,ii}(2,1:240:end,tt),2))','color',[.6 .6 .6]);
        hold on
        I=polarplot([0 data_phase_mean{2,ii}(2,tt)],[0 data_PLV_mean{2,ii}(2,tt)]);
        set(I,'linewidth',6,'color',color)
        thetaticks([0 90 180 270])
        title({['\fontsize{14} PLV = ', num2str(data_PLV_mean{2,ii}(2,tt))];[]},'color',color)

        % mPFC IL <--> dHPC
        subplot(3,2,tt+5)
        polarplot([zeros(1,size(data_phase{2,ii}(3,1:240:end,tt),2)) data_phase{2,ii}(3,1:240:end,tt)]',repmat([0 1],1,size(data_phase{2,ii}(3,1:240:end,tt),2))','color',[.6 .6 .6]);
        hold on
        I=polarplot([0 data_phase_mean{2,ii}(3,tt)],[0 data_PLV_mean{2,ii}(3,tt)]);
        set(I,'linewidth',6,'color',color)
        thetaticks([0 90 180 270])
        title({['\fontsize{14} PLV = ', num2str(data_PLV_mean{2,ii}(3,tt))];[]},'color',color)

    end
end

%clear('ii','I','data_phase','data_phase_mean','data_PLV_mean')

%% Save Figures

newStr = id(1:end-8);

path = '/Users/flavio/Desktop';
%path = files.FilesLoaded{1, 1}.folder;

name = strcat(path,'/',newStr,'_PLV_retrieval');

%get all figures
figlist = findobj(allchild(0), 'flat', 'Type', 'figure');

name_figs = {'_3_5Hz','_6_8Hz','_7_10Hz'};
set(gcf,'renderer','Painters')

% Loop through figure
for ii = 1:numel(figlist)
    name_loop = strcat(name,name_figs{ii});
    FigHandle = figlist(ii);
    saveas(FigHandle,name_loop,'png')
    exportgraphics(gcf,strcat(name_loop,'.eps'),'Resolution', 300)

end

close all
clear('name','newStr','path')

%% Delta phase from Euler representation of angles. Time bins NO OVERLAP
% % Extracts relative phase and length of circular variance (PLV)
%
% time_bins     = 0.02; % seconds
% time_bins_idx = time_bins * parameters.decimated_srate;
%
%
% hilb.delta_phase_bins = cell(size(hilb.phase));
% hilb.euler_differences_bins = cell(size(hilb.phase));
% hilb.PLV_bins = cell(size(hilb.phase));
% hilb.mean_phase_bins = cell(size(hilb.phase));
%
% for ii = 1:size(hilb.phase,2)
%     for hh = 1:size(hilb.phase,1)
%         for jj = 1:size(hilb.phase{hh,ii},1)
%
%             bins = 2:time_bins_idx:size(hilb.phase{hh,ii},2);
%
%             for bb = 2:size(bins,2)
%
%                 if hh == 1
%                     hilb.delta_phase_bins{hh,ii}(jj,bb-1) = circ_mean(hilb.delta_phase{hh,ii}(jj,bins(bb-1):bins(bb)),[],2);
%                 else
%                     hilb.delta_phase_bins{hh,ii}(jj,bb-1,:) = circ_mean(hilb.delta_phase{hh,ii}(jj,bins(bb-1):bins(bb),:),[],2);
%                 end
%
%             end
%
%
%             if hh == 1
%                 hilb.euler_differences_bins{hh,ii}(jj,:) = exp(1i*(hilb.delta_phase_bins{hh,ii}(jj,:)));
%                 hilb.PLV_bins{hh,ii}(jj,:)               = abs(mean(hilb.euler_differences_bins{hh,ii}(jj,:),2));
%                 hilb.mean_phase_bins{hh,ii}(jj,:)        = angle(mean(hilb.euler_differences_bins{hh,ii}(jj,:),2));
%             else
%                 hilb.euler_differences_bins{hh,ii}(jj,:,:) = exp(1i*(hilb.delta_phase_bins{hh,ii}(jj,:,:)));
%                 hilb.PLV_bins{hh,ii}(jj,:,:)               = abs(mean(hilb.euler_differences_bins{hh,ii}(jj,:,:),2));
%                 hilb.mean_phase_bins{hh,ii}(jj,:,:)        = angle(mean(hilb.euler_differences_bins{hh,ii}(jj,:,:),2));
%             end
%
%         end
%     end
% end
%
% clear('time_bins','time_bins_idx','hh','jj','ii','bb')


%% Amplitude Sliding window. Time bins and overlap.
% Extracts relative amplitude

% Define time bins
time_bins     = 3; % seconds
time_bins_idx = time_bins * parameters.decimated_srate;

% Overlap
timeoverlap    = .95; % percentage
overlap = round((time_bins_idx)-(timeoverlap*time_bins_idx));

% All time series
hilb.amplitude_bins = cell(size(hilb.phase));

% Average to get trials
hilb.amplitude_mean_trial_bins   = cell(size(hilb.phase));

% Average to get session
hilb.amplitude_mean_session_bins = cell(size(hilb.phase));


for ii = 1:size(hilb.phase,2)
    for hh = 1:size(hilb.phase,1)
        for jj = 1:size(hilb.phase{hh,ii},1)

            bins = (2:overlap:size(hilb.phase{hh,ii},2)-time_bins_idx);

            % loop and overlap over time bins
            for bb = 2:size(bins,2)

                if hh == 1 % condition for baseline
                    hilb.amplitude_bins{hh,ii}(jj,bb-1)   = mean(hilb.amplitude{hh,ii}(jj,bins(bb-1):bins(bb-1) + time_bins_idx -1),2);

                else % condition for 4 Hz and 8 Hz
                    hilb.amplitude_bins{hh,ii}(jj,bb-1,:) = mean(hilb.amplitude{hh,ii}(jj,bins(bb-1):bins(bb-1) + time_bins_idx -1,:),2);

                end
            end

            if hh == 1 % condition for baseline

                hilb.amplitude_bins{hh,ii}(jj,:)                  = normalize(hilb.amplitude_bins{hh,ii}(jj,:),'range');
                hilb.amplitude_mean_trial_bins{hh,ii}(jj,:)       = mean(hilb.amplitude_bins{hh,ii}(jj,:),2);
                hilb.amplitude_mean_session_bins{hh,ii}(jj,:)     = mean(hilb.amplitude_mean_trial_bins{hh,ii}(jj,:),2);

            else % condition for 4 Hz and 8 Hz

                hilb.amplitude_bins{hh,ii}(jj,:,:)                = normalize(hilb.amplitude_bins{hh,ii}(jj,:,:),'range');
                hilb.amplitude_mean_trial_bins{hh,ii}(jj,:)       = mean(hilb.amplitude_bins{hh,ii}(jj,:,:),2);
                hilb.amplitude_mean_session_bins{hh,ii}(jj,:)     = mean(mean(hilb.amplitude_mean_trial_bins{hh,ii}(jj,:,:),2),3);

            end

        end
    end
end

hilb.parameters.timev_bins = linspace(-hilb.parameters.time_pre, 10 + hilb.parameters.time_pre, size(hilb.amplitude_bins{4,1},2));


clear('time_bins','time_bins_idx','hh','jj','ii','bb','timeoverlap','overlap','bins')

%% Save
% newStr = id(1:end-8);
% 
% path = '/Users/flavio/Desktop';
% %path = files.FilesLoaded{1, 1}.folder;
% 
% name = strcat(path,'/',newStr,'_PLV_extinction_Trials_TimeCourse');
% 
% %get all figures
% figlist = findobj(allchild(0), 'flat', 'Type', 'figure');
% 
% name_figs = {'_3_5Hz','_6_8Hz','_8_10Hz'};
% 
% % Loop through figure
% for ii = 1:numel(figlist)
%     name_loop = strcat(name,name_figs{ii});
%     FigHandle = figlist(ii);
%     saveas(FigHandle,name_loop,'png')
% end
% 
% close all
% clear('name','newStr','path')

%% Save

%newStr = regexprep(files.id.name,'.mat','_');
%newStr = files.id(ms).name(1:end-8);
newStr = id(1:end-8);

path = '/Users/flavio/Desktop';
%path = files.FilesLoaded{1, 1}.folder;

name = strcat(path,'/',newStr,'_PLV_retrieval');

% save data
save(name,'hilb','-v7.3')

%% last update 23/05/2024 - 17:00
%  listening:
