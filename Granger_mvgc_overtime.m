%% Granger Causality 

% 1) MVGC Multivariate Granger Causality MATLAB toolbox
% hosted at http://www.sussex.ac.uk/sackler/mvgc
% Current version is mvgc_v1.3, last updated March 2022

% L. Barnett and A. K. Seth, Granger causality for state-space models, Phys. Rev. E 91(4) Rapid Communication, 2015.
% L. Barnett and A. K. Seth, "The MVGC Multivariate Granger Causality Toolbox: A new approach to Granger-causal inference", J. Neurosci. Methods 223, pp 50-68, 2014.

% ----------------------------------

% 2) Time Spectrum and Time vector.
% VAR model estimation regression mode (LWR) by Mike X Cohen (mikexcohen@gmail.com)

% - The code relies on the following functions : - grangerX.m 
%                                                - BSMART ToolBox (http://www.sahs.uth.tmc.edu/hliang/software)


% ----------------------------------

% Part of this code was adapted from "mvgc_demo_statespace.m"
% by Flavio Mourao.
% email: mourao.fg@tamu.edu
% Maren Lab - Department of Psychological and Brain Sciences
% Texas A&M University

% Started in:  01/2024
% Last update: 02/2024

%% Start ToolBox
startup_mvgc1 % path /Users/flavio/Documents/MATLAB/MVGC1-master
clc

clear('have_genvar_mex','mvgc_root','mvgc_version')

%% Selected data


data.lfp{9,1} = [];

tpre = 10 * parameters.decimated_srate;
tpos = 10 * parameters.decimated_srate;

events = data.events{2, 1}(CSIT{ms},:);

for ii = 1:size(events,1)
    data.lfp{9,1}(:,:,ii) = data.lfp{5,1}(:,events(ii,1) - tpre : events(ii,2) + tpos);
end

%%
mvgc = [];


% The data will be downsampled to 250 and z-scored.
mvgc.parameters.decimate_factor = 10; % 
mvgc.data = [];



for ii = 1:size(data.lfp{9,1},1)
    for jj = 1:size(data.lfp{9,1},3)
        mvgc.data(ii,:,jj) = zscore(decimate(data.lfp{9,1}(ii,:,jj),mvgc.parameters.decimate_factor));
    end
end



clear('ii','jj')

%% Granger over time - Spectrum and Time vector.
%  VAR model estimation regression mode (LWR) by Mike X Cohen (mikexcohen@gmail.com)

% - The code relies on the following functions : - grangerX.m 
%                                                - BSMART ToolBox (http://www.sahs.uth.tmc.edu/hliang/software)



tic

fprintf('Granger over time - Spectrum and Time vector. VAR model estimation regression mode (LWR) by Mike X Cohen\n')

% Choose channels combinations
Combinations_ = nchoosek(1:size(mvgc.data,1),2); % all possible combinations

% Define band frequency
% 2 - 12 Hertz
mvgc.parameters.frex_idx_2_12Hz = linspace(2,12,length(short_fft.freq2plot1));

win = 2048; % in seconds.
mvgc.parameters.time_winx  = round(win/(1000/100)); % in samples


% Baseline - Comment out this for-loop and make the changes in the following for-loop as suggested above to calculate everything simultaneously.
for cc = 1:size(Combinations_,1)

    [mvgc.Time_x2y_Spec{1,1}(:,:,cc),mvgc.Time_y2x_Spec{1,1}(:,:,cc),mvgc.Time_x2y_FInt{1,1}(:,:,cc),mvgc.Time_y2x_FInt{1,1}(:,:,cc)] = ...
        GrangerX(mvgc.data(Combinations_(cc,:),:),100,mvgc.parameters.frex_idx_2_12Hz,mvgc.parameters.time_winx,30);
    
end


clear('tt','ii','cc','tt','steps ','win')

fprintf('Done.')

toc

%% PLot Granger Time Spec. 

Session = 2;

Combinations_1 = nchoosek(1:size(mvgc.data,1),2); % all possible combinations
Combinations_2 = flip(Combinations_1,2);

% Define band frequency
% 2 - 12 Hertz
steps         = diff(mvgc.parameters.freqs); % according to the fft time window
mvgc.parameters.frex_2_12Hz     = 2:steps(1):12;
mvgc.parameters.frex_idx_2_12Hz = dsearchn(mvgc.parameters.freqs,mvgc.parameters.frex_2_12Hz');

% Time window
win = 1000; % in seconds.
mvgc.parameters.time_winx  = round(win/(1000/mvgc.parameters.fs)); % in samples


% Time vector_Granger
time_v_spec = linspace(0,10,size(mvgc.Time_x2y_Spec{Session,1},2));

% Behavior data
data_2_plot_Behav =[];
for jj = 1:length(CSIT{ms})
    data_2_plot_Behav(:,:,jj) = decimate(data.behavior{3,CSIT{ms}(jj)*2}(1,1:end-1),4);
end

% Normalize behavior vector to the Granger time window
for jj = 1:size(data_2_plot_Behav,3)
    for ii = 1:length(data_2_plot_Behav) - mvgc.parameters.time_winx
        data_2_plot_Behav_m(1,ii,jj) = mean(data_2_plot_Behav(1,ii:ii+mvgc.parameters.time_winx,jj),2);
    end
end

% Time vector_Behavior
time_v_Behav = linspace(0,10,size(data_2_plot_Behav_m,2));


%sgtitle({'Amplitude Spectrum via short-window FFT';['Hamming window = ' num2str(short_fft.timewin./1000) 's' ' - ' 'overlap = ' num2str(short_fft.overlap) '%']})


for jj = 1:length(CSIT{ms})
    
    figure 
    set(gcf, 'Position', get(0, 'Screensize'));
    set(gcf,'color','w');

    for cc = 1:size(Combinations_1,1)
        
        subplot(3,size(Combinations_1,1),cc)
        contourf(time_v_spec, mvgc.parameters.freqs(mvgc.parameters.frex_idx_2_12Hz),smoothdata(mvgc.Time_x2y_Spec{Session,jj}(:,:,cc),2,'gaussian',35),150,'linecolor','none');
        xlabel('Time (s)','FontSize',14), ylabel('Frequency (Hz)','FontSize',14)
        clim([0 3])


        if cc == 1
            title('mPFC PL --> mPFC IL')
        elseif cc == 2
            title('mPFC PL --> HPC')
        else
            title('mPFC IL --> HPC')
            c = colorbar;
            c.Label.String = 'Granger (A.U.)';
        end


        subplot(3,size(Combinations_1,1),cc+size(Combinations_1,1))
        contourf(time_v_spec, mvgc.parameters.freqs(mvgc.parameters.frex_idx_2_12Hz),smoothdata(mvgc.Time_y2x_Spec{Session,jj}(:,:,cc),2,'gaussian',35),150,'linecolor','none');
        xlabel('Time (s)','FontSize',14), ylabel('Frequency (Hz)','FontSize',14)
        clim([0 3])


        if cc == 1
            title('mPFC IL --> mPFC PL')
        elseif cc == 2
            title('HPC --> mPFC PL')
        else
            title('HPC --> mPFC IL')
            c = colorbar;
            c.Label.String = 'Granger (A.U.)';
        end


        subplot(3,size(Combinations_1,1),cc+size(Combinations_1,1)+3)
        plot(time_v_Behav,data_2_plot_Behav_m(:,:,jj),'linew',2,'color',[0.7 0.7 0.7])
        set(gca,'ytick',[])
        set(gca,'ycolor',[1 1 1])
        box off
        hold on
        ylim([0 100])
        
        bb = data_2_plot_Behav_m(:,:,jj);
        ff = (data_2_plot_Behav_m(:,:,jj) < 5);

        plot(time_v_Behav(ff), bb(ff) ,'k.','linew', 2, 'color', [.6 0 0])
        xlabel('Time (s)','FontSize',14)
        legend('movement','freezing')
    
    end
,
% Save
newStr1 = id(1:end-20);
Path    = files.FilesLoaded{1,1}(ms).folder;
name_1  = strcat(Path,'/',newStr1,'_Granger_over_time_CS_Trial_',num2str(CSIT{ms}(jj)));

saveas(gcf,name_1,'png') % save figure

close all

end

clear('ii','ch','steps','win','newStr1','name_1','bb','ff','c','cc','jj','timev','time_v_Behav','data_2_plot_Behav','data_2_plot_Behav_m','time_v_spec','win','steps','Session')
%% Save data

% Settings
%ms = 1;
%Path    = files.FilesLoaded{1,1}(ms).folder;
Path = '/Users/flavio/Desktop';

newStr1 = id(1:end-8);
name_1 = strcat(Path,'/',newStr1,'_MVGC_Granger_Context');

% Save data
save(name_1,'mvgc','-v7.3')

clear('name','newStr1','path') 



%% last update 13/10/2024
%  listening:
