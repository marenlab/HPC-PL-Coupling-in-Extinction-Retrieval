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

%ms=1

% Extinction ans Retrieval Sessions
temp_data{1,1} = data.lfp{6,1}(:,B_clean{ms}(1):B_clean{ms}(2)-1); % baseline
temp_data{2,1} = data.lfp{7,1}(:,:,CSIT{ms}(1,:)); % CS or Time win Contextual Exposure
temp_data{3,1} = data.lfp{8,1}(:,:,CSIT{ms}(1,:)); % ITI or Time win Contextual Exposure. Be carefful with the index here, usualliy ITI is idx 8. 7 will be just for context exposure

% temp_data{2,2} = data.lfp{7,2}(:,:,CSIT{ms}(2,:)); % CS or Time win Contextual Exposure
% temp_data{3,2} = data.lfp{8,2}(:,:,CSIT{ms}(2,:)); % ITI or Time win Contextual Exposure. Be carefful with the index here, usualliy ITI is idx 8. 7 will be just for context exposure
% 

mvgc = [];


% The data will be downsampled to 250 and z-scored.
mvgc.parameters.decimate_factor = 4; % 
mvgc.data = [];

for tt = 1:size(temp_data,1)
    for tt_1 = 1:size(temp_data,2)

        if isempty(temp_data{tt,tt_1})
            continue
        end

        if tt == 1 % case for baseline period

            for ii = 1:size(temp_data{tt,tt_1},1)
                for jj = 1:size(temp_data{tt,tt_1},3)
                    mvgc.data{tt,tt_1}(ii,:,jj) = normalize(decimate(temp_data{tt,tt_1}(ii,:,jj),mvgc.parameters.decimate_factor),2,'zscore');
                end
            end

        else

            for ii = 1:size(temp_data{tt,tt_1},1)
                for jj = 1:size(temp_data{tt,tt_1},3)
                    mvgc.data{tt,tt_1}(ii,:,jj) = normalize(decimate(temp_data{tt,tt_1}(ii,:,jj),mvgc.parameters.decimate_factor),2,'zscore');
                end
            end
        end
    end
end


clear('temp_data','tt','ii','jj')

%% Parameters

mvgc.parameters.regmode   = 'LWR';                                                            % VAR model estimation regression mode ('OLS', 'LWR' or empty for default)
mvgc.parameters.icregmode = 'LWR';                                                            % information criteria regression mode ('OLS', 'LWR' or empty for default)
mvgc.parameters.number_channels = 3;

mvgc.parameters.morder    = 'BIC';                                                          % model order to use ('actual', 'AIC', 'BIC' or supplied numerical value)
mvgc.parameters.momax     = 55;                                                               % maximum model order for model order estimation

mvgc.parameters.acmaxlags = [];                                                               % maximum autocovariance lags (empty for automatic calculation)

mvgc.parameters.tstat     = 'F';                                                              % statistical test for MVGC:  'F' for Granger's F-test (default) or 'chi2' for Geweke's chi2 test
mvgc.parameters.alpha     = 0.05;                                                             % significance level for significance test
mvgc.parameters.mhtc      = 'FDRD';                                                           % multiple hypothesis test correction (see routine 'significance')

mvgc.parameters.fs        = parameters.decimated_srate./mvgc.parameters.decimate_factor;      % sample rate (Hz)
mvgc.parameters.fres      = 4096;                                                             % frequency resolution (empty for automatic calculation)

mvgc.parameters.freqs     = sfreqs(mvgc.parameters.fres,mvgc.parameters.fs);                  % frequency vector based on frequency resolution

%% Model order estimation (<mvgc_schema.html#3 |A2|>)
% The model's order will be defined only once, based on the baseline of
% each animal, or defined cell

data_2_estimate = mvgc.data{1,1};

% Calculate information criteria up to specified maximum model order.

ptic('\n*** tsdata_to_infocrit\n');
[mvgc.parameters.AIC,mvgc.parameters.BIC,mvgc.parameters.moAIC,mvgc.parameters.moBIC] = tsdata_to_infocrit(data_2_estimate,mvgc.parameters.momax,mvgc.parameters.icregmode);
ptoc('*** tsdata_to_infocrit took ');

% Plot information criteria.

mvgc.parameters.amo = 5; % actual model order

fprintf('\nbest model order (AIC) = %d\n',mvgc.parameters.moAIC);
fprintf('best model order (BIC) = %d\n',mvgc.parameters.moBIC);
fprintf('actual model order     = %d\n',mvgc.parameters.amo);

% Select model order.

if     strcmpi(mvgc.parameters.morder,'actual')
       mvgc.parameters.morder = mvgc.parameters.amo;
       fprintf('\nusing actual model order = %d\n',mvgc.parameters.morder);

       figure;
       plot_tsdata([mvgc.parameters.AIC mvgc.parameters.BIC]',{'AIC','BIC'},1/mvgc.parameters.fs);
       title(['Model order estimation - using the best actual model order = ' num2str(mvgc.parameters.morder)]);

elseif strcmpi(mvgc.parameters.morder,'AIC')
       mvgc.parameters.morder = mvgc.parameters.moAIC;
       fprintf('\nusing AIC best model order = %d\n',mvgc.parameters.morder);

       figure;
       plot_tsdata([mvgc.parameters.AIC mvgc.parameters.BIC]',{'AIC','BIC'},1/mvgc.parameters.fs);
       title(['Model order estimation - using AIC best model order = ' num2str(mvgc.parameters.morder)]);
       xline(mvgc.parameters.morder/ mvgc.parameters.fs, '--r', ...
       ['\leftarrow Order = ' num2str(mvgc.parameters.morder)], ...
       'LabelOrientation', 'horizontal', 'FontSize', 10);

       xline(25/ mvgc.parameters.fs, '--r', ...
       ['\leftarrow Order = ' num2str(25)], ...
       'LabelOrientation', 'horizontal', 'FontSize', 10);


elseif strcmpi(mvgc.parameters.morder,'BIC')
       mvgc.parameters.morder = mvgc.parameters.moBIC;
       fprintf('\nusing BIC best model order = %d\n',mvgc.parameters.morder);

       figure;
       plot_tsdata([mvgc.parameters.AIC mvgc.parameters.BIC]',{'AIC','BIC'},1/mvgc.parameters.fs);
       title(['Model order estimation - using BIC best model order = ' num2str(mvgc.parameters.morder)]);
       xline(mvgc.parameters.morder/ mvgc.parameters.fs, '--r', ...
       ['\leftarrow Order = ' num2str(mvgc.parameters.morder)], ...
       'LabelOrientation', 'horizontal', 'FontSize', 10);

       xline(25/ mvgc.parameters.fs, '--r', ...
       ['\leftarrow Order = ' num2str(25)], ...
       'LabelOrientation', 'horizontal', 'FontSize', 10);


else
       fprintf('\nusing specified model order = %d\n',mvgc.parameters.morder);

       figure;
       plot_tsdata([mvgc.parameters.AIC mvgc.parameters.BIC]',{'AIC','BIC'},1/mvgc.parameters.fs);
       title(['Model order estimation - using specified best model order = ' num2str(mvgc.parameters.morder)]);
       xline(mvgc.parameters.morder/ mvgc.parameters.fs, '--r', ...
       ['\leftarrow Order = ' num2str(mvgc.parameters.morder)], ...
       'LabelOrientation', 'horizontal', 'FontSize', 10);

       xline(25/ mvgc.parameters.fs, '--r', ...
       ['\leftarrow Order = ' num2str(25)], ...
       'LabelOrientation', 'horizontal', 'FontSize', 10);


end

clear('data_2_estimate');

%% Save Granger model estimation figure

% Settings
% ms = 1;
newStr1 = id(1:end-8);
%Path        = files.FilesLoaded{1,1}(ms).folder;
Path = '/Users/flavio/Desktop';

name_1 = strcat(Path,'/',newStr1,'_Granger_model_estimation');

% save figure
saveas(gcf,name_1,'png')

close all

clear('name_1','newStr1','path') 

%% VAR model estimation (<mvgc_schema.html#3 |A2|>)
% % each trial
% 
% % Estimate VAR model of selected order from data.
% 
% for tt = 1:size(mvgc.data,1)
%     for tt_1 = 1:size(temp_data,2)
% 
%         if isempty(mvgc.data{tt,tt_1})
%             continue
%         end
% 
%         for ii = 1:size(mvgc.data{tt,tt_1},3)
% 
%             ptic('\n*** tsdata_to_var... ');
%             [mvgc.parameters.A{tt,tt_1},mvgc.parameters.SIG{tt,tt_1}] = tsdata_to_var(mvgc.data{tt,tt_1}(:,:,ii),mvgc.parameters.morder,[]);
%             assert(~isbad(mvgc.parameters.A{tt,tt_1}),'VAR estimation failed - bailing out');
%             ptoc;
% 
%             % Report information on the estimated VAR, and check for errors.
%             %
%             % _IMPORTANT:_ We check the VAR model for stability and symmetric
%             % positive-definite residuals covariance matrix. _THIS CHECK SHOULD ALWAYS BE
%             % PERFORMED!_ - subsequent routines may fail if there are errors here. If there
%             % are problems with the data (e.g. non-stationarity, colinearity, etc.) there's
%             % also a good chance they'll show up at this point - and the diagnostics may
%             % supply useful information as to what went wrong.
% 
%             mvgc.info{tt,tt_1} = var_info(mvgc.parameters.A{tt,tt_1},mvgc.parameters.SIG{tt,tt_1});
%             assert(~mvgc.info{tt,tt_1}.error,'VAR error(s) found - bailing out');
% 
%         end
%     end
% end
% 
% clear('tt','ii')

%% VAR model estimation (<mvgc_schema.html#3 |A2|>)
% Full model

% Estimate VAR model of selected order from data.

for tt = 1:size(mvgc.data,1)
    for tt_1 = 1:size(mvgc.data,2)

        if isempty(mvgc.data{tt,tt_1})
            continue
        end

        ptic('\n*** tsdata_to_var... ');
        [mvgc.parameters.A{tt,tt_1},mvgc.parameters.SIG{tt,tt_1}] = tsdata_to_var(mvgc.data{tt,tt_1},mvgc.parameters.morder,[]);
        assert(~isbad(mvgc.parameters.A{tt,1}),'VAR estimation failed - bailing out');
        ptoc;

        % Report information on the estimated VAR, and check for errors.
        %
        % _IMPORTANT:_ We check the VAR model for stability and symmetric
        % positive-definite residuals covariance matrix. _THIS CHECK SHOULD ALWAYS BE
        % PERFORMED!_ - subsequent routines may fail if there are errors here. If there
        % are problems with the data (e.g. non-stationarity, colinearity, etc.) there's
        % also a good chance they'll show up at this point - and the diagnostics may
        % supply useful information as to what went wrong.

        mvgc.info{tt,tt_1} = var_info(mvgc.parameters.A{tt,tt_1},mvgc.parameters.SIG{tt,tt_1});
        assert(~mvgc.info{tt,1}.error,'VAR error(s) found - bailing out');

    end
end

clear('tt','ii')


%% Granger causality calculation: frequency domain  (<mvgc_schema.html#3 |A14|>)

% If not specified, we set the frequency resolution to something sensible. Warn if
% resolution is very large, as this may lead to excessively long computation times,
% and/or out-of-memory issues.

if isempty(mvgc.parameters.fres)
    mvgc.parameters.fres = 2^nextpow2(mvgc.info.acdec); % based on autocorrelation decay; alternatively, you could try fres = 2^nextpow2(nobs);
	fprintf('\nfrequency resolution auto-calculated as %d (increments ~ %.2gHz)\n',mvgc.parameters.fres,mvgc.parameters.fs/2/mvgc.parameters.fres);
end

if mvgc.parameters.fres > 20000 % adjust to taste
	fprintf(2,'\nWARNING: large frequency resolution = %d - may cause computation time/memory usage problems\nAre you sure you wish to continue [y/n]? ',mvgc.parameters.fres);
	istr = input(' ','s'); if isempty(istr) || ~strcmpi(istr,'y'); fprintf(2,'Aborting...\n'); return; end
end

% Calculate spectral pairwise-conditional causalities at given frequency
% resolution by state-space method.

for tt = 1:size(mvgc.data,1)
    for tt_1 = 1:size(mvgc.data,2)

        if isempty(mvgc.data{tt,tt_1})
            continue
        end

        %for ii = 1:size(mvgc.data{tt,tt_1},3) % in case to perform for each trial...need to add index below

            ptic('\n*** var_to_spwcgc... ');
            mvgc.F_spect{tt,tt_1} = var_to_spwcgc(mvgc.parameters.A{tt,tt_1},mvgc.parameters.SIG{tt,tt_1},mvgc.parameters.fres);
            assert(~isbad(mvgc.F_spect{tt,tt_1},false),'spectral GC calculation failed - bailing out');
            ptoc;

        %end
    end
end

clear('tt','ii')

%% Plot
% % Plot spectral causal graph.
% % Plot pairwise spectral quantities in |P|, a 3-dim numerical matrix with
% % first index representing target ("to"), second index source ("from")
% % and third index frequency range - typically spectral causalities
% 
% % Considering all trials/time epochs. Decide this for freezing analysis
% CSIT{ms} = 1:1:length(mvgc.F_spect);
% 
% % ms = 1;
% 
% % Baseline
% figure;
% set(gcf, 'Position', get(0, 'Screensize'));
% 
% %sgtitlex('Baseline - Pairwise-conditional Granger causality - frequency domain.');
% plot_spw(mvgc.F_spect{1,1},mvgc.parameters.fs,[2 12]);
% 
% newStr1 = id(1:end-8);
% Path = '/Users/flavio/Desktop';
% name_1 = strcat(Path,'/',newStr1,'_Granger_Spec_baseline');
% saveas(gcf,name_1,'png') % save figure
% 
% close all
% 
% 
% % CS and/or ITI trials
% %CSIT{ms} = [1:1:3]; % Trial number
% session = 2; % Type 2 to CS and 3 to ITI
% 
%  data2plot = [];
% for ii = 1:length(CSIT{ms})
%     if isempty(mvgc.F_spect{session,CSIT{ms}(ii)})
%         continue
%     end
% 
%     data2plot(:,:,:,ii) = mvgc.F_spect{session,CSIT{ms}(ii)};
% 
% end
% 
% % All chosed Trials - saving separately
% % for jj = 1:length(CSIT{ms})
% %     figure;
% %     set(gcf, 'Position', get(0, 'Screensize'));
% % 
% %     sgtitlex(['CS - Pairwise-conditional Granger causality - frequency domain. Trial = ' num2str(CSIT{ms}(jj))]);
% %     plot_spw(data2plot(:,:,:,jj),mvgc.parameters.fs,[2 12]);
% %     
% %     newStr1 = id(1:end-8);
% %     Path    = files.FilesLoaded{1,1}(ms).folder;
% %     name_1 = strcat(Path,'/',newStr1,'_Granger_Spec_CS_trial_',num2str(CSIT{ms}(jj)));
% %     saveas(gcf,name_1,'png') % save figure
% % 
% %     close all
% % 
% % end
% 
% % Mean Trials
% figure;
% set(gcf, 'Position', get(0, 'Screensize'));
% 
% if isempty(data2plot)
%     data2plot = ones(size(mvgc.F_spect{3,1}));
% end
% 
% %sgtitlex('CS - Pairwise-conditional Granger causality - frequency domain. Mean trials');
% sgtitlex('Freezing - Pairwise-conditional Granger causality - frequency domain. Mean trials');
% 
% plot_spw(mean(data2plot,4),mvgc.parameters.fs,[2 12]);
% 
% % Save
% newStr1 = id(1:end-8);
% %Path    = files.FilesLoaded{1,1}(ms).folder;
% Path = '/Users/flavio/Desktop';
% 
% %name_1 = strcat(Path,'/',newStr1,'_Granger_Spec_Mean_CS_Trials_1');
% % name_1 = strcat(Path,'/',newStr1,'_Granger_Spec_Mean_Freezing');
% name_1 = strcat(Path,'/',newStr1,'_Granger_Spec_Mean_Context');
% 
% saveas(gcf,name_1,'png') % save figure
% 
% %close all
% 
% 
% 
% session = 3; % Type 2 to CS and 3 to ITI
% 
% data2plot = [];
% for ii = 1:length(CSIT{ms})
%     if isempty(mvgc.F_spect{session,CSIT{ms}(ii)})
%         continue
%     end
% 
%     data2plot(:,:,:,ii) = mvgc.F_spect{session,CSIT{ms}(ii)};
% 
% end
% 
% % Mean Trials
% figure;
% set(gcf, 'Position', get(0, 'Screensize'));
% 
% if isempty(data2plot)
%     data2plot = ones(size(mvgc.F_spect{3,1}));
% end
% 
% %sgtitlex('CS - Pairwise-conditional Granger causality - frequency domain. Mean trials');
% sgtitlex('non-Freezing - Pairwise-conditional Granger causality - frequency domain. Mean trials');
% 
% plot_spw(mean(data2plot,4),mvgc.parameters.fs,[2 12]);
% 
% % Save
% newStr1 = id(1:end-8);
% %Path    = files.FilesLoaded{1,1}(ms).folder;
% Path    = '/Users/flavio/Desktop';
% 
% %name_1 = strcat(Path,'/',newStr1,'_Granger_Spec_Mean_CS_Trials_1');
% %name_1 = strcat(Path,'/',newStr1,'_Granger_Spec_Mean_non_Freezing');
% name_1 = strcat(Path,'/',newStr1,'_Granger_Spec_Mean_Context');
% 
% saveas(gcf,name_1,'png') % save figure
% 
% %close all
% 
% clear('ii','jj')
% clear('session','name_1','newStr1','path','data2plot') 


%% Plot
% % Plot spectral causal graph. Version 2.
% % Plot pairwise spectral quantities in |P|, a 3-dim numerical matrix with
% % first index representing target ("to"), second index source ("from")
% % and third index frequency range - typically spectral causalities
% 
% % Define parameters
% % all possible combinations
% Combinations_2 = nchoosek(1:size(mvgc.data{1,1},1),2);
% Combinations_1 = flip(Combinations_2,2); % all possible combinations
% 
% %CSIT{ms} = [1:1:10];
% sub_idx = reshape((1:length(CSIT{ms})*3),[],3);
% 
% steps                          = diff(mvgc.parameters.freqs); % according to the fft time window
% mvgc.parameters.frex_2_12Hz     = 0:steps(1):12;
% mvgc.parameters.frex_idx_2_12Hz = dsearchn(mvgc.parameters.freqs,mvgc.parameters.frex_2_12Hz');
% 
% 
% 
% 
% figure
% set(gcf, 'Position', get(0, 'Screensize'));
% %sgtitlex('Baseline - Pairwise-conditional Granger causality - frequency domain.');
% %sgtitlex('Pairwise-conditional Granger causality - frequency domain.');
% 
% % Baseline
% for cc = 1:size(Combinations_1,1)
% 
%         subplot(4,3,cc)
%         plot(mvgc.parameters.freqs(mvgc.parameters.frex_idx_2_12Hz),squeeze(mvgc.F_spect{1,1}(Combinations_1(cc,1),Combinations_1(cc,2),mvgc.parameters.frex_idx_2_12Hz)))
%         hold
%         plot(mvgc.parameters.freqs(mvgc.parameters.frex_idx_2_12Hz),squeeze(mvgc.F_spect{1,1}(Combinations_2(cc,1),Combinations_2(cc,2),mvgc.parameters.frex_idx_2_12Hz)))
% 
%         if          cc == 1
%                     legend('PL --> IL', 'IL --> PL')
%                     ylabel({'Baseline';'Granger causality'},'FontWeight','bold','FontSize',14)
% 
%         elseif      cc == 2
%                     legend('PL --> HPC', 'HPC --> PL')
%         else
%                     legend('IL --> HPC', 'HPC --> IL')
%         end
% 
% xlabel('(Hz)','FontSize',12)
% xlim([2 12])
% %ylim([0 0.9])
% 
% end
% 
% 
% % Type 2 to CS or freezing and 3 to ITI or non-freezing
% session = 2; % Type 2 to CS and 3 to ITI
% 
% data2plot = [];
% 
% for ii = 1:length(CSIT{ms})
% 
%     if isempty(mvgc.F_spect{session,CSIT{ms}(ii)})
%         continue
%     end
% 
%     data2plot(:,:,:,ii) = mvgc.F_spect{session,CSIT{ms}(ii)};
% 
% end
% 
% data2plot_mean = mean(data2plot,4);
% 
% if isempty(data2plot_mean)
%     data2plot_mean = nan(size(mvgc.F_spect{3,1}));
% end
% 
% for cc = 1:size(Combinations_1,1)
% 
%         subplot(3,3,cc+3)
%         plot(mvgc.parameters.freqs(mvgc.parameters.frex_idx_2_12Hz),squeeze(data2plot_mean(Combinations_1(cc,1),Combinations_1(cc,2),mvgc.parameters.frex_idx_2_12Hz)))
%         hold
%         plot(mvgc.parameters.freqs(mvgc.parameters.frex_idx_2_12Hz),squeeze(data2plot_mean(Combinations_2(cc,1),Combinations_2(cc,2),mvgc.parameters.frex_idx_2_12Hz)))
% 
%         if          cc == 1
%                     legend('PL --> IL', 'IL --> PL')
%                     ylabel({'First Time Epochs';'Granger causality'},'FontWeight','bold','FontSize',14)
%         elseif      cc == 2
%                     legend('PL --> HPC', 'HPC --> PL')
%         else
%                     legend('IL --> HPC', 'HPC --> IL')
%         end
% 
% xlabel('(Hz)','FontSize',12)
% xlim([2 12])
% %ylim([0 0.9])
% end
% 
% 
% % Type 2 to CS or freezing and 3 to ITI or non-freezing
% session = 3; 
% 
% data2plot = [];
% 
% for ii = 1:length(CSIT{ms})
% 
%     if isempty(mvgc.F_spect{session,CSIT{ms}(ii)})
%         continue
%     end
% 
%     data2plot(:,:,:,ii) = mvgc.F_spect{session,CSIT{ms}(ii)};
% 
% end
% 
% data2plot_mean = mean(data2plot,4);
% 
% if isempty(data2plot_mean)
%     data2plot_mean = nan(size(mvgc.F_spect{3,1}));
% end
% 
% for cc = 1:size(Combinations_1,1)
% 
%         subplot(3,3,cc+6)
%         plot(mvgc.parameters.freqs(mvgc.parameters.frex_idx_2_12Hz),squeeze(data2plot_mean(Combinations_1(cc,1),Combinations_1(cc,2),mvgc.parameters.frex_idx_2_12Hz)))
%         hold
%         plot(mvgc.parameters.freqs(mvgc.parameters.frex_idx_2_12Hz),squeeze(data2plot_mean(Combinations_2(cc,1),Combinations_2(cc,2),mvgc.parameters.frex_idx_2_12Hz)))
% 
%         if          cc == 1
%                     legend('PL --> IL', 'IL --> PL')
%                     ylabel({'Last Time Epochs';'Granger causality'},'FontWeight','bold','FontSize',14)
% 
%         elseif      cc == 2
%                     legend('PL --> HPC', 'HPC --> PL')
%         else
%                     legend('IL --> HPC', 'HPC --> IL')
%         end
% 
% xlabel('(Hz)','FontSize',12)
% xlim([2 12])
% %ylim([0 0.9])
% 
% end
% 
% 
% % Save
% newStr1 = id(1:end-8);
% Path    = '/Users/flavio/Desktop';
% 
% %name_1  = strcat(Path,'/',newStr1,'_Granger_Spec_Mean_CS_Trials_2');
% %name_1  = strcat(Path,'/',newStr1,'_Granger_ALL_Spec_Mean_Freezing');
% name_1  = strcat(Path,'/',newStr1,'_Granger_ALL_Spec_Mean_Context');
% 
% %Path    = files.FilesLoaded{1,1}(ms).folder;
% 
% saveas(gcf,name_1,'png') % save figure
% 
% close all
% 
% clear('ii','jj')
% clear('Combinations_2','Combinations_1','sub_idx','steps','data2plot_mean','cc','session','name_1','newStr1','path','data2plot') 

%% Granger causality calculation: time domain  (<mvgc_schema.html#3 |A13|>)

% Calculate time-domain pairwise-conditional causalities from VAR model parameters
% by state-space method [4]. The VAR model is transformed into an equivalent state-
% space model for computation. Also return p-values for specified test (F-test or
% likelihood-ratio test; this is optional - if p-values are not required, then it
% is not necessary to supply time series |X|, regression mode |regmode|, or test
% specification |tstat|).

for tt = 1:size(mvgc.data,1)
    for tt_1 = 1:size(mvgc.data,2)

        if isempty(mvgc.data{tt,tt_1})
            continue
        end

        % for ii = 1:size(mvgc.data{tt,1},3)

            ptic('*** var_to_pwcgc... ');
            [mvgc.stats.F{tt,tt_1},mvgc.stats.pval{tt,tt_1}] = var_to_pwcgc(mvgc.parameters.A{tt,tt_1},mvgc.parameters.SIG{tt,tt_1},mvgc.data{tt,tt_1},mvgc.parameters.regmode,mvgc.parameters.tstat);
            ptoc;

        % end
    end
end

% Check for failed GC calculation
for tt = 1:size(mvgc.data,1)
    for tt_1 = 1:size(mvgc.data,2)

        if isempty(mvgc.data{tt,tt_1})
            continue
        end

        % for ii = 1:size(mvgc.data{tt,tt_1},3)

            assert(~isbad(mvgc.stats.F{tt,tt_1},false),'GC calculation failed - bailing out');

            % Significance-test p-values, correcting for multiple hypotheses.

            mvgc.stats.sig{tt,tt_1} = significance(mvgc.stats.pval{tt,tt_1},mvgc.parameters.alpha,mvgc.parameters.mhtc);

        % end
    end
end

clear('tt','ii')

%% Plot
% Plot time-domain causal graph, p-values and significance.

% ms = 1;

% % Baseline
% figure
% set(gcf, 'Position', get(0, 'Screensize'));
% 
% sgtitlex('Baseline - Pairwise-conditional Granger causality - time domain');
% 
% subplot(1,3,1);
% plot_pw(mvgc.stats.F{1,1});
% title('Pairwise-conditional GC');
% hb = colorbar('location','southoutside');
% 
% subplot(1,3,2);
% plot_pw(mvgc.stats.pval{1,1});
% title(['p-values (' mvgc.parameters.tstat   '-test)']);
% hb = colorbar('location','southoutside');
% 
% subplot(1,3,3);
% plot_pw(mvgc.stats.sig{1,1});
% title(['Significant at \alpha = ' num2str(mvgc.parameters.alpha)]);
% hb = colorbar('location','southoutside');
% 
% newStr1 = id(1:end-20);
% name_1 = strcat(Path,'/',newStr1,'_Granger_Stats_baseline');
% saveas(gcf,name_1,'png') % save figure
% close all


% CS or ITI trials
% CSIT{ms} = [1 2 4]; % Choose CS and or ITI
% Session = 2; % Type 2 to CS and 3 to ITI
% sub_idx = reshape((1:length(CSIT{ms})*3),[],3);
% 
% figure
% set(gcf, 'Position', get(0, 'Screensize'));
% 
% sgtitlex('Pairwise-conditional Granger causality - time domain - CS Trials');
% 
% 
% for ii = 1:size(sub_idx,1)
%     s = subplot(3,size(sub_idx,1),sub_idx(ii,1));
%     plot_pw(mvgc.stats.F{Session,CSIT{ms}(ii)});
%     title({['F-test'],['Trial: ' num2str(CSIT{ms}(ii))]});
%     %clim([0 1])
% 
%     if ii == size(sub_idx,1)
%         s2Pos = get(s,'position');
%         hb = colorbar('location','eastoutside');
%         set(s,'position',s2Pos);
%     end
% 
% end
% 
% for ii = 1:size(sub_idx,1)
%     s = subplot(3,size(sub_idx,1),sub_idx(ii,2));
%     plot_pw(mvgc.stats.pval{Session,CSIT{ms}(ii)});
%     title(['p-values (' mvgc.parameters.tstat   '-test)']);
% 
%     if ii == size(sub_idx,1)
%         s2Pos = get(s,'position');
%         hb = colorbar('location','eastoutside');
%         set(s,'position',s2Pos);
%     end
% 
% end
% 
% for ii = 1:size(sub_idx,1)
%     s = subplot(3,size(sub_idx,1),sub_idx(ii,3));
%     plot_pw(mvgc.stats.sig{Session,CSIT{ms}(ii)});
%     title(['Significant at \alpha = ' num2str(mvgc.parameters.alpha)]);
% 
%     if ii == size(sub_idx,1)
%         colorbar
%         s2Pos = get(s,'position');
%         hb = colorbar('location','eastoutside');
%         set(s,'position',s2Pos);
%     end
% 
% end
% 
% % Save
% newStr1 = id(1:end-8);
% Path    = files.FilesLoaded{1,1}(ms).folder;
% name_1  = strcat(Path,'/',newStr1,'_Granger_Stats_CSTrials');
% 
% saveas(gcf,name_1,'png') % save figure
% 
% close all
% 
% clear('ii','s','s2Pos','hb','sub_idx','name_1','newStr1')

%% Granger causality calculation: frequency domain -> time-domain  (<mvgc_schema.html#3 |A15|>)

% Check that spectral causalities average (integrate) to time-domain
% causalities, as they should according to theory.

% Define band frequency
% 3 - 6 Hertz
steps                          = diff(mvgc.parameters.freqs); % according to the fft time window
mvgc.parameters.frex_3_6Hz     = 3:steps(1):6;
mvgc.parameters.frex_idx_3_6Hz = dsearchn(mvgc.parameters.freqs,mvgc.parameters.frex_3_6Hz');

for tt = 1:size(mvgc.data,1)
    for tt_1 = 1:size(mvgc.data,2)

        if isempty(mvgc.data{tt,tt_1})
            continue
        end

        %for ii = 1:size(mvgc.data{tt,tt_1},3)


            fprintf('\nfrequency-domain GC integration check... ');
            mvgc.Fint_3_6Hz{tt,tt_1} = smvgc_to_mvgc(mvgc.F_spect{tt,tt_1}(:,:,mvgc.parameters.frex_idx_3_6Hz)); % integrate spectral MVGCs


            mvgc.parameters.amax = maxabs(mvgc.stats.F{tt,tt_1}+mvgc.Fint_3_6Hz{tt,tt_1})/2;

            if mvgc.parameters.amax < 1e-5; mvgc.parameters.amax = 1; end % in case all GCs very small
            mvgc.parameters.mre = maxabs(mvgc.stats.F{tt,tt_1}-mvgc.Fint_3_6Hz{tt,tt_1})/mvgc.parameters.amax;
            if mvgc.parameters.mre < 1e-5
                fprintf('OK (maximum relative error ~ %.0e)\n',mvgc.parameters.mre);
            else
                fprintf(2,'WARNING: high maximum relative error ~ %.0e\n',mvgc.parameters.mre);
            end

        %end
    end
end


% 6 - 8 Hertz
steps                          = diff(mvgc.parameters.freqs); % according to the fft time window
mvgc.parameters.frex_6_8Hz     = 6:steps(1):8;
mvgc.parameters.frex_idx_6_8Hz = dsearchn(mvgc.parameters.freqs,mvgc.parameters.frex_6_8Hz');

for tt = 1:size(mvgc.data,1)
    for tt_1 = 1:size(mvgc.data,2)

        if isempty(mvgc.data{tt,tt_1})
            continue
        end

        %for ii = 1:size(mvgc.data{tt,tt_1},3)

            fprintf('\nfrequency-domain GC integration check... ');
            mvgc.Fint_6_8Hz{tt,tt_1} = smvgc_to_mvgc(mvgc.F_spect{tt,tt_1}(:,:,mvgc.parameters.frex_idx_6_8Hz)); % integrate spectral MVGCs


            mvgc.parameters.amax = maxabs(mvgc.stats.F{tt,tt_1}+mvgc.Fint_6_8Hz{tt,tt_1})/2;

            if mvgc.parameters.amax < 1e-5; mvgc.parameters.amax = 1; end % in case all GCs very small
            mvgc.parameters.mre = maxabs(mvgc.stats.F{tt,tt_1}-mvgc.Fint_6_8Hz{tt,tt_1})/mvgc.parameters.amax;
            if mvgc.parameters.mre < 1e-5
                fprintf('OK (maximum relative error ~ %.0e)\n',mvgc.parameters.mre);
            else
                fprintf(2,'WARNING: high maximum relative error ~ %.0e\n',mvgc.parameters.mre);
            end

        %end
    end
end

% 8 - 10 Hertz
steps                          = diff(mvgc.parameters.freqs); % according to the fft time window
mvgc.parameters.frex_8_10Hz     = 8:steps(1):10;
mvgc.parameters.frex_idx_8_10Hz = dsearchn(mvgc.parameters.freqs,mvgc.parameters.frex_8_10Hz');

for tt = 1:size(mvgc.data,1)
    for tt_1 = 1:size(mvgc.data,2)

        if isempty(mvgc.data{tt,tt_1})
            continue
        end

        %for ii = 1:size(mvgc.data{tt,1},3)


            fprintf('\nfrequency-domain GC integration check... ');
            mvgc.Fint_8_10Hz{tt,tt_1} = smvgc_to_mvgc(mvgc.F_spect{tt,tt_1}(:,:,mvgc.parameters.frex_idx_8_10Hz)); % integrate spectral MVGCs


            mvgc.parameters.amax = maxabs(mvgc.stats.F{tt,tt_1}+mvgc.Fint_8_10Hz{tt,tt_1})/2;

            if mvgc.parameters.amax < 1e-5; mvgc.parameters.amax = 1; end % in case all GCs very small
            mvgc.parameters.mre = maxabs(mvgc.stats.F{tt,tt_1}-mvgc.Fint_8_10Hz{tt,tt_1})/mvgc.parameters.amax;
            if mvgc.parameters.mre < 1e-5
                fprintf('OK (maximum relative error ~ %.0e)\n',mvgc.parameters.mre);
            else
                fprintf(2,'WARNING: high maximum relative error ~ %.0e\n',mvgc.parameters.mre);
            end

        %end
    end
end

clear('tt','ii','steps')



%% Granger over time - Spectrum and Time vector.
% %  VAR model estimation regression mode (LWR) by Mike X Cohen (mikexcohen@gmail.com)
% 
% % - The code relies on the following functions : - grangerX.m 
% %                                                - BSMART ToolBox (http://www.sahs.uth.tmc.edu/hliang/software)
% 
% % -------------------------
% % This session takes up a lot of time. It has been slightly modified (lines 597 to 602 and 605 and 606) to calculate only baseline the trials of interest 
% % -------------------------
% 
% tic
% 
% fprintf('Granger over time - Spectrum and Time vector. VAR model estimation regression mode (LWR) by Mike X Cohen\n')
% 
% % Choose channels combinations
% Combinations_ = nchoosek(1:size(mvgc.data,1),2); % all possible combinations
% 
% % Define band frequency
% % 2 - 12 Hertz
% steps         = diff(mvgc.parameters.freqs); % according to the fft time window
% mvgc.parameters.frex_2_12Hz     = 2:steps(1):12;
% mvgc.parameters.frex_idx_2_12Hz = dsearchn(mvgc.parameters.freqs,mvgc.parameters.frex_2_12Hz');
% 
% win = 1000; % in seconds.
% mvgc.parameters.time_winx  = round(win/(1000/mvgc.parameters.fs)); % in samples
% 
% 
% % Baseline - Comment out this for-loop and make the changes in the following for-loop as suggested above to calculate everything simultaneously.
% for cc = 1:size(Combinations_,1)
% 
%     [mvgc.Time_x2y_Spec{1,1}(:,:,cc),mvgc.Time_y2x_Spec{1,1}(:,:,cc),mvgc.Time_x2y_FInt{1,1}(:,:,cc),mvgc.Time_y2x_FInt{1,1}(:,:,cc)] = ...
%         grangerX(mvgc.data{1,1}(Combinations_(cc,:),:),mvgc.parameters.fs,mvgc.parameters.freqs(mvgc.parameters.frex_idx_2_12Hz),mvgc.parameters.time_winx,mvgc.parameters.morder);
% 
% end
% 
% % CS-Trials or/and ITI trials
% for tt = 2% 1:size(mvgc.data,1) -> 2 for CS-trials
%     for ii = 1:length(CSIT{ms})% size(mvgc.data{tt,1},3) -> selected trials in Analysis.m
% 
%         for cc = 1:size(Combinations_,1)
% 
%             [mvgc.Time_x2y_Spec{tt,ii}(:,:,cc),mvgc.Time_y2x_Spec{tt,ii}(:,:,cc),mvgc.Time_x2y_FInt{tt,ii}(:,:,cc),mvgc.Time_y2x_FInt{tt,ii}(:,:,cc)] = ... 
%                 grangerX(mvgc.data{tt,1}(Combinations_(cc,:),:,CSIT{ms}(ii)),mvgc.parameters.fs,mvgc.parameters.freqs(mvgc.parameters.frex_idx_2_12Hz),mvgc.parameters.time_winx,mvgc.parameters.morder);
% 
%         end
%    end
% end
% 
% clear('tt','ii','cc','steps ','win')
% 
% fprintf('Done.')
% 
% toc
% 
% %% PLot Granger Time Spec. 
% 
% Session = 2;
% 
% Combinations_1 = nchoosek(1:size(mvgc.data,1),2); % all possible combinations
% Combinations_2 = flip(Combinations_1,2);
% 
% % Define band frequency
% % 2 - 12 Hertz
% steps         = diff(mvgc.parameters.freqs); % according to the fft time window
% mvgc.parameters.frex_2_12Hz     = 2:steps(1):12;
% mvgc.parameters.frex_idx_2_12Hz = dsearchn(mvgc.parameters.freqs,mvgc.parameters.frex_2_12Hz');
% 
% % Time window
% win = 1000; % in seconds.
% mvgc.parameters.time_winx  = round(win/(1000/mvgc.parameters.fs)); % in samples
% 
% 
% % Time vector_Granger
% time_v_spec = linspace(0,10,size(mvgc.Time_x2y_Spec{Session,1},2));
% 
% % Behavior data
% data_2_plot_Behav =[];
% for jj = 1:length(CSIT{ms})
%     data_2_plot_Behav(:,:,jj) = decimate(data.behavior{3,CSIT{ms}(jj)*2}(1,1:end-1),4);
% end
% 
% % Normalize behavior vector to the Granger time window
% for jj = 1:size(data_2_plot_Behav,3)
%     for ii = 1:length(data_2_plot_Behav) - mvgc.parameters.time_winx
%         data_2_plot_Behav_m(1,ii,jj) = mean(data_2_plot_Behav(1,ii:ii+mvgc.parameters.time_winx,jj),2);
%     end
% end
% 
% % Time vector_Behavior
% time_v_Behav = linspace(0,10,size(data_2_plot_Behav_m,2));
% 
% 
% %sgtitle({'Amplitude Spectrum via short-window FFT';['Hamming window = ' num2str(short_fft.timewin./1000) 's' ' - ' 'overlap = ' num2str(short_fft.overlap) '%']})
% 
% 
% for jj = 1:length(CSIT{ms})
%     
%     figure 
%     set(gcf, 'Position', get(0, 'Screensize'));
%     set(gcf,'color','w');
% 
%     for cc = 1:size(Combinations_1,1)
%         
%         subplot(3,size(Combinations_1,1),cc)
%         contourf(time_v_spec, mvgc.parameters.freqs(mvgc.parameters.frex_idx_2_12Hz),smoothdata(mvgc.Time_x2y_Spec{Session,jj}(:,:,cc),2,'gaussian',35),150,'linecolor','none');
%         xlabel('Time (s)','FontSize',14), ylabel('Frequency (Hz)','FontSize',14)
%         clim([0 3])
% 
% 
%         if cc == 1
%             title('mPFC PL --> mPFC IL')
%         elseif cc == 2
%             title('mPFC PL --> HPC')
%         else
%             title('mPFC IL --> HPC')
%             c = colorbar;
%             c.Label.String = 'Granger (A.U.)';
%         end
% 
% 
%         subplot(3,size(Combinations_1,1),cc+size(Combinations_1,1))
%         contourf(time_v_spec, mvgc.parameters.freqs(mvgc.parameters.frex_idx_2_12Hz),smoothdata(mvgc.Time_y2x_Spec{Session,jj}(:,:,cc),2,'gaussian',35),150,'linecolor','none');
%         xlabel('Time (s)','FontSize',14), ylabel('Frequency (Hz)','FontSize',14)
%         clim([0 3])
% 
% 
%         if cc == 1
%             title('mPFC IL --> mPFC PL')
%         elseif cc == 2
%             title('HPC --> mPFC PL')
%         else
%             title('HPC --> mPFC IL')
%             c = colorbar;
%             c.Label.String = 'Granger (A.U.)';
%         end
% 
% 
%         subplot(3,size(Combinations_1,1),cc+size(Combinations_1,1)+3)
%         plot(time_v_Behav,data_2_plot_Behav_m(:,:,jj),'linew',2,'color',[0.7 0.7 0.7])
%         set(gca,'ytick',[])
%         set(gca,'ycolor',[1 1 1])
%         box off
%         hold on
%         ylim([0 100])
%         
%         bb = data_2_plot_Behav_m(:,:,jj);
%         ff = (data_2_plot_Behav_m(:,:,jj) < 5);
% 
%         plot(time_v_Behav(ff), bb(ff) ,'k.','linew', 2, 'color', [.6 0 0])
%         xlabel('Time (s)','FontSize',14)
%         legend('movement','freezing')
%     
%     end
% ,
% % Save
% newStr1 = id(1:end-20);
% Path    = files.FilesLoaded{1,1}(ms).folder;
% name_1  = strcat(Path,'/',newStr1,'_Granger_over_time_CS_Trial_',num2str(CSIT{ms}(jj)));
% 
% saveas(gcf,name_1,'png') % save figure
% 
% close all
% 
% end
% 
% clear('ii','ch','steps','win','newStr1','name_1','bb','ff','c','cc','jj','timev','time_v_Behav','data_2_plot_Behav','data_2_plot_Behav_m','time_v_spec','win','steps','Session')
%% Save data

% Settings
%ms = 1;
%Path    = files.FilesLoaded{1,1}(ms).folder;
Path = '/Users/flavio/Desktop';

newStr1 = id(1:end-8);
name_1 = strcat(Path,'/',newStr1,'_MVGC_Granger_extinction');

% Save data
save(name_1,'mvgc','-v7.3')

clear('name','newStr1','path') 


%% last update 13/10/2024
%  listening:
