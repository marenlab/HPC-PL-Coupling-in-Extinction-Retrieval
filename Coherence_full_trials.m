
%% Coherence

% 1)
% Magnitute Square Coherence using Welch’s overlapped averaged
% Matlab based welch function --> mscoher

% 2)
% LFP coherogram by multi-taper estimation
% Adapted from --> bz_MTCoherogram.m
% Copyright (C) 2010-2014 by Michaël Zugaro

% The code relies on the following functions:
% --> bz_MTCoherogram -> https://github.com/buzsakilab/buzcode/blob/master/analysis/SpectralAnalyses/bz_MTCoherogram.m
% --> Chronux analysis software -> http://chronux.org
% --> Freely Moving Animal (FMA) Toolbox -> http://fmatoolbox.sourceforge.net/FMAToolbox


% Flavio Mourao.
% email: mourao.fg@tamu.edu
% Maren Lab - Department of Psychological and Brain Sciences
% Texas A&M University

% Started in:  03/2024
% Last update: 03/2024

%%
fprintf('\n Magnitude-squared coherence estimate... \n');

%%  mscoher

% variable: coher_ :

% First row: baseline
% Second row: CS tones
% Third row: IT
% 
% In each cell:
%     - Rows: combinations
%     - Columns: frequencies

coher_.mscohr   = []; % mscoher from matlab built function


% All possible channels combinations
coher_.mscohr.params.combinations  = nchoosek(1:size(data.lfp{5,1},1),2);
% -> row 1 --> mPFC PL <--> mPFC IL
% -> row 2 --> mPFC PL <--> dHPC
% -> row 2 --> mPFC IL <--> dHPC

% Parameters
% Time window
coher_.mscohr.params.timewin    = 2048; % in ms

% Convert time window to points
coher_.mscohr.params.timewinpnts  = hamming(round(coher_.mscohr.params.timewin/(1000/parameters.decimated_srate)));

% Number of overlap samples
coher_.mscohr.params.overlap = 90;
coher_.mscohr.params.noverlap = floor(coher_.mscohr.params.overlap*0.01*(round(coher_.mscohr.params.timewin/(1000/parameters.decimated_srate))));

% nFFT
coher_.mscohr.params.nFFT = 2^nextpow2(round(coher_.mscohr.params.timewin/(1000/parameters.decimated_srate)));
%coher_.mscohr.params.nFFT = 2^15;

% Baseline
not1 = 6;

for ii = 1:length(coher_.mscohr.params.combinations)
    temp_1 = zscore(data.lfp{not1,1}(coher_.mscohr.params.combinations(ii,1),10000:170000),[],2);
    temp_2 = zscore(data.lfp{not1,1}(coher_.mscohr.params.combinations(ii,2),10000:170000),[],2);

    if ii == 1
        [coher_.mscohr.data{1,1}(ii,:),coher_.mscohr.params.freq] = mscohere(temp_1,temp_2,coher_.mscohr.params.timewinpnts,coher_.mscohr.params.overlap,coher_.mscohr.params.nFFT,parameters.decimated_srate,'mimo');
    else
        coher_.mscohr.data{1,1}(ii,:) = mscohere(temp_1,temp_2,coher_.mscohr.params.timewinpnts,coher_.mscohr.params.overlap,coher_.mscohr.params.nFFT,parameters.decimated_srate,'mimo');

    end
end

% CS-tones
% not1 = 7;
% 
% for ii = 1:length(coher_.mscohr.params.combinations)
%     for jj = 1:size(data.lfp{not1,1},3)
% 
%         temp_1 = zscore(data.lfp{not1,1}(coher_.mscohr.params.combinations(ii,1),:,jj),[],2);
%         temp_2 = zscore(data.lfp{not1,1}(coher_.mscohr.params.combinations(ii,2),:,jj),[],2);
% 
%         if ii == 1
%             [coher_.mscohr.data{2,1}(ii,:,jj),coher_.mscohr.params.freq] = mscohere(temp_1,temp_2,coher_.mscohr.params.timewinpnts,coher_.mscohr.params.overlap,coher_.mscohr.params.nFFT,parameters.decimated_srate);
%         else
%             coher_.mscohr.data{2,1}(ii,:,jj) = mscohere(temp_1,temp_2,coher_.mscohr.params.timewinpnts,coher_.mscohr.params.overlap,coher_.mscohr.params.nFFT,parameters.decimated_srate);
%         end
% 
%     end
% end
% 
% % ITI
% not1 = 8;
% 
% for ii = 1:length(coher_.mscohr.params.combinations)
%     for jj = 1:size(data.lfp{not1,1},3)
% 
%         temp_1 = zscore(data.lfp{not1,1}(coher_.mscohr.params.combinations(ii,1),:,jj),[],2);
%         temp_2 = zscore(data.lfp{not1,1}(coher_.mscohr.params.combinations(ii,2),:,jj),[],2);
%         
%         if ii == 1
%             [coher_.mscohr.data{3,1}(ii,:,jj),coher_.mscohr.params.freq] = mscohere(temp_1,temp_2,coher_.mscohr.params.timewinpnts,coher_.mscohr.params.overlap,coher_.mscohr.params.nFFT,parameters.decimated_srate);
%         else
%             coher_.mscohr.data{3,1}(ii,:,jj) = mscohere(temp_1,temp_2,coher_.mscohr.params.timewinpnts,coher_.mscohr.params.overlap,coher_.mscohr.params.nFFT,parameters.decimated_srate);
%         end
% 
%     end
% end

clear ('ii','jj','temp_1','temp_2','not1')

%% Coherogram by Magnitute Square Coherence using Welch.s overlapped averaged. Time blocks

% coher_.coherogram = [];
% 
% % All possible channels combinations
% coher_.coherogram.params.combinations  = nchoosek(1:size(data.lfp{5,1},1),2);
% % -> row 1 --> mPFC PL <--> mPFC IL
% % -> row 2 --> mPFC PL <--> dHPC
% % -> row 2 --> mPFC IL <--> dHPC
% 
% 
% % Welch Parameters
% % Time window
% coher_.coherogram.params.timewin    = 1000; % in ms
% 
% % Convert time window to points
% coher_.coherogram.params.timewinpnts  = hamming(round(coher_.coherogram.params.timewin/(1000/parameters.decimated_srate)));
% 
% % Number of overlap samples
% coher_.coherogram.params.overlap_1 = 95;
% coher_.coherogram.params.noverlap = floor(coher_.coherogram.params.overlap_1*0.01*(round(coher_.coherogram.params.timewin/(1000/parameters.decimated_srate))));
% 
% % nFFT
% coher_.coherogram.params.nFFT = 2^15; %2^nextpow2(round(coher_.mscohr.params.timewin/(1000/parameters.srate)));
% 
% 
% % Time course  Parameters
% % Window length
% coher_.coherogram.params.window = 4000; %s ;
% 
% overlap = .9375; % (%)
% coher_.coherogram.params.overlap_2 = overlap * coher_.coherogram.params.window;
% coher_.coherogram.params.step = coher_.coherogram.params.window - coher_.coherogram.params.overlap_2;
%    
% 
% % Baseline
% not1 = 6;
% 
% Nwindow = (length(data.lfp{not1,1})-coher_.coherogram.params.window)/coher_.coherogram.params.step+1;
% 
% for nwin = 1:Nwindow
%     win = (1:coher_.coherogram.params.window) + (nwin-1)* coher_.coherogram.params.step ;
% 
%     for ii = 1:length(coher_.coherogram.params.combinations)
% 
%         temp_1 = data.lfp{not1,1}(coher_.coherogram.params.combinations(ii,1),win);
%         temp_2 = data.lfp{not1,1}(coher_.coherogram.params.combinations(ii,2),win);
% 
%         [coher_.coherogram.data{1,1}(ii,:,nwin),coher_.coherogram.params.freq] = mscohere(temp_1, temp_2, coher_.coherogram.params.timewinpnts, coher_.coherogram.params.overlap_1, coher_.coherogram.params.nFFT, parameters.decimated_srate);
%         coher_.coherogram.params.timev{1,1}(nwin) = median(win)/parameters.decimated_srate;
% 
%     end
% 
% end
% 
% % CS-tones
% not1 = 7;
% 
% Nwindow = (length(data.lfp{not1,1})-coher_.coherogram.params.window)/coher_.coherogram.params.step+1;
% 
% for nwin = 1:Nwindow
%     win = (1:coher_.coherogram.params.window) + (nwin-1)* coher_.coherogram.params.step ;
% 
%     for ii = 1:length(coher_.coherogram.params.combinations)
%         for jj = 1:size(data.lfp{not1,1},3)
% 
%         temp_1 = data.lfp{not1,1}(coher_.coherogram.params.combinations(ii,1),win,jj);
%         temp_2 = data.lfp{not1,1}(coher_.coherogram.params.combinations(ii,2),win,jj);
% 
%         [coher_.coherogram.data{2,1}(ii,:,nwin,jj),coher_.coherogram.params.freq] = mscohere(temp_1, temp_2, coher_.coherogram.params.timewinpnts, coher_.coherogram.params.overlap_1, coher_.coherogram.params.nFFT, parameters.decimated_srate);
%         coher_.coherogram.params.timev{2,1}(nwin) = median(win)/parameters.decimated_srate;
% 
%         end
%     end
% end
% 
% % ITI
% not1 = 8;
% 
% Nwindow = (length(data.lfp{not1,1})-coher_.coherogram.params.window)/coher_.coherogram.params.step+1;
% 
% for nwin = 1:Nwindow
%     win = (1:coher_.coherogram.params.window) + (nwin-1)* coher_.coherogram.params.step ;
% 
%     for ii = 1:length(coher_.coherogram.params.combinations)
%         for jj = 1:size(data.lfp{not1,1},3)
% 
%         temp_1 = data.lfp{not1,1}(coher_.coherogram.params.combinations(ii,1),win,jj);
%         temp_2 = data.lfp{not1,1}(coher_.coherogram.params.combinations(ii,2),win,jj);
% 
%         [coher_.coherogram.data{3,1}(ii,:,nwin,jj),coher_.coherogram.params.freq] = mscohere(temp_1, temp_2, coher_.coherogram.params.timewinpnts, coher_.coherogram.params.overlap_1, coher_.coherogram.params.nFFT, parameters.decimated_srate);
%         coher_.coherogram.params.timev{3,1}(nwin) = median(win)/parameters.decimated_srate;
% 
%         end
%     end
% end
% 
% clear ('nwin','ii','jj','temp_1','temp_2','not1','nwindow','win','overlap')

%% Coherogram by multi-taper estimation
%
%    [coherogram,phase,t,f] = MTCoherogram(lfp1,lfp2,<options>)
%
%    lfp1,lfp2      wide-band LFPs (one channel each).
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'srate'       sampling rate (in Hz) (default = from timestamps if
%                   available, otherwise 1250Hz)
%     'range'       frequency range (in Hz) (default = all)
%     'window'      duration (in s) of the time window (default = 5)
%     'overlap'     overlap between successive windows (default = window/2)
%     'step'        step between successive windows (default = window/2)
%     'tapers'      relative resolution and order of the tapers [NW K]
%                   (default = [3 5])
%     'pad'         FFT padding (see help for <a href="matlab:help cohgramc">cohgramc</a>) (default = 0)
%     'show'        plot results (default = 'off')
%     'cutoffs'     cutoff values for color plot (default = [0 1])
%    =========================================================================
%
%  NOTES
%
%    The LFP can be provided either as a time stamped matrix (list of time-voltage
%    pairs), or as a voltage vector - in which case the frequency must be specified.
%
%    The time displacement between successive short time coherences can be supplied
%    either as a 'step' (explicit time difference) or as an 'overlap' (between
%    successive time windows).
%
%  OUTPUT
%
%    coherogram     coherogram magnitude
%    phase          coherogram phase
%    t              time bins
%    f              frequency bins
%
%%
fprintf('\n Coherogram by multi-taper estimation... \n');

startup_chronux
startup_FMA

clear("chronux_root",'startup_FMA','FMA_root')

%% Parameters for coherogram following instructions above

% First row: baseline
% In each cell:
%     - Rows: combinations
%     - Columns: frequencies
%     - 3th dimention: time

% Second row: CS tones
% Third row: ITI - Inter trials

% In each cell:
%     - Rows: combinations
%     - Columns: frequencies
%     - 3th dimention: time
%     - 4th dimention: trials


 coher_.coherogram_MT = [];

coher_.coherogram_MT.params.combinations  = nchoosek(1:size(data.lfp{5,1},1),2);
coher_.coherogram_MT.params.srate = parameters.decimated_srate;
coher_.coherogram_MT.params.range = [1 100];
coher_.coherogram_MT.params.window = [18];

overlap = .95; % (%)

coher_.coherogram_MT.params.overlap = overlap * coher_.coherogram_MT.params.window;
coher_.coherogram_MT.params.step = coher_.coherogram_MT.params.window - coher_.coherogram_MT.params.overlap;
coher_.coherogram_MT.params.tapers = [3 5];
coher_.coherogram_MT.params.pad = 0;
coher_.coherogram_MT.params.show = 'on';
coher_.coherogram_MT.params.cutoffs = [0.1 .6];

%%

% params = struct('tapers',[3 5],'pad',[0],'fs',1000,'fpass',[1 100],'err',[0],'trialave',[0]);

%% LFP coherogram by multi-taper estimation

% All possible channels combinations following mscoherence parameters defined above

% baseline

not1 = 6;

for ii = 1:length(coher_.coherogram_MT.params.combinations)
    temp_1 = data.lfp{not1,1}(coher_.coherogram_MT.params.combinations(ii,1),B_clean{ms}(1):B_clean{ms}(2));
    temp_2 = data.lfp{not1,1}(coher_.coherogram_MT.params.combinations(ii,2),B_clean{ms}(1):B_clean{ms}(2));


    [coher_.coherogram_MT.data{1,1}(ii,:,:),coher_.coherogram_MT.phase{1,1}(ii,:,:),coher_.coherogram_MT.params.timev{1,1},coher_.coherogram_MT.params.freq] = ...
        bz_MTCoherogram(temp_1',temp_2','frequency',coher_.coherogram_MT.params.srate,'range',coher_.coherogram_MT.params.range,'window',coher_.coherogram_MT.params.window,...
        'overlap',coher_.coherogram_MT.params.overlap,'step',coher_.coherogram_MT.params.step,'tapers',coher_.coherogram_MT.params.tapers,'pad',coher_.coherogram_MT.params.pad,'show',coher_.coherogram_MT.params.show,'cutoffs',coher_.coherogram_MT.params.cutoffs);
  
%     [C(ii,:,:),phi(ii,:,:),S12(ii,:,:,:),S1(ii,:,:),S2(ii,:,:),t,f] = ...
%         cohgramcpb(temp_1', temp_2',[18000 9000],params);

end

% 
% % CS-tones
% not1 = 7;
% 
% for ii = 1:length(coher_.coherogram_MT.params.combinations)
%     for jj = 1:size(data.lfp{not1,1},3)
% 
%         temp_1 = zscore(data.lfp{not1,1}(coher_.coherogram_MT.params.combinations(ii,1),:,jj),[],2);
%         temp_2 = zscore(data.lfp{not1,1}(coher_.coherogram_MT.params.combinations(ii,2),:,jj),[],2);
% 
%         [coher_.coherogram_MT.data{2,1}(ii,:,:,jj),coher_.coherogram_MT.phase{2,1}(ii,:,:,jj),coher_.coherogram_MT.params.timev{2,1},coher_.coherogram_MT.params.freq] = ...
%             bz_MTCoherogram(temp_1',temp_2','frequency',coher_.coherogram_MT.params.srate,'range',coher_.coherogram_MT.params.range,'window',coher_.coherogram_MT.params.window,...
%             'overlap',coher_.coherogram_MT.params.overlap,'step',coher_.coherogram_MT.params.step,'tapers',coher_.coherogram_MT.params.tapers,'pad',coher_.coherogram_MT.params.pad,'show',coher_.coherogram_MT.params.show,'cutoffs',coher_.coherogram_MT.params.cutoffs);
% 
%     end
% end
% 
% % ITI - ( Inter trials)
% not1 = 8;
% 
% for ii = 1:length(coher_.coherogram_MT.params.combinations)
%     for jj = 1:size(data.lfp{not1,1},3)
% 
%         temp_1 = zscore(data.lfp{not1,1}(coher_.coherogram_MT.params.combinations(ii,1),:,jj),[],2);
%         temp_2 = zscore(data.lfp{not1,1}(coher_.coherogram_MT.params.combinations(ii,2),:,jj),[],2);
% 
%         [coher_.coherogram_MT.data{3,1}(ii,:,:,jj),coher_.coherogram_MT.phase{3,1}(ii,:,:,jj),coher_.coherogram_MT.params.timev{3,1},coher_.coherogram_MT.params.freq] = ...
%             bz_MTCoherogram(temp_1',temp_2','frequency',coher_.coherogram_MT.params.srate,'range',coher_.coherogram_MT.params.range,'window',coher_.coherogram_MT.params.window,...
%             'overlap',coher_.coherogram_MT.params.overlap,'step',coher_.coherogram_MT.params.step,'tapers',coher_.coherogram_MT.params.tapers,'pad',coher_.coherogram_MT.params.pad,'show',coher_.coherogram_MT.params.show,'cutoffs',coher_.coherogram_MT.params.cutoffs);
% 
%     end
% end


clear ('ii','jj','temp_1','temp_2','not1','overlap')

%% Save

%newStr = regexprep(files.id.name,'.mat','_');
%newStr = files.id(ms).name(1:end-8);
newStr = id(1:end-5);

path = '/Users/flavio/Desktop';
%path = files.FilesLoaded{1, 1}.folder;

name = strcat(path,'/',newStr,'_coherence_nfft_32768_fullTrial_2048TimeW');

% save data
save(name,'coher_','-v7.3')

clear('name','newStr','path') 

%% last update 27/03/2024 - 17:56
% listening: Godspeedyou! Black emperor - https://www.youtube.com/watch?v=LViR6liUK2k 
