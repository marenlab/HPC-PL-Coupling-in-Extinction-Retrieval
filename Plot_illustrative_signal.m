


signal_filter = [];

hilb_extinction.parameters.filters_not = [58 62];
hilb_extinction.parameters.filters = [0 100; 30 100];

params.bandstop = 1;

% Baseline - NOT
for ii = 1:size(data.lfp{6,1},1)
    for jj = 1:size(hilb_extinction.parameters.filters_not,1)

        signal_filter{1,1}(ii,:) = fun_myfilters(data.lfp{6,1}(ii,10000:170000),parameters.decimated_srate,hilb_extinction.parameters.filters_not(jj,:),'iir',params); % just filtering

    end

end

params.bandstop = 0;

% Baseline
for ii = 1:size(data.lfp{6,1},1)
    for jj = 1:size(hilb_extinction.parameters.filters,1)

        signal_filter{2,jj}(ii,:) = fun_myfilters(signal_filter{1,1}(ii,:),parameters.decimated_srate,hilb_extinction.parameters.filters(jj,:),'eegfilt',params); % just filtering

    end

end

%%
timev = linspace(0, 160, length(signal_filter{2,1}));

%%

xlim_ = [1 180];

figure

subplot(1,3,1)
plot(timev,signal_filter{2,1}(1,:),'color',[.8 .8 .8])
hold on
plot(timev,signal_filter{2,2}(1,:),'linew',1)
ylim([-500 500])
xlim(xlim_)

subplot(1,3,2)
plot(timev,signal_filter{2,1}(2,:),'color',[.8 .8 .8])
hold on
plot(timev,signal_filter{2,2}(2,:),'linew',1)
ylim([-500 500])
xlim(xlim_)

subplot(1,3,3)
plot(timev,signal_filter{2,1}(3,:),'color',[.8 .8 .8])
hold on
plot(timev,signal_filter{2,2}(3,:),'linew',1)
ylim([-500 500])
xlim(xlim_)