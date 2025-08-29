movement_exposureB_raw(ms,:)   = data.behavior{1, 1}(1,1:180000);
movement_extinction_raw(ms,:)  = data.behavior{3, 1}(1,1:180000);
movement_retrieval_raw(ms,:)   = data.behavior{3, 1}(1,1:180000);


%%

movement_exposureB = [];
movement_extinction = [];
movement_retrieval = [];

for ii = 1:size(movement_extinction_raw,1)
    movement_exposureB(ii,:)  = decimate(movement_exposureB_raw(ii,:),100);
    movement_extinction(ii,:) = decimate(movement_extinction_raw(ii,:),100);
    movement_retrieval(ii,:)  = decimate(movement_retrieval_raw(ii,:),100);
end

%%
movement_exposureB_filt  = [];
movement_extinction_filt = [];
movement_retrieval_filt  = [];

params.bandstop = 0;
for ii = 1:size(movement_extinction_raw,1)
    movement_exposureB_filt(ii,:)  = abs(hilbert(fun_myfilters(movement_exposureB(ii,:),100,[1 5],'eegfilt',params)));
    movement_extinction_filt(ii,:) = abs(hilbert(fun_myfilters(movement_extinction(ii,:),100,[1 5],'eegfilt',params)));
    movement_retrieval_filt(ii,:)  = abs(hilbert(fun_myfilters(movement_retrieval(ii,:),100,[1 5],'eegfilt',params)));
end

for ii = 1:size(movement_extinction_raw,1)
    movement_exposureB_filt(ii,:)  = normalize(movement_exposureB_filt(ii,:)./sum(movement_exposureB_filt(ii,:)),'range');
    movement_extinction_filt(ii,:) = normalize(movement_extinction_filt(ii,:)./sum(movement_extinction_filt(ii,:)),'range');
    movement_retrieval_filt(ii,:)  = normalize(movement_retrieval_filt(ii,:)./sum(movement_retrieval_filt(ii,:)),'range');
end

%     data.lfp{5,2}(jj,:)  = fun_myfilters(data.lfp{5,1}(jj,:),parameters.decimated_srate,parameters.filter.longcutoff_1,'eegfilt',params);


time_v = linspace(1,170,size(movement_extinction,2));


%% Cumulative time series

cumulative_sum_exposureB = nan(size(movement_exposureB_filt));
cumulative_sum_exposureC = nan(size(movement_exposureC_filt));
cumulative_sum_extinction = nan(size(movement_extinction_filt));
cumulative_sum_retrieval = nan(size(movement_retrieval_filt));


for ii = 1:size(movement_extinction,1)
    cumulative_sum_exposureB(ii,:)        = cumsum(movement_exposureB(ii,:));
    cumulative_sum_exposureB_total(ii,1)  = cumulative_sum_exposureB(ii,end);
    cumulative_sum_exposureB_area(ii,:)   = trapz(movement_exposureB(ii,:),2);

    cumulative_sum_extinction(ii,:)        = cumsum(movement_extinction(ii,:));
    cumulative_sum_extinction_total(ii,1)  = cumulative_sum_extinction(ii,end);
    cumulative_sum_extinction_area(ii,:)   = trapz(movement_extinction(ii,:),2);

    cumulative_sum_retrieval(ii,:)        = cumsum(movement_retrieval(ii,:));
    cumulative_sum_retrieval_total(ii,1)  = cumulative_sum_retrieval(ii,end);
    cumulative_sum_retrieval_area(ii,:)   = trapz(movement_retrieval(ii,:),2);


end

%%
cumulative_sum_exposureB_stat{1,1} = mean(cumulative_sum_exposureB,1,'omitnan');
cumulative_sum_exposureB_stat{1,2} = prctile(cumulative_sum_exposureB,25,1);
cumulative_sum_exposureB_stat{1,3} = prctile(cumulative_sum_exposureB,75,1);

cumulative_sum_extinction_stat{1,1} = mean(cumulative_sum_extinction,1,'omitnan');
cumulative_sum_extinction_stat{1,2} = prctile(cumulative_sum_extinction,25,1);
cumulative_sum_extinction_stat{1,3} = prctile(cumulative_sum_extinction,75,1);

cumulative_sum_retrieval_stat{1,1} = mean(cumulative_sum_retrieval,1,'omitnan');
cumulative_sum_retrieval_stat{1,2} = prctile(cumulative_sum_retrieval,25,1);
cumulative_sum_retrieval_stat{1,3} = prctile(cumulative_sum_retrieval,75,1);

%%
figure;
plot(time_v, cumulative_sum_exposureB_stat{1,1}, 'LineWidth', 2,'Color',[.4 .4 .4]);
hold on;
plot(time_v, cumulative_sum_extinction_stat{1,1}, 'LineWidth', 2);
plot(time_v, cumulative_sum_retrieval_stat{1,1}, 'LineWidth', 2);

xlabel('amplitude');
ylabel('CDF');
legend('exposure', 'extinction','retrieval');
title('Comparação das funções de distribuição acumulada (CDF)');
grid on;

%%

movement_exposureB_filt_conc   = movement_exposureB_(:);
movement_extinction_filt_conc  = movement_extinction_filt(:);
movement_retrieval_filt_conc   = movement_retrieval_filt(:);

[f1, x1] = ecdf(movement_exposureB_filt_conc);
[f2, x2] = ecdf(movement_extinction_filt_conc);
[f3, x3] = ecdf(movement_retrieval_filt_conc);

% CDFs
figure;
plot(x1, f1, 'LineWidth', 2,'Color',[.4 .4 .4]);
hold on;
plot(x2, f2, 'LineWidth', 2);
plot(x3, f3, 'LineWidth', 2);

xlabel('amplitude');
ylabel('CDF');
legend('exposure', 'extinction','retrieval');
title('Comparação das funções de distribuição acumulada (CDF)');
grid on;


%% Teste K-S 
[h, p, ksstat] = kstest2(movement_retrieval_filt_conc, movement_extinction_filt_conc);

fprintf('h: %.4f\n', h);
fprintf('KS: %.4f\n', ksstat);
fprintf('Valor-p: %.4e\n', p); % Formato para visualizar valores muito pequenos

