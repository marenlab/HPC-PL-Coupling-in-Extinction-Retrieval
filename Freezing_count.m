%% loop over all animals

% Set freezing bouts
% Context exposure -> 1 and 3
% Extinction - pre tone -> 3 and 5
% Extinction - fear retrieval -> 3 and 5
% Extinction retrieval - pre tone -> 1
% Extinction - fear retrieval -> 3 and 5

idx = [3 5];

% Set time to %
% - 180 s for pre tone
% - (10 s CS + 30s ITI) * number of trials 
time = 180;

% number of events
time_freezing_{ms,1} = cellfun(@(x) diff(x,[],2),data.events_behavior(idx,CSIT{ms}),'UniformOutput',0);
time_freezing_{ms,2} = vertcat(time_freezing_{ms,1}{:})./1000;
time_freezing_{ms,3} = sum(time_freezing_{ms,2});
time_freezing_{ms,4} = (time_freezing_{ms,3}./time).*100;

%% Concatenated all animals

time_freezing_exposure{1,5}           = vertcat(time_freezing_exposure{1:3,2});
time_freezing_extinction_pretone{1,5} = vertcat(time_freezing_extinction_pretone{1:3,2});
time_freezing_extinction_CS_ITI{1,5}  = vertcat(time_freezing_extinction_CS_ITI{1:3,2});
time_freezing_retrieval_pre_tone{1,5} = vertcat(time_freezing_retrieval_pre_tone{1:3,2});
time_freezing_retrieval_CS_ITI{1,5}   = vertcat(time_freezing_retrieval_CS_ITI{1:3,2});

%% averaging

time_freezing_exposure{2,5}           = vertcat(time_freezing_exposure{1:3,2});
time_freezing_extinction_pretone{2,5} = vertcat(time_freezing_extinction_pretone{1:3,2});
time_freezing_extinction_CS_ITI{2,5}  = vertcat(time_freezing_extinction_CS_ITI{1:3,2});
time_freezing_retrieval_pre_tone{2,5} = vertcat(time_freezing_retrieval_pre_tone{1:3,2});
time_freezing_retrieval_CS_ITI{2,5}   = vertcat(time_freezing_retrieval_CS_ITI{1:3,2});

%%

num_bins = 20;
edges = linspace(0, 20, num_bins + 1);
lim_ = [0 0.1];

% Contar ocorrências em cada bin
counts1 = histcounts(time_freezing_exposure{1,5}, edges);
counts2 = histcounts(time_freezing_extinction_pretone{1,5},edges);
counts3 = histcounts(time_freezing_retrieval_pre_tone{1,5},edges);
counts4 = histcounts(time_freezing_extinction_CS_ITI{1,5},edges);
counts5 = histcounts(time_freezing_retrieval_CS_ITI{1,5},edges);

% Normalizar pela duração total da gravação para obter probabilidade por segundo
prob1 = counts1 / 180;  % Divide pelo tempo total do vetor 1
prob2 = counts2 / 180;  % Divide pelo tempo total do vetor 2
prob3 = counts3 / 180;  % Divide pelo tempo total do vetor 3
prob4 = counts4 / 200;  % Divide pelo tempo total do vetor 4
prob5 = counts5 / 200;  % Divide pelo tempo total do vetor 5



% Obter centros dos bins para alinhamento no plot
bin_centers = edges(1:end-1) + diff(edges)/2;

figure;

subplot(251)
bar(bin_centers, prob1, 'FaceAlpha', 0.5, 'FaceColor', 'b', 'DisplayName', 'exposure');
xlabel('Time (s)');
ylabel('Probability');
ylim(lim_)
title('Context Exposure')

subplot(252)
bar(bin_centers, prob2, 'FaceAlpha', 0.5, 'FaceColor', 'r', 'DisplayName', 'extinction pretone');
xlabel('Time (s)');
ylim(lim_)
title('Extinction pre Tone')

subplot(253)
bar(bin_centers, prob3, 'FaceAlpha', 0.5, 'FaceColor', 'g', 'DisplayName', 'retrieval pre_tone');
xlabel('Time (s)');
ylim(lim_)
title('Retrieval pre Tone')

subplot(254)
bar(bin_centers, prob4, 'FaceColor', 'r', 'DisplayName', 'extinction CS-ITI');
xlabel('Time (s)');
ylim(lim_)
title('Fear retrieval')

subplot(255)
bar(bin_centers, prob5,'FaceColor', 'g', 'DisplayName', 'retrieval CS-ITI');
xlabel('Time (s)');
ylim(lim_)
title('Extinction retrieval')


% ECDF
bandw = 0.3;

[f1, xi1] = ksdensity(counts1, 'Bandwidth', bandw);
cdf_values1 = cumsum(f1) * (xi1(2) - xi1(1));

[f2, xi2] = ksdensity(counts2, 'Bandwidth', bandw);
cdf_values2 = cumsum(f2) * (xi2(2) - xi2(1));

[f3, xi3] = ksdensity(counts3, 'Bandwidth', bandw);
cdf_values3 = cumsum(f3) * (xi3(2) - xi3(1));

[f4, xi4] = ksdensity(counts4, 'Bandwidth', bandw);
cdf_values4 = cumsum(f4) * (xi4(2) - xi4(1));

[f5, xi5] = ksdensity(counts5, 'Bandwidth', bandw);
cdf_values5 = cumsum(f5) * (xi5(2) - xi1(1));


subplot(2,5,[6 7])
ecdf(counts1)
hold on
plot(xi1, cdf_values1, 'LineWidth', 2);

ecdf(counts2)
plot(xi2, cdf_values2, 'LineWidth', 2);

ecdf(counts3)
plot(xi3, cdf_values3, 'LineWidth', 2);

ecdf(counts4)
plot(xi4, cdf_values4, 'LineWidth', 2);

ecdf(counts5)
plot(xi5, cdf_values5, 'LineWidth', 2);

xlabel('Time');
ylabel('ECDF');

ylim([0 1])
xlim([-1 20])

legend(' ','exposure',' ','extinction-pretone',' ','retrieval-pretone',' ','fear',' ','retrieval')


% KSTest2
ks_matriz_p = nan(5,5);
ks_matriz_h = nan(5,5);


[ks_matriz_h(2,1),ks_matriz_p(2,1)] = kstest2(time_freezing_exposure{1,5},time_freezing_extinction_pretone{1,5});
[ks_matriz_h(3,1),ks_matriz_p(3,1)] = kstest2(time_freezing_exposure{1,5},time_freezing_retrieval_pre_tone{1,5});
[ks_matriz_h(4,1),ks_matriz_p(4,1)] = kstest2(time_freezing_exposure{1,5},time_freezing_extinction_CS_ITI{1,5});
[ks_matriz_h(5,1),ks_matriz_p(5,1)] = kstest2(time_freezing_exposure{1,5},time_freezing_retrieval_CS_ITI{1,5});

[ks_matriz_h(3,2),ks_matriz_p(3,2)] = kstest2(time_freezing_extinction_pretone{1,5},time_freezing_retrieval_pre_tone{1,5});
[ks_matriz_h(4,2),ks_matriz_p(4,2)] = kstest2(time_freezing_extinction_pretone{1,5},time_freezing_extinction_CS_ITI{1,5});
[ks_matriz_h(5,2),ks_matriz_p(5,2)] = kstest2(time_freezing_extinction_pretone{1,5},time_freezing_retrieval_CS_ITI{1,5});

[ks_matriz_h(4,3),ks_matriz_p(4,3)] = kstest2(time_freezing_retrieval_pre_tone{1,5},time_freezing_extinction_CS_ITI{1,5});
[ks_matriz_h(5,3),ks_matriz_p(5,3)] = kstest2(time_freezing_retrieval_pre_tone{1,5},time_freezing_retrieval_CS_ITI{1,5});

[ks_matriz_h(5,4),ks_matriz_p(5,4)] = kstest2(time_freezing_extinction_CS_ITI{1,5},time_freezing_retrieval_CS_ITI{1,5});


py_path = "/opt/anaconda3/bin/python";


subplot(2,5,[8 9])
p = imagesc(ks_matriz_p);
set(p,'AlphaData',~isnan(ks_matriz_p))
c = colorbar;
c.Label.String = 'P values';
Py_map = getPyPlot_cMap('RdBu_r', [], [], py_path);
colormap(Py_map)
clim([0 .1])

hold on
for ii = 1:size(ks_matriz_p, 1)
    for jj = 1:size(ks_matriz_p, 2)
        text(jj, ii, sprintf('%.4f', ks_matriz_p(ii, jj)), ...
            'FontSize', 10, 'Color', 'w', 'HorizontalAlignment', 'center');
    end
end

%% Concatenated all animals

% time_freezing_exposure{1,5}           = vertcat(time_freezing_exposure{1:3,2});
% time_freezing_extinction_pretone{1,5} = vertcat(time_freezing_extinction_pretone{1:3,2});
% time_freezing_extinction_CS_ITI{1,5}  = vertcat(time_freezing_extinction_CS_ITI{1:3,2});

time_freezing_retrieval_pre_tone_male = [];
time_freezing_retrieval_CS_ITI_male = [];
time_freezing_retrieval_pre_tone_female = [];
time_freezing_retrieval_CS_ITI_female = [];

time_freezing_retrieval_exposure_male{1,5} = vertcat(time_freezing_exposure{1:3,2});
time_freezing_retrieval_exposure_female{1,5}   = vertcat(time_freezing_exposure{4:6,2});


time_freezing_extinction_pre_tone_male{1,5} = vertcat(time_freezing_extinction_pretone{1:3,2});
time_freezing_extinction_CS_ITI_male{1,5}   = vertcat(time_freezing_extinction_CS_ITI{1:3,2});

time_freezing_extinction_pre_tone_female{1,5} = vertcat(time_freezing_extinction_pretone{4:6,2});
time_freezing_extinction_CS_ITI_female{1,5}   = vertcat(time_freezing_extinction_CS_ITI{4:6,2});


time_freezing_retrieval_pre_tone_male{1,5} = vertcat(time_freezing_retrieval_pre_tone{1:3,2});
time_freezing_retrieval_CS_ITI_male{1,5}   = vertcat(time_freezing_retrieval_CS_ITI{1:3,2});

time_freezing_retrieval_pre_tone_female{1,5} = vertcat(time_freezing_retrieval_pre_tone{4:6,2});
time_freezing_retrieval_CS_ITI_female{1,5}   = vertcat(time_freezing_retrieval_CS_ITI{4:6,2});

%% Plot

num_bins = 10;
edges = linspace(0, 20, num_bins + 1);
lim_ = [0 0.1];

% Contar ocorrências em cada bin
counts1 = histcounts(time_freezing_retrieval_pre_tone_male{1,5},edges);
counts2 = histcounts(time_freezing_retrieval_pre_tone_female{1,5},edges);
counts3 = histcounts(time_freezing_retrieval_CS_ITI_male{1,5},edges);
counts4 = histcounts(time_freezing_retrieval_CS_ITI_female{1,5},edges);

% Normalizar pela duração total da gravação para obter probabilidade por segundo
prob1 = counts1 / 180;  % Divide pelo tempo total do vetor 1
prob2 = counts2 / 180;  % Divide pelo tempo total do vetor 2
prob3 = counts3 / 200;  % Divide pelo tempo total do vetor 3
prob4 = counts4 / 200;  % Divide pelo tempo total do vetor 4

% Obter centros dos bins para alinhamento no plot
bin_centers = edges(1:end-1) + diff(edges)/2;

figure;

subplot(251)
bar(bin_centers, prob1, 'FaceAlpha', 0.5, 'FaceColor', 'b', 'DisplayName', 'preCS Male');
xlabel('Time (s)');
ylabel('Probability');
ylim(lim_)
title('preCS Male')

subplot(252)
bar(bin_centers, prob2, 'FaceAlpha', 0.5, 'FaceColor', 'r', 'DisplayName', 'preCS Famale');
xlabel('Time (s)');
ylim(lim_)
title('preCS Famale')

subplot(253)
bar(bin_centers, prob3, 'FaceAlpha', 0.5, 'FaceColor', 'g', 'DisplayName', 'CS-ITI Male');
xlabel('Time (s)');
ylim(lim_)
title('CS-ITI Male')

subplot(254)
bar(bin_centers, prob4, 'FaceColor', 'r', 'DisplayName', 'CS-ITI Female');
xlabel('Time (s)');
ylim(lim_)
title('CS-ITI Female')


% ECDF
bandw = 0.3;

[f1, xi1] = ksdensity(counts1, 'Bandwidth', bandw);
cdf_values1 = cumsum(f1) * (xi1(2) - xi1(1));

[f2, xi2] = ksdensity(counts2, 'Bandwidth', bandw);
cdf_values2 = cumsum(f2) * (xi2(2) - xi2(1));

[f3, xi3] = ksdensity(counts3, 'Bandwidth', bandw);
cdf_values3 = cumsum(f3) * (xi3(2) - xi3(1));

[f4, xi4] = ksdensity(counts4, 'Bandwidth', bandw);
cdf_values4 = cumsum(f4) * (xi4(2) - xi4(1));



subplot(2,5,[6 7])
ecdf(counts1)
hold on
plot(xi1, cdf_values1, 'LineWidth', 2);

ecdf(counts2)
plot(xi2, cdf_values2, 'LineWidth', 2);

ecdf(counts3)
plot(xi3, cdf_values3, 'LineWidth', 2);

ecdf(counts4)
plot(xi4, cdf_values4, 'LineWidth', 2);

% ecdf(counts5)
% plot(xi5, cdf_values5, 'LineWidth', 2);

xlabel('Time');
ylabel('ECDF');

ylim([0 1])
xlim([-1 20])

legend(' ','preCS Male',' ','preCS Female',' ','Retrieval Male',' ','Retrieval Female')


% KSTest2
ks_matriz_p = nan(4,4);
ks_matriz_h = nan(4,4);

[ks_matriz_h(2,1),ks_matriz_p(2,1)] = kstest2(time_freezing_retrieval_pre_tone_female{1,5},time_freezing_retrieval_pre_tone_male{1,5});
[ks_matriz_h(3,1),ks_matriz_p(3,1)] = kstest2(time_freezing_retrieval_pre_tone_male{1,5},time_freezing_retrieval_CS_ITI_male{1,5});
[ks_matriz_h(4,1),ks_matriz_p(4,1)] = kstest2(time_freezing_retrieval_CS_ITI_female{1,5},time_freezing_retrieval_pre_tone_male{1,5});

[ks_matriz_h(3,2),ks_matriz_p(3,2)] = kstest2(time_freezing_retrieval_pre_tone_female{1,5},time_freezing_retrieval_CS_ITI_male{1,5});
[ks_matriz_h(4,2),ks_matriz_p(4,2)] = kstest2(time_freezing_retrieval_CS_ITI_female{1,5},time_freezing_retrieval_pre_tone_female{1,5});

[ks_matriz_h(4,3),ks_matriz_p(4,3)] = kstest2(time_freezing_retrieval_CS_ITI_female{1,5},time_freezing_retrieval_CS_ITI_male{1,5});


py_path = "/opt/anaconda3/bin/python";


subplot(2,5,[8 9])
p = imagesc(ks_matriz_p);
set(p,'AlphaData',~isnan(ks_matriz_p))
c = colorbar;
c.Label.String = 'P values';
Py_map = getPyPlot_cMap('RdBu_r', [], [], py_path);
colormap(Py_map)
clim([0 .1])

hold on
for ii = 1:size(ks_matriz_p, 1)
    for jj = 1:size(ks_matriz_p, 2)
        text(jj, ii, sprintf('%.4f', ks_matriz_p(ii, jj)), ...
            'FontSize', 10, 'Color', 'w', 'HorizontalAlignment', 'center');
    end
end

