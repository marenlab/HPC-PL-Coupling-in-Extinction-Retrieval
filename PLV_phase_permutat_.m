
%all_data{6,3} = data.lfp{6,1};

plv_real{6,2} = hilb.PLV_mean_session{1,2}(2:end,1);


%%
% Filtering and getting angles from Hilbert transform

filters = [6 8];
fs = 1000;

data_angle = [];

for ss = 1:size(all_data,1)
    for ii = 1:size(all_data,2)

        for ll = 1:size(all_data{ss,ii},1)

            for jj = 1:size(filters,1)

                temp = eegfilt(all_data{ss,ii}(ll,10000:170000),fs,filters(jj,1),filters(jj,2));
                data_angle{ss,ii}(ll,:) = angle(hilbert(temp));

            end
        end

    end

end


clear ('ss','ii','jj','ll','temp','filters')

%%

nsurrog = 1000; % number of rearrangements
trials = 1; % select trials

% plv_null = zeros(2,nsurrog);
% phase_diff_shuffled_all = zeros(2,nsurrog);


for tt = 1:size(trials,1)

    for ss = 1:size(data_angle,1)
        for gg = 1:size(data_angle,2)

            % Baseline
            phase1 = data_angle{ss,gg}(1,:);
            phase2 = data_angle{ss,gg}(2,:);
            phase3 = data_angle{ss,gg}(3,:);


            for ii = 1:nsurrog

                % shift_amount = randi([1, length(phase1)-1]);
                % shuffled_phase1 = circshift(phase1,shift_amount,2);  % permute signal - phase1


                % dHPC - mPFC PL
                shuffled_phase1 = shuffle_esc(unwrap(phase1),fs,200)'; % modified function for Shift array circularly
                phase_diff_shuffled = angle(exp(1i * (unwrap(phase3) - shuffled_phase1)));      
                phase_diff_shuffled_all{ss,gg}(1,ii,tt) = circ_mean(phase_diff_shuffled,[],2);  

                z_shuff = mean(exp(1i * (unwrap(phase3) - shuffled_phase1)));
                plv_null{ss,gg}(1,ii,tt) = abs(z_shuff)';                     

                % dHPC - mPFC IL
                shuffled_phase2 = shuffle_esc(unwrap(phase1),fs,200)'; % modified function for Shift array circularly
                phase_diff_shuffled = angle(exp(1i * (unwrap(phase3) - shuffled_phase2)));      
                phase_diff_shuffled_all{ss,gg}(2,ii,tt) = circ_mean(phase_diff_shuffled,[],2); 

                z_shuff = mean(exp(1i * (unwrap(phase3) - shuffled_phase2))); 
                plv_null{ss,gg}(2,ii,tt) = abs(z_shuff)';                      



            end
        end
    end
end



%% Real values

for ii = 1:size(plv_null,2)
    plv_real{6,ii} = median(plv_real{5,ii},2);
end

for ii = 1:size(plv_null,1)
    for jj = 1:size(plv_null,2)
        plv_real_mean{ii,jj} = median(plv_real{ii,jj},2);
    end
end

for ii = 1:size(plv_null,2)
        % plv_real{5,ii} = cat(3,plv_real{:, ii});
        % plv_real{5,ii} = median(plv_real{5,ii},3);

        % plv_real_mean{5,ii} = cat(3,plv_real_mean{:, ii});
        % plv_real_mean{5,ii} = median(plv_real_mean{5,ii},3);

        % plv_null{5,ii} = cat(4,plv_null{:, ii});
        % plv_null{5,ii} = median(plv_null{5,ii},4);
        % 
        % plv_z{5,ii} = cat(4,plv_z{:, ii});
        % plv_z{5,ii} = median(plv_z{5,ii},4);
        % 
        plv_z_mean{5,ii} = cat(3,plv_z_mean{:, ii});
        plv_z_mean{5,ii} = median(plv_z_mean{5,ii},3);

        % plv_null_mean{5,ii} = cat(4,plv_null_mean{:, ii});
        % plv_null_mean{5,ii} = median(plv_null_mean{5,ii},4);

        % phase_diff_shuffled_all{5,ii} = cat(4,phase_diff_shuffled_all{:,ii});
        % phase_diff_shuffled_all{5,ii} = median(phase_diff_shuffled_all{5,ii},4);

end



        % plv_real_mean{5,ii} = cat(3,plv_real{:, ii});
        % plv_real_mean{5,ii} = median(plv_real{5,ii},3);


clear ('ii','tt','ss','gg','phase1','phase2','phase3','shuffled_phase1','shuffled_phase2','phase_diff_shuffled','z_shuff','nsurrog','trials')


%% Comparing z values

for ii = 1:size(plv_null,1)
    for jj = 1:size(plv_null,2)
        plv_z{ii,jj} = (plv_real{ii,jj} - squeeze(mean(plv_null{ii,jj},2)))./squeeze(std(plv_null{ii,jj},[],2));
    end
end

% for ii = 1:size(plv_real_mean,1)
%     for jj = 1:size(plv_real_mean,2)
%         plv_z_mean{ii,jj} = (plv_real_mean{ii,jj} - squeeze(mean(plv_null_mean{ii,jj},2)))./squeeze(std(plv_null_mean{ii,jj},[],2));
%     end
% end

for ii = 1:size(plv_null,1)
    for jj = 1:size(plv_null,2)
        plv_p{ii,jj} = erfc(plv_z{ii,jj} / sqrt(2)) / 2;
    end
end


%% Plots trials
 
blk_trial = 1;

group = 1; % 3 Extinction Retreival / 2 Extinction / 1 Exposure

if group == 3
    figure
    set(gcf,'color','w');
    set(gcf, 'Position', [1522, -330, 1925, 680]);

    color_ = [1, .4, 0];
elseif group == 2
    color_ = [.8, 0, 0]; 
    hold on
else
    color_ = [.4, .4, .4]; %
    hold on
end


hold on
titles_ = {'Animal 1','Animal 2','Animal 3','Animal 4','Animal 5','Animal 6'};

sgtitle({'\fontsize{18} \bf PLV surrogated ',[]})

% dHPC <-> mPFC-PL
for cc = 1:size(plv_null,1)

    subplot(2,6,cc)
    h = histogram(normalize(plv_null{cc, group}(1,:,blk_trial),2),50,'Normalization','pdf');
    h.FaceColor = color_;
    h.FaceAlpha = .5;
    h.EdgeColor = 'none';
    xlim([-4 30])
    ylim([0 .5])
    xline(1.645)
    ylabel({'\fontsize{14} dHPC <-> mPFC-IL','\fontsize{12} pdf'})
    box off


    if cc~=1
        set(gca, 'YTick', [],'YColor','none')
    end
    
    hold on
    plot([plv_z{cc, group}(1,blk_trial) plv_z{cc, group}(1,blk_trial)],[0 .1],'color',color_, 'linew',4);
    xlabel(['p = ' num2str(plv_p{cc, 1}(1,blk_trial))])
    ylabel({'\fontsize{14} dHPC <-> mPFC-PL','\fontsize{12} pdf'})
    box off

    % if cc~=5
        title({['\fontsize{14}' titles_{cc}],[]})
    % else
    %     title({['\fontsize{14}' 'all animals'],[]})
    % end

end

% dHPC <-> mPFC-IL
for cc = 1:size(plv_null,1)

    subplot(2,6,cc+6)
    h = histogram(normalize(plv_null{cc, group}(2,:,blk_trial),2),50,'Normalization','pdf');
    h.FaceColor = color_;
    h.FaceAlpha = .5;
    h.EdgeColor = 'none';
    xlim([-4 30])
    ylim([0 .5])
    xline(1.645)
    ylabel({'\fontsize{14} dHPC <-> mPFC-IL','\fontsize{12} pdf'})
    box off


    if cc~=1
        set(gca, 'YTick', [],'YColor','none')
    end
    
    hold on
    plot([plv_z{cc, group}(2,blk_trial) plv_z{cc, group}(2,blk_trial)],[0 .1],'color',color_, 'linew',4);
    xlabel(['p = ' num2str(plv_p{cc, 1}(2,blk_trial))])
    ylabel({'\fontsize{14} dHPC <-> mPFC-IL','\fontsize{12} pdf'})
    box off

    if cc~=5
        title({['\fontsize{14}' titles_{cc}],[]})
    else
        title({['\fontsize{14}' 'all animals'],[]})
    end

end


%% Plots averaging

blk_trial = 1; % select trials block

group = 2; % 1 chR2 / 2 mCherry

if group == 1
    figure
    set(gcf,'color','w');
    set(gcf, 'Position', [1522, -330, 1925, 680]);

    color_ = [.2, .6, 1]; % chR2
else
    color_ = [.8, 0, 0]; % mCherry
    hold on

end

hold on
titles_ = {'Animal 1','Animal 2','Animal 3','Animal 4'};

sgtitle({'\fontsize{18} \bf Habituation - 8 Hz',[]})

% dHPC <-> mPFC-IL
for cc = 1:size(plv_null_mean,1)

    subplot(3,5,cc)
    h = histogram(normalize(plv_null_mean{cc, group}(1,:,blk_trial),2),50,'Normalization','pdf');
    h.FaceColor = color_;
    h.FaceAlpha = .5;
    h.EdgeColor = 'none';
    xlim([-4 30])
    ylim([0 .5])
    xline(1.645)
    ylabel({'\fontsize{14} dHPC <-> mPFC-IL','\fontsize{12} pdf'})
    box off


    if cc~=1
        set(gca, 'YTick', [],'YColor','none')
    end
    
    hold on
    plot([plv_z_mean{cc, group}(1,blk_trial) plv_z_mean{cc, group}(1,blk_trial)],[0 .1],'color',color_, 'linew',4);
    xlabel(['p = ' num2str(plv_p{cc, 1}(1,blk_trial))])
    ylabel({'\fontsize{14} mPFC-IL','\fontsize{12} pdf'})
    box off

    if cc~=5
        title({['\fontsize{14}' titles_{cc}],[]})
    else
        title({['\fontsize{14}' 'all animals'],[]})
    end

end

% dHPC <-> mPFC-PL
for cc = 1:size(plv_null_mean,1)

    subplot(3,5,cc+10)
    h = histogram(normalize(plv_null_mean{cc, group}(2,:,blk_trial),2),50,'Normalization','pdf');
    h.FaceColor = color_;
    h.FaceAlpha = .5;
    h.EdgeColor = 'none';
    xlim([-4 30])
    ylim([0 .5])
    xline(1.645)
    ylabel({'\fontsize{14} dHPC <-> mPFC-IL','\fontsize{12} pdf'})
    box off


    if cc~=1
        set(gca, 'YTick', [],'YColor','none')
    end
    
    hold on
    plot([plv_z_mean{cc, group}(2,blk_trial) plv_z_mean{cc, group}(2,blk_trial)],[0 .1],'color',color_, 'linew',4);
    xlabel(['p = ' num2str(plv_p{cc, 1}(2,blk_trial))])
    ylabel({'\fontsize{14} mPFC-IL','\fontsize{12} pdf'})
    box off

    if cc~=5
        title({['\fontsize{14}' titles_{cc}],[]})
    else
        title({['\fontsize{14}' 'all animals'],[]})
    end

end

%% Save

newStr = 'permut_chR2_mcherry_blk_';
path = '/Users/flavio/Desktop';
name = strcat(path,'/',newStr,num2str(blk_trial),'_PLV_extinction');

saveas(gcf,name,'png')
exportgraphics(gcf,strcat(name,'.eps'),'ContentType','vector')

close all

clear('cc','color_1','h','titles_','blk_trial')

%%

