%%

function [data, parameters] = Extracting_LFPs_and_events(id,FilesLoaded,Path)

% Extracting LFPs and Events from Intan/Open Ephys or Plexon *.nex files

% ------------------------------------------------------
% Important note
% Obs.: For now, I am analyzing the lab repository. The LFP records were decimated from 40k to 1k. There were 82 channels in Plexon *.nex file;
%       however, officially, there were recorded 32 channels: 16 in the mPFC and 16 in the HPC, according to the paper: https://doi.org/10.1038/s41467-023-42315-1

% It seems, according to the labels, FP = field potential. Channel 49 to channel 81. Channel 82 is the signal related to vibration for freezing analysis.
% -> Ask Steve or Tucge if any other synchronization signals were saved. For example, the start/end of CS.

% ------------------------------------------------------

% - Extract, organize and save data from Intan/Open Ephys and Plexon files: *.continuous, *.events, *.nex

% - The code relies on the following functions : -> open_files.m
%                                                -> load_open_ephys_data.m                  (https://github.com/open-ephys/analysis-tools)
%                                                -> load_open_ephys_binary.m                (https://github.com/open-ephys/analysis-tools)
%                                                -> npy-matlab                              (https://github.com/open-ephys/analysis-tools)
%                                                -> OmniPlex and MAP Offline SDK Bundle     (https://plexon.com/software-downloads/#software-downloads-SDKs)
%                                                -> FieldTrip                               (https://www.fieldtriptoolbox.org)

%                                                -> xml2struct (https://www.mathworks.com/matlabcentral/fileexchange/28518-xml2struct)

% - Option: down sampling data

% - Outputs from OpenEphys files:

%   "data"
%   -> data.raw{1,1}  -> Cell -> First cell column: original signal.
%                        Each cell: Rows: Channels x Columns: Time
%   -> data.timev_raw -> time vector. Original sample rate

%   -> data.data{1,1} -> Cell -> First cell column: signal decimated.
%                        Each cell: Rows: Channels x Columns: Time
%   -> data.timev     -> time vector. Signal decimated

%   -> events         -> External TTls. Events that are detected in a continuous data stream
%                        _ Supports up to 8 inputs. For digital inputs - labels: 0 - 7.
%                        _ ts -> All timestamps
%                        _ ts_sort .... -> sorted according to the labels
%                                          each cell column corresponds to a recorded event type
%                                          within the cell -> lines -> timestamps(seconds)


%   "parameters"
%   ->  Record informations and parameters
%   ->  Impedance test made by OpenEphys
%       Columns: number of impedance measures


% - Outputs from Plexon files:


% by Flavio Mourao.


%% Define and save number of channels loaded in case of *.continuous (Open Ephys extension files)
if contains(FilesLoaded{1}, 'continuous') == 1
    parameters.nch = sum(contains(FilesLoaded, 'continuous'));
else
    parameters.nch = 1; % Case for Binary files from RHD 2132 or *.nex files; this will change accordingly later
end

%% Choose factor to LFP down sampling
% - Manually -
% parameters.downsampling = 6;
parameters.downsampling = 1;

% - Request from user -
%         prompt        = {'Decimation Factor:'};
%         dlgtitle      = 'Please enter';
%         dims          = [1 30];
%         default_input = {'6'};
%
%         input = inputdlg(prompt,dlgtitle,dims,default_input); %gui
%
%         parameters.downsampling = str2double(input{1,1});
%
%         clear ('prompt','dlgtitle','dims','default_input','input')

%% Loop to extract data


for jj = 1:length(FilesLoaded)
    baseFileName = FilesLoaded{jj};

    if contains(baseFileName,{'settings','messages'})
        continue
    end

    fullFileName = fullfile(Path,baseFileName);

    %Identify the file extension
    [~, ~, fExt] = fileparts(baseFileName);


    switch lower(fExt)

        % Case for load *.continuous (Open Ephys extension files)

        case '.continuous'

            if jj == 1
                fprintf(1, '\nExtracting Signal from *.continuous (Open Ephys extension)\n');
            end

            % Identify the channel number and print the name on the Command Window:
            % channels   1 to 16 and/or 17 to 32

            channel = split(baseFileName,{'100_CH','.continuous'});
            fprintf(1, '\nExtracting Signal from Channel %s\n', channel{2, 1});


            if      jj == 1 && parameters.downsampling == 1

                % Load datafiles (*.continuous), timestamps e record info.
                % Raw data - Rows: Time  x Columns: Channels
                % Remove linear trend (Matlab build function 'detrend'/ Method: 'constant' - subtract the mean from the data)
                [data_temp, data.timev_raw, info] = load_open_ephys_data(fullFileName);
                data(1,1).raw{1,1}  = zeros(parameters.nch, ceil(length(data_temp)));
                data.raw{1,1}       = (detrend(data_temp, 'constant'))';  % Raw data
                parameters.header   = info.header;                        % Data File Header
                parameters.srate    = info.header.sampleRate;

                % Normalizing time vector according to start record
                parameters.const_time = min(data.timev_raw); % exact time when the record was started after play viewing
                data.timev_raw  = data.timev_raw - parameters.const_time;  % Time stamp (sec)

            elseif  jj == 1 && parameters.downsampling > 1

                % Load datafiles (*.continuous), timestamps e record info.
                % Data - Rows: Time  x Columns: Channels
                % Remove linear trend (Matlab build function 'detrend'/ Method: 'constant' - subtract the mean from the data)
                [data_temp, data_timev, info] = load_open_ephys_data(fullFileName);
                data_temp       = (detrend(data_temp, 'constant'))';  % Raw data
                parameters.header   = info.header;

                % Downsampling with Matlab decimate function
                % Data - Rows: Channels x Columns: Time
                data(1,1).data{1,1}  = zeros(parameters.nch, ceil(length(data_temp)/parameters.downsampling));
                data.data{1,1}       = decimate(data_temp,parameters.downsampling);

                % Organize parameters according to the downsampling information
                parameters.srate  = info.header.sampleRate./parameters.downsampling;  % Sampling frequency after downsamplig(Hz)
                parameters.header.downsampling = parameters.downsampling;

                % Normalizing time vector according to start record
                parameters.const_time = min(data_timev); % exact time when the record was started after play viewing
                data.timev      = (data_timev(1:parameters.downsampling:end)) - parameters.const_time;  % Time stamp (sec)

            elseif  jj > 1 && parameters.downsampling == 1

                % Load datafiles (*.continuous).
                % Rows: Time x Columns: Channels
                % Remove linear trend (function 'detrend')
                data.raw{1,1}(jj,:)  = (detrend(load_open_ephys_data(fullFileName),'constant'))'; % Raw data

            elseif  jj > 1 && parameters.downsampling > 1

                % Load datafiles (*.continuous).
                % Rows: Time x Columns: Channels
                % Remove linear trend (function 'detrend')
                data_temp(jj,:)  = (detrend(load_open_ephys_data(fullFileName),'constant'))'; % Raw data

                % Downsampling with Matlab decimate function
                % data - Rows: Channels x Columns: Time
                data.data{1,1}(jj,:) = decimate(data_temp(jj,:),parameters.downsampling); % parameters.downsampling with Matlab decimate function

            end


        % Case for load events

        case '.events'

            % Identify TTL events file and print the name on the Command Window:
            fprintf(1, '\nExtracting Events \n');

            % Load datafiles (*.continuous), timestamps e record info.
            [data.events.labels, data.events.ts, parameters.events.info] = load_open_ephys_data(fullFileName);

            % Sort Events
            % Trigger/events labels according to the digital inputs
            labels = 0:7;

            % data.events.ts_sort -> each cell column corresponds to a recorded event type
            %                        within the cell -> lines -> timestamps(seconds)

            for ii = 1:length(labels)
                data.events.ts_sort{1,ii} = data.events.ts(data.events.labels(:,1) == labels(ii));
            end

            % Normalizing events according to start record

            for ii = 1:length(data.events.ts_sort)
                if isempty(data.events.ts_sort{ii})
                    continue;
                else
                    data.events.ts_sort{ii} = data.events.ts_sort{ii} - parameters.const_time;
                end
            end


        % Case for load Binary files *.dat.

        case '.oebin'

            fprintf(1, '\nExtracting Signal from Binary File - %s\n', id);
            %fprintf(1, '\nExtracting Signal from Binary File \n')

            if  parameters.downsampling == 1

                % Load datafiles (*.continuous), timestamps e record info.
                % Raw data - Rows: Time  x Columns: Channels
                % Remove linear trend (Matlab build function 'detrend'/ Method: 'constant' - subtract the mean from the data)
                data_temp = load_open_ephys_binary(fullFileName,'continuous',1);
                data.raw{1,1}  = zeros(size(data_temp.Data));
                data.raw{1,1} = (detrend(data_temp.Data', 'constant'))';  % Raw data

                parameters.header = data_temp.Header;                     % Data File Header
                parameters.srate  = data_temp.Header.sample_rate;
                parameters.header.sampleRate  = data_temp.Header.sample_rate; % Redundancy to normalize the variable name
                clear('parameters.Header.sample_rate');

                % Normalizing time vector according to start record
                data.timev_raw = double(data_temp.Timestamps)./data_temp.Header.sample_rate;
                parameters.const_time = min(data.timev_raw); % exact time when the record was started after play viewing
                data.timev_raw  = data.timev_raw - parameters.const_time;  % Time stamp (sec)

            elseif  parameters.downsampling > 1

                % Load datafiles (*.continuous), timestamps e record info.
                %Data - Rows: Time  x Columns: Channels
                % Remove linear trend (Matlab build function 'detrend'/ Method: 'constant' - subtract the mean from the data)
                data_temp = load_open_ephys_binary(fullFileName,'continuous',1);


                data_detrended  = zeros(size(data_temp.Data));
                data_detrended = (detrend(data_temp.Data', 'constant'))';     % Raw data

                parameters.header = data_temp.Header;                         % Data File Header
                parameters.header.sampleRate  = data_temp.Header.sample_rate; % Redundancy to normalize the variable name
                clear('parameters.Header.sample_rate');

                % Downsampling with Matlab decimate function
                % data - Rows: Channels x Columns: Time
                data(1,1).data{1,1}  = zeros(size(data_temp.Data,1), ceil(length(data_temp.Data)/parameters.downsampling));

                for kk = 1:size(data.raw{1,1},1)
                    data.data{1,1}(kk,:) = decimate(data_detrended(kk,:),parameters.downsampling);
                end

                % Organize parameters according to the downsampling information
                parameters.srate  = data_temp.Header.sample_rate./parameters.downsampling;  % Sampling frequency after downsamplig(Hz)
                parameters.header.downsampling = parameters.downsampling;

                % Number of channels recorded
                parameters.nch = size(data_temp.Data,1);

                % Normalizing time vector according to start record
                data.timev            = double(data_temp.Timestamps)./data_temp.Header.sample_rate;
                parameters.const_time = min(data.timev); % exact time when the record was started after play viewing
                data.timev            = (data.timev(1:parameters.downsampling:end)) - parameters.const_time;  % Time stamp (sec)

            end

            % load events from Binary File
            fprintf(1, '\nExtracting Events\n');

            % Load datafiles
            data_temp_events = load_open_ephys_binary(fullFileName,'events',1);

            % Rename variable
            data.events.labels = data_temp_events.ChannelIndex;

            % Sort Events
            % Trigger/events labels according to the digital inputs
            labels = 1:8;

            % Convert samples to seconds
            data.events.ts = double(data_temp_events.Timestamps)./data_temp.Header.sample_rate;

            % data.events.ts_sort -> each cell column corresponds to a recorded event type
            %                        within the cell -> lines -> timestamps(seconds)

            for ii = 1:length(labels)
                data.events.ts_sort{1,ii} = data.events.ts(data_temp_events.ChannelIndex(:,1) == labels(ii));
            end

            % Normalizing events according to start record

            for ii = 1:length(data.events.ts_sort)
                if isempty(data.events.ts_sort{ii})
                    continue;
                else
                    data.events.ts_sort{ii} = data.events.ts_sort{ii} - parameters.const_time;
                end
            end



            % Case for load *.nex files from Plexon

        case {'.nex5','.nex'}% or .nex5 files
 
            % ------------------------------------------------------

            % Important Notes:
            % It seems, according to the labels, FP = field potential and WB = wide band. Channel 49 to channel 81. Channel 82 is the signal related to vibration for freezing analysis.
            % -> Ask Steve or Tucge if any other synchronization signals were saved. For example, the start/end of CS.

            % You may need to adjust the scale in Plexon under the Preamp/Total gain to have a gain of 153. 
            % This is because Plexon imports DDT files as 16-bit files with an assumed max resolution of 5 V, which makes the resolution 0.153 uV/bit, but you will probably 
            % want to adjust that resolution to be 1 uV / bit. However, import your data without this change first to see if the scaling is off, then adjust accordingly and check again.
            
            % To import *.nex5 the original code from fieldtrip was adapted. The sample conouter was not right....

            % ------------------------------------------------------


            fprintf(1, '\nExtracting Signal from Plexon File - %s\n', id);

            % For Plexon records, I chose to always perform the conversion without downsampling the data. Decimation will be an option during data processing

            if  parameters.downsampling == 1 || parameters.downsampling > 1

                % Load datafiles (*.nex), timestamps e record info.
                % Raw data - Rows: Channels x Columns: Time
                % Remove linear trend (Matlab build function 'detrend'/ Method: 'constant' - subtract the mean from the data)

                parameters.header = ft_read_header(fullFileName); % Data File Header

                data_temp         = ft_read_data(fullFileName);   % Raw data
                data_temp         = data_temp(contains(parameters.header.label, ["FP","WB"]) == 1,:); % remove empty channels according to the labels

                % Remove NaN`s
                %data_temp(sum(isnan(data_temp),1) == size(data_temp,2),:) = []; % remove rows with all NaN's in a row
                data_temp(:,sum(isnan(data_temp),1) == size(data_temp,1)) = []; % remove columns with all NaN's in a column

                data.lfp{1,1}  = zeros(size(data_temp));
                data.lfp{1,1}  = (detrend(data_temp', 'constant','omitnan'))';  % Raw data
                
                % Scaling in uV
                data.lfp{1,1} = data.lfp{1,1} .* 1/0.153;

                % Behavior
                data.behavior{1,1} = abs(data.lfp{1,1}(end,:));                            % Isolating the behavior variable. Configured on the last line of the record. Usually last channel from analog.
                data.behavior{1,1} = (data.behavior{1,1}./max(data.behavior{1,1})).*100;   % normalize values to %

                data.lfp{1,1}(end,:) = [];                           % Deleting the channel related to the variable's behavior along with the field recording.

                parameters.original_srate  = parameters.header.Fs;   % Redundancy to normalize the variable name to further analysis
                parameters.nch             = size(data.lfp{1,1},1);  % number of channels

                % Normalizing time vector according to start record
                data.timev = linspace(0,length(data.lfp{1,1})/parameters.original_srate,length(data.lfp{1,1}));


                % load events from Plexon File
                fprintf(1, '\nExtracting Events from Plexon File\n');

                % According to the events recorded by Nature Communication:

                % The variable 'events' is a struct.  The 'sample' field is at the decimated sampling rate.
                %                                     The 'timestamp' field is at the original sampling rate

                %   - The events begin approximately in 180 seconds (baseline period).
                %   - 'EVT05' is the label where CS is applied but depends on the plexons settings from old MedPC chamber
                %   - 'EVT08' is the label where CS is applied but depends on the plexons settings from new MedPC chamber video freeze system
                %   - 'Segment' refers to the interval between each CS (ITI).
                %   -  Each timestamp separated by 1 sec (Totty and Steve was using tone pips).

                % CS and the ITI varied. CS varied between 9 and 10 secs, and the ITI between 20 and 30 secs.
                % According to the paper, CS = 10 sec and ITI = 30 sec."


                % Records made by Tugce: 

                % The variable 'events' is a struct.  Both fields have the original sampling rate.
                %                                     
                %   - The events begin approximately in 180 seconds (baseline period).
                %   - 'EVT05' is the label where CS is applied but depends on the plexons settings.

                data.events{1,1}  = ft_read_event(fullFileName);   % Data events

                % elseif  parameters.downsampling

            end


        % Case for Load *.mat. Behavior events analyzed from Guide_Video_Track

        case '.mat'

            fprintf(1, '\nExtracting %s\n', 'Behavior events analyzed from Guide_Video_Track');

            data.events.behavior = load(fullFileName);
            % correcting TS_LFPsec to zero if the record has started after a viewing time
            data.events.behavior.TS_LFPsec = data.events.behavior.TS_LFPsec - parameters.const_time;

            
        % Case for load *.xml. Impedance values

        case '.xml'

            fprintf(1, '\nExtracting %s\n', 'Impedance Values');

            temp = xml2struct(fullFileName);

            for ii = 1:length(temp.CHANNEL_IMPEDANCES.CHANNEL)
                parameters.impedance(ii,jj) = str2double(temp.CHANNEL_IMPEDANCES.CHANNEL{1, ii}.Attributes.magnitude);

            end

            % Delete zero columns
            parameters.impedance(:, ~any(parameters.impedance,1)) = [];

    end

end

%% Save
% newStr1 = files.id(ms).name(1:2);
% newStr2 = files.id(ms).name(end-2:end);

% name = strcat('E:\Projetos 2\Flavio\Samir\Analysis\Terceiro dia\',newStr1,newStr2,'_sFFT_FullTrial_stats.mat');
% name = strcat('/Users/flavio/Google Drive/MATLAB/my_functions/PNPD/',newStr1,newStr2,'_sFFT_FullTrial_stats.mat');
%
% save(name,'data','parameters','-v7.3')
%
% clear('name','variable','newStr1','newStr2')
%%
fprintf('\n Done. \n');

end

%% last update:
%  listening: 
