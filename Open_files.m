function [files] = Open_files()

% Organizes files path from Intan/Open Ephys/Plexon:  *.continuous, *.events, *.oebin, *.mat, *.xml, *.nex

% Two possibilities: - Load individual and selected files
%                    - Load all files in selected folders

% - Outputs:

%   "files"      ->  id : General folders and/or files information and settings
%                ->  FilesLoads: Cell with files information and settings


% by Flavio Mourao.

%% - Request from user -

prompt        = {'To load individual files type (1). To load a folder with subfolders type (2)';'To load bevaviour timestamps *.raw file type (1) or not type (0)'};
dlgtitle      = 'Please enter';
dims          = [1 80];
default_input = {'1';'0'};

input = inputdlg(prompt,dlgtitle,dims,default_input); %gui

%% Load files
    % *.continuous -> Signal
    % *.events -> Events
    % *.mat. Behavior events analyzed from Guide_Video_Track

switch input{1, 1} % string
    
    case '1'
        
        % Load individual selected files
        [FilesLoaded,Path] = uigetfile({'*.continuous; *.events; *.oebin; *.mat; *.xml; *.nex; *.nex5'},...
            'Select individual channels from an experimental animal','MultiSelect', 'on');

        if isequal(FilesLoaded,0)
           fprintf(1, '\nSelected Canceled. At least one file must be selected\n'); 
        end
        
        % Define a struct with files informations from dir organization'
        % BEWARE ! This organization changes according to the operating system.
        % files.FilesLoaded = repmat(struct('name',[],'folder',[],'date',[],'bytes',[],'isdir',[],'datenum',[]), 1, length(FilesLoaded));

        % Filename can be organize as a single char or a group char in a cell depending on the number os files selected
        if ischar(FilesLoaded)
            files.FilesLoaded = {dir(fullfile(Path, FilesLoaded))}; % condition for a single file selected
            files.id = struct('name',files.FilesLoaded{1,1}.name);
            
        else
            for ff = 1:length(FilesLoaded) % loop over multiple files selected
                files.FilesLoaded{1,1}(ff) = dir(fullfile(Path, char(FilesLoaded(ff))));
                files.id(ff) = struct('name',files.FilesLoaded{1,1}(ff).name);
            end
        end

        % Optional - Uncomment the line below for sort data. Channels based on a specific file properties.
        % data.Channels = nestedSortStruct(files.FilesLoaded,'name',1); % Perform a nested sort of a struct array based on multiple fields.
        % >>> https://uk.mathworks.com/matlabcentral/fileexchange/28573-nested-sort-of-structure-arrays?focused=5166120&tab=function
    
        if str2num(input{2, 1}) == 1
            % Load individual selected files
            [FilesLoaded_behav,Path_behav] = uigetfile({'*.mat;'},...
                'Select individual bevaviour file(s) from an experimental animal','MultiSelect', 'on');
            
            if isequal(FilesLoaded_behav,0)
                fprintf(1, '\nSelected Canceled. At least one file must be selected\n');
            end
            
            % Define a struct with files informations from dir organization'
            % BEWARE ! This organization changes according to the operating system.
            % files.FilesLoaded = repmat(struct('name',[],'folder',[],'date',[],'bytes',[],'isdir',[],'datenum',[]), 1, length(FilesLoaded));
            
            % Filename can be organize as a single char or a group char in a cell depending on the number os files selected
            if ischar(FilesLoaded_behav)
                files.FilesLoaded_behav = {dir(fullfile(Path_behav, FilesLoaded_behav))}; % condition for a single file selected
                files.id_behav = struct('name',files.FilesLoaded_behav{1,1}.name);
                
            else
                for ff = 1:length(FilesLoaded_behav) % loop over multiple files selected
                    files.FilesLoaded_behav{1,1}(ff) = dir(fullfile(Path_behav, char(FilesLoaded_behav(ff))));
                    files.id_behav(ff) = struct('name',files.FilesLoaded_behav{1,1}(ff).name);
                end
            end
            
            % Optional - Uncomment the line below for sort data. Channels based on a specific file properties.
            % data.Channels = nestedSortStruct(files.FilesLoaded,'name',1); % Perform a nested sort of a struct array based on multiple fields.
            % >>> https://uk.mathworks.com/matlabcentral/fileexchange/28573-nested-sort-of-structure-arrays?focused=5166120&tab=function
        end
        
        
        
    case '2'      
        
        % Load current folder
        Path = uigetdir;        
        folders = dir(Path);
        folders = folders(~ismember({folders.name}, {'.', '..'})); % Just removing the . and .. entries. Correspond to the current folder and the parent folder 
        
        % Load selected folders ans files from the list
        [indx,tf] = listdlg('PromptString',{'Select Folders and Files.',''},...
            'SelectionMode','multiple','ListString',{folders.name});        
        
        if tf == 0
           fprintf(1, '\nAt least one file or folder must be selected\n'); 
        
        else   
            files.id = folders(indx);

            files.FilesLoaded     = cell(1,length(files.id));        
            for kk = 1:length(files.id)
                
                if files.id(kk).isdir == 1
                    AllFiles_params = dir(fullfile(Path,files.id(kk).name,'**'));  % dir to folders and subfolders
                    AllFiles_params = AllFiles_params([AllFiles_params.isdir]==0); % set only files from the dir list
                    AllFiles_params = AllFiles_params(contains({AllFiles_params.name},{'.continuous','all_channels.events','oebin','mat','xml','nex','nex5'})); % select desired files 
                else
                    AllFiles_params = dir(fullfile(Path,files.id(kk).name)); % dir to single files
                
                end
                
                files.FilesLoaded{1,kk} = AllFiles_params;
                
            end
        end
        
        if str2num(input{2, 1}) == 1
            
            % Load current folder
            Path_behav = uigetdir;
            folders_behav = dir(Path_behav);
            folders_behav = folders_behav(~ismember({folders_behav.name}, {'.', '..'})); % Just removing the . and .. entries. Correspond to the current folder and the parent folder
            
            % Load selected folders ans files from the list
            [indx_behav,tf_behav] = listdlg('PromptString',{'Select Behaviour Folders and Files.',''},...
                'SelectionMode','multiple','ListString',{folders_behav.name});
            
            if tf_behav == 0
                fprintf(1, '\nAt least one file or folder must be selected\n');
                
            else
                files.id_behav = folders_behav(indx_behav);
                
                files.FilesLoaded_behav     = cell(1,length(files.id_behav));
                for kk = 1:length(files.id_behav)
                    
                    if files.id_behav(kk).isdir == 1
                        AllFiles_params_behav = dir(fullfile(Path_behav,files.id_behav(kk).name,'**'));  % dir to folders and subfolders
                        AllFiles_params_behav = AllFiles_params_behav([AllFiles_params_behav.isdir]==0); % set only files from the dir list
                        AllFiles_params_behav = AllFiles_params_behav(contains({AllFiles_params_behav.name},{'*.mat'})); % select desired files
                    else
                        AllFiles_params_behav = dir(fullfile(Path_behav,files.id_behav(kk).name)); % dir to single files
                        
                    end
                    
                    files.FilesLoaded_behav{1,kk} = AllFiles_params_behav;
                    
                end
            end
        end
            
        
fprintf('\n Load Path Done. \n');

end

%% last update 10/12/2023 - 
%  listening: 