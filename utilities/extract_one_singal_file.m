%@Homls 2016-11-08
%Extract the same file in different folder, e.g   R_Torques.dat
%    ParentDir = 'SIM_HIRO_ONE_SA_SUCCESS';   the successful demonstrations
%   ParentDir = 'SIM_HIRO_ONE_SA_FAILURE';      the failure demonstrations
%the pwd = 'F:\Matlab\Emily B. Fox\BPARHMMtoolbox\logs'
function [DataCell, R_State, folders_name] = extract_one_singal_file(ParentDir, dataType, ndemos, idata) 
    DataCell = {};
    R_State = {};
    ParentDir = strcat('logs/',ParentDir);
    cd(ParentDir);
    
    %SIM_HIRO_ONE_SA_ERROR_CHARAC_Prob
    specificFailure = '+x';   
    %folder = {'FC008+x0.0090','FC011+x0.0105','FC150+x0.0105'};
    if strcmp(ParentDir,'logs/SIM_HIRO_ONE_SA_ERROR_CHARAC_Prob')
        cd (specificFailure);
        ParentDir = specificFailure;
    end
   folders = dir;

%id = linspace(6,6+ndemos-1,ndemos); %+y extract special demonstrations
%id = 4 + idata;
 %ntrain = 10; %for sim_success
 ntrain = 5; %for sim_failure
 while length(DataCell) ~= ntrain
    DataCell = {};
    R_State = {};
    id = randi([3,size(folders,1)],1,ntrain); % randomly extract the demonstrations
    folders_name = {};
    for i=1:size(id, 2) % a folder per demonsration      
        if(folders(id(i)).isdir == 1 && ~strcmp(folders(id(i)).name,'..') && ~strcmp(folders(id(i)).name,'.') && ~strcmp(folders(id(i)).name,'.git') && ~strcmp(folders(id(i)).name,'ERROR'))
            folders_name = [folders_name, {folders(id(i)).name}];
            foldname = folders(id(i)).name;
            cd(foldname); 
            data = [];
            iState = [];
            for j = 1: length(dataType)
                file = dir([char(dataType(j)) '.dat*']);
                filename = strcat('/',foldname,'/',file(1).name);
                raw_data = load(filename);         
                d = raw_data(:,2:end);  %delete the time column
                if strcmp(dataType(j), 'R_Torques')             
                    d = smoothWrenchFilter(d(:,1)); % wrenchVec = d;
                 %   d = unique(d,'rows','stable'); %first-order detective
                 %   d = [d,[d(1,:);diff(d)]];
                end
                data = [data, d];    
            end
            DataCell = [DataCell; {data}];
            
            file = dir([char('R_State') '.dat*']);
            filename = strcat('/',foldname,'/',file(1).name);
            tState = load(filename);
            for in = 1:length(tState)
                idex= find(raw_data(:,1) == tState(in));
                iState = [iState,idex];
            end    
            R_State = [R_State; {iState}];   
            
            cd ..
        end
    end
 end
%  nbData = [];
%  for i = 1: length(DataCell) 
%      DataCell{i} = DataCell{i}';
%      nbData = [nbData, length(DataCell{i})]; % so as to find out the max length
%  end  
%  data = [];
%  for n=1:length(DataCell)
%      DataCell{n} = spline(1:size(DataCell{n},2), DataCell{n}, linspace(1,size(DataCell{n},2), max(nbData))); %Resampling alignment
%      data = [data, DataCell{n}(:)]; 
%  end
%     DataCell = {};
%     DataCell = {data};
    cd ../..
    
    if strcmp(ParentDir,specificFailure)
        cd ..
    end
end