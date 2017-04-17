%@Homls 2016-11-08
%Extract the same file in different folder, e.g   R_Torques.dat
%    ParentDir = 'SIM_HIRO_ONE_SA_SUCCESS';   the successful demonstrations
%   ParentDir = 'SIM_HIRO_ONE_SA_FAILURE';      the failure demonstrations
%the pwd = 'F:\Matlab\Emily B. Fox\BPARHMMtoolbox\logs'
function [DataCell, R_State, folders_name] = extract_file(ParentDir, dataType, ndemos, idata) 
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

 %  randsel = randperm(size(folders,1)); % randomly extract the demonstrations
 %  id = randsel(1:ndemos);
 % id = [6];
 %id = (1:length(folders)); %extract all the demonstrations
 %id = [1, 3, 5, 6, 7, 10] + 2; %+r extract special demonstrations
 %id = [5, 6, 8, 9, 10, 11] + 2; %+x extract special demonstrations
 %id = [1, 2, 3, 4, 6, 10, 12] + 2; %+x+r extract special demonstrations 
 %id = [4, 5, 6, 7, 8, 9, 10, 12] + 2; %+x+y extract special demonstrations 
 %id = [                        ] + 2; %+x+y+r extract special demonstrations
 %id = [1, 2, 3, 4] + 2; %+x+y-r extract special demonstrations
 %id = [2, 3, 4, 5, 8, 9] + 2; %+x+y-r extract special demonstrations
 
 
%id = linspace(6,6+ndemos-1,ndemos); %+y extract special demonstrations
id = 12 + idata;
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
%                     d = smoothWrenchFilter(d); % wrenchVec = d;
                    d = unique(d,'rows','stable');
                    d = [d,[d(1,:);diff(d)]];
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
            iState = [iState,size(data,1)];
            R_State = [R_State; {iState}];
        
            cd ..
        end
    end
    cd ../..
    if strcmp(ParentDir,specificFailure)
        cd ..
    end
end