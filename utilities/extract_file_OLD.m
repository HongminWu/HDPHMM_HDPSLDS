%@Homls 2016-11-08
%Extract the same file in different folder, e.g   R_Torques.dat
%    ParentDir = 'SIM_HIRO_ONE_SA_SUCCESS';   the successful demonstrations
%   ParentDir = 'SIM_HIRO_ONE_SA_FAILURE';      the failure demonstrations
%the pwd = 'F:\Matlab\Emily B. Fox\BPARHMMtoolbox\logs'
function DataCell = extract_file(ParentDir, FILE_NAME, nTrail)
    DataCell = {};
    cd(ParentDir);
    folders = dir;
    randsel = randperm(size(folders,1));
    id = randsel(1:nTrail);
    for i=1:size(id, 2) % a folder per demonsration      
        if(folders(id(i)).isdir == 1 && ~strcmp(folders(id(i)).name,'..') && ~strcmp(folders(id(i)).name,'.') && ~strcmp(folders(id(i)).name,'.git'))
            foldname = strcat(ParentDir,'/',folders(id(i)).name);
            cd ..
            cd(foldname); 
            file = dir(FILE_NAME);
            filename = strcat(foldname,'/',file.name);
            d = load(filename);
            d = d(:,2:end);  %delete the time column
            DataCell = [DataCell; {d}];
            cd ..
        end
    end
    cd ..
end

%     for i=1:size(folders) % a folder per demonsration      
%         if(folders(i).isdir == 1 && ~strcmp(folders(i).name,'..') && ~strcmp(folders(i).name,'.') && ~strcmp(folders(i).name,'.git'))
%             foldname = strcat(ParentDir,'/',folders(i).name);
%             cd ..
%             cd(foldname); 
%             file = dir(FILE_NAME);
%             filename = strcat(foldname,'/',file.name);
%             d = load(filename);
%             d = d(:,2:end);  %delete the time column
%             DataCell = [DataCell; {d}];
%             cd ..
%         end
%     end
%     cd ..
% end