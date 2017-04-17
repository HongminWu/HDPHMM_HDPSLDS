function DataCell = pick_n_place_file()
    DataCell = {};
    csv_file_name = 'wrench-filtered-right.csv';
    ParentDir = './logs/pick_n_place/pick_n_place_box_random_smach/sucess';
    cd(ParentDir);
    folders = dir;
  % randsel = randperm(size(folders,1)); % randomly extract the demonstrations
  % id = randsel(1:nTrail);
  id = [3];
    for i=1:size(id, 2) % a folder per demonsration      
        if(folders(id(i)).isdir == 1 && ~strcmp(folders(id(i)).name,'..') && ~strcmp(folders(id(i)).name,'.') && ~strcmp(folders(id(i)).name,'.git') && ~strcmp(folders(id(i)).name,'ERROR'))
            cd(folders(id(i)).name);   
            csv_data = csvread(csv_file_name,1,5);
            smooth_data = smoothWrenchFilter(csv_data); % wrenchVec = d;   
            DataCell = [DataCell; {smooth_data}];
            cd ..
        end
    end
    cd ../../../..
end