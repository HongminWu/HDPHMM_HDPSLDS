    clc;
    clear global;
    global_variables;
    cd(rootPath); 
    
    LearnedModelPath  = strcat(modelPath, 'TASK');
    if (exist(LearnedModelPath,'dir') == 0)
        mkdir(LearnedModelPath);
    end
    cd (LearnedModelPath);
    delete *;
    %trialIdx = {'02','03','04','05','06','07','08','09','10','11','12','13','14','15','16','17','18','19','20','21'};
    trialIdx = {'05','06','07','08','09','10'};
    for  idx=1:length(trialIdx) % for testing each state
        train_task(trialIdx{idx}, TRAINING_DATASET_PATH, LearnedModelPath);
    end
   cd (rootPath)
