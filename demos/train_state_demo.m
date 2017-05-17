function train_state_demo()
    clc;
    clear global;
    global_variables;
    global STATE TRAINING_DATASET_PATH
    cd(rootPath); 
    
    [APPROACH, ROTATION, INSERTION, MATING] = deal(1,2,3,4);
    for training_state = ROTATION : MATING
        LearnedModelPath  = strcat(modelPath, char(STATE(training_state)));
        if (exist(LearnedModelPath,'dir') == 0)
            mkdir(LearnedModelPath);
        end
        cd (LearnedModelPath);
        delete *;
        trialIdx = {'02','03','04','05','06','07','08','09','10','11','12','13','14','15','16','17','18','19','20','21'};
        for  idx=1:length(trialIdx) % for testing each state
            train_state(trialIdx{idx}, training_state, TRAINING_DATASET_PATH, LearnedModelPath);
        end
    end
    
    calculate_state_threshold; % Optimal model selection and Calculate the expected log-likelihood
    
    test_state_demo;
    
   cd (rootPath)
end
