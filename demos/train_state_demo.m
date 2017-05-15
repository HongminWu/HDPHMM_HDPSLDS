function train_state_demo()
    clc;
    clear global;
    global_variables;
    global STATE TRAINING_DATASET_PATH
    cd(rootPath); 
    
    [APPROACH, ROTATION, INSERTION, MATING] = deal(1,2,3,4);
    for training_state = APPROACH : MATING
        LearnedModelPath  = strcat(modelPath, char(STATE(training_state)));
        if (exist(LearnedModelPath,'dir') == 0)
            mkdir(LearnedModelPath);
        end
        cd (LearnedModelPath);
        delete *;
        for trialID = 7:14
            train_state(trialID, training_state, TRAINING_DATASET_PATH, LearnedModelPath);
        end
    end
   cd (rootPath)
end
