function train_state_demo()
    clc;
    clear global;
    global_variables;
    cd(rootPath); 
    
    [APPROACH, ROTATION, INSERTION, MATING] = deal(1,2,3,4);
    for training_state = APPROACH : MATING
        LearnedModelPath  = strcat(modelPath, char(STATE(training_state)));
        cd (LearnedModelPath);
        delete *;
        for trialID = 20:20
            train_state(trialID, training_state, TRAINING_DATASET_PATH, LearnedModelPath);
        end
    end
   cd (rootPath)
end
