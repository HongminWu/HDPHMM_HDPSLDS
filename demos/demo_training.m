function demo_training()
    clc;
%     clear;
    clear global;
    global_variables;
    cd(rootPath); 
    disp('....................................WARNING.......................................................')
    disp('....................................WARNING.......................................................')
    disp('....................................WARNING.......................................................')
    disp('....................................WARNING.......................................................')
    disp('.........................Are you sure to re-training the model?...................................')
    
    [APPROACH, ROTATION, INSERTION, MATING] = deal(1,2,3,4);
    for training_state = APPROACH : MATING
        LearnedModelPath  = strcat(modelPath, char(STATE(training_state)));
        cd (LearnedModelPath);
        delete *;
        for trialID = 5:7 
            run(trialID, training_state, TRAINING_DATASET_PATH, LearnedModelPath);
        end
    end
end
