function train_task_demo()
    clc;
    clear global;
    global_variables;
    cd(rootPath); 
    
    LearnedModelPath  = strcat(modelPath, 'TASK');
        cd (LearnedModelPath);
        delete *;
        for train_task_id = 20:20
            train_task(train_task_id, TRAINING_DATASET_PATH, LearnedModelPath);
        end
   cd (rootPath)
end
