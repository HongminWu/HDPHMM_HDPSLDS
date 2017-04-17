function vis_neglog_c(neglog_c,dataset,idata,path)
        neglog_c = cumsum(neglog_c);
        gHandle = figure;
        plot(neglog_c); 
        title('Cumulative class-specific log-likelihood');
        if strcmp(dataset, 'SIM_HIRO_ONE_SA_SUCCESS') || strcmp(dataset,'REAL_HIRO_ONE_SA_SUCCESS')
            fname = strcat('success_',num2str(idata));
            cd(path);
        elseif strcmp(dataset, 'SIM_HIRO_ONE_SA_FAILURE')
            fname = strcat('failure_',num2str(idata));
            cd(path);
        end
        saveas(gHandle, fname,'jpg');
        cd ../..
end