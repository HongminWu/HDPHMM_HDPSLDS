function [likelihood log_normalizer] = compute_likelihood(data_struct,theta,obsModelType,Kz,Ks)

switch obsModelType
   
    case 'Gaussian'
        
        invSigma = theta.invSigma;
        mu = theta.mu;
        
        T = size(data_struct.obs,2);
        dimu = size(data_struct.obs,1);
        
        log_likelihood = zeros(Kz,Ks,T);
        for kz=1:Kz  %hidden state
            for ks=1:Ks %nComponents
                
                cholinvSigma = chol(invSigma(:,:,kz,ks));
                dcholinvSigma = diag(cholinvSigma); %variance for each dimension
                
                u = cholinvSigma*(data_struct.obs - mu(:,kz*ones(1,T),ks));
                
                log_likelihood(kz,ks,:) = -0.5*sum(u.^2,1) + sum(log(dcholinvSigma));
            end
        end
        
        %max(log_likelihood,[],1) find the maximum values of each column
        %max(log_likelihood,[],2) find the maximum values of each row
         log_normalizer = max(max(log_likelihood,[],1),[],2);    
         log_likelihood = log_likelihood - log_normalizer(ones(Kz,1),ones(Ks,1),:); % minus the maximum value, for the sake of log_likelihood are negative
         likelihood = exp(log_likelihood);
         log_normalizer = log_normalizer - (dimu/2)*log(2*pi);
        
    case {'AR','SLDS'}
        
        invSigma = theta.invSigma;
        A = theta.A;
        X = data_struct.X;
        
        T = size(data_struct.obs,2);
        dimu = size(data_struct.obs,1);
        
        log_likelihood = zeros(Kz,Ks,T);
        if isfield(theta,'mu')
            
            mu = theta.mu;
            
            for kz=1:Kz
                for ks=1:Ks
                    
                    cholinvSigma = chol(invSigma(:,:,kz,ks));
                    dcholinvSigma = diag(cholinvSigma);
                    
                    u = cholinvSigma*(data_struct.obs - A(:,:,kz,ks)*X-mu(:,kz*ones(1,T),ks));
                    
                    log_likelihood(kz,ks,:) = -0.5*sum(u.^2,1) + sum(log(dcholinvSigma));
                    
                end
            end
        else
            
            for kz=1:Kz
                for ks=1:Ks
                    
                    cholinvSigma = chol(invSigma(:,:,kz,ks));
                    dcholinvSigma = diag(cholinvSigma);
                    
                    u = cholinvSigma*(data_struct.obs - A(:,:,kz,ks)*X);
                    
                    log_likelihood(kz,ks,:) = -0.5*sum(u.^2,1) + sum(log(dcholinvSigma));
                    
                end
            end
            
        end
        
       
        log_normalizer = max(max(log_likelihood,[],1),[],2);
        log_likelihood = log_likelihood - log_normalizer(ones(Kz,1),ones(Ks,1),:);
        likelihood = exp(log_likelihood);
       
        log_normalizer = log_normalizer - (dimu/2)*log(2*pi);
       
    case 'Multinomial'

        likelihood = theta.p(:,:,data_struct.obs);
        log_normalizer = zeros(1,size(data_struct.obs,2));
       
end