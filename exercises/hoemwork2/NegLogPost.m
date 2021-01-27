function NegLogPost=NeckerError(Sguess)

        %This is the error function that must return the negative logarithm
        %of the posterior probability
        
        %You can use all the variables from the calling function
        %NeckerExercise here to help you       

        %-- Here is help code for computing the angles between all edges
        anglist=[];
        for k=1:size(Sguess,2),
            con=[edg(find(edg(:,1)==k),2); edg(find(edg(:,2)==k),1)]';
            anglist=[anglist; [ k con([1 2]); k con([1 3]); k con([2 3])]];
        end


        u1=Sguess(:,anglist(:,2))-Sguess(:,anglist(:,1));
        u2=Sguess(:,anglist(:,3))-Sguess(:,anglist(:,1));

        u1=u1./repmat(sqrt(sum(u1.^2)),3,1);
        u2=u2./repmat(sqrt(sum(u2.^2)),3,1);
        Angles=acos(sum(u1.*u2))*180/pi;
        %%-- end of help code
        
        % arg_min[((I-I_proj).^2/(sigma_N^2))+((S-S_hat).^2/(sigma_S^2))] 
        % scene S = normal distribution N (mean_S,sigma_S)
        % noise N = normal distribution N (mean_N = 0, sigma_N=1)
        % I=M*[S; ones(1,size(S,2))];
        % I_guess = M*[Sguess; ones(1,size(Sguess,2))];
        % 1st part
        
        I_proj = M*[Sguess; ones(1,size(Sguess,2))]; %projection of S_guess initial on the 2-D image plane
        I_proj_r = Mr*[Sguess; ones(1,size(Sguess,2))];
        sigma_N = 1; % initial assumption for the noise N(0,1)
        sigma_S = 1;% initial assumption, usually given by tests
         % i have 8 vertices and i want to calculate the distance between
         % each pair
        for i = 1 : 8
            for j = i+1 : 8
            %L = pdist2(S(:,i),S(:,i))
            %l(i) = norm(L)
            l(i,j) = norm(S(:,i)-S(:,j));
            end
        end
        
        for i = 1 : 8
            for j = i+1 : 8
            %L = pdist2(S(:,i),S(:,i))
            %l(i) = norm(L)
            lguess(i,j) = norm(Sguess(:,i)-Sguess(:,j));
            end
        end
        %l = norm(S(:,1)-S(:,2)) %lenght of first side of the real cube (2)       
        
        
        %Likelihood = norm(I-I_proj)/sigma_N.^2;
        Prior_angles = norm(90-anglist)/sigma_S.^2;
        Prior_lenght = norm (l-lguess)/sigma_S.^2;
        Prior_position = norm (S-Sguess)/sigma_S.^2;
        %NegLogPost= Likelihood+ Prior_position;
        %NegLogPost= Likelihood+ Prior_angles +Prior_lenght+ Prior_position;
        %NegLogPost= Likelihood+ Prior_angles +Prior_lenght;
        %NegLogPost= Likelihood+ Prior_angles;
        %NegLogPost= Likelihood;
        
        Likelihood = norm(I-I_proj)/sigma_N.^2+norm(I-I_proj_r)/sigma_N.^2;
        NegLogPost= Likelihood;
    end