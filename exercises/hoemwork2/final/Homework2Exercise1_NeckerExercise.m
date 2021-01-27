function NeckerExercise

%----------------------------------------------------------------------
% 02458 Cognitive Modelling - Necker Exercise
% Orignial version by Jouko Lampinen, Helsinki University of Technology
% Modified by Tobias Andersen, Technical University of Denmark
%
% Variables:
%  S   - 3D vertex list of the true structure, as 3xN matrix
%  edg - edgea list as Mx2 matrix, each row contains the indices of
%        the vertices linked by the edge
%  I   - projection of S to the image plane, as 2xN matrix
%  M   - projection matrix to project S to I,
%          I=M*[S; ones(1,size(S,2))];
%        M can be constructed by
%          M=viewmtx(AZ,EL);
%	 where Az is the azimuth and EL elevation angle of the observer.
%        Note: for orthographic projection only 3 columns of M are needed,
%        so that I=M(1:2,1:3)*S.
%        For perspective projection (M=viewmtx(AZ,EL,PHI)) M is 2 by 4 and
%        the homogenous coordinates [S;ones(1,size(S,2))] are used.
%        For generality we use full projection matrix here.
%
% The task is to estimate Shat below given the observed image I

%-- vertices of the true underlying scene
S=[ 0 1 1 0 0 1 1 0 ;0 0 1 1 1 1 0 0;0 0 0 0 1 1 1 1]*2-1;

%-- edges as start point, end point index pairs
edg=[ 1 2; 1 4 ; 1 8 ; 2 3 ; 2 7 ; 3 4 ; 3 6 ; 4 5 ; 5 6 ; 5 8 ; 6 7 ; 7 8];

%-- select the viewing angle
AZ=-32; EL=25;

%-- construct the projection matrix

M=viewmtx(AZ,EL); M=M(1:2,:);

%-- compute the 2D image I_proj
I=M*[S; ones(1,size(S,2))];

% Make an initial random guess, Sinit, of the true underlying scene
Sguess=rand(size(S)); %S_guess? 

% Use an optimization routine to find a better guess minimizing
% the error in form of the negative logarithm of the posterior
options = optimset('MaxFunEvals',1000000,'TolFun',1e-3,'TolX',1e-3);
Shat = fminunc(@NeckerError,Sguess,options);

% Draw the true underlying scene, S, and the best fit, Shat.
figure(1); clf;
plotscene(S,edg,'ro-');
hold on; plotscene(Shat,edg,'bo-'); hold off; axis equal; rot(AZ,EL);
title('3D Scene, best fit');axis off;


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
        sigma_N = 2; % initial assumption
        sigma_S = 2;% initial assumption
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
        
        
        Likelihood = norm(I-I_proj)/sigma_N.^2;
        Prior_angles = norm(90-anglist)/sigma_S.^2;
        Prior_lenght = norm (l-lguess)/sigma_S.^2;
        Prior_position = norm (S-Sguess)/sigma_S.^2;
        NegLogPost= Likelihood+ Prior_position;
        %NegLogPost= Likelihood+ Prior_angles +Prior_lenght+ Prior_position;
        %NegLogPost= Likelihood+ Prior_angles +Prior_lenght;
        %NegLogPost= Likelihood+ Prior_angles;
        %NegLogPost= Likelihood;
    end

end