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

% Vertices and edges of a hexagonal cylinder
th=(0:60:359)/180*pi;
S2=[sin(th);cos(th)];
n=length(th);
edg2=[ [1:n]' [2:n 1]'];
S=[ [S2 S2]; kron([-1 1],ones(1,n))];
edg=[edg2 ; edg2+n; [1:n; (1:n)+n]'];

%-- select the viewing angle
AZ=-20; EL=26;

%-- construct the projection matrix

M=viewmtx(AZ,EL); M=M(1:2,:);

%-- compute the 2D image
I=M*[S; ones(1,size(S,2))];

% Make an initial random guess, Sinit, of the true underlying scene
Sinit=rand(size(S));

% Use an optimization routine to find a better guess minimizing
% the error in form of the negative logarithm of the posterior
options = optimset('MaxFunEvals',1000000,'TolFun',1e-3,'TolX',1e-3);
Shat = fminunc(@HexaError,Sinit,options);

% Draw the true underlying scene, S, and the best fit, Shat.
figure(1); clf;
plotscene(S,edg,'ro-');
hold on; plotscene(Shat,edg,'bo-'); hold off; axis equal; rot(AZ,EL);
title('3D Scene, maximum likelihood fit');axis off;



    function NegLogPost=HexaError(Sguess)

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
        
        
        I_proj = M*[Sguess; ones(1,size(Sguess,2))]; %projection of S_guess initial on the 2-D image plane
        sigma_N = 2; % initial assumption
        sigma_S = 2;% initial assumption
        n = 12;
         % i have 12 vertices and i want to calculate the distance between
         % each pair
        for i = 1 : n
            for j = i+1 : n
            %L = pdist2(S(:,i),S(:,i))
            %l(i) = norm(L)
            l(i,j) = norm(S(:,i)-S(:,j));
            end
        end
        
        for i = 1 : n
            for j = i+1 : n
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
        %NegLogPost= Likelihood+ Prior_position; % easy game 
        
        %NegLogPost= Likelihood+ Prior_angles +Prior_lenght+ Prior_position;
        %redundant
        
        %NegLogPost= Likelihood+ Prior_lenght; %not enough
        
        NegLogPost= Likelihood+ Prior_angles +Prior_lenght; 
         %with 8 vertexes not enough
         %with 12 vertexes, much better
         
        %NegLogPost= Likelihood+ Prior_angles; not enough
        %NegLogPost= Likelihood; %not enough info to have a correct Shat
          
    end
end