%ex.3

function NeckerExercise

%----------------------------------------------------------------------
% 02458 Cognitive Modelling - Necker Exercise WITH TEO EYES
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

%-- select the viewing angles
AZ1=-32; EL1=25; %left eye
AZ2=+32; EL2=25; %right eye

%-- construct the projection matrixs
M1=viewmtx(AZ1,EL1); M1=M1(1:2,:); %left eye
M2=viewmtx(AZ2,EL2); M2=M2(1:2,:); %right eye

sigma_N = 1; % initial assumption for the noise N(0,1)
sigma_S_1 = 1;% initial assumption, usually given by tests
sigma_S_2 = 1;% initial assumption, usually given by tests
r1 = 1/sigma_S_1; %reliability of view as opposite to the standard deviation
r2 = 1/sigma_S_2;
%-- compute the 2D image real I_proj

%I=M*[S; ones(1,size(S,2))];
I1=M1*[S; ones(1,size(S,2))]; % from left eye
I2=M2*[S; ones(1,size(S,2))]; % from right eye


% Make an initial random guess, Sinit, of the true underlying scene
Sguess_1=rand(size(S));%S_guess_left_eye
Sguess_2=rand(size(S));%S_guess_right_eye
Sguess = Sguess_1*(r1/(r1+r2))+Sguess_2*(r2/(r1+r2))



% Use an optimization routine to find a better guess minimizing
% the error in form of the negative logarithm of the posterior
options = optimset('MaxFunEvals',1000000,'TolFun',1e-3,'TolX',1e-3);
Shat = fminunc(@NeckerError,Sguess,options);

% Draw the true underlying scene, S, and the best fit, Shat.
figure
subplot(1,2,1)
%figure(1); clf;
plotscene(S,edg,'ro-');
hold on; plotscene(Shat,edg,'bo-'); hold off; axis equal; rot(AZ1,EL1);
title('3D Scene, left eye prespective tot,best fit');axis off;
subplot(1,2,2)
%figure(2); clf;
plotscene(S,edg,'ro-');
hold on; plotscene(Shat,edg,'bo-'); hold off; axis equal; rot(AZ2,EL2);
title('3D Scene, right eye prespective tot, best fit');axis off;


    function NegLogPost=NeckerError(Sguess)
        
        I_proj_1 = M1*[Sguess; ones(1,size(Sguess,2))]; %projection of S_guess initial on the 2-D image plane left eye
        I_proj_2 = M2*[Sguess; ones(1,size(Sguess,2))]; %projection of S_guess initial on the 2-D image plane right eye
        
        Prior_position_1 = norm (S-Sguess_1)/sigma_S_1.^2;
        Prior_position_2 = norm (S-Sguess_2)/sigma_S_2.^2;

        
        Likelihood = (norm(I1-I_proj_1)/sigma_N.^2)+norm(I2-I_proj_2)/sigma_N.^2; %we dont need priors to have a well defined sense of space and depth
        NegLogPost= Likelihood;
        %NegLogPost= Likelihood + Prior_position_1+ Prior_position_r;
    end

end