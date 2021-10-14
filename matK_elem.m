function [Kel] = matK_elem(S1, S2, S3)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mat_elem :
% calcul la matrices de raideur elementaire en P1 lagrange
%
% SYNOPSIS [Kel] = mat_elem(S1, S2, S3)
%          
% INPUT * S1, S2, S3 : les 2 coordonnees des 3 sommets du triangle 
%                      (vecteurs reels 1x2)
%
% OUTPUT - Kel matrice de raideur elementaire (matrice 3x3)
%
% NOTE (1) le calcul est exacte (pas de condensation de masse)
%      (2) calcul direct a partir des formules donnees par 
%          les coordonnees barycentriques 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global epsilon;

% preliminaires, pour faciliter la lecture:
x1 = S1(1); y1 = S1(2);
x2 = S2(1); y2 = S2(2);
x3 = S3(1); y3 = S3(2);

% les 3 normales a l'arete opposees (de la longueur de l'arete)
norm = zeros(3, 2);
norm(1, :) = [y2-y3, x3-x2];
norm(2, :) = [y3-y1, x1-x3];
norm(3, :) = [y1-y2, x2-x1];

% D est, au signe pres, deux fois l'aire du triangle
D = ((x2-x1)*(y3-y1) - (y2-y1)*(x3-x1));
if (abs(D) <= eps) 
  error('l aire d un triangle est nulle!!!'); 
end;

w_hat = [-9/32,25/96,25/96,25/96];
M_hat = [1/3,1/3;1/5,1/5;1/5,3/5;3/5,1/5];
dF = [x2-x1,x3-x1;y2-y1,y3-y1];
gradW = [-1,1,0;-1,0,1];
% calcul de la matrice de raideur
% -------------------------------
Kel = zeros(3,3);
for i=1:3
  for j=1:3
    % Kel(i,j) = norm(i,:)*norm(j,:)'/(2*D);
    for k=1:4
        Kel(i,j) = Kel(i,j) + w_hat(k)*(A((x1+M_hat(k,1)*(x2-x1)+M_hat(k,2)*(x3-x1))/eps,(y1+M_hat(k,1)*(y2-y1)+M_hat(k,2)*(y3-y1))/eps)*((dF')^(-1)*gradW(:,i)))'*((dF')^(-1)*gradW(:,j))*abs(D);
    end
  end; % j
end; % i

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                        fin de la routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
