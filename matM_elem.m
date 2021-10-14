function Mel = matM_elem(S1, S2, S3)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matM_elem :
% calcul la matrices de masse elementaire en P1 lagrange
%
% SYNOPSIS Mel = matM_elem(S1, S2, S3)
%          
% INPUT * S1, S2, S3 : les 2 coordonnees des 3 sommets du triangle 
%                      (vecteurs reels 1x2)
%
% OUTPUT - Mel matrice de masse elementaire (matrice 3x3)
%
% NOTE (1) le calcul est exacte (pas de condensation de masse)
%      (2) calcul direct a partir des formules donnees par 
%          les coordonnees barycentriques 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% preliminaires, pour faciliter la lecture:
x1 = S1(1); y1 = S1(2);
x2 = S2(1); y2 = S2(2);
x3 = S3(1); y3 = S3(2);

xc1 = (x1+x2)/2; yc1 = (y1+y2)/2;
xc2 = (x1+x3)/2; yc2 = (y1+y3)/2;
xc3 = (x3+x2)/2; yc3 = (y3+y2)/2;

lambda = [(y2-y3)*(xc1-x3)-(x2-x3)*(yc1-y3),(y2-y3)*(xc2-x3)-(x2-x3)*(yc2-y3),(y2-y3)*(xc3-x3)-(x2-x3)*(yc3-y3);(y3-y1)*(xc1-x1)-(x3-x1)*(yc1-y1),(y3-y1)*(xc2-x1)-(x3-x1)*(yc2-y1),(y3-y1)*(xc3-x1)-(x3-x1)*(yc3-y1);(y1-y2)*(xc1-x2)-(x1-x2)*(yc1-y2),(y1-y2)*(xc2-x2)-(x1-x2)*(yc2-y2),(y1-y2)*(xc3-x2)-(x1-x2)*(yc3-y2)];


% D est, au signe pres, deux fois l'aire du triangle
D = ((x2-x1)*(y3-y1) - (y2-y1)*(x3-x1));
if (abs(D) <= eps) 
  error('l aire d un triangle est nulle!!!'); 
end;

% calcul de la matrice de masse
% -----------------------------
Mel = zeros(3,3);
for i=1:3
	for j=1:3
		% Mel(i,j) = (1/(6*D))*(lambda(i,1)*lambda(j,1)+lambda(i,2)*lambda(j,2)+lambda(i,3)*lambda(j,3));
        if i==j
            Mel(i,j) = D/12;
        else
            Mel(i,j) = D/24;
        end
    end; % j
end; % i

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                        fin de la routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
