% =====================================================
%
%
% une routine pour la mise en oeuvre des EF P1 Lagrange
% pour l'equation de Laplace suivante, avec conditions periodiques
% sur le maillage nom_maillage.msh
%
% | -div(A grad u ) u + u= f,   dans \Omega
% |         u periodique,   sur le bord
%
% =====================================================

clear();

global epsilon; epsilon = 5e-1;
global eta; eta = 1e-1;
global ast; ast = 'non';

% lecture du maillage et affichage
% ---------------------------------
nom_maillage = 'geomCarre_cell.msh';
[Nbpt,Nbtri,Coorneu,Refneu,Numtri,~,~,~,~]=lecture_msh(nom_maillage);
% Refneu = Refneu~=0;
% ----------------------
% calcul des matrices EF
% ----------------------

% declarations
% ------------
KK = sparse(Nbpt,Nbpt); % matrice de rigidite
KKb = sparse(Nbpt,Nbpt); % matrice de rigidite avec A = 1
MM = sparse(Nbpt,Nbpt); % matrice de masse
LL1 = zeros(Nbpt,1);     % vecteur second membre
LL2 = zeros(Nbpt,1);     % vecteur second membre

% boucle sur les triangles
% ------------------------
for l=1:Nbtri
  % Coordonnees des sommets du triangles
  S1=Coorneu(Numtri(l,1),:);
  S2=Coorneu(Numtri(l,2),:);
  S3=Coorneu(Numtri(l,3),:);
  % calcul des matrices elementaires du triangle l 
  
   Kel=matK_elem(S1, S2, S3);
   Mel=matM_elem(S1, S2, S3);
   Kbel=matKb_elem(S1, S2, S3);
    
  % On fait l'assemmblage de la matrice globale et du second membre
  for i=[1:3]
      for j=[1:3]
          MM(Numtri(l,i),Numtri(l,j))=MM(Numtri(l,i),Numtri(l,j))+Mel(i,j);
          KK(Numtri(l,i),Numtri(l,j))=KK(Numtri(l,i),Numtri(l,j))+Kel(i,j);
          KKb(Numtri(l,i),Numtri(l,j))=KKb(Numtri(l,i),Numtri(l,j))+Kbel(i,j);
      end
  end
end % for l

% Second Membre
LL1 = -KK*Coorneu(:,1);
LL2 = -KK*Coorneu(:,2);

% Projection sur l espace V_p
% matrice de projection 
N0 = Nbpt-sum(Refneu~=0);
NS = sum(Refneu==5);
X = eye(sum(Refneu==3));
X = X(end:-1:1,:);
PP = sparse([ones(1,4) zeros(1,sum(Refneu~=5));zeros(sum(Refneu==1),4) eye(sum(Refneu==1)) zeros(sum(Refneu==1),sum(Refneu==1)) X zeros(sum(Refneu==1),sum(Refneu==1)+N0);zeros(sum(Refneu==2),4+sum(Refneu==1)) eye(sum(Refneu==2)) zeros(sum(Refneu==2),sum(Refneu==3)) X zeros(sum(Refneu==2),N0);zeros(N0,sum(Refneu~=0)) eye(N0)]);
% spy(PP)
AA = eta*MM+KK;
AAp = PP*AA*PP';
LLp1 = PP*LL1;
LLp2 = PP*LL2;



% inversion
% ----------
WWp1 = AAp\LLp1;
WWp2 = AAp\LLp2;

% Expression de la solution dans toute la base

WW1 = PP'*WWp1;
WW2 = PP'*WWp2;

global A_ast;A_ast=[(Coorneu(:,1)+WW1)'*KK*(Coorneu(:,1)+WW1),(Coorneu(:,1)+WW1)'*KK*(Coorneu(:,2)+WW2);(Coorneu(:,2)+WW2)'*KK*(Coorneu(:,1)+WW1),(Coorneu(:,2)+WW2)'*KK*(Coorneu(:,2)+WW2)];
global ast; ast = 'oui';

nom_maillage = 'geomCarre01.msh';
[Nbpt,Nbtri,Coorneu,Refneu,Numtri,~,~,~,~]=lecture_msh(nom_maillage);



% declarations
% ------------
KK = sparse(Nbpt,Nbpt); % matrice de rigidite
KKb = sparse(Nbpt,Nbpt); % matrice de rigidite avec A = 1
MM = sparse(Nbpt,Nbpt); % matrice de masse
LL = zeros(Nbpt,1);     % vecteur second membre

% boucle sur les triangles
% ------------------------
for l=1:Nbtri
  % Coordonnees des sommets du triangles
  S1=Coorneu(Numtri(l,1),:);
  S2=Coorneu(Numtri(l,2),:);
  S3=Coorneu(Numtri(l,3),:);
  % calcul des matrices elementaires du triangle l 
  
   Kel=matK_elem(S1, S2, S3);
   Mel=matM_elem(S1, S2, S3);
   Kbel=matKb_elem(S1, S2, S3);
    
  % On fait l'assemmblage de la matrice globale et du second membre
  for i=[1:3]
      for j=[1:3]
          MM(Numtri(l,i),Numtri(l,j))=MM(Numtri(l,i),Numtri(l,j))+Mel(i,j);
          KK(Numtri(l,i),Numtri(l,j))=KK(Numtri(l,i),Numtri(l,j))+Kel(i,j);
          KKb(Numtri(l,i),Numtri(l,j))=KKb(Numtri(l,i),Numtri(l,j))+Kbel(i,j);
      end
  end
end % for l

% Calcul du second membre L
% -------------------------
FF = f(Coorneu(:,1),Coorneu(:,2));
LL = MM*FF;

% Projection sur l espace V_0
% matrice de projection 
N0 = Nbpt-sum(Refneu~=0);
PP = sparse([zeros(N0,sum(Refneu~=0)) eye(N0)]);
AA = KK;
AA0 = PP*AA*PP';
LL0 = PP*LL;

% inversion
% ----------
UU0 = AA0\LL0;

% Expression de la solution dans toute la base
UU = PP'*UU0;

% visualisation
% -------------
affiche(UU, Numtri, Coorneu, sprintf('Cellule - %s', 'geomCarre\_per.msh'));

validation = 'non';
% validation
% ----------
if strcmp(validation,'oui')
    UU_exact = sin(pi*Coorneu(:,1)+pi/2).*sin(pi*Coorneu(:,2)+pi/2);
    affiche(UU_exact, Numtri, Coorneu, sprintf('Cellule exact - %s', 'geomCarre\_per.msh'));
    % Calcul de l erreur L2
    EE_L2 = (UU_exact-UU)'*MM*(UU_exact-UU)
    log(sqrt(EE_L2/(UU_exact'*MM*UU_exact)));
    % Calcul de l erreur H1
    EE_H1 = (UU_exact-UU)'*(KKb)*(UU_exact-UU)
    log(sqrt(EE_H1/(UU_exact'*(KKb)*UU_exact)));
    % attention de bien changer le terme source (dans FF)
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                        fin de la routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%