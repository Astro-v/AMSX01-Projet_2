% =====================================================
%
%
% une routine pour la mise en oeuvre des EF P1 Lagrange
% pour l'equation de Laplace suivante, avec A 1-periodique
% et u_eps avec des condition de Dirichlet
% 
% =====================================================

clear();

global eps; eps = 5e-1;

% lecture du maillage et affichage
% ---------------------------------
nom_maillage = 'geomCarre_per.msh';
[Nbpt,Nbtri,Coorneu,Refneu,Numtri,Reftri,Nbaretes,Numaretes,Refaretes]=lecture_msh(nom_maillage);

% ----------------------
% calcul des matrices EF
% ----------------------

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

% Projection sur l espace V_p
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
affiche(UU, Numtri, Coorneu, sprintf('Micro periodique - %s', 'geomCarre\_per.msh'));

validation = 'non';
% validation
% ----------
if strcmp(validation,'oui')
UU_exact = sin(pi*Coorneu(:,1)).*sin(pi*Coorneu(:,2));
affiche(UU_exact, Numtri, Coorneu, sprintf('Micro periodique exact - %s', 'geomCarre\_per.msh'));
% Calcul de l erreur L2
EE_L2 = (UU_exact-UU)'*MM*(UU_exact-UU);
log(sqrt(EE_L2/(UU_exact'*MM*UU_exact)))
% Calcul de l erreur H1
EE_H1 = (UU_exact-UU)'*(KKb)*(UU_exact-UU);
log(sqrt(EE_H1/(UU_exact'*(KKb)*UU_exact)))
% attention de bien changer le terme source (dans FF)
end