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


% lecture du maillage et affichage
% ---------------------------------
nom_maillage = 'geomCarre_per01.msh';
[Nbpt,Nbtri,Coorneu,Refneu,Numtri,Reftri,Nbaretes,Numaretes,Refaretes]=lecture_msh(nom_maillage);
% Refneu = Refneu~=0;
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
NS = sum(Refneu==5);
X = eye(sum(Refneu==3));
X = X(end:-1:1,:);
PP = [ones(1,4) zeros(1,sum(Refneu~=5));zeros(sum(Refneu==1),4) eye(sum(Refneu==1)) zeros(sum(Refneu==1),sum(Refneu==1)) X zeros(sum(Refneu==1),sum(Refneu==1)+N0);zeros(sum(Refneu==2),4+sum(Refneu==1)) eye(sum(Refneu==2)) zeros(sum(Refneu==2),sum(Refneu==3)) X zeros(sum(Refneu==2),N0);zeros(N0,sum(Refneu~=0)) eye(N0)];
AA = MM+KK;
AAp = PP*AA*PP';
LLp = PP*LL;



% inversion
% ----------
UUp = AAp\LLp;

% Expression de la solution dans toute la base

UU = PP'*UUp;

% visualisation
% -------------
affiche(UU, Numtri, Coorneu, sprintf('Periodique - %s', 'geomCarre\_per.msh'));

validation = 'oui';
% validation
% ----------
if strcmp(validation,'oui')
    UU_exact = sin(pi*Coorneu(:,1)+pi/2).*sin(pi*Coorneu(:,2)+pi/2);
    affiche(UU_exact, Numtri, Coorneu, sprintf('Periodique exact - %s', 'geomCarre\_per.msh'));
    % Calcul de l erreur L2
    EE_L2 = (UU_exact-UU)'*MM*(UU_exact-UU)
    log(sqrt(EE_L2/(UU_exact'*MM*UU_exact)));
    % Calcul de l erreur H1
    EE_H1 = (UU_exact-UU)'*(KKb)*(UU_exact-UU)
    log(sqrt(EE_H1/(UU_exact'*(KKb)*UU_exact)));
    % attention de bien changer le terme source (dans FF)
end
if 0
    h = [0.2;0.1;0.05;0.025];
    err_L2 = [-3.0991;-4.4263;-5.8235;-7.2133];
    err_H1 = [-3.0204;-4.2248;-5.4403;-6.6234];
    h = log(1./h);

    figure()
    plot(h,err_L2,h,err_H1)
    legend('Norme L^2','Norme H^1')
    xlabel({'$\log(1/h)$'},'Interpreter','latex')
    ylabel({'$\log(\Vert u-u_h\Vert/\Vert u\Vert)$'},'Interpreter','latex')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                        fin de la routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

