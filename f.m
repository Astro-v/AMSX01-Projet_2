function val = f(x,y)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% f :
% Evaluation de la fonction second membre.
%
% SYNOPSIS val = f(x,y)
%          
% INPUT * x,y : les 2 coordonnees du point ou on veut evaluer la fonction.
%
% OUTPUT - val: valeur de la fonction sur ce point.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% val = 2*pi^2*sin(pi*x).*sin(pi*y); % Pour A = I
% val = 3*pi^2*sin(pi*x).*sin(pi*y); % Pour A = [1,0;0,2]
% val = pi^2*(-2*cos(2*pi*x).*cos(pi*x)+sin(2*pi*x).*sin(pi*x)+6*sin(pi*x)).*sin(pi*y); % Pour A = [2 + sin(2*pi*x),0;0,4]
val = pi^2*(-2*cos(3*pi*x)+6*sin(pi*x)).*sin(pi*y); % Pour A = [2 + sin(2*pi*x),0;0,4 + sin(2*pi*x)]
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                     fin de la fonction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
