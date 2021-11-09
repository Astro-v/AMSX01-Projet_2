function val = mat_A(x,y)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mat_A :
% Evaluation de la matrice A .
%
% SYNOPSIS val = mat_A(x,y)
%          
% INPUT * x,y : les 2 coordonnees du point ou on veut evaluer la fonction.
%
% OUTPUT - val: valeur de la matrice sur ce point.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global ast;
if (ast=='non')
    %val = 1;

    %val = [1 0;0 1];

    %val = [1 0;0 2];

    %val = [2+sin(2*pi*x),0;0,4];

    val = [2+sin(2*pi*x),0;0,4+sin(2*pi*x)];

    %val = (2+sin(2*pi*x))*(4+sin(2*pi*y))*[1,0;0,1];
else
    global A_ast;
    val = A_ast;
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                     fin de la fonction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
