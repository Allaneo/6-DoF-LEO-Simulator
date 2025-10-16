%% Modelo de estados para un cuerpo rigido con x_lineal = [r(3x1);v(3x1)] y a_angular = [q(4x1);omega(3x1)]
function [xdot_lineal, xdot_angular] = MD_modelo_de_estados(x_lineal,F_tot,x_angular,M_tot,m,I)
%% Parametros 

%% Calculo lineal
A_lineal = [zeros(3),eye(3);
    zeros(3,6)];

xdot_lineal = A_lineal*x_lineal+[zeros(3,1);F_tot]/m;

%% Calculo angular
q=x_angular(1:4);
q = q/norm(q);          %normalizacion de q

wx=x_angular(5);
wy=x_angular(6);
wz=x_angular(7);
omega = [wx;wy;wz];

Omega = [ 0, -wx, -wy, -wz;
          wx,  0,  wz, -wy;
          wy, -wz, 0,   wx;
          wz,  wy, -wx,  0];

qdot = 0.5 * Omega * q;

omegadot = I\( M_tot - cross(omega, I*omega));

xdot_angular = [qdot;omegadot];
end