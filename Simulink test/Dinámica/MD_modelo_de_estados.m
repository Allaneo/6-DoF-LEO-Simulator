function [xdot_lineal, xdot_angular] = MD_modelo_de_estados(x_lineal,F_tot,x_angular,M_tot,m,I)
% Modelo de estados para un cuerpo rigido con x_lineal = [r(3x1);v(3x1)] y x_angular = [q(4x1);omega(3x1)]

% Implementa el modelo dinámico 6-DOF del satélite separando la dinámica 
% traslacional y rotacional: integra posición y velocidad a partir de la 
% fuerza total, e integra la actitud (cuaternión) y la velocidad angular a
% partir del torque total y la inercia. Es el bloque central de evolución 
% de estados dentro del simulador, alimentado por los módulos de fuerzas y 
% momentos perturbadores.

% Entradas:
% x_lineal: estado traslacional [r; v] (6x1)
% F_tot: fuerza total en el mismo marco que x_lineal (3x1) [N]
% x_angular: estado rotacional [q; omega] (7x1), cuaternión scalar-first
% M_tot: torque total en Body (3x1) [N*m]
% m: masa [kg]
% I: matriz de inercia en Body (3x3) [kg*m^2]

% Salidas:
% xdot_lineal: derivada del estado traslacional (6x1)
% xdot_angular: derivada del estado rotacional (7x1)
%% Dinámica lineal
A_lineal = [zeros(3),eye(3);
    zeros(3,6)];

xdot_lineal = A_lineal*x_lineal+[zeros(3,1);F_tot]/m;

%% Dinámica angular
q  = x_angular(1:4);
wx = x_angular(5);
wy = x_angular(6);
wz = x_angular(7);
omega = [wx;wy;wz];

Omega = [ 0, -wx, -wy, -wz;
          wx,  0,  wz, -wy;
          wy, -wz, 0,   wx;
          wz,  wy, -wx,  0];

qdot = 0.5 * Omega * q;

omegadot = I\( M_tot - cross(omega, I*omega));

xdot_angular = [qdot;omegadot];
end