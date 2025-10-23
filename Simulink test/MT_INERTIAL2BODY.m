function  r_body= MT_INERTIAL2BODY(q,r_inertial)
% Inputs:
% q: cuaternión [4x1]
% vector a transformar r_inertial de [3x1]

% salida: vector r_body en terna body [3x1]

q0=q(1); q1=q(2); q2=q(3); q3=q(4);
C_IB = [  q0^2+q1^2 - q2^2 - q3^2,   2*(q1*q2 - q0*q3),         2*(q1*q3 + q0*q2);
          2*(q1*q2 + q0*q3),         q0^2 - q1^2 + q2^2 - q3^2, 2*(q2*q3 - q0*q1);
          2*(q1*q3 - q0*q2),         2*(q2*q3 + q0*q1),         q0^2 - q1^2 - q2^2 + q3^2 ];

% inertial -> body
C_BI = C_IB.';

% direction in body frame
r_body = C_BI * r_inertial;
end
