function  r_inertial= MT_BODY2INERTIAL(q,r_inertial)
% Inputs:
% q: cuaternión [4x1]
% vector a transformar r_body de [3x1]

% salida: vector r_inertial en terna inertial [3x1]

q0=q(1); q1=q(2); q2=q(3); q3=q(4);
C_IB = [  q0^2+q1^2 - q2^2 - q3^2,   2*(q1*q2 - q0*q3),         2*(q1*q3 + q0*q2);
          2*(q1*q2 + q0*q3),         q0^2 - q1^2 + q2^2 - q3^2, 2*(q2*q3 - q0*q1);
          2*(q1*q3 - q0*q2),         2*(q2*q3 + q0*q1),         q0^2 - q1^2 - q2^2 + q3^2 ];

% direction in intertial frame
r_inertial = C_IB * r_inertial;
end
