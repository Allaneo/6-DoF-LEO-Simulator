function Ja_body = MT_ECEF2BODY_matrix(q, era,Ja_ecef)
    c = cos(era); s = sin(era);
    C_IE = [c, -s, 0;
          s, c, 0; 
          0, 0, 1];
    
    q0=q(1); q1=q(2); q2=q(3); q3=q(4);
    C_IB = [  q0^2+q1^2 - q2^2 - q3^2,   2*(q1*q2 - q0*q3),         2*(q1*q3 + q0*q2);
          2*(q1*q2 + q0*q3),         q0^2 - q1^2 + q2^2 - q3^2, 2*(q2*q3 - q0*q1);
          2*(q1*q3 - q0*q2),         2*(q2*q3 + q0*q1),         q0^2 - q1^2 - q2^2 + q3^2 ];
    C_BI = C_IB';

    C_BE = C_BI * C_IE ; 
    C_EB = C_BE';

    Ja_body = C_BE * Ja_ecef * C_EB;
end