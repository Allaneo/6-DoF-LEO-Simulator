function B_ecef_T = MT_NED2ECEF(B_ned_nT, phi, lambda)
% B_ned_nT: [Bn; Be; Bd] en nT
% phi, lambda: geodésicos [deg]
% salida: B_ecef en Tesla

sphi = sind(phi); cphi = cosd(phi);
sl   = sind(lambda); cl = cosd(lambda);

C_ecef_ned = [ ...
    -sphi*cl,  -sl,   -cphi*cl;
    -sphi*sl,   cl,   -cphi*sl;
     cphi,      0,    -sphi ];

B_ecef_T = C_ecef_ned * (B_ned_nT(:) * 1e-9);  % nT -> T
end
