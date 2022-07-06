Refinement = 6;
[xx,yy] = getMesh (Refinement);
E = 0.03;
Disp = 0.3;
Power = 2;
LengthScale = 0.1;
Cov = CovarianceKernel (xx,yy,Refinement,LengthScale,E,Disp,Power);
[U,S,Vt] = svd(Cov);

