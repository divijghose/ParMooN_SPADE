pkg load statistics
Refinement = 4;
N_R = 500;
[xx,yy] = getMesh (Refinement);
E = 0.03;
Disp = 0.3;
Power = 2;
LengthScale = 0.1;
Cov = CovarianceKernel (xx,yy,Refinement,LengthScale,E,Disp,Power);
[U,S,Vt] = svd(Cov);
dS = diag(S);
i=1;
while(cumsum(dS)(i)/sum(dS)<0.999)
  i++;
endwhile
Unew = U(:,1:i);
SNew = dS(1:i);
GaussD = stdnormal_rnd(i,N_R);
Z = GaussD.*SNew;
RealznVector = Unew*Z;

[doMean,doCoeff,doMode] = DO_Init(RealznVector, N_R);
