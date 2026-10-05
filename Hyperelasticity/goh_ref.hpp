// Hand-written reference kernel (same signature as the generated one)
#include <cmath>
inline void goh_ref(const double* F, double p, double c, double k1, double k2,
                    double k, double gam, double* psi, double* P, double* A)
{
  const double a[2][3] = {{0., cos(gam), sin(gam)}, {0., cos(gam), -sin(gam)}};  // fibres in the y-z plane
  double H[2][3][3];
  for (int s=0;s<2;s++) for (int i=0;i<3;i++) for (int j=0;j<3;j++)
    H[s][i][j] = k*(i==j) + (1-3*k)*a[s][i]*a[s][j];
  double C[3][3]; double I1=0;
  for (int i=0;i<3;i++) for (int j=0;j<3;j++){ C[i][j]=0; for(int m=0;m<3;m++) C[i][j]+=F[3*m+i]*F[3*m+j]; }
  I1=C[0][0]+C[1][1]+C[2][2];
  double cof[3][3];
  for (int i=0;i<3;i++) for (int j=0;j<3;j++){
    int i1=(i+1)%3,i2=(i+2)%3,j1=(j+1)%3,j2=(j+2)%3;
    cof[i][j]=F[3*i1+j1]*F[3*i2+j2]-F[3*i1+j2]*F[3*i2+j1]; }
  double J=0; for (int j=0;j<3;j++) J+=F[j]*cof[0][j];
  double E[2],e[2],FH[2][3][3];
  for (int s=0;s<2;s++){ double t=0; for(int i=0;i<3;i++)for(int j=0;j<3;j++) t+=C[i][j]*H[s][i][j];
    E[s]=t-1; e[s]=exp(k2*E[s]*E[s]);
    for(int i=0;i<3;i++)for(int j=0;j<3;j++){ FH[s][i][j]=0; for(int m=0;m<3;m++) FH[s][i][j]+=F[3*i+m]*H[s][m][j]; } }
  *psi = 0.5*c*(I1-3)+p*(J-1)+k1/(2*k2)*((e[0]-1)+(e[1]-1));
  for (int i=0;i<3;i++) for (int j=0;j<3;j++){
    double v=c*F[3*i+j]+p*cof[i][j];
    for(int s=0;s<2;s++) v+=2*k1*E[s]*e[s]*FH[s][i][j];
    P[3*i+j]=v; }
  const int eps[3][3][3]={{{0,0,0},{0,0,1},{0,-1,0}},{{0,0,-1},{0,0,0},{1,0,0}},{{0,1,0},{-1,0,0},{0,0,0}}};
  for(int i=0;i<3;i++)for(int j=0;j<3;j++)for(int kk=0;kk<3;kk++)for(int l=0;l<3;l++){
    double v=c*(i==kk)*(j==l);
    for(int n=0;n<3;n++)for(int m=0;m<3;m++) v+=p*eps[i][kk][n]*eps[j][l][m]*F[3*n+m];
    for(int s=0;s<2;s++){
      double d1=k1*E[s]*e[s], d2=k1*(1+2*k2*E[s]*E[s])*e[s];
      v+=4*d2*FH[s][i][j]*FH[s][kk][l]+2*d1*(i==kk)*H[s][j][l]; }
    A[9*(3*i+j)+3*kk+l]=v; }
}
