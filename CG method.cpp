//Question 3
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#define pi 3.14159265 
int main(){
int N=101,i,j,k;	
double length=N,dx=(1/(length-1)),dy=(1/(length-1));
double Tanalytical[N][N]={0}, b[N*N]={0};
double Ap[N*N]={0},An[N*N]={0},As[N*N]={0},Ae[N*N]={0},Aw[N*N]={0}; 
double phi_new[N*N]={0},phi_old[N*N]={0},Temp[N][N]={0}; 
double sum=0,Y[N*N]={0}; 

//Analytical Solution

for(i=0;i<N;i++){
Tanalytical[N-1][i]=1;
}
for(j=1;j<N-1;j++){
for(i=1;i<N-1;i++){ 
    Tanalytical[j][i]= sin(2*pi*(i*dx)) *sin(2*pi*(j*dy));
}
}
 
 //Computational solution

for(i=0;i<N*N;i++){
Ap[i]=1;
}
for(i=N;i<((N*N)-N);i++){
if(((i%N)!=0)&&((i%N)!=(N-1))){
Ap[i]=-4;
An[i]=1;
As[i]=1;
Ae[i]=1;
Aw[i]=1;
}
}
  
 int l;l=0;
 for(i=0;i<(N);i++){
	for(j=0;j<(N);j++){
		l=(i)*(N)+j;
		if((i==0)||(i==N-1)||(j==0)||(j==N-1))
			b[l]=0.0;
		else	
		 b[l]= (-8)*(pi*pi)*sin(2*pi*(i*dx))*sin(2*pi*(j*dy))*dx*dx;
	
	//	printf("b[%d]  %lf\n",l,b[l]);
	}
}
 
 
 double R[N*N]={0},P[N*N]={0},R_old[N*N]={0},AP[N*N]={0},X[N*N]={0}; 
 double RTR=0,RoTRo=0,PTAP=0,err;
 double alpha,beta;
k=0;
err=1.0;
for(i=0;i<N*N;i++){
	R[i]=b[i];
	P[i]=R[i];
}
while(fabs(err>1e-6)){
    RTR=0;
	for(i=0;i<N*N;i++){
	RTR=RTR+R[i]*R[i];
    }
   for(i=0;i<N*N;i++){
   	if(i==0)
   	   AP[i]=Ap[i]*P[i] +An[i]*P[i+1] +Ae[i]*P[i+N];
   	else if(i>=1&&i<N)
	     AP[i]=Ap[i]*P[i] +An[i]*P[i+1] +Ae[i]*P[i+N]  +As[i]*P[i-1];
	else     
	 AP[i]=Ap[i]*P[i] +An[i]*P[i+1] +Ae[i]*P[i+N] +As[i]*P[i-1] +Aw[i]*P[i-N];
   }
   PTAP=0;
   for(i=0;i<N*N;i++){
   PTAP=PTAP+P[i]*AP[i];
  }
  alpha=RTR/PTAP;
  for(i=0;i<N*N;i++){
   X[i]=X[i] +alpha*P[i];
}
for(i=0;i<N*N;i++){
    R_old[i]=R[i];
}
for(i=0;i<N*N;i++){
    R[i]=R[i]-alpha*AP[i];
}
 RTR=0;
	for(i=0;i<N*N;i++){
	RTR=RTR+R[i]*R[i];
    }
RoTRo=0;    
for(i=0;i<N*N;i++){
    RoTRo=RoTRo+R_old[i]*R_old[i];
}
beta=RTR/RoTRo;
for(i=0;i<N*N;i++){
    P[i]=R[i]+beta*P[i];
}
	k++;
	err=RTR;
}
 
k=0;
for(i=0;i<N;i++){
   for(j=0;j<N;j++){
       Temp[j][i]=X[k];
       k++;
   }
}
//Results 
FILE *fp;
fp=fopen("CGTemperature.dat","w");
for(i=0;i<N;i++){
for(j=0;j<N;j++){
fprintf(fp,"%lf\t%lf\t%lf\t%lf\n",(j*dx),(i*dy),Temp[i][j],Tanalytical[i][j]);
}
}
fclose(fp);
fp=fopen("X_CGMid line temperatures.dat","w");
for(j=0;j<N;j++){
fprintf(fp,"%lf\t%lf\t%lf\n",(j*dx),Temp[50][j],Tanalytical[50][j]);
}
fclose(fp);
fp=fopen("Y_CGMid line temperatures.dat","w");
for(j=0;j<N;j++){
fprintf(fp,"%lf\t%lf\t%lf\n",(j*dy),Temp[j][50],Tanalytical[j][50]);
}
fclose(fp);
return 0;
}
