 //Question 2
 #include<stdio.h>
 #include<stdlib.h>
 #include<math.h>
 #define pi 3.14159265 
 int main(){
 int N=101,i,j,k;	
 double length=N,L=2.0,W=2.0,dx=(2.0/(length-1)),dy=(2.0/(length-1)); 
 double Tanalytical[N][N]={0},x[N*N]={0},b[N*N]={0};
 double Ap[N*N]={0},An[N*N]={0},As[N*N]={0},Ae[N*N]={0},Aw[N*N]={0};
 double Lw[N*N]={0},Ls[N*N]={0},Lp[N*N]={0},Un[N*N]={0},Ue[N*N]={0},Mnw[N*N]={0},Mse[N*N]={0};
 double phi_new[N*N]={0},phi_old[N*N]={0},Temp[N][N]={0},Tm=2.0;
 double x1[N*N]={0}, x2[N*N]={0}, X[N*N]={0};
 double sum=0,Y[N*N]={0},error=1.0,alpha=0;  

 //Analytical Solution

  for(i=0;i<N;i++){
     Tanalytical[N-1][i]=1;
 }
 for(j=1;j<N-1;j++){
     for(i=1;i<N-1;i++){ 
        Tanalytical[j][i]=Tm*((sinh((pi*(dy*j))/L))/(sinh((pi*W)/L)) *(sin((pi*(i*dx))/L)) );     
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
    k=0;
   for(i=0;i<((N*N));i++){
         if(i%N==(N-1)){
	
    b[i]=Tm*sin((pi*(k*dx))/L);
    k++;
   }
 }


 for(i=0;i<((N*N));i++){
      Lw[i]=Aw[i]/(1+alpha*Un[i-N]);
      Ls[i]=As[i]/(1+alpha*Ue[i-1]);
      if(i==0){
        Lp[i]=Ap[i];
  }
    else if((i>=1)&&(i<N)){
        Lp[i]=Ap[i]-(Ls[i]*Un[i-1]) +(alpha*Ls[i]*Ue[i-1]);
     }
    else{
       Lp[i]=Ap[i]-(Lw[i]*Ue[i-N])-(Ls[i]*Un[i-1]) +(alpha*Ls[i]*Ue[i-1]) +(alpha*Lw[i]*Un[i-N]);
   }
       Un[i]=(An[i] -alpha*Lw[i]*Un[i-N])/Lp[i];
       Ue[i]=(Ae[i] -alpha*Ls[i]*Ue[i-1])/Lp[i];
  }
 for(i=0;i<((N*N));i++){
   if(i>=N){
     Mnw[i]=Lw[i]*Un[i-N];
  }
  if(i>=1){
    Mse[i]=Ls[i]*Ue[i-1];
  }
 }
 while(error>1e-6){ 
  error=0.0;
  for(i=0;i<(N*N)-N+1;i++){	
      x1[i]=Mse[i]*phi_old[i+N-1] -Mse[i]*alpha*(phi_old[i-1]+phi_old[i+N]-phi_old[i]);
   }
  for(i=N-1;i<N*N;i++){
    x2[i]=Mnw[i]*phi_old[i-N+1] -Mnw[i]*alpha*(phi_old[i-N] +phi_old[i+1]-phi_old[i]); 
   }
 for(i=0;i<N*N;i++){
   X[i]=b[i]+x1[i]+x2[i];
 } 

//Forward Substitution LY=X

 for(i=0;i<N*N;i++){
   if((i>=1)&&(i<N)){
     sum=(Ls[i]*Y[i-1]);
   }
  else if(i>=N){
     sum=(Ls[i]*Y[i-1])+(Lw[i]*Y[i-N]);
   }

    Y[i]=(X[i]-sum)/Lp[i];
    sum=0;
  }

//Backward substitution; U phi=Y

  sum=0.0;
  for(i=(N*N)-1;i>=0;i--){ 
     if((i<(N*N)-1)&&(i>(N*N)-N)){
        sum=(Un[i]*phi_new[i+1]);
   }
    else if(i<=((N*N)-N)){
       sum=(Un[i]*phi_new[i+1])+(Ue[i]*phi_new[i+N]);
   }

    phi_new[i]=(Y[i]-sum)/1;
    sum=0;
  }

  for(i=0;i<N*N;i++){
      error=error+fabs(phi_new[i]-phi_old[i]);
   }
 
  for(i=0;i<N*N;i++){
     phi_old[i]=phi_new[i];
   }
  } 

  k=0;
  for(i=0;i<N;i++){
    for(j=0;j<N;j++){
        Temp[j][i]=phi_new[k];
        k++;
    }
  }

//Results 
 FILE *fp;
 fp=fopen("SIP Temperature.dat","w");
 for(i=0;i<N;i++){
   for(j=0;j<N;j++){
    fprintf(fp,"%lf\t%lf\t%lf\t%lf\n",(j*dx),(i*dy),Temp[i][j],Tanalytical[i][j]);
  }
  }
 fclose(fp);
 fp=fopen("x_SIP Mid line temperatures.dat","w");
  for(j=0;j<N;j++){
    fprintf(fp,"%lf\t%lf\t%lf\n",(j*dx),Temp[50][j],Tanalytical[50][j]);
  }
  fclose(fp);
  fp=fopen("y_SIP Mid line temperatures.dat","w");
     for(j=0;j<N;j++){
       fprintf(fp,"%lf\t%lf\t%lf\n",(j*dx),Temp[j][50],Tanalytical[j][50]);
 }
 fclose(fp);
 return 0;
}
