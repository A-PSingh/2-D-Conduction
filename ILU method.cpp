  //Question 1
 #include<stdio.h>
 #include<stdlib.h>
 #include<math.h>
 #define pi 3.14159265 
 int main(){
 int N=101,i,j,k;	
 double length=N,dx=(1/(length-1)),dy=(1/(length-1));
 double Tanalytical[N][N]={0},x[N*N]={0},b[N*N]={0};
 double Ap[N*N]={0},An[N*N]={0},As[N*N]={0},Ae[N*N]={0},Aw[N*N]={0};
 double Lw[N*N]={0},Ls[N*N]={0},Lp[N*N]={0},Un[N*N]={0},Ue[N*N]={0},Mnw[N*N]={0},Mse[N*N]={0};
 double phi_new[N*N]={0},phi_old[N*N]={0},Temp[N][N]={0};;
 double x1[N*N]={0}, x2[N*N]={0}, X[N*N]={0};
 double sum=0,Y[N*N]={0},error=1.0; 

 //Analytical Solution

    for(i=0;i<N;i++){
         Tanalytical[N-1][i]=1;
           }
    for(j=1;j<N-1;j++){
        for(i=1;i<N-1;i++){
           for (k=1;k<111;k++)
                    { 
                 sum=sum+(((1+pow(-1,k+1))/(k)) *((sinh(k*pi*(dy*j)))/(sinh(k*pi))) *(sin(k*pi*(i*dx))) ); 	   
                 }
                 Tanalytical[j][i]=(2/pi)*sum;
                 sum=0;
        }}
 
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
   }}
   for(i=0;i<((N*N));i++){
        if(i%N==(N-1)){
        b[i]=1;
    }}

    for(i=0;i<((N*N));i++){
          Lw[i]=Aw[i];
          Ls[i]=As[i];
       if(i==0){
         Lp[i]=Ap[i];
    }
   else if((i>=1)&&(i<N)){
          Lp[i]=Ap[i]-(Ls[i]*Un[i-1]);
   }
   else{
        Lp[i]=Ap[i]-(Lw[i]*Ue[i-N])-(Ls[i]*Un[i-1]);
      }
        Un[i]=An[i]/Lp[i];
        Ue[i]=Ae[i]/Lp[i];
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
            x1[i]=Mse[i]*phi_old[i+N-1];
          }
        for(i=N-1;i<N*N;i++){
            x2[i]=Mnw[i]*phi_old[i-N+1];
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
   fp=fopen("TemperatureILU.dat","w");
   for(i=0;i<N;i++){
       for(j=0;j<N;j++){
            fprintf(fp,"%lf\t%lf\t%lf\t%lf\n",(j*dx),(i*dy),Temp[i][j],Tanalytical[i][j]);
    }
  }
  fclose(fp);
  fp=fopen("x_Mid line temperatures.dat","w");
  for(j=0;j<N;j++){
           fprintf(fp,"%lf\t%lf\t%lf\n",(j*dx),Temp[50][j],Tanalytical[50][j]);
  }
  fclose(fp);
  fp=fopen("y_Mid line temperatures.dat","w");
  for(j=0;j<N;j++){
  fprintf(fp,"%lf\t%lf\t%lf\n",(j*dx),Temp[j][50],Tanalytical[j][50]);
  }
  fclose(fp);
  return 0;
 }
