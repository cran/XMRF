#include<math.h>
#include<string.h>
#include<R.h>
#include<stdlib.h>

int sign(double x){
  int ans = 1;
  if(x<0) ans= -1;
  if(x==0) ans = 0;
  return ans;
}

double absolute(double x){
  double ans = x;
  if(x<0) ans= -1.0*x;
  return ans;
}

//void foo2(int *nin, double *x, double *y, double *xa, double *ya)

//void LPGM_neighborhood(double *X, double *Y, int *startb, int *intercept, int *lambda, int *nin, int *pin, int *nlamsin, double *alphas, double *Bmat)
void LPGM_neighborhood(double *Xt, double *Y, double *startb, double *lambda, int *nin, int *pin, int *nlamsin, double *alphas, double *Bmat)
{
  unsigned int n = nin[0];
  unsigned int p = pin[0];
  unsigned int nlams = nlamsin[0];
  unsigned int goon = 1;
  double thr = 0.001;
  double maxit = 1000;
  double *Bhat, *tt1, *smthgrad, *oldb, *newb, *tempc, *tempco; 
  int i, j, l, iter;
  double ind, temp, tempo, temp2, tempo2, temp3, diff, tk, bs, smthobj_newb, smthobj_oldb;
  //double ninv;
  
  Bhat = (double*)malloc(p*sizeof(double));
  tt1 = (double*)malloc(p*sizeof(double));
  smthgrad = (double*)malloc(p*sizeof(double));
  oldb = (double*)malloc(p*sizeof(double));
  newb = (double*)malloc(p*sizeof(double));
  tempc = (double*)malloc(n*sizeof(double));
  tempco = (double*)malloc(n*sizeof(double));
  
  //ninv = 1.0/(double)n;
  
  //printf("start!!!\n\n");
  
  if(Bhat == NULL || tt1 ==  NULL || smthgrad == NULL || oldb == NULL || newb == NULL || tempc == NULL){
	//printf("Error. Out of memory. \n");
	//exit(0);
	goon = 0;
	Bmat = NULL;
	alphas = NULL;
  }
  
  //printf("start!!!\n\n");
  
  //printf("Here 1\n");
  
  if(goon == 1){
  	  
	  // Initialize the Bhat
	  temp = 0;
	  for(i=0; i<p; i++) {
		Bhat[i]=0;
		temp += startb[i];
	  }
	  
	  if(temp != 0){memcpy(Bhat, startb, p*sizeof(double));}
	  //if(temp != 0){for(j=0; j<p; j++){Bhat[j] = startb[j];}}
	  
	  //printf("Here 2\n");

	  
	  // Do : tt1 <- t(Y)%*%Xt
	  for(j=0; j<p; j++){
		temp = 0;
		for(i=0; i<n; i++){
		  // equivalent to Y[i]*Xt[i][j]
		  temp += (Y[i]*Xt[j+p*i]);
		}
		tt1[j] = temp;
	  }
	  
	  // check
	//  printf("tt1\n");
	//  for(j=0; j<p; j++) printf("%.2f\t", tt1[j]);
	  
	  //printf("Here 3\n");  


	  // Loop through each lambda
	  for(l=0; l<nlams; l++){
	  //check
	//	printf("\n-------------------------------------------------\n");
	//	printf("\n-------------------------------------------------\n");
	//	printf("\nLambda: %f\n", lambda[l]);
		iter =1;
		ind = 1.0;
		
		while(thr<ind && iter<maxit){
	//	  printf("\n-------------------------------------------------\n");
	//	  printf("Iter : %d\n", iter);
		
		  //printf("Here 4\n");
		  memcpy(oldb, Bhat, p*sizeof(double));
		  //for(j=0;j<p;j++) {oldb[j] = Bhat[j];}
		  tk = 1.0;
		  bs = 0.5;
		  
		  // Compute: smthgrad = -(t(Xt)%*%(Y-exp(Xt%*%Bhat)))/n
		  // 1. get exp(Xt %*% Bhat)
		  for(i=0; i<n; i++){
			tempc[i] = 0;
			//for(j=0; j<p; j++){tempc[i] += Xt[i+n*j]*Bhat[j];}
			for(j=0; j<p; j++){tempc[i] += (Xt[j+p*i]*Bhat[j]);}
			tempc[i] = exp(tempc[i]);
		  }
		  // 2. get all
		  for(j=0; j<p; j++){
			smthgrad[j] = 0;
			for(i=0; i<n; i++){smthgrad[j] += (-1.0*Xt[j+p*i]*(Y[i]-tempc[i]));}
			//smthgrad[j] = smthgrad[j]/n;
			smthgrad[j] = smthgrad[j]/(double)n;
			//smthgrad[j] = smthgrad[j]*ninv;
		  }
		  
		  // check
	//	  printf("\t\tsmthgrad:\t");
	//	  for(j=0; j<p; j++) printf("%10.2f", smthgrad[j]);
		  
		  // check
	//	  printf("\n\t\toldb:\t");
	//	  for(j=0; j<p; j++) printf("%10.2f", oldb[j]);

		  // Compute: newb
		  //    tmp = oldb-tk*smthgrad
		  //    newb = matrix(sign(tmp)*sapply(abs(tmp) - tk*lambda[i],max,0),p+1,1);
		  //    newb[1] = tmp[1];
		  newb[0] = oldb[0] - tk*smthgrad[0];
		  for(j=1; j<p; j++){
			temp = oldb[j] - tk*smthgrad[j];
			//temp2 = abs(temp) - tk*(double)lambda[l];
			temp2 = absolute(temp) - tk*lambda[l];
			//printf("temp: %10.2f, abs.temp: %10.2f, lambda: %5.2f, temp2 : %10.2f, ", temp, abs(temp), (double)lambda[l], temp2);
			//printf("temp: %10.2f, temp2 : %10.2f\n", temp, temp2);
			newb[j] = sign(temp)* temp2;
			if(temp2 < 0) {newb[j] = 0;}
		  } 
		  
		  // check
	//	  printf("\n\t\tnewb:\t");
	//	  for(j=0; j<p; j++) printf("%10.2f", newb[j]);
		  
		  // Get smthobj_newb and smthobj_oldb, compute:
		  // 		smthobj_newb = -(1/n)*(tt1%*%newb-sum(exp(Xt%*%newb)))
		  //		smthobj_oldb = -(1/n)*(tt1%*%oldb-sum(exp(Xt%*%oldb)))
		  // 1. get sum(exp(Xt%*%newb)) and sum(exp(Xt%*%oldb))
		  temp = 0;
		  tempo = 0;
		  for(i=0; i<n; i++){
			tempc[i] = 0;
			tempco[i] = 0;
			for(j=0; j<p; j++){
				tempc[i] += (Xt[j+p*i]*newb[j]);
				tempco[i] += (Xt[j+p*i]*oldb[j]);
			}
			tempc[i] = exp(tempc[i]);
			tempco[i] = exp(tempco[i]);
			
			temp += tempc[i];
			tempo += tempco[i];
		  }
		  
		  // 2. get (tt1%*%newb) and (tt1%*%oldb)
		  temp2 = 0;
		  tempo2 = 0;
		  for(j=0; j<p; j++){
			temp2 += (tt1[j]*newb[j]);
			tempo2 += (tt1[j]*oldb[j]);
		  }
		  
		  // 3. get smthobj_newb and get smthobj_oldb
		  smthobj_newb = (temp -temp2)/(double)n;
		  smthobj_oldb = (tempo -tempo2)/(double)n;
		  //smthobj_newb = (temp -temp2)*ninv;
		  //smthobj_oldb = (tempo -tempo2)*ninv;
		  
		  
		  //check
		  //printf("\n\tsmthobj_newb : %g \t smthobj_oldb: %g\n", smthobj_newb, smthobj_oldb);
		  
		  // Begin line search
		  // Loop if smthobj_newb > (smthobj_oldb -t(smthgrad)%*%(oldb-newb)+(1/(2*tk))*sum((oldb-newb)^2)) 
		  diff = 0;
		  temp = 0;
		  temp2 = 0;
		  for(j=0; j<p; j++){
			temp3 = oldb[j] - newb[j];
			temp += (smthgrad[j]*temp3);
			temp2 += (temp3*temp3);
		  }
		  diff = smthobj_oldb - temp + (temp2/(2.0*tk));

		  //printf("here 5!\n");
	//	  printf("\n\tsmthobj_newb: %10g \t smthobj_oldb: %10g \t diff: %10g\n", smthobj_newb, smthobj_oldb, diff);
		  
		  while(smthobj_newb > diff || isinf(smthobj_newb)){
	//      while(smthobj_newb > diff){
			tk = bs * tk;
			//    tmp = oldb-tk*smthgrad
			//    newb = matrix(sign(tmp)*sapply(abs(tmp) - tk*lambda[i],max,0),p+1,1);
			//    newb[1] = tmp[1];
			newb[0] = oldb[0] - tk*smthgrad[0];
			for(j=1; j<p; j++){
			  temp = oldb[j] - tk*smthgrad[j];
			  //temp2 = abs(temp) - tk*lambda[l];
			  temp2 = absolute(temp) - tk*lambda[l];
			  newb[j] = sign(temp)* temp2;
			  if(temp2 < 0) {newb[j] = 0;}
			}
		  
			//Compute: smthobj_newb = -(1/n)*(tt1%*%newb-sum(exp(Xt%*%newb)))
			// 1. get sum(exp(Xt%*%newb))
			temp = 0;
			for(i=0; i<n; i++){
			  tempc[i] = 0;
			  //for(j=0; j<p; j++){tempc[i] += Xt[i+n*j]*newb[j];}
			  for(j=0; j<p; j++){tempc[i] += (Xt[j+p*i]*newb[j]);}
			  tempc[i] = exp(tempc[i]);
			  temp += tempc[i];
			}
		  
	/*        // 2. get (tt1%*%newb) 
			temp2 = 0;
			for(j=0; j<p; j++){
			  temp2 += (tt1[j]*newb[j]);
			}
		
			// 3. get smthobj_newb
			smthobj_newb = (temp -temp2)/(double)n;
			
			// Get new diff
			diff = 0;
			temp = 0;
			temp2 = 0;
			for(j=0; j<p; j++){
				temp3 = oldb[j] - newb[j];
				temp += (smthgrad[j]*temp3);
				temp2 += (temp3*temp3);
			}
			diff = smthobj_oldb - temp + (temp2/(2.0*tk));
	*/		
			// 2. get (tt1%*%newb) and part that compute for diff
			temp2 = 0;
			tempo = 0;
			tempo2 = 0;
			
			for(j=0; j<p; j++){
			  temp2 += (tt1[j]*newb[j]);
			  
			  temp3 = oldb[j] - newb[j];
			  tempo += (smthgrad[j]*temp3);
			  tempo2 += (temp3*temp3);
			}
			
			// 3. get smthobj_newb and diff
			smthobj_newb = (temp -temp2)/(double)n;
			//smthobj_newb = (temp -temp2)*ninv;
			diff = smthobj_oldb - tempo + (tempo2/(2.0*tk));
			

			//printf("smthobj_newb: %g\t diff: %g\n", smthobj_newb, diff);
	//		printf("\n\tsmthobj_newb: %10g \t smthobj_oldb: %10g \t diff: %10g\n", smthobj_newb, smthobj_oldb, diff);
			
		  } // end of line search loop
	 
		  memcpy(Bhat, newb, p*sizeof(double));
		  //for(j=0; j<p; j++){ Bhat[j] = newb[j];}
		  iter += 1;
		  
		  //check
	//	  printf("Bhat: \t");
	//	  for(j=0; j<p; j++){ printf("%10g ", Bhat[j]);}
	//	  printf("\noldb: \t");
	//	  for(j=0; j<p; j++){ printf("%10g", oldb[j]);}
		  
		  // Compute: ind = vec.norm(oldb - Bhat)/vec.norm(oldb)
		  //  vec.norm sqrt(t(A)%*%A)
		  temp = 0.0;
		  temp2 = 0.0;
		  for(j=0; j<p; j++){
			temp += ((oldb[j]-Bhat[j])*(oldb[j]-Bhat[j]));
			temp2 += (oldb[j]*oldb[j]);
		  }
		  ind = sqrt(temp)/sqrt(temp2);
		  
		  if(temp2 == 0){ind = 1;}
	//	  printf("ind: %g\n", ind);
		  
		  // Check fo stupid values!!!!
		  //  if (is.na(ind)) { ind = 0}
		  //  if (!is.finite(ind)) {ind = 1}
		  
		  //--
		} // end while for iter loop
		
		// Store the final residual & estimated Bhat
		alphas[l] = Bhat[0];
		//printf("p: %d\n", p);
		for(j=1; j<p; j++){
		  Bmat[l+nlams*(j-1)] = Bhat[j];
		}
		//--
	  } // end nlams loop
	  

  }
  
  if(Bhat != NULL) {free(Bhat);}
  if(tt1 != NULL) {free(tt1);}
  if(smthgrad != NULL) {free(smthgrad);}
  if(oldb != NULL) {free(oldb);}
  if(newb != NULL) {free(newb);}
  if(tempc != NULL) {free(tempc);}
  if(tempco != NULL) {free(tempco);}
}
