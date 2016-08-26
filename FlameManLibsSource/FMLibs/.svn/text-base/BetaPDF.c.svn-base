#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "BetaPDF.h"

void BetaPDF(int nx, double *x, double mean, double var, double *pdf,
                     int *pdfBound)
{
  int j, jMean, j0, j1;
  double dx, sigma;
  double alpha, beta, tmp;
  double sumPDF, meanPDF;
  double truncate = 1.0e-4;
  
  /* ensure that the mean and variance are properly bounded */
  mean = fmin(1.0, fmax(0.0, mean));
  var  = fmin(fmax(0.0, var), mean*(1.0 - mean));
  
  sigma = sqrt(var);    /*  compute the standard deviation */
  
  /*  compute the shape parameters for the beta distribution */
  tmp   = mean*(1-mean)/var - 1;
  alpha = mean*tmp;
  beta  = (1-mean)*tmp;
  if (fabs(alpha - 1) < 1.0e-6) alpha = 1.0;
  if (fabs(beta - 1)  < 1.0e-6) beta = 1.0;
  
/*   printf("mean  = %15.6e, var  = %15.6e\n", mean, var); */
/*   printf("alpha = %15.6e, beta = %15.6e\n", alpha, beta); */

  /*  find the location of the mesh point before the mean */
  jMean = -1;
  do { jMean++; } while (x[jMean] < mean);
  jMean--;
  
  /*  set the pdf to zero */
  for (j=0; j<nx; ++j)
    pdf[j] = 0.0;
  
  /*  if the mesh size about the mean is less than 1.5 standard deviations */
  /*  compute at two points about the mean such that the mean is conserved */
  if (x[jMean+1] - x[jMean] > 1.5*sigma)
    {
      j0 = jMean;
      j1 = jMean + 1;
      
      /*  compute pdf at points about the mean */
      dx = (x[j1]-x[j0]);
      pdf[j0] = fabs(x[j1] - mean)/dx;
      pdf[j1] = fabs(mean - x[j0])/dx;
      
      /*  normalise the PDF */
      sumPDF  = pdf[j0]+pdf[j1];
      pdf[j0] = pdf[j0]/sumPDF;
      pdf[j1] = pdf[j1]/sumPDF;
      
      /*  store bounds of non-zero PDF values */
      if (pdfBound != NULL) {
        pdfBound[0] = j0;
        pdfBound[1] = j1;
      }
      return;
    }
  
  /* ----- compute a beta distribution based on alpha and beta ----- */
  sumPDF  = 0.0;
  if (alpha == 1.0 && beta == 1.0) /*  uniform distribution (not likely) */
    {
/* HP      fprintf(stderr, "Beta distribution is uniform\n"); */
      j0 = 0;
      j1 = nx-1;
      pdf[j0] = pdf[j1] = 0.5;
      sumPDF += 1.0;
      for (j=j0+1; j < j1; j++)
        {
          dx = 0.5*x[j+1]-0.5*x[j-1];
          pdf[j] = 1./(x[j1] - x[j0]);
          sumPDF += pdf[j];
        }
    }
  /*  u shaped distribution */
  else if ( (alpha < 1.0 && beta < 1.0) )
    {
/* HP      fprintf(stderr, "Beta distribution is uniform or u-shaped\n"); */
      /*  need to calculate over whole domain, ignore the ends */
      j0 = 0;
      j1 = nx-1;
      /*  compute the PDF over all the interior points */
      for (j=j0+1; j < j1; j++)
        {
          dx = 0.5*(x[j+1]-x[j-1]);
          pdf[j] = scaledBetaDistribution(x[j], dx, alpha, beta);
          sumPDF += pdf[j];
        }
      /*  set the boundary points to a value which preserves the mean */
      /*  in theory they are undefined (infinite) */
      pdf[j0] = pdf[j1] = 0.5*(1.0 - sumPDF);
      sumPDF = 1.0;
    }
  /*  decreasing distribution */
  else if ( (alpha < 1.0 && beta >= 1.0) || (alpha == 1.0 && beta > 1.0) )
    {
/* HP      fprintf(stderr, "Beta distribution is decreasing\n"); */
      /*  start at the left boundary and march till less than truncate */
      j0 = 0;
      /* dx = 0.5*(x[1]-x[0]); */
      /* pdf[j0] = scaledBetaDistribution(x[j0], dx, alpha, beta); */
      /* sumPDF += pdf[j0]; */
      for (j = j0+1; j<nx-1; j++)
        {
          dx = 0.5*(x[j+1]-x[j-1]);
          pdf[j] = scaledBetaDistribution(x[j], dx, alpha, beta);
          
          if (pdf[j] < truncate)
            {
              pdf[j]  = 0.0;
              pdf[j0] = 1.0 - sumPDF;
              break;
            }
          
          j1 = j;
          sumPDF += pdf[j];
          pdf[j0] = 1.0 - sumPDF;
        }
      sumPDF = 1.0;
    }
  /*  increasing distribution */
  else if ( (alpha > 1.0 && beta <= 1.0) || (alpha == 1.0 && beta < 1.0) )
    {
/* HP      fprintf(stderr, "Beta distribution is increasing\n"); */
      /*  start at the right boundary and march till less than truncate */
      j1 = nx-1;
      /* dx = 0.5*(x[j1] - x[j1-1]); */
      /* pdf[j1] = scaledBetaDistribution(x[j1], dx, alpha, beta); */
      /* sumPDF += pdf[j1]; */
      for (j = j1-1; j>0; j--)
        {
          dx = 0.5*(x[j+1]-x[j-1]);
          pdf[j] = scaledBetaDistribution(x[j], dx, alpha, beta);
          
          if (pdf[j] < truncate)
            {
              pdf[j]  = 0.0;
              pdf[j1] = 1.0 - sumPDF;
              break;
            }
          
          j0 = j;
          sumPDF += pdf[j];
          pdf[j1] = 1.0 - sumPDF;
        }
      sumPDF = 1.0;
    }
  /*  unimodal distribution */
  else
    {
/* HP      fprintf(stderr, "Beta distribution is unimodal\n"); */
      /*  start from the mean and decrease x till the pdf is negligible */
      for (j=jMean; j > 0; j--)
        {
          dx = 0.5*(x[j+1]-x[j-1]);
          pdf[j] = scaledBetaDistribution(x[j], dx, alpha, beta);
          
         if (pdf[j] < truncate)
            {
              pdf[j] = 0.0;
              break;
            }
          
          j0 = j;
          sumPDF += pdf[j];
        }
      
      /*  start from the mean and decrease x till the pdf is negligible */
      for (j=jMean+1; j < nx - 1; j++)
        {
          dx = 0.5*(x[j+1]-x[j-1]);
          pdf[j] = scaledBetaDistribution(x[j], dx, alpha, beta);
          
          if (pdf[j] < truncate)
            {
              pdf[j] = 0.0;
              break;
            }
          
          j1 = j;
          sumPDF += pdf[j];
        }
      
      if (j0 == 1)
        {
          pdf[0] = 1.0 - sumPDF;
          sumPDF = 1.0;
          j0 = 0;
        }
      
      if (j1 == nx - 2)
        {
          pdf[nx-1] = 1.0 - sumPDF;
          sumPDF = 1.0;
          j1 = nx-1;
        }
    }
  /* ----- end beta distribution computation ----- */
  
  /* ----- normalize the PDF ----- */
  meanPDF = 0.0;
  for (j = j0; j <= j1; j++)
    {
      pdf[j]   = pdf[j]/sumPDF;
      meanPDF += pdf[j]*x[j];
    }
  
/*   if (fabs(meanPDF - mean) > 1.0e-3) */
/*     { */
/*      printf("WARNING: mean of beta distribution not conserved\n"); */
/*      printf("alpha = %5.4f  beta = %5.4f  j0 = %d j1 = %d\n", alpha, beta, j0, j1 ); */
/*      printf("mean = %5.4f  meanPDF = %5.4f  error %6.4f percent\n", */
/*             mean, meanPDF, fabs(mean-meanPDF)/mean*100); */
/*      exit(2); */
/*     } */
  
  /*  store the bounding indices of the non-zero PDF */
  if (pdfBound != NULL) {
    pdfBound[0] = j0;
    pdfBound[1] = j1;
  }
}

double scaledBetaDistribution(double x, double dx, double alpha, double beta)
{
  double tmp, fgam;
  
  fgam = lnGamma(alpha+beta) - lnGamma(alpha) - lnGamma(beta);
  tmp = (alpha-1.0)*log(x) + (beta-1.0)*log(1.0-x) + fgam;

  return exp(tmp)*dx;
}

/* This function led to overflow*/
/* double scaledBetaDistribution(double x, double dx, double alpha, double beta) */
/* { */
/*   double pdf; */
/*    */
   /*  compute value of beta pdf at x[j] */ 
/*   pdf = pow(x, alpha-1)*pow(1-x, beta-1)*exp(lnGamma(alpha+beta)-lnGamma(alpha)-lnGamma(beta)); */
/*    */
/*   pdf = pdf*dx;*/          /*  scale for the mesh size */ 
/*   return fmin(pdf, 1.);*/  /*  ensure that the scaled PDF is less than 1 */ 
/* } */

double lnGamma(double xx)
{
  int j;
  double x, y, tmp, ser = 1.000000000190015;
  double cof[6] = {76.18009172947146, -86.50532032941677, 24.01409824083091, 
  -1.231739572450155, 0.1228650973866179E-2, -.53395239384953E-5};
  
  y = x = xx;
  tmp = x + 5.5;
  tmp -= (x+0.5)*log(tmp);
  for (j = 0; j<=5; j++) ser += cof[j]/++y;
  return -tmp + log(2.5066282746310005*ser/x);
}
