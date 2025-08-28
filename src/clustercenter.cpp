#include <Rcpp.h>
#include "clustercenter.h"
#include "distcomp.h"

// Here we can specify functions that compute "centers" for clusters
// Currently it is not possible for the code to use information about happy
// points (and so on) with respect to the previous center, because the caller
// will not pass this information

// [[Rcpp::export]]
void optimClusterCenterEuclid2(NumericVector clustx, NumericVector clusty, double& centerx, double& centery) {
  centerx = mean(clustx);
  centery = mean(clusty);
  return;
}

// computes the exact barycenter of a cluster for squared euclidean distances
// alorithm from Drezner 1991, Algorithm1
void exactClusterCenterEuclid2(NumericVector clustx, NumericVector clusty, double& centerx, double& centery, double penp) {
  int n = clustx.length();
  double bestcenterx = clustx[0];
  double bestcentery = clusty[0];
  double bestdist= 2*n*penp; //replace this with the cost of an empty point if you want to save that step
  
  for (int i = 0; i < (n-1); i++)
  {
    for (int j = (i+1); j < n; j++)
    {
      if (ISNA(clustx[i]) || ISNA(clustx[j]))
      {
        /* virtual points fulfill the distance condition but are useless for this */
      }else
      {
        double tempdist = (clustx[i]-clustx[j]) * (clustx[i]-clustx[j]) + (clusty[i]-clusty[j]) * (clusty[i]-clusty[j]); //same as dprime2 but without the special cases of NA points as they are already left out
        
        if (tempdist <= 8*penp) //maximum euclidean distance is 2* sqrt(2*penp); square this to get max dprime2 between points = 8*penp
        {
          //two intersection points based on vectors
          double a = sqrt(2.0*penp/tempdist -1.0/4.0);
          double x = (clustx[i] + clustx[j])/2.0 + a*(clusty[j] - clusty[i]);
          double y = (clusty[i] + clusty[j])/2.0 + a*(clustx[i] - clustx[j]);

          NumericVector res = bestPoint(x, y, clustx, clusty, penp, i, j);
          if (res[0] < bestdist)
          {
            bestdist = res[0];
            bestcenterx = res[1];
            bestcentery = res[2];
          }

          x = (clustx[i] + clustx[j])/2.0 - a*(clusty[j] - clusty[i]);
          y = (clusty[i] + clusty[j])/2.0 - a*(clustx[i] - clustx[j]);
          
          res = bestPoint(x, y, clustx, clusty, penp, i, j);
          if (res[0] < bestdist)
          {
            bestdist = res[0];
            bestcenterx = res[1];
            bestcentery = res[2];
          }

        } // both points are close enough for circles to intersect
      } // both points are "real" 
      
    }
  }
  
  centerx = bestcenterx;
  centery = bestcentery;
  
  return;
}

//Drezners Algorithm called from outside, barycenter and cost are returned
List exactClusterCenterEuclid2(NumericVector clustx, NumericVector clusty, double penp, bool aleph) {
  int n = clustx.length();
  int skip = 0; //how many times was an index skipped
  //double nclose; //number of active points NOT IMPLEMENTED YET
  ///*
  double bestcenterx = clustx[0];
  double bestcentery = clusty[0];
  double bestdist= 2*(n-1)*penp;//*/

  bool cont; //check if the new solution can be better than the current/aleph
  if (aleph)
  {
    bestcenterx = NA_REAL;
    bestcentery = NA_REAL; 
    bestdist= n*penp; //cost of an empty point
  }
 
  NumericMatrix dmat(n,n);
  dmat = cross_dmat(clustx,clusty);
  LogicalVector candidates = (dmat(_,0) > -1); //candidates(n,1)??
  
  for (int i = 0; i < (n-1); i++)
  {
    //test if enough points are close to be better than the empty bary
    NumericVector testdist = dmat(_,i);
    cont = true;
    
    skippoints(bestdist, penp, n, skip, i, cont, testdist, candidates);
 
    if(cont){ //current set has potential to get better solution
      for (int j = (i+1); j < n; j++)
      {
        if (ISNA(clustx[i]) || ISNA(clustx[j]))
        {
          /* virtual points fulfill the distance condition but are useless for this */
        }else
        {
          //double tempdist = (clustx[i]-clustx[j]) * (clustx[i]-clustx[j]) + (clusty[i]-clusty[j]) * (clusty[i]-clusty[j]); //same as dprime2 but without the special cases of NA points as they are already left out
          double tempdist = testdist[j];
          
          if (tempdist <= 8*penp) //maximum euclidean distance is 2* sqrt(2*penp); square this to get max dprime2 between points = 8*penp
          {
            //two intersection points based on vectors
            double a = sqrt(2.0*penp/tempdist -1.0/4.0);
            double x = (clustx[i] + clustx[j])/2.0 + a*(clusty[j] - clusty[i]);
            double y = (clusty[i] + clusty[j])/2.0 + a*(clustx[i] - clustx[j]);

            NumericVector res = bestPoint(x, y, clustx, clusty, penp, i, j);
            if (res[0] < bestdist)
            {
              bestdist = res[0];
              bestcenterx = res[1];
              bestcentery = res[2];
            }

            x = (clustx[i] + clustx[j])/2.0 - a*(clusty[j] - clusty[i]);
            y = (clusty[i] + clusty[j])/2.0 - a*(clustx[i] - clustx[j]);
            
            res = bestPoint(x, y, clustx, clusty, penp, i, j);
            if (res[0] < bestdist)
            {
              bestdist = res[0];
              bestcenterx = res[1];
              bestcentery = res[2];
            }

          } // both points are close enough for circles to intersect
        } // both points are "real" 
      } //for j; check all intersections
    } // at least half the points are close
  } // for i

  double cost = sum(cross_dprime2(bestcenterx,bestcentery,clustx,clusty,penp)); //should be the same as bestdist
  //Rcout << "Cost: " << cost << " bestdist " << bestdist << std::endl;
  List res = List::create(Named("barycenterx") = bestcenterx, _["barycentery"] = bestcentery, _["cost"] = cost, _["skipped"] = skip);
  return(res);
}

//extra function for checking if a point can be skipped
void skippoints(double bestdist, double penp, int n, int& skip, int i, bool& cont, NumericVector testdist, LogicalVector& candidates){
  NumericVector relevantdist = testdist[candidates]; //this is just for checking the lower bound!
  LogicalVector isclose = (relevantdist <= 8*penp); // in the calculation of the barycenters all points are included
  
  if ((n-sum(isclose))*2*penp >= bestdist) 
  {
    candidates[i] = 0;
    cont = false;
    //IntegerVector sequence = seq(i+1,n-1);
    //NumericVector pointsafter = testdist[sequence];
    NumericVector pointsafter = testdist[seq(i+1,n-1)];
    //Rcout << "testdist " << testdist << " i " << i << " skipped " << sum(pointsafter <= 8*penp ) << std::endl;
    skip+= sum(pointsafter <= 8*penp); //how many points would the algorithm check in the for(j...) loop
  }
  return;
}

//Drezners Algorithm called from outside, barycenter and cost are returned
List vanillaClusterCenterEuclid2(NumericVector clustx, NumericVector clusty, double penp) {
  int n = clustx.length();
  int count = 0; //how many times was an index skipped
  //double nclose; //number of active points NOT IMPLEMENTED YET
  double bestcenterx = clustx[0];
  double bestcentery = clusty[0];
  double bestdist= 2*(n-1)*penp;

  NumericMatrix dmat(n,n);
  dmat = cross_dmat(clustx,clusty);
  for (int i = 0; i < (n-1); i++)
  {
    //test if enough points are close to be better than the empty bary
    NumericVector testdist = dmat(_,i);
    for (int j = (i+1); j < n; j++)
    {
      if (ISNA(clustx[i]) || ISNA(clustx[j]))
      {
        /* virtual points fulfill the distance condition but are useless for this */
      }else
      {
        //double tempdist = (clustx[i]-clustx[j]) * (clustx[i]-clustx[j]) + (clusty[i]-clusty[j]) * (clusty[i]-clusty[j]); //same as dprime2 but without the special cases of NA points as they are already left out
        double tempdist = testdist[j];
        
        if (tempdist <= 8*penp) //maximum euclidean distance is 2* sqrt(2*penp); square this to get max dprime2 between points = 8*penp
        {
          count++;
          //two intersection points based on vectors
          double a = sqrt(2.0*penp/tempdist -1.0/4.0);
          double x = (clustx[i] + clustx[j])/2.0 + a*(clusty[j] - clusty[i]);
          double y = (clusty[i] + clusty[j])/2.0 + a*(clustx[i] - clustx[j]);

          NumericVector res = bestPoint(x, y, clustx, clusty, penp, i, j);
          if (res[0] < bestdist)
          {
            bestdist = res[0];
            bestcenterx = res[1];
            bestcentery = res[2];
          }

          x = (clustx[i] + clustx[j])/2.0 - a*(clusty[j] - clusty[i]);
          y = (clusty[i] + clusty[j])/2.0 - a*(clustx[i] - clustx[j]);
          
          res = bestPoint(x, y, clustx, clusty, penp, i, j);
          if (res[0] < bestdist)
          {
            bestdist = res[0];
            bestcenterx = res[1];
            bestcentery = res[2];
          }

        } // both points are close enough for circles to intersect
      } // both points are "real" 
    } //for j; check all intersections
  } // for i

  double cost = sum(cross_dprime2(bestcenterx,bestcentery,clustx,clusty,penp)); //should be the same as bestdist
  List res = List::create(Named("barycenterx") = bestcenterx, _["barycentery"] = bestcentery, _["cost"] = cost, _["calculations"] = count);
  return(res);
}

//this function calculates the best subset (out of 4 candidates) for calculating the barycenter from
//
NumericVector bestPoint(double x, double y, NumericVector clustx, NumericVector clusty, double penp, int i, int j){
  int n = clustx.length();
  double nhappy = 0.0;
  double meanx = 0.0;
  double meany = 0.0;
  for (int k = 0; k < n; k++)
  {
    if (!ISNA(clustx[k]))
    {
      double dist = dprime2(x,y,clustx[k],clusty[k],2*penp);
      if (dist < 2*penp) //in theory we need also points that are = 2*penp away, but this has probability 0
      {
        nhappy+=1.0;
        meanx += clustx[k];
        meany += clusty[k];
      }
    }
  }

  double candx = (meanx + clustx[i] + clustx[j])/(nhappy+2.0); //these include clustx[i] and clustx[j]
  double candy = (meany + clusty[i] + clusty[j])/(nhappy+2.0); //these include clusty[i] and clusty[j]
  double d = sum(cross_dprime2(candx,candy,clustx,clusty,penp));

  double bestdist = d;
  double bestcenterx = candx;
  double bestcentery = candy;

  candx = (meanx + clustx[j])/(nhappy+1.0);
  candy = (meany + clusty[j])/(nhappy+1.0);
  d = sum(cross_dprime2(candx,candy,clustx,clusty,penp));
  if (d<bestdist)
  {
    bestdist = d;
    bestcenterx = candx;
    bestcentery = candy;
  }

  candx = (meanx + clustx[i])/(nhappy+1.0);
  candy = (meany + clusty[i])/(nhappy+1.0);
  d = sum(cross_dprime2(candx,candy,clustx,clusty,penp));
  if (d<bestdist)
  {
    bestdist = d;
    bestcenterx = candx;
    bestcentery = candy;
  }
  
  if (nhappy >0)
  {
    candx = meanx/nhappy;
    candy = meany/nhappy;
    d = sum(cross_dprime2(candx,candy,clustx,clusty,penp));
    if (d<bestdist)
    {
      bestdist = d;
      bestcenterx = candx;
      bestcentery = candy;
    }
  }
  
  NumericVector res = {bestdist, bestcenterx, bestcentery};
  //res = {bestdist, bestcenterx, bestcentery};

  return(res);
}

//iterated version of the heuristic used in MSM2020
List heuristicCenterEuclid2(NumericVector clustx, NumericVector clusty, double& centerx, double& centery, double penp, bool bounds){

  NumericVector foox;
  NumericVector fooy;
  double cx;
  double cy;
  LogicalVector isclose;
  
  bool change = true;
  int it = 1;

  NumericVector distvec = cross_dprime2(centerx,centery,clustx,clusty,2*penp);
  if (bounds)
  {
    isclose = (distvec <= 2*penp);
  }
  else{
    isclose = (distvec < 2*penp);
  }
  
  
  // if no point is close sample one from the clusterpoints
  if (sum(isclose) == 0)
  {
    int ind = sample(clustx.length(), 1, false)(0)-1;
    cx = clustx[ind];
    cy = clusty[ind];
  }
  else{
    foox = clustx[isclose];
    fooy = clusty[isclose];
    cx = mean(foox);
    cy = mean(fooy);
  }
  //optimClusterCenterEuclid2(clustx[isclose], clusty[isclose], centerx, centery);
  
  
  if (cx == centerx && cy == centery)
  {
    change = false;
  }
  
  while(change && it < 100)
  {
    it++;
    
    double ocx = cx;
    double ocy = cy;
    
    distvec = cross_dprime2(cx,cy,clustx,clusty,2*penp);
    if (bounds)
    {
      isclose = (distvec <= 2*penp);
    }
    else{
      isclose = (distvec < 2*penp);
    }
    foox = clustx[isclose];
    fooy = clusty[isclose];
    cx = mean(foox);
    cy = mean(fooy);
    
    if (cx == ocx && cy == ocy)
    {
      change = false;
    }//if nothing changed
  }//while 
  
  /*Rcout << "dist " << distvec << std::endl; //distance does not change
  distvec = cross_dprime2(cx,cy,clustx,clusty,penp);
  Rcout << "dist " << distvec << std::endl;*/
  double nclose = sum(isclose);
  double cost = sum(distvec);
  List res = List::create(Named("barycenterx") = cx, _["barycentery"] = cy, _["cost"] = cost, _["iterations"] = it, _["closePoints"] = nclose);
  return(res);
}
