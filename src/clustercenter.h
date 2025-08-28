#ifndef CENTERING_H
#define CENTERING_H

#include <Rcpp.h>
using namespace Rcpp;

// Functions for computing "centers" for clusters (on various spaces, for various p, and so on)
void optimClusterCenterEuclid2(NumericVector clustx, NumericVector clusty, double& centerx, double& centery);

// Algorithm1 of Drezner1991
void exactClusterCenterEuclid2(NumericVector clustx, NumericVector clusty, double& centerx, double& centery, double penp);
//Drezner1991 from outside, barycenter and cost are returned; with some features to make it faster
List exactClusterCenterEuclid2(NumericVector clustx, NumericVector clusty, double penp, bool aleph);
//Drezner1991 from outside, barycenter and cost are returned; just the algorithm, no enhancements on the candidates

void skippoints(double bestdist, double penp, int n, int& skip, int i, bool& cont, NumericVector testdist, LogicalVector& candidates); //subroutine of cheking whether a point can be skipped; is lower bound on the objective function larger than the current best; adds number of skipped points to int skip

List vanillaClusterCenterEuclid2(NumericVector clustx, NumericVector clusty, double penp);

//this function calculates the best subset (out of 8 candidates) for calculating the barycenter from
NumericVector bestPoint(double x, double y, NumericVector clustx, NumericVector clusty, double penp, int i, int j);

//subroutine of iterating the heuristic clustercenter
List heuristicCenterEuclid2(NumericVector clustx, NumericVector clusty, double& centerx, double& centery, double penp, bool bounds);

  
#endif
