#include "mex.h"
#include <math.h>
#include <time.h>
#include <sys/time.h>

#define DEBUG 0

/*
 * Find the neighbours in a 3-D cube around each voxel. 
 *
 * (fpereira@princeton.edu)
 
 in:
 
 out:

 assumes:
 - coordinate indexing starts at 1

%   This file is part of Simitar
%
%   Simitar is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%   Simitar is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%
%   You should have received a copy of the GNU General Public License
%   along with Simitar.  If not, see <http://www.gnu.org/licenses/>.
%
*/

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
  /* input */
  double *examples1;
  double *examples2;
  int *labels;             double *plabels; /* original labels */
  int *labelsGroup;        double *plabelsGroup;
  double *voxelsToNeighbours; double *pvoxelsToNeighbours;
  double *numberOfNeighbours; double *pnumberOfNeighbours;
  int nPermutations;       double *pnPermutations;

  /* output */
  double *accuracy;
  double *accuracyCount;
  double *fraction;

  double *similarity;

  /* new ones (used to keep the patterns in each neighbourhood) */

  double *colToCoord;
  double *vptr; double *nptr; double *cptr; double *vnptr;
  double total;
  int neighbourRadius; double *pneighbourRadius; int neighbourMax;
  int pos;

  double *eptr1; double *eptr2;
  double *tptr1; double *tptr2;
  double *mptr1; double *mptr2;
  double *sptr1; double *sptr2;
  int e1; int e2;
  int c1; int c2;
  double correlation; double *sptr;
  int moffset;


  /* rest of variables pertaining to input/output contents */
  int nGroups;
  int itmp;
  int iprev,iseen;
  int group;
  int *groupStarts; double *pgroupStarts; mxArray *groupStarts_mx;
  int *groupEnds;   double *pgroupEnds;   mxArray *groupEnds_mx;
  int *groupSize;   double *pgroupSize;   mxArray *groupSize_mx;
  int *groupNclass; double *pgroupNclass; mxArray *groupNclass_mx;
  int **groupPermutations;
  int p;
  int nClasses;
  int n,m;
  int n1,m1,n2,m2;
  int g,c,gidx;
  int *gptr;    /* group pointer */
  double *eptr, *e2ptr; /* example pointer */
  int ngc;
  double *sumX1; mxArray *sumX1_mx; double *x1ptr;
  double *sumX2; mxArray *sumX2_mx; double *x2ptr;
  double *mptr,*m2ptr;
  double *exampleLargestValue;
  double *examplesCorrect; mxArray *examplesCorrect_mx;
  double *examplesSquared;
  double *means;
  double *meansSquared;
  double *vars;
  double *sumX1train;
  double *sumX2train;
  int gtest, gtrain;
  int *nc;
  int ntest, ntrain, ntrainmo;
  int nn;
  int neighbour;
  double voxelValue;

  
  /* everything else */
  int     e,i,ig,j,k,v;

  double *X,*M,*Ez,*Vz,*sigma;
  double *pEz,*pVz,*pX,*pM,*pL;
  double *pEzM, *pEzM2;
  double **mEz,**mVz,**mX,**mM,**mL;
  double fixed;
  double *log_px;
  double *log_px_pd;
  mxArray *tmp1_mx,*tmp2_mx,*sigma2_mx,*logsigma_mx,*pEzM_mx,*pEzM2_mx;
  double  *tmp1,*tmp2,*sigma2,*logsigma;
  int     status,mrows,ncols;
  int     cX,rX,cM,rM,cEz,rEz,cVz,rVz,cS,rS;
  int     idx;
  char    buffer[10];

  /*  check for proper number of arguments */
  /* NOTE: You do not need an else statement when using mexErrMsgTxt
     within an if statement, because it will never get to the else
     statement if mexErrMsgTxt is executed. (mexErrMsgTxt breaks you out of
     the MEX-file) */

  if(nrhs!=2) 
    mexErrMsgTxt("two inputs required");

  if(nlhs!=2) 
    mexErrMsgTxt("two outputs required");

  /* check to make sure the input arguments are double matrices */
  if( !mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) ) {
    mexErrMsgTxt("inputs must be matrices");
  }

  /*  create a pointer to the input matrices */
  colToCoord       = mxGetPr(prhs[0]);
  pneighbourRadius = mxGetPr(prhs[1]);
  neighbourRadius = (int) round(*pneighbourRadius);
  m = mxGetN(prhs[0]);

  neighbourMax = (int) round(pow(2*neighbourRadius+1,3))-1;
  
  mexPrintf("finding neighbours within radius %d\n",neighbourRadius);fflush(stdout);

  /*
    mexPrintf("#voxels=%d radius=%d max=%d\n",m,neighbourRadius,neighbourMax);fflush(stdout);
    return;
  */

  /* mexPrintf("pvn = %d,%d pnn = %d,%d\n",mxGetN(prhs[3]),mxGetM(prhs[3]),mxGetN(prhs[4]),mxGetM(prhs[4]));fflush(stdout);*/

  /*  create output matrices */

  plhs[0] = mxCreateDoubleMatrix(m,1,mxREAL);
  /*
    plhs[1] = mxCreateDoubleMatrix(m,neighbourMax,mxREAL);
  */
  plhs[1] = mxCreateDoubleMatrix(neighbourMax,m,mxREAL);

  numberOfNeighbours = mxGetPr(plhs[0]);
  voxelsToNeighbours = mxGetPr(plhs[1]);

  /* mexPrintf("%p\t%p\n",numberOfNeighbours,voxelsToNeighbours);fflush(stdout);return; */

  /*
   * Main loop
   */

  /*
    vptr  = colToCoord;
    for (i=0; i<m; i++)
    {
    mexPrintf("%1.2f,%1.2f,%1.2f\n",vptr[0],vptr[1],vptr[2]);fflush(stdout);
    vptr = vptr + 3;
    }
    return;
  */

  cptr  = numberOfNeighbours;
  vptr  = colToCoord;

  for (i=0; i<m; i++)
    {
      nptr = colToCoord;
      cptr[i] = 0;

      vnptr = voxelsToNeighbours + i*neighbourMax;

      if (remainder((double)i,1000) == 0) {
	mexPrintf(" %d",i); fflush(stdout);
      }

      /*
	for (j=0;j<neighbourMax;j++) {*vnptr=0;}
	vnptr = voxelsToNeighbours + (i-1)*neighbourMax;;
      */

      /*
	mexPrintf("%p\t%p\t%2.0f\t%p\t%p\n",vptr,nptr,cptr[i],voxelsToNeighbours,vnptr);fflush(stdout);
      */

      /* mexPrintf("voxel %d (%2.0f,%2.0f,%2.0f) adding\n",i,vptr[0],vptr[1],vptr[2]);fflush(stdout);*/

      for (j=0; j<m; j++) {

	/* detect whether this voxel is within the radius */
	total = 0;
	if (abs(vptr[0]-nptr[0]) > neighbourRadius) {
	  total++;
	} else {
	  if (abs(vptr[1]-nptr[1]) > neighbourRadius) {
	    total++;
	  } else {
	    if (abs(vptr[2]-nptr[2]) > neighbourRadius) {
	      total++;
	    }
	  }
	}

	/* mexPrintf("%2.0f,%2.0f,%2.0f\t%2.0f,%2.0f,%2.0f\t%2.0f\t%2.0f\n",vptr[0],vptr[1],vptr[2],nptr[0],nptr[1],nptr[2],total,cptr[i]);fflush(stdout); */

	if ((total == 0) && (i!=j)) {
	  /*mexPrintf("\t%d (%2.0f,%2.0f,%2.0f)\n",j,nptr[0],nptr[1],nptr[2]);fflush(stdout);*/
	  /* if it is, add to the neigbours */
	  /*
	   *vnptr = j + 1;
	   vnptr++;
	  */
	  pos = (int) round(cptr[i]);
	  vnptr[pos] = j+1;
	  cptr[i] = cptr[i] + 1;
	}

	nptr = nptr + 3;
      }


      vnptr = voxelsToNeighbours + i*neighbourMax;

      /*
	mexPrintf("%d\t",i); fflush(stdout);
	for (j=0; j<neighbourMax; j++) {
	  mexPrintf("%2.0f ",vnptr[j]);fflush(stdout);
	}; mexPrintf("\n"); fflush(stdout);
      */
      vptr  = vptr + 3;;

      /*
	mexPrintf("%p\t%p\t%2.0f\t%p\t%p\n",vptr,nptr,cptr[i],voxelsToNeighbours,vnptr);fflush(stdout);
      */

      /* return; */

      /* mexPrintf("%d\t%p\t%p\t%d\n",i,vptr,nptr,cptr[i]);fflush(stdout); */
      /*return;*/

      /* mexPrintf("aqui!\n");fflush(stdout); return; */
    }

  mexPrintf("\nfound all neighbours\n");fflush(stdout);

  /*
  for (i=0; i<m; i++) {
    vnptr = voxelsToNeighbours + i*neighbourMax;
    mexPrintf("%d\t",i); fflush(stdout);
    for (j=0; j<neighbourMax; j++) {
      mexPrintf("%2.0f ",vnptr[j]);fflush(stdout);
    }; mexPrintf("\n"); fflush(stdout);
  }
  */

  return;
}
