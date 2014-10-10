/**
 Rectangular Histogram of Gradients (HoG) descriptor mex implementation

 See Histograms of Oriented Gradients for Human Detection
 [Dalal, Triggs, 2005]
 for more details.

 Usage:

 hog(image, thetabins, [cellx celly], [gridx gridy], boolean signed_gradient)

 image size is the patch size (ex. 36)

 image should be a single valued grayscale image.

 theta bins should be the number of angular bins (typically 9)

 cellx and celly should represent the number of subdivions of the image where
 the number of pixels for each cell is given by
 width / cellx x height / celly
 This should probably be a factor of the patch size
 (ie. if patch is 36, this should be [6 6] , [9 9], etc.

 gridx and gridy should represent the number of neighboring cells to
 accumulate when counting the histogram ex. [1, 1]

 signed gradient, if true respects the sign when binning ie. [0..360],
 else the sign is dropped [0..180].

 returns a feature that is size 1 x thetabins*cellx*celly
 (default is 9x6x6=324)

 Default:

 hog(image, 36, 9, [6 6], [1 1], false)

 Paul Sastrasinh (psastras) 2011
 */

#include <mex.h>
#include <matrix.h>
#include <math.h>
#include "hog.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  // Read the input arguments
  if(nrhs < 1)
  {
    mexErrMsgTxt("Expected at least one input, the single channel image.");
  }

  mxClassID imgid = mxGetClassID(prhs[0]);
  if (imgid != mxSINGLE_CLASS)
    mexErrMsgTxt("Expected first argument to be single (grayscale image).");

  float *imageData = (float *)mxGetData(prhs[0]);
  const int *dims = mxGetDimensions(prhs[0]);
  int iWidth = dims[1]; int iHeight = dims[0];

  // number of angle bins
  int ntheta = 9;

  // cellsize (image size / pixels)
  int ncellsx = 6;
  int ncellsy = 6;

  // block size (cells)
  int ngridx = 1;
  int ngridy = 1;

  // controls the interpolation level
  // 0 - no interpolation
  // 1 - bilinear (not complete)
  // 2 - trilinear (not implemented yet)
  // Note: To be honest, I don't think the interpolation is even necessary -
  // overlapping bins (blocks) coupled with gaussian weighting does
  // pretty much the exact same thing as interpolation...getting some hard
  // numbers on this would be interesting.
  int interplevel = 0;

  bool signedgradient = false;

  if(nrhs >= 2)
  {
    ntheta = (int)mxGetScalar(prhs[1]);
  }

  if(nrhs >= 3)
  {
    double *ncells = (double *)mxGetData(prhs[2]);
    ncellsx = (int)ncells[0];
    ncellsy = (int)ncells[1];
  }

  if(nrhs >= 4)
  {
    double *ngrids = (double *)mxGetData(prhs[3]);
    ngridx = (int)ngrids[0];
    ngridy = (int)ngrids[1];
  }

  if(nrhs >= 5) {
    signedgradient = (int)mxGetScalar(prhs[4]);
  }

  if(nrhs >= 6) {
    mexWarnMsgTxt("There appear to be extra arguments...proceeding anyway.");
  }

  if(interplevel == 1)
    mexWarnMsgTxt("Bilinear interpolation only partially implemented.\n");
  if(interplevel > 2) {
    mexWarnMsgTxt("Invalid interpolation level.  Selecting 0 instead.\n");
    interplevel = 0;
  }

  plhs[0] = mxCreateNumericMatrix(1, ntheta * ncellsx * ncellsy,
                                  mxSINGLE_CLASS, mxREAL);
  float *H = (float *)mxGetData(plhs[0]);
  // grad = new GradientResponse[iWidth * iHeight];
  normalize(imageData, imageData, iWidth, iHeight);
  hog(imageData, H, iWidth, iHeight, ntheta, ncellsx, ncellsy,
      ngridx, ngridy, signedgradient, interplevel);
  // delete [] grad;
}

