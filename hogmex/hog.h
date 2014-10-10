#ifndef __HOG_H__
#define __HOG_H__

#include <cmath>
#define MIN(x, y) (x) < (y) ? (x) : (y)
#define MAX(x, y) (x) > (y) ? (x) : (y)

// struct to hold the gradient data during the convolution step
struct GradientResponse
{
  float theta, G; // x gradient, y gradient, angle, magnitude
};

GradientResponse *grad = 0;

// image variance normalization
void normalize(const float *imageDataIn, float *imageDataOut,
               const int iWidth, const int iHeight)
{
  float mu = 0.f, m2 = 0.f;

  // one pass variance and mean computation 
  for(int i=0;i<iWidth*iHeight;i++)
  {
    float delta = imageDataIn[i] - mu;
    mu += delta / (float)(i+1);
    m2 += delta * (imageDataIn[i] - mu);
  }

  const float stdev = sqrtf(m2 / (float)(iWidth*iHeight)) + 0.01f;

  for(int i=0;i<iWidth*iHeight;i++)
  {
    imageDataOut[i] = (imageDataIn[i] - mu) / stdev;
  }
}

/* hog */
void hog(const float *imageDataIn, float * descriptorOut,
         const int iWidth, const int iHeight,
		 const int ntheta,
         const int ncellsx, const int ncellsy,
         const int ngridx,  const int ngridy,
		 const bool signedgradient, const int interplevel)
{
  const int celldimy = (int)(iHeight / (float)(ncellsy + ngridx));
  const int celldimx = (int)(iWidth / (float)(ncellsx + ngridy));
  // gaussian kernel
  const int ydim = celldimy*(ngridy+1);
  const int xdim = celldimx*(ngridx+1);
  float gauss[xdim*ydim];

  float sigma = MAX(xdim, ydim) / 1.25f;
  float s = 0.f;
  for(int y=0, i=0; y<ydim; y++) {
    for(int x=0; x<xdim; x++, i++) {
      gauss[i]=(1./(2*M_PI*sigma*sigma))*exp(-(x*x+y*y)/(2*sigma*sigma));
      s += gauss[i];
    }
  }
  for(int i=0; i<xdim*ydim; i++) { gauss[i] /= s; }
  
  // convolve gradients
  for(int y=0; y<iHeight; y++)
  {
    for(int x=0; x<iWidth; x++)
    {
      const float &xl = x-1<0       ?0.:imageDataIn[(x-1)*iHeight+y];
      const float &xr = x+1>=iWidth ?0.:imageDataIn[(x+1)*iHeight+y];
      const float &yt = y-1<0       ?0.:imageDataIn[x*iHeight+y-1];
      const float &yb = y+1>=iHeight?0.:imageDataIn[x*iHeight+y+1];

      const float &gx = xr - xl;
      const float &gy = yb - yt;

      GradientResponse *r = &grad[x * iHeight + y];
      r->G = sqrtf(gx*gx + gy*gy);
      r->theta = atan2(gy, gx);
      if(!signedgradient && r->theta < 0) {
        r->theta += M_PI;
      }
    }
  }

  // contrast over overlapping blocks
  for(int y=0,i=0; y<ncellsy; y++) {
  const int yloc = celldimy * y;
  for(int x=0; x<ncellsx; x++, i++)
  {
    const int xloc = celldimx * x;
    // current cell in the histogram
    float *hist = &descriptorOut[i * ntheta]; 
    // For each pixel in each cell
    for(int yy=yloc,z=0; yy < celldimy*(ngridy+1)+yloc; yy++) {
    for(int xx = xloc; xx < celldimx*(ngridx+1)+xloc; xx++, z++)
    {
      const GradientResponse *p = &grad[xx * iHeight + yy];
      const float &theta =
          signedgradient ? (p->theta + 3.14159265358979f) * 0.5f :
          p->theta; //[0...PI]
      const float &w = gauss[z];

      if(interplevel == 0) {  // no interpolation
        const int bin =
            (int)(MAX(MIN(theta / (3.14159265358979f / ntheta),
                          ntheta - 1), 0));
        hist[bin] += w*p->G;
      // } else if(interplevel == 1) {
      //   // bilinear interpolation on the angles
      //   float t = theta / (3.14159265358979f / ntheta);
      //   const int z1 = ((int)(t)) % ntheta;
      //   const int z2 = ((int)(t+1)) % ntheta;
      //   t -= z1;
      //   hist[z1] += (1.f-t) * p->G * w;
      //   hist[z2] += t * p->G * w;
      //   // missing spatial interpolation ...
      // }
      // else if(interplevel == 2)
      //   // trilinear interpolation on cells(x,y) and angles (theta)
      // {
      //   mexErrMsgTxt("Trilinear interpolation not yet implemented.");
      }
    }
    }

    // normalize the cell
    float sum = 0.0;
    for(int i=0; i<ntheta; i++){ sum += hist[i]*hist[i]; }
    sum = sqrt(sum) + 0.01f;
    for(int i=0; i<ntheta; i++){ hist[i] /= sum; }
  }
  }
}

#endif // __HOG_H__


