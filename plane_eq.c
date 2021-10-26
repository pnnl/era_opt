#include "system_includes.h"
#include "plane_eq.h"
void plane_eq(double *points,double *abcd, int *colinear_p) {
  /*
    given three points in 3 dimensional space stored in points 
    in (x1,y1,z1) (x2,y2,z2) (x3,y3,z3) oder
    find the coefficients of the plane passing through the 
    three points
    Adapted from my findabcd.f routine from geosim.
    We avoid using determiniants as these points will come from
    long narrow triangles and thus the system may have a high condition 
    number.
      !
      !     let the points be represented by (x_i, y_i, z_i) then 
      !     solving the following system of equations for a,b,c,d
      !     yeilds the necessary coefficients.
      !
      !     x1 y1 z1 1    a      0
      !    (x2 y2 z2 1)  (b)  = (0)
      !     x3 y3 z3 1    c      0
      !      ?  ?  ? ?    d      1
      !    
      ! The general challange is what to use for the scaling coefficients
      ! listed as ? in the last row. While it is tempting to use
      ! all 1's, this breaks down in some case, e.g. the event that the 
      ! solution should be  x = y with c and d being 0
      ! but it is sufficient to take the last row to be some unit vector
      ! of all zeros except one element. Which element we can determine
      ! as we solve and factor the system with total pivoting.
      !
      ! on success (the points are not colinear) we return  colinear = 0, 
      ! If the points are co-linear we return a value of 1 in colinear
      !
      ! jcol and irow are used to keep pivoting information.
      ! we need to make an adjustment as d should be exactly 0.
  */
  double coeffs[12];
  double xp[4];
  double vmax;
  double vmax2;
  double absvmax;
  double absvmax2;
  double tv;
  double abstv;
  double tv2;
  double abstv2;
  double det;
  double emult;
  double ac;
  double recip_vmax;
  double recip_vmax2;
  double a;
  double b;
  double c;
  double d;
  double pscale;
  double p13;
  double p12;
  double p22;
  double p23;

  int    jcol[4];
  int    irow[4];

  int    i;
  int    j;

  int    k;
  int    kk;

  int    ip0;
  int    jp0;

  int    ip1;
  int    jp1;

  int    ip2;
  int    jp2;

  int    ic;
  int    jc;

  int    iprowpos;
  int    ipcolpos;

  int    ie_rowpos;
  int    ie_rowcol;

  int    colinear;
      
  colinear = 0;
  a = 0.0;
  b = 0.0;
  c = 0.0;
  d = 0.0;
  for (i=0;i<4;i++) {
    irow[i] = i;
    jcol[i] = i;
  }
  /*
  ! spread out the coefficients out points into linear 3 x 4 array
  ! inserting a fourth column of 1's (stored in row major order.
  */
  k = 0;
  kk = 0;
  for (j=0;j<3;j++) {
    for (i=0;i<3;i++) {
      coeffs[k] = points[kk];
      kk += 1;
      k  += 1;
    }
    coeffs[k] = 1.0;
    k += 1;
  }
  /*
    Find the maximal pivotal element in the first three columns and
    rows ( we don't let column 3 of all 1's participate in the pivoting)
  */
  vmax = coeffs[0];
  absvmax = vmax;
  if (absvmax < 0.0) {
    absvmax = -absvmax;
  }
  /*
    ip0,jp0 will record to the row and column number of vmax,
    To be swapped in to the 0,0 row and column.
    ic, and jc will be used to track the row and column of 
    the current element we are examining.
  */
  ip0  = 0;
  jp0  = 0;
  ic   = 0;
  jc   = 0;
  for (k=1;k<12;k++) {
    tv = coeffs[k];
    abstv = tv;
    if (abstv < 0.0) abstv = - abstv;
    jc = jc + 1;
    if (jc == 4) {
      jc = 0;
      ic = ic + 1;
    }
    /*
      we don't let column 3 participate in the pivoting
    */
    if (jc != 3) {
      if (abstv > absvmax) {
	/*
	  New largest element, record its position in ic,jc
        */
	absvmax = abstv;
	vmax    = tv;
	ip0     = ic;
	jp0     = jc;
      }
    }
  }
  /*
    Pivot maximum element into 0,0 position.
  */
  if (ip0 != 0) {
    /*
      if maximum was not in row 0, swap row pointers 
      of row 0 and row ip0.
    */
    irow[0]   = ip0;
    irow[ip0] = 0;
  }
  if (jp0 != 0) {
    jcol[0]    = jp0;
    jcol[jp0]  = 0;
  }
  det = vmax;
  if (absvmax == 0.0) {
    colinear = 1;
  } else {
    recip_vmax = 1.0/vmax;
    /*
      Scale row ip0 by recip_vmax
    */ 
    iprowpos = ip0 << 2;
    /*
      First scale the 1 in column 3.
    */
    coeffs[iprowpos + 3] = recip_vmax;
    ipcolpos  = 0;
    for (k = iprowpos;k <iprowpos+3;k++) {
      if (ipcolpos != jp0) {
	coeffs[k] = coeffs[k] * recip_vmax;
      }
      ipcolpos = ipcolpos + 1;
    }	
    /*
      now eliminate the elements in column jp0 in rows other 
      than ip0 tracking their resultant maximum absolute value
      as we go - as that will be the next pivot element
    */
    absvmax2 = 0.0;
    ip1   = irow[1];
    jp1   = jcol[1];
    for (i=1;i<3;i++) {
      ie_rowpos = irow[i] << 2;
      emult    = -coeffs[ie_rowpos + jp0];
      /*
	subtract emult * row ip0 from row irow[i].
	note as both rhs elements are 0 the rhs values do not change
	when subtracting a multiplie of the pivot row from the
	elimination row.
      */
      for (j=0;j<4;j++) {
	if (j != jp0) {
	  tv2 = coeffs[ie_rowpos + j];
	  tv2 = tv2 + (emult * coeffs[iprowpos+j]);
	  coeffs[ie_rowpos + j] = tv2;
	  /*
	    don't let column 3 be a pivot column
	  */
	  if (j != 3) {
	    abstv2 = tv2;
	    if (abstv2 < 0.0)abstv2 = -abstv2;
	    if (abstv2 > absvmax2) {
	      absvmax2 = abstv2;
	      vmax2 = tv2;
	      ip1   = irow[i];
	      jp1   = j;
	    }
	  }
	}
      } /* end for j */
    } /* end for i */
    /*
      now if absvmax2 == 0.0 we have colinearity
    */
    if (absvmax2 == 0.0) {
      colinear = 1;
    } else {
      if (ip1 != irow[1]) {
	irow[2] = irow[1];
	irow[1] = ip1;
      }
      if (jp1 != jcol[1]) {
	jcol[2] = jcol[1];
	jcol[1] = jp1;
      }
      det = det * vmax2;
      /*
	scale the new pivot row by recip_vmax2
      */
      recip_vmax2 = 1.0/vmax2;
      iprowpos = ip1 << 2;
      /*
	columns to be scaled are column 3 and column 3 - jp0 - jp1
      */
      p13 = coeffs[iprowpos+3]; 
      p13 = p13 * recip_vmax2;
      coeffs[iprowpos+3] = p13;
      jp2 = 3 - jp0 - jp1;
      p12 = coeffs[iprowpos + jp2];
      p12 = p12 * recip_vmax2;
      coeffs[iprowpos+jp2] = p12;
      /*
	now subtract a multiple of the new pivot row, ip1, 
	from row ip2 = 3-ip0-ip1 , the last of the first
	three rows. to be eliminated.
      */
      ie_rowpos = (3 - ip0 - ip1) << 2;
      emult    =  - coeffs[ie_rowpos + jp1];
      p22 = coeffs[ie_rowpos + jp2];
      p22 = p22 + emult * p12;
      coeffs[ie_rowpos + jp2] = p22;
      p23 = coeffs[ie_rowpos + 3];
      p23 = p23 + emult * p13;
      coeffs[ie_rowpos+3] = p23;
      /*
	If p22 is large then the correct fourth row to use
	is (0,0,0,1) 
	But for small p22, say p22 < 1, (if p22 is 0 then 
	d is actually 0) then use a vector of zeros
	and a 1 in column jp2.
      */
      /*
      det = p22 * det;
	d should always be zero.
      */
      /* 
	 THe plane contains the origin, so d is 0.
	 0  0  p22  p23   0
	 0  0    1    0   1  use this as the third row and rhs.
	*/
      xp[3] = 0.0;
      xp[2] = 1.0;
      det   = 0.0;
      /*
      if (p22 == 0.0) { 
	// 
	//   THe plane contains the origin, so d is 0.
	//   0  0  p22  p23   0
	//   0  0    1    0   1  use this as the third row and rhs.
	//
      	xp[3] = 0.0;
	xp[2] = 1.0;
      } else {
	xp[3] = det;
	//
	//  0  0  p22  p23   0
	//  0  0    0    1   det  use this as the third row and rhs,
	//                        so as to match scaling of the 
	//			determinant method.               
	//
	xp[2] = -p23*det/p22;
      }
      */
      /*
	finish back solve
      */             
      iprowpos = ip1 << 2;
      xp[1] = 0.0 - (xp[3] * coeffs[iprowpos+3]);
      xp[1] = xp[1] - (xp[2] *  coeffs[iprowpos + jp2]);
      iprowpos = ip0 << 2;
      xp[0] = 0.0 - (xp[3] * coeffs[iprowpos+3]);
      xp[0] = xp[0] - (xp[2] * coeffs[iprowpos + jp2]);
      xp[0] = xp[0] - (xp[1] * coeffs[iprowpos + jp1]);
      abcd[jp0] = xp[0];
      abcd[jp1] = xp[1];
      abcd[jp2] = xp[2];
      abcd[3]   = xp[3];
      a = abcd[0];
      b = abcd[1];
      c = abcd[2];
      d = abcd[3];
      /*
	Now we want to normalize the coefficients.
      */
      pscale = (a * a) + (b * b) + (c * c);
      pscale = sqrt(pscale);
      /*
      if (d > 0.0) {
	pscale = -pscale;
      } else {
        if (d == 0.0) {
      */
      if (c < 0.0) {
	pscale = - pscale;
      } else {
	if (c == 0.0) {
	  if (b <0.0) {
	    pscale = -pscale;
	  }
	}
      }
      /*
	}
      }
      */
      pscale = 1.0/pscale;
      for (i=0;i<4;i++) {
	abcd[i] = abcd[i] * pscale;
      }
    }
  }
  *colinear_p = colinear;
}
      
            
            

