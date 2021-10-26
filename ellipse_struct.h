#ifndef _ELLIPSE_STRUCT_H_
#define _ELLIPSE_STRUCT_H_ 1
struct ellipse_struct {
  /*
%  ellipse_struct - structure that defines the best fit to an ellipse                       
%                       a           - sub axis (radius) of the X axis of the non-tilt ellipse   
%                       b           - sub axis (radius) of the Y axis of the non-tilt ellipse   
%                       phi         - orientation in radians of the ellipse (tilt)              
%                       x0          - center at the X axis of the non-tilt ellipse              
%                       y0          - center at the Y axis of the non-tilt ellipse              
%                       x0_in       - center at the X axis of the tilted ellipse                
%                       y0_in       - center at the Y axis of the tilted ellipse                
%                       long_axis   - size of the long axis of the ellipse                      
%                       short_axis  - size of the short axis of the ellipse                     
%                       status      - status of detection of an ellipse 1 = found,
%                                     0 = not an ellipse.
  */
  double a; 
  double b;
  double phi;
  double x0;
  double y0;
  double x0_in;
  double y0_in;
  double long_axis;
  double short_axis;
  int    status;
}
;
#endif
