#include <stdio.h>
#include <stdlib.h>
#include <math.h>

static const uint64_t MATRIX_FIELDS=26;
#define _btw 0
#define _bth 1
#define _bbw 2
#define _bbh 3
#define _brw 4
#define _brh 5
#define _bfw 6
#define _bfh 7
#define _baw 8
#define _bah 9
#define _blen 10
#define _bh 11
#define _nw 12
#define _nht 13
#define _nhb 14
#define _fnw 15
#define _fnh 16
#define _rad 17
#define _acenter 18
#define _astart 19
#define _aend 20
#define _xend 17
#define _ycenter 18
#define _ystart 19
#define _yend 20
#define _parent 21
#define _max_leaf_idx 22
#define _is_leaf 23
#define _visible 24
#define _visible_faces 25

//void *test(uint32_t rows, uint32_t cols, double *a, void *callback(uint32_t index));

void *test(uint64_t rows, uint64_t cols, double *a, void *callback(uint64_t index)){
  uint64_t r = 0;
  uint64_t c = 0;\
  uint64_t i = 0;

  callback(r); 

  
  for (r=0; r<rows; r++){
    i = (cols * r) + 2;
    a[i] += 1;
    if (r<10){
      callback(r);
        }
    /* for (c=0; c<cols; c++){ */
    /*   i = (cols * r) + c; */
    /*   //printf("%d %d %f \n", r, c, a[i]); */
    /*   a[i] += 1.0; */
    /* }  */
  }  
}


void draw_circular_region(double cx, double cy,
                          double rx, double ry,
                          double rw, double rh,
                          double zoom_factor,
                          uint64_t nnodes,
                          double *imgdata,
                          int update_visibility,
                          void *py_draw_node(uint64_t _)){
 
  uint64_t nid, idx, parent_idx;
  uint64_t collapse_until = 0;
  
  for (nid=0; nid<nnodes; nid++){
    idx =  (MATRIX_FIELDS * (uint64_t) nid);
    parent_idx = (uint64_t) imgdata[idx + _parent];
    
    if (update_visibility == 1){
      imgdata[idx + _visible] = 0;
      
      if (collapse_until !=0 && nid <= collapse_until){        
        continue;
      }
      /* else if (collapse_until > nid){ */
      /*   collapse_until = 0; */
      /* } */
    
      if (imgdata[idx + _is_leaf] == 0 && (imgdata[idx + _fnh] * zoom_factor) <= 2.0){
        collapse_until = imgdata[idx + _max_leaf_idx];
      }      
      imgdata[idx + _visible] = 1;      
    }
    if (imgdata[idx + _visible] == 1){
      if (RectAnnulusIntersect(cx, cy,
                               imgdata[parent_idx + _rad],
                               imgdata[idx + _rad],
                               imgdata[idx + _astart],
                               imgdata[idx + _aend],
                               rx, ry, rw, rh) > 1){
        py_draw_node(nid);        
      }
    }
    /* if intersects(rx, ry, rw, rh, r1, r2, a1, a2){         */
    /*   } */     
  }
}



/// Check Intersection Functions
/// To test, call:
///
/// rx = rectangle top-left X
/// ry = rectangle top-left Y
/// rw = rectangle width 
/// rh = rectangle height
/// r1 = annulus inner radius
/// r2 = annulus outer radius
/// a1 = annulus arc angle start radians
/// a2 = annulus arc angle end in radians 
///
/// i = RectAnnulusIntersect(rx, ry, rw, rh, r1, r2, a1, a2)
///
/// intersection if i > 0


void cart2angle(double x, double y, double *a){  
  //  angles in qt4 are clockwise. 3 o'clock is 0 
  //      270
  // 180        0
  //      90

  if (x>0 && y<0){
    *a = - atan(y/x);
  }
  else if (x<0 && y<0){
    *a = M_PI - atan(y/x);
  }
  else if (x<0 && y>0){
    *a = M_PI - atan(y/x);
  }
  else if (x>0 && y>0){
    *a = (2*M_PI) - atan(y/x);
  }

}

double to_degrees(double radians) {
  return radians * (180.0 / M_PI);  
}

double to_radians(double degrees) {
  return (degrees * M_PI)/ 180.0;
}


void tocart(double cx, double cy, double rad, double angle, double *x, double *y){
  *x = (rad * cos(angle));
  *y = (rad * sin(angle));
  //printf("TOCART %f %f %f %f \n",rad, angle, *x, *y);  
}

void topolar(double x, double y, double *rad, double *angle){
  double kk;
  *rad = sqrt(x*x + y*y);
  cart2angle(x, y, angle);
  //printf("TOPOLAR %f %f %f %f \n",x, y, *rad, kk);  
}


int cartPointInAnnulus(double rad1, double rad2, double astart, double aend,
                       double px, double py){
  double pr, pa;
  topolar(px, py, &pr, &pa);

  if (pr > rad2 || pr < rad1){
    return 0;
  }
  if (pa < astart || pa > aend){
    return 0;
  }
  return 1;
}

int polarPointInRect(double cx, double cy,
                     double r1x, double r1y,
                     double r4x, double r4y,
                     double rad, double angle){
  double x, y; 
  tocart(cx, cy, rad, angle, &x, &y);
  //printf("%f,%f %f,%f \n\n", r1x, r1y, r4x, r4y);
  
  if (x < r1x || x > r4x){
    return 0;
  }
  if (y > r1y ||  y < r4y){
    return 0;
  }
  return 1;
}


int RectAnnulusIntersect(double rad1, double rad2,
                         double astart, double aend, double rx, double ry,
                         double rw, double rh, double cx, double cy){
  
  double r1x, r1y, r2x, r2y, r3x, r3y, r4x, r4y;
  //astart = to_radians(astart);
  //aend = to_radians(aend);
  
  // translate rect coords to main circle center
  r1x = rx - cx;
  r1y = cy - ry;
  
  r2x = r1x + rw;
  r2y = r1y;
  
  r3x = r1x;
  r3y = r1y - rh;
  
  r4x = r1x + rw;
  r4y = r1y - rh;

  //printf("RECT %f,%f   %f,%f   %f,%f   %f,%f\n", r1x, r1y, r2x, r2y, r3x, r3y, r4x, r4y);
  
  // Check if any vertex of rect is wihin annulus
  if      (cartPointInAnnulus(rad1, rad2, astart, aend, r1x, r1y) == 1){
    return 11;
      }
  else if (cartPointInAnnulus(rad1, rad2, astart, aend, r2x, r2y) == 1){
    return 12;
      }
  else if (cartPointInAnnulus(rad1, rad2, astart, aend, r3x, r3y) == 1){
    return 13;
      }
  else if (cartPointInAnnulus(rad1, rad2, astart, aend, r4x, r4y) == 1){
    return 14;
      }
    
  // Check if any vertex of annulus is within rect
  else if (polarPointInRect(cx, cy, r1x, r1y, r4x, r4y, rad1, astart)){
    return 21;
  }
  else if (polarPointInRect(cx, cy, r1x, r1y, r4x, r4y, rad2, astart)){
    return 22;
      }
  else if (polarPointInRect(cx, cy, r1x, r1y, r4x, r4y, rad1, aend)){
    return 23;
  }
  else if (polarPointInRect(cx, cy, r1x, r1y, r4x, r4y, rad2, aend)){
    return 24;
  }
  
  // Check if any rect segment intersect with the annulus arcs
  else if (SegmentArcIntersect(rad1, astart, aend, r1x, r1y, r2x, r2y)==1){
    return 31;
  }
  else if (SegmentArcIntersect(rad1, astart, aend, r3x, r3y, r4x, r4y)==1){
    return 32;
  }
  else if (SegmentArcIntersect(rad1, astart, aend, r1x, r1y, r3x, r3y)==1){
    return 33;
  }
  else if (SegmentArcIntersect(rad1, astart, aend, r2x, r2y, r4x, r4y)==1){
    return 34;
  }
  else if (SegmentArcIntersect(rad2, astart, aend, r1x, r1y, r2x, r2y)==1){
    return 35;
  }
  else if (SegmentArcIntersect(rad2, astart, aend, r3x, r3y, r4x, r4y)==1){
    return 36;
  }
  else if (SegmentArcIntersect(rad2, astart, aend, r1x, r1y, r3x, r3y)==1){
    return 37;
  }
  else if (SegmentArcIntersect(rad2, astart, aend, r2x, r2y, r4x, r4y)==1){
    return 38;
  }
  // check if any rect segment intersect with the side segments of the annulus


  return 0;
}
  
void intersects(double cx, double cy,
                double rx, double ry,
                double rw, double rh,  
                uint64_t nnodes, double *matrix){
  
  uint64_t i, idx, parent_idx;
  uint64_t matches=0;
  int match=0;
  
  for (i=0; i<nnodes; i++){
    idx =  (MATRIX_FIELDS * (uint64_t) i); 
    parent_idx = (uint64_t) matrix[idx + _parent];

    match = RectAnnulusIntersect(cx, cy,
                         matrix[parent_idx + _rad],
                         matrix[idx + _rad],
                         matrix[idx + _astart],
                         matrix[idx + _aend],
                         rx, ry, rw, rh);
    if (match > 0){
      matches = matches +1;
    }
    
  }
  printf("MATCHES xxxxxxxxxxxxxx: %d\n", matches);
  
}
  

int SegmentArcIntersect(double R,
                        double astart, double aend,
                        double aX, double aY,
                        double bX, double bY){

  double a, iX, iY, t, t1, t2, dt, dist, dX, dY, dl, nearestX, nearestY, cX, cY;
  cX = 0;
  cY = 0;
  
  dX = bX - aX;
  dY = bY - aY;
  if ((dX == 0) && (dY == 0))
    {
      // A and B are the same points, no way to calculate intersection
      return;        
    }
  
  dl = (dX * dX + dY * dY);
  t = ((cX - aX) * dX + (cY - aY) * dY) / dl;
  
  // point on a line nearest to circle center
  nearestX = aX + t * dX;
  nearestY = aY + t * dY;

  dist = (double) sqrt(pow((cX-nearestX),2)+ pow((cY-nearestY),2));
    
  // dist = point_dist(nearestX, nearestY, cX, cY);
  //  printf("DIst %f %f\n", dist, R);
  
  if (dist == R)
    {
      // line segment touches circle; one intersection point
      iX = nearestX;
      iY = nearestY;

      if (t < 0 || t > 1)
        {
          // intersection point is not actually within line segment
          return 0;
        }else
        {
          cart2angle(iX, iY, &a);
          if (a >= astart && a <= aend){                      
            return 1;
          }
        }
    }
  else if (dist < R)
    {
      // two possible intersection points

      dt = sqrt(R * R - dist * dist) / sqrt(dl);

      // intersection point nearest to A
      t1 = t - dt;
      iX = aX + t1 * dX;
      iY = aY + t1 * dY;
      if (t1 < 0 || t1 > 1)
        {

          // intersection point is not actually within line segment              
        }else{

          cart2angle(iX, iY, &a);

          if (a >= astart && a <= aend){
            return 1;
          }
      }

      // intersection point farthest from A
      t2 = t + dt;
      iX = aX + t2 * dX;
      iY = aY + t2 * dY;
      if (t2 < 0 || t2 > 1)
        {
          // intersection point is not actually within line segment          
        } else {
        cart2angle(iX, iY, &a);
        if (a >= astart && a <= aend){                      
          return 1;
        }

      }        
    } else{
    return 0;
  }

}
