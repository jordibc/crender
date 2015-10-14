import cPickle
import gzip
import sys
import time
import math

from ctypes import *
from numpy.ctypeslib import ndpointer
import numpy as np

from PyQt4 import QtCore
from PyQt4.QtGui import *
from PyQt4.QtCore import *

R90 = math.pi/2.0
R180 = math.pi
R270 = 3 * R90
R360 = 2 * math.pi 

def get_arc_path(inner_r, outter_r, arc_start, arc_end, cx, cy):
    "Returns a Qt Path representing an Annulus"
    arc_span = arc_end - arc_start
    arc_start, arc_end, arc_span = map(math.degrees, (arc_start, arc_end, arc_span))
    
    inner_diam = inner_r * 2.0
    outter_diam = outter_r * 2.0
    path = QPainterPath()
    rect1 = QRectF(cx - inner_r, cy - inner_r, inner_diam, inner_diam)
    rect2 = QRectF(cx - outter_r, cy - outter_r, outter_diam, outter_diam)
    path.arcMoveTo(rect1, -arc_start)
    i2 = path.currentPosition()
    path.arcMoveTo(rect1, -arc_end)
    o1 = path.currentPosition()
    path.moveTo(i2)
    path.arcTo(rect1, -arc_start, -arc_span)
    path.lineTo(o1)
    path.arcTo(rect2, -arc_end, arc_span)
    path.closeSubpath()
    return path


def check_intersect(sample_paths):
    """Check what annulus intersect a giver region of the scene represented by a rectangle"""

    # /// each sample is a list with the following values in the following order
    # ///
    # /// rx = rectangle top-left X
    # /// ry = rectangle top-left Y
    # /// rw = rectangle width 
    # /// rh = rectangle height
    # /// r1 = annulus inner radius
    # /// r2 = annulus outer radius
    # /// a1 = annulus arc angle start radians
    # /// a2 = annulus arc angle end in radians 
    # /// cx = circle center X for the referred annulus
    # /// cy = circle center Y for the referred annulus

    # Register C method
    crender = cdll.LoadLibrary("crender.so")    
    C_intersects = crender.RectAnnulusIntersect
    C_intersects.restype = c_int
    C_intersects.argtypes = [c_double, c_double,
                            c_double, c_double,
                            c_double, c_double,
                            c_double, c_double,
                            c_double, c_double]

    qt, native = [], [] 

    qt_paths = []
    for rx, ry, rw, rh, r1, r2, a1, a2, cx, cy in sample_paths:
        annulus = get_arc_path(r1, r2, a1, a2, cx, cy)
        rect = QRectF(rx, ry, rw, rh)
        qt_paths.append([annulus, rect])
        
    t1 = time.time()
    for annulus, rect in qt_paths:
        s = annulus.intersects(rect)
        qt.append(bool(s))
    print "time qt", time.time()-t1

    t1 = time.time()
    for rx, ry, rw, rh, r1, r2, a1, a2, cx, cy in sample_paths:
        s = C_intersects(r1, r2, a1, a2, rx, ry, rw, rh, cx, cy)
        native.append(bool(s))
    print "time native", time.time()-t1
    
    print 'total:', len(qt)
    print "correct:", len([1 for a in (np.array(qt) == np.array(native)) if a == True])
    print "fail:", len([1 for a in (np.array(qt) == np.array(native)) if a == False ])
    print qt
    print native
    
if __name__ == "__main__":
    samples = [
        #1
        (0, 0, 5, 5, 1, 20, 0, R270, 10, 10), # rect vertex within annulus
        #2
        (0, 0, 30, 30, 1, 10, 0, R90/2, 10, 10), # annulus within rect
        #3
        (4, 0, 3, 40, 1, 8, R180, R270, 10, 10), # rect side intersects annulus but vertexes are outside
        #4
        (0, 0, 60, 4, 1, 8, R180, R360, 10, 10), # rect side subtly intersects annulus outer arc
        #5
        (0, 0, 60, 4, 1, 5, R180, R360, 10, 10), # rect side subtly intersects annulus outer arc        
    ]
    check_intersect(samples)
