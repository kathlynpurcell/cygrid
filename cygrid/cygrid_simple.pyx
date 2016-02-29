import numpy as np
import traceback,sys
cimport numpy as np
from kapteyn import wcs
cimport cython
from cython.parallel import prange
from cython cimport view

from libc.math cimport exp, cos, sin, sqrt, asin


def eprint():
    print sys.exc_info()[1]
    print traceback.print_tb(sys.exc_info()[2])
#

# Build command
#/vol/software/software/tools/epd/epd-x86_64/bin/python setup.py build_ext --inplace

cdef class pygrid:
    
    cdef np.ndarray datacube, weights, xwcs, ywcs
    cdef object header, proj, proj2d, yxshape
    cdef int zsize, ysize, xsize
    
    def __init__(self, header, **kwargs):
        
        self.header = header
        self.proj = wcs.Projection(header)
        self.proj2d  = self.proj.sub(axes=[self.proj.lonaxnum, self.proj.lataxnum] )
        
        self.datacube = np.zeros((header['NAXIS3'], header['NAXIS2'], header['NAXIS1']), dtype = np.float32)
        self.weights = self.datacube.copy()
        print "self.datacube",(header['NAXIS3'], header['NAXIS2'], header['NAXIS1'])
        
        self.zsize = header['NAXIS3']
        self.ysize = header['NAXIS2']
        self.xsize = header['NAXIS1']
        self.yxshape = (self.ysize,self.xsize)
        self.ywcs, self.xwcs = np.indices(self.yxshape)
        #zdummy = np.ones_like(self.ywcs)
        try:
            #self.xwcs, self.ywcs, zdummy = self.proj2d.toworld((self.xwcs.flatten()+1,self.ywcs.flatten()+1,zdummy.flatten()))
            self.xwcs, self.ywcs = self.proj2d.toworld((self.xwcs.flatten()+1,self.ywcs.flatten()+1))
        except:
            # this will set invalid coords to nan, in case, something
            # was wrong with coord trafo
            print "something was wrong with coord trafo"
            print "will try to set wrong coords to nan"
            eprint()
            #self.xwcs, self.ywcs, zdummy = self.proj2d.toworld()
            self.xwcs, self.ywcs = self.proj2d.toworld()
        self.xwcs = self.xwcs.reshape(self.yxshape).astype(np.float64)
        self.ywcs = self.ywcs.reshape(self.yxshape).astype(np.float64)
    
    @cython.boundscheck(False)  
    @cython.wraparound(False)
    @cython.cdivision(True)
    def grid(self, np.ndarray[double, ndim = 1] lons, np.ndarray[double, ndim = 1] lats, float[:,:] fluxes, float[:,:] weights, float kr, float ksq):
        print "Gridding",len(fluxes),"spectra in datacube..."
        
        cdef double pi = 3.141592653589793
        cdef double degtorad = pi/180.
        
        cdef int z,y,x,i
        cdef int speccount = len(fluxes)
        cdef double sdist, sweight, lonmin, lonmax, latmin, latmax, l1, l2, b1, b2, sinbdiff, sinldiff
        cdef double mynan = np.nan
        
        cdef float[:,:,:] datacubeview = self.datacube
        cdef float[:,:,:] weightsview = self.weights
        cdef double[:,:] ywcsview = self.ywcs
        cdef double[:,:] xwcsview = self.xwcs
        
        cdef double kernelradius = kr
        cdef double kernelsigmasq = ksq
        
        lonmin = np.min(lons)
        lonmax = np.max(lons)
        latmin = np.min(lats)
        latmax = np.max(lats)
        print "old lonmin,lonmax,latmin,latmax: %.3f %.3f %.3f %.3f"%(lonmin,lonmax,latmin,latmax)
        #isInField=lambda x,y: ((xwcsview[y,x] + kernelradius) > lonmin) & ((xwcsview[y,x] - kernelradius) < lonmax) & ((ywcsview[y,x] + kernelradius) > latmin) & ((ywcsview[y,x] - kernelradius) < latmax)
        zeroLongField=False
        # testing for Ra=0 problem
        lons2=lons.copy()
        lons2[lons2<180.]+=360.
        if abs(np.max(lons2)-np.min(lons2))<abs(lonmax-lonmin)-1e-4:
            lonmin,lonmax = np.min(lons2),np.max(lons2)
            print "new lonmin,lonmax,latmin,latmax: %.3f %.3f %.3f %.3f"%(lonmin,lonmax,latmin,latmax)
            zeroLongField=True
            #isInField=lambda x,y: (((xwcsview[y,x] + kernelradius) >0 ) & ((xwcsview[y,x] - kernelradius) < lonmax-360.) | ((xwcsview[y,x] + kernelradius) > lonmin) & ((xwcsview[y,x] - kernelradius) < 360.)) & ((ywcsview[y,x] + kernelradius) > latmin) & ((ywcsview[y,x] - kernelradius) < latmax)
        
        for y in prange(self.ysize, nogil = True, schedule = 'dynamic'):
            for x in range(self.xsize):
                #if (mynan==xwcsview[y,x]): continue
                #if (mynan==ywcsview[y,x]): continue
                if (xwcsview[y,x]!=xwcsview[y,x]): continue
                if (ywcsview[y,x]!=ywcsview[y,x]): continue
                if ( ((not zeroLongField) & (xwcsview[y,x]+kernelradius>lonmin) & (xwcsview[y,x]-kernelradius<lonmax)) | (zeroLongField & (((xwcsview[y,x]+kernelradius>0 ) & (xwcsview[y,x]-kernelradius<lonmax-360.)) | ((xwcsview[y,x]+kernelradius>lonmin) & (xwcsview[y,x]-kernelradius<360.)))) ) & (ywcsview[y,x]+kernelradius>latmin) & (ywcsview[y,x]-kernelradius<latmax):
                    
                    l1 = xwcsview[y,x] * degtorad
                    b1 = ywcsview[y,x] * degtorad
                    
                    for i in range(speccount):
                        
                        l2 = lons[i] * degtorad
                        b2 = lats[i] * degtorad
                        
                        sinbdiff = sin( (b1 - b2)/ 2.)
                        sinldiff = sin( (l1 - l2)/ 2.)
                        
                        sdist = 2. * asin( sqrt( sinbdiff * sinbdiff + cos(b1) * cos(b2) * sinldiff * sinldiff) ) / degtorad
                        
                        if sdist < kernelradius:
                            sweight = exp(sdist*sdist/-2./kernelsigmasq)
                            for z in range(self.zsize):
                                datacubeview[z,y,x] += fluxes[i,z] * weights[i,z] * sweight
                                weightsview[z,y,x] += weights[i,z] * sweight
    
    def getDatacube(self):
        return self.datacube/self.weights
        
    def getWeights(self):
        return self.weights
        
    def getProjection(self):
        return self.proj
        
    def getProjection2D(self):
        return self.proj2d
        
    def setXWCS(self, np.ndarray[double, ndim = 2] wcs):
        self.xwcs = wcs
        
    def setYWCS(self, np.ndarray[double, ndim = 2] wcs):
        self.ywcs = wcs




























