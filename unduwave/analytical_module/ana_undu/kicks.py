"""
Contains the functionality for loading and processing b-field data
"""
from unduwave.unduwave_incl import *
import unduwave.helpers.file_folder_helpers as ff_h
# import unduwave.quantities.quantities as quantities
from unduwave.attribute_classes.attributes import _attribute
from unduwave.attribute_classes.attributes import _attribute_collection
from unduwave.analytical_module.ana_undu.grids import grid_interpolator

import multiprocessing

def calculate_Be2_intgrl(self,
		gamma,
		nhor=None,
		dhor=None,
		horPnts=None,
		horIntrvl=None,
		nvert=None,
		dvert=None,
		vertPnts=None,
		vertIntrvl=None,
		longPnts=None,
		nlong=None,
		dlong=None,
		longIntrvl=None,
	) :

	alphaFac=q_el/(gamma*m_el*v_c)
	facKicks=-1e-6*1/2*(v_c*1e-9)**2
	facPot=1e-3

	vertPnts=self.createPointList(
		vals=self._bmap.yvals,npnts=nvert,dpnts=dvert,pnts=vertPnts,intrvl=vertIntrvl)
	horPnts=self.createPointList(
		vals=self._bmap.zvals,npnts=nhor,dpnts=dhor,pnts=horPnts,intrvl=horIntrvl)
	longPnts=self.createPointList(
		vals=self._bmap.xvals,npnts=nlong,dpnts=dlong,pnts=longPnts,intrvl=longIntrvl)
	if longPnts is None :
		# longIntrvl=[self._bmap.xvals[0],self._bmap.xvals[-1]]
		longPnts=self._bmap.xvals

	xvals=longPnts

	coordLimits=[
		[longPnts[0],longPnts[-1]],
		[vertPnts[0],vertPnts[-1]],
		[horPnts[0],horPnts[-1]],
		]
	YY, ZZ = np.meshgrid(vertPnts,horPnts,indexing='ij')
	screen0=np.zeros((len(vertPnts),len(horPnts),3))
	screen0[:,:,0] = xvals[0]
	screen0[:,:,1] = YY
	screen0[:,:,2] = ZZ
	pntList=screen0.reshape( (screen0.shape[0]*screen0.shape[1],3) )

	potFun=np.zeros( (len(xvals),len(pntList),1) )

	def run_integral(p0) :
		pnts, ds, intgrls=self.integrate_along_line(
			p0=p0,
			n0=np.array([1.0,0.0,0.0]),
			limit=400,
			epsrel=1e-3,
			epsabs=1e-5,
			which='ByBz',
			square=True,
			coordLimits=coordLimits
			)
		return intgrls

	res=map(run_integral,pntList)

	for ind,splineVec in enumerate(res) : 
		vertSpline=splineVec[0]
		horSpline=splineVec[1]
		potFun[:,ind,0]=(vertSpline(xvals)+horSpline(xvals))*facPot

	integralsByBzSq=potFun.reshape( (len(xvals),screen0.shape[0],screen0.shape[1],1) )
	potentialGrid=grid_interpolator(
		funs=integralsByBzSq,
		xvals=xvals,
		yvals=vertPnts,
		zvals=horPnts,
		)
	return potentialGrid

def calculate_kick_map(self,
		gamma,
		nhor=None,
		dhor=None,
		horPnts=None,
		horIntrvl=None,
		nvert=None,
		dvert=None,
		vertPnts=None,
		vertIntrvl=None,
		longPnts=None,
		nlong=None,
		dlong=None,
		longIntrvl=None,
	) :

	alphaFac=q_el/(gamma*m_el*v_c)
	facKicks=-1e-6*1/2*(v_c*1e-9)**2
	facPot=1e-3

	vertPnts=self.createPointList(
		vals=self._bmap.yvals,npnts=nvert,dpnts=dvert,pnts=vertPnts,intrvl=vertIntrvl)
	horPnts=self.createPointList(
		vals=self._bmap.zvals,npnts=nhor,dpnts=dhor,pnts=horPnts,intrvl=horIntrvl)
	longPnts=self.createPointList(
		vals=self._bmap.xvals,npnts=nlong,dpnts=dlong,pnts=longPnts,intrvl=longIntrvl)
	if longPnts is None :
		# longIntrvl=[self._bmap.xvals[0],self._bmap.xvals[-1]]
		longPnts=self._bmap.xvals

	xvals=longPnts

	coordLimits=[
		[longPnts[0],longPnts[-1]],
		[vertPnts[0],vertPnts[-1]],
		[horPnts[0],horPnts[-1]],
		]
	YY, ZZ = np.meshgrid(vertPnts,horPnts,indexing='ij')
	screen0=np.zeros((len(vertPnts),len(horPnts),3))
	screen0[:,:,0] = xvals[0]
	screen0[:,:,1] = YY
	screen0[:,:,2] = ZZ
	pntList=screen0.reshape( (screen0.shape[0]*screen0.shape[1],3) )

	potFun=np.zeros( (len(xvals),len(pntList),2) )

	def run_integral(p0) :
		pnts, ds, intgrls=self.integrate_along_line(
			p0=p0,
			n0=np.array([1.0,0.0,0.0]),
			limit=400,
			epsrel=1e-3,
			epsabs=1e-5,
			square=False,
			coordLimits=coordLimits
			)
		return intgrls

	res=map(run_integral,pntList)

	for ind,splineVec in enumerate(res) : 
		vertSpline=splineVec[0]
		horSpline=splineVec[1]
		potFun[:,ind,0]=vertSpline(xvals)
		potFun[:,ind,1]=horSpline(xvals)

	integralsByBzSq=potFun.reshape( (len(xvals),screen0.shape[0],screen0.shape[1],2) )
	integralsByBzSq=integralsByBzSq**2
	integralSumByBzSq=np.sum(integralsByBzSq,axis=3)
	integralSumByBzSq_t=np.zeros( (integralSumByBzSq.shape[0],integralSumByBzSq.shape[1],integralSumByBzSq.shape[2],1) )
	integralSumByBzSq_t[:,:,:,0]=integralSumByBzSq
	potentialGrid=grid_interpolator(
		funs=integralSumByBzSq_t,
		xvals=xvals,
		yvals=vertPnts,
		zvals=horPnts,
		)

	"""
	################################
	"""

	# def run_integral_focus_pot(p0) :

	# 	pntsIntPot, dsPot, intgrlPot=potentialGrid.integrate_along_line(
	# 		p0=p0,
	# 		n0=np.array([1.0,0.0,0.0]),
	# 		limit=400,
	# 		epsrel=1e-3,
	# 		epsabs=1e-5,
	# 		coordLimits=coordLimits,
	# 		)
	# 	return np.array([intgrlPot[-1]])

	# res=map(run_integral_focus_pot,pntList)

	# focusPotFun=np.zeros( (len(xvals),len(pntList),1) )
	# for ind,splineVec in enumerate(res) : 
	# 	potSpline=splineVec[0]
	# 	focusPotFun[:,ind,0]=potSpline(xvals)

	# focusPotFun=focusPotFun.reshape((len(xvals),screen0.shape[0],screen0.shape[1],1))		
	# focusPotentialGrid=grid_interpolator(
	# 	funs=focusPotFun,
	# 	xvals=xvals,
	# 	yvals=vertPnts,
	# 	zvals=horPnts,
	# 	)

	# gX,gY,gZ=focusPotentialGrid.calcGradient()

	# kickFun=np.zeros( (len(xvals),screen0.shape[0],screen0.shape[1],2) )

	# kickFun[:,:,:,0] = (gY.funs*facKicks).reshape( (len(xvals),screen0.shape[0],screen0.shape[1]) )
	# kickFun[:,:,:,1] = (gZ.funs*facKicks).reshape( (len(xvals),screen0.shape[0],screen0.shape[1]) )
	# kickMap=grid_interpolator(
	# 	funs=kickFun,
	# 	xvals=xvals,
	# 	yvals=vertPnts,
	# 	zvals=horPnts,
	# 	)
	# return focusPotentialGrid, kickMap

	gX,gY,gZ=potentialGrid.calcGradient()
	def run_integral_kick(p0) :
		pntsIntY, dsY, intgrlY=gY.integrate_along_line(
			p0=p0,
			n0=np.array([1.0,0.0,0.0]),
			limit=400,
			epsrel=1e-3,
			epsabs=1e-5,
			coordLimits=coordLimits,
			)
		pntsIntZ, dsZ, intgrlZ=gZ.integrate_along_line(
			p0=p0,
			n0=np.array([1.0,0.0,0.0]),
			limit=400,
			epsrel=1e-3,
			epsabs=1e-5,
			coordLimits=coordLimits,
			)
		pntsIntPot, dsP, intgrlP=potentialGrid.integrate_along_line(
			p0=p0,
			n0=np.array([1.0,0.0,0.0]),
			limit=400,
			epsrel=1e-3,
			epsabs=1e-5,
			coordLimits=coordLimits,
			)
		return np.array([intgrlY[-1],intgrlZ[-1],intgrlP[-1]])

	res=map(run_integral_kick,pntList)
	kickFun=np.zeros( (1,len(pntList),2) )
	potFun=np.zeros( (1,len(pntList),1) )
	for ind,splineVec in enumerate(res) : 
		vertSpline=splineVec[0]
		horSpline=splineVec[1]
		potSpline=splineVec[2]
		kickFun[0,ind,0]=vertSpline(xvals[-1])*facKicks
		kickFun[0,ind,1]=horSpline(xvals[-1])*facKicks
		potFun[0,ind,0]=potSpline(xvals[-1])*facPot

	kickFun=kickFun.reshape((1,screen0.shape[0],screen0.shape[1],2))
	kickMap=grid_interpolator(
		funs=kickFun,
		xvals=[xvals[-1]],
		yvals=vertPnts,
		zvals=horPnts,
		)

	potFun=potFun.reshape((1,screen0.shape[0],screen0.shape[1],1))
	potMap=grid_interpolator(
		funs=potFun,
		xvals=[xvals[-1]],
		yvals=vertPnts,
		zvals=horPnts,
		)
	return potentialGrid, potMap, kickMap

def integrate_along_line(
		vals,
		limit=400,
		epsrel=1e-3,
		epsabs=1e-5,
	) : 

	nvals=3
	ds=[i for i in range(len(vals))]
	firstInts=[]
	for i in range(nvals) :
		cs_integr = CubicSpline(ds, vals[:,i])
		firstInt=np.zeros((len(ds),))
		firstInt[0]=0.0
		for ind, x in enumerate(ds[0:-1]) :
			integral=quad( cs_integr, x, ds[ind+1], 
				limit = limit,
				epsrel=epsrel,
				epsabs=epsabs,
				)[0]
			firstInt[ind+1] = integral
			firstInt[ind+1]=firstInt[ind]+firstInt[ind+1]
		firstInts.append(CubicSpline(ds,firstInt))

	return firstInts

def _init_shared(shared_data) :
	global SHARED
	SHARED=shared_data

def worker(arg):
	myData=SHARED['data']
	myVals=myData[:,arg[0],arg[1]]
	integral=integrate_along_line(
		vals=myVals,
		limit=400,
		epsrel=1e-3,
		epsabs=1e-5,
	)
	return integral

def calc_beff_mp(inputs, shared_data=None, processes=4):
	# Must guard pool creation when module may be imported on Windows
	with multiprocessing.Pool(processes=processes, initializer=_init_shared, initargs=(shared_data,)) as p:
		return p.map(worker, inputs)