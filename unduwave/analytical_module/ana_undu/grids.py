"""
Contains the functionality for loading and processing b-field data
"""
from unduwave.unduwave_incl import *
import unduwave.helpers.file_folder_helpers as ff_h
from unduwave.constants import *
# import unduwave.quantities.quantities as quantities
from unduwave.attribute_classes.attributes import _attribute
from unduwave.attribute_classes.attributes import _attribute_collection
import unduwave.analytical_module.ana_undu.bfield as bfield

def get_line_through_point(p0,n0,xlims,ylims,zlims) :

	norm=np.linalg.norm(n0)
	n0=n0/norm
	xmax=xlims[-1]
	xmin=xlims[0]
	ymax=ylims[-1]
	ymin=ylims[0]
	zmax=zlims[-1]
	zmin=zlims[0]

	maxs=np.array([xmax,ymax,zmax])
	mins=np.array([xmin,ymin,zmin])
	with np.errstate(divide='ignore', invalid='ignore'):
		tsMaxT=(maxs-p0)/n0
		tsMinT=(mins-p0)/n0
	tsMax=[]
	tsMin=[]
	for t in tsMaxT : 
		if (not np.isinf(t)) and ( not np.isnan(t) ) :
			tsMax.append(t)
	for t in tsMinT : 
		if (not np.isinf(t)) and ( not np.isnan(t) ) :
			tsMin.append(t)
	minT=tsMin[np.argmin(np.abs(tsMin))]
	maxT=tsMax[np.argmin(np.abs(tsMax))]
	pMin=p0+n0*minT
	pMax=p0+n0*maxT

	return pMin, pMax

def grid_get_lines_init_shared(
		funs,
		xvals,
		yvals,
		zvals, 
		nvec,
		interpolators,
		nlinepnts,
		) :
	global SHARED
	SHARED={
		'funs':funs,
		'xvals' : xvals, 
		'yvals' : yvals, 
		'zvals' : zvals, 
		'nvec' : nvec,
		'interpolators' : interpolators,
		'nlinepnts' : nlinepnts,
		}

def grid_get_lines_worker(p0):
	funs=SHARED['funs']
	xvals=SHARED['xvals']
	yvals=SHARED['yvals']
	zvals=SHARED['zvals']
	nvec=SHARED['nvec']
	nlinepnts=SHARED['nlinepnts']
	interpolators=SHARED['interpolators']

	pMin, pMax = get_line_through_point(
		p0=p0,
		n0=nvec,
		xlims=[xvals[0],xvals[-1]],
		ylims=[yvals[0],yvals[-1]],
		zlims=[zvals[0],zvals[-1]],
		)

	n0N=nvec/np.linalg.norm(nvec)
	distP=np.linalg.norm(p0-pMin)
	tVal0=distP
	sVal0=np.dot(p0,n0N)
	ds0=sVal0-tVal0

	nTs=np.linspace(0,1,nlinepnts)
	pnts=np.zeros( (nlinepnts,3) )
	for indt, t in enumerate(nTs) :
		pnts[indt,:]=pMin+t*(pMax-pMin)

	vals=np.zeros( (nlinepnts,len(interpolators)) )

	for ind_interpol, interpolator in enumerate(interpolators) :
		valstmp=interpolator(pnts)
		vals[:,ind_interpol]=valstmp

	ds=np.zeros((pnts.shape[0],))
	for ind,pnt in enumerate(pnts[1:]) :
		diff=pnt-pnts[ind]
		diff=np.linalg.norm(diff)
		ds[ind+1]=ds[ind]+diff
	ds=ds+ds0

	return [pnts, vals, ds	]

	# return {'pnts':pnts, 'ds' : ds, 'firstInts' : firstInts}

def grid_get_lines_calc(
		funs, 
		xvals,
		yvals,
		zvals,
		p0s,
		nvec,
		interpolators,
		nlinepnts=100,
		processes=3,
		):

	initargs=(funs,xvals,yvals,zvals,nvec,interpolators,nlinepnts)
	# Must guard pool creation when module may be imported on Windows
	with multiprocessing.Pool(processes=processes, initializer=grid_get_lines_init_shared, initargs=initargs) as p:
		res_all=p.map(grid_get_lines_worker, p0s)

	return res_all

def grid_line_integrals_init_shared(
		funs,
		xvals,
		yvals,
		zvals, 
		limit=400,
		epsrel=1e-3,
		epsabs=1e-5,
		) :
	global SHARED
	SHARED={
		'funs':funs,
		'xvals' : xvals, 
		'yvals' : yvals, 
		'zvals' : zvals, 
		'limit' : limit, 
		'epsrel':epsrel,
		'epsabs':epsabs
		}

def grid_line_integrals_worker(cur_ind):
	funs=SHARED['funs']
	xvals=SHARED['xvals']
	yvals=SHARED['yvals']
	zvals=SHARED['zvals']
	limit=SHARED['limit']
	epsrel=SHARED['epsrel']
	epsabs=SHARED['epsabs']

	pnts=[ np.array([xval,yvals[cur_ind[0]],zvals[cur_ind[1]]]) for xval in xvals ]
	pnts=np.array(pnts)

	vals=funs[:,cur_ind[0],cur_ind[1],:]

	ds=np.zeros((pnts.shape[0],))
	for ind,pnt in enumerate(pnts[1:]) :
		diff=pnt-pnts[ind]
		diff=np.linalg.norm(diff)
		ds[ind+1]=ds[ind]+diff
	ds=ds+xvals[0]

	firstInts=[]
	for i in range(funs.shape[-1]) :
		valsIntgrt=vals[:,i]
		# if square :
		# 	valsIntgrt=np.square(valsIntgrt)
		cs_integr = CubicSpline(ds, valsIntgrt)
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

	return np.array(firstInts)

def grid_line_integrals_calc(
		funs, 
		xvals,
		yvals,
		zvals,
		limit=400,
		epsrel=1e-3,
		epsabs=1e-5,
		processes=3,
		):

	shape=funs.shape
	allindices=[]
	for indy in range(shape[1]) :
		for indz in range(shape[2]) :
			allindices.append([indy,indz])

	initargs=(funs,xvals,yvals,zvals,limit,epsrel,epsabs)
	# Must guard pool creation when module may be imported on Windows
	with multiprocessing.Pool(processes=processes, initializer=grid_line_integrals_init_shared, initargs=initargs) as p:
		res_all=p.map(grid_line_integrals_worker, allindices)
	nres_all=np.array(res_all)
	# nres_all=nres_all.reshape( (1,shape[1],shape[2],shape[3]) )

	return nres_all

def grid_pnts_init_shared(pnts,interpolators) :
	global SHARED
	SHARED={'pnts':pnts, 'interpolators' : interpolators}

def grid_pnts_worker(cur_ind):
	pnts=SHARED['pnts']
	interpolators=SHARED['interpolators']
	res=interpolators[cur_ind](np.round(pnts,5))
	return res

def grid_pnts_calc(pnts, interpolators, processes=3):
	nvals=len(interpolators)
	# Must guard pool creation when module may be imported on Windows
	with multiprocessing.Pool(processes=processes, initializer=grid_pnts_init_shared, initargs=(pnts,interpolators,)) as p:
		res_all=p.map(grid_pnts_worker, [*range(nvals)])

	res=np.zeros( (len(pnts),nvals) )
	for i in range(nvals) :
		res[:,i]=res_all[i]

	return res

class unduGridError(Exception) :
	pass

def createPointList(vals,npnts=None,dpnts=None,pnts=None,intrvl=None) :
	try:
		if npnts is None :
			npnts=len(vals)

		if not (pnts is None) :
			if True in (pnts > vals[-1]) :
				raise unduGridError("bmap_interpolation: createPointList: pnts outside of limit.")
			if True in (pnts < vals[0]) :
				raise unduGridError("bmap_interpolation: createPointList: pnts outside of limit.")
		elif not (dpnts is None) :
			npnts=int((vals[-1]-vals[0])/dpnts)
			if npnts < 1 :
				npnts=1
			pnts=np.linspace(vals[0],vals[-1],npnts)
		elif not (intrvl is None) :
			if intrvl[-1] > vals[-1] :
				intrvl[-1]=vals[-1]
			if intrvl[0] < vals[0] :
				intrvl[0]=vals[0]
			pnts = np.linspace(intrvl[0],intrvl[-1],npnts)
		else :
			if npnts == 1 :
				pnts = [(vals[0]+vals[-1])/2]
			else :
				pnts = np.linspace(vals[0],vals[-1],npnts)
	except unduGridError as e :
		print(e)

	return pnts

class grid_interpolator : 

	def __init__(self,
		funs,
		xvals,
		yvals=None,
		zvals=None,
		) : 
		self._g_xvals=np.array(xvals)
		self._g_yvals=np.array(yvals)
		self._g_zvals=np.array(zvals)
		self.setfuns(funs=funs)
		self.set_limits()
		self.processes=6

	def setfuns(self,funs) :
		shapeFuns=funs.shape
		self._g_nvals=1
		if len(shapeFuns) > 3 :
			self._g_nvals=shapeFuns[3]
		self._g_funs=funs
		self.create_grid_interpolation()

	def create_difference_grid(self,grid,nx=None,ny=None,nz=None) :

		xvals1=self._g_xvals
		xvals2=grid._g_xvals
		yvals1=self._g_yvals
		yvals2=grid._g_yvals
		zvals1=self._g_zvals
		zvals2=grid._g_zvals		

		nx2=[]
		for xval in xvals1 : 
			if (xval >= xvals2[0]) and (xval <= xvals2[-1]) :
				nx2.append(xval)
		xvals=[]
		for xval in xvals2 : 
			if (xval >= nx2[0]) and (xval <= nx2[-1]) :
				xvals.append(xval)
		ny2=[]
		for yval in yvals1 : 
			if (yval >= yvals2[0]) and (yval <= yvals2[-1]) :
				ny2.append(yval)
		yvals=[]
		for yval in yvals2 : 
			if (yval >= ny2[0]) and (yval <= ny2[-1]) :
				yvals.append(yval)
		nz2=[]
		for zval in zvals1 : 
			if (zval >= zvals2[0]) and (zval <= zvals2[-1]) :
				nz2.append(zval)
		zvals=[]
		for zval in zvals2 : 
			if (zval >= nz2[0]) and (zval <= nz2[-1]) :
				zvals.append(zval)


		myfuns, xs, ys, zs=self.get_values_on_grid(
			xs=xvals,
			ys=yvals,
			zs=zvals,
			nx=nx,
			ny=ny,
			nz=nz
			)

		otherfuns, xs, ys, zs=grid.get_values_on_grid(
			xs=xvals,
			ys=yvals,
			zs=zvals,
			nx=nx,
			ny=ny,
			nz=nz
			)
		myfuns=myfuns.reshape((len(xs),len(ys),len(zs),myfuns.shape[-1]))
		otherfuns=otherfuns.reshape((len(xs),len(ys),len(zs),myfuns.shape[-1]))
		new_grid=np.zeros( myfuns.shape )
		new_grid=np.abs(myfuns-otherfuns)

		diff_drig=grid_interpolator(
			funs=new_grid,
			xvals=xs,
			yvals=ys,
			zvals=zs,
			)
		return diff_drig

	def plot_grid_map(self,
			indPlot=0,
			zlab='',
			xPos=None,
			yPos=None,
			zPos=None,
			nHor=None,
			nVert=None,
			nfig=None,
			save=False,
			filename=None,
			title=None,
			) : 
		try:
			bInd=indPlot
			notNone=0
			if not (xPos is None) :

				if nHor is None :
					nHor=len(self._g_zvals)
				if nVert is None :
					nVert=len(self._g_yvals)

				pnts, dsY, dsZ=self.get_values_on_plane(
					p0=np.array([xPos,0.0,0.0]),
					v0=np.array([0.0,1.0,0.0]),
					h0=np.array([0.0,0.0,1.0]),
					nh=nHor,
					nv=nVert,
					)
				xPlotVals=dsZ
				yPlotVals=dsY
				bPlotVals=pnts[:,:,bInd]
				bPlotVals=bPlotVals.T
				xLab="z [mm]"
				yLab="y [mm]"
				notNone=notNone+1
				pos=f"x={xPos:.2f} mm"
			elif not (yPos is None) :

				if nHor is None :
					nHor=len(self._g_xvals)
				if nVert is None :
					nVert=len(self._g_zvals)
				pnts, dsZ, dsX=self.get_values_on_plane(
					p0=np.array([0.0,yPos,0.0]),
					h0=np.array([1.0,0.0,0.0]),
					v0=np.array([0.0,0.0,1.0]),
					nh=nHor,
					nv=nVert,
					)
				xPlotVals=dsX
				yPlotVals=dsZ
				bPlotVals=pnts[:,:,bInd]
				bPlotVals=bPlotVals.T
				xLab="x [mm]"
				yLab="z [mm]"
				notNone=notNone+1
				pos=f"y={yPos:.2f} mm"

				# if filename.find("undu_By_interpol") >= 0:
				# 	pdb.set_trace()
			elif not (zPos is None) :

				if nHor is None :
					nHor=len(self._g_xvals)
				if nVert is None :
					nVert=len(self._g_yvals)

				pnts, dsY, dsX=self.get_values_on_plane(
						p0=np.array([0.0,0.0,zPos]),
					h0=np.array([1.0,0.0,0.0]),
					v0=np.array([0.0,1.0,0.0]),
					nh=nHor,
					nv=nVert,
					)
				xPlotVals=dsX
				yPlotVals=dsY
				bPlotVals=pnts[:,:,bInd]
				bPlotVals=bPlotVals.T
				xLab="x [mm]"
				yLab="y [mm]"
				notNone=notNone+1
				pos=f"z={zPos:.2f} mm"
			if (notNone < 1) or (notNone>1) :
				raise unduBfieldError("grid_interpolator.plot_grid_map: Only one of xPos, yPos or zPos must be given.")
			if title is None :
				title=f"B-Field Map at\n {pos}"

			X_data,Z_data = np.meshgrid(xPlotVals,yPlotVals,indexing='ij')

			clf=False
			if nfig is None :
				clf=True
				nfig=0

			# if title.find("Radia Kickmap") >= 0 :
			# 	pdb.set_trace()
			fig = plt.figure(nfig,figsize=(2*13*cm_inch, 2*6.5*cm_inch), dpi=150)
			if clf:
				plt.clf()
			fig.suptitle(title, fontsize=12)

			ax = plt.gca()

			plt.tight_layout()

			cmap = plt.colormaps["plasma"]
			cmap = cmap.with_extremes(bad=cmap(0))

			pcm = ax.pcolormesh(X_data,Z_data,bPlotVals, cmap=cmap)
			cb=fig.colorbar(pcm, ax=ax, label=zlab)                    

			axCb = cb.ax
			text = axCb.yaxis.label
			axCb.yaxis.get_offset_text().set_fontsize(8)
			axCb.tick_params(axis='both', which='major', labelsize=8)
			font = matplotlib.font_manager.FontProperties(size=8)
			text.set_font_properties(font)

			ax.set_ylabel(yLab, fontsize=8)
			ax.set_xlabel(xLab, fontsize=8)
			ax.tick_params(axis='both', which='major', labelsize=8)

			if save :
				if filename is None :
					filename = "colormap.png"
				plt.savefig(filename , bbox_inches='tight')

			plt.draw()
			plt.ion()
		except unduGridError as e :
			print(e)
			pdb.set_trace()

	def write_grid_data(self,file,cols=None) :
		path, fileN = os.path.split(file)
		os.makedirs(path, exist_ok=True)

		x_vals=self._g_xvals
		y_vals=self._g_yvals
		z_vals=self._g_zvals
		bvals = self._g_funs
		shapeFuns=self._g_nvals
		if cols is None :
			cols = ['x','y','z']
			for i in range(shapeFuns) :
				cols.append(f'f_{i}')

		with open( file, 'w') as o_f:
			for col in cols :
				o_f.write(f"{col} ")
			o_f.write(f'\n')
			for indx, xval in enumerate(np.unique(x_vals)) :
				for indy, yval in enumerate(np.unique(y_vals)) :
					for indz, zval in enumerate(np.unique(z_vals)) :
						o_f.write(f"{xval}  {yval}  {zval}")		
						for i in range(shapeFuns) :
							o_f.write(f"  {bvals[indx,indy,indz,i]}")		
						o_f.write(f"\n")

	@staticmethod
	def load_radia_kicks(file,xval=0.0) :

		try:
			with open(file, 'r') as o_f:
				lines = o_f.readlines()

			horKickSt=None
			vertKickSt=None
			potSt=None
			for indl,line in enumerate(lines) : 
				if line.find("# Horizontal Kick [T2m2]")>=0 : 
					horKickSt=indl+2
				if line.find("# Vertical Kick [T2m2]")>=0 : 
					vertKickSt=indl+2
				if line.find("# Longitudinally Integrated Squared Transverse Magnetic Field [T2m]")>=0 : 
					potSt=indl+2
			if horKickSt is None :
				raise unduGridError("grid_interpolator: load_radia_kicks: horizontal kicks not found")
			if vertKickSt is None :
				raise unduGridError("grid_interpolator: load_radia_kicks: vertical kicks not found")
			if potSt is None :
				raise unduGridError("grid_interpolator: load_radia_kicks: potential not found")

			yvals=[]
			zvals=[]
			zvalstmp=lines[horKickSt]
			zvalstmp=zvalstmp.strip().split(' ')
			for val in zvalstmp :
				if len(val) > 0 :
					zvals.append(float(val))
			zvals=np.array(zvals)

			#load hor kicks
			horKicks=[]
			indS=horKickSt+1
			while True :
				if lines[indS].find("# Vertical Kick [T2m2]") >= 0 :
					break
				if indS >= (len(lines)) : 
					break
				lineVals=[]
				lineValsTmp=lines[indS].strip().split(' ')
				yval=None
				for indv,val in enumerate(lineValsTmp) :
					if len(val) > 0 :
						if yval is None :
							yval=float(val)
						else :
							lineVals.append(float(val))
				yvals.append(yval)
				horKicks.append(lineVals)
				indS=indS+1
			horKicksTmp=np.array(horKicks)
			horKicks=np.zeros( (1,len(yvals),len(zvals),1) )
			horKicks[:,:,:,0]=horKicksTmp
			yvals=np.array(yvals)

			#load vert kicks
			vertKicks=[]
			indS=vertKickSt+1
			while True :
				if lines[indS].find("# Longitudinally Integrated Squared Transverse Magnetic Field [T2m]") >= 0 :
					break
				if indS >= (len(lines)) : 
					break
				lineVals=[]
				lineValsTmp=lines[indS].strip().split(' ')
				yval=None
				for indv,val in enumerate(lineValsTmp) :
					if len(val) > 0 :
						if yval is None :
							yval=float(val)
						else :
							lineVals.append(float(val))
				vertKicks.append(lineVals)
				indS=indS+1
			vertKicksTmp=np.array(vertKicks)
			vertKicks=np.zeros( (1,len(yvals),len(zvals),1) )
			vertKicks[:,:,:,0]=vertKicksTmp

			#load potential
			potential=[]
			indS=potSt+1
			while True :
				if indS >= (len(lines)) : 
					break
				lineVals=[]
				lineValsTmp=lines[indS].strip().split(' ')
				yval=None
				for indv,val in enumerate(lineValsTmp) :
					if len(val) > 0 :
						if yval is None :
							yval=float(val)
						else :
							lineVals.append(float(val))
				potential.append(lineVals)
				indS=indS+1
			potentialTmp=np.array(potential)
			potential=np.zeros( (1,len(yvals),len(zvals),1) )
			potential[:,:,:,0]=potentialTmp

			horKickGrid=grid_interpolator(
				funs=np.flip(horKicks,axis=1),
				xvals=np.array([xval]),
				yvals=np.flip(yvals)*1e3,
				zvals=zvals*1e3,
				)

			vertKickGrid=grid_interpolator(
				funs=np.flip(vertKicks,axis=1),
				xvals=np.array([xval]),
				yvals=np.flip(yvals)*1e3,
				zvals=zvals*1e3,
				)

			potentialGrid=grid_interpolator(
				funs=np.flip(potential,axis=1),
				xvals=np.array([xval]),
				yvals=np.flip(yvals)*1e3,
				zvals=zvals*1e3,
				)

		except unduGridError as e :
			print(e)
		return horKickGrid, vertKickGrid, potentialGrid

	@staticmethod
	def load_grid_data(
			file, 
			header=0,
			skiprows=None,
			) :
		"""
		"""

		data = pd.read_csv( file, skiprows=skiprows, dtype=object, delimiter=r"\s+",header=header)

		try:
			if (not ( 'x' in data.columns)) or (not ( 'y' in data.columns)) or (not ( 'z' in data.columns)) : 
				raise unduGridError("grid_interpolator: load_grid_data: not all coordinates in data")
		except unduGridError as e :
			print(e)

		not_xyz=[]
		for col in data.columns:
			if (not (col == 'x')) :
				if (not (col == 'y')) :
					if (not (col == 'z')) :
						not_xyz.append(col)
			data[col] = data[col].astype(float)

		unique_x=data['x'].unique()
		unique_y=data['y'].unique()
		unique_z=data['z'].unique()

		nxmap=len(unique_x)
		nymap=len(unique_y)
		nzmap=len(unique_z)
		funs=np.zeros((nxmap,nymap,nzmap,len(not_xyz)))
		for ind,col in enumerate(not_xyz) :
			dataC=data[col]
			dataC=np.array(dataC)
			dataC=dataC.reshape((nxmap,nymap,nzmap,))
			funs[:,:,:,ind]=dataC

		grid=grid_interpolator(
			funs=funs,
			xvals=unique_x,
			yvals=unique_y,
			zvals=unique_z,
			)
		return grid

	@staticmethod
	def crop_grid(data, xlims=None,ylims=None,zlims=None) :

		if xlims is None :
			indx0=0
			indxE=len(data._g_xvals)
		else :
			indx0=None
			indxE=None
			for indx,valx in enumerate(data._g_xvals) :
				if (indx0 is None ) and (valx>=xlims[0]) :
					indx0=indx
				if (indxE is None ) and (valx>xlims[1]) :
					indxE=indx
					break
			if indxE is None :
				indxE=len(xvals)
			elif indxE < 0 :
				indxE=0

		if ylims is None :
			indy0=0
			indyE=len(data._g_yvals)
		else :
			indy0=None
			indyE=None
			for indy,valy in enumerate(data._g_yvals) :
				if (indy0 is None ) and (valy>=ylims[0]) :
					indy0=indy
				if (indyE is None ) and (valy>ylims[1]) :
					indyE=indy
					break
			if indyE is None :
				indyE=len(yvals)
			elif indyE < 0 :
				indyE=0

		if zlims is None :
			indz0=0
			indzE=len(data._g_zvals)
		else :
			indz0=None
			indzE=None
			for indz,valz in enumerate(data._g_zvals) :
				if (indz0 is None ) and (valz>=zlims[0]) :
					indz0=indz
				if (indzE is None ) and (valz>zlims[1]) :
					indzE=indz
					break
			if indzE is None :
				indzE=len(zvals)

		nfuns=data._g_funs[indx0:indxE,indy0:indyE,indz0:indzE,:]

		newgrid=grid_interpolator(
			funs=nfuns,
			xvals=data._g_xvals[indx0:indxE],
			yvals=data._g_yvals[indy0:indyE],
			zvals=data._g_zvals[indz0:indzE],
			)
		return newgrid
			
	def set_limits(self) :
		self._g_xmin=round(min(self._g_xvals),5)
		self._g_xmax=round(max(self._g_xvals),5)
		self._g_ymin=round(min(self._g_yvals),5)
		self._g_ymax=round(max(self._g_yvals),5)
		self._g_zmin=round(min(self._g_zvals),5)
		self._g_zmax=round(max(self._g_zvals),5)

	def is_in_limits(self,pnt) : 
		if round(pnt[0],5) >= round(self._g_xmin,5) :
			if round(pnt[0],5) <= round(self._g_xmax,5) :
				if round(pnt[1],5) >= round(self._g_ymin,5) :
					if round(pnt[1],5) <= round(self._g_ymax,5) :
						if round(pnt[2],5) >= round(self._g_zmin,5) :
							if round(pnt[2],5) <= round(self._g_zmax,5) :
								return True
		return False

	def create_grid_interpolation(self) :
		"""
		create interpolation from bvals, xvals,...
		"""
		pnts=( np.round(self._g_xvals,6),np.round(self._g_yvals,6),np.round(self._g_zvals,6) )
		# self._interp = RegularGridInterpolator(pnts, self.funs)
		self.interps=[]
		for i in range(self._g_nvals) :
			self.interps.append(RegularGridInterpolator(pnts, self._g_funs[:,:,:,i]))

	def get_values(self,pnts) :

		try:
			for pnt in pnts :
				if not self.is_in_limits(pnt=pnt) :
					raise unduGridError(f"grid_interpolator: get_values: pnt {pnt} lies outside.")
		except unduGridError as e :
			print(e)
			return
		npnts=len(pnts)

		res=grid_pnts_calc(
			pnts=pnts, 
			interpolators=self.interps, 
			processes=self.processes,
			)

		# res=np.zeros( (npnts,self._g_nvals) )
		# for i in range(self._g_nvals) :
		# 	res[:,i]=self.interps[i](pnts)

		return res

	def __call__(self,pnts) :		
		return self.get_values(pnts=pnts)

	def get_values_on_plane(self,p0,v0,h0,nv=None,nh=None) : 
		"""
		Returns the values on a plane including point p0, with n0 and v0 spanning the plane
		"""
		if nv is None :
			nv=len(self.yvals)
		if nh is None :
			nh=len(self.zvals)
		pntsHor0, valsHor0, dsHor = self.get_values_along_line(p0=p0,n0=h0,n=nh)
		# shiftHor=p0-pntsHor0[0]
		# normHor=np.linalg.norm(shiftHor)
		# dsHor=dsHor-normHor
		multVals=False
		if len(valsHor0.shape) > 1 :
			multVals=True
			nvals=valsHor0.shape[-1]
			valsGrid=np.zeros( (nv,nh,nvals) )
		else :
			valsGrid=np.zeros( (nv,nh) )

		all_line_vals = grid_get_lines_calc(
				funs=self._g_funs, 
				xvals=self._g_xvals,
				yvals=self._g_yvals,
				zvals=self._g_zvals,
				p0s=pntsHor0,
				nvec=v0,
				nlinepnts=nv,
				interpolators=self.interps,
				processes=self.processes,
				)

		for indh, one_line in enumerate(all_line_vals) :
			[pntsVert, valVert, dsVert] = one_line
			# normVert=p0[1]-pntsVert[0][1]
			# dsVert=dsVert-normVert
			if multVals :
				for indV in range(nvals) : 
					valsGrid[:,indh,indV] = valVert[:,indV]
			else :
				valsGrid[:,indh] = valVert[:]

		return valsGrid, dsVert, dsHor

	def get_values_on_grid(self,
			xs=None,
			ys=None,
			zs=None,
			nx=None,
			ny=None,
			nz=None,
			) :
		"""
		"""
		try:
			if xs is None :
				if nx is None :
					nx=len(self._g_xvals)
				xs=np.linspace(self._g_xvals.min(),self._g_xvals.max(),nx)

			if ys is None :
				if ny is None :
					ny=len(self._g_yvals)
				ys=np.linspace(self._g_yvals.min(),self._g_yvals.max(),ny)

			if zs is None :
				if nz is None :
					nz=len(self._g_zvals)
				zs=np.linspace(self._g_zvals.min(),self._g_zvals.max(),nz)
			npnts=len(xs)*len(ys)*len(zs)
			pnts=np.zeros( (npnts,3) )
			indPnts=0
			for indx, xval in enumerate(xs) :
				for indy, yval in enumerate(ys) :
					for indz, zval in enumerate(zs) :
						pnts[indPnts,:]=np.array([xval,yval,zval])
						indPnts=indPnts+1

			interpValues=self.get_values(pnts=pnts)
		except unduGridError as e :
			print(e)
		return interpValues, xs, ys, zs

	def calcGradient(self) : 

		gradients=[]
		for i in range(self._g_nvals) :
			funs=self._g_funs[:,:,:,i]
			if len(self._g_xvals) < 2 :
				myGradient=np.gradient(funs[0,:,:],self._g_yvals,self._g_zvals,edge_order=2)
				grad0=np.zeros( (1,funs.shape[1],funs.shape[2]) ) 
				grad0[0,:,:]=myGradient[0]
				grad1=np.zeros( (1,funs.shape[1],funs.shape[2]) ) 
				grad1[0,:,:]=myGradient[1]

				gradients.append(grid_interpolator(		
					funs=grad0.reshape((1,len(self._g_yvals),len(self._g_zvals),1)),
					xvals=copy.deepcopy(self._g_xvals),
					yvals=copy.deepcopy(self._g_yvals),
					zvals=copy.deepcopy(self._g_zvals),
					))
				gradients.append(grid_interpolator(		
					funs=grad1.reshape((1,len(self._g_yvals),len(self._g_zvals),1)),
					xvals=copy.deepcopy(self._g_xvals),
					yvals=copy.deepcopy(self._g_yvals),
					zvals=copy.deepcopy(self._g_zvals),
					))

			else :
				myGradient=np.gradient(funs,self._g_xvals,self._g_yvals,self._g_zvals,edge_order=2)

				gradients.append(grid_interpolator(		
					funs=myGradient[0].reshape((len(self._g_xvals),len(self._g_yvals),len(self._g_zvals),1)),
					xvals=copy.deepcopy(self._g_xvals),
					yvals=copy.deepcopy(self._g_yvals),
					zvals=copy.deepcopy(self._g_zvals),
					))
				gradients.append(grid_interpolator(		
					funs=myGradient[1].reshape((len(self._g_xvals),len(self._g_yvals),len(self._g_zvals),1)),
					xvals=copy.deepcopy(self._g_xvals),
					yvals=copy.deepcopy(self._g_yvals),
					zvals=copy.deepcopy(self._g_zvals),
					))
				gradients.append(grid_interpolator(		
					funs=myGradient[2].reshape((len(self._g_xvals),len(self._g_yvals),len(self._g_zvals),1)),
					xvals=copy.deepcopy(self._g_xvals),
					yvals=copy.deepcopy(self._g_yvals),
					zvals=copy.deepcopy(self._g_zvals),
					))
		return gradients

	def get_line_through_point(self,p0,n0,coordLimits=None) :
		try:
			self.set_limits()
			if not self.is_in_limits(pnt=p0) : 
				raise unduGridError("grid_interpolator: get_values_along_line: p0 not in map volume.")
			norm=np.linalg.norm(n0)
			n0=n0/norm
			xmax=self._g_xmax
			xmin=self._g_xmin
			ymax=self._g_ymax
			ymin=self._g_ymin
			zmax=self._g_zmax
			zmin=self._g_zmin
			if not (coordLimits is None) :
				if len(coordLimits) > 2 :
					if not (coordLimits[0] is None) :
						xmin=coordLimits[0][0]
						xmax=coordLimits[0][1]
					if not (coordLimits[1] is None) :
						ymin=coordLimits[1][0]
						ymax=coordLimits[1][1]
					if not (coordLimits[2] is None) :
						zmin=coordLimits[2][0]
						zmax=coordLimits[2][1]
			maxs=np.array([xmax,ymax,zmax])
			mins=np.array([xmin,ymin,zmin])
			with np.errstate(divide='ignore', invalid='ignore'):
				tsMaxT=(maxs-p0)/n0
				tsMinT=(mins-p0)/n0
			tsMax=[]
			tsMin=[]
			for t in tsMaxT : 
				if (not np.isinf(t)) and ( not np.isnan(t) ) :
					tsMax.append(t)
			for t in tsMinT : 
				if (not np.isinf(t)) and ( not np.isnan(t) ) :
					tsMin.append(t)
			minT=tsMin[np.argmin(np.abs(tsMin))]
			maxT=tsMax[np.argmin(np.abs(tsMax))]
			pMin=p0+n0*minT
			pMax=p0+n0*maxT
		except unduGridError as e :
			print(e)
		return pMin, pMax

	def get_values_along_line(self,
			p0,
			n0,
			dx=None,
			dy=None,
			dz=None,
			n=100,
			coordLimits=None,
			) :

		pMin, pMax = self.get_line_through_point(
			p0=p0,
			n0=n0,
			coordLimits=coordLimits,
			)

		n0N=n0/np.linalg.norm(n0)
		distP=np.linalg.norm(p0-pMin)
		totDist=np.linalg.norm(pMax-pMin)
		tVal0=distP
		sVal0=np.dot(p0,n0N)
		ds0=sVal0-tVal0
		nTs=np.linspace(0,1,n)
		pnts=np.zeros( (n,3) )
		for indt, t in enumerate(nTs) :
			pnts[indt,:]=pMin+t*(pMax-pMin)

		vals=self(pnts)

		ds=np.zeros((pnts.shape[0],))
		for ind,pnt in enumerate(pnts[1:]) :
			diff=pnt-pnts[ind]
			diff=np.linalg.norm(diff)
			ds[ind+1]=ds[ind]+diff
		ds=ds+ds0
		return pnts, vals, ds

	def integrate_along_line(self,
			p0,
			n0,
			limit=400,
			epsrel=1e-3,
			epsabs=1e-5,
			coordLimits=None,
			square=False,
		) : 

		pnts,vals,ds = self.get_values_along_line(
			p0=p0,
			n0=n0,
			dx=None,
			dy=None,
			dz=None,
			n=limit,
			coordLimits=coordLimits,
			)

		firstInts=[]
		for i in range(self._g_nvals) :
			valsIntgrt=vals[:,i]
			if square :
				valsIntgrt=np.square(valsIntgrt)
			cs_integr = CubicSpline(ds, valsIntgrt)
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

		return pnts, ds, firstInts

class bmap_hallbach_parametrization :

	def __init__(self,
			**kwargs
			) : 

		super(bmap_hallbach_parametrization, self).__init__(
			**kwargs
			)

	def create_hallbach_parametrization(self) :
		"""
		create interpolation from bvals, xvals,...
		"""

		return bfield

	def get_field(self,pos) :
		"""
		also set 0 outside of area?
		y>?
		x
		"""
		field=getFieldFromRepresentation(pos)

		return field

