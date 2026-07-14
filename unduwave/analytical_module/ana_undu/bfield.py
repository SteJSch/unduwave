"""
Contains the functionality for loading and processing b-field data
"""
from unduwave.unduwave_incl import *
import unduwave.helpers.file_folder_helpers as ff_h
import unduwave.helpers.numerical_helpers as n_h
import unduwave.quantities.quantities as quantities
from unduwave.attribute_classes.attributes import _attribute
from unduwave.attribute_classes.attributes import _attribute_collection
import unduwave.analytical_module.ana_undu.grids as grid
import unduwave.analytical_module.ana_undu.kicks as kicks
import unduwave.constants as uc

class unduBfieldError(Exception) :
	pass

def beffs_init_shared(all_pnts,xs,prd_lngth,n_max) :
	global SHARED
	SHARED={'pnts':all_pnts, 'prd_lngth' : prd_lngth, 'n_max' : n_max,'xs':xs}

def beffs_worker(cur_inds):
	indy=cur_inds[0]
	indz=cur_inds[1]
	pnts=SHARED['pnts']
	data_x=SHARED['xs']
	prd_lngth=SHARED['prd_lngth']
	n_max=SHARED['n_max']

	beffs=np.zeros( (pnts.shape[3]) )
	for indDim in range(pnts.shape[3]) :

		dataB=pnts[:,indy,indz,indDim]

		in_spline = CubicSpline(data_x , dataB)

		b0, bc, bs = n_h.fourier_series_coeff_numpy(
			f=in_spline, 
			T=prd_lngth, 
			d0=data_x[0],
			N=n_max, 
			return_complex=False,
			)

		beff = (b0/2.0)**2
		for ind,tmp in enumerate(bs) :
			beff = beff + (tmp**2+bc[ind]**2)/(ind+1)**2

		beffs[indDim]=math.sqrt(beff)

	return beffs

def beff_grid_calc(all_pnts, xs, prd_lngth, n_max=10, processes=6):
	# Must guard pool creation when module may be imported on Windows
	try:
		if len(all_pnts.shape)<4:
			raise unduBfieldError("beff_grid_calc: Wrong shape of all_pnts")
	except unduBfieldError as e: 
		print(e)
		return
	shape=all_pnts.shape
	allindices=[]
	for indy in range(shape[1]) :
		for indz in range(shape[2]) :
			allindices.append([indy,indz])
	with multiprocessing.Pool(processes=processes, initializer=beffs_init_shared, initargs=(all_pnts,xs,prd_lngth,n_max,)) as p:
		beffsall=p.map(beffs_worker, allindices)

	npbeffsall=np.array(beffsall)

	beffsalln = np.zeros( (1,)+npbeffsall.shape )
	beffsalln[0,:,:]=npbeffsall

	return beffsalln

def create_indefinite_integral_cs( spline, xs, integr_limit) : 
	xmin = min(xs)
	y_1st = []
	for ind, x in enumerate(xs) : 
		if ind == 0 :
			y_1st.append( quad( spline, xmin, x, limit = integr_limit )[0] )
		else :
			y_1st.append( y_1st[-1] + quad( spline, xs[ind-1], x, limit = integr_limit )[0] )
	cs_integr = CubicSpline(xs, y_1st)
	return cs_integr

class bfield(quantities.quantity) : 
	"""
	Holds the bfield data
	"""

	def __init__(self, 
			unitsXB=[0.001,1.0],
			**kwargs,
			) :

		super(bfield, self).__init__(**kwargs)
		self._unitsXB=unitsXB
		self._data=np.ones((0,0))

		self.bvals=np.array([])
		self.xvals=np.array([0.0])
		self.zvals=np.array([0.0])
		self.yvals=np.array([0.0])
		self._interpolator=None

	"""
	bfield class

	has x,y,z,bx,by,bz

	if has x
	find new x before, with right interval
	find new x after

	loop new x+oldx+newx make harm bm, add to bm already there

	loop new xs
	add 0 for all other bs
	"""

	def write_field(self,
			file,
			cols=None,
			outType='std',
			unitsXB=None,
			filez=None,
			filex=None,
			data=None,
			) :
		path, fileN = os.path.split(file)
		os.makedirs(path, exist_ok=True)
		if outType == 'std' : 
			# if field == 'map' :
			# 	pass
			# else :
			self.write_field_std(file=file,unitsXB=unitsXB,whatStr=cols)
		elif outType == 'waveBy' : 
			self.write_field_std(file=file,unitsXB=unitsXB)
		elif outType == 'waveByz' : 
			self.write_field_waveByz(filey=file,filez=filez,unitsXB=None)
		elif outType == 'waveBxyz' : 
			self.write_field_waveBxyz(filey=file,filex=filex,filez=filez,unitsXB=None)
		elif outType == 'unduOut' : 
			self.write_field_unduOut(file=file,unitsXB=unitsXB)
		elif outType == 'bmap' : 
			self.write_field_map(file=file,unitsXB=unitsXB,wave=False)
		elif outType == 'bmapWAVE' : 
			self.write_field_map(file=file,unitsXB=unitsXB,wave=True)
		elif outType == 'bmapUNDU' : 
			self.write_field_map_undu(file=file,unduMap=data,unitsXB=unitsXB)

		# convert_x_mm_b_T_file_to_wave_std(file_in, out_path )

	def write_field_unduOut(self,file,unitsXB=None,colx='x') :
		"""
		UndumagOut file Units are mm for the length and T for B
		"""
		path, fileN = os.path.split(file)
		os.makedirs(path, exist_ok=True)
		if unitsXB is None :
			unitsXB=self._unitsXB
		unitConvX=unitsXB[0]/0.001
		unitConvB=unitsXB[1]/1
		x_vals=self.get_xvals(colx=colx)
		bvals = self.bvals
		with open( file, 'w') as o_f:
			for ind, xval in enumerate(x_vals) :
				o_f.write(f"  {xval*unitConvX}   {bvals[ind,1]*unitConvB}   {bvals[ind,2]*unitConvB}   0.0000000E+000   0.0000000E+000   0.0000000E+000   0.0000000E+000          0\n")		

	def write_field_std(self,file,unitsXB=None, whatStr='', colx='x') :
		"""
		UndumagOut file Units are mm for the length and T for B
		"""
		path, fileN = os.path.split(file)
		os.makedirs(path, exist_ok=True)
		if unitsXB is None :
			unitsXB=self._unitsXB
		unitConvX=self._unitsXB[0]/unitsXB[0] # here the x-coord has to be in m?!!!
		unitConvB=self._unitsXB[1]/unitsXB[1]
		# unitConvX=unitsConv[0]/0.001 # x in mm
		# unitConvB=unitsConv[1]
		x_vals=self.get_xvals(colx=colx)
		bvals=self.bvals
		doBx=True
		doBy=True
		doBz=True
		if len(whatStr) > 0 :
			if not ('Bx' in whatStr) :
				doBx=False
			if not ('By' in whatStr) :
				doBx=False
			if not ('Bz' in whatStr) :
				doBz=False
		with open( file, 'w') as o_f:
			for ind, xval in enumerate(x_vals) :
				lineStr=f"{xval*unitConvX}"
				if doBx :
					lineStr=lineStr + f" {bvals[ind,0]*unitConvB}"
				if doBy :
					lineStr=lineStr + f" {bvals[ind,1]*unitConvB}"
				if doBz :
					lineStr=lineStr + f" {bvals[ind,2]*unitConvB}"
				lineStr=lineStr+'\n'
				o_f.write(lineStr)

	def write_field_waveByz(self,filey,filez,unitsXB=None,colx='x') :
		"""
		UndumagOut file Units are mm for the length and T for B
		"""
		path, fileN = os.path.split(file)
		os.makedirs(path, exist_ok=True)
		if unitsXB is None :
			unitsXB=self._unitsXB
		unitConvX=unitsXB[0]/1 # here the x-coord has to be in m?!!!
		unitConvB=unitsXB[1]/1
		x_vals=self.get_xvals(colx=colx)
		b_vals = self.bvals
		with open( filey, 'w') as o_f:
			o_f.write('Comment\n')
			o_f.write('1.0 1.0\n')
			o_f.write(str(len(x_vals))+'\n')
			for ind, xval in enumerate(x_vals) :
				o_f.write(f"{xval*unitConvX} {b_vals[ind,1]*unitConvB}\n")
		with open( filez, 'w') as o_f:
			o_f.write('Comment\n')
			o_f.write('1.0 1.0\n')
			o_f.write(str(len(x_vals))+'\n')
			for ind, xval in enumerate(x_vals) :
				o_f.write(f"{xval*unitConvX} {b_vals[ind,2]*unitConvB}\n")

	def write_field_waveBxyz(self,filex,filey,filez,unitsXB=None,colx='x') :
		"""
		Outputs 3 separate files with 1D-data
		"""
		path, fileN = os.path.split(filex)
		os.makedirs(path, exist_ok=True)
		path, fileN = os.path.split(filey)
		os.makedirs(path, exist_ok=True)
		if unitsXB is None :
			unitsXB=self._unitsXB
		unitConvX=unitsXB[0]/1 # here the x-coord has to be in m?!!!
		unitConvB=unitsXB[1]/1
		x_vals=self.xvals
		if len(x_vals.shape) == 1 :
			xvalList=x_vals
		elif len(x_vals.shape) == 2 :
			xvalList=x_vals[:,0]
		elif len(x_vals.shape) == 3 :
			xvalList=x_vals[:,0,0]
		bvals = self.bvals
		with open( filex, 'w') as o_f:
			o_f.write('Comment\n')
			o_f.write('1.0 1.0\n')
			o_f.write(str(len(xvalList))+'\n')
			for ind, xval in enumerate(xvalList) :
				o_f.write(f"{xval*unitConvX} {bvals[ind,0]*unitConvB}\n")
		with open( filey, 'w') as o_f:
			o_f.write('Comment\n')
			o_f.write('1.0 1.0\n')
			o_f.write(str(len(xvalList))+'\n')
			for ind, xval in enumerate(xvalList) :
				o_f.write(f"{xval*unitConvX} {bvals[ind,1]*unitConvB}\n")
		with open( filez, 'w') as o_f:
			o_f.write('Comment\n')
			o_f.write('1.0 1.0\n')
			o_f.write(str(len(xvalList))+'\n')
			for ind, xval in enumerate(xvalList) :
				o_f.write(f"{xval*unitConvX} {bvals[ind,2]*unitConvB}\n")

	def write_field_map(self,file,unitsXB=None,wave=False) :
		path, fileN = os.path.split(file)
		os.makedirs(path, exist_ok=True)
		if unitsXB is None :
			unitsXB=self._unitsXB
		unitConvX=1
		unitConvB=1
		if wave :
			unitConvX=unitsXB[0]/1 # xyz in wave field map are in m
			unitConvB=unitsXB[1]/1
		else :
			unitConvX=self._unitsXB[0]/unitsXB[0] # xyz in wave field map are in m
			unitConvB=self._unitsXB[1]/unitsXB[1]

		x_vals=self.xvals
		y_vals=self.yvals
		z_vals=self.zvals
		bvals = self.bvals
		# shape=x_vals.shape

		# nx=len(x_vals)
		# ny=len(y_vals)
		# nz=len(z_vals)
		# nlines=nx*ny*nz
		# xn=x_vals.reshape(nlines)
		# yn=y_vals.reshape(nlines)
		# zn=z_vals.reshape(nlines)
		# bComp=bvals.reshape(nlines,3)
		with open( file, 'w') as o_f:
			if wave :
				o_f.write('! WAVE: x y z Bx By Bz with x as long. beam axis\n')
				o_f.write('@ date (yyyy.month.day) and time =  2025.06.25 18:59:58\n')
				o_f.write('@ run =           1\n')
				o_f.write('@ comment = WAVE.EXAMPLE\n')
				o_f.write('@ scaling = 1.0 1.0 1.0 1.0 1.0 1.0\n')
				o_f.write('@ offset = 0.0, 0.0, 0.0 0.0 0.0 0.0\n')
			else :
				o_f.write('x y z Bx By Bz\n')
			for indx, xval in enumerate(np.unique(x_vals)) :
				for indy, yval in enumerate(np.unique(y_vals)) :
					for indz, zval in enumerate(np.unique(z_vals)) :
						o_f.write(f"{xval*unitConvX}  {yval*unitConvX}  {zval*unitConvX}  {bvals[indx,indy,indz,0]*unitConvB}  {bvals[indx,indy,indz,1]*unitConvB}  {bvals[indx,indy,indz,2]*unitConvB}\n")		

	def write_field_map_undu(self,file,unduMap,unitsXB=None) :
		path, fileN = os.path.split(file)
		os.makedirs(path, exist_ok=True)
		if unitsXB is None :
			unitsXB=self._unitsXB
		unitConvX=unitsXB[0]/1 # xyz in wave field map are in m
		unitConvB=unitsXB[1]/1

		with open( file, 'w') as o_f:

			o_f.write('* Some Undulator, Run: RunNr\n')
			o_f.write('* Some Date\n')
			o_f.write('* imoth imag mat ityp matmod x/mm y/mm z/mm Bx/T By/T Bz/T B/T Hx/T Hy/T Hz/T H/T Mx/T My/T Mz/T M/T BxDip/T ByDip/T BzDip/T ifail kfail cmag cmoth\n')

			for line in unduMap.to_dict('records') :
				for key, el in line.items() :
					o_f.write(f"{el}	")
				o_f.write("\n")

	def fromQuants(self,
			xQuant,
			byQuant,
			bzQuant=None,
			) :
		self.xvals=np.array( xQuant._data )
		nvals=len(self.xvals)
		self.bvals=np.zeros((nvals,3))
		self.bvals[:,1]=np.array( byQuant._data )
		if not (bzQuant is None) :
			self.bvals[:,2]=np.array( bzQuant._data )

	def clear_undumag_bmap_exact(self,undu_bmap,vHor,vVert,fileOut=None) :

		if len(vVert) < 2 :
			vVert=[0.0]
		if len(vHor) < 2 :
			vHor=[0.0]
		newData=[]
		for ind, mapLine in enumerate(undu_bmap.to_dict('records')) :
			xVal=mapLine['x']
			vv=mapLine['y']
			vh=mapLine['z']
			if True in np.isclose(vHor, vh, atol=1e-3, rtol=1e-3) : 
				if True in np.isclose(vVert, vv, atol=1e-3, rtol=1e-3) : 
					newData.append(mapLine)

		newData=pd.DataFrame(newData)
		unique_x=newData['x'].unique()
		unique_y=newData['y'].unique()
		unique_z=newData['z'].unique()

		numYZ=[]
		for indx, xval in enumerate(unique_x) :

			if indx < (len(unique_x)-1) :
				nextX=unique_x[indx+1]
			else :
				nextX=None
			if indx > 0 :
				lastX=unique_x[indx-1]
			else :
				lastX=None
			if abs(xval) > 0 :
				if (nextX is None) and (lastX is None) :
					fac=0.01
				elif nextX is None :
					distNext = abs(xval-lastX)
					fac=abs(distNext/(2*xval))
				elif lastX is None :
					distNext = abs(xval-nextX)
					fac=abs(distNext/(2*xval))
				else :
					distNext=abs(xval-nextX)
					distLast=abs(xval-lastX)
					minDist=min(distNext,distLast)
					fac=abs(minDist/(2*xval))
			else :
				fac=0.0
			if xval >= 0 :
				small=xval*(1-fac)
				big=xval*(1+fac)
			else :
				big=xval*(1-fac)
				small=xval*(1+fac)

			cut_map=newData[ (newData['x'] >= small) & (newData['x'] <= big) ]
			y_vals=cut_map['y'].unique()
			z_vals=cut_map['z'].unique()
			numYZ.append({ 'xval' : xval, 'ny' : len(y_vals), 'nz' : len(z_vals), 'map' : cut_map })

		numYZ=pd.DataFrame(numYZ)
		zValueCounts=numYZ['nz'].value_counts().to_dict()
		yValueCounts=numYZ['ny'].value_counts().to_dict()

		maxZVal=None
		numMaxZ=0
		for key, el in zValueCounts.items() :
			if el > numMaxZ :
				numMaxZ=el
				maxZVal=key

		maxYVal=None
		numMaxY=0
		for key, el in yValueCounts.items() :
			if el > numMaxY :
				numMaxY=el
				maxYVal=key

		newMap=[]
		for el in numYZ.to_dict("records") :
			if el['ny'] == maxYVal :
				if el['nz'] == maxZVal :
					cut=el['map'].to_dict("records")
					newMap=newMap+cut
		newMap=pd.DataFrame(newMap)
		if not fileOut is None :
			self.write_field_map_undu(file=fileOut,unduMap=newMap,unitsXB=None)
		return newMap

	def load_clear_undu_map(self,
			mapFile,
			vHor,
			vVert, 
			mapOutName,
			):
		"""
		if mapOut None the mapFile will be replaced
		"""
		cols=[ 'imoth', 'imag', 'mat', 'ityp', 'matmod', 'x', 'y', 'z', 'Bx', 'By', 'Bz', 'B', 'Hx', 'Hy', 'Hz', 'H', 'Mx', 'My', 'Mz', 'M', 'BxDip', 'ByDip', 'BzDip', 'ifail', 'kfail', 'cmag', 'cmoth' ]
		skiprows=range(0, 3)
		data = pd.read_csv( mapFile, skiprows=skiprows, dtype=object, delimiter=r"\s+",header=None)
		data.columns = cols
		cols_float = ['x', 'y', 'z', 'Bx', 'By', 'Bz']
		for col in cols_float:
			data[col] = data[col].astype(float)

		data=self.clear_undumag_bmap_exact(
			undu_bmap=data,
			vHor=vHor,
			vVert=vVert,
			fileOut=mapOutName
			)

		unique_x=data['x'].unique()
		unique_y=data['y'].unique()
		unique_z=data['z'].unique()

		nxmap=len(unique_x)
		nymap=len(unique_y)
		nzmap=len(unique_z)

		self.xvals=unique_x
		self.yvals=unique_y
		self.zvals=unique_z

		self.bvals=np.zeros((nxmap,nymap,nzmap,3))
		try:
			bxGrid=np.array(data['Bx'].to_list()).reshape((nxmap,nymap,nzmap))
			byGrid=np.array(data['By'].to_list()).reshape((nxmap,nymap,nzmap))
			bzGrid=np.array(data['Bz'].to_list()).reshape((nxmap,nymap,nzmap))
		except:
			print("bfield: load_field_from_file: Grid inconsistent.")
			return
		self.bvals[:,:,:,0]=bxGrid
		self.bvals[:,:,:,1]=byGrid
		self.bvals[:,:,:,2]=bzGrid
		return

	def load_field_from_file(self,
			file, 
			fieldMap=False,
			cols=None,
			unduFile = False, 
			radiaFile=False,
			header=None,
			skiprows=None,
			) :
		"""
		If radiaFile, unduFile and fieldMap are false, the file loaded is supposed
		to contain 1 space-dimension (ordered) and up to 3 B-components

		Undufile=True : The on-axis field file from an undumag simulation is loaded,
		that is 1-space dimension data with by,bz components.

		fieldMap=True, 
		unduFile=False
		unduFile=True
		Undumag field-map file is loaded
		general 3D map (x,y,z) with full B-components
		fileOut : If a field-map is loaded from undumag and cleared, it can be put out into this

		"""

		if unduFile and (not fieldMap) :
			data = pd.read_csv( file, dtype=object, sep='\\s+', header = header)
			data.columns = ['x','By','Bz','intBy','intBz','int2By','int2Bz','quark']
			for col in data.columns:
				data[col] = data[col].astype(float)

			self.xvals=np.round(np.array(data['x'].to_list()),8)
			nvals=len(self.xvals)
			self.bvals=np.zeros((nvals,3))
			self.bvals[:,1]=data['By'].to_list()
			self.bvals[:,2]=data['Bz'].to_list()
			return
		# field map
		if fieldMap :
			if unduFile : 
				cols=[ 'imoth', 'imag', 'mat', 'ityp', 'matmod', 'x', 'y', 'z', 'Bx', 'By', 'Bz', 'B', 'Hx', 'Hy', 'Hz', 'H', 'Mx', 'My', 'Mz', 'M', 'BxDip', 'ByDip', 'BzDip', 'ifail', 'kfail', 'cmag', 'cmoth' ]
				skiprows=range(0, 3)
			elif radiaFile :
				cols=[ "z", "x", "y", "Bz", "Bx", "By", "B", "Hz", "Hx", "Hy", "H", "Mz", "Mx", "My", "M" ]
				skiprows=0
			elif cols is None :
				print("bfield: load_field_from_file: cols has to be set for general map file")
				return
			data = pd.read_csv( file, skiprows=skiprows, dtype=object, delimiter=r"\s+",header=header)
			data.columns = cols
			cols_float = ['x', 'y', 'z', 'Bx', 'By', 'Bz']
			for col in cols_float:
				data[col] = data[col].astype(float)

			if radiaFile :
				data['x']=np.round(data['x'],8)
				unique_x=data['x'].unique()
				unique_y=data['y'].unique()
				unique_z=data['z'].unique()
				nxmap=len(unique_x)
				nymap=len(unique_y)
				nzmap=len(unique_z)
				dataNP=np.zeros( ( len(unique_y), len(unique_z), len(unique_x), 3 ) )			
				dataNP[:,:,:,0]=np.array(data['Bx']).reshape( ( len(unique_y), len(unique_z), len(unique_x) ) )
				dataNP[:,:,:,1]=np.array(data['By']).reshape( ( len(unique_y), len(unique_z), len(unique_x) ) )
				dataNP[:,:,:,2]=np.array(data['Bz']).reshape( ( len(unique_y), len(unique_z), len(unique_x) ) )
				dataNP=dataNP.swapaxes(0,2).swapaxes(1,2)

				try:
					bxGrid=dataNP[:,:,:,0]
					byGrid=dataNP[:,:,:,1]
					bzGrid=dataNP[:,:,:,2]
				except:
					print("bfield: load_field_from_file: Grid inconsistent.")
					return

			else :
				unique_x=data['x'].unique()
				unique_y=data['y'].unique()
				unique_z=data['z'].unique()
				nxmap=len(unique_x)
				nymap=len(unique_y)
				nzmap=len(unique_z)
				try:
					bxGrid=np.array(data['Bx'].to_list()).reshape((nxmap,nymap,nzmap,))
					byGrid=np.array(data['By'].to_list()).reshape((nxmap,nymap,nzmap,))
					bzGrid=np.array(data['Bz'].to_list()).reshape((nxmap,nymap,nzmap,))
				except:
					print("bfield: load_field_from_file: Grid inconsistent.")
					return

			self.xvals=np.round(unique_x,8)
			self.yvals=np.round(unique_y,8)
			self.zvals=np.round(unique_z,8)

			self.bvals=np.zeros((nxmap,nymap,nzmap,3))

			self.bvals[:,:,:,0]=bxGrid
			self.bvals[:,:,:,1]=byGrid
			self.bvals[:,:,:,2]=bzGrid

			return
			
		if cols is None :
			print("bfield: load_field_from_file: cols has to be set for general file")
			return
		data = pd.read_csv( file, dtype=object, sep='\\s+', header = header, skiprows=skiprows)
		data.columns = cols
		for col in data.columns:
			data[col] = data[col].astype(float)
		xvals=[]
		yvals=[]
		zvals=[]
		bxvals=[]
		byvals=[]
		bzvals=[]
		nvalsx=1
		nvalsy=1
		nvalsz=1
		if 'x' in cols :
			self.xvals=np.round(np.array(data['x'].to_list()),8)
			nvals=len(self.xvals)
		elif 'y' in cols :
			self.yvals=np.round(np.array(data['y'].to_list()),8)
			nvals=len(self.yvals)
		elif 'z' in cols :
			self.zvals=np.round(np.array(data['z'].to_list()),8)
			nvals=len(self.zvals)
		self.bvals=np.zeros((nvals,3))
		if 'Bx' in cols :
			self.bvals[:,0]=data['Bx'].to_list()
		if 'By' in cols :
			self.bvals[:,1]=data['By'].to_list()
		if 'Bz' in cols :
			self.bvals[:,2]=data['Bz'].to_list()

	def move(self,dist,colx) : 
		bfieldO=self.clone_me()
		xvals=bfieldO.get_xvals(colx=colx)
		xvals=xvals+dist
		bfieldO.xvals=xvals
		return bfieldO

	def center(self,colx = 'x'):
		bfieldO=self.clone_me()
		xvals=bfieldO.get_xvals(colx=colx)
		x0=xvals[0]
		xE=xvals[-1]
		delta=-0.5*(x0+xE)
		bfieldO=bfieldO.move(dist=delta,colx=colx)
		return bfieldO

	def gauge(self, gauge_fac ) : 
		bf = self.clone_me()
		bf.bvals=bf.bvals*gauge_fac
		return bf

	def expandPeriodic(self, periodLength, nperiodsAdd, start=0.0 ) : 
		bf = self.clone_me()
		xvals=bf.xvals
		bvals=bf.bvals

		# bspline=CubicSpline(xvals,bvals[:,1])

		# xvalsIntrp=np.linspace(start-periodLength/2.0,start+periodLength*3.0/2.0,1000)
		# bvalsIntrp=bspline(xvalsIntrp)

		# xvalPer=[]
		# bvalPer=[]
		# for indx,xval in enumerate(xvalsIntrp) :
		# 	if xval >= start : 
		# 		if xval < (start+periodLength) :
		# 			xvalPer.append(xval)
		# 			bvalPer.append([0.0,bvalsIntrp[indx],0.0])
		# bf.xvals=np.array(xvalPer)
		# bf.bvals=np.array(bvalPer)
		# n_bint1, n_bint2, bint1z, bint2z = bf.getTotalFieldIntegrals()
		# print(f"{n_bint1}, {n_bint2}, {bint1z}, {bint2z}")

		# return

		xvalsAdd=[]
		bvalsAdd=[]
		xvalsN=[]
		bvalsN=[]
		xvalsEnd=[]
		bvalsEnd=[]
		indStart=None
		for indx, xval in enumerate(xvals) :
			if xval < (start+periodLength) :
				xvalsN.append(xval)
				bvalsN.append(bvals[indx])
			if xval>=start :
				if indStart is None :
					indStart=indx-1
				if xval < (start+periodLength) :
					xvalsAdd.append(xval)
					bvalsAdd.append(bvals[indx,:])
				else :
					xvalsEnd.append(xval)
					bvalsEnd.append(bvals[indx,:])

		addX=periodLength
		for period in range(nperiodsAdd) :
			for indAdd, xval in enumerate(xvalsAdd) :
				xvalsN.append(xval+addX)
				bvalsN.append(bvalsAdd[indAdd])
			addX=addX+periodLength
		addX=periodLength*nperiodsAdd
		xvalsEnd=np.array(xvalsEnd)
		xvalsEnd=xvalsEnd+addX
		xvalsEnd=xvalsEnd.tolist()
		xvalsN=xvalsN+xvalsEnd
		xvalsN=np.array(xvalsN)

		bvalsN=bvalsN+bvalsEnd
		bvalsN=np.array(bvalsN)

		bf.xvals=xvalsN
		bf.bvals=bvalsN

		return bf

	def get_pos_ind(self,colx,valx) :
		vals=self.get_xvals(colx=colx)
		for ind, val in enumerate(vals) :
			if math.isclose(val, valx, rel_tol=1e-5, abs_tol=0.01) :
				return ind

		differences = np.abs(vals - valx)
		closest_index = differences.argmin()
		# closest_value = vals[closest_index]

		return closest_index

	def get_xvals(self,colx) :
		if colx == 'x' :
			xvals=self.xvals
		elif colx == 'y' :
			xvals=self.yvals
		elif colx == 'z' :
			xvals=self.zvals
		return xvals

	def set_xvals(self,xvals,colx) :
		if colx == 'x' :
			self.xvals=xvals
		elif colx == 'y' :
			self.yvals=xvals
		elif colx == 'z' :
			self.zvals=xvals

	def add_field(self, 
			bfieldA, 
			pos = None,#'front', 'back'
			colx = 'x',
			pntsOut=None,
			) : 
		"""
		Glues the field bfieldA to this field.
		pos=None - The fields are overlayed as they are
		'front' - start of bfieldA is shifted onto end of this data
		'back' - end of bfieldA is shifted to start of this data

		"""
		bfieldO=self.clone_me()
		xvalsO=bfieldO.get_xvals(colx=colx)
		xvalsA=bfieldA.get_xvals(colx=colx)

		deltaX=0.0
		if pos=='back' :
			deltaX=xvalsO[-1]-xvalsA[0]
		if pos=='front' :
			deltaX=(xvalsO[0]-xvalsA[-1])
		xvalsA=xvalsA+deltaX

		f_spline_b_orig=CubicSpline(xvalsO , bfieldO.bvals)
		f_spline_b_add=CubicSpline(xvalsA , bfieldA.bvals)

		xO0=xvalsO[0]
		xOE=xvalsO[-1]
		xA0=xvalsA[0]
		xAE=xvalsA[-1]

		xStart=min(xO0,xA0)
		xEnd=max(xOE,xAE)

		if pntsOut is None :
			dnstO=len(xvalsO)/(xOE-xO0)
			dnstA=len(xvalsA)/(xAE-xA0)
			dnst=max(dnstO,dnstA)
			npnts=int(dnst*(xEnd-xStart))
			newX=np.linspace(xStart,xEnd,npnts)
		else :
			newX=pntsOut

		nBs=[]
		for xV in newX :
			nB=np.array([0.0,0.0,0.0])
			if (xV>=xO0) and (xV<=xOE) :
				nB=nB+f_spline_b_orig(xV)
			if (xV>=xA0) and (xV<=xAE) :
				nB=nB+f_spline_b_add(xV)
			nBs.append(nB)
		bvals=np.zeros((len(newX),3))
		bvals[:]=nBs

		bfieldO.set_xvals(
			xvals=newX,
			colx=colx
			)
		bfieldO.bvals=bvals
		return bfieldO


	def create_harm_field(
			self,
			periodLength, 
			amplitude, 
			nperiods=None, 
			phase_shift = 0, 
			num_pnts_per_period = 100, 
			colx = 'x', 
			coly = 'By',
			deltaX=None,
			) : 
		"""
		Creates a sine field: 
		amplitude*sin( k_per * x + phase_shift ) for all x lying in deltaX
		or if deltaX=None, it is set to 
		[ -nperiods*periodLength/2.0, nperiods*periodLength/2.0 ]
		creates bfield class filled with data, with num_pnts_per_period pnts per period in deltaX

		"""
		if (deltaX is None) :
			deltaX=[ -nperiods*periodLength/2.0, nperiods*periodLength/2.0 ]
		length = deltaX[-1] - deltaX[0]
		numPerAct=(length/periodLength)
		xpnts=np.linspace(deltaX[0],deltaX[1],int(numPerAct*num_pnts_per_period))

		bvalN=[]
		for pnt in xpnts : 
			val= amplitude*math.sin( 2*math.pi/periodLength * pnt + phase_shift ) 
			if coly == 'Bx' :
				bvalN.append([val,0.0,0.0])
			elif coly == 'By' :
				bvalN.append([0.0,val,0.0])
			elif coly == 'Bz' :
				bvalN.append([0.0,0.0,val])
		bvals=np.zeros((len(xpnts),3))
		bvals[:]=bvalN
		newFld=bfield(unitsXB=self._unitsXB)
		if colx == 'x' :
			newFld.xvals=xpnts
		elif colx == 'y' :
			newFld.yvals=xpnts
		elif colx == 'z' :
			newFld.zvals=xpnts
		newFld.bvals=bvals
		return newFld

	def create_harm_field_with_ends(
			self,
			periodLength, 
			amplitude, 
			nperiods,
			num_pnts_per_period = 100, 
			colx = 'x', 
			coly = 'By',
			shift=0,
			) :

		periodField=self.create_harm_field(
			periodLength=periodLength, 
			amplitude=amplitude, 
			phase_shift = math.pi, 
			num_pnts_per_period = num_pnts_per_period, 
			colx = colx, 
			coly = coly,
			deltaX=[0,periodLength*nperiods], 
			)

		be1=self.create_harm_field(
			periodLength=periodLength, 
			amplitude=amplitude*0.75, 
			phase_shift = 0, 
			num_pnts_per_period = num_pnts_per_period/2.0, 
			colx = colx, 
			coly = coly,
			deltaX=[0,periodLength/2.0], 
			)
		be2=self.create_harm_field(
			periodLength=periodLength, 
			amplitude=amplitude*0.25, 
			phase_shift = math.pi, 
			num_pnts_per_period = num_pnts_per_period/2.0, 
			colx = colx, 
			coly = coly,
			deltaX=[0,periodLength/2.0], 
			)
		be3=self.create_harm_field(
			periodLength=periodLength, 
			amplitude=amplitude*0.75, 
			phase_shift = math.pi, 
			num_pnts_per_period = num_pnts_per_period/2.0, 
			colx = colx, 
			coly = coly,
			deltaX=[0,periodLength/2.0], 
			)
		be4=self.create_harm_field(
			periodLength=periodLength, 
			amplitude=amplitude*0.25, 
			phase_shift = 0.0,
			num_pnts_per_period = num_pnts_per_period/2.0, 
			colx = colx, 
			coly = coly,
			deltaX=[0,periodLength/2.0], 
			)

		bfully=be1.add_field(
			bfieldA=be2, 
			pos = 'front',#'front', 'back', None
			colx=colx,
			pntsOut=None,
			)
		bfully=bfully.add_field(
			bfieldA=periodField, 
			pos = 'back',#'front', 'back', None
			colx=colx,
			pntsOut=None,
			)
		bfully=bfully.add_field(
			bfieldA=be3, 
			pos = 'back',#'front', 'back', None
			colx=colx,
			pntsOut=None,
			)
		bfully=bfully.add_field(
			bfieldA=be4, 
			pos = 'back',#'front', 'back', None
			colx=colx,
			pntsOut=None,
			)
		bfully=bfully.center(colx = 'x')
		return bfully

	def clone_me(self) : 
		bf = bfield()
		bf = copy.deepcopy(self)
		return bf

	def calc_integrals_grid(self,
			longIntrvl=None,
			nlongs=None,
			nhor=None,
			dhor=None,
			horPnts=None,
			horIntrvl=None,
			nvert=None,
			dvert=None,
			vertPnts=None,
			vertIntrvl=None,
			intrgnt_limit=400, 
			epsrel=1e-3,
			epsabs=1e-5,
			processes=1,
			) :

		myInterpol=self.create_grid_interpolation()
		myInterpol.processes=processes

		vertPnts=grid.createPointList(
			vals=self.yvals,npnts=nvert,dpnts=dvert,pnts=vertPnts,intrvl=vertIntrvl)
		horPnts=grid.createPointList(
			vals=self.zvals,npnts=nhor,dpnts=dhor,pnts=horPnts,intrvl=horIntrvl)
		longPnts=grid.createPointList(
			vals=self.xvals,npnts=nlongs,dpnts=None,pnts=None,intrvl=longIntrvl)
		if longPnts is None :
			longPnts=self.xvals

		xvals=longPnts

		# print("here")#???
		nfuns, xs, ys, zs=myInterpol.get_values_on_grid(
			xs=longPnts,
			ys=vertPnts,
			zs=horPnts,
			nx=None,
			ny=None,
			nz=None,
			)
		nfuns=nfuns.reshape((len(xs),len(ys),len(zs),nfuns.shape[-1]))

		# print("here2")
		line_intgrls=grid.grid_line_integrals_calc(
			funs=nfuns, 
			xvals=longPnts,
			yvals=vertPnts,
			zvals=horPnts,
			limit=intrgnt_limit,
			epsrel=epsrel,
			epsabs=epsabs,
			processes=processes,
			)

		firstIntFuns=np.zeros( (len(xvals),len(line_intgrls),3) )

		# print("here3")
		for ind,splineVec in enumerate(line_intgrls) : 
			longSpline=splineVec[0]
			vertSpline=splineVec[1]
			horSpline=splineVec[2]
			firstIntFuns[:,ind,0]=longSpline(xvals)
			firstIntFuns[:,ind,1]=vertSpline(xvals)
			firstIntFuns[:,ind,2]=horSpline(xvals)

		integralsBxByBz=firstIntFuns.reshape( (len(xvals),len(vertPnts),len(horPnts),3) )

		# print("here4")
		firstIntGrid=grid.grid_interpolator(
			funs=integralsBxByBz,
			xvals=xvals,
			yvals=vertPnts,
			zvals=horPnts,
			)

		# print("here5")
		scndLineInts=grid.grid_line_integrals_calc(
			funs=firstIntGrid._g_funs, 
			xvals=longPnts,
			yvals=vertPnts,
			zvals=horPnts,
			limit=intrgnt_limit,
			epsrel=epsrel,
			epsabs=epsabs,
			processes=processes,
			)

		scndIntFuns=np.zeros( (len(xvals),len(line_intgrls),3) )

		# print("here6")
		for ind,splineVec in enumerate(scndLineInts) : 
			longSpline=splineVec[0]
			vertSpline=splineVec[1]
			horSpline=splineVec[2]
			scndIntFuns[:,ind,0]=longSpline(xvals)
			scndIntFuns[:,ind,1]=vertSpline(xvals)
			scndIntFuns[:,ind,2]=horSpline(xvals)

		scndIntegralsBxByBz=scndIntFuns.reshape( (len(xvals),len(vertPnts),len(horPnts),3) )

		# print("here7")
		scndIntGrid=grid.grid_interpolator(
			funs=scndIntegralsBxByBz,
			xvals=xvals,
			yvals=vertPnts,
			zvals=horPnts,
			)
		firstIntGrid.processes=processes
		scndIntGrid.processes=processes

		return firstIntGrid, scndIntGrid

	def calc_kicks_grid(self,
			prd_lngth,
			longIntrvl=None,
			nlongs=None,
			nhor=None,
			dhor=None,
			horPnts=None,
			horIntrvl=None,
			nvert=None,
			dvert=None,
			vertPnts=None,
			vertIntrvl=None,
			intrgnt_limit=400, 
			n_max = 10, 
			epsrel=1e-3,
			epsabs=1e-5,
			processes=1,
			method='harmonic',# 'full', 'harmonic',
			beamEnGeV=1, # Beam Energy in GeV
			) :


		facKicks=-1e-6*1/2*(uc.v_c*1e-9/beamEnGeV)**2

		if method=='harmonic' :

			beffs, ys, zs, longIntrvl=self.calc_beff_grid(
				prd_lngth=prd_lngth,
				nlongs=nlongs,
				longIntrvl=longIntrvl,
				nhor=nhor,
				dhor=dhor,
				horPnts=horPnts,
				horIntrvl=horIntrvl,
				nvert=nvert,
				dvert=dvert,
				vertPnts=vertPnts,
				vertIntrvl=vertIntrvl,
				intrgnt_limit=intrgnt_limit, 
				n_max = n_max, 
				epsrel=epsrel,
				epsabs=epsabs,
				processes=processes,
				)

			nshape=beffs._g_funs.shape

			facFocusPot=1/(8*(math.pi**2))*prd_lngth**3
			focusPot = np.zeros( (nshape[0],nshape[1],nshape[2],1) )
			# focusPot=copy.deepcopy(beffs.funs)
			tmp=(beffs._g_funs**2)*facFocusPot
			focusPot[:,:,:,0]=0*tmp[:,:,:,0]+tmp[:,:,:,1]+tmp[:,:,:,2]
			focusPotGrid=grid.grid_interpolator(
					funs=focusPot,
					xvals=np.array([longIntrvl[-1]]),
					yvals=ys,
					zvals=zs,
				)
			gradients=focusPotGrid.calcGradient()
			kickY=gradients[0]
			kickZ=gradients[1]
			kickY.setfuns(funs=kickY._g_funs*facKicks)
			kickZ.setfuns(funs=kickZ._g_funs*facKicks)

			return kickY, kickZ, focusPotGrid

		elif method=='full' :

			facPot=1e-3
			myInterpol=self.create_grid_interpolation()

			vertPnts=grid.createPointList(
				vals=self.yvals,npnts=nvert,dpnts=dvert,pnts=vertPnts,intrvl=vertIntrvl)
			horPnts=grid.createPointList(
				vals=self.zvals,npnts=nhor,dpnts=dhor,pnts=horPnts,intrvl=horIntrvl)
			longPnts=grid.createPointList(
				vals=self.xvals,npnts=nlongs,dpnts=None,pnts=None,intrvl=longIntrvl)
			if longPnts is None :
				longPnts=self.xvals

			xvals=longPnts

			nfuns, xs, ys, zs=myInterpol.get_values_on_grid(
				xs=longPnts,
				ys=vertPnts,
				zs=horPnts,
				nx=None,
				ny=None,
				nz=None,
				)
			nfuns=nfuns.reshape((len(xs),len(ys),len(zs),nfuns.shape[-1]))

			line_intgrls=grid.grid_line_integrals_calc(
				funs=nfuns, 
				xvals=longPnts,
				yvals=vertPnts,
				zvals=horPnts,
				limit=intrgnt_limit,
				epsrel=epsrel,
				epsabs=epsabs,
				processes=processes,
				)

			potFun=np.zeros( (len(xvals),len(line_intgrls),2) )

			for ind,splineVec in enumerate(line_intgrls) : 
				vertSpline=splineVec[1]
				horSpline=splineVec[2]
				potFun[:,ind,0]=vertSpline(xvals)
				potFun[:,ind,1]=horSpline(xvals)

			integralsByBzSq=potFun.reshape( (len(xvals),len(vertPnts),len(horPnts),2) )
			integralsByBzSq=integralsByBzSq**2
			integralSumByBzSq=np.sum(integralsByBzSq,axis=3)
			integralSumByBzSq_t=np.zeros( (integralSumByBzSq.shape[0],integralSumByBzSq.shape[1],integralSumByBzSq.shape[2],1) )
			integralSumByBzSq_t[:,:,:,0]=integralSumByBzSq*facKicks
			potentialGrid=grid.grid_interpolator(
				funs=integralSumByBzSq_t,
				xvals=xvals,
				yvals=vertPnts,
				zvals=horPnts,
				)
			gX,gY,gZ=potentialGrid.calcGradient()

			line_intgrls_gY=grid.grid_line_integrals_calc(
				funs=gY._g_funs, 
				xvals=longPnts,
				yvals=vertPnts,
				zvals=horPnts,
				limit=intrgnt_limit,
				epsrel=epsrel,
				epsabs=epsabs,
				processes=processes,
				)

			line_intgrls_gZ=grid.grid_line_integrals_calc(
				funs=gZ._g_funs, 
				xvals=longPnts,
				yvals=vertPnts,
				zvals=horPnts,
				limit=intrgnt_limit,
				epsrel=epsrel,
				epsabs=epsabs,
				processes=processes,
				)

			line_intgrls_pot=grid.grid_line_integrals_calc(
				funs=potentialGrid._g_funs, 
				xvals=longPnts,
				yvals=vertPnts,
				zvals=horPnts,
				limit=intrgnt_limit,
				epsrel=epsrel,
				epsabs=epsabs,
				processes=processes,
				)

			kickFun=np.zeros( (1,len(vertPnts)*len(horPnts),2) )
			potFun=np.zeros( (1,len(vertPnts)*len(horPnts),1) )

			for ind,line_intgrl_gY in enumerate(line_intgrls_gY) : 
				vertSpline=line_intgrl_gY[0]
				horSpline=line_intgrls_gZ[ind][0]
				potSpline=line_intgrls_pot[ind][0]
				kickFun[0,ind,0]=vertSpline(xvals[-1])
				kickFun[0,ind,1]=horSpline(xvals[-1])
				potFun[0,ind,0]=potSpline(xvals[-1])*facPot

			kickFun=kickFun.reshape((1,len(vertPnts),len(horPnts),2))
			kickMapY=grid.grid_interpolator(
				funs=kickFun[:,:,:,0:1],
				xvals=[xvals[-1]],
				yvals=vertPnts,
				zvals=horPnts,
				)
			kickMapZ=grid.grid_interpolator(
				funs=kickFun[:,:,:,1:2],
				xvals=[xvals[-1]],
				yvals=vertPnts,
				zvals=horPnts,
				)

			potFun=potFun.reshape((1,len(vertPnts),len(horPnts),1))
			potMap=grid.grid_interpolator(
				funs=potFun,
				xvals=[xvals[-1]],
				yvals=vertPnts,
				zvals=horPnts,
				)
			return kickMapY, kickMapZ, potMap, gY,gZ

		return kickY, kickZ

	def calc_beff_grid(self,
			prd_lngth,
			nlongs=None,
			longIntrvl=None,
			nhor=None,
			dhor=None,
			horPnts=None,
			horIntrvl=None,
			nvert=None,
			dvert=None,
			vertPnts=None,
			vertIntrvl=None,
			intrgnt_limit=None, 
			n_max = 10, 
			epsrel=1e-3,
			epsabs=1e-5,
			processes=1,
			) :

		if longIntrvl is None :
			longIntrvl = [0,prd_lngth]

		data_y=grid.createPointList(
			vals=self.yvals,
			npnts=nvert,
			dpnts=dvert,
			pnts=vertPnts,
			intrvl=vertIntrvl
		)
		data_z=grid.createPointList(
			vals=self.zvals,
			npnts=nhor,
			dpnts=dhor,
			pnts=horPnts,
			intrvl=horIntrvl,
		)
		data_x=grid.createPointList(
			vals=self.xvals,
			npnts=nlongs,
			intrvl=longIntrvl,
		)

		interpolation=self.create_grid_interpolation()
		allpnts, xs, ys, zs=interpolation.get_values_on_grid(
			xs=data_x,
			ys=data_y,
			zs=data_z,
			nx=None,
			ny=None,
			nz=None,
			)

		allpnts=allpnts.reshape( (len(data_x),len(data_y),len(data_z),3) )

		beffs=beff_grid_calc(
			all_pnts=allpnts, 
			xs=xs,
			prd_lngth=prd_lngth, 
			n_max=n_max, 
			processes=processes,
			)
		beffs=beffs.reshape((1,len(data_y),len(data_z),3))

		beffGrid=grid.grid_interpolator(
				funs=beffs,
				xvals=[prd_lngth],
				yvals=data_y,
				zvals=data_z,
			)

		return beffGrid, data_y, data_z, longIntrvl

	def find_zero_in_array(self,array,xs) :
		lastby=array[0]
		yzeros=[]
		for ind, by in enumerate(array[1:]) :
			ind=ind+1
			some0=False
			if round(lastby,8) >= 0.0 :
				if round(by,8) < 0.0 :
					some0=True
			else :
				if round(by,8) >= 0.0 :
					some0=True
			if some0 :
				xzero=0.5*(xs[ind]+xs[ind-1])
				yzeros.append(xzero)
			lastby=by
		return yzeros

	def find_zero_crossings(self, nperiods) :
		byS,bzS = self.getSplines()
		xs=np.linspace( self.xvals[0], self.xvals[-1], 50*nperiods )
		bys=byS(xs)
		bzs=bzS(xs)
		yzeros=self.find_zero_in_array(array=self.bvals[:,1],xs=self.xvals)
		zzeros=self.find_zero_in_array(array=self.bvals[:,2],xs=self.xvals)
		return yzeros, zzeros

	def getSplines(self) : 
		by = CubicSpline(self.xvals , self.bvals[:,1] )
		bz = CubicSpline(self.xvals , self.bvals[:,2] )
		return by,bz
		# [ fst_int, snd_int ] = self.field_integrals( spline = cspline, xs = xs )

	def field_integrals( self) : 
		by,bz=self.getSplines()
		integr_limit =  10000
		bint1y = create_indefinite_integral_cs( spline = by, xs = self.xvals, integr_limit = integr_limit) 
		bint2y = create_indefinite_integral_cs( spline = bint1y, xs = self.xvals, integr_limit = integr_limit)
		bint1z = create_indefinite_integral_cs( spline = bz, xs = self.xvals, integr_limit = integr_limit) 
		bint2z = create_indefinite_integral_cs( spline = bint1z, xs = self.xvals, integr_limit = integr_limit)
		return bint1y, bint2y, bint1z, bint2z

	def getTotalFieldIntegrals(self) :
		bint1y, bint2y, bint1z, bint2z = self.field_integrals()
		lastX=self.xvals.max()
		return bint1y(lastX), bint2y(lastX), bint1z(lastX), bint2z(lastX)

	def getFieldIntegrals(self, colx = 'x') : 

		xvals=self.get_xvals(colx=colx)
		bvals=self.bvals
		f_spline_b=CubicSpline(xvals , bvals)

		xs = np.linspace( xvals[0], xvals, spline_pnts )
		fst_int = para['frst_int_spline']
		snd_int = para['scnd_int_spline']
		fst_int_y = fst_int(xs)
		snd_int_y = snd_int(xs)
		fst_data = pd.DataFrame( { colx : xs, coly : fst_int_y } )
		snd_data = pd.DataFrame( { colx : xs, coly : snd_int_y } )
		fst_field = self.create_field_from_data( data = fst_data, para = None )
		fst_field.make_ana()
		snd_field = self.create_field_from_data( data = snd_data, para = None )
		snd_field.make_ana()
		return fst_field, snd_field

	def get_integral_data(self, colx = 'x', coly = 'By', intrvl = None) : 
		para = self.get_para()
		numEx = len(para['extrema'][0])
		if para['numPer'] is None :
			if ana.is_odd(num = numEx) : 
				if numEx <= 4 :
					numPer = 1
				else :
					numPer = int((numEx-5)/2)
			else : 
				numPer = int((numEx)/2)
		else :
			numPer = para['numPer']
		if numPer > 0 :
			spline_pnts = numPer*300
		else :
			spline_pnts = 300

		frst_int_spline = para['frst_int_spline']
		snd_int_spline = para['scnd_int_spline']
		xs = np.linspace( self.data[colx].to_list()[0], self.data[colx].to_list()[-1], spline_pnts )

		if not (intrvl is None) : 
			if len(intrvl) >= 2 :
				ind_start = None
				ind_end = None		
				for ind,elem in enumerate(xs) :	
					if (elem >= intrvl[0]) and (ind_start is None) :
						ind_start = ind
					if (elem >= intrvl[1]) and (ind_end is None) :
						ind_end = ind
						break
				xs = xs[ind_start:ind_end]
			# max_diff = max(diff_vals)
			# min_diff = min(diff_vals)

			# intrvl = 0.2*(max_diff - min_diff)
			# intrvl_plt = [min_diff-intrvl,max_diff+intrvl]

		frst_int_data = frst_int_spline(xs)
		snd_int_data = snd_int_spline(xs)
		first_data = { 'x' : xs, 'intB1' : frst_int_data }
		second_data = { 'x' : xs, 'intB2' : snd_int_data }
		return first_data, second_data

	def integrate_fld( self, xs = [None,None], colx = 'x', coly = 'By' ) : 
		spline = self.get_para()['spline']
		data_x = self.data[colx].to_list()
		x_min = data_x[0]
		x_max = data_x[-1]
		if xs[0] is None : 
			xs[0] = x_min
		if xs[1] is None : 
			xs[1] = x_max
		numPer = self.get_para()['numPer']
		if not ( numPer is None ) :
			integr_limit =  400 * numPer
		else :
			integr_limit =  1000
		return quad( spline, xs[0], xs[-1], limit = integr_limit )[0]

	def compare_fields(self, bf, colx = 'x', coly = 'By' ) :
		numPer = self.get_para()['numPer']
		if numPer is None :
			numPer = 10
		integr_limit = 100*numPer
		[diff_spl, x_vals] = ana.calc_metric_diff( data_simu = self.data, data_meas = bf.data, integr_limit = integr_limit,colx = colx, coly = coly)
		return [diff_spl, x_vals]

	def get_field_interval(self, intrvl, colx = 'x', coly = 'By') :
		data =self.data[ (self.data['x'] >= intrvl[0]) & (self.data['x'] <= intrvl[1]) ]
		field = self.create_field_from_data( data = data, para = None )
		return field

	def plot_fld_map(self,
			bWhat, # "Bx", "By" or "Bz"
			xPos=None,
			yPos=None,
			zPos=None,
			nfig=None,
			filename=None,
			title=None,
			) : 
		try:
			bInd=None
			if bWhat == 'Bx' :
				bInd=0
				zLab="B$_x$ [T]"
			elif bWhat == 'By' :
				bInd=1
				zLab="B$_y$ [T]"
			elif bWhat == 'Bz' :
				bInd=2
				zLab="B$_z$ [T]"
			else :
				raise unduBfieldError("bfield.plot_fld_map: Wrong bWhat Parameter")
			notNone=0
			if not (xPos is None) :
				ind=self.get_pos_ind(
					colx='x',
					valx=xPos,
					)
				if ind is None :
					raise unduBfieldError(f"bfield.plot_fld_map: xPos={xPos:.4f} not found.")
				xPlotVals=self.zvals
				yPlotVals=self.yvals
				bPlotVals=self.bvals[ind,:,:,bInd]
				bPlotVals=bPlotVals.T
				xLab="z [mm]"
				yLab="y [mm]"
				notNone=notNone+1
				pos=f"x={xPos:.2f} mm"
			elif not (yPos is None) :
				ind=self.get_pos_ind(
					colx='y',
					valx=yPos,
					)
				xPlotVals=self.xvals
				yPlotVals=self.zvals
				xLab="x [mm]"
				yLab="z [mm]"
				bPlotVals=self.bvals[:,ind,:,bInd]
				notNone=notNone+1
				if ind is None :
					raise unduBfieldError(f"bfield.plot_fld_map: yPos={yPos:.4f} not found.")
				pos=f"y={yPos:.2f} mm"
			elif not (zPos is None) :
				ind=self.get_pos_ind(
					colx='z',
					valx=zPos,
					)
				xPlotVals=self.xvals
				yPlotVals=self.yvals
				xLab="x [mm]"
				yLab="y [mm]"
				bPlotVals=self.bvals[:,:,ind,bInd]
				notNone=notNone+1
				if ind is None :
					raise unduBfieldError(f"bfield.plot_fld_map: zPos={zPos:.4f} not found.")
				pos=f"z={zPos:.2f} mm"
			if (notNone < 1) or (notNone>1) :
				raise unduBfieldError("bfield.plot_fld_map: Only one of xPos, yPos or zPos must be given.")
			if title is None :
				title=f"B-Field Map at\n {pos}"

			X_data,Z_data = np.meshgrid(xPlotVals,yPlotVals,indexing='ij')

			if nfig is None :
				nfig=0
				plt.clf()

			fig = plt.figure(nfig,figsize=(2*13*cm_inch, 2*6.5*cm_inch), dpi=150)

			fig.suptitle(title, fontsize=12)

			ax = plt.gca()

			plt.tight_layout()

			cmap = plt.colormaps["plasma"]
			cmap = cmap.with_extremes(bad=cmap(0))

			pcm = ax.pcolormesh(X_data,Z_data,bPlotVals, cmap=cmap)
			cb=fig.colorbar(pcm, ax=ax, label=zLab)                    

			axCb = cb.ax
			text = axCb.yaxis.label
			axCb.yaxis.get_offset_text().set_fontsize(8)
			axCb.tick_params(axis='both', which='major', labelsize=8)
			font = matplotlib.font_manager.FontProperties(size=8)
			text.set_font_properties(font)

			ax.set_ylabel(yLab, fontsize=8)
			ax.set_xlabel(xLab, fontsize=8)
			ax.tick_params(axis='both', which='major', labelsize=8)

			if not (filename is None) :
				plt.savefig(filename , bbox_inches='tight')

			plt.draw()
			plt.ion()

		except unduBfieldError as e: 
			print(e)
			pdb.set_trace()

		return

	"""
	Grid - Things
	"""

	def create_difference_grid(self,bmap,nx=None,ny=None,nz=None) :

		selfBfieldMap=self.get_field_map_interpol(
			nx=nx,
			ny=ny,
			nz=nz,
			)
		selfBfieldRepr=selfBfieldMap.create_grid_interpolation()

		otherBfieldMap=bmap.get_field_map_interpol(
			nx=nx,
			ny=ny,
			nz=nz,
			)

		otherBfieldRepr=otherBfieldMap.create_grid_interpolation()
		diff_grid=selfBfieldRepr.create_difference_grid(
			grid=otherBfieldRepr,
			nx=nx,
			ny=ny,
			nz=nz,
			)

		return diff_grid

	def create_grid_interpolation(self) :
		"""
		create interpolation from bvals, xvals,...
		"""
		pnts=( self.xvals,self.yvals,self.zvals )
		return grid.grid_interpolator(
				funs=self.bvals,
				xvals=self.xvals,
				yvals=self.yvals,
				zvals=self.zvals,
			)

	def get_field_map_interpol(self,
			interpolator=None,
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
			if interpolator is None :
				interpolator=self.create_grid_interpolation()

			bs, xs, ys, zs=interpolator.get_values_on_grid(xs=xs,ys=ys,zs=zs,nx=nx,ny=ny,nz=nz)

			bfieldMap=self.clone_me()
			bfieldMap.xvals=xs
			bfieldMap.yvals=ys
			bfieldMap.zvals=zs

			nxmap=len(xs)
			nymap=len(ys)
			nzmap=len(zs)

			bfieldMap.bvals=bs.reshape( (nxmap,nymap,nzmap,interpolator._g_nvals) )

			return bfieldMap

		except unduBfieldError as e :
			print(e)

	def plot_fld(self, 
			colx = 'x', 
			coly = 'By', 
			title=None, 
			plt_extrm = False, 
			folder = '', 
			add = '',
			nfig=None,
			plot=True,
			save=False,
			filename=None
			) : 

		if nfig is None :
			fig = plt.figure(figsize=(13*cm_inch, 6.5*cm_inch), dpi=150)
		else :
			fig = plt.figure(nfig,figsize=(13*cm_inch, 6.5*cm_inch), dpi=150)
		if colx == 'x' :
			xvals=self.xvals
		elif colx == 'y' :
			xvals=self.yvals
		if colx == 'z' :
			xvals=self.zvals

		plt.plot(xvals,self.bvals, label = add)

		if title is None :
			fig.suptitle(f"{coly}-Field", fontsize=14)
		else :
			fig.suptitle(title, fontsize=14)
		plt.axhline(y=0.0,color = 'black')
		plt.xlabel(f'x [{self._unitsXB[0]}*m]', fontsize=12),
		plt.ylabel(f'{coly} [{self._unitsXB[1]}*T]', fontsize=12)
		ax = plt.gca()
		ax.legend(loc='upper right')
		if save :
			if not (filename is None) :
				plt.savefig(folder+filename+".png", bbox_inches='tight')
		if plot :
			# plt.draw()
			# plt.ion()
			plt.show()
		# if save :
		# 	plt.clf()
		return fig

	def calc_regression_B(self,nfig = None,color='black',linestyle='--', nperiods = None, add = '',intrvl = None, plot_it=True) : 

		if nperiods is None :
			nperiods = self.get_para()['numPer']

		if not (intrvl is None) :
			tmp_fld = self.get_field_interval(intrvl=intrvl, colx = 'x', coly = 'By')
			data_tmp = tmp_fld.data
		else : 
			data_tmp = self.data

		m,y0, avrg = ana.calculate_least_square(data=self.data)

		xs = np.linspace( data_tmp['x'].to_list()[0], data_tmp['x'].to_list()[-1], 300*nperiods )
		ys = [ x*m+y0 for x in xs ]

		if plot_it :
			if nfig is None :
				plt.clf()
				fig = plt.figure(figsize=(13*cm_inch, 6.5*cm_inch), dpi=150)
			else :
				fig = plt.figure(nfig,figsize=(13*cm_inch, 6.5*cm_inch), dpi=150)

			ax = plt.gca()
			plt.plot(xs,ys, label = f'{add} Fit:\n x*{m:.2E}+{y0:.2E};Av:{avrg:.2E}',color=color,linestyle=linestyle)
			box1 = ax.get_position()
			ax.set_position([box1.x0, box1.y0, box1.width * 0.8, box1.height])
			ax.legend(loc='center left', bbox_to_anchor=(1.0, 0.5))
			plt.savefig(f"lsq_regr_{add}.png", bbox_inches='tight')
			plt.draw()

		return m, y0, avrg

	def plot_diffs(self, bf, colx = 'x', coly = 'By', xlims = [],folder = '', add = '',fun=None)  : 

		# self.make_ana()
		# bf.make_ana()
		[diff, xs] = self.compare_fields(bf = bf, colx = 'x', coly = 'By' )

		if len(xlims) >= 2 :
			ind_start = None
			ind_end = None		
			for ind,elem in enumerate(xs) :	
				if (elem >= xlims[0]) and (ind_start is None) :
					ind_start = ind
				if (elem >= xlims[1]) and (ind_end is None) :
					ind_end = ind
					break
			diff_vals = diff(xs[ind_start:ind_end])
		else : 
			diff_vals = diff(xs)
			xvals = self.data[colx].to_list()
			xlims = [xvals[0],xvals[-1]]
		max_diff = max(diff_vals)
		min_diff = min(diff_vals)

		intrvl = 0.2*(max_diff - min_diff)
		intrvl_plt = [min_diff-intrvl,max_diff+intrvl]

		# intgrl_prdc = quad( diff, xlims[0], xlims[-1], limit = 300 )[0] /( xlims[-1] - xlims[0] )

		fig = plt.figure(figsize=(2*6.5*cm_inch, 6.5*cm_inch), dpi=150)
		fig.suptitle(f"Differences in {xlims}", fontsize=14)
		plt.plot(xs, diff(xs),label = 'Diff',color='black')
		plt.xlabel('x [mm]', fontsize=12),
		ax = plt.gca()
		ax.set_xlim(xlims)
		ax.set_ylim(intrvl_plt)
		plt.axvline(x=-20,color='red')
		plt.axvline(x=20,color='red')
		ax.set_ylabel('Abs. Loc. Difference', fontsize=12)

		ax2 = ax.twinx() 
		ax2.plot(self.data['x'], self.data['By'], label = f"Field",color='red',linestyle='-.')

		if not (fun is None) : 
			fun(nfig=0)

		ax2.tick_params(axis ='y', labelcolor = 'red')
		ax2.set_ylabel('B [a.u.]', fontsize=12,color='black')
		box1 = ax.get_position()
		ax.set_position([box1.x0, box1.y0, box1.width * 0.8, box1.height])
		ax.legend(loc='center left', bbox_to_anchor=(1.5, 0.5))
		box2 = ax2.get_position()
		ax2.set_position([box2.x0, box2.y0, box2.width * 0.8, box2.height])
		ax2.legend(loc='center left', bbox_to_anchor=(1.5, 0.3))

		if self.get_para()['file'] is None :
			plt.savefig(folder+f"diff_zoom_{xlims}"+add+".png", bbox_inches='tight')
		else :
			plt.savefig(folder+self.get_para()['file']+f"_diff_zoom_{xlims}"+add+".png", bbox_inches='tight')

		fig = plt.figure(figsize=(6.5*cm_inch, 6.5*cm_inch), dpi=150)
		fig.suptitle(f"Fields in {xlims}", fontsize=14)
		plt.plot(self.data['x'], self.data['By'], label = f"Original",color='black')
		ax = plt.gca()
		ax.plot(bf.data['x'], bf.data['By'], label = f"Minimized", linestyle = 'dashed',color='red')
		if not (fun is None) : 
			fun(nfig=0)
		ax.set_xlim(xlims)
		plt.axvline(x=-20,color='red')
		plt.axvline(x=20,color='red')
		ax.set_ylabel('B [a.u.]', fontsize=12)
		box1 = ax.get_position()
		ax.set_position([box1.x0, box1.y0, box1.width * 0.8, box1.height])
		ax.legend(loc='center left', bbox_to_anchor=(1.0, 0.5))

		if self.get_para()['file'] is None :
			plt.savefig(folder+f"fields_zoom_{xlims}"+add+".png", bbox_inches='tight')
		else :
			plt.savefig(folder+self.get_para()['file']+f"_fields_zoom_{xlims}"+add+".png", bbox_inches='tight')

	def find_best_gauge_factor(self,bf) : 
		self.make_ana()
		bf.make_ana()

		val0 =  self.data['By'].max()
		bfmax_g = self.gauge_b_field_data(col='By', gauge_fac=1/val0 )
		val =  bf.data['By'].max()
		bfmin_g = bf.gauge_b_field_data(col='By', gauge_fac=1/val )

		facs = np.linspace(0.8,1.2,15)

		res = []
		for fac in facs :
			sml_tmp = bfmin_g.gauge_b_field_data(col='By', gauge_fac=fac )
			xlims = bfmax_g.data['x'].to_list()
			[diff, xs] = bfmax_g.compare_fields(bf = sml_tmp, colx = 'x', coly = 'By' )
			intgrl = quad( diff, xlims[0], xlims[-1], limit = 300 )[0] /( xlims[-1] - xlims[0] )
			res.append( {'fac':fac,'intgrl':intgrl} )
		res = pd.DataFrame(res)

		csfac = CubicSpline(res['fac'],res['intgrl'])
		simu_x_lin = np.linspace( res['fac'].to_list()[0], res['fac'].to_list()[-1], num=1000)
		cs_vals = csfac(simu_x_lin)
		min_ind = min( (v, i) for i, v in enumerate(cs_vals) )[1]
		min_fac = simu_x_lin[min_ind]

		fig = plt.figure(num=0,figsize=(2*6.5*cm_inch, 6.5*cm_inch), dpi=150)
		fig.suptitle(f"Diff. vs. Factor - min {min_fac:.4f}", fontsize=14)
		plt.plot(res['fac'], res['intgrl'],marker='x',color='black')
		plt.plot(simu_x_lin, cs_vals,color='blue')

		plt.axvline(x=min_fac,color='red')
		plt.xlabel('factor', fontsize=12)
		ax = plt.gca()
		ax.set_ylabel('Abs. Loc. Difference', fontsize=12)
		plt.savefig(f"bst_fac.png", bbox_inches='tight')
		plt.clf()

		return min_fac*val0/val

	def plot_intgrls(self, colx = 'x', coly = 'By',folder='', add = '', nfig = None,save=False,plot=True,title=None,filename=None)  : 
		xs = self.data[colx].to_list()
		xslin = np.linspace( xs[0], xs[-1], num=10000)

		if nfig is None :
			plt.clf()
			fig = plt.figure(figsize=(13*cm_inch, 6.5*cm_inch), dpi=150)
		else :
			fig = plt.figure(nfig,figsize=(13*cm_inch, 6.5*cm_inch), dpi=150)
		if coly == 'By' :
			fst_int_spline = self.get_para()['frst_int_spline']
			scnd_int_spline = self.get_para()['scnd_int_spline']
		elif coly == 'Bz' :
			fst_int_spline = self.get_para()['frst_int_spline_z']
			scnd_int_spline = self.get_para()['scnd_int_spline_z']
		res_1 = fst_int_spline(xslin)
		res_2 = scnd_int_spline(xslin)
		maxv = max(res_1)
		minv = abs(min(res_1))
		val1 = maxv
		if maxv < minv : 
			val1 = minv

		maxv = max(res_2)
		minv = abs(min(res_2))
		val2 = maxv
		if maxv < minv : 
			val2 = minv

		plt.plot(xslin,res_1, label = 'First Integral')
		ax = plt.gca()
		ax2 = ax.twinx()
		ax2.plot(xslin,res_2, label = 'Second Integral',color='r')

		# plt.axhline(y=0.0, color='r', linestyle='-')
		fig.suptitle(f"First and Second Integral {add}", fontsize=14)
		ax.set_xlabel('x [mm]', fontsize=12),
		ax.set_ylabel('First Integral [Tmm]', fontsize=12)
		# ax.set_ylim([-val1,val1])
		ax2.set_ylabel('Second Integral [T$mm^2$]', fontsize=12,color='red')
		ax2.tick_params(axis ='y', labelcolor = 'red') 
		# ax2.set_ylim([-val2,val2])
		# ax.set_xlim([-150,150])
		# ax.legend(loc='best', bbox_to_anchor=(0.4, 0.4, 0.0, 0.0))
		# ax.axhline(y=0.0,color = 'black')
		ax2.axhline(y=0.0,color = 'black',linestyle='-.',linewidth=1)
		# ax2.legend(loc='best', bbox_to_anchor=(0.4, 0.3, 0.0, 0.0))
		# plt.xlim( [-0.3,0.3] )
		if save:
			if self.get_para()['file'] is None :
				plt.savefig(folder+"field_integral_"+add+".png", bbox_inches='tight')
			else :
				if filename is None :
					plt.savefig(folder+self.get_para()['file']+"_field_integral_"+add+".png", bbox_inches='tight')
				else :
					plt.savefig(folder+filename+".png", bbox_inches='tight')
		print("Field integrals: 1st: ", res_1[-1], " and 2nd: ", res_2[-1])
		if plot :
			plt.draw()
		return [ res_1, res_2 ]

	def plot_zeros(self, colx = 'x', coly = 'By')  : 
		# # # Plotting Zero Crossings in Simu and Meas
		fig = plt.figure(figsize=(6.5*cm_inch, 6.5*cm_inch), dpi=150)
		fig.suptitle("Zero Crossings", fontsize=14)
		zeros = self.get_para()['zeros']
		plt.plot(self.data[colx],self.data[coly], label = 'Simu')
		plt.plot(zeros, [0 for i in [*range(0,len(zeros))] ], "x", label = 'Zeros Simu')
		plt.axhline(y=0.0, color='r', linestyle='-')
		plt.xlabel('x [m]', fontsize=12),
		plt.ylabel('B [a.u.]', fontsize=12)
		if self.get_para()['file'] is None :
			plt.savefig("zero_crossings.png", bbox_inches='tight')
		else :
			plt.savefig(self.get_para()['file']+"_zero_crossings.png", bbox_inches='tight')
		plt.draw()

	def plot_rel_diffs(self)  : 
		# # # # Plotting relative differences between neighbouring Extrema
		fig = plt.figure(figsize=(13*cm_inch, 6.5*cm_inch), dpi=150)
		para = self.get_para()
		xmin = para['minima'][0]
		xmax = para['maxima'][0]
		diffs_min = para['diffs_min']
		diffs_max = para['diffs_max']
		plt.plot(xmin, diffs_min, "-x", linewidth=2, label =  "Minima" )
		plt.plot(xmax, diffs_max, "-x", linewidth=2, label = "Maxima" )
		fig.suptitle("Relative Peak Variations", fontsize=14)
		plt.xlabel('x [m]', fontsize=12),
		plt.ylabel('Relative Variation in %', fontsize=12)
		ax = plt.gca()
		ax.set_yscale('log')
		ax.legend(loc='best', bbox_to_anchor=(0.5, 0.5, 0.5, 0.5))
		if self.get_para()['file'] is None :
			plt.savefig("relative_extrema_differences.png", bbox_inches='tight')
		else :
			plt.savefig(self.get_para()['file']+"_relative_extrema_differences.png", bbox_inches='tight')
		plt.draw()
