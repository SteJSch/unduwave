"""
Contains the functionality for loading and processing b-field data
"""
from unduwave.unduwave_incl import *
import unduwave.helpers.file_folder_helpers as ff_h
import unduwave.quantities.quantities as quantities
from unduwave.attribute_classes.attributes import _attribute
from unduwave.attribute_classes.attributes import _attribute_collection

class bfield() : 
	"""
	Holds the bfield data
	"""

	def __init__(self, 
			unitsXB=[0.001,1.0],
			) :
		super().__init__()
		self._unitsXB=unitsXB
		self._data=np.ones((0,0))

		self.bvals=np.array([])
		self.xvals=np.array([0.0])
		self.zvals=np.array([0.0])
		self.yvals=np.array([0.0])

	"""
	bfield class

	has x,y,z,bx,by,bz

	if has x: 
		find new x before, with right interval
		find new x after

		loop new x+oldx+newx make harm bm, add to bm already there

		loop new xs :
			add 0 for all other bs
	"""

	def write_field(self,
			file,
			cols,
			outType='std',
			unitsXB=None,
			filez=None,
			filex=None,
			) :
		if outType == 'std' : 
			if field == 'map' :
				pass
			else :
				self.write_field_std(file=file,unitsXB=unitsXB)
		elif outType == 'waveBy' : 
			self.write_field_std(file=file,unitsXB=unitsXB)
		elif outType == 'waveByz' : 
			self.write_field_waveByz(filey=file,filez=filez,unitsXB=None)
		elif outType == 'waveBxyz' : 
			self.write_field_waveBxyz(filey=file,filex=filex,filez=filez,unitsXB=None)
		elif outType == 'unduOut' : 
			self.write_field_unduOut(file=file,unitsXB=unitsXB)
		elif outType == 'mapWave' : 
			pass
		# convert_x_mm_b_T_file_to_wave_std(file_in, out_path )

	def write_field_unduOut(self,file,unitsXB=None,colx='x') :
		"""
		UndumagOut file Units are mm for the length and T for B
		"""
		if unitsXB is None :
			unitsXB=self._unitsXB
		unitConvX=unitsXB[0]/0.001
		unitConvB=unitsXB[1]/1
		x_vals=self.get_xvals(colx=colx)
		bvals = self.bvals
		with open( file, 'w') as o_f:
			for ind, xval in enumerate(x_vals) :
				o_f.write(f"  {xval*unitConvX}   {bvals[ind,1]*unitConvB}   {bvals[ind,2]*unitConvB}   0.0000000E+000   0.0000000E+000   0.0000000E+000   0.0000000E+000          0\n")		

	def write_field_std(self,file,unitsConv=None, whatStr='', colx='x') :
		"""
		UndumagOut file Units are mm for the length and T for B
		"""
		if unitsConv is None :
			unitsConv=self._unitsXB
		unitConvX=unitsConv[0]/0.001 # x in mm
		unitConvB=unitsConv[1]
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

	def write_field_map_wave(self,file,unitsXB=None) :
		if unitsXB is None :
			unitsXB=self._unitsXB
		unitConvX=unitsXB[0]/1 # xyz in wave field map are in m
		unitConvB=unitsXB[1]/1

		x_vals=self.xvals
		y_vals=self.yvals
		z_vals=self.zvals
		bvals = self.bvals
		shape=x_vals.shape

		nx=shape[0]
		ny=shape[1]
		nz=shape[2]
		nlines=nx*ny*nz
		xn=x_vals.reshape(nlines)
		yn=y_vals.reshape(nlines)
		zn=z_vals.reshape(nlines)
		# bComp=bvals.reshape(nlines,3)
		with open( file, 'w') as o_f:

			o_f.write('! WAVE: x y z Bx By Bz with x as long. beam axis\n')
			o_f.write('@ date (yyyy.month.day) and time =  2025.06.25 18:59:58\n')
			o_f.write('@ run =           1\n')
			o_f.write('@ comment = WAVE.EXAMPLE\n')
			o_f.write('@ scaling = 1.0 1.0 1.0 1.0 1.0 1.0\n')
			o_f.write('@ offset = 0.0, 0.0, 0.0 0.0 0.0 0.0\n')
			for indx, xval in enumerate(np.unique(x_vals)) :
				for indy, yval in enumerate(np.unique(y_vals)) :
					for indz, zval in enumerate(np.unique(z_vals)) :
						o_f.write(f"  {xval*unitConvX}  {yval*unitConvX}  {zval*unitConvX}  {bvals[indx,indy,indz,0]*unitConvB}  {bvals[indx,indy,indz,1]*unitConvB}  {bvals[indx,indy,indz,2]*unitConvB}\n")		

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
			unduFile=False :
			unduFile=True :
				Undumag field-map file is loaded
				general 3D map (x,y,z) with full B-components
		"""

		if unduFile and (not fieldMap) :
			data = pd.read_csv( file, dtype=object, sep='\\s+', header = header)
			data.columns = ['x','By','Bz','intBy','intBz','int2By','int2Bz','quark']
			for col in data.columns:
				data[col] = data[col].astype(float)

			self.xvals=np.array(data['x'].to_list())
			nvals=len(self.xvals)
			self.bvals=np.zeros((nvals,3))
			self.bvals[:,1]=data['By'].to_list()
			self.bvals[:,2]=data['Bz'].to_list()
			return
		"""
		field map
		"""
		if fieldMap :
			if unduFile : 
				cols=[ 'imoth', 'imag', 'mat', 'ityp', 'matmod', 'x', 'y', 'z', 'Bx', 'By', 'Bz', 'B', 'Hx', 'Hy', 'Hz', 'H', 'Mx', 'My', 'Mz', 'M', 'BxDip', 'ByDip', 'BzDip', 'ifail', 'kfail', 'cmag', 'cmoth' ]
				skiprows=range(0, 3)
			elif cols is None :
				print("bfield: load_field_from_file: cols has to be set for general map file")
				return
			data = pd.read_csv( file, skiprows=skiprows, dtype=object, delimiter=r"\s+",header=None)
			data.columns = cols
			cols_float = ['x', 'y', 'z', 'Bx', 'By', 'Bz']
			for col in cols_float:
				data[col] = data[col].astype(float)

			unique_x=data['x'].unique()
			unique_y=data['y'].unique()
			unique_z=data['z'].unique()

			xGrid,yGrid,zGrid = np.meshgrid(
				unique_x,
				unique_y,
				unique_z,
				indexing='ij'
				)
			nxmap=len(unique_x)
			nymap=len(unique_y)
			nzmap=len(unique_z)

			self.bvals=np.zeros((nxmap,nymap,nzmap,3))
			bxGrid=np.array(data['Bx'].to_list()).reshape((nxmap,nymap,nzmap))
			byGrid=np.array(data['By'].to_list()).reshape((nxmap,nymap,nzmap))
			bzGrid=np.array(data['Bz'].to_list()).reshape((nxmap,nymap,nzmap))
			
			self.bvals[:,:,:,0]=bxGrid
			self.bvals[:,:,:,1]=byGrid
			self.bvals[:,:,:,2]=bzGrid

			self.xvals=xGrid
			self.yvals=yGrid
			self.zvals=zGrid
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
			self.xvals=np.array(data['x'].to_list())
			nvals=len(self.xvals)
		elif 'y' in cols :
			self.yvals=np.array(data['y'].to_list())
			nvals=len(self.yvals)
		elif 'z' in cols :
			self.zvals=np.array(data['z'].to_list())
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

	def calc_beff(self,prd_lngth, n_max = 10,colx = 'x') : 
		"""
		This Function only works for very special use cases! In general we would have to add full fourier analysis
		and determination of the phase of the reconstructed function..
		"""

		# fig = plt.figure(figsize=(13*cm_inch, 6.5*cm_inch), dpi=150)

		numPer = 4
		if not ( numPer is None ) :
			integr_limit =  400 * numPer
		else :
			integr_limit =  1000

		data_x = self.get_xvals(colx=colx)
		data_y = self.bvals[:,1]
		data_z = self.bvals[:,2]
		datas=[data_y,data_z]
		beffs=[]

		for dataB in datas :
			in_spline = CubicSpline(data_x , dataB)
			xs = np.linspace( -1/2*prd_lngth, 1/2*prd_lngth, 1000 )
			ys = in_spline(xs)
			# plt.plot(xs,ys, label = 'Original Data Cut')
			# plt.plot(data_x,dataB, label = 'Original Data-By')
			b0=quad( in_spline, -1/2*prd_lngth, 1/2*prd_lngth, limit = integr_limit )[0]
			bs=[]
			bc=[]
			for i in range(1,n_max+1) :

				xs_sin = np.sin(2*math.pi*xs/prd_lngth*i) #+math.pi/2.0
				xs_cos = np.cos(2*math.pi*xs/prd_lngth*i) #+math.pi/2.0
				new_data_s = []
				new_data_c = []
				for ind_by, by in enumerate(ys) : 
					new_val_s = 2/prd_lngth*by*xs_sin[ind_by]
					new_val_c = 2/prd_lngth*by*xs_cos[ind_by]
					new_data_s.append(new_val_s)
					new_data_c.append(new_val_c)
				# pdb.set_trace()
				# plt.plot(xs,new_data, label = f'{i}')

				cspline_sine = CubicSpline(xs , new_data_s)
				cspline_cos = CubicSpline(xs , new_data_c)
				int_val_s=quad( cspline_sine, xs[0], xs[-1], limit = integr_limit )[0]
				int_val_c=quad( cspline_cos, xs[0], xs[-1], limit = integr_limit )[0]
				bs.append(int_val_s)
				bc.append(int_val_c)
			beff = b0**2
			for ind,tmp in enumerate(bs) :
				beff = beff + (tmp**2+bc[ind]**2)/(ind+1)**2
			beff=math.sqrt(beff)
			beffs.append(beff)

		# re-calc fourier data
		# re_data = []
		# eff_field = []
		# for ind_x, x_val in enumerate(xs) : 
		# 	val = b0
		# 	for i in range(1,n_max+1) :
		# 		ampl=math.sqrt(bs[i-1]**2+bc[i-1]**2)
		# 		phs=math.atan2(bs[i-1],bc[i-1])
		# 		val = val + ampl*np.cos((i)*2*math.pi/prd_lngth*x_val+phs)
		# 	re_data.append(val)
		# 	eff_field.append( beff*np.sin(2*math.pi/prd_lngth*x_val+phs) )
		# plt.plot(xs,re_data, label = f'Reconstructed from Fourier')
		# plt.plot(xs,eff_field, label = f'Effective Field')

		# fig.suptitle("Fourier-Coefficients Bn", fontsize=14)
		# plt.xlabel('x [mm]', fontsize=12),
		# plt.ylabel('B [T]', fontsize=12)
		# ax = plt.gca()
		# ax.legend(loc='best')
		# plt.ion()
		# if coly == 'Bz' :
		# 	pdb.set_trace()
		beffFull=math.sqrt(beffs[0]**2+beffs[1]**2)
		return beffFull, beffs

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

	def field_integrals( self,xs, spline = None ) : 
		integr_limit =  10000
		cs_1st_int = ana.create_indefinite_integral_cs( spline = spline, xs = xs, integr_limit = integr_limit) 
		cs_2nd_int = ana.create_indefinite_integral_cs( spline = cs_1st_int, xs = xs, integr_limit = integr_limit)
		return [ cs_1st_int, cs_2nd_int ]

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
