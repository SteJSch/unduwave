"""
Contains the functionality for loading and processing b-field data
"""
from unduwave.unduwave_incl import *
import unduwave.helpers.file_folder_helpers as ff_h
import unduwave.quantities.quantities as quantities
from unduwave.attribute_classes.attributes import _attribute
from unduwave.attribute_classes.attributes import _attribute_collection

def convert_x_mm_b_T_file_to_wave_std( folder_in, file_in, out_path ) : 
	"""
	Loads a file in folder_in called file_in with two cols: x[mm] and B[T] - no header to separator
	and converts, depending on b_type, to wave std and copies to out_path (path+filename)
	"""
	data = pd.read_csv( folder_in+file_in, dtype=object, header = None,delimiter=r"\s+" )
	data.columns = [ "x", "B" ]
	for col in data.columns:
		data[col] = data[col].astype(float)

	data['x'] = data['x'].apply( lambda x : x*1e-3 )
	data.to_csv(out_path,index=False, header = None, sep = ' ')

	lines = []
	with open(out_path, 'r') as o_f:
		# read an store all lines into list
		lines = o_f.readlines()
	lines.insert(0,str(len(data))+'\n')
	lines.insert(0,'1.0 1.0\n')
	lines.insert(0,'Comment\n')
	with open( out_path, 'w') as o_f:
		for ind, line in enumerate(lines) :
			o_f.write(line)

def interpolate_bfield_data(b_fields, val, colInt, coly = 'By',num_supports=1000) :

	xmin = b_fields[0]['data']['x'].min()
	xmax = b_fields[0]['data']['x'].max()
	for field in b_fields : 
		if field['data']['x'].min() < xmin :
			xmin = field['data']['x'].min()
		if field['data']['x'].max() > xmax :
			xmax = field['data']['x'].max()
	support_x = np.linspace( xmin, xmax, num= int((xmax-xmin)+1) )
	interp_data = []
	for field in b_fields : 
		x_vals = field['data']['x'].to_list()
		y_vals = field['data'][coly].to_list()
		spline = CubicSpline(x_vals, y_vals)
		spline_supp_vals = spline( support_x )
		interp_data.append( { colInt : field[colInt], 'spline' : spline, 'supp_vals' : spline_supp_vals } )

	supp_list = []
	for ind, supp in enumerate(support_x) :
		supp_list.append( { 'ind' : ind, 'x' : supp, 'gB' : [] } )
		for interp in interp_data : 
			supp_list[-1]['gB'].append( { colInt : interp[colInt], 'Bg' : interp['supp_vals'][ind] } )

	num_gaps = len(b_fields)
	num_gaps_spline = num_gaps * 10
	for supp in supp_list : 
		g_data = pd.DataFrame(supp['gB'])
		g_vals = g_data[colInt].to_list()
		b_vals = g_data['Bg'].to_list()
		g_spline = CubicSpline(g_vals, b_vals)
		supp.update( { 'g_spline' : g_spline } )

	data_for_g = []
	for supp in supp_list : 
		data_for_g.append( { 'x' : supp['x'], coly : supp['g_spline'](val) } )
	data_for_g = pd.DataFrame(data_for_g)
	g_spline = CubicSpline( support_x, data_for_g[coly] )

	y_vals = g_spline(support_x)
	data_df = pd.DataFrame({'x':support_x,coly:y_vals})
	field = bfield(data=data_df)
	ret = { 
		colInt : val, 
		'bfield' : field,
		'data' : data_df
	} 
	return ret

def load_b_fields(folder, hints = [], undu_file = False) : 
	"""
	Load field files from "folder" containing files with fields for different gaps
	The file format should be two rows, the first one 'x' in [mm], the second one 'By' [T]
	The file-name format should be: file_name_front + 'gap_' + str(gap) + '_' + file_name_back + '.' + file_ending, eg: foo_gap_4.5_v1.csv
	Returns list of dictionaries: { 'gap' : gap,'data' : pd.DataFrame(data), 'file_name' : file_name }, where data is a panda dataframe with rows
	"x" and "By"
	"""
	data_ret = []
	files = f_h.find_files_exptn(folder = folder, hints = hints, exptns = [])
	fields = []
	for file in files :
		bf = bfield(file = folder + '/' + file, undu_file = undu_file)
		fields.append({ 'field' : bf, 'gap' : bf.get_para()['gap'] })
	sort_fields = make_bf_para_list( bfields = fields, para_name = 'gap'  )
	return ana.sort_data_list( data_l = sort_fields, key = 'gap'  )

def load_undu_b_fields(folder, hints = [], exptns = []) : 
	data_ret = []
	files = f_h.find_files_exptn(folder = folder, hints = hints, exptns = exptns)
	fields = []
	for file in files :
		bf = bfield(file = folder + '/' + file, undu_file = True)
		fields.append(bf)
	return fields

def make_bf_para_list( bfields, para_name = 'gap'  ) :
	"""
	Creates a dictionary containing the bfields under 'field' and the parameter para_name under the key para_name
	para_name needs to be defined in para of bf's!
	"""
	fields = []
	for field in bfields :
		if field is dict :
			fields.append({ 'field' : field['field'], para_name : field['field'].get_para()[para_name] })
		else :
			fields.append({ 'field' : field, para_name : field.get_para()[para_name] })
	return fields

def plot_b_field_data( b_fields, col_names = 'gap', col_val_plt = [] ) :
	"""
	Plots all the fields in list b_fields, if col_names and col_val_plt are given, then only the fields which parameter
	col_names is equal to some value in col_val_plt list.
	e.g. : plot all fields with gap 15 or 20
	"""
	fig = plt.figure(figsize=(13*cm_inch, 6.5*cm_inch), dpi=150)
	for b_field in b_fields :
		b_field = b_field['field']
		plt_it = True
		if len(col_val_plt) > 0 :
			plt_it = False
			for val in col_val_plt : 
				if b_field.get_para()[col_names] == val :
					plt_it = True
					break
		if plt_it :
			b_field.plot_data(colX = 'x', colY = 'By', fig = fig, label_para = col_names)
	fig.suptitle("Magnetic Field", fontsize=14)
	plt.xlabel('x [mm]', fontsize=12),
	plt.ylabel('Magnetic Field [T]', fontsize=12)
	ax = plt.gca()
	ax.legend(loc='best', bbox_to_anchor=(0.8, 0.5, 0.0, 0.0))	
	plt.draw()

def interpolate_b_data(b_fields, gap, lim_peak, num_support_per_extrema = 20, colx = 'x', coly = 'By') :
	"""
	interpolates b-field data for a given gap using already present dataframes for different gaps
	takes b-field data list loaded with load_b_fields_gap, a gap number, the lim_peak value 
	used for findind the extrema in the data (for determination of the number of periods)
	and the number of support positions per extrema for the calculation of splines from the data
	Returns a list containing a dictionary: { 'gap' : gap,'data' : pd.DataFrame(data), 'file_name' : file_name }, where data is a panda dataframe 
	containing the interpolated data and file_name contains the value of the gap at which this data is determined
	"""
	fields_cut = []
	for field in b_fields :
		field = field['field']
		[field_cut, x_min, x_max] = field.cut_data_support(col_cut = colx)
		fields_cut.append(field_cut)

	num_ex = []
	for field in fields_cut : 
		[ex_x, ex_y] = ana.find_extrema(data = field.data, lim = lim_peak, colx = colx, coly = coly)
		num_ex.append( len(ex_x) )
	max_num_ex = max(num_ex)
	num_supports = max_num_ex * num_support_per_extrema
	support_x = np.linspace( x_min, x_max, num=num_supports)
	interp_data = []
	for field in fields_cut : 
		x_vals = field.data[colx].to_list()
		y_vals = field.data[coly].to_list()
		spline = CubicSpline(x_vals, y_vals)
		spline_supp_vals = spline( support_x )
		interp_data.append( { 'gap' : field.get_para()['gap'], 'spline' : spline, 'supp_vals' : spline_supp_vals } )

	supp_list = []
	for ind, supp in enumerate(support_x) :
		supp_list.append( { 'ind' : ind, 'x' : supp, 'gB' : [] } )
		for interp in interp_data : 
			supp_list[-1]['gB'].append( { 'gap' : interp['gap'], 'Bg' : interp['supp_vals'][ind] } )

	num_gaps = len(b_fields)
	num_gaps_spline = num_gaps * 10
	for supp in supp_list : 
		g_data = pd.DataFrame(supp['gB'])
		g_vals = g_data['gap'].to_list()
		b_vals = g_data['Bg'].to_list()
		g_spline = CubicSpline(g_vals, b_vals)
		supp.update( { 'g_spline' : g_spline } )

	data_for_g = []
	for supp in supp_list : 
		data_for_g.append( { colx : supp['x'], coly : supp['g_spline'](gap) } )
	data_for_g = pd.DataFrame(data_for_g)
	g_spline = CubicSpline( support_x, data_for_g['By'] )
	intrp_f = b_fields[-1]['field'].clone_me()
	para = intrp_f.get_para()
	para.update( { 'gap' : gap } )
	para.update( { 'file' : 'interp_b_field_gap_'+str(gap)+'_.dat' } )
	intrp_f.set_para( para )
	intrp_f.data = data_for_g
	return intrp_f

def cut_data_support(bfs, col_cut = 'x') : 
	"""
	takes b-field data list loaded with load_b_fields_gap and
	cuts the fields to the smallest common support in the column
	col_cut - afterwards all fiel-data is defined on the same col-values
	"""
	b_field_x_small = None
	b_field_x_high = None
	bfs_r = []
	for bf in bfs :
		bf = bf.clone_me()
		x_vals = bf.data[col_cut].to_list()
		if b_field_x_high is None :
			b_field_x_high = x_vals[-1]
		else : 
			if x_vals[-1] < b_field_x_high : 
				b_field_x_high = x_vals[-1]
		if b_field_x_small is None :
			b_field_x_small = x_vals[0]
		else : 
			if x_vals[0] > b_field_x_small : 
				b_field_x_small = x_vals[0]

		bf.data = bf.data[ (bf.data[col_cut] >= b_field_x_small) & \
					(bf.data[col_cut] <= b_field_x_high) ]
		bfs_r.append(bf)
	return [ bfs_r, b_field_x_small, b_field_x_high]

def create_harm_field( period, amplitude, deltaX, phase_shift = 0, num_pnts_per_period = 100, colx = 'x', coly = 'By' ) : 
	"""
	Creates a sine field: amplitude*sin( k_per * x + phase_shift ) for all x lying in interval deltaX, 
	creates bfield class filled with data, with num_pnts_per_period pnts per period in deltaX
	"""
	length = deltaX[-1] - deltaX[0]
	numPer = length / period
	xpnts = np.linspace( deltaX[0], deltaX[-1], (int)(numPer*num_pnts_per_period) )
	data = []
	for pnt in xpnts : 
		data.append( { colx : pnt, coly : amplitude*math.sin( 2*math.pi/period * pnt + phase_shift ) } )
	bf_r = bfield()
	bf_r = bf_r.create_field_from_data( data = pd.DataFrame(data) )
	return bf_r

def create_field_with_ends(period_length,amplY,nperiods,num_pnts_per_period=100,amplZ=0,shift=0,gaugeY=1.0,gaugeZ=-1.0) :
	bfield=create_harm_field( 
		period=period_length, 
		amplitude=amplY, 
		deltaX=[0,period_length*nperiods], 
		phase_shift = 0, num_pnts_per_period = num_pnts_per_period, 
		colx = 'x', coly = 'By' 
		)
	be1=create_harm_field( 
		period=period_length,
		amplitude=amplY*0.25, 
		deltaX=[0,period_length/2.0], 
		phase_shift = 0, num_pnts_per_period = num_pnts_per_period/2.0, 
		colx = 'x', coly = 'By' 
		)
	be2=create_harm_field( 
		period=period_length, 
		amplitude=amplY*0.75, 
		deltaX=[0,period_length/2.0], 
		phase_shift = math.pi, num_pnts_per_period = num_pnts_per_period/2.0, 
		colx = 'x', coly = 'By' 
		)
	be3=create_harm_field( 
		period=period_length, 
		amplitude=amplY*0.75, 
		deltaX=[0,period_length/2.0], 
		phase_shift = 0, num_pnts_per_period = num_pnts_per_period/2.0, 
		colx = 'x', coly = 'By' 
		)
	be4=create_harm_field( 
		period=period_length, 
		amplitude=amplY*0.25, 
		deltaX=[0,period_length/2.0], 
		phase_shift = math.pi, num_pnts_per_period = num_pnts_per_period/2.0, 
		colx = 'x', coly = 'By' 
		)

	bfully=be1.glue_bf( bf=be2, pos = 'back', colx = 'x', coly = 'By') 
	bfully=bfully.glue_bf( bf=bfield, pos = 'back', colx = 'x', coly = 'By') 
	bfully=bfully.glue_bf( bf=be3, pos = 'back', colx = 'x', coly = 'By') 
	bfully=bfully.glue_bf( bf=be4, pos = 'back', colx = 'x', coly = 'By') 
	bfully=bfully.get_cntrd_bf(lim_peak = None, colx = 'x', coly = 'By' )
	bfully=bfully.gauge_b_field_data(col='By',gauge_fac=gaugeY)

	if not (amplZ == 0) :
		bfieldz=create_harm_field( 
			period=period_length, 
			amplitude=amplZ, 
			deltaX=[0,period_length*nperiods], 
			phase_shift = 0, num_pnts_per_period = num_pnts_per_period, 
			colx = 'x', coly = 'Bz' 
			)
		be1z=create_harm_field( 
			period=period_length, 
			amplitude=amplZ*0.25, 
			deltaX=[0,period_length/2.0], 
			phase_shift = 0, num_pnts_per_period = num_pnts_per_period/2.0, 
			colx = 'x', coly = 'Bz' 
			)
		be2z=create_harm_field( 
			period=period_length, 
			amplitude=amplZ*0.75, 
			deltaX=[0,period_length/2.0], 
			phase_shift = math.pi, num_pnts_per_period = num_pnts_per_period/2.0, 
			colx = 'x', coly = 'Bz' 
			)
		be3z=create_harm_field( 
			period=period_length, 
			amplitude=amplZ*0.75, 
			deltaX=[0,period_length/2.0], 
			phase_shift = 0, num_pnts_per_period = num_pnts_per_period/2.0, 
			colx = 'x', coly = 'Bz' 
			)
		be4z=create_harm_field( 
			period=period_length, 
			amplitude=amplZ*0.25, 
			deltaX=[0,period_length/2.0], 
			phase_shift = math.pi, num_pnts_per_period = num_pnts_per_period/2.0, 
			colx = 'x', coly = 'Bz' 
			)

		bfullz=be1z.glue_bf( bf=be2z, pos = 'back', colx = 'x', coly = 'Bz') 
		bfullz=bfullz.glue_bf( bf=bfieldz, pos = 'back', colx = 'x', coly = 'Bz') 
		bfullz=bfullz.glue_bf( bf=be3z, pos = 'back', colx = 'x', coly = 'Bz') 
		bfullz=bfullz.glue_bf( bf=be4z, pos = 'back', colx = 'x', coly = 'Bz') 
		bfullz=bfullz.get_cntrd_bf(lim_peak = None, colx = 'x', coly = 'Bz' )
		bfullz.move(dist=shift)
		bfullz=bfullz.gauge_b_field_data(col='Bz',gauge_fac=gaugeZ)
	return bfully, bfullz

# funs to manipulate bfield lists
# interpolate, e.g., min bfield, 

class bfield(_attribute_collection) : 
	"""
	Holds the bfield data
	"""

	bx = quantities.quantity(
		api=None,
		data=None,
		name=f"B$_x$",
		plot_name=f'Bx',
		description=f'Magnetic Field Bx',
		unit='T',
		)
	by = quantities.quantity(
		api=None,
		data=None,
		name=f"B$_y$",
		plot_name=f'By',
		description=f'Magnetic Field By',
		unit='T',
		)
	bz = quantities.quantity(
		api=None,
		data=None,
		name=f"B$_z$",
		plot_name=f'Bz',
		description=f'Magnetic Field Bz',
		unit='T',
		)
	xvals = quantities.quantity(
		api=None,
		data=None,
		name=f"x",
		plot_name=f'x',
		description=f'x coordinate',
		unit='m',
		)
	yvals = quantities.quantity(
		api=None,
		data=None,
		name=f"y",
		plot_name=f'y',
		description=f'y coordinate',
		unit='m',
		)
	zvals = quantities.quantity(
		api=None,
		data=None,
		name=f"z",
		plot_name=f'z',
		description=f'z coordinate',
		unit='m',
		)

	def __init__(self, 
			data = None,
			unitsXB=[0.001,1.0],
			) :
		super().__init__()
		self._data=data
		self._unitsXB=unitsXB

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

	def write_field_unduOut(self,file,unitsXB=None) :
		"""

		UndumagOut file Units are mm for the length and T for B
		"""
		if unitsXB is None :
			unitsXB=self._unitsXB
		unitConvX=unitsXB[0]/0.001
		unitConvB=unitsXB[1]/1
		x_vals=self.xvals._data
		b_y = self.by._data
		b_z=self.bz._data
		with open( file, 'w') as o_f:
			for ind, xval in enumerate(x_vals) :
				o_f.write(f"  {xval*unitConvX:.4g}   {b_y[ind]*unitConvB:.4g}   {b_z[ind]*unitConvB:.4g}   0.0000000E+000   0.0000000E+000   0.0000000E+000   0.0000000E+000          0\n")		

	def write_field_std(self,file,unitsXB=None) :
		"""

		UndumagOut file Units are mm for the length and T for B
		"""
		if unitsXB is None :
			unitsXB=self._unitsXB
		unitConvX=unitsXB[0]/0.001
		unitConvB=unitsXB[1]/1
		x_vals=self.xvals._data
		b_y = self.by._data
		with open( file, 'w') as o_f:
			for ind, xval in enumerate(x_vals) :
				o_f.write(f"{xval*unitConvX:.4g} {b_y[ind]*unitConvB:.4g}\n")

	def write_field_waveByz(self,filey,filez,unitsXB=None) :
		"""

		UndumagOut file Units are mm for the length and T for B
		"""
		if unitsXB is None :
			unitsXB=self._unitsXB
		unitConvX=unitsXB[0]/1 # here the x-coord has to be in m?!!!
		unitConvB=unitsXB[1]/1
		x_vals=self.xvals._data
		b_y = self.by._data
		b_z = self.bz._data
		with open( filey, 'w') as o_f:
			o_f.write(str(len(x_vals))+'\n')
			o_f.write('1.0 1.0\n')
			o_f.write('Comment\n')
			for ind, xval in enumerate(x_vals) :
				o_f.write(f"{xval*unitConvX:.4g} {b_y[ind]*unitConvB:.4g}\n")
		with open( filez, 'w') as o_f:
			o_f.write(str(len(x_vals))+'\n')
			o_f.write('1.0 1.0\n')
			o_f.write('Comment\n')
			for ind, xval in enumerate(x_vals) :
				o_f.write(f"{xval*unitConvX:.4g} {b_z[ind]*unitConvB:.4g}\n")

	def write_field_waveBxyz(self,filex,filey,filez,unitsXB=None) :
		"""

		UndumagOut file Units are mm for the length and T for B
		"""
		if unitsXB is None :
			unitsXB=self._unitsXB
		unitConvX=unitsXB[0]/1 # here the x-coord has to be in m?!!!
		unitConvB=unitsXB[1]/1
		x_vals=self.xvals._data
		b_x = self.bx._data
		b_y = self.by._data
		b_z = self.bz._data
		with open( filex, 'w') as o_f:
			o_f.write(str(len(x_vals))+'\n')
			o_f.write('1.0 1.0\n')
			o_f.write('Comment\n')
			for ind, xval in enumerate(x_vals) :
				o_f.write(f"{xval*unitConvX:.4g} {b_x[ind]*unitConvB:.4g}\n")
		with open( filey, 'w') as o_f:
			o_f.write(str(len(x_vals))+'\n')
			o_f.write('1.0 1.0\n')
			o_f.write('Comment\n')
			for ind, xval in enumerate(x_vals) :
				o_f.write(f"{xval*unitConvX:.4g} {b_y[ind]*unitConvB:.4g}\n")
		with open( filez, 'w') as o_f:
			o_f.write(str(len(x_vals))+'\n')
			o_f.write('1.0 1.0\n')
			o_f.write('Comment\n')
			for ind, xval in enumerate(x_vals) :
				o_f.write(f"{xval*unitConvX:.4g} {b_z[ind]*unitConvB:.4g}\n")

	def write_field_map_wave(file) :

		with open( f"{folder}field_out.map", 'w') as o_f:

			o_f.write('! WAVE: x y z Bx By Bz with x as long. beam axis\n')
			o_f.write('@ date (yyyy.month.day) and time = 2025.06.25 18:59:58\n')
			o_f.write('@ run =           960\n')
			o_f.write('@ comment = WAVE.EXAMPLE\n')
			o_f.write('@ scaling = 1.0 1.0 1.0 1.0 1.0 1.0\n')
			o_f.write('@ offset = 0.0, 0.0, 0.0 0.0 0.0 0.0\n')
			for ind, line in enumerate(bMap) :
				o_f.write(f"  {line['x']:.2f}  {line['y']:.2f}  {line['z']:.2f}  {line['Bx']:.2f}  {line['By']:.2f}  {line['Bz']:.2f}\n")		

	def load_field_from_file(self,
			file, 
			fieldMap=False,
			cols=None,
			unduFile = False, 
			radiaFile=False,
			header=None,
			skiprows=None,
			) :

		if unduFile and (not fieldMap) :
			data = pd.read_csv( file, dtype=object, sep='\\s+', header = header)
			data.columns = ['x','By','Bz','intBy','intBz','int2By','int2Bz','quark']
			for col in data.columns:
				data[col] = data[col].astype(float)

			self.xvals._data=data['x'].to_list()
			self.by._data=data['By'].to_list()
			self.bz._data=data['Bz'].to_list()
			return
		"""
		field map
		"""
		if fieldMap :
			if unduFile : 
				cols=[ 'imoth', 'imag', 'mat', 'ityp', 'matmod', 'x', 'y', 'z', 'Bx', 'By', 'Bz', 'B', 'Hx', 'Hy', 'Hz', 'H', 'Mx', 'My', 'Mz', 'M', 'BxDip', 'ByDip', 'BzDip', 'ifail', 'kfail', 'cmag', 'cmoth' ]
				skiprows=range(0, 3)
			elif cols is None :
				print("bfield: load_field_from_file: cols has to be set for random map")
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
				)
			nxmap=len(unique_x)
			nymap=len(unique_y)
			nzmap=len(unique_z)

			bxGrid=np.array(data['Bx'].to_list()).reshape((nxmap,nymap,nzmap))
			byGrid=np.array(data['By'].to_list()).reshape((nxmap,nymap,nzmap))
			bzGrid=np.array(data['Bz'].to_list()).reshape((nxmap,nymap,nzmap))
			self.xvals._data=xGrid
			self.yvals._data=yGrid
			self.zvals._data=zGrid
			self.bx._data=bxGrid
			self.by._data=byGrid
			self.bz._data=bzGrid
			return
			
		if cols is None :
			print("bfield: load_field_from_file: cols has to be set for random file")
			return
		data = pd.read_csv( file, dtype=object, sep='\\s+', header = header, skiprows=skiprows)
		data.columns = cols
		for col in data.columns:
			data[col] = data[col].astype(float)
		if 'x' in cols :
			self.xvals._data=data['x'].to_list()
		if 'y' in cols :
			self.yvals._data=data['y'].to_list()
		if 'z' in cols :
			self.zvals._data=data['z'].to_list()
		if 'Bx' in cols :
			self.bx._data=data['Bx'].to_list()
		if 'By' in cols :
			self.by._data=data['By'].to_list()
		if 'Bz' in cols :
			self.bz._data=data['Bz'].to_list()

	def clone_me(self) : 
		bf = bfield()
		bf.data = copy.deepcopy(self.data)
		bf.set_para(self.get_para())
		return bf

	def load_field_from_radia(self, file) :
		self.load_field_from_file(file=file)

		#data in radia coordinates
		# self.data.columns = ['x','y','z','Bx','By','Bz','B','Hx','Hy','Hz','H','Mx',
		# 		'My','Mz','M']

		# undu_x->radia_y,...,y->z,z->x
		#data in undu coordinates
		self.data.columns = ['z','x','y','Bz','Bx','By','B','Hz','Hx','Hy','H','Mz',
				'Mx','My','M']

# 		'By','Bz','intBy','intBz','int2By','int2Bz','quark']
# str(x)+ " " + str(y)+ " " + str(z)+ " " \
#         + str(Bx)+ " " + str(By)+ " " + str(Bz)+ " " + str(B)+ " " \
#         + str(Hx)+ " " + str(Hy)+ " " + str(Hz)+ " " + str(H)+ " " \
#         + str(Mx)+ " " + str(My)+ " " + str(Mz)+ " " + str(M) + NL

	def get_para(self) : 
		return self.para

	def set_para(self, para) : 
		self.para = dict(para)
	
	def get_std_para(self, b_type = None ) : 
		if b_type is None : 
			b_type = self.type_On_Axis_By
		para = { 'type' : b_type, 'spline_num_pnts' : 1000, 'gap' : None, 'file' : None, 'undu_file' : False, 'extrema' : None, 'numPer' : None, \
				'minima' : None, 'maxima' : None , 'zeros' : None, 'diffs_min' : None, 'diffs_max' : None, 'spline' : None, 'frst_int_spline' : None, 'scnd_int_spline' : None, \
				'frst_int' : None, 'scnd_int' : None, 'prd_lngth' : None, 'strt_end_prdc' : None, 'frst_int_end_strct' : None,\
				'frst_int_prdc' : None }
		return para

	def save_data(self, file) : 
		self.data.to_csv( file , index = False, sep = ' ' )

	def make_child_para(self) : 
		para = self.get_para()
		new_para = self.get_std_para()
		new_para['type'] = para['type']
		new_para['spline_num_pnts'] = para['spline_num_pnts']
		new_para['gap'] = para['gap']
		new_para['file'] = para['file']
		new_para['undu_file'] = para['undu_file']
		return new_para

	def get_gap_from_filename(self, file) : 
		gap_str = file.split( 'gap_' )
		gap = None
		if len(gap_str) > 1 : 
			gap = float(gap_str[-1].split('_')[0])
		return gap

	def plot_data(self, colX = 'x', colY = 'By', fig = None, label_para = None) :
		end_s = False
		if fig is None :  
			end_s = True
			fig = plt.figure(figsize=(13*cm_inch, 6.5*cm_inch), dpi=150)			
		data = self.data
		if not (label_para is None ) :
			plt.plot(data[colX],data[colY], '-', label = label_para+'='+str(self.get_para()[label_para]))
		else :
			plt.plot(data[colX],data[colY], '-')
		if end_s :
			fig.suptitle("Magnetic Field", fontsize=14)
			plt.xlabel('x [mm]', fontsize=12)
			plt.ylabel('Magnetic Field [T]', fontsize=12)
			plt.pause(0.25)
			plt.ion()

	def get_period_length_frm_zeros(self, zeros=None) :
		if zeros is None :
			zeros=self.get_para()['zeros']
		if len(zeros) < 3 :
			return None
		ind_fst_zero = 0
		for ind, zero in enumerate(zeros) :
			if zero > 0 :
				ind_fst_zero = ind 
				break
		if ind_fst_zero + 1 >= len(zeros) :
			prd_lngth = zeros[ind_fst_zero] - zeros[ind_fst_zero-2]
		else : 
			prd_lngth = zeros[ind_fst_zero + 1] - zeros[ind_fst_zero-1]
		return [prd_lngth, zeros[ind_fst_zero]]

	def make_ana(self, lim_peak = None, colx = 'x', coly = 'By') : 
		# count max,min, extr, 1+2.int., num periods, end-structs?, 
		ret_data = []
		# data_cntrd = self.get_cntrd_bf( lim_peak, colx, coly)
		data_cntrd = self
		para = self.get_para()

		[xEx, yEx] = ana.find_extrema(data = data_cntrd.data, lim = lim_peak, colx = colx, coly = coly)
		numEx = len(xEx)
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
		para['extrema'] = [ xEx, yEx ]
		para['numPer'] = numPer

		[xMin, yMin] = ana.find_minima(data = data_cntrd.data, lim = lim_peak, colx = colx, coly = coly)
		[xMax, yMax] = ana.find_maxima(data = data_cntrd.data, lim = lim_peak, colx = colx, coly = coly)
		diffs_min = ana.rel_difference( points = yMin )
		diffs_max = ana.rel_difference( points = yMax )
		para['minima'] = [xMin, yMin]
		para['maxima'] = [xMax, yMax]
		para['diffs_min'] = diffs_min
		para['diffs_max'] = diffs_max
		try:
			cspline = CubicSpline(data_cntrd.data['x'] , data_cntrd.data['By'])
		except:
			lx=data_cntrd.data['x'].to_list()
			for inde, el in enumerate(lx):
				if inde < (len(lx)-1):
					if el >= lx[inde+1] :
						pdb.set_trace()
						print("stop")
			pdb.set_trace()
		xs = np.linspace( self.data[colx].to_list()[0], self.data[colx].to_list()[-1], spline_pnts )
		by = cspline(xs)

		if ('Bz' in data_cntrd.data.keys()) :
			csplinez = CubicSpline(data_cntrd.data['x'] , data_cntrd.data['Bz'])
			[ fst_int_z, snd_int_z ] = self.field_integrals( spline = csplinez, xs = xs )
			para.update( { 'spline_z' : csplinez } )
			para.update( { 'frst_int_spline_z' : fst_int_z } )
			para.update( { 'scnd_int_spline_z' : snd_int_z } )
			para['frst_int_z'] = fst_int_z(xs[-1])
			para['scnd_int_z'] = snd_int_z(xs[-1])

		zeros = ana.find_zero_crossings( data = pd.DataFrame( { colx : xs, coly : by } ), colx = colx, coly = coly)
		para['zeros'] = zeros

		if para['prd_lngth'] is None :
			res_prd = self.get_period_length_frm_zeros(zeros = zeros)
			if res_prd is None : 
				para['prd_lngth'] = None
			else : 
				para['prd_lngth'] = res_prd[0]
		[ fst_int, snd_int ] = self.field_integrals( spline = cspline, xs = xs )
		para['spline'] = cspline
		para['frst_int_spline'] = fst_int
		para['scnd_int_spline'] = snd_int
		para['frst_int'] = fst_int(xs[-1])
		para['scnd_int'] = snd_int(xs[-1])

		if not (para['prd_lngth'] is None) :
			para['strt_end_prdc'] = [ -numPer/2.0 * para['prd_lngth'], numPer/2.0 * para['prd_lngth'] ]
			# find the periodic and end-integrals
			fst_int_begin = fst_int( para['strt_end_prdc'][0] )
			fst_int_prdc = fst_int( para['strt_end_prdc'][1] ) - fst_int_begin
			fst_int_end = para['frst_int'] - fst_int_prdc - fst_int_begin
			fst_int_end_strct = fst_int_begin + fst_int_end
			para['frst_int_end_strct'] = fst_int_end_strct
			para['frst_int_prdc'] = fst_int_prdc

	def calc_beff(self,prd_lngth, n_max = 10,colx = 'x',coly='By') : 
		"""
		This Function only works for very special use cases! In general we would have to add full fourier analysis
		and determination of the phase of the reconstructed function..
		"""

		fig = plt.figure(figsize=(13*cm_inch, 6.5*cm_inch), dpi=150)

		numPer = 4
		if not ( numPer is None ) :
			integr_limit =  400 * numPer
		else :
			integr_limit =  1000

		data_x = self.data[colx].to_list()
		data_y = self.data[coly].to_list()
		in_spline = CubicSpline(data_x , data_y)

		xs = np.linspace( -1/2*prd_lngth, 1/2*prd_lngth, 1000 )

		ys = in_spline(xs)
		plt.plot(xs,ys, label = 'Original Data Cut')
		plt.plot(data_x,data_y, label = 'Original Data')
		# plt.show()
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

		# re-calc fourier data
		re_data = []
		eff_field = []
		for ind_x, x_val in enumerate(xs) : 
			val = b0
			for i in range(1,n_max+1) :
				ampl=math.sqrt(bs[i-1]**2+bc[i-1]**2)
				phs=math.atan2(bs[i-1],bc[i-1])
				val = val + ampl*np.cos((i)*2*math.pi/prd_lngth*x_val+phs)
			re_data.append(val)
			eff_field.append( beff*np.sin(2*math.pi/prd_lngth*x_val+phs) )
		plt.plot(xs,re_data, label = f'Reconstructed from Fourier')
		plt.plot(xs,eff_field, label = f'Effective Field')

		fig.suptitle("Fourier-Coefficients Bn", fontsize=14)
		plt.xlabel('x [mm]', fontsize=12),
		plt.ylabel('B [T]', fontsize=12)
		ax = plt.gca()
		ax.legend(loc='best')
		# plt.ion()
		# if coly == 'Bz' :
		# 	pdb.set_trace()

		return beff

	def get_field_integrals_obj(self, colx = 'x', coly = 'By') : 
		para = self.get_para()

		numPer = para['numPer']
		if not (numPer is None) :
			if numPer > 0 :
				spline_pnts = numPer*300
			else :
				spline_pnts = 300
		else: 
			spline_pnts = 3000

		xs = np.linspace( self.data[colx].to_list()[0], self.data[colx].to_list()[-1], spline_pnts )
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

	def move(self,dist) : 
		new_x = []
		for el in self.data['x'].to_list() :
			new_x.append(el+dist)
		self.data['x'] = new_x

	def gauge_b_field_data(self, col, gauge_fac ) : 
		"""
		takes a list of b-fields as returned by load_b_fields_gap and col, the name of the 
		data column of the data objects to gauge and the number gauge_fac
		each element in the column col is then multiplied by gauge_fac
		"""
		bf = self.clone_me()
		bf.data[col] = bf.data[col].apply( lambda x: x*gauge_fac )
		return bf

	def get_cntrd_bf(self, lim_peak = None, colx = 'x', coly = 'By' ) : 
		"""
		takes b-field data list loaded with load_b_fields_gap and centers each field 
		according to the position of the first and last peak identified for which |peak| >= lim_peak
		"""
		bf = self.clone_me()
		bf.data = ana.center_data(data = self.data, strat = { 'name' : 'extrema', 'lim' : lim_peak } , colx = colx, coly = coly)
		return bf

	def convert_x_mm_b_T_file_to_wave_std(self, file_in, out_path ) : 
		"""
		Loads a file in folder_in called file_in with two cols: x[mm] and B[T] - no header to separator
		and converts, depending on b_type, to wave std and copies to out_path (path+filename)
		"""
		data = pd.read_csv( file_in, dtype=object, delim_whitespace=True )
		data.columns = [ "x", "B" ]
		for col in data.columns:
			data[col] = data[col].astype(float)

		data['x'] = data['x'].apply( lambda x : x*1e-3 )
		data.to_csv(out_path,index=False, header = None, sep = ' ')

		lines = []
		with open(out_path, 'r') as o_f:
			# read an store all lines into list
			lines = o_f.readlines()
		lines.insert(0,str(len(data))+'\n')
		lines.insert(0,'1.0 1.0\n')
		lines.insert(0,'Comment\n')
		with open( out_path, 'w') as o_f:
			for ind, line in enumerate(lines) :
				o_f.write(line)

	def save_prepared_b_data(self, folder, name_add = '' ) :
		"""
		takes b-field data list loaded with load_b_fields_gap, a folder and a string name_add 
		and saves all b-field data into the folder, the files are named according ot the 
		file_name propertie in the b_fields dictionaries with name_add added to the file names
		"""
		filename = self.get_para()['file']
		filename_split = filename.split('.')
		filename_ending = filename_split[-1]
		file_name_no_ending = filename.replace( filename_ending,'' )
		filename = file_name_no_ending + name_add + filename_ending
		self.data.to_csv(folder+filename,index=False, sep = ' ')

	def normalize_data( self, strat = 'extrema', colx = 'x', coly = 'By' ) : 
		bf_nrmlzd = self.clone_me()
		bf_nrmlzd.data = ana.normalize_data( data = bf_nrmlzd.data, strat=strat, \
				colx = colx, coly = coly )
		return bf_nrmlzd

	def field_integrals( self,xs, spline = None ) : 
		numPer = self.get_para()['numPer']
		if spline is None : 
			spline = self.get_para()['spline']
		if not ( numPer is None ) :
			if numPer > 0 :
				integr_limit =  400 * numPer
			else :
				integr_limit = 400
		else :
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

	def extract_avrg_period(self, prds_to_avrg = 3, min_rel_chng = 0.01, lim_peak = None, colx = 'x', coly = 'By' ) :
		"""
		Extracts prds_to_avrg number of periods btwn zero-crossings in region where min_rel_chng is not exceeded by the 
		neighbouring extrema. Averages the extracted periods.
		Analyze-field should have been performed beforehand
		"""
		zeros = self.get_para()['zeros']
		ind_zero = 0
		# look for closest zero to 0 (on x-axis)
		for ind, zero in enumerate(zeros) :
			if zero >= 0 : 
				ind_zero = ind 
				break
		# partition the data around this pnt such that the extracted periods lie around
		# equally distributed
		numzerosfrwd = len(zeros) - ind_zero
		np_frwd = int(prds_to_avrg/2)
		np_bckwd = prds_to_avrg - np_frwd
		frwd_chng = False
		if 2 * np_frwd >= numzerosfrwd : 
			np_frwd = int(numzerosfrwd/2)
			np_bckwd = prds_to_avrg - np_frwd
			frwd_chng = True
		if 2 * np_bckwd > ind_zero : 
			np_bckwd = int(ind_zero/2)
			prds_to_avrg = np_bckwd + np_frwd
		if prds_to_avrg < 1 :
			print("Not enough periods to average.")
			return None
		print("Averaging ",prds_to_avrg," periods between x=",zeros[ind_zero-2*np_bckwd]," and x=",\
				zeros[ind_zero+2*np_frwd])
		zeros_use = zeros[ ind_zero - 2*np_bckwd : ind_zero + 2*np_frwd + 1]
		# Check accurracy of the maxima contained in this region
		diffs_min = self.get_para()['diffs_min']
		xmins = self.get_para()['minima'][0]
		start = False
		ind_intrvl_chng = []
		for ind, elem in enumerate(diffs_min) : 
			if ( elem <= min_rel_chng ) and (not start) : 
				ind_intrvl_chng.append(xmins[ind])
				start = True
			if ( ( elem >= min_rel_chng ) and start ) or ((ind+1) >= len(diffs_min)) : 
				ind_intrvl_chng.append(xmins[ind])
				break
		if len(ind_intrvl_chng) < 2 : 
			print("No region with sufficient accurracy found.")
			return None
		for zero in zeros_use :
			if ((zero < ind_intrvl_chng[0]) or (zero > ind_intrvl_chng[1])) : 
				print("accurracy condition not met.")
				return None
		ind_zero_strt = 0
		para_new = self.make_child_para()
		para_new['numPer'] = 1
		x_intrv_btwn_z = zeros_use[1] - zeros_use[0]
		bfs_1per = []
		# extract data 
		while True : 
			# extract slightly more data before and after each intervall
			# to perform a cubicspline approx. - otherwise 
			# the period is too short
			x_strt = zeros_use[ind_zero_strt] - 0.1 * x_intrv_btwn_z
			x_end = zeros_use[ind_zero_strt+2] + 0.1 * x_intrv_btwn_z
			data_prd = self.data[ (self.data[colx] >= x_strt) & (self.data[colx] <= x_end) ]
			x_vals = data_prd[colx].to_list()
			x_vals_spline = np.linspace( x_vals[0], x_vals[-1], num=1000)
			spline_tmp = CubicSpline(data_prd[colx],data_prd[coly])
			data_spln = spline_tmp(x_vals_spline)
			spline_df = pd.DataFrame( { colx : x_vals_spline, coly : data_spln } )
			zero_spln = ana.find_zero_crossings( data = spline_df, colx = colx, coly = coly, strat = { 'name' : 'brute_force', 'lim' : None })
			spline_df = spline_df[ (spline_df[colx] >= zero_spln[0]) & (spline_df[colx] <= zero_spln[2]) ]
			ind_zero_strt = ind_zero_strt + 2
			bf = self.create_field_from_data( data = spline_df, para = para_new )
			bfs_1per.append(bf.get_cntrd_bf(lim_peak = lim_peak, colx = colx, coly = coly ) )
			if (ind_zero_strt + 2) >= len(zeros_use) : 
				break

		[bfs_r, x_min, x_max] = cut_data_support(bfs = bfs_1per, col_cut = 'x')
		x_vals_spline = None
		data_avrgd = None
		num_fields = len(bfs_r)
		for bf in bfs_r : 
			data = bf.data
			x_vals = data[colx].to_list()
			if x_vals_spline is None :
				x_vals_spline = np.linspace( x_vals[0], x_vals[-1], num=1000)
			spline_tmp = CubicSpline(data[colx],data[coly])
			data_spln = spline_tmp(x_vals_spline)
			if data_avrgd is None :
				data_avrgd = data_spln
			else :
				for ind, bval in enumerate(data_spln) : 
					data_avrgd[ind] = data_avrgd[ind] + bval
		for ind, b_val in enumerate(data_avrgd) : 
			data_avrgd[ind] = data_avrgd[ind] / num_fields
		bf_1prd = self.create_field_from_data( data = pd.DataFrame( { colx : x_vals_spline, coly : data_avrgd } ), para = para_new )
		return [ bf_1prd, bfs_r, zeros_use, 2*np_bckwd ]

	def compare_fields(self, bf, colx = 'x', coly = 'By' ) :
		numPer = self.get_para()['numPer']
		if numPer is None :
			numPer = 10
		integr_limit = 100*numPer
		[diff_spl, x_vals] = ana.calc_metric_diff( data_simu = self.data, data_meas = bf.data, integr_limit = integr_limit,colx = colx, coly = coly)
		return [diff_spl, x_vals]

	def extract_amplitude(self, lim_peak = None, colx = 'x', coly = 'By') : 
		[xExs, yExs] = ana.find_extrema(data = self.data, lim = lim_peak, colx = colx, coly = coly)
		maxEx = max(yExs)
		mainExX = []
		mainExY = []
		avrgEx = 0
		for ind, yEx in enumerate(yExs) : 
			if (abs(yEx) >= 0.9*maxEx) and ( abs(yEx) <= 1.1*maxEx ) :
				avrgEx = avrgEx + abs(yEx)
				mainExX.append( xExs[ind] )
				mainExY.append( yEx )
		avrgEx = avrgEx / len(mainExX)
		minEx = avrgEx/20
		sideExX = []
		sideExY = []
		for ind, yEx in enumerate(yExs) : 
			if (abs(yEx) > minEx) and (abs(yEx) < 0.9*maxEx) : 
				sideExX.append( xExs[ind] )
				sideExY.append( yEx )
		return [ avrgEx, mainExX, mainExY, sideExX, sideExY ]

	def integrate_zeros(self, colx = 'x', coly = 'By') : 
		"""
		Check if mirror symmetrie there - center, calc area, compare
		"""
		ars = []
		zeros = self.get_para()['zeros'].copy()
		zeros.append(None)
		for ind, zero in enumerate(zeros) : 
			if ind == 0 : 
				intrvl = [None,zero]
			else :
				intrvl = [zeros[ind-1],zero]
			ar = self.integrate_fld(xs = intrvl, colx = colx, coly = coly)
			ars.append(ar)
		return ars

	def glue_bf(self, bf, pos = 'front', colx = 'x', coly = 'By') : 
		"""
		Glues the field bf to this field adding it either to the front or back (=pos)
		"""
		mydata = self.data.to_dict('records')
		bfdata = bf.data.to_dict('records')
		mydataBeg=mydata[0]['x']
		mydataEnd=mydata[-1]['x']
		bfdataBeg=bfdata[0]['x']
		bfdataEnd=bfdata[-1]['x']
		# if bfdataBeg >= mydataBeg :
		# 	if bfdataBeg <= mydataEnd :

		# 		if bfdataEnd <= mydataEnd :
		# 			return self

		# 		for bfdata_i, bfda in enumerate(bfdata) :
		# 			if bfda['x'] > mydataEnd :
		# 				bfdata=bfdata[bfdata_i:]
		# 				break
		# elif bfdataEnd >= mydataBeg :
		# 	if bfdataEnd <= mydataEnd :

		# 		if bfdataBeg >= mydataBeg :
		# 			return self

		# 		for bfdata_i, bfda in enumerate(bfdata) :
		# 			if bfda['x'] == mydataBeg :
		# 				bfdata=bfdata[:bfdata_i]
		# 				break
		# pdb.set_trace()
		if pos == 'front' : 
			my_pos = mydata[0][colx]
			bf_pos = bfdata[-1][colx]
			my_ind = 0
			bf_ind = -1
			delta = my_pos - bf_pos #- 1
		elif pos == 'back' : 
			my_pos = mydata[-1][colx]
			bf_pos = bfdata[0][colx]
			my_ind = -1
			bf_ind = 0
			delta = my_pos - bf_pos #+ 1
		for line in bfdata : 
			line[colx] = line[colx] + delta
		# average
		avrg = 0.5*(mydata[my_ind][coly] + bfdata[bf_ind][coly] )
		mydata.pop( my_ind )
		bfdata[bf_ind][coly] = avrg
		full_data = mydata+bfdata
		full_data_df = pd.DataFrame( full_data )
		full_data_df = full_data_df.sort_values(by=[colx])
		bf_r_pa = self.make_child_para()
		bf_r = self.create_field_from_data( data = full_data_df, para = bf_r_pa)
		return bf_r

	def reverse_field(self, colx = 'x', coly = 'By') : 
		data = self.data.copy()
		para_new = self.make_child_para()
		data[colx] = data[colx].apply( lambda x : -x )
		sorted_df = data.sort_values(by=[colx], ascending=True).reset_index(drop=True)
		bf_1 = self.create_field_from_data( data = sorted_df, para = para_new )
		return bf_1

	def get_field_interval(self, intrvl, colx = 'x', coly = 'By') :

		data =self.data[ (self.data['x'] >= intrvl[0]) & (self.data['x'] <= intrvl[1]) ]
		field = self.create_field_from_data( data = data, para = None )
		return field

	def cut_bf(self, pos, lim_peak = None, colx = 'x', coly = 'By') : 
		"""
		Cuts the field at pos, returns both parts as separate bfields with interpolated data vals
		"""
		spline = self.get_para()['spline']
		prd_lngth = self.get_para()['prd_lngth']
		xdata = self.data[colx].to_list()

		# num_per_space = (int)(( xdata[-1] - xdata[0] ) / prd_lngth )
		# xlin = np.linspace( xdata[0], xdata[-1], 100*num_per_space )
		# blin = spline(xlin)
		# spline_df = pd.DataFrame( { colx : xlin, coly : blin } )
		# data_half1 = spline_df[ spline_df[colx] <= pos ]
		# data_half2 = spline_df[ spline_df[colx] >= pos ]

		data_half1 =self.data[ (self.data['x'] >= xdata[0]) & (self.data['x'] <= pos) ]
		data_half2 =self.data[ (self.data['x'] >= pos) & (self.data['x'] <= xdata[-1]) ]

		para_new = self.make_child_para()
		bf_1 = self.create_field_from_data( data = data_half1, para = para_new )
		bf_2 = self.create_field_from_data( data = data_half2, para = para_new )
		bf_1.make_ana(lim_peak = lim_peak, colx = colx, coly = coly)
		bf_2.make_ana(lim_peak = lim_peak, colx = colx, coly = coly)
		return [bf_1, bf_2]

	def plot_fld(self, colx = 'x', coly = 'By', title=None, plt_extrm = False, folder = '', add = '',fun=None,nfig=None,plot=True,save=False,filename=None) : 

		if nfig is None :
			fig = plt.figure(figsize=(13*cm_inch, 6.5*cm_inch), dpi=150)
		else :
			fig = plt.figure(nfig,figsize=(13*cm_inch, 6.5*cm_inch), dpi=150)
		data = self.data
		plt.plot(data[colx],data[coly], label = add)

		if plt_extrm :
			para = self.get_para()
			mins = para['minima']
			maxs = para['maxima']
			plt.plot(mins[0], mins[1], "o", linewidth=2, label =  "Minima" )
			plt.plot(maxs[0], maxs[1], "o", linewidth=2, label = "Maxima" )

		if title is None :
			fig.suptitle("B-Field with Extrema", fontsize=14)
		else :
			fig.suptitle(title, fontsize=14)
		plt.axhline(y=0.0,color = 'black')
		# plt.axvline(x=-20.0,color = 'red')
		# plt.axvline(x=20.0,color = 'red')
		plt.xlabel('x [mm]', fontsize=12),
		plt.ylabel(f'{coly} [T]', fontsize=12)
		ax = plt.gca()
		if not (fun is None) : 
			fun(nfig=0)
		ax.legend(loc='upper right')
		if save :
			if filename is None :
				if self.get_para()['file'] is None :
					plt.savefig(folder+"field_extrema_"+add+".png", bbox_inches='tight')
				else:
					plt.savefig(folder+self.get_para()['file']+"_field_extrema_"+add+".png", bbox_inches='tight')
			else :				
				plt.savefig(folder+filename+".png", bbox_inches='tight')
		if plot :
			plt.draw()
		if save :
			plt.clf()
		return fig

	def local_averages(self, intrvl = None, nperiods = None, 
		perlength = None, colx = 'x', coly = 'By', plot_it=True, add = '',nfig=None) : 
		if nperiods is None :
			nperiods = self.get_para()['numPer']
		if perlength is None :
			perlength = self.get_para()['prd_lngth']

		data = self.data
		xs = self.data[colx].to_list()
		if intrvl is None : 
			intrvl = [xs[0],xs[-1]]

		frst_int_spline = self.get_para()['frst_int_spline']
		res = { 'pos' : [], 'avrg' : [] }

		iper = 0
		ind_data = 0
		avrg_avrg = 0
		num_avrgs = 0
		while True : 
			avrg = 0.0
			num_pnts = 0
			istart = intrvl[0] + iper*perlength
			iend = istart + perlength
			if iend > intrvl[-1] :
				break
			res['pos'].append(0.5*(istart+iend))

			intgrl = frst_int_spline(iend) - frst_int_spline(istart)
			avrg = intgrl/perlength

			iper = iper + 1
			res['avrg'].append(avrg)
			avrg_avrg = avrg_avrg + avrg
			num_avrgs = num_avrgs + 1
		avrg_avrg = avrg_avrg / num_avrgs
		if plot_it :
			plt.clf()
			fig = plt.figure(nfig,figsize=(13*cm_inch, 6.5*cm_inch), dpi=150)
			plt.plot(res['pos'],res['avrg'], label = 'Average')
			fig.suptitle("Local Average over x", fontsize=14)
			plt.axhline(y=0.0,color = 'black')
			plt.xlabel('x [mm]', fontsize=12),
			plt.ylabel('Average <F>', fontsize=12)
			plt.savefig("avrg_"+add+".png", bbox_inches='tight')
			plt.draw()

		return res, avrg_avrg

	def add_field(self,field, colx = 'x', middle = None) :
		data = self.data
		data_add = field.data

		my_spline_by=None
		f_spline_by=None
		my_spline_bz=None
		f_spline_bz=None
		if ('By' in self.data.keys()) :
			my_spline_by = CubicSpline(self.data['x'] , self.data['By'])
		if ('By' in field.data.keys()) :
			f_spline_by = CubicSpline(field.data['x'] , field.data['By'])

		if ('Bz' in self.data.keys()) :
			my_spline_bz = CubicSpline(self.data['x'] , self.data['Bz'])
		if ('Bz' in field.data.keys()) :
			f_spline_bz = CubicSpline(field.data['x'] , field.data['Bz'])

		x_me_min = data[colx].to_list()[0]
		x_me_max = data[colx].to_list()[-1]

		x_f_min = data_add[colx].to_list()[0]
		x_f_max = data_add[colx].to_list()[-1]

		minx = x_me_min
		if x_f_min < x_me_min : 
			minx = x_f_min
		maxx = x_me_max
		if x_f_max > x_me_max :
			maxx = x_f_max

		num_pnts = int((len(data) + len(data_add))*1.5)
		unit_x = np.linspace( minx, maxx, num_pnts )
		unit_x_r = []
		unit_x_l = []
		for xval in unit_x :
			if xval < middle : 
				unit_x_l.append(xval)
			else:
				unit_x_r.append(xval)

		new_vals_y=None
		new_vals_z=None
		if ('By' in self.data.keys()) and ('By' in field.data.keys()) :
			new_vals_y_l = [f_spline_by(x) for x in unit_x_l ]
			new_vals_y_r = [my_spline_by(x) for x in unit_x_r ]
			new_vals_y = new_vals_y_l+new_vals_y_r
		elif ('By' in self.data.keys()) : 
			new_vals_y = [my_spline_by(x) for x in unit_x ]
		elif ('By' in field.data.keys()) : 
			new_vals_y = [f_spline_by(x) for x in unit_x ]
		if ('Bz' in self.data.keys()) and ('Bz' in field.data.keys()) :
			new_vals_z_l = [f_spline_bz(x) for x in unit_x_l ]
			new_vals_z_r = [my_spline_bz(x) for x in unit_x_r ]
			new_vals_z = new_vals_z_l+new_vals_z_r
		elif ('Bz' in self.data.keys()) : 
			new_vals_z = [my_spline_bz(x) for x in unit_x ]
		elif ('Bz' in field.data.keys()) : 
			new_vals_z = [f_spline_bz(x) for x in unit_x ]

		if ( not ( new_vals_y is None )) and ( not ( new_vals_z is None )) :
			new_data = pd.DataFrame( { 'x' : unit_x, 'By' : new_vals_y, 'Bz' : new_vals_z } )
		elif ( not ( new_vals_y is None ))  :
			new_data = pd.DataFrame( { 'x' : unit_x, 'By' : new_vals_y } )
		elif ( not ( new_vals_z is None ))  :
			new_data = pd.DataFrame( { 'x' : unit_x, 'Bz' : new_vals_z } )
		pdb.set_trace()
		add_field = self.create_field_from_data( data = new_data )
		return add_field

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
