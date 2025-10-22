"""
Contains functionality for handling files.

Module to handle various file operations including finding, moving, copying, and deleting files.
"""

from unduwave.unduwave_incl import *
from unduwave.quantities.quantities import quantity
from unduwave.attribute_classes.attributes import _attribute
from unduwave.attribute_classes.attributes import _attribute_collection
import unduwave.constants as uc

# compares the results in the lists new_res and old_res 
# and returns absolute relative difference for each list-entry
# if elem in old_res is zero while new_res not: returns -2 
# if both elem zero: returns -1
def compare_sum_res( new_res, old_res  ) :
	res = []
	for ind, elem in enumerate( old_res ) :
		if not ( elem == 0 ) :
			res.append( abs(  ( new_res[ind] - elem ) / elem ) )
		else : 
			if elem == new_res[ind] :
				res.append(-1)
			else : 
				res.append(-2)
	return res

# see Hofmann eq. 8.30, p164
def bessel_sum( type, order_m, arg_a, arg_b, accur = 0.001 ) :
	max_num_runs = 100
	ind = -1
	a = order_m*arg_a
	b = order_m*arg_b
	bessel = ss.jv( 0, a ) * ss.jv( order_m, b )
	while True :
		ind = ind + 1
		j_l_a = ss.jv( ind, a )
		j_ml_a = (-1)**(ind)*j_l_a
		if type == 1 :
			bessel_tmp = j_l_a * ss.jv( order_m + 2*ind, b )
			if not (ind == 0) :
				bessel_tmp = bessel_tmp + j_ml_a * ss.jv( order_m - 2*ind, b )
		elif type == 2 :
			# if ( not ( order_m == 0 ) ) & ( not ( arg_b == 0 ) ) :
			# 	fac = 2*( order_m + 2 * ind )/( order_m * arg_b )
			# 	fac2 = 2*( order_m - 2 * ind )/( order_m * arg_b )
			# 	bessel_tmp = fac*j_l_a * ss.jv( order_m + 2*ind, b )
			# 	bessel_tmp = bessel_tmp + fac2 * j_ml_a * ss.jv( order_m - 2*ind, b )
			# else : 
			bessel_tmp = j_l_a * ( ss.jv( order_m + 2*ind + 1, b ) + ss.jv( order_m + 2*ind - 1, b ) )
			if not (ind==0) :
				bessel_tmp = bessel_tmp + j_ml_a * ( ss.jv( order_m - 2*ind + 1, b ) + ss.jv( order_m - 2*ind - 1, b ) )
		bessel_tmp = bessel + bessel_tmp
		[comp_res] = compare_sum_res( new_res = [ bessel_tmp ], old_res = [bessel] )
		if ind > (order_m+10) :
			if comp_res == -1 :
				break
			elif comp_res >= 0 :
				if comp_res < accur : 
					bessel = bessel_tmp
					break
		bessel = bessel_tmp
		if ind > max_num_runs : 
			break
	# print("Num runs: ", ind)
	return bessel

def xyz_from_spherical(rs,thetas,phis) : 

	shape_tensor = ( len(rs), len(thetas), len(phis), 3 )
	coord_tensor = np.zeros(shape_tensor)

	for i_r, r in enumerate(rs) :
		for i_theta, theta in enumerate(thetas) :
			for i_phi, phi in enumerate(phis) :

				x=r*math.sin(theta)*math.cos(phi)
				y=r*math.sin(theta)*math.sin(phi)
				z=r*math.cos(theta)
				coord_tensor[i_r,i_theta,i_phi,:] = [x,y,z]
	return coord_tensor

def spherical_coord_from_xyz(xs,ys,zs) : 

	shape_tensor = ( len(xs), len(ys), len(zs), 3 )
	coord_tensor = np.zeros(shape_tensor)

	for i_x, x in enumerate(xs) :
		for i_y, y in enumerate(ys) :

			for i_z, z in enumerate(zs) :

				r = math.sqrt( x**2 + y**2 + z**2 )
				if z > 0:
					theta = math.atan(math.sqrt(x**2+y**2)/z) 
				elif z < 0 :
					theta = math.pi + math.atan(math.sqrt(x**2+y**2)/z) 
				elif (z==0.0) and (not (x*y == 0.0)) :
					theta = math.pi/2.0 
				else:
					theta = 0.0

				if x > 0:
					phi = math.atan(y/x)
				elif (x < 0) and (y >= 0) :
					phi = math.atan(y/x) + math.pi
				elif (x < 0) and (y < 0) :
					phi = math.atan(y/x) - math.pi
				elif (x == 0.0) and (y > 0) :
					phi = math.pi/2.0 
				elif (x == 0.0) and (y < 0) :
					phi = -math.pi/2.0
				else: 
					phi=0.0
				coord_tensor[i_x,i_y,i_z,:] = [r,theta,phi]

	return coord_tensor

def make_coords(fun,coords1, coords2) :
	"""
	Takes 2 lists of coordinates for which fun is evaluated. Fun is expected to give back cartesian coordinates of a point.
	Returns 3 2d-arrays containing the grid coordinates
	"""
	coordr = []
	zg = []
	xg = []
	yg = []
	for ind2, coord2 in enumerate(coords2) :
		zg.append([])
		xg.append([])
		yg.append([])
		for ind1, coord1 in enumerate(coords1) :
			x,y,z = fun(coord1,coord2)
			xg[-1].append(x)
			yg[-1].append(y)
			zg[-1].append( z )
	# xg,yg = np.meshgrid(coords1, coords2)
	return [zg,xg,yg]

class aspectrum(_attribute_collection) :

	def __init__(self,ebeam,undulator) :
		self._ebeam=ebeam
		self._undulator=undulator

	# def Fm_sgm_strong(phi,theta,order,undu_K=None,
	# 		undu_K_red=None,beam_gamma_red=None,bessel_arg_a=None,bessel_arg_b=None,
	# 		bessel_sum_1=None,bessel_sum_2=None,Fm_denom=None,Fm_fac=None) : 
	# 	if Fm_fac is None:
	# 		Fm_fac = 3*order**2/(math.pi*(1+undu_K**2/2.0)**2*undu_K_red**2)
	# 	if bessel_arg_a is None :
	# 		bessel_arg_a = undu_K_red**2/(4*(1+beam_gamma_red**2*theta**2))
	# 	if bessel_arg_b is None :
	# 		bessel_arg_b = 2*undu_K_red*beam_gamma_red*theta*math.cos(phi)/(1+beam_gamma_red**2*theta**2)
	# 	if bessel_sum_1 is None :
	# 		bessel_sum_1 = bessel_sum( type=1, order_m=order, arg_a=bessel_arg_a, arg_b=bessel_arg_b, accur = 0.001 )
	# 	if bessel_sum_2 is None :
	# 		bessel_sum_2 = bessel_sum( type=2, order_m=order, arg_a=bessel_arg_a, arg_b=bessel_arg_b, accur = 0.001 )
	# 	if Fm_denom is None: 
	# 		Fm_denom = (1+beam_gamma_red**2*theta**2)**3
	# 	return Fm_fac*(2*bessel_sum_1*beam_gamma_red*theta*math.cos(phi)-bessel_sum_2*undu_K_red)**2/Fm_denom

	def spectral_fun(self,omega,order) : 
		fund_freq=self._undulator.properFrequencyStrong.get()
		N_periods=self._undulator.numPeriods.get()
		freq_arg = (omega-order*fund_freq)/fund_freq*math.pi*N_periods
		if freq_arg == 0.0 :
			freq_fun = N_periods/fund_freq
		else :
			freq_fun = N_periods/fund_freq*( math.sin(freq_arg)/freq_arg )**2
		return freq_fun

	def fun_transv(self,theta,phi) : 
		beam_energy=self._ebeam.beam_en.get()*1e9*uc.q_el
		gamma = beam_energy/(uc.m_el*uc.v_c**2)
		beta = math.sqrt( 1 - 1 / gamma )

		r = 1/gamma**4 * ((1-beta*math.cos(theta))**2-(1-beta**2)*math.sin(theta)**2*math.cos(phi)**2)/(1-beta*math.cos(theta))**5
		coords = xyz_from_spherical([r],[theta],[phi])
		return [coords[0,0,0,0],coords[0,0,0,1],coords[0,0,0,2]]

	def fun_transv_ultra_rel(self,theta,phi) : 
		beam_energy=self._ebeam.beam_en.get()*1e9*uc.q_el
		gamma = beam_energy/(uc.m_el*uc.v_c**2)
		beta = math.sqrt( 1 - 1 / gamma )

		r = gamma**2*(1-2*gamma**2*theta**2*math.cos(2*phi)+gamma**4*theta**4)/(1+gamma**2*theta**2)**5
		coords = xyz_from_spherical([r],[theta],[phi])
		return [coords[0,0,0,0],coords[0,0,0,1],coords[0,0,0,2]]

	def flux_strong_undu(self,x,y,z,omega,order) : 

		current=self._ebeam.current.get()
		beam_energy=self._ebeam.beam_en.get()*1e9*uc.q_el

		undu_lambda=self._undulator.periodLength.get()

		undu_B_eff=self._undulator.bEff.get()
		N_periods=self._undulator.numPeriods.get()

		beam_gamma = beam_energy/(uc.m_el*uc.v_c**2)
		beam_beta = math.sqrt( 1 - 1 / beam_gamma )
		v_el = uc.v_c * beam_beta

		undu_k = 2 * math.pi / ( undu_lambda )
		undu_K = uc.q_el * undu_B_eff / ( uc.m_el * uc.v_c * undu_k )
		# undu_K = math.sqrt(2)
		undu_Omega = undu_k*beam_beta*uc.v_c

		beam_beta_red = beam_beta * ( 1 - undu_K**2/(4*beam_beta**2 * beam_gamma**2) )
		beam_gamma_red=beam_gamma/math.sqrt(1+undu_K**2/2)
		undu_K_red = undu_K / math.sqrt( 1 +  undu_K**2/2 )

		avrg_total_power = uc.R_el*uc.v_c**3*uc.m_el*undu_k**2*undu_K**2*beam_gamma**2/3.0
		const_int_fac = avrg_total_power/uc.hbar*beam_gamma_red**2
		def flux_strong_undu_intgrnd(Theta,phi) : 
			bessel_arg_a = undu_K_red**2/(4*(1+beam_gamma_red**2*Theta**2))
			bessel_arg_b = 2*undu_K_red*beam_gamma_red*Theta*math.cos(phi)/(1+beam_gamma_red**2*Theta**2)
			bessel_sum_1 = bessel_sum( type=1, order_m=order, arg_a=bessel_arg_a, arg_b=bessel_arg_b, accur = 0.001 )
			bessel_sum_2 = bessel_sum( type=2, order_m=order, arg_a=bessel_arg_a, arg_b=bessel_arg_b, accur = 0.001 )
			Fm_fac = 3*order**2/(math.pi*(1+undu_K**2/2.0)**2*undu_K_red**2)
			Fm_denom = (1+beam_gamma_red**2*Theta**2)**3
			Fm_sgm = Fm_fac*(2*bessel_sum_1*beam_gamma_red*Theta*math.cos(phi)-bessel_sum_2*undu_K_red)**2/Fm_denom
			Fm_pi=3*Fm_fac*(2*bessel_sum_1*beam_gamma_red*Theta*math.sin(phi))**2/Fm_denom

			fund_freq = 2*beam_gamma_red**2*undu_Omega/(1+beam_gamma_red**2*Theta**2)
			freq_fun = self.spectral_fun(omega=omega,order=order)

			integrand = (Fm_sgm + Fm_pi)*freq_fun
			# integrand = (Fm_sgm + Fm_pi)
			# integrand = (Fm_sgm)

			# integrand=math.exp(-Theta)
			return integrand

		# print(f'Fund freq {2*beam_gamma_red**2*undu_Omega*uc.hbar/uc.q_el}')
		coord_tensor = spherical_coord_from_xyz([x],[y],[z])
		# pdb.set_trace()
		res = flux_strong_undu_intgrnd(coord_tensor[0,0,0,1],coord_tensor[0,0,0,2])
		return res*const_int_fac

		on_ax_fac=uc.fein_const*order**2*beam_gamma**2*current*undu_K**2*N_periods**2
		on_ax_fac = on_ax_fac/(q_el*(1+undu_K**2/2.0)**2)
		j1=ss.jv( int((order -1)/2.0), order*undu_K_red**2/4.0 )
		j2=ss.jv( int((order +1)/2.0), order*undu_K_red**2/4.0 )
		on_ax_fac=on_ax_fac*(j1-j2)**2
		on_ax_hof_orig=on_ax_fac

		on_ax_kim=uc.fein_const*N_periods**2*beam_gamma**2*current/uc.q_el
		kim_arg=order*undu_K**2/(4*(1+undu_K**2/2.0))
		j1_kim=ss.jv( int((order -1)/2.0), kim_arg )
		j2_kim=ss.jv( int((order +1)/2.0), kim_arg )
		fn_kim=undu_K**2*order**2/(1+undu_K**2/2.0)**2*(j1_kim-j2_kim)**2
		on_ax_kim_orig=on_ax_kim*fn_kim

		fac_fun = 1.744*1e14
		fac_fun_real = uc.fein_const/(m_el*v_c**2/(uc.q_el*1e9))**2*1/uc.q_el
		on_ax_kim_flux_orig_f = (N_periods**2)*(beam_energy/(uc.q_el*1e9))**2*current*fn_kim
		on_ax_kim_flux_orig = on_ax_kim_flux_orig_f * fac_fun/36
		on_ax_kim_flux_real = on_ax_kim_flux_orig_f * fac_fun_real

		print(f"Original Flux-Densities: Hof: {on_ax_hof_orig:.4e} Nphot/(s*mrad^2) and kim: {on_ax_kim_orig:.4e}")
		print(f"Flux-Density-Kim natural units with Kim-Fac: {on_ax_kim_flux_orig:.4e} and with RealFactor: {on_ax_kim_flux_real:.4e}")
		print(f"Flux-Factor natural units Comparisson: Kim: {fac_fun:.4e} and ana: {fac_fun_real:.4e}")

		unit_fac_mrad=1e6/(z*1e3)**2
		unit_fac_rad=1/(z*1e3)**2
		print(f"Unit-Change-mrad: {on_ax_hof_orig*unit_fac_mrad:.4e} Npho/(s mm^2) and Unit-Change-rad: {on_ax_hof_orig*unit_fac_rad:.4e} Npho/(s mm^2)")

		fac_flux = 1.431e14
		flux_kim = fac_flux * N_periods * (1+undu_K**2/2.0)*fn_kim/order * current

		print(f"The kim-flux then is: {flux_kim:.4e}")

		pdb.set_trace()

		return res*const_int_fac

		num_integration_pnts = 100

		r = 6 # m
		dX_scrn = 0.04 # m
		dY_scrn = 0.04 # m

		xs = np.linspace(-dX_scrn/2.0,dX_scrn/2.0,num_integration_pnts)
		ys = np.linspace(-dY_scrn/2.0,dY_scrn/2.0,num_integration_pnts)
		zs = [r]
		coord_tensor = spherical_coord_from_xyz(xs,ys,zs)

		x_plot = []
		y_plot = []
		distro_plot = []
		shape_tensor = ( len(xs), len(ys) )
		distro_flux = np.zeros(shape_tensor)

		for ind_x, x in enumerate(xs) : 
			for ind_y, y in enumerate(ys) : 
				theta = coord_tensor[ind_x,ind_y,0,1]
				phi = coord_tensor[ind_x,ind_y,0,2]
				x_plot.append(x)
				y_plot.append(y)
				# distro_plot.append(flux_strong_undu_intgrnd(theta,phi))
				distro_flux[ind_x,ind_y] = flux_strong_undu_intgrnd(theta,phi)

		# pdb.set_trace()

		# shape_tensor = ( len(theta_vals), len(phi_vals) )
		# distro_flux = np.zeros(shape_tensor)

		# theta_plot = []
		# phi_plot = []
		# distro_plot = []
		# for ind_theta, theta in enumerate(theta_vals) : 
		# 	for ind_phi, phi in enumerate(phi_vals) : 
		# 		theta_plot.append(theta)
		# 		phi_plot.append(phi)
		# 		distro_plot.append(flux_strong_undu_intgrnd(theta,phi))
		# 		# distro_flux[ind_theta,ind_phi] = flux_strong_undu_intgrnd(theta,phi)

		integral = 0.0
		# integral = const_int_fac * scipy.integrate.simps(scipy.integrate.simps(distro_flux, phi_vals), theta_vals)
		# print(f"Integrated for freq {omega/q_el*hbar}")

		if plot : 
			fig = plt.figure(figsize=(13*cm_inch, 6.5*cm_inch), dpi=150)
			fig.suptitle("Flux Distribution", fontsize=14)
			ax = plt.gca()

			# ax.contour(theta_vals, phi_vals, distro_flux)
			ax = fig.add_subplot(111, projection='3d')

			X, Y = np.meshgrid(xs,ys)
			surf = ax.plot_surface(X, Y, distro_flux, cmap=cm.coolwarm,
			                       linewidth=0, antialiased=False)

			# ax.plot_trisurf(x_plot, y_plot, distro_plot, cmap=cm.jet, linewidth=0.2)
			# ax.plot_trisurf(theta_plot, phi_plot, distro_plot, cmap=cm.jet, linewidth=0.2)
			ax.set_xlabel('x', fontsize=12)
			ax.set_ylabel('y', fontsize=12)
			# ax.set_zlabel('Flux', fontsize=12)
			plt.savefig(f"flux_distro_{omega:.2e}_order_{order}.png"   , bbox_inches='tight')
			plt.show()

		return integral