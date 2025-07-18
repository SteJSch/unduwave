"""
Contains functionality for handling files.

Module to handle various file operations including finding, moving, copying, and deleting files.
"""

from unduwave.unduwave_incl import *

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
	ind = 0
	a = order_m*arg_a
	b = order_m*arg_b
	bessel = ss.jv( 0, a ) * ss.jv( order_m, b )
	while True :
		ind = ind + 1
		j_l_a = ss.jv( ind, a )
		j_ml_a = (-1)**(ind)*j_l_a
		if type == 1 :		
			bessel_tmp = j_l_a * ss.jv( order_m + 2*ind, b )
			bessel_tmp = bessel_tmp + j_ml_a * ss.jv( order_m - 2*ind, b )
		elif type == 2 :
			if ( not ( order_m == 0 ) ) & ( not ( arg_b == 0 ) ) :
				fac = 2*( order_m + 2 * ind )/( order_m * arg_b )
				fac2 = 2*( order_m - 2 * ind )/( order_m * arg_b )
				bessel_tmp = fac*j_l_a * ss.jv( order_m + 2*ind, b )
				bessel_tmp = bessel_tmp + fac2 * j_ml_a * ss.jv( order_m - 2*ind, b )
			else : 
				bessel_tmp = j_l_a * ( ss.jv( order_m + 2*ind + 1, b ) + ss.jv( order_m + 2*ind - 1, b ) )
				bessel_tmp = bessel_tmp + j_ml_a * ( ss.jv( order_m - 2*ind + 1, b ) + ss.jv( order_m - 2*ind - 1, b ) )
		bessel_tmp = bessel + bessel_tmp
		[comp_res] = compare_sum_res( new_res = [ bessel_tmp ], old_res = [bessel] )
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

def electr_field_full( times, obs_pnt, undu_paras, beam_paras ) :
	res_e = []
	x_c = undu_paras.max_deflection
	z_c = undu_paras.unduPara**2 / (8*undu_paras.beta**2 *undu_paras.gamma**2 * undu_paras.waveNum)
	z_cL = undu_paras.red_beta * velC
	for time in times :
		r_el = [ x_c * math.cos( undu_paras.unduFreq * time ), 0, \
				z_cL * time + z_c * math.sin( 2 * undu_paras.unduFreq * time ) ]
		res_e.append(r_el)

	return res_e

def flux_at_order_m( harmOrder, undu_paras, beam_paras ) : 
	fac = harmOrder * pc.feinC * pc.velC * beam_paras.gamma**2 * undu_paras.waveNum * undu_paras.unduPara**2 * \
		undu_paras.nPeriods / (2*math.pi*(1+undu_paras.unduPara**2/2)**2)
	bessel1 = ss.jv( (harmOrder-1)/2, harmOrder*undu_paras.red_unduPara**2/4 )
	bessel2 = ss.jv( (harmOrder+1)/2, harmOrder*undu_paras.red_unduPara**2/4 )
	maxFlux = fac*( bessel1 - bessel2 )**2 * ( 1/(4*math.pi)*0.1 )
	return maxFlux

# see Hofmann eq. 8.32, p165
def electric_field_strong_time_fourier(time, arg_a, arg_b, order, theta, phi, undu_paras, beam_paras ) : 
	res = - pc.elCharge * undu_paras.waveNum * undu_paras.red_gamma**3 / ( math.pi * pc.vacPermittivity * pc.elRadius * ( 1 + undu_paras.red_gamma**2*theta**2 ) )
	res = res * math.sin( order * ( undu_paras.get_first_harm( theta = theta ) * time + math.pi / 2 ) )
	bessel_sum1 = bessel_sum( type = 1, order_m = order, arg_a = arg_a, arg_b = arg_b, accur = 0.001 )
	bessel_sum2 = bessel_sum( type = 2, order_m = order, arg_a = arg_a, arg_b = arg_b, accur = 0.001 )
	x_comp = res * ( 2 * undu_paras.red_gamma * theta * math.cos( phi ) * bessel_sum1 - undu_paras.red_unduPara * bessel_sum2 )
	y_comp = res * ( 2 * undu_paras.red_gamma*theta*math.sin(phi)*bessel_sum1 )
	return [ x_comp, y_comp, 0 ]

# see Hofmann eq. 8.34, p165
def electric_field_strongU(time, theta, phi, undu_paras, beam_paras, accur = 0.001 ) : 
	e_fields = [ [0, 0, 0] ]
	num_avrg = 4
	arg_a = undu_paras.red_unduPara**2 / ( 4 * ( 1 + undu_paras.red_gamma**2 * theta**2 ) )
	arg_b = 2*undu_paras.red_unduPara*undu_paras.red_gamma*theta*math.cos(phi) / ( 1 + undu_paras.red_gamma**2 * theta**2 )	
	max_num_runs = 100
	ind = 1
	while True :		
		new_field = electric_field_strong_time_fourier(time = time, \
				arg_a = arg_a, arg_b = arg_b, order = ind, theta = theta, phi = phi, undu_paras = undu_paras, beam_paras = beam_paras )
		zeros = []
		for ind_f, elem in enumerate( new_field ) :
			if not ( elem == 0 ) :
				new_field[ind_f] = elem + e_fields[-1][ind_f]
			else : 
				zeros.append(1)
		if len(zeros) < 3 : 
			comp_res = compare_sum_res( new_res = new_field, old_res= e_fields[-1] )
			e_fields.append( new_field )
			stop = []
			for elem in comp_res : 
				if elem == -1 : 
					stop.append(1)
				elif elem >= 0:
					if elem < accur :
						stop.append(1)
			if len(stop) == 3 :
				break
		ind = ind + 1
		if ind > max_num_runs : 
			break
	return e_fields[-1]
