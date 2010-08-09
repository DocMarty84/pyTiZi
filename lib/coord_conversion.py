def Cartesian_To_Fractional(x, y, z, box):
	z = (z/box.temp_gamma_term) / box.c
	y = ((y-z*box.c*box.temp_beta_term)/box.temp_gamma_sin) / box.b
	x = (x-y*box.b*box.temp_gamma_cos-z*box.c*box.temp_beta_cos) / box.a
	
	return [x,y,z]

def Fractional_To_Cartesian(x, y, z, box):
	x = x*box.a + y*box.b*box.temp_gamma_cos + z*box.c*box.temp_beta_cos
	y = y*box.b*box.temp_gamma_sin + z*box.c*box.temp_beta_term
	z = z*box.c*box.temp_gamma_term
	
	return [x,y,z]
