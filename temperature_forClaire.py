from numpy import loadtxt, pi, array, sqrt, meshgrid, linspace, zeros, atleast_2d
from matplotlib import pyplot as plt
from sys import exit

def dipole_temperatures(cs_single_list, Evec_file, np_radius, inc_irrad, kappa, save_T_data=False, save_Q_data=False):
	'''
	Input Arguments:
	- cs_single_list: array of SLRpy_finite-style cross section output file names for individual NPs. 2nd column is absorption cross section in [nm^2].
	- Evec_file: File containing total local E-field at each dipole location. 
	- np_radius: array-type object of NP radii [nm].
	- inc_irrad: incident irradiance [W * nm^-2].
	- kappa: Thermal conductivity of embedding medium [W * nm^-1 * K^-1]
	- save_T_data: switch to enable writing temp. data to file -- see NOTE_2 below.
	- save_Q_data: switch to enable writing heat power absorbed data to file -- see NOTE_3 below. 

	Returns:
	- dT_array: array of internal temperatures -- see NOTE_2 below.
	- pos_array: [N x 3] array of cartesian vector components for all dipoles
	- dipole_rad_array: [N x] array of sphere radii -- ordering matches that of pos_array

	NOTE_1: 
		- (CAUTION): At the moment care must be taken to order np_radius and cs_single_list elements identically! 
	
	NOTE_2:
		- Whether it is the TOTAL or just the direct OPTICAL absorption induced temp. increase
			depends on if save_T_data is enabled or disabled, respectively. Furthermore, if save_T_data is 
			enabled, then a file named "T_internal.out" is produced with N lines for N dipoles and the format:

			col_1: 		col_2:		col_3: 			col_4: 			col_5:
			x_n [nm]	y_n [nm] 	z_n [n]		T_(total) [K] 	T_(optical) [K]

	NOTE_3:
		- If save_Q_data is enabled, then the heat power absorbed data is written to a file named 
			"Q_heat_power.out". This file will have N lines for N dipoles and the format:

			col_1: 		col_2:		col_3: 		col_4: 
			x_n [nm]	y_n [nm] 	z_n [n]		Q_n [W]
	'''

	# read-in absorption cross section data from specifiedfiles:
	abs_cs_list = []
	for cs_file in cs_single_list:
	 	cs_single_data = loadtxt(cs_file)
	 	abs_cs_list.append(cs_single_data[1])

	# All field data contained in single file: 
	Evec_data = atleast_2d(loadtxt(Evec_file))

	num_dipoles = Evec_data.shape[0]
	#print("\nNumber of dipoles = {0:d}".format(num_dipoles))

	dT_array  = []
	pos_array = []
	dipole_rad_array = []
	Q_array = []
	# loop over target dipoles:
	for i in range(num_dipoles):

		type_id = int(Evec_data[i,10])

		# log the dipole position and radius [nm]:
		pos_array.append(array([Evec_data[i,0], Evec_data[i,1], Evec_data[i,2]]))
		dipole_rad_array.append(np_radius[type_id])

		# evaluate local electromagnetic enhancement factor:
		EF = Evec_data[i,9]

		# heat power absorbed [W]:
		Q = abs_cs_list[type_id]*inc_irrad*(EF**2)
		Q_array.append(Q)
		print("Q_{0:d} = {1:e} W ".format(type_id, Q))

		# temperature change due to optical absorption [K].
		dT = Q/(4*pi*kappa*np_radius[type_id])

		dT_array.append(dT)
	#for

	if save_Q_data:
		out_file_name = "Q_heat_power.out"
		outfile = open(out_file_name, "w")
		for i in range(num_dipoles):
			output_line = "{0:.2f}\t{1:.2f}\t{2:.2f}\t{3:e}\n".format(Evec_data[i,0], Evec_data[i,1], Evec_data[i,2], Q_array[i])
			outfile.write(output_line)
		outfile.close()
		print("\nHeat power (Q) logged to file: {0:s}".format(out_file_name))
		print("Units = [Watts]\n")
	#if

	if save_T_data:
		'''
		Up to this point, the internal temperature increase calculated is that
		due only to the direct optical absorption. If this block is enabled, 
		then the TOTAL internal field is evaluated and then written to the 
		file: T_internal.out. 
		'''

		T_int_total = zeros(num_dipoles)

		for i in range(num_dipoles):
			for j in range(num_dipoles):
				if i == j:
					T_int_total[i] += dT_array[i]
				else:
					R = sqrt( (pos_array[i][0] - pos_array[j][0])**2 + (pos_array[i][1] - pos_array[j][1])**2 + (pos_array[i][2] - pos_array[j][2])**2 )
					T_int_total[i]+= dT_array[j]*np_radius[int(Evec_data[j,10])]/R


		out_file_name = "T_internal.out"
		outfile = open(out_file_name, "w")
		outfile.write("## Format: col0=x[nm]\tcol1=y[nm]\tcol2=z[nm]\tcol3=T_total[K]\tcol4=T_optical[K]\n")
		for i in range(num_dipoles):
			output_line = "{0:.2f}\t{1:.2f}\t{2:.2f}\t{3:.2f}\t{4:.2f}\n".format(Evec_data[i,0], Evec_data[i,1], Evec_data[i,2], T_int_total[i], dT_array[i])
			outfile.write(output_line)
		outfile.close()
	#if

	return array(dT_array), array(pos_array), array(dipole_rad_array)



if __name__ == "__main__":


	## DMREF PROPOSAL INPUTS:
	
	# Incident laser irradiance[W * nm^-2]
	I =  10**(-10) # this corresponds to 10^8 W/m^2

	# Thermal conductivity of embedding medium [W * nm^-1 * K^-1]   
	kappa = 0.6*(10**-9) 	# currently water

	'''
	## Single 80 nmD:
	dir_base = "/Users/marcbourgeois/Desktop/Work/UWashington/THERMOPLAS/SLRpy_FINITE_HYAK_Copy/single_80nmD_n1.0/"
	cs_file_list = [dir_base + "CS_379.5nm.out"]
	Evec_file = dir_base + "Evec_379.50.out"
	np_radius_list = [40] # [nm]
	'''

	'''
	## Single 90 nmD:
	dir_base = "/Users/marcbourgeois/Desktop/Work/UWashington/THERMOPLAS/SLRpy_FINITE_HYAK_Copy/single_90nmD_n1.0/"
	cs_file_list = [dir_base + "CS_379.5nm.out"]
	Evec_file = dir_base + "Evec_379.50.out"
	np_radius_list = [45] # [nm]
	'''

	
	## Single 100 nmD:
	#dir_base = "/Users/marcbourgeois/Desktop/Work/UWashington/THERMOPLAS/SLRpy_FINITE_HYAK_Copy/single_100nmD_n1.0/"
	#cs_file_list = [dir_base + "CS_379.5nm.out"]
	#Evec_file = dir_base + "Evec_379.50.out"
	#np_radius_list = [50] # [nm]
	

	'''
	## Dimer HOMO:
	dir_base = "/Users/marcbourgeois/Desktop/Work/UWashington/THERMOPLAS/SLRpy_FINITE_HYAK_Copy/"
	cs_file_list = [dir_base + "single_90nmD_n1.0/CS_379.5nm.out", dir_base + "single_90nmD_n1.0/CS_379.5nm.out"]
	Evec_file = dir_base + "DIMER_HOMO/Evec_379.50.out"
	#cs_file_list = [dir_base + "single_90nmD_n1.0/CS_416nm.out", dir_base + "single_90nmD_n1.0/CS_416nm.out"]
	#Evec_file = dir_base + "DIMER_HOMO/Evec_416.00.out"
	#cs_file_list = [dir_base + "single_90nmD_n1.0/CS_431nm.out", dir_base + "single_90nmD_n1.0/CS_431nm.out"]
	#Evec_file = dir_base + "DIMER_HOMO/Evec_431.00.out"
	np_radius_list = [45, 45] # [nm]
	'''

	'''
	## Dimer HETERO:
	dir_base = "/Users/marcbourgeois/Desktop/Work/UWashington/THERMOPLAS/SLRpy_FINITE_HYAK_Copy/"
	cs_file_list = [dir_base + "single_100nmD_n1.0/CS_379.5nm.out", dir_base + "single_80nmD_n1.0/CS_379.5nm.out"]
	Evec_file = dir_base + "DIMER_HETERO/Evec_379.50.out"
	#cs_file_list = [dir_base + "single_100nmD_n1.0/CS_416nm.out", dir_base + "single_80nmD_n1.0/CS_416nm.out"]
	#Evec_file = dir_base + "DIMER_HETERO/Evec_416.00.out"
	#cs_file_list = [dir_base + "single_100nmD_n1.0/CS_431nm.out", dir_base + "single_80nmD_n1.0/CS_431nm.out"]
	#Evec_file = dir_base + "DIMER_HETERO/Evec_431.00.out"
	np_radius_list = [50, 40] # [nm]
	'''

	'''
	## Dimer HOMO 1D:
	dir_base = "/Users/marcbourgeois/Desktop/Work/UWashington/THERMOPLAS/SLRpy_FINITE_HYAK_Copy/"
	np_radius_list = [45, 45] # [nm]
	#cs_file_list = [dir_base + "single_90nmD_n1.0/CS_379.5nm.out", dir_base + "single_90nmD_n1.0/CS_379.5nm.out"]
	#Evec_file = dir_base + "1D_DIMER_HOMO/Evec_379.50.out"
	#cs_file_list = [dir_base + "single_90nmD_n1.0/CS_416nm.out", dir_base + "single_90nmD_n1.0/CS_416nm.out"]
	#vec_file = dir_base + "1D_DIMER_HOMO/Evec_416.00.out"
	cs_file_list = [dir_base + "single_90nmD_n1.0/CS_431nm.out", dir_base + "single_90nmD_n1.0/CS_431nm.out"]
	Evec_file = dir_base + "1D_DIMER_HOMO/Evec_431.00.out"
	'''


	## Dimer HETERO 1D:
	#dir_base = "/Users/marcbourgeois/Desktop/Work/UWashington/THERMOPLAS/SLRpy_FINITE_HYAK_Copy/"
	#np_radius_list = [50, 40] # [nm]
	#cs_file_list = [dir_base + "single_100nmD_n1.0/CS_379.5nm.out", dir_base + "single_80nmD_n1.0/CS_379.5nm.out"]
	#Evec_file = dir_base + "1D_DIMER_HETERO/Evec_379.50.out"
	#cs_file_list = [dir_base + "single_100nmD_n1.0/CS_416nm.out", dir_base + "single_80nmD_n1.0/CS_416nm.out"]
	#Evec_file = dir_base + "1D_DIMER_HETERO/Evec_416.00.out"
	#cs_file_list = [dir_base + "single_100nmD_n1.0/CS_431nm.out", dir_base + "single_80nmD_n1.0/CS_431nm.out"]
	#Evec_file = dir_base + "1D_DIMER_HETERO/Evec_431.00.out"


	## 1D MONOMER -- Windowed Source:
	dir_base = "/Users/marcbourgeois/Desktop/Work/UWashington/THERMOPLAS/SLRpy_FINITE_HYAK_Copy/"
	np_radius_list = [45] # [nm]
	cs_file_list = [dir_base + "single_90nmD_n1.0/CS_421nm.out"]
	#Evec_file = dir_base + "1D_MONO/N_201/WINDOWED/Evec_379.50_full_center.out"
	#Evec_file = dir_base + "1D_MONO/N_201/WINDOWED/Evec_379.50_m40um_0um_LHS.out"
	#cs_file_list = [dir_base + "single_90nmD_n1.0/CS_421nm.out"]
	Evec_file = dir_base + "1D_MONO/N_201/WINDOWED/Evec_421.00_4000nm_center.out"
	#Evec_file = dir_base + "1D_MONO/N_201/WINDOWED/Evec_421.00_10000nm_center.out"
	#Evec_file = dir_base + "1D_MONO/N_201/WINDOWED/Evec_421.00_20000nm_center.out"
	#Evec_file = dir_base + "1D_MONO/N_201/WINDOWED/Evec_421.00_40000nm_center.out"
	#Evec_file = dir_base + "1D_MONO/N_201/WINDOWED/Evec_421.00_full_center.out"
	#Evec_file = dir_base + "1D_MONO/N_201/WINDOWED/Evec_421.00_m40um_m20um_LHS.out"
	#Evec_file = dir_base + "1D_MONO/N_201/WINDOWED/Evec_421.00_m40um_0um_LHS.out"
	#Evec_file = dir_base + "1D_MONO/N_201/WINDOWED/Evec_421.00_m40um_20um_LHS.out"
	#cs_file_list = [dir_base +"single_90nmD_n1.0/CS_434nm.out"]
	#Evec_file = dir_base + "1D_MONO/N_201/WINDOWED_2.5AOI/Evec_434.00_m2.5AOI.out"


	T_array, r_array, a_array  = dipole_temperatures(cs_file_list, Evec_file, np_radius_list, I, kappa, save_T_data=True, save_Q_data=False)

	
