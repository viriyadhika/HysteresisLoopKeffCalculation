import numpy as np
import matplotlib.pyplot as plt
from dataclasses import dataclass
from scipy.integrate import simps
import tkinter
from tkinter import filedialog

"""
This script calculates the Keff from the field vs. magnetization curve
The input format is a .dat file exported from Origin Pro 
Format is as following:

Row
1. 1st row Header (to be thrown out)
2. 2nd row onwards (data)

Column: according to the follwoing variables in input options
1. field_column_num_op
2. mag_column_num_op
3. field_column_num_ip
4. mag_column_num_ip

External field is in Oe and Magnetization in emu/cc
"""

@dataclass
class InputOptions:
	#data_path: str = r''
	#results_path: str = r''
	forc_path: str = None
	labels: str = None
	experimental_results_path_filename: str = ''
	title: str = ''
	input_file_ext: str = '.dat'
	legend: bool = True
	field_column_num_op: int = 0 #column number of external Field_op
	mag_column_num_op: int = 1 #column number of op data
	field_column_num_ip: int = 2 #column number of external Field_ip
	mag_column_num_ip: int = 3 #column number of ip data
	linspace_data_point: int = 1000 #no. of data point in the linspace
	extrapolate_point: int = 10 #how many data points to extrapolate

def normalize(base_data, normalize_data):
	average_normalized = np.mean(normalize_data[0:10])
	average_base = np.mean(base_data[0:10])
	
	return average_normalized / average_base * base_data

def test():
	in_opt = InputOptions()
	tk_root = tkinter.Tk()
	files = filedialog.askopenfilename(title = 'Select .dat files')
	readfile = np.genfromtxt(files, delimiter='\t', comments='#',skip_header=1, invalid_raise=False)# put the file address here
	tk_root.destroy()
  
	i = 0
	while i < np.size(readfile,1)/4:
		analysis_main(readfile, 4*i + in_opt.field_column_num_ip, 4*i + 
			in_opt.field_column_num_op, 4*i + in_opt.mag_column_num_ip, 4*i + in_opt.mag_column_num_op)
		i += 1
	
def analysis_main(file, field_column_num_ip1, field_column_num_op1, mag_column_num_ip1, mag_column_num_op1):
	print('\nTitle: Data no.', int(field_column_num_ip1/4) + 1)
	
	in_opt = InputOptions()
	
#Filter the data from NaN
	mag_data_op = file[:,mag_column_num_op1][~np.isnan(file[:, mag_column_num_op1])]
	mag_data_ip= file[:, mag_column_num_ip1][~np.isnan(file[:, mag_column_num_ip1])]
	field_data_op = file[:, field_column_num_op1][~np.isnan(file[:, field_column_num_op1])]
	field_data_ip = file[:, field_column_num_ip1][~np.isnan(file[:, field_column_num_ip1])]	

#Normalize the data
	mag_data_ip = normalize(mag_data_ip, mag_data_op)

#Dividing plots in to sections
#Section one is a negative sweep, section two is to return to positive H
	i = 1
	field_data_op_diff = np.sign(np.diff(field_data_op))
	while i < field_data_op.size:
		#The point where H flips direction
		if field_data_op_diff[i] != field_data_op_diff[i-1]: 
			break
		i += 1
	op_flipping = i
	
	i = 1
	field_data_ip_diff = np.sign(np.diff(field_data_ip))
	while i < field_data_ip.size:
		#The point where H flips direction
		if field_data_ip_diff[i] != field_data_ip_diff[i-1]: 
			break
		i += 1
	ip_flipping = i

	mag_data_op1 = np.absolute(mag_data_op[:op_flipping])
	mag_data_op2 = np.absolute(mag_data_op[op_flipping:])
	mag_data_ip1 = np.absolute(mag_data_ip[:ip_flipping])
	mag_data_ip2 = np.absolute(mag_data_ip[ip_flipping:])
	
	field_data_op1 = field_data_op[:op_flipping]
	field_data_op2 = field_data_op[op_flipping:]
	field_data_ip1 = field_data_ip[:ip_flipping]
	field_data_ip2 = field_data_ip[ip_flipping:]

#Sort the mag data based on the field data 
	mag_data_op1 = np.array(mag_data_op1)[np.array(sorted(range(len(field_data_op1)), key = lambda k:field_data_op1[k]))]
	mag_data_ip1 = np.array(mag_data_ip1)[np.array(sorted(range(len(field_data_ip1)), key = lambda k:field_data_ip1[k]))]
	mag_data_op2 = np.array(mag_data_op2)[np.array(sorted(range(len(field_data_op2)), key = lambda k:field_data_op2[k]))]
	mag_data_ip2 = np.array(mag_data_ip2)[np.array(sorted(range(len(field_data_ip2)), key = lambda k:field_data_ip2[k]))]

#Sort the field data
	field_data_op1 = sorted(field_data_op1)
	field_data_ip1 = sorted(field_data_ip1)
	field_data_op2 = sorted(field_data_op2)
	field_data_ip2 = sorted(field_data_ip2)

#To detect which is the highest field
	min_field1_low = min(field_data_op1[0], field_data_ip1[0])
	min_field1_high = max(field_data_op1[-1], field_data_ip1[-1])
	
	min_field2_low = min(field_data_op2[0], field_data_ip2[0])
	min_field2_high = max(field_data_op2[-1], field_data_ip2[-1])

	
#Convert the curve to the same basis
	op_field1 = np.linspace(min_field1_low, min_field1_high, num = in_opt.linspace_data_point)
	ip_field1 = np.linspace(min_field1_low, min_field1_high, num = in_opt.linspace_data_point)
	op_field2 = np.linspace(min_field2_low, min_field2_high, num = in_opt.linspace_data_point)
	ip_field2 = np.linspace(min_field2_low, min_field2_high, num = in_opt.linspace_data_point)

#Interpolate the magnetization data to the basis	
	ip_mag1 = np.interp(ip_field1, field_data_ip1, mag_data_ip1, 
					 left = np.mean(mag_data_ip1[:in_opt.extrapolate_point]), 
					 right = np.mean(mag_data_ip1[-in_opt.extrapolate_point:]))
	op_mag1 = np.interp(op_field1, field_data_op1, mag_data_op1, 
					 left = np.mean(mag_data_op1[:in_opt.extrapolate_point]), 
					 right = np.mean(mag_data_op1[-in_opt.extrapolate_point:]))
	ip_mag2 = np.interp(ip_field2, field_data_ip2, mag_data_ip2, 
					 left = np.mean(mag_data_ip2[:in_opt.extrapolate_point]), 
					 right = np.mean(mag_data_ip2[-in_opt.extrapolate_point:]))
	op_mag2 = np.interp(op_field2, field_data_op2, mag_data_op2, 
					 left = np.mean(mag_data_op2[:in_opt.extrapolate_point]), 
					 right = np.mean(mag_data_op2[-in_opt.extrapolate_point:]))
	
#Average of two op curves
	tmp_fig = plt.figure(figsize=(10, 8))
	plt.plot(op_field1, op_mag1, '-b', label = "OP1" )
	plt.plot(op_field2, op_mag2, '-r',label="OP2")
	plt.legend(loc='upper left',fontsize=12.5)
	plt.xlabel('External Field (Oe)', {'color': 'k', 'fontsize': 25})
	plt.ylabel('Magnetization (emu/cc)', {'color': 'k', 'fontsize': 25})
	plt.tick_params(axis='both', which='major', labelsize=12.5)
	area_op1 = simps(op_mag1, op_field1)
	area_op2 = simps(op_mag2, op_field2)
	print(area_op1)
	print(area_op2)
	op_average = (area_op1 + area_op2) / 2
	print("Avg area op(ergs/cm^3) is equal to ", op_average, "\n")
	
	plt.plot(ip_field1, ip_mag1, '-k', label = "IP1" )
	plt.plot(ip_field2, ip_mag2, '-g',label="IP2")
	plt.legend(loc='upper left',fontsize=12.5)
	plt.xlabel('External Field (Oe)', {'color': 'k', 'fontsize': 25})
	plt.ylabel('Magnetization (emu/cc)', {'color': 'k', 'fontsize': 25})
	plt.tick_params(axis='both', which='major', labelsize=12.5)
	area_ip1 = simps(ip_mag1, ip_field1)
	area_ip2 = simps(ip_mag2, ip_field2)
	print(area_ip1)
	print(area_ip2)
	ip_average = (area_ip1 + area_ip2) / 2
	print("Avg area ip(ergs/cm^3) is equal to ", ip_average)
	
	plt.show(tmp_fig)
	#plt.savefig(os.path.join(in_opt.results_path, in_opt.title + 'Magnetization_Plot_%d_'%i), dpi=300)
	plt.close(tmp_fig)

#Keff is the difference in area of the curves/2 because for one quadrant only
	Avg_Keffective = (op_average - ip_average) / 2
	print("K effective(ergs/cm^3)is equal to ", Avg_Keffective)
	print("K effective(MJ/m^3)is equal to ", Avg_Keffective/10000000)

# run the main function
if __name__ == '__main__':
	test()
	# test()
