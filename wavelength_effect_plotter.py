import matplotlib.pyplot as plt
import numpy as np
import csv
import argparse
from matplotlib.backends.backend_pdf import PdfPages

def get_data(filename):
    arr=[]
    with open(filename) as csvfile:
        csvreader=csv.reader(csvfile,delimiter="|")
        for row in csvreader:
            arr.append(row)
    C_l=np.array(arr[0]).astype(np.float)
    C_d=np.array(arr[1]).astype(np.float)
    Cl_upon_Cd=np.array(C_l)/np.array(C_d)
    return(C_l,C_d,Cl_upon_Cd)


parser=argparse.ArgumentParser()
parser.add_argument("pdf_output_name",type=str)
parser.add_argument("function",type=str,choices=['Sine','Triangle','Square'],default='Sine',help='Select the (Sine/Square/Triangle Wave function')
parser.add_argument("amplitude", type=str)
parser.add_argument("wavelength",type=str,nargs='+',choices=['08','16','32','64'])
parser.add_argument("-b","--baseline",action='store_false',default=True)
args=parser.parse_args()
print(args.pdf_output_name)
print(args.function)
print(args.wavelength)
print(args.amplitude)
print(args.baseline)



function=args.function
wavelength=args.wavelength
amplitude=args.amplitude

filenames=[]

for w in wavelength:
	filenames.append('Simulations_'+function+'/'+function+'_A_0.00'+amplitude+'_w_'+w+'.csv')
print(filenames)
C_l=[];C_d=[];Cl_upon_Cd=[]
for filename in filenames:
	cl,cd,clcd=get_data(filename)
	C_l.append(cl)
	C_d.append(cd)
	Cl_upon_Cd.append(clcd)

bcl,bcd,bclcd=get_data('Major_Project_Simulations/A_0.00_w_s_upon_00.csv')
print(bcl)
x=np.linspace(1,12,12)
with PdfPages(args.pdf_output_name) as pdf:
	plt.figure()
	if args.baseline:
		plt.plot(x,bcl,linestyle='dashed',marker='.',label='Baseline Wing',alpha=0.8,linewidth=0.8)
		print('Cl baseline:',bcl[0],'-------',bcl[-1])
	for i,filename in enumerate(filenames):
		plt.plot(x,C_l[i],linestyle='dashed',marker='.',label=function+' A'+amplitude+r'$ \lambda$'+wavelength[i],alpha=0.8,linewidth=0.8)
	print('Cl :',C_l[-1][0],'-------',C_l[-1][-1])
	plt.xlabel('AoA (degrees)')
	plt.ylabel('$C_{L}$')
	# plt.title('Effect of Wavelength on $C_{l}$')
	plt.legend()
	pdf.savefig(bbox_inches = 'tight')
	plt.close()
	plt.figure()
	if args.baseline:
		plt.plot(x,bcd,linestyle='dashed',marker='.',label='Baseline Wing',alpha=0.8,linewidth=0.8)
		print('Cd baseline:',bcd[0],'-------',bcd[-1])
	for i,filename in enumerate(filenames):
		plt.plot(x,C_d[i],linestyle='dashed',marker='.',label=function+' A'+amplitude+r'$ \lambda$'+wavelength[i],alpha=0.8,linewidth=0.8)
	print('Cd :',C_d[-1][0],'-------',C_d[-1][-1])
	plt.xlabel('AoA (degrees)')
	plt.ylabel('$C_{Di}$')
	# plt.title('Effect of Wavelength on $C_{d}$')
	plt.legend()
	pdf.savefig(bbox_inches = 'tight')
	plt.close()
	plt.figure()
	if args.baseline:
		plt.plot(x,bclcd,linestyle='dashed',marker='.',label='Baseline Wing',alpha=0.8,linewidth=0.8)
		print('ClCd baseline:',bclcd[0],'-------',bclcd[-1])
	for i,filename in enumerate(filenames):
		plt.plot(x,Cl_upon_Cd[i],linestyle='dashed',marker='.',label=function+' A'+amplitude+r'$ \lambda$'+wavelength[i],alpha=0.8,linewidth=0.8)
	print('Cd :',Cl_upon_Cd[-1][0],'-------',Cl_upon_Cd[-1][-1])
	plt.xlabel('AoA (degrees)')
	plt.ylabel('$C_{L}/C_{Di}$')
	# plt.title('Effect of Wavelength on $C_{l}/C_{d}$')
	plt.legend()
	pdf.savefig(bbox_inches = 'tight')
	plt.close()