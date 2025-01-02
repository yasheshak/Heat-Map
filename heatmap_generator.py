import matplotlib.pyplot as plt
import matplotlib.patches as mplpatches
import numpy as np
import argparse

parser = argparse.ArgumentParser()

#Use assignment 3 argparse and change the two input files
parser.add_argument('-p', '--phase', default='peak_phases.phase', type=str, action='store', help='input phase file')
parser.add_argument('-c', '--expression', default='FPKM_expression_data.exp', type=str, action='store', help='input expression file')
parser.add_argument('-o', '--output', default='genes_CT_heatmap.png', type=str, action='store', help='output png file')

args = parser.parse_args()

phaseFile = args.phase
expFile = args.expression
outFile = args.output

print(phaseFile, expFile, outFile)


plt.style.use('BME163')

#Figure dimensions
figureWidth=5
figureHeight=3 

plt.figure(figsize=(figureWidth,figureHeight))

#main panel dimensions
panel1Width=0.75
panel1Height=2.5
panel1 = plt.axes([0.5/figureWidth,0.3/figureHeight,panel1Width/figureWidth,panel1Height/figureHeight])

#Mini color map guide dimesnions
panel2Width=0.1
panel2Height=0.2
panel2 = plt.axes([1.3/figureWidth,1.45/figureHeight,panel2Width/figureWidth,panel2Height/figureHeight])

#Main panel labels
panel1.set_xlabel('CT')
panel1.set_ylabel('Number of genes')
panel1.set_xticks([ 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5])
panel1.set_xticklabels(['0','','6', '', '12', '', '18', ''])
panel1.set_yticks(range(0,1201,200) )

#colormap labels
panel2.set_xticks([])
panel2.yaxis.tick_right()
panel2.set_yticks([0,1])
panel2.set_yticklabels(['Min', 'Max'])

#PLOTTING COLOR MAP
#Following code was given by Professor Vollmers in Lecture 9
viridis5 = (253/255, 231/255, 37/255)
viridis4 = (94/255, 201/255, 98/255)
viridis3 = (33/255, 145/255, 140/255)
viridis2 = (59/255, 82/255, 139/255)
viridis1 = (68/255, 1/255, 84/255)

R1=np.linspace(viridis1[0],viridis2[0],26)
G1=np.linspace(viridis1[1],viridis2[1],26)
B1=np.linspace(viridis1[2],viridis2[2],26)

R2=np.linspace(viridis2[0],viridis3[0],26)
G2=np.linspace(viridis2[1],viridis3[1],26)
B2=np.linspace(viridis2[2],viridis3[2],26)

R3=np.linspace(viridis3[0],viridis4[0],26)
G3=np.linspace(viridis3[1],viridis4[1],26)
B3=np.linspace(viridis3[2],viridis4[2],26)

R4=np.linspace(viridis4[0],viridis5[0],26)
G4=np.linspace(viridis4[1],viridis5[1],26)
B4=np.linspace(viridis4[2],viridis5[2],26)

R=np.concatenate((R1[:-1],R2[:-1],R3[:-1],R4),axis=None)
G=np.concatenate((G1[:-1],G2[:-1],G3[:-1],G4),axis=None)
B=np.concatenate((B1[:-1],B2[:-1],B3[:-1],B4),axis=None) #Professor Vollmer's code ends here

#Plotting color rectangles on panel2
color_rects = []
num_rects = len(R)

for i in range(num_rects):
    rect = mplpatches.Rectangle((0, i / num_rects), 1, 1 / num_rects,
                             facecolor=(R[i], G[i], B[i]), edgecolor='none')
    color_rects.append(rect)

for rect in color_rects:
    panel2.add_patch(rect)


#Parse through exp file and put into dictionary
exp_dict = {}
with open("BME163_Input_Data_Assignment6.exp", 'r') as exp_file:
    next(exp_file)

    for line in exp_file:
        columns = line.strip().split('\t')

        ensembl_id = columns[1]
        fpkm_values = [float(columns[i]) for i in range(4, 12)]
        exp_dict[ensembl_id] = fpkm_values 

# Normalize the values in the dictionary
normalized_dict = {}
for ensembl_id, values in exp_dict.items():
    min_value = min(values)
    max_value = max(values)
    
    normalized_values = [int(((v - min_value) / (max_value - min_value)) * 100) for v in values]
    normalized_dict[ensembl_id] = normalized_values

      

#Parse through the phase file and put into dictionary
phase_dict = {}
with open("BME163_Input_Data_Assignment6.phase", 'r') as phase_file:
    next(phase_file)

    for line in phase_file:
        columns = line.strip().split('\t')

        ensembl_id = columns[0]
        peak_phase = float(columns[1])
        phase_dict[ensembl_id] = peak_phase


#Sort the normalized dict and only take out sorted values for the next steps
sorted_exp_values = [i for _, i in sorted(normalized_dict.items(), key=lambda x: phase_dict.get(x[0], 0))]
 
#convert phase dict into a list and sort the list and only take our sorted values for the next steps
CT_list = list(phase_dict.items())
sorted_phase_values = [value for _, value in sorted(CT_list, key=lambda x: x[1])]


####PLOTTING NOW
#set the x and y limits based on the values and sorted values
#print(len(values))
#print(len(sorted_exp_values))

panel1.set_xlim(0, 8) #length of values
panel1.set_ylim(0, 1262) #length of values in sorted exp dict

#used lecture 18 code but looped through the sorted lists for value, already normalized earlier, use enumerate for loop inside loop
for ypos, values in enumerate(reversed(sorted_exp_values)):
    for xpos, x in enumerate(values):
        rectangles = mplpatches.Rectangle([xpos, ypos], 1, 1, 
                                         facecolor=(R[x], G[x], B[x]), 
                                         edgecolor='black', 
                                         linewidth=0)
        panel1.add_patch(rectangles)


plt.savefig(outFile,dpi=600)
