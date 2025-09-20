from smoldyn import Simulation as S
import numpy as np
from pathlib import Path
import math
import os
import sys

def Count_SPEN(t, args):
    global active_receptor_count, curr_t, count1, count2, count3
    curr_t = t # Current time point

    # chromosome centers
    chr1_Xist_center = np.array([-5,0,0])
    chr2_Xist_center = np.array([0,0,0])
    chr3_Xist_center = np.array([5,0,0])
    
    # We counted bound SPEN within 2.5 um to the chromosome center to be considered on the SAME chromosome
    dis_threshold = 2.5

    # Count the number of bound SPEN in each output files.

    # Output files containing Xist bound to one SPEN
    with open(cwd+"/"+outputfilename_XS,'rb') as f:
        f.seek(-2,os.SEEK_END)
        while f.tell() > 0 and f.read(1) != b'\n':
            f.seek(-2, os.SEEK_CUR) # Read the last line of the file
        float_line = f.readline().decode('utf-8').strip().split() # Extract numbers
        float_line = list(map(float, float_line))
        XS_x = float_line[1:len(float_line):3] # x coordinates
        XS_y = float_line[2:len(float_line):3] # y coordinates
        XS_z = float_line[3:len(float_line):3] # z coordinates

    # Output files containing Xist bound to two SPENs
    with open(cwd+"/"+outputfilename_XS2,'rb') as f:
        f.seek(-2,os.SEEK_END)
        while f.tell() > 0 and f.read(1) != b'\n':
            f.seek(-2, os.SEEK_CUR)
        float_line = f.readline().decode('utf-8').strip().split()
        float_line = list(map(float, float_line))
        XS2_x = float_line[1:len(float_line):3]
        XS2_y = float_line[2:len(float_line):3]
        XS2_z = float_line[3:len(float_line):3]

    # Output files containing Xist bound to three SPENs
    with open(cwd+"/"+outputfilename_XS3,'rb') as f:
        f.seek(-2,os.SEEK_END)
        while f.tell() > 0 and f.read(1) != b'\n':
            f.seek(-2, os.SEEK_CUR)
        float_line = f.readline().decode('utf-8').strip().split()
        float_line = list(map(float, float_line))
        XS3_x = float_line[1:len(float_line):3]
        XS3_y = float_line[2:len(float_line):3]
        XS3_z = float_line[3:len(float_line):3]

    # Output files containing Xist bound to four SPENs
    with open(cwd+"/"+outputfilename_XS4,'rb') as f:
        f.seek(-2,os.SEEK_END)
        while f.tell() > 0 and f.read(1) != b'\n':
            f.seek(-2, os.SEEK_CUR)
        float_line = f.readline().decode('utf-8').strip().split()
        float_line = list(map(float, float_line))
        XS4_x = float_line[1:len(float_line):3]
        XS4_y = float_line[2:len(float_line):3]
        XS4_z = float_line[3:len(float_line):3]

    # Output files containing Xist bound to five SPENs
    with open(cwd+"/"+outputfilename_XS5,'rb') as f:
        f.seek(-2,os.SEEK_END)
        while f.tell() > 0 and f.read(1) != b'\n':
            f.seek(-2, os.SEEK_CUR)
        float_line = f.readline().decode('utf-8').strip().split()
        float_line = list(map(float, float_line))
        XS5_x = float_line[1:len(float_line):3]
        XS5_y = float_line[2:len(float_line):3]
        XS5_z = float_line[3:len(float_line):3]

    # Output files containing Xist bound to six SPENs
    with open(cwd+"/"+outputfilename_XS6,'rb') as f:
        f.seek(-2,os.SEEK_END)
        while f.tell() > 0 and f.read(1) != b'\n':
            f.seek(-2, os.SEEK_CUR)
        float_line = f.readline().decode('utf-8').strip().split()
        float_line = list(map(float, float_line))
        XS6_x = float_line[1:len(float_line):3]
        XS6_y = float_line[2:len(float_line):3]
        XS6_z = float_line[3:len(float_line):3]

    # Output files containing Xist bound to seven SPENs
    with open(cwd+"/"+outputfilename_XS7,'rb') as f:
        f.seek(-2,os.SEEK_END)
        while f.tell() > 0 and f.read(1) != b'\n':
            f.seek(-2, os.SEEK_CUR)
        float_line = f.readline().decode('utf-8').strip().split()
        float_line = list(map(float, float_line))
        XS7_x = float_line[1:len(float_line):3]
        XS7_y = float_line[2:len(float_line):3]
        XS7_z = float_line[3:len(float_line):3]

    # Output files containing Xist bound to eight SPENs
    with open(cwd+"/"+outputfilename_XS8,'rb') as f:
        f.seek(-2,os.SEEK_END)
        while f.tell() > 0 and f.read(1) != b'\n':
            f.seek(-2, os.SEEK_CUR)
        float_line = f.readline().decode('utf-8').strip().split()
        float_line = list(map(float, float_line))
        XS8_x = float_line[1:len(float_line):3]
        XS8_y = float_line[2:len(float_line):3]
        XS8_z = float_line[3:len(float_line):3]

    # Output files containing Xist bound to nine SPENs
    with open(cwd+"/"+outputfilename_XS9,'rb') as f:
        f.seek(-2,os.SEEK_END)
        while f.tell() > 0 and f.read(1) != b'\n':
            f.seek(-2, os.SEEK_CUR)
        float_line = f.readline().decode('utf-8').strip().split()
        float_line = list(map(float, float_line))
        XS9_x = float_line[1:len(float_line):3]
        XS9_y = float_line[2:len(float_line):3]
        XS9_z = float_line[3:len(float_line):3]

    # Output files containing Xist bound to ten SPENs
    with open(cwd+"/"+outputfilename_XS10,'rb') as f:
        f.seek(-2,os.SEEK_END)
        while f.tell() > 0 and f.read(1) != b'\n':
            f.seek(-2, os.SEEK_CUR)
        float_line = f.readline().decode('utf-8').strip().split()
        float_line = list(map(float, float_line))
        XS10_x = float_line[1:len(float_line):3]
        XS10_y = float_line[2:len(float_line):3]
        XS10_z = float_line[3:len(float_line):3]

    # Calculate distances of SPEN molecules from the centers of chromosomes

    all_coordinates_XS = np.column_stack((XS_x, XS_y, XS_z))
    XS_distance1 = np.linalg.norm(all_coordinates_XS - chr1_Xist_center, axis=1)
    XS_distance2 = np.linalg.norm(all_coordinates_XS - chr2_Xist_center, axis=1)
    XS_distance3 = np.linalg.norm(all_coordinates_XS - chr3_Xist_center, axis=1)

    all_coordinates_XS2 = np.column_stack((XS2_x, XS2_y, XS2_z))
    XS2_distance1 = np.linalg.norm(all_coordinates_XS2 - chr1_Xist_center, axis=1)
    XS2_distance2 = np.linalg.norm(all_coordinates_XS2 - chr2_Xist_center, axis=1)
    XS2_distance3 = np.linalg.norm(all_coordinates_XS2 - chr3_Xist_center, axis=1)

    all_coordinates_XS3 = np.column_stack((XS3_x, XS3_y, XS3_z))
    XS3_distance1 = np.linalg.norm(all_coordinates_XS3 - chr1_Xist_center, axis=1)
    XS3_distance2 = np.linalg.norm(all_coordinates_XS3 - chr2_Xist_center, axis=1)
    XS3_distance3 = np.linalg.norm(all_coordinates_XS3 - chr3_Xist_center, axis=1)

    all_coordinates_XS4 = np.column_stack((XS4_x, XS4_y, XS4_z))
    XS4_distance1 = np.linalg.norm(all_coordinates_XS4 - chr1_Xist_center, axis=1)
    XS4_distance2 = np.linalg.norm(all_coordinates_XS4 - chr2_Xist_center, axis=1)
    XS4_distance3 = np.linalg.norm(all_coordinates_XS4 - chr3_Xist_center, axis=1)

    all_coordinates_XS5 = np.column_stack((XS5_x, XS5_y, XS5_z))
    XS5_distance1 = np.linalg.norm(all_coordinates_XS5 - chr1_Xist_center, axis=1)
    XS5_distance2 = np.linalg.norm(all_coordinates_XS5 - chr2_Xist_center, axis=1)
    XS5_distance3 = np.linalg.norm(all_coordinates_XS5 - chr3_Xist_center, axis=1)

    all_coordinates_XS6 = np.column_stack((XS6_x, XS6_y, XS6_z))
    XS6_distance1 = np.linalg.norm(all_coordinates_XS6 - chr1_Xist_center, axis=1)
    XS6_distance2 = np.linalg.norm(all_coordinates_XS6 - chr2_Xist_center, axis=1)
    XS6_distance3 = np.linalg.norm(all_coordinates_XS6 - chr3_Xist_center, axis=1)

    all_coordinates_XS7 = np.column_stack((XS7_x, XS7_y, XS7_z))
    XS7_distance1 = np.linalg.norm(all_coordinates_XS7 - chr1_Xist_center, axis=1)
    XS7_distance2 = np.linalg.norm(all_coordinates_XS7 - chr2_Xist_center, axis=1)
    XS7_distance3 = np.linalg.norm(all_coordinates_XS7 - chr3_Xist_center, axis=1)

    all_coordinates_XS8 = np.column_stack((XS8_x, XS8_y, XS8_z))
    XS8_distance1 = np.linalg.norm(all_coordinates_XS8 - chr1_Xist_center, axis=1)
    XS8_distance2 = np.linalg.norm(all_coordinates_XS8 - chr2_Xist_center, axis=1)
    XS8_distance3 = np.linalg.norm(all_coordinates_XS8 - chr3_Xist_center, axis=1)

    all_coordinates_XS9 = np.column_stack((XS9_x, XS9_y, XS9_z))
    XS9_distance1 = np.linalg.norm(all_coordinates_XS9 - chr1_Xist_center, axis=1)
    XS9_distance2 = np.linalg.norm(all_coordinates_XS9 - chr2_Xist_center, axis=1)
    XS9_distance3 = np.linalg.norm(all_coordinates_XS9 - chr3_Xist_center, axis=1)

    all_coordinates_XS10 = np.column_stack((XS10_x, XS10_y, XS10_z))
    XS10_distance1 = np.linalg.norm(all_coordinates_XS10 - chr1_Xist_center, axis=1)
    XS10_distance2 = np.linalg.norm(all_coordinates_XS10 - chr2_Xist_center, axis=1)
    XS10_distance3 = np.linalg.norm(all_coordinates_XS10 - chr3_Xist_center, axis=1)

    # Count the number of SPEN on chromosome X1
    count1 = (np.sum(XS_distance1 < dis_threshold) + 2*np.sum(XS2_distance1 < dis_threshold) + 3*np.sum(XS3_distance1 < dis_threshold)
    + 4*np.sum(XS4_distance1 < dis_threshold) + 5*np.sum(XS5_distance1 < dis_threshold) + 6*np.sum(XS6_distance1 < dis_threshold)
    + 7*np.sum(XS7_distance1 < dis_threshold)  + 8*np.sum(XS8_distance1 < dis_threshold) + 9*np.sum(XS9_distance1 < dis_threshold)
    + 10*np.sum(XS10_distance1 < dis_threshold))

    # Count the number of SPEN on chromosome X2
    count2 = (np.sum(XS_distance2 < dis_threshold) + 2*np.sum(XS2_distance2 < dis_threshold) + 3*np.sum(XS3_distance2 < dis_threshold)
    + 4*np.sum(XS4_distance2 < dis_threshold) + 5*np.sum(XS5_distance2 < dis_threshold) + 6*np.sum(XS6_distance2 < dis_threshold)
    + 7*np.sum(XS7_distance2 < dis_threshold)  + 8*np.sum(XS8_distance2 < dis_threshold) + 9*np.sum(XS9_distance2 < dis_threshold)
    + 10*np.sum(XS10_distance2 < dis_threshold))

    # Count the number of SPEN on chromosome X3
    count3 = (np.sum(XS_distance3 < dis_threshold) + 2*np.sum(XS2_distance3 < dis_threshold) + 3*np.sum(XS3_distance3 < dis_threshold)
              + 4*np.sum(XS4_distance3 < dis_threshold) + 5*np.sum(XS5_distance3 < dis_threshold) + 6*np.sum(XS6_distance3 < dis_threshold)
              + 7*np.sum(XS7_distance3 < dis_threshold)  + 8*np.sum(XS8_distance3 < dis_threshold) + 9*np.sum(XS9_distance3 < dis_threshold)
              + 10*np.sum(XS10_distance3 < dis_threshold))

    # Activator synthesis rate, calculated using number of bound SPEN. See model equations in paper.
    act_syn = a_act/(1 +(count1/K_n)**n) + a_act/(1 +(count2/K_n)**n) + a_act/(1 +(count3/K_n)**n)
    
    # Add the command that controls activator production rate regulated by SPEN at each time step
    add_command = "set reaction_production Activator_production " + str(act_syn*dt)

    # The time when the added command takes effect
    ontime = curr_t + dt
    s.addCommand(add_command,cmd_type="@",on=ontime)

    # Calculate Xist dissociation rate from the chromosome regulated by bound SPEN
    Xist_degrad1 = k2/(1+(count1/K_S)**m)
    Xist_degrad2 = k2/(1+(count2/K_S)**m)
    Xist_degrad3 = k2/(1+(count3/K_S)**m)

    # Similarly, add the command that controls Xist dissociation rate on Chromosome X1, X2, and X3
    add_command = "set reaction_rate Xistb_dissociation_1 " + str(Xist_degrad1)
    s.addCommand(add_command,cmd_type="@",on=ontime)
    
    add_command = "set reaction_rate Xistb_dissociation_2 " + str(Xist_degrad2)
    s.addCommand(add_command,cmd_type="@",on=ontime)
    
    add_command = "set reaction_rate Xistb_dissociation_3 " + str(Xist_degrad3)
    s.addCommand(add_command,cmd_type="@",on=ontime)

    # Read all activator coordinates and calculate the number of total activators
    with open(cwd+"/"+outputfilename_Activator,'rb') as f:
        f.seek(-2,os.SEEK_END)
        while f.tell() > 0 and f.read(1) != b'\n':
            f.seek(-2, os.SEEK_CUR)
        float_line = f.readline().decode('utf-8').strip().split()
        float_line = list(map(float, float_line))
        act_x = float_line[1:len(float_line):3]

    n_act = len(act_x)

    # Calculate Xist synthesis rate, regulated by activator    
    Xist_synthesis_rate = a_x*n_act/(K_a + n_act)
    add_command = "set reaction_rate XIST_transcription " + str(Xist_synthesis_rate)
    s.addCommand(add_command,cmd_type="@",on=ontime)

    return 0

# Not used. Just for the required format of s.connect
def Rate_change(var):
    return 0


def build_smoldyn():
    global a_act, K_n, n, k2, K_S, m, a_x, K_a, dt, t, cwd, outputfilename_Activator, outputfilename_XS, outputfilename_XS2, outputfilename_XS3, outputfilename_XS4, outputfilename_XS5, outputfilename_XS6, outputfilename_XS7, outputfilename_XS8, outputfilename_XS9, outputfilename_XS10, activator_filename, s, count1, count2, count3
    dt = float(sys.argv[1]) # time step
    filename = sys.argv[2] + ".cfg" # Smoldyn file name in which we will add additional commands in each time step
    K_S = float(sys.argv[3]) # quantity of SPEN at which Xist dissociation is half max
    a_act = float(sys.argv[4]) # activator synthesis rate from single X chromosome
    K_n = float(sys.argv[5]) # quantity of bound SPEN at which activator synthesis rate is half max.
    n = float(sys.argv[6]) # Hill coefficient for SPEN supressing activator synthesis rate. 
    k2 = float(sys.argv[7]) # maximum dissociation rate for Xist
    m = float(sys.argv[8]) # Hill coefficient for bound SPEN reducing dissociation rate Xist
    a_x = float(sys.argv[9]) # xist synthesis rate from single X chromosome
    K_a = float(sys.argv[10]) # quantity of activator at which Xist transcription is half max

    # Output file name for activator coordinates
    outputfilename_Activator = "Activator_" + sys.argv[2] + ".txt"

    # Output file names for SPEN coordinates (in complexes with Xist)
    outputfilename_XS = "XS_" + sys.argv[2] + ".txt"
    outputfilename_XS2 = "XS2_" + sys.argv[2] + ".txt"
    outputfilename_XS3 = "XS3_" + sys.argv[2] + ".txt"
    outputfilename_XS4 = "XS4_" + sys.argv[2] + ".txt"
    outputfilename_XS5 = "XS5_" + sys.argv[2] + ".txt"
    outputfilename_XS6 = "XS6_" + sys.argv[2] + ".txt"
    outputfilename_XS7 = "XS7_" + sys.argv[2] + ".txt"
    outputfilename_XS8 = "XS8_" + sys.argv[2] + ".txt"
    outputfilename_XS9 = "XS9_" + sys.argv[2] + ".txt"
    outputfilename_XS10 = "XS10_" + sys.argv[2] + ".txt"

    # Execute Smoldyn functions at every 500 time steps until final simulation time point = 20000 min
    cwd = os.getcwd()
    s=S.fromFile(Path(__file__).parent / filename,"")
    s.connect(Count_SPEN, Rate_change, step=500)
    s.run(stop=20000,dt=dt,start=0, display=False, overwrite=True)
    
def main():
    build_smoldyn()

if __name__ == "__main__":
    main()
