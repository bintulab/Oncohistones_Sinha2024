#this script crawls into dedupscaled folder and runs SEACR peak caller on all processed .bedgraph CUT&RUN data 
#(base) bintulab@BintuDigitalStorm-PC:~/win_e/Joydeb/20210700_Joydeb_CUTRUN_pilot_H33$ python scripts/run_seacr.py


#note, improve code vby adding arguments for filepaths directly and add heatma


import os
import subprocess


inputdir= "/home/bintulab/win_f/Joydeb/20230124_CNR_6_0/SEACR_input/"
print(inputdir)

for files in os.listdir(inputdir) :
    #process files command line
    if files.endswith('.bedgraph') :
    #run seacr peak calling
        fname=print(files)
        type(fname)
        os.chdir("..")
        print("Running SEACR on file:", files)
        subprocess.run(['/home/bintulab/win_f/Joydeb/20230124_CNR_6_0/scripts/SEACR/SEACR_1.3.sh', inputdir+files, '0.01', 'non', 'stringent', 'PEAKS_'+files]) 
        print("Done!")
        print()
        print()
        os.chdir(inputdir)




