import multiprocessing
import subprocess
import sys
#from datetime import datetime, timedelta
import pandas
import numpy
#import requests
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import Align

'''THIS SCRIPT IMPLEMENTS PARALELL PROCESSING (CURRENTLY HARD SET TO 6 THREADS) 
TO PROCESS LARGE BATCHES OF SEQUENCES. IT IS CURRENTLY SET TO RUN WITH 'ncbi_seqs321.csv'
WHICH CONTAINS METADATA FROM NCBI VIRUS THAT HAS BEEN MERGED WITH A FATSA FILE 
OF THE SAME SEQUENCES. THE SCRIPT SPLITS THE SEQUENCE FILE AND CALLS A PROCESS FOR EACH PART
THE PROCESSES THEN OUTPUT DATA INTO THEIR OWN FILES.
'''


cpu_count = multiprocessing.cpu_count()
#print(cpu_count)

cpu_count = 6 # why does the above think my computer has 12 cores?


'''UPDATE WITH DATA SOURCE!!'''
#GETS THE SEQUENCES TABLE
sequences = pandas.read_csv("ncbi_seqs321.csv")
pandas.set_option("precision",0)



size = len(sequences.index)
end = -1

# OPENS SUBPROCESSES
for x in range(1, cpu_count+1):
    start = end + 1
    end =  start +(size // cpu_count)
    end = min(end, size)
    s = sequences.loc[start:end]
    s.to_csv("seq"+str(x)+".csv")
    subprocess.Popen(args=["python",'process.py',str(x),str(start),str(end)])


print("processes have all started")
