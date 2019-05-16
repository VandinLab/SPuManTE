import math
import numpy as np
import sys
import os
import re
import time


datasets = ["T10I4D100K","a9a","bms-pos","bms-web2","chess","connect","mushroom","pumb-star","susy","T40I10D100K","accidents","bms-web1","breast-cancer","cod-rna","covtype","ijcnn1","phishing","retail","svmguide3"]

datasets = ["a9a","mushroom","breast-cancer","retail","ijcnn1","svmguide3","cod-rna","covtype"]

sample_sizes = [1000, 10000, 100000, 1000000]
small_wait = 5
medium_wait = 10
long_wait = 30
alpha = 0.05
gamma = 0.01

run_fisher = 1
parallel = 1

runs = 1

def run_experiment(db , sz):
    for run_id in range(runs):
        time.sleep(np.random.randint(small_wait)+small_wait)
        cmd1 = "python run_unconditional.py -db "+str(db)+" -sz "+str(sz)+" && "
        cmd2 = "python run_fisher.py -db "+str(db)+" -sz "+str(sz)
        cmd = cmd1+cmd2
        print cmd
        os.system(cmd)
        time.sleep(np.random.randint(medium_wait)+medium_wait)

if parallel == 1:
    ids = list()
    from multiprocessing import Pool
    pool = Pool()
    for sz in sample_sizes:
        for db in datasets:
            ids.append( pool.apply_async(run_experiment , [db , sz]))

    todo_total = len(ids)
    done_ = 0
    for id in ids:
        id.get()
        done_ = done_ + 1
        print "DONE "+str(done_)+"/"+str(todo_total)
else:
    for sz in sample_sizes:
        for db in datasets:
            cmd = "python run_unconditional.py -db "+str(db)+" -sz "+str(sz)
            print cmd
            os.system(cmd)
            time.sleep(medium_wait)
            if run_fisher == 1:
                cmd = "python run_fisher.py -db "+str(db)+" -sz "+str(sz)
                print cmd
                os.system(cmd)
                time.sleep(medium_wait)
