import math
import os
import numpy as np
import sys
import time
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-db","--dataset", help="dataset name")
parser.add_argument("-sz","--samplesize", type=float ,help="sample size (>0).",default=10000.)
parser.add_argument("-a","--alpha", type=float ,help="Upper bound on FWER (in (0,1))",default=0.05)
parser.add_argument("-g","--gamma", type=float ,help="confidence for Amira (in (0,alpha))",default=0.01)
parser.add_argument("-r","--results",help="path where to write results",default="results_unconditional.csv")
parser.add_argument("-v","--verbose", help="increase output verbosity (def. false)")
args = parser.parse_args()

wait_time = 1

def get_result(pattern , path ,  verbose=1):
    fin = open(path,'r')
    for line in fin:
        if pattern in line:
            line = line.replace('\n','')
            if verbose == 1:
                print line
            return line[len(pattern):]
    fin.close()

if not args.samplesize:
    print "Argument samplesize is needed"
    parser.print_help(sys.stderr)
    exit()
if not args.dataset:
    print "dataset name is needed!"
    parser.print_help(sys.stderr)
    exit()

print "dataset = "+str(args.dataset)
print "sample size = "+str(args.samplesize)
print "alpha = "+str(args.alpha)
print "gamma = "+str(args.gamma)

def checkinput_transactions(path):
    fin = open(path , 'r')
    t_id = -1
    for line in fin:
        line = line.replace('\n','')
        t_id = t_id + 1
        items = line.split(' ')
        trans_items = set()
        for item in items:
            if len(item) > 0:
                item = int(item)
                if item not in trans_items:
                    trans_items.add(item)
                else:
                    print "Error! item "+str(item)+" is at least doubled in transaction "+str(t_id)
                    exit()
    print "checked "+path


def run_unconditional():
    global sample_size
    correct_path = "../unconditional/correct/fim_closed"
    enumerate_path = "../unconditional/enumerate/fim_closed"
    work_dir_path = "work_dir/"
    if not os.path.exists(work_dir_path):
        os.system("mkdir "+str(work_dir_path))
    temp_file_path = work_dir_path+"out_"+str(args.dataset)+"_sz"+str(args.samplesize)+"_a"+str(args.alpha)+"_g"+str(args.gamma)

    # run amira to create the sample

    amira_file_path = work_dir_path+"out_amira_"+str(args.dataset)+"_sz"+str(args.samplesize)+"_g"+str(args.gamma)+".txt"
    cmd = "python run_amira.py -db "+str(args.dataset)+" -sz "+str(args.samplesize)
    print cmd
    os.system(cmd)
    time.sleep(wait_time)

    results_patterns = ("eps: " , "total: ","create_sample: ")
    amira_out_path = work_dir_path+"out_amira_"+str(args.dataset)+"_sz"+str(args.samplesize)+"_g"+str(args.gamma)+".txt"
    results = list()
    for pattern in results_patterns:
        results.append(float(get_result(pattern , amira_out_path)))

    epsilon = results[0]
    running_time_amira = (results[1]-results[2]) / 1000.

    trans_path = work_dir_path+"sample_amira_"+str(args.dataset)+"_sz"+str(args.samplesize)+"_g"+str(args.gamma)+".dat"
    labels_path = work_dir_path+"sample_amira_"+str(args.dataset)+"_sz"+str(args.samplesize)+"_g"+str(args.gamma)+".labels"

    checkinput_transactions(trans_path)
    time.sleep(wait_time)

    if not os.path.isfile(labels_path):
        print "path to dataset not correct! "+labels_path
        exit()

    if not os.path.isfile(trans_path):
        print "path to dataset not correct! "+trans_path
        exit()

    # correct
    cmd = correct_path+" "+temp_file_path+" "+str(args.alpha-args.gamma)+" "+str(labels_path)+" "+str(trans_path)+" "+str(epsilon)
    print "correcting "+str(args.dataset)+"..."
    if args.verbose:
        print cmd
        time.sleep(wait_time)
    os.system(cmd)
    time.sleep(wait_time)

    results_patterns = ("\t Corrected significance threshold: " , "\t Final LCM support: ")
    out_path = temp_file_path+"_results.txt"
    res = list()
    for pattern in results_patterns:
        res.append(float(get_result(pattern , out_path)))

    results_patterns = ["Total execution time (s): "]
    out_path = temp_file_path+"_timing.txt"
    for pattern in results_patterns:
        res.append(float(get_result(pattern , out_path)))

    corr_sign_thr = res[0]
    testability_support = res[1]
    time_to_correct = res[2]

    out_path = temp_file_path+"temp_out.txt"
    cmd = enumerate_path+" "+temp_file_path+" "+str(corr_sign_thr)+" "+str(int(testability_support))+" "+str(labels_path)+" "+str(trans_path)+" "+str(epsilon)+" > "+out_path+" 2>&1"
    print "enumerating "+str(args.dataset)+"..."
    if args.verbose:
        print cmd
        time.sleep(wait_time)
    os.system(cmd)
    time.sleep(wait_time)
    patterns = list()
    results_patterns = ["N = ", "n1 = ","Number of significant patterns found: "]
    results_patterns.append("Tested patterns ")
    results_patterns.append("explored_contingency_tables ")
    results_patterns.append("count_number_of_ci_rejects ")
    results_patterns.append("count_number_of_simpleub_rejects ")
    results_patterns.append("count_number_of_enumerate_rejects ")
    results_patterns.append("count_number_of_enumerate_nonrejects ")
    results_patterns.append("Total execution time (s): ")
    res = list()
    for pattern in results_patterns:
        #print "checking "+pattern
        try:
            returned = float(get_result(pattern , out_path))
            res.append(returned)
        except ValueError:
            print "Problem in unconditional when checking for "+pattern
            print "DB = "+str(args.db)+" sz = "+str(args.sz)
            res.append(-1.0)
    res.append(running_time_amira)
    res.append(epsilon)
    res.append(time_to_correct)

    print args.results
    fout = open(args.results , 'a')
    out = args.dataset+";"+str(args.samplesize)+";"+str(args.alpha)+";"+str(args.gamma)+";"
    out = out+str(corr_sign_thr)+";"+str(testability_support)+";"
    for res_ in res:
        out = out+str(res_)+";"
    out = out+"\n"
    fout.write(out)


run_unconditional()
