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
parser.add_argument("-r","--results",help="path where to write results",default="results_fisher.csv")
parser.add_argument("-v","--verbose", help="increase output verbosity (def. false)")
args = parser.parse_args()

wait_time = 1
run_amira = 0

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


def run_fisher():
    global sample_size
    correct_path = "../fisher/correct/fim_closed"
    enumerate_path = "../fisher/enumerate/fim_closed"
    work_dir_path = "work_dir/"
    if not os.path.exists(work_dir_path):
        os.system("mkdir "+str(work_dir_path))
    temp_file_path = work_dir_path+"out_"+str(args.dataset)+"_sz"+str(args.samplesize)+"_a"+str(args.alpha)

    # create random sample

    if run_amira == 1:
        cmd = "python run_amira.py -db "+str(args.dataset)+" -sz "+str(args.samplesize)+" -wd "+str(work_dir_path)
        print cmd
        os.system(cmd)
        time.sleep(wait_time)

    trans_path = work_dir_path+"sample_amira_"+str(args.dataset)+"_sz"+str(args.samplesize)+"_g"+str(0.01)+".dat"
    labels_path = work_dir_path+"sample_amira_"+str(args.dataset)+"_sz"+str(args.samplesize)+"_g"+str(0.01)+".labels"

    if not os.path.isfile(labels_path):
        print "path to dataset not correct! "+labels_path
        exit()

    if not os.path.isfile(trans_path):
        print "path to dataset not correct! "+trans_path
        exit()

    checkinput_transactions(trans_path)
    time.sleep(wait_time)

    # correct
    cmd = correct_path+" "+temp_file_path+" "+str(args.alpha)+" "+str(labels_path)+" "+str(trans_path)
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
    cmd = enumerate_path+" "+temp_file_path+" "+str(corr_sign_thr)+" "+str(int(testability_support))+" "+str(labels_path)+" "+str(trans_path)+" > "+out_path+" 2>&1"
    print "enumerating "+str(args.dataset)+"..."
    if args.verbose:
        print cmd
        time.sleep(wait_time)
    os.system(cmd)
    time.sleep(wait_time)
    patterns = list()
    results_patterns = ["N = ", "n1 = " ,"Number of significant patterns found: "]
    results_patterns.append("Tested patterns ")
    results_patterns.append("Total execution time (s): ")
    res = list()
    for pattern in results_patterns:
        print "checking "+pattern
        try:
            returned = float(get_result(pattern , out_path))
            res.append(returned)
        except ValueError:
            print "Problem in fisher when checking for "+pattern
            print "DB = "+str(args.db)+" sz = "+str(args.sz)
            res.append(-1.0)
    res.append(time_to_correct)

    print args.results
    fout = open(args.results , 'a')
    out = args.dataset+";"+str(args.samplesize)+";"+str(args.alpha)+";"
    out = out+str(corr_sign_thr)+";"+str(testability_support)+";"
    for res_ in res:
        out = out+str(res_)+";"
    out = out+"\n"
    fout.write(out)


run_fisher()
