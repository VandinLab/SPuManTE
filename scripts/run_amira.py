import math
import os
import numpy as np
import sys
import time
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-db","--dataset", help="dataset name")
parser.add_argument("-sz","--samplesize", type=float ,help="sample size (>0).",default=10000.)
parser.add_argument("-g","--gamma", type=float ,help="confidence for Amira (in (0,alpha))",default=0.01)
parser.add_argument("-v","--verbose", help="increase output verbosity (def. false)")
parser.add_argument("-wd", help="path to the work directory ",default="work_dir/")
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
print "gamma = "+str(args.gamma)


def dataset_add_label(transactions_path , labels_path , transactions_with_label_path):
    transactions = open(transactions_path , 'r')
    labels = open(labels_path , 'r')
    out_transactions = open(transactions_with_label_path , 'w')

    maxitem = 0
    for trans in transactions:
        items = trans.split(' ')
        for item in items:
            try:
                maxitem = max(maxitem , int(item))
            except ValueError:
                pass
    #print "maxitem "+str(maxitem)
    label0id = maxitem + 1
    label1id = maxitem + 2

    transactions.close()
    transactions = open(transactions_path , 'r')
    for trans in transactions:
        label = labels.readline()
        if int(label) == 1:
            trans = trans.replace("\n","")
            trans = trans.replace("\r","")
            trans = trans+" "+str(label1id)+"\n"
        else:
            trans = trans.replace("\n","")
            trans = trans.replace("\r","")
            trans = trans+" "+str(label0id)+"\n"
        out_transactions.write(trans)

    return label0id , label1id

def dataset_remove_label(label0id , label1id , sample_path , sample_labels , transactions_with_label_path):
    transactions_with_label = open(transactions_with_label_path , 'r')
    out_transactions = open(sample_path , 'w')
    out_labels = open(sample_labels , 'w')

    for trans in transactions_with_label:
        label = 0
        trans_ = ""
        trans = trans.replace("\n","")
        trans = trans.replace("\r","")
        items = trans.split(' ')
        #print items
        for item in items:
            if len(item) > 0:
                if int(item) == label1id:
                    label = 1
                else:
                    trans_ = trans_ + str(item) + " "
        trans_ = trans_ + "\n"
        label = str(label)+"\n"
        out_transactions.write(trans_)
        out_labels.write(label)



def run_amira():
    global sample_size
    amira_path = "../amira/amira"
    work_dir_path = args.wd
    if not os.path.exists(work_dir_path):
        os.system("mkdir "+str(work_dir_path))
    temp_file_path = work_dir_path+"out_amira_"+str(args.dataset)+"_sz"+str(args.samplesize)+"_g"+str(args.gamma)+".txt"

    labels_path = "../datasets/"+args.dataset+"/"+args.dataset+".labels"
    trans_path = "../datasets/"+args.dataset+"/"+args.dataset+".dat"
    transactions_with_label_path = work_dir_path+"dataset_"+str(args.dataset)+"_sz"+str(args.samplesize)+"_g"+str(args.gamma)+"_with_labels.dat"

    if not os.path.isfile(trans_path):
        print "path to dataset not correct! "+trans_path
        exit()

    if not os.path.isfile(labels_path):
        print "path to dataset not correct! "+labels_path
        exit()

    label0id , label1id = dataset_add_label(trans_path , labels_path , transactions_with_label_path)

    #checkinput_transactions(trans_path)
    time.sleep(wait_time)

    if not os.path.isfile(transactions_with_label_path):
        print "path to dataset not correct! "+transactions_with_label_path
        exit()

    # correct
    min_freq = max(0.1 , math.sqrt(1.0 / float(args.samplesize)))
    sample_path_with_labels = work_dir_path+"sample_amira_"+str(args.dataset)+"_sz"+str(args.samplesize)+"_g"+str(args.gamma)+"_with_labels.dat"
    sample_path = work_dir_path+"sample_amira_"+str(args.dataset)+"_sz"+str(args.samplesize)+"_g"+str(args.gamma)+".dat"
    sample_labels_path = work_dir_path+"sample_amira_"+str(args.dataset)+"_sz"+str(args.samplesize)+"_g"+str(args.gamma)+".labels"
    cmd = amira_path+" -f -n -p -s "+sample_path_with_labels+" "+str(args.gamma)+" "+str(min_freq)+" "+str(args.samplesize)+" "+str(transactions_with_label_path)+" > "+temp_file_path
    print "creating sample and computing epsilon for "+str(args.dataset)+"..."
    if args.verbose:
        print cmd
        time.sleep(wait_time)
    os.system(cmd)
    time.sleep(wait_time)

    # now remove labels form sample
    dataset_remove_label(label0id , label1id , sample_path , sample_labels_path , sample_path_with_labels)

    results_patterns = ("eps: " , "total: ","create_sample: ")
    results = list()
    for pattern in results_patterns:
        results.append(float(get_result(pattern , temp_file_path)))

    epsilon = results[0]
    running_time = results[1]
    print epsilon
    fout = open(temp_file_path , 'w')
    i = 0
    for pattern in results_patterns:
        fout.write(pattern+str(results[i])+"\n")
        i = i + 1


run_amira()
