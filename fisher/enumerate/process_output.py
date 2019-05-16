# Import most important libraries
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import hypergeom


# ----------------------------- FILE I/O FUNCTIONS -----------------------------

# This function takes as input the file containing all significant
# itemsets and returns a list with each significant itemset as a Python
# set data type. Also, a list with the corresponding frequencies is returned,
# for consistency checks
def load_significant_itemsets(significant_itemsets_filename):
    f = open(significant_itemsets_filename,'r')
    significant_itemsets = list()
    frequencies = list()
    for line in f:
        # Itemsets are stored in the file as lines of the form:
        # item1 item2 ... itemN (x)
        # First I split the string by spaces, getting a list of strings:
        # ['item1','item2',...,'itemN','(x)']
        split_line = line.split()
        # List comprehensions in Python are super useful. If you don't know
        # about them check them out in the Python doc! 
        # Here I use them to convert each string 'itemJ' to an integer,
        # discarding the last element of the list, which is the frequency
        # of the itemset, x, and not an item of the itemset
        itemset = [int(i) for i in split_line[:-1]]
        # Add itemset to list of significant itemsets
        significant_itemsets.append(set(itemset))
        # Parse frequency to an integer, removing the first and last
        # characters, which are parentheses
        frequencies.append(int(split_line[-1][1:-1]))
    # Close the file
    f.close()
    # Return both lists
    return (significant_itemsets,frequencies)

# This function takes as input the file containing the cell counts and pvalues
# corresponding to all significant itemsets (see function above). It returns a
# Numpy matrix with 3 columns: the first one is the cell count a, the second
# one the margin x, and the last column the corresponding p-value.
# Note that the second column should coincide with the list called frequencies
# in the function above. If not, something went wrong!
def load_significant_itemset_pvalues(pvalues_filename):
    # This file is easier to read, since it consists of a set of rows all
    # with 3 columns, separated by commas. Therefore, one can use the
    # Numpy built-in function loadtxt to read the file comfortably
    X = np.loadtxt(pvalues_filename,delimiter=',',skiprows=1);
    return X

# This function loads a transaction database in memory. It takes as input the
# file containing the transactions. The inner code has no comments, since it
# is almost the same as the function load_significant_itemsets above
def load_transaction_database(transaction_database_file):
    f = open(transaction_database_file,'r')
    transactions = list()
    for line in f:
        transactions.append(set([int(i) for i in line.split()]))
    f.close()
    return transactions

# This function loads the class labels, returning a Numpy vector
def load_class_labels(labels_file):
    return np.loadtxt(labels_file)

# ------------------ FUNCTIONS TO COMPUTE PVALUES ---------

# Function to evaluate Fisher's exact test
def fisher_exact_test_pval(a,x,n,N):
    if int(x)==0 or int(x)==N:
        pval = 1
    else:
        prob_mass = hypergeom.pmf(np.arange(N+1),N,n,int(x))
        left_tail = prob_mass[:(a+1)].sum()
        right_tail = prob_mass[a:].sum()
        pval = np.minimum(left_tail,right_tail)
    return pval

# -----------------FUNCTIONS TO CHECK RESULTS ---------------

# This function checks if a given itemset is present in the set of significant
# itemsets (returns 1) or not (returns 0)
# The inputs are the itemset to be checked (encoded as a Python set) and a list
# with all significant itemsets (each itemset also encoded as a Python list)
def check_itemset(itemset,significant_itemsets):
    return itemset in significant_itemsets

# This function computes the support x of a given itemset within a transaction
# database, the cell count a and the corresponding p-value using Fisher's exact
# test. It takes as input the itemset, the transaction database and the class
# labels, all with the format returned by the FILE I/O functions above
def eval_itemset(itemset,transaction_database,class_labels):
    # This could also be done with a for loop, but functional programming is
    # cool! Some functions, such as map, filter or reduce come really handy
    # sometimes
    occurrences = np.array(map(lambda t: itemset.issubset(t),transaction_database))
    x = occurrences.sum()
    a = occurrences[class_labels==1].sum()
    pval = fisher_exact_test_pval(a,x,int(class_labels.sum()),len(class_labels))
    return (a,x,pval)

if __name__ == "__main__":
    path_to_transaction_database = '../datasets/transactions/mushroom.dat'
    path_to_class_labels = '../datasets/labels/mushroom.dat'
    transaction_database = load_transaction_database(path_to_transaction_database)
    class_labels = load_class_labels(path_to_class_labels)
    # Assume we have executed the C code to find significantly associated
    # patterns, generating a list of significant itemsets in a file itemsets.dat
    # and the corresponding pvalues in a file pvalues.dat
    (sig_itemsets,frequencies) = load_significant_itemsets('itemsets.dat')
    pvalues = load_significant_itemset_pvalues('pvalues.dat')
    # Check if all values reported are correct
    X = np.zeros(pvalues.shape)
    i = 0
    for itemset in sig_itemsets:
        X[i,:] = np.array(eval_itemset(itemset,transaction_database,class_labels))
        i = i + 1
    # X and pvalues should be identical (up to some numerical precision
    # differences)
    print np.abs((X-pvalues)).max(axis=0)
