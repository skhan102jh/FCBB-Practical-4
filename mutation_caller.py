import pysam
import sys 
import argparse
import math 

# METHODS TO CALCULATE PROBABILITIES
def p_base_given_base(observed_base, assumed_base):
    e = 0.1
    if (observed_base != assumed_base):
        return e / 3
    else:
        return 1 - e

def p_base_given_allele(observed_base, allele):
    return 0.5 * p_base_given_base(observed_base, allele[0]) + 0.5 * p_base_given_base(observed_base, allele[1])

# SET UP PARSING FOR COMMAND LINE
parser = argparse.ArgumentParser(description='Bringing in normal and cancer bam')
parser.add_argument('-n','--normal',required=True, help='Normal BAM file')
parser.add_argument('-c','--cancer',required=True, help='Cancer BAM file')
args = parser.parse_args()

# Open files of interest 
normal_bam = pysam.AlignmentFile(args.normal, "rb")
cancer_bam = pysam.AlignmentFile(args.cancer, "rb")

# Calculate max number of columns in the .bam files 
normal_columns = 0
for pileupcolumn in normal_bam.pileup():
    normal_columns += 1
cancer_columns = 0
for pileupcolumn in cancer_bam.pileup():
    cancer_columns += 1
num_positions = max(normal_columns, cancer_columns)
valid_position = [True for i in range(num_positions)]

# Indicate invalid positions in cancer.bam (coverage < 20)
for pileupcolumn in cancer_bam.pileup():
    if (pileupcolumn.n < 20):
        valid_position[pileupcolumn.pos] = False


# Create a list of "truth" alleles    
normal_alleles = ["empty" for i in range(num_positions)]
allele_list =['AA','CC','GG','TT','AC','AG','AT','CG','CT','GT']
# FILL THE DICTIONARY FOR NORMAL ALLELES
for pileupcolumn in normal_bam.pileup():
    # IF COLUMN IS INVALID, INDICATE INSUFFICIENT COVERAGE
    if (pileupcolumn.n < 20 or not valid_position[pileupcolumn.pos]):
        valid_position[pileupcolumn.pos] = False
        sys.stdout.write("Insufficient coverage at position " + str(pileupcolumn.pos) + "\n")
    # IF COLUMN IS VALID, DETERMINE ALLELE THAT MAXIMIZES LIKELIHOOD
    else: 
        # Create a dictionary that holds the likelihood of the observed data 
        # given a specific allele
        max_likelihood = 0
        max_allele = ""
        for allele in allele_list:
            prob_data_given_allele = 1
            for pileupread in pileupcolumn.pileups:
                observed_base = pileupread.alignment.query_sequence[pileupread.query_position]
                prob_data_given_allele *= p_base_given_allele(observed_base, allele)
            if (prob_data_given_allele > max_likelihood):
                max_likelihood = prob_data_given_allele
                max_allele = allele
        # Return the allele that yields maximum likelihood of observed data
        if (math.log(max_likelihood) < -50):
            sys.stdout.write("Position " + str(pileupcolumn.pos) + " has ambiguous genotype\n")
            valid_position[pileupcolumn.pos] = False
        else: 
            normal_alleles[pileupcolumn.pos] = max_allele
        # Find the allele that yields the highest likelihood of observing data

for pileupcolumn in cancer_bam.pileup():
    # IF COLUMN IS INVALID, INDICATE INSUFFICIENT COVERAGE
    if (pileupcolumn.n < 20 or not valid_position[pileupcolumn.pos]):
        continue 
    # IF COLUMN IS VALID, SEE IF THERE IS A CANDIDATE MUTATION 
    else: 
        likelihood = 1 
        normal_allele = normal_alleles[pileupcolumn.pos]
        for pileupread in pileupcolumn.pileups:
            observed_base = pileupread.alignment.query_sequence[pileupread.query_position]
            likelihood *= p_base_given_allele(observed_base, normal_allele)
        # Determine if the probability indicates a mutation 
        if (math.log(likelihood) < -75):
            sys.stdout.write("Position " + str(pileupcolumn.pos) 
                + " has a candidate somatic mutation (Log-likelihood="
                + str(math.log(likelihood)) + ")\n")
