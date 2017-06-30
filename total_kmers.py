import sys
import matplotlib.pyplot as plt
from subprocess import call


def move_fastq_files(phenotype):
    i = 0
    phenotypeFile = open(phenotype, 'r')
    phenotypeFile.readline()  # skip the first line
    for line in phenotypeFile:
        splitLine = line.split()
        sampleName = splitLine[0] + "_chrY.fastq.gz"
        call(["mv", "./data/" + sampleName, "../chrY_compressed-fastq/"])
        i += 1
        if (i >= 190):
            break


def add_samples_to_config_file(n, phenotype):
    """
    Add the n last samples described in the file 'phenotype' into the config file of de-kupl.

    :param n: the number of sample to add
    :param phenotype: a file containing informations about the different
    samples
    :type n: int
    :type phenotype: str
    """
    i = 0
    phenotypeFile = open(phenotype, 'r')
    configFile = open("config.json", 'a')
    phenotypeFile.readline()  # skip the first line
    for line in phenotypeFile:
        splitLine = line.split()
        sampleName = splitLine[0] + "_chrY"
        if (int(splitLine[4]) == 1):
            condition = 'A'
        else:
            condition = 'B'
        configFile.write(", {\n\t\t\t\"name\" : \"" + sampleName + "\",\n\t")
        configFile.write("\t\t\"condition\" : \"" + condition + "\"\n\t\t}")
        i += 1
        if (i >= int(n)):
            break
    phenotypeFile.close()
    configFile.write("\n\t]\n}")
    configFile.close()


def total_kmers(n, countsFile):
    """
    Count the total of k-mers for one sample in the tabulated file 'countsFile' and writes it into a file called 'total_kmers.txt'.

    The first column in this file is the tags for each k-mer, and the other
    columns are the occurrences of each k-mers for each samples.

    :param n: the n-th sample in the file
    :param countsFile: the name of the k-mers file
    :type n: int
    :type countsFile: str
    """
    i = 0
    countsFile = open(countsFile, 'r')
    line = countsFile.readline()  # skip the first line
    splitLine = line.split('\t')
    sampleName = splitLine[n]
    print("Counting the total of k-mers for " + sampleName)
    for line in countsFile:
        splitLine = line.split('\t')
        i += int(splitLine[n])
    countsFile.close()
    resultsFile = open("total_kmers.txt", 'a')
    resultsFile.write(sampleName + '\t' + str(i) + '\n')
    resultsFile.close()


def total_kmers_for_one_sample(sampleName, countsFile):
    """
    Count the total of k-mers for one sample in the tabulated file 'countsFile' and writes it into a file called 'total_kmers.txt'.

    The first column in this file is the tags for each k-mer, and the other
    columns are the occurrences of each k-mers for each samples.

    :param sampleName: the sample in the file
    :param countsFile: the name of the k-mers file
    :type sampleName: str
    :type countsFile: str
    """
    i, j = 0, 0
    countsFile = open(countsFile, 'r')
    line = countsFile.readline()
    splitLine = line.split('\t')
    while (j < len(splitLine) - 1):
        if (splitLine[j] == sampleName):
            break
        j += 1
    for line in countsFile:
        splitLine = line.split('\t')
        i += int(splitLine[j])
    countsFile.close()
    resultsFile = open("total_kmers.txt", 'a')
    resultsFile.write(sampleName + '\t' + str(i) + '\n')
    resultsFile.close()


def total_kmers_histogram(totalKmerFile):
    """
    Build a histogram with the file containing the total of k-mers for each samples.

    :param totalKmerFile: the name of the file
    :type totalKmerFile: str
    """
    totalKmerList = []
    with open(totalKmerFile, 'r') as tkf:
        for line in tkf:
            splitLine = line.split()
            if (splitLine[0] == "Average"):
                break
            totalKmerList += [int(splitLine[2])]
    plt.hist(totalKmerList, bins=50, color='green', edgecolor='blue')
    plt.axis([200000000, 1200000000, 0, 30])
    plt.xlabel('Nombre total de k-mers')
    plt.ylabel("Nombre d'individus")
    plt.title("Répartition des individus selon le nombre total de k-mers comptés.")
    plt.savefig('total_kmers.png')


def total_kmers_average():
    """Give the average of the total k-mers of all the samples in the file 'total_kmers.txt'."""
    totals_list = []
    nbLines = 0
    totalKmerFile = open("total_kmers.txt", 'r')
    resultsFile = open("total_kmers.txt", 'a')
    for line in totalKmerFile:
        splitLine = line.split()
        totals_list += [int(splitLine[1])]
        nbLines += 1
    totalKmerFile.close()
    average = sum(totals_list) // nbLines
    resultsFile.write("Average\t" + str(average))
    resultsFile.close()


def get_condition(phenotype, sample):
    cpt = 0
    phenotypeFile = open(phenotype, 'r')
    for line in phenotypeFile:
        splitLine = line.split()
        if (splitLine[0] == sample[:-5]):
            print(int(splitLine[4]))
            cpt = 1
            if (cpt != 1):
                print(3)


def genomecov_mean(mergedGenomecovFile, finalFileName, condition):
    finalGenomecovFile = open(finalFileName, 'a')
    with open(mergedGenomecovFile) as mgf:
        for line in mgf:
            splitLine = line.split('\t')
            genomecovSum = 0
            for i in range(2, len(splitLine)):
                genomecovSum += int(splitLine[i])
            genomecovMean = genomecovSum // (len(splitLine) -2)
            finalGenomecovFile.write(splitLine[0] + condition + '\t' + splitLine[1] + '\t' + str(genomecovMean) + '\n')
    finalGenomecovFile.close()

def add_conditions(phenotype, textFile):
    """
    Add the condition of each sample into the file containing the k-mer totals.

    :param phenotype: the file containing the conditions
    :param textFile: the file with the k-mer totals
    :type phenotype: str
    :type textFile: str
    """
    phenotypeFile = open(phenotype, 'r')
    totalKmersFile = open(textFile, 'r')
    resultsFile = open(textFile + "_with_conditions.txt", 'w')
    for line in totalKmersFile:
        totalSplitLine = line.split()
        phenotypeLine = phenotypeFile.readline()
        phenotypeSplitLine = phenotypeLine.split()
        if (totalSplitLine[0] == "Average"):
            resultsFile.write(line)
            lastLine = totalKmersFile.readline()
            resultsFile.write(lastLine)
            break
        while (totalSplitLine[0][:-5] != phenotypeSplitLine[0]):
            oLine = phenotypeFile.readline()
            phenotypeSplitLine = oLine.split()
        phenotypeFile.seek(0, 0)
        resultsFile.write(line[:-1] + ' ' + phenotypeSplitLine[4] + '\n')
    phenotypeFile.close()
    totalKmersFile.close()
    resultsFile.close()


def add_normalization_factor(sampleConditionFile):
    """
    Create a new file from the 'sampleConditionFile' and add the normalization factor for each samples.

    :param sampleConditionFile: the file with the sample conditions
    :type sampleConditionFile: str
    """
    totalKmersList = []
    nfList = []
    i = 0
    sampleConditionFile = open(sampleConditionFile, 'r')
    sampleConditionFile.readline()
    nfFile = open("sample_conditions_full.tsv", 'w')
    nfFile.write("sample\tcondition\tnormalization_factor\n")
    totalKmersFile = open("total_kmers.txt", 'r')
    for line in totalKmersFile:
        splitLine = line.split()
        if (splitLine[0] == "Average"):
            average = int(splitLine[2])
            break
        totalKmersList += [int(splitLine[2])]
    for elt in totalKmersList:
        nfList += [average / elt]
    for line in sampleConditionFile:
        nfFile.write(line[:-1] + '\t' + str(nfList[i]) + '\n')
        i += 1
    sampleConditionFile.close()
    nfFile.close()
    totalKmersFile.close()


def select_contigs(threshold, pvalue, contigFile):
    """
    Select the contigs in the file "contigFile" according to the two other parameters : threshold is the minimum length of the contig, pvalue the maximum pvalue associated.

    :param threshold: the minimum length
    :param pvalue: the maximum pvalue
    :param contigFile: the file containing the contigFile
    :type threshold: int
    :type pvalue: float
    :type contigFile: str
    """
    i = 1
    resultsFile = open("contigs.fasta", 'w')
    with open(contigFile, 'r') as ctgFile:
        ctgFile.readline()  # skip the header
        for line in ctgFile:
            splitLine = line.split('\t')
            if (len(splitLine[1]) > int(threshold) and (int(splitLine[3] < pvalue))):
                resultsFile.write(">contig " + str(i) + " length=" + str(len(splitLine[1])) + " pvalue=" + splitLine[3] + " nb_kmers=" + str(splitLine[0]) + '\n')
                resultsFile.write(splitLine[1] + '\n')
                i += 1
    resultsFile.close()


def cut_contig(contig):
    """
    Cut the contig in k-mers of 31 nucleotids.

    :param contig: the contig to cut
    :type contig: str
    """
    kmerList = []
    for i in range(0, len(contig)-30):
        kmerList += [contig[i:i+31]]
    return kmerList


def revcomp(kmer):
    """
    Give the reverse complement of the kmer.

    :param kmer: the kmer to reverse complement
    :type kmer: str
    """
    length = len(kmer) - 1
    revcompResult = list(kmer)
    for i in range(len(kmer)):
        revcompResult[length - i] = revcomp_letter(kmer[i])
    return "".join(revcompResult)


def revcomp_letter(nucleotid):
    """
    Give the complement of the nucleotid.

    :param nucleotid: a nucleotid
    :type letter: str
    """
    if (nucleotid == 'A'):
        return 'T'
    elif (nucleotid == 'T'):
        return 'A'
    elif (nucleotid == 'C'):
        return 'G'
    elif (nucleotid == 'G'):
        return 'C'
    else:
        return 'N'


def kmer_counts(kmer, countsFile, filename, cpt):
    revcompResult = revcomp(kmer)
    resultsFile = open(filename, 'a')
    with open(countsFile, 'r') as cFile:
        headerLine = cFile.readline()
        if (cpt == 0):
            resultsFile.write("orientation\t" + headerLine)
        for line in cFile:
            splitLine = line.split('\t')
            if (kmer == splitLine[0]):
                resultsFile.write("forward\t"+ line)
                # break
            if (revcompResult == splitLine[0]):
                resultsFile.write("reverse\t" + line)
                # break
            else:
                continue


def ratios(contig, phenotype):
    zHealthy, zSick, sHealthy, sSick = 0, 0, 0, 0
    i = 6 # skip the first columns as they are not used in this function
    with open(contig, 'r') as ctFile:
        headerLine = ctFile.readline()
        splitHeader = headerLine.split('\t')
        countsLine = ctFile.readline()
        splitCounts = countsLine.split('\t')
        while (i < len(splitCounts)):
            sample = [splitHeader[i], splitCounts[i]]
            with open(phenotype, 'r') as phenotypeFile:
                for line in phenotypeFile:
                    splitPhenotype = line.split()
                    if (splitPhenotype[0] == sample[0][:-5]):
                        break
                if (int(splitPhenotype[4]) == 1):
                    if (float(sample[1]) == 0):
                        zHealthy += 1
                    else:
                        sHealthy += 1
                else:
                    if (float(sample[1]) == 0):
                        zSick += 1
                    else:
                        sSick += 1
            i += 1
        print([zHealthy, sHealthy, zSick, sSick])
        print([round(((zHealthy/len(splitCounts)) * 100), 2), round(((sHealthy/len(splitCounts)) * 100), 2), round(((zSick/len(splitCounts)) * 100), 2), round(((sSick/len(splitCounts)) * 100), 2)])


def contig_or_not(contig):
    resultsList = []
    dr = int(contig)
    with open("cTetA.tsv", 'r') as contigFile:
        headerLine = contigFile.readline()
        splitHeader = headerLine.split('\t')
        countsLine = contigFile.readline()
        splitCounts = countsLine.split('\t')
        for i in range(6, len(splitHeader) - 1):
            if(dr == 0 and round(float(splitCounts[i])) == 0):
                resultsList += [splitHeader[i]]
            if(dr == 1 and round(float(splitCounts[i])) != 0):
                resultsList += [splitHeader[i]]
    print(resultsList)


def below_nb_kmers(nb_kmers):
    nbKmers = int(nb_kmers)
    with open("total_kmers_with_conditions.txt", 'r') as totalKmersFile:
        for line in totalKmersFile:
            splitLine = line.split()
            if (int(splitLine[2]) < nbKmers and splitLine[3] == "1"):
                print(splitLine[3])


def main():
    if (sys.argv[1] == "total-kmers"):
        for i in range(246, 340):
            total_kmers(i, sys.argv[2])
    elif (sys.argv[1] == "total-kmers-average"):
        total_kmers_average()
    elif (sys.argv[1] == "add-samples-to-config-file"):
        add_samples_to_config_file(sys.argv[2], sys.argv[3])
    elif (sys.argv[1] == "move-fastq-files"):
        move_fastq_files(sys.argv[2])
    # significantly_low_total_kmers(sys.argv[1])
    elif (sys.argv[1] == "add-conditions"):
        add_conditions(sys.argv[2], sys.argv[3])
    elif(sys.argv[1] == "add-normalization-factor"):
        add_normalization_factor(sys.argv[2])
    # total_kmers_histogram(sys.argv[1])
    elif (sys.argv[1] == "select-contigs"):
        select_contigs(sys.argv[2], sys.argv[3], sys.argv[4])
    elif (sys.argv[1] == "total-kmers-for-one-sample"):
        total_kmers_for_one_sample(sys.argv[2], sys.argv[3])
    elif (sys.argv[1] == "cut-contig"):
        kmerList = cut_contig(sys.argv[2])
        cpt = 0
        for i in range(len(kmerList)):
            print(kmerList[i])
            kmer_counts(kmerList[i], sys.argv[3], sys.argv[4], cpt)
            cpt += 1
    elif (sys.argv[1] == "kmer-counts"):
        kmer_counts(sys.argv[2], sys.argv[3])
    elif (sys.argv[1] == "ratios"):
        ratios(sys.argv[2], sys.argv[3])
    elif (sys.argv[1] == "below-nb-kmers"):
        below_nb_kmers(sys.argv[2])
    elif (sys.argv[1] == "get-condition"):
        get_condition(sys.argv[2], sys.argv[3])
    elif (sys.argv[1] == "genomecov-mean"):
        genomecov_mean(sys.argv[2], sys.argv[3], sys.argv[4])
    elif (sys.argv[1] == "contig-or-not"):
        contig_or_not(sys.argv[2])
    else:
        print("This function doesn't exist.")


if __name__ == "__main__":
    main()
