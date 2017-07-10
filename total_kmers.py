import sys
from math import sqrt


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
    """
    Get the condition associated to the sample.

    :param phentoype: the name of a file containing informations about the condition of each samples of a population
    :param sample: a sample
    :type phenotype: str
    :type sample: str
    """
    cpt = 0
    phenotypeFile = open(phenotype, 'r')
    for line in phenotypeFile:
        splitLine = line.split()
        if (splitLine[0] == sample[:-5]):
            print(int(splitLine[4]))
            cpt = 1
            if (cpt != 1):
                print(3)


def genomecov_mean_and_sd(mergedGenomecovFile, finalFileName, condition):
    """
    Compute the mean and the standard deviation for each line of a file containing the coverage per position for several samples.

    :param mergedGenomecovFile: the name of a file containing the coverage for each samples
    :param finalFileName: the name of the results file
    :param condition: a suffix used to distinguish two different conditions
    :type mergedGenomecovFile: str
    :type finalFileName: str
    :type condition: str
    """
    finalGenomecovFile = open(finalFileName, 'a')
    with open(mergedGenomecovFile) as mgf:
        for line in mgf:
            splitLine = line.split('\t')
            genomecovSum = 0
            individualCoverage = []
            for i in range(2, len(splitLine)):
                genomecovSum += int(splitLine[i])
                individualCoverage += [int(splitLine[i])]
            genomecovMean = genomecovSum // (len(splitLine) -2)
            genomecovVariance = sum([(x - genomecovMean)**2 for x in individualCoverage]) // (len(splitLine) - 2)
            genomecovSD = sqrt(genomecovVariance)
            finalGenomecovFile.write(splitLine[0] + condition + '\t' + splitLine[1] + '\t' + str(genomecovMean) + '\t' + str(genomecovSD) +'\n')
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
            average = int(splitLine[1])
            break
        totalKmersList += [int(splitLine[1])]
    for elt in totalKmersList:
        nfList += [average / elt]
    for line in sampleConditionFile:
        nfFile.write(line[:-1] + '\t' + str(nfList[i]) + '\n')
        i += 1
    sampleConditionFile.close()
    nfFile.close()
    totalKmersFile.close()


def below_nb_kmers(nb_kmers):
    """
    Select all the samples which have a total of k-mers below the threshold.
    
    :param nb_kmers: the threshold
    :type nb_kmers: int
    """
    with open("total_kmers_with_conditions.txt", 'r') as totalKmersFile:
        for line in totalKmersFile:
            splitLine = line.split()
            if (int(splitLine[2]) < nb_kmers and splitLine[3] == "1"):
                print(splitLine[3])


def main():
    """Allow to use this file like an executable."""
    if (sys.argv[1] == "total-kmers"):
        for i in range(0, 340):
            total_kmers(i, sys.argv[2])
    elif (sys.argv[1] == "total-kmers-average"):
        total_kmers_average()
    elif (sys.argv[1] == "add-samples-to-config-file"):
        add_samples_to_config_file(sys.argv[2], sys.argv[3])
    elif (sys.argv[1] == "add-conditions"):
        add_conditions(sys.argv[2], sys.argv[3])
    elif(sys.argv[1] == "add-normalization-factor"):
        add_normalization_factor(sys.argv[2])
    elif (sys.argv[1] == "get-condition"):
        get_condition(sys.argv[2], sys.argv[3])
    elif (sys.argv[1] == "genomecov-mean-and-sd"):
        genomecov_mean_and_sd(sys.argv[2], sys.argv[3], sys.argv[4])
    elif (sys.argv[1] == "below-nb-kmers"):
        below_nb_kmers(int(sys.argv[2]))
    else:
        print("This function doesn't exist.")


if __name__ == "__main__":
    main()