import sys

def select_contigs(lengthThreshold, pValueThreshold, contigFile):
    """
    Select contigs from a file and write them into another file.

    The contigs will be selected according to a minimum length and a maximum pvalue.

    :param lengthThreshold: the minimum length
    :param pValueThreshold: the maximum pvalue
    :param contigFile: the file containing the contigFile
    :type lengthThreshold: int
    :type pValueThreshold: float
    :type contigFile: str
    """
    i = 1
    resultsFile = open("contigs.fasta", 'w')
    with open(contigFile, 'r') as ctgFile:
        ctgFile.readline()  # skip the header
        for line in ctgFile:
            splitLine = line.split('\t')
            if (len(splitLine[1]) > int(lengthThreshold) and (int(splitLine[3] < pValueThreshold))):
                resultsFile.write(">contig " + str(i) + " length=" + str(len(splitLine[1])) + " pvalue=" + splitLine[3] + " nb_kmers=" + str(splitLine[0]) + '\n')
                resultsFile.write(splitLine[1] + '\n')
                i += 1
    resultsFile.close()


def cut_contig(contig):
    """
    Cut the contig in k-mers of 31 nucleotids.

    :param contig: the contig to cut
    :type contig: str
    :return: a list of k-mers
    :rtype: list
    """
    kmerList = []
    for i in range(0, len(contig)-30):
        kmerList += [contig[i:i+31]]
    return kmerList


def revcomp(kmer):
    """
    Give the reverse complement of the kmer.

    :param kmer: a k-mer
    :type kmer: str
    :return: the reverse complement
    :rtype: str
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
    :return: the complement
    :rtype: str
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
    """
    Select the lines within 'countsFile' corresponding to the k-mer and write it in a new file called 'filename'.

    The cpt parameters allows to write the header only once in the results file.

    :param kmer: a k-mer
    :param countsFile: a file containing informations about the k-mer on a certain line
    :param filename: the name given to the results file
    :param cpt: a counter
    :type kmer: str
    :type countsFile: str
    :type filename: str
    :type cpt: int
    """
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
                break
            if (revcompResult == splitLine[0]):
                resultsFile.write("reverse\t" + line)
                break
            else:
                continue


def ratios(contig, phenotype):
    """
    Compute ratios corresponding to the presence/absence of differential repeatition of the contig within the population associated with the condition of the samples.

    :param contig: the contig
    :param phenotype: the name of a file containing informations about the condition of each samples of the population
    :type contig: str
    :type phenotype: str
    """
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


def contig_or_not(diffRepeat, contigFile):
    """
    Check for each samples if the contig is differentially repeated or not.

    The file contains two lines of merged-diff-counts.tsv, the final output of DE-kupl. The first one is the header,
    the second contains the normalized counts for each samples.
    'boolean' allows to select the samples with a differential repetition (boolean=1) or the other ones (boolean=0).
    :param diffRepeat: selection criteria
    :param contigFile: the name of the file
    :type boolean: int
    :type contigFile: str
    """
    resultsList = []
    dr = int(diffRepeat)
    with open(contigFile, 'r') as ctgFile:
        headerLine = ctgFile.readline()
        splitHeader = headerLine.split('\t')
        countsLine = ctgFile.readline()
        splitCounts = countsLine.split('\t')
        for i in range(6, len(splitHeader) - 1):
            if(dr == 0 and round(float(splitCounts[i])) == 0):
                resultsList += [splitHeader[i]]
            if(dr == 1 and round(float(splitCounts[i])) != 0):
                resultsList += [splitHeader[i]]
    print(resultsList)


def main():
    """Allow to use this file like an executable."""
    if (sys.argv[1] == "select-contigs"):
        select_contigs(sys.argv[2], sys.argv[3], sys.argv[4])
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
    elif (sys.argv[1] == "contig-or-not"):
        contig_or_not(sys.argv[2])
    else:
        print("This function doesn't exist.")


if __name__ == "__main__":
    main()