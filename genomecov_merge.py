import sys


def genomecov_merge(filename, resultsFilename, nbLines):
    """
    Merge every 'nbLines' lines of the filename into a single line which is the average of the 'nbLines' lines, then write the result in "resultsFilename".

    :param filename: a file containing the coverage for each position on the chromosome Y
    :param resultsFilename: the name the user wants for the output of this function
    :param n: the number of line to merge
    :type filename: str
    :type resultsFile: str
    :type n: int
    """
    resultsFile = open(resultsFilename, 'w')
    with open(filename, 'r') as genomecovFile:
        while (True):
            i = 1
            linesList = []
            while (i <= nbLines):
                genomecovLine = genomecovFile.readline()
                if (genomecovLine == ""):
                    break
                else:
                    splitLine = genomecovLine.split()
                    linesList += [splitLine]
                i += 1
            if (len(linesList) == 0):
                break
            locusmin = int(linesList[0][1])
            locusmax = int(linesList[-1][1])
            locus = (locusmin + locusmax) // 2
            depthMean = 0
            for j in range(0, len(linesList) - 1):
                depthMean += float(linesList[j][2])
            if (i != 2):
                depthMean = depthMean / len(linesList)
            resultsFile.write(linesList[0][0] + " " + str(locus) + " " + str(round(depthMean)) + '\n')
    resultsFile.close()


if __name__ == "__main__":
    genomecov_merge(sys.argv[1], sys.argv[2], int(sys.argv[3]))