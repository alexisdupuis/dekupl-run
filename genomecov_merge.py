import sys


def genomecov_merge(filename, resultsFilename, n):
    """
    Merge every n lines of the filename into a single line which is the average of the n lines, then write the result in "resultsFilename".

    :param filename: a file containing the coverage for each position on the chromosome Y
    :param resultsFilename: the name the user want for the output of this function
    :param n: the number of line to merge
    :type filename: str
    :type resultsFile: str
    :type n: int
    """
    nbLines = int(n)
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
            locusmax = int(linesList[i-2][1])
            locus = (locusmin + locusmax) // 2
            depthMean = 0
            for j in range(0, len(linesList)):
                depthMean += float(linesList[j][2])
            if (i != 2):
                depthMean = depthMean / len(linesList)
            resultsFile.write(linesList[0][0] + " " + str(locus) + " " + str(round(depthMean)) + '\n')
    resultsFile.close()


"""
def genomecov_merge_bedgraph(filename, resultsFilename, n):
    nbLines = int(n)
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
                    if (float(splitLine[3]) > 50):
                        continue
                    linesList += [splitLine]
                i += 1
            if (len(linesList) == 0):
                break
            locusmin = linesList[0][1]
            locusmax = linesList[i-2][2]
            depthMean = 0
            for j in range(0, len(linesList)):
                depthMean += float(linesList[j][3])
            if (i != 2):
                depthMean = depthMean / (i-2)
            resultsFile.write(linesList[0][0] + " " + locusmin + " " + locusmax + " " + str(round(depthMean)) + '\n')
    resultsFile.close()
"""


if __name__ == "__main__":
    genomecov_merge(sys.argv[1], sys.argv[2], sys.argv[3])