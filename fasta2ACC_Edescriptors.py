"""Usage: fasta2ACC -a <antigen fasta filename> -na <non-antigen fasta filename> -o <output filename>
-l <lag variable> -n <number of records from the fastafile>
If you want all the records from the fasta file to be in output file set -n <all>"""
import sys, os, csv
from Bio import SeqIO
import getopt

d = {'A': (0.008, 0.134, -0.475, -0.039, 0.181),
     'R': (0.171, -0.361, 0.107, -0.258, -0.364),
     'N': (0.255, 0.038, 0.117, 0.118, -0.055),
     'D': (0.303, -0.057, -0.014, 0.225, 0.156),
     'C': (-0.132, 0.174, 0.07, 0.565, -0.374),
     'Q': (0.149, -0.184, 0.03, 0.035, -0.112),
     'E': (0.221, -0.28, -0.315, 0.157, 0.303),
     'G': (0.218, 0.562, -0.024, 0.018, 0.106),
     'H': (0.023, -0.177, 0.041, 0.28, -0.021),
     'I': (-0.353, 0.071, -0.088, -0.195, -0.107),
     'L': (-0.267, 0.018, -0.265, -0.274, 0.206),
     'K': (0.243, -0.339, -0.044, -0.325, -0.027),
     'M': (-0.239, -0.141, -0.155, 0.321, 0.077),
     'F': (-0.329, -0.023, 0.072, -0.002, 0.208),
     'P': (0.173, 0.286, 0.407, -0.215, 0.384),
     'S': (0.199, 0.238, -0.015, -0.068, -0.196),
     'T': (0.068, 0.147, -0.015, -0.132, -0.274),
     'W': (-0.296, -0.186, 0.389, 0.083, 0.297),
     'Y': (-0.141, -0.057, 0.425, -0.096, -0.091),
     'V': (-0.274, 0.136, -0.187, -0.196, -0.299)}


def makedata(fname):
    handle = open(fname)
    data = {}
    for line in handle:
        lst = line[:-1].split('\t')
        data[lst[0]] = [float(x) for x in lst[1:]]
    handle.close()
    return data


def chck(sequence):
    if len(sequence) < 5:
        return False
    letters = {}
    for letter in sequence:
        if letter in ['X', 'U', 'B', 'J', 'O', 'Z']:
            return False
        if letter not in letters.keys():
            letters[letter] = 1
        else:
            letters[letter] += 1

    if len(letters.keys()) <= 4 and letters.keys() == ['A', 'C', 'T', 'G']:
        return False
    else:
        return True


def Edescriptor_import(seq, descriptors):
    # Returns the E descriptor values for every aminoacid in a sequence
    # as a dictionary with keyword the index of the AA
    z = {}
    lst = list(seq)
    i = 0
    for aa in lst:
        z[i] = descriptors[aa]
        i += 1
    return z


def acc_calculate(z, lag):
    # Calculates ACC of E descriptors of a sequence for a given lag value
    # and returns a dict with ACC
    acc = {}
    dnumber = len(z[list(z)[0]])
    for k in range(dnumber):
        for j in range(dnumber):
            for l in range(1, lag+1):
                f = 0.0
                for i in range(len(z.keys())-l):

                    f = f+(z[i][k]*z[i+l][j]/(len(z.keys())-l))

                acc['ACC'+str(k+1)+str(j+1)+str(l)] = f

    return acc


def makeGInumber(record_id):
    g = record_id.split('|')

    if len(g) > 1:
        ginum = g[1]
    else:
        ginum = g[0]
    return ginum


def readFastaFile(fname, recnum):
    # Reads a fale in fasta format. Returns list with records.
    # Before writing checks if the sequence is not DNA and
    # for unknown characters for aminoacids.
    if recnum == 'all':
        recnum = 0
    data = []
    sequence = []
    handle = open(fname, "r")
    h = open("log.log", "w")
    i = 0
    for record in SeqIO.parse(handle, "fasta"):
        if not chck(str(record.seq)):
            ginum = makeGInumber(record.id)
            print("Problem with sequence %s.\nSequence %s was skipped." % (ginum, ginum))
            r = "Record with ginumber "+ str(ginum)+" is either not a protein or there is an unidentified character in its sequence"
            h.write(r)
            continue
        if str(record.seq) not in sequence:
            sequence.append(str(record.seq))
            data.append(record)
            i += 1
            if i == int(recnum):
                break
        else:

            ginum = makeGInumber(record.id)
            print("There is a record with sequence equal to the sequence with gi number %s.\nRecord with gi number %s would be skipped.\n" %(ginum,ginum))
            continue
    h.close()
    handle.close()
    return data


def write2csvFile(fname, records1, records0, descriptors, lag):
    # Writes the ACC values to Tab separated file with legend.
    with open(fname, "w") as csvfile:
        i = 0
        writer = csv.writer(csvfile, delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL)

        for record in records1:
            ginum = makeGInumber(record.id)
            sequence = str(record.seq)
            if len(sequence) < lag:
                continue
            zvalues = Edescriptor_import(sequence, descriptors)
            acc_values = acc_calculate(zvalues, lag)
            acc = [str(acc_values[key]) for key in sorted(acc_values.keys())]
            if not i:
                legend = ['protein_id']+sorted(acc_values.keys())+['class']
                writer.writerow(legend)
                i += 1
            row = [str(ginum)] + acc+['1']
            writer.writerow(row)
            i += 1
        for record in records0:
            ginum = makeGInumber(record.id)
            sequence = str(record.seq)
            if len(sequence) < lag:
                continue
            zvalues = Edescriptor_import(sequence, descriptors)
            acc_values = acc_calculate(zvalues, lag)
            acc = [str(acc_values[key]) for key in sorted(acc_values.keys())]
            if not i:
                legend = ['protein_id']+sorted(acc_values.keys())+['class']
                writer.writerow(legend)
                i += 1
            row = [str(ginum)] + acc+['0']
            writer.writerow(row)
            i += 1

    return True


def process(fasta1, fasta0, output, lag, nrec):
    data1 = readFastaFile(cwd + '/' + fasta1, nrec)
    data0 = readFastaFile(cwd + '/' + fasta0, nrec)
    write2csvFile(cwd + '/' + output, data1, data0, d, int(lag))
    return True


def main(argv=None):
    if argv is None:
        argv = sys.argv

    try:
        opts, args = getopt.getopt(argv[1:], 'h:a:u:o:l:n:', ["help", "iafile=", "inafile=", "ofile=", "lag=", "nrec="])
    except getopt.error as msg:
        print(msg)
        print("for help use --help")
        sys.exit(2)
    for o, a in opts:
        if o in ("-h", "--help"):
            print(__doc__)
            sys.exit(0)
        elif o in ("-a", "--iafile"):
            fasta1 = a
        elif o in ("-u", "--inafile"):
            fasta2 = a
        elif o in ("-o", "--ofile"):
            output = a
        elif o in ("-l", "--lag"):
            lag = a
        elif o in ("-n", "--nrec"):
            nrec = a
        else:
            print(o)

    process(fasta1, fasta2, output, lag, nrec)


if __name__ == '__main__':
    cwd = os.getcwd()
    print("Current working directory: {0}".format(cwd))

    if len(sys.argv) == 11:
        sys.exit(main())
    else:
        print("Sorry. Too few arguments", len(sys.argv))
        print(__doc__)
        sys.exit(0)
