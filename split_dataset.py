"""Usage: split_dataset.py -a <path to immunogens fasta filename> -u <path to non-immunogens fasta filename> """


import getopt
from Bio import SeqIO
import random, os, sys


def makeTestandTrain(fname, interval):
    counter = 0
    records = {}
    for i in interval:
        records[i] = []
    handle = open(fname, "r")
    for record in SeqIO.parse(handle, "fasta"):
        counter += 1
        seq = str(record.seq)
        for l in sorted(records.keys()):
            if len(seq) < l:
                records[l].append(record)
                break

    handle.close()
    for key in sorted(records.keys()):
        test = []
        train = []
        for i in range(len(records[key])):
            if i and i % 5:
                train.append(records[key][i])
            else:
                test.append(records[key][i])

    return (test, train)


def process(fasta1, fasta0):
    f1 = fasta1
    f2 = fasta0
    intrvl = [20, 50, 100, 150, 200, 250, 300, 350, 400, 500, 600, 750, 1000, 2750]
    result = makeTestandTrain(f1, intrvl)
    with open(cwd + "/testset_immunogens.fasta", "w") as output_handle:
        SeqIO.write(result[0], output_handle, "fasta")
    with open(cwd + "/trainset_immunogens.fasta", "w") as output_handle:
        SeqIO.write(result[1], output_handle, "fasta")
    intrvl = [45, 70, 100, 150, 200, 250, 300, 350, 400, 500, 600, 750, 1000, 2750]
    result = makeTestandTrain(f2, intrvl)
    random.shuffle(result[0])
    result[1].append(result[0].pop())
    with open(cwd + "/testset_nonimmunogens.fasta", "w") as output_handle:
        SeqIO.write(result[0], output_handle, "fasta")
    with open(cwd + "/trainset_nonimmunogens.fasta", "w") as output_handle:
        SeqIO.write(result[1], output_handle, "fasta")
    return True


def main(argv=None):
    if argv is None:
        argv = sys.argv

    try:
        opts, args = getopt.getopt(argv[1:], 'h:a:u:l:n:', ["help", "iafile=", "inafile="])
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
        else:
            print(o)

    process(fasta1, fasta2)


if __name__ == '__main__':
    cwd = os.getcwd()

    if len(sys.argv) == 5:
        sys.exit(main())
    else:
        print("Sorry. Too few arguments", len(sys.argv))
        print(__doc__)
        sys.exit(0)
