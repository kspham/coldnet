__author__ = 'snow'


def simon_parsing(csv_file):
    """
    parsing SIMON csv file from 10.000 ASD
    :param csv_file:
    :return: dictionary with key = patient and value = patient mutation
    """
    csv = open(csv_file, 'r')
    patient_mutations = {}
    csv.readline()
    for l in csv:
        lsp = l.split(',')
        if len(lsp) > 7:
            lsp = l.split('"')
            chr, start, end, ref = lsp[0].strip(',').split(',')
            var = lsp[1]
            tech, patient = lsp[2].lstrip(',').split(',')
        else:
            chr, start, end, ref, var, tech, patient = l.split(',')
        if patient not in patient_mutations.keys():
            patient_mutations[patient] = []
        patient_mutations[patient].append((chr, start, end, ref, var))
    return patient_mutations


def main():
    d = simon_parsing("/home/snow/autism/somatic.csv")
    for k in d.keys():
        print k
        print len(d[k])


if __name__ == "__main__":
    main()
