__author__ = 'snow'
from intervaltree import Interval, IntervalTree


class GeneConversion:
    """
    Convert position to refseq, refseq to gene symbol.
    """
    def __init__(self):
        gene_file = "/home/snow/pro/data/ref_gene"
        refseq_symbol = "/home/snow/pro/data/refseq2symbol.txt"
        self.map_refseq2symbol = {}
        self.map_symbol2refseq = {}
        self.intervals_dict = {}
        self.read_to_interval_tree(gene_file)
        self.read_refseq_symbol(refseq_symbol)

    def read_refseq_symbol(self, refseq_symbol):
        f = open(refseq_symbol)
        for l in f:
            lsp = l.strip().split()
            self.map_refseq2symbol[lsp[0]] = lsp[1]
            self.map_symbol2refseq[lsp[1]] = lsp[0]

    def read_to_interval_tree(self, gene_file):
        genes = open(gene_file)
        genes.readline()
        for l in genes:
            refseq, chromosome, strand, start_s, end_s = l.strip().split(',')
            start, end = int(start_s), int(end_s)
            if chromosome not in self.intervals_dict.keys():
                self.intervals_dict[chromosome] = IntervalTree()
            self.intervals_dict[chromosome][start:end] = refseq

    def refseq2symbol(self):
        pass

    def symbol2refseq(self):
        pass

    def pos2refseq(self, chr, start, end):
        results = []
        if chr not in self.intervals_dict.keys():
            print ("Chromosome not in dictionary --- probably error!")
            return None
        ivs = self.intervals_dict[chr][start:end]
        if len(ivs) == 0:
            return None
        for iv in ivs:
            results.append(iv.data)
        return results

    def pos2symbol(self, chr, start, end):
        results = []
        refseqs = self.pos2refseq(chr, start, end)
        for r in refseqs:
            s = self.map_refseq2symbol[r]
            if s not in results:
                results.append(s)
        return results


def main():
    a = GeneConversion()
    print (a.pos2refseq('11',89057522, 89057523))
    print (a.pos2symbol('11',89057522, 89057523))

if __name__ == "__main__":
    main()