from Bio import SeqIO
from Bio import Seq

class blacklist():
    def __init__(self, path: str):
        self.blacklist_sequences = []
        try:
            with open(path) as handle:
                for record in SeqIO.parse(handle, "fasta"):
                    self.blacklist_sequences.append(record.seq)
            self.is_empty = len(self.blacklist_sequences) == 0
        except IOError:
            Warning("Invalid path to blacklist or invalid blacklist file. Blacklist will be empty")
            self.is_empty = True
    def is_blacklisted(self, nucl: Seq):
        if self.is_empty:
            return False
        if nucl.sequence in self.blacklist_sequences:
            return True
        return False