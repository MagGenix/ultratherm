from Bio import SeqIO

class blacklist():
    def __init__(self, path: str):
        self.blacklist_sequences = []
        try:
            with open(path) as handle:
                for record in SeqIO.parse(handle, "fasta"):
                    self.blacklist_sequences.append(record.seq)
            self.is_empty = len(self.blacklist_sequences) == 0
            self.blacklist_path = path
        except IOError:
            Warning("Invalid path to blacklist or invalid blacklist file. Blacklist will be empty")
            self.is_empty = True
            self.blacklist_path = ""
