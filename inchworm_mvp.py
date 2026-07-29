"""
This module contains functions for assembling contigs from reads using k-mers.
By Anders Ekergård, 2026
For more information, see the README.md file.
"""


def find_kmers(reads: list[str], k: int) -> list[str]:
    """
    Find all k-mers in a list of reads.
    args:
        reads: list[str] - list of reads
        k: int - the length of each k-mer
    returns:
        list[str] - a list of all k-mers in the reads
    """
    kmers: list[str] = []
    for read in reads:
        for i in range(len(read) - k + 1):
            kmer = read[i:i+k]
            kmers.append(kmer)
    return kmers

def kmers_in_contig(contig: str, k: int) -> set[str]:
    """
    Find all k-mers in a contig.
    args:
        contig: str - the contig sequence
        k: int - the length of each k-mer
    returns:
        set[str] - a set of all k-mers in the contig
    """
    result: set[str] = set()
    for i in range(len(contig) - k + 1):
        kmer = contig[i:i+k]
        result.add(kmer)
    return result


def most_common_kmers(kmers: list[str])-> dict[str, int]:
    """    
    Find the most common kemers from a list of kmers
    args:
        kmers: list[str] - list of kmers
    
    returns:
        dict[str, int] - dictionary of kmers and their counts
    """

    most_common_kmers ={}
    for kmer in kmers:
        if kmer in most_common_kmers:
            most_common_kmers[kmer] += 1
        else:
            most_common_kmers[kmer] = 1
    return most_common_kmers


def expand_forward(kmers_count: dict [str, int], contig: str, k: int)-> str:
    """

    Expand in the contig by finding the next kmer that overlaps with the current kmer.
    Each k-mer is used at most once so the assembly cannot get stuck in a cycle.
    args:
        kmers_count: dict[str, int] - dictionary of kmers and their counts
        contig: str - the current contig sequence
        k: int - the length of each k-mer
    returns:
        str - the expanded contig sequence
    """
    if not kmers_count:
        return ""

    used_kmers: set[str] = set()

    
    while True:
        candidates: list[str] = []
        for kmer in kmers_count:
            if kmer in used_kmers:
                continue
            if kmer[:-1] == contig[-(k-1):]:
                candidates.append(kmer)
                used_kmers.add(kmer)
        if not candidates:
            break

        best = max(candidates, key=kmers_count.get)
        
        contig += best[-1]
    return contig
def expand_backward(kmers_count: dict [str, int], contig: str, k: int)->str:
    """
    Expand in the contig by finding the next kmer that overlaps with the current kmer.
    Each k-mer is used at most once so the assembly cannot get stuck in a cycle.
    args:
        kmers_count: dict[str, int] - dictionary of kmers and their counts
        contig: str - the current contig sequence
        k: int - the length of each k-mer
    returns:
        str - the expanded contig sequence
    """
    used_kmers: set[str] = set()
    while True:
        candidates: list[str] = []
        for kmer in kmers_count:
            if kmer in used_kmers:
                continue
            if kmer[1:] == contig[:k-1]:
                candidates.append(kmer)
                used_kmers.add(kmer)
        if not candidates:
            break       
        best = max(candidates, key=kmers_count.get)
        contig =  best[0]+contig
    return contig
def expand(kmers_count: dict[str, int], k: int) -> list[str]:
    """
    Expand the contig by finding the next kmer that overlaps with the current kmer.
    args:
        kmers_count: dict[str, int] - dictionary of kmers and their counts
        k: int - the length of each k-mer
    returns:
        list[str] - a list of expanded contig sequences
    """
    used_kmers: set[str] = set()
    contigs: list[str] = []

    available = {kmer: c for kmer, c in kmers_count.items() if kmer not in used_kmers}
    while available:
        start = max(available, key=available.get)
        contig = expand_backward(available, start, k)
        contig = expand_forward(available, contig, k) 
        contigs.append(contig)
        used_kmers |= kmers_in_contig(contig, k)
        available = {kmer: c for kmer, c in kmers_count.items() if kmer not in used_kmers}

    return contigs
def load_reads(filepath: str) -> list[str]:
    """
    Load reads from a FASTQ file.
    """
    reads: list[str] = []
    with open(filepath) as f:
        for index, seq in enumerate(f):
            if index % 4 == 1:
                reads.append(seq.strip())
    return reads

def save_contigs(contigs: list[str], filepath: str) -> None:
    """
    Save contigs to a FASTA file.
    """
    with open(filepath, "w") as f: 
        for index, seq in enumerate(contigs):
            f.write(f">contig_{index}\n{seq}\n")
    return None
if __name__ == "__main__":
    k = 3  # Example k-mer length
    kmers_list = find_kmers(load_reads("reads.fastq"), k)

    results = expand(most_common_kmers(kmers_list), k)
    print(f"Results: {results}")
    save_contigs(results, "contigs.fasta")
