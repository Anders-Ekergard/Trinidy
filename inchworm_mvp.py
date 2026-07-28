reads = ["ATGGC", "TGGCA", "GGCAT"]
k = 3
kmers: list[str] =[]
for read in reads:
    for i in range(len(read)-k +1):
        kmers.append(read[i:i+k])
def kmers_in_contig(contig: str, k: int) -> set[str]:
    result: set[str] = set()
    for i in range(len(contig) - k + 1):
        kmer = contig[i:i+k]
        result.add(kmer)
    return result
def most_common_kmers(kmers: list[str])-> dict[str, int]: # type: ignore
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
if __name__ == "__main__":
    kmers_list = most_common_kmers(kmers)
    print(f"Most common kmers: {kmers_list}")
    results = expand(kmers_list, k)
    print(f"Results: {results}")