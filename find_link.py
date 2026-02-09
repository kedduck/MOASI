import sys
from collections import defaultdict

class UnionFind:
    def __init__(self):
        self.parent = {}
    
    def find(self, x):
        if x not in self.parent:
            self.parent[x] = x
        if self.parent[x] != x:
            self.parent[x] = self.find(self.parent[x])
        return self.parent[x]
    
    def union(self, x, y):
        rx, ry = self.find(x), self.find(y)
        if rx != ry:
            self.parent[ry] = rx

def interval_distance(s_start, s_end, p_start, p_end):
    if s_end < p_start:
        return p_start - s_end
    elif p_end < s_start:
        return s_start - p_end
    else:
        return 0

def main(bed_file, max_dist=87000):
    style_by_chrom = defaultdict(list)
    pollen_by_chrom = defaultdict(list)
    gene_to_line = {}

    with open(bed_file, 'r') as f:
        for line in f:
            stripped = line.strip()
            if not stripped:
                continue
            parts = stripped.split('\t')
            if len(parts) < 5:
                continue
            chrom, start, end, gene_id, gene_type = parts[0], int(parts[1]), int(parts[2]), parts[3], parts[4]
            if gene_type == "Style":
                style_by_chrom[chrom].append((start, end, gene_id, stripped))
            elif gene_type == "Pollen":
                pollen_by_chrom[chrom].append((start, end, gene_id, stripped))
            gene_to_line[gene_id] = stripped

    all_style_genes = set(gid for chrom in style_by_chrom.values() for _, _, gid, _ in chrom)
    all_pollen_genes = set(gid for chrom in pollen_by_chrom.values() for _, _, gid, _ in chrom)

    first_cluster = True

    for chrom in sorted(set(style_by_chrom.keys()) | set(pollen_by_chrom.keys())):
        styles = style_by_chrom[chrom]
        pollens = pollen_by_chrom[chrom]
        if not styles or not pollens:
            continue

        uf = UnionFind()
        all_genes = set()

        for _, _, gid, _ in styles:
            all_genes.add(gid)
        for _, _, gid, _ in pollens:
            all_genes.add(gid)

        for s_start, s_end, s_id, _ in styles:
            for p_start, p_end, p_id, _ in pollens:
                dist = interval_distance(s_start, s_end, p_start, p_end)
                if dist <= max_dist:
                    uf.union(s_id, p_id)

        components = defaultdict(list)
        for gene in all_genes:
            root = uf.find(gene)
            components[root].append(gene)

        valid_clusters = []
        style_gene_set = all_style_genes
        pollen_gene_set = all_pollen_genes

        for comp_genes in components.values():
            if len(comp_genes) < 2:
                continue
            has_style = any(g in style_gene_set for g in comp_genes)
            has_pollen = any(g in pollen_gene_set for g in comp_genes)
            if not (has_style and has_pollen):
                continue

            lines = [gene_to_line[gid] for gid in comp_genes]
            lines.sort(key=lambda x: int(x.split('\t')[1]))  # 按 start 排序行
            min_start = int(lines[0].split('\t')[1])
            valid_clusters.append((min_start, lines))

        valid_clusters.sort(key=lambda x: x[0])

        for min_start, lines in valid_clusters:
            if not first_cluster:
                print("###")
            for ln in lines:
                print(ln)
            first_cluster = False

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python cluster_style_pollen_pairs.py <bed_file>", file=sys.stderr)
        sys.exit(1)
    main(sys.argv[1])
