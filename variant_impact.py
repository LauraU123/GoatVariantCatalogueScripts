import pysam

# chromosomes in T2Tgoat1.0 ref
chrs_dir = {"1": 166528736, "2": 147555794, "3": 132669774, "4": 126865584, "5": 127915634, "6": 126033241, "7": 114753091, "8": 117817743, "9": 100793790, "10": 109114799, "11": 113240263, "12" : 98214097, "13": 90012720, "14": 100201825, "15": 87982261, "16": 90863288, "17": 80799350, "18": 71942568, "19": 71331175, "20": 80554990, "21": 80325541, "22" : 70628230, "23" : 55774162, "24" : 70861467, "25": 51558643, "26": 55986131, "27": 50558960, "28": 50812434, "29" : 56632579, "X": 144688805, "Y" : 20906000}


#calculating total bases affected by large INS in VCF files
def get_ins(vcf_path):
    """Sums up affected bases by INS in SV vcf files"""
    total_nr =0
    total_sum = 0
    with pysam.VariantFile(vcf_path) as vcf:
        for rec in vcf:
            if rec.info.get("SVTYPE") == "INS": 
                start = rec.pos
                total_nr+=1
                svlen = rec.info.get("SVLEN")
                if type(svlen) == tuple:
                    total_sum+=svlen[0]
                else:
                    total_sum+=svlen
    print(total_nr)
    print(total_sum)
print("Sawfish ins")
get_ins("../../UPLOAD_VARIANT/2026april_sawfish2.vcf.gz")
print("Sniffles ins")
get_ins("../../UPLOAD_VARIANT/2026april_sniffles2.vcf.gz")


def parse_vcf(vcf_path, chrom, svtype):
    """Parses VCF by chromosome and SVtype, returning intervals affected across it"""
    intervals = []
    with pysam.VariantFile(vcf_path) as vcf:
        for rec in vcf.fetch(chrom):
            if rec.info.get("SVTYPE") == svtype: 
                start = rec.pos
                end = rec.stop
                intervals.append((start, end))
    return intervals

def merge_intervals(intervals):
    """Takes as input intervals from parse_vcf, and merges them if they overlap"""
    intervals.sort()
    merged = []
    for start, end in intervals:
        if not merged or merged[-1][1] < start: merged.append([start, end])
        else: merged[-1][1] = max(merged[-1][1], end)
    return merged

# finding nr of affected bases
def affected_bases(intervals):
    """Sums up affected bases by calculating """
    merged = merge_intervals(intervals)
    return sum(end - start + 1 for start, end in merged)


# calculating fraction total affected for Sawfish2
for type_ in ["DEL", "DUP", "INV"]:
    total_length = 0
    total_affected = 0
    for chr_, length in chrs_dir.items():
        vars_ = parse_vcf("../../UPLOAD_VARIANT/2026april_sawfish2.vcf.gz", chr_, type_)
        interval = merge_intervals(vars_)
        affected = affected_bases(interval)
        total_length+=length
        total_affected+=affected
    print("Sawfish2")
    print(type_, total_length, total_affected, total_affected/total_length)


# calculating fraction total for Sniffles2

for type_ in ["DEL", "DUP", "INV"]:
    total_length = 0
    total_affected = 0
    for chr_, length in chrs_dir.items():
        vars_ = parse_vcf("../../UPLOAD_VARIANT/2026april_sniffles2.vcf.gz", chr_, type_)
        interval = merge_intervals(vars_)
        affected = affected_bases(interval)
        total_length+=length
        total_affected+=affected
    print("Sniffles2")
    print(type_, total_length, total_affected, total_affected/total_length)

### FOR SMALL VARIANT PROPORTION FINDING 
def snv_counts(vcf_path):
    """Finding counts of INS and DEL in small variant VCF files"""
    sums = {"INS": 0, "DEL": 0}
    with pysam.VariantFile(vcf_path) as vcf:
        for rec in vcf:
            REF = rec.ref
            ref_len = len(REF)
            for alt in rec.alts:
                diff = len(alt)-ref_len
                if diff > 0:
                    sums["INS"] +=diff
                elif diff < 0: 
                    sums["DEL"] += -(diff)
    print(sums)
