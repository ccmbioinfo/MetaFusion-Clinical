#!/usr/bin/env python

import sys
import os
import argparse


class FusionResult:
    def __init__(self, tool_name, line, idxes):
        tmp = line.split()
        (idx_chr1, idx_pos1, idx_strand1, idx_chr2, idx_pos2, idx_strand2, idx_frame, idx_split_cnt, idx_pair_cnt,
         idx_gene1, idx_gene2, idx_gene_location1, idx_gene_location2) = idxes
        self.tool = tool_name
        # DRAGEN
        # MAF--RP11-679B19.2  0.687111    chr16:79619740:-    chr16:79246798:-    ENSG00000178573.6
        # ENSG00000261722.1   1   0   1
        # M70292:226:000000000-JG5LD:1:2101:16221:13209;M70292:226:000000000-JG5LD:1:1118:25958:17836;
        if self.tool == "DRAGEN":  # chr17:61906923:+
            # need to remove "chr" prefix 
            self.chr1 = tmp[idx_chr1].split(":")[0][3:]
            self.pos1 = tmp[idx_pos1].split(":")[1]
            self.strand1 = tmp[idx_strand1].split(":")[2]
            self.chr2 = tmp[idx_chr2].split(":")[0][3:]
            self.pos2 = tmp[idx_pos2].split(":")[1]
            self.strand2 = tmp[idx_strand2].split(":")[2]
            self.gene1 = tmp[idx_gene1].split("--")[0]
            self.gene2 = tmp[idx_gene2].split("--")[1]
            self.pair_cnt = tmp[idx_pair_cnt] if idx_pair_cnt != "NA" else -1
            # assign frame based on whether --examine_coding_effect flag was used
            self.frame = tmp[idx_frame] if idx_frame != "NA" else "NA"

        # STAR-FUSION NEEDS FIELDS ASSIGNED DIFFERENTLY
        elif self.tool == "STAR-Fusion":  # chr17:61906923:+
            # need to remove "chr" prefix 
            self.chr1 = tmp[idx_chr1].split(":")[0][3:]
            self.pos1 = tmp[idx_pos1].split(":")[1]
            self.strand1 = tmp[idx_strand1].split(":")[2]
            self.chr2 = tmp[idx_chr2].split(":")[0][3:]
            self.pos2 = tmp[idx_pos2].split(":")[1]
            self.strand2 = tmp[idx_strand2].split(":")[2]
            self.gene1 = tmp[idx_gene1].split("--")[0]
            self.gene2 = tmp[idx_gene2].split("--")[1]
            self.pair_cnt = tmp[idx_pair_cnt] if idx_pair_cnt != "NA" else -1
            # assign frame based on whether --examine_coding_effect flag was used
            self.frame = tmp[idx_frame] if idx_frame != "NA" else "NA"

        elif self.tool == "STAR-SEQR":  # chr10:124152694:-
            # need to remove "chr" prefix 
            self.chr1 = tmp[idx_chr1].split(":")[0][3:]
            self.pos1 = tmp[idx_pos1].split(":")[1]
            self.strand1 = tmp[idx_strand1].split(":")[2]
            self.chr2 = tmp[idx_chr2].split(":")[0][3:]
            self.pos2 = tmp[idx_pos2].split(":")[1]
            self.strand2 = tmp[idx_strand2].split(":")[2]
            self.gene1 = tmp[idx_gene1]
            self.gene2 = tmp[idx_gene2]
            # left and right paired reads separate, take maximum value
            self.pair_cnt = max(tmp[idx_pair_cnt], tmp[idx_pair_cnt + 1])
            self.frame = "NA"

        elif self.tool == "arriba":
            self.chr1 = tmp[idx_chr1].split(":")[0]
            self.pos1 = tmp[idx_pos1].split(":")[1]
            self.strand1 = tmp[idx_strand1].split("/")[0]
            self.chr2 = tmp[idx_chr2].split(":")[0]
            self.pos2 = tmp[idx_pos2].split(":")[1]
            self.strand2 = tmp[idx_strand2].split("/")[0]
            self.gene1 = tmp[idx_gene1]
            self.gene2 = tmp[idx_gene2]
            self.pair_cnt = tmp[idx_pair_cnt] if idx_pair_cnt != "NA" else -1
            self.frame = tmp[idx_frame]

        elif self.tool == "FusionMap":
            self.chr1 = tmp[idx_chr1]
            self.pos1 = tmp[idx_pos1]
            self.strand1 = tmp[idx_strand1] if idx_strand1 != "NA" else "NA"
            self.chr2 = tmp[idx_chr2]
            self.pos2 = tmp[idx_pos2]
            self.strand2 = tmp[idx_strand2] if idx_strand2 != "NA" else "NA"
            self.gene1 = tmp[idx_gene1]
            self.gene2 = tmp[idx_gene2]
            self.pair_cnt = tmp[idx_pair_cnt] if idx_pair_cnt != "NA" else -1
            self.frame = tmp[idx_frame]

        elif self.tool == "CICERO":
            tmp = line.split("\t")
            self.chr1 = tmp[idx_chr1][3:] if tmp[idx_chr1].startswith(("chr", "Chr")) else tmp[idx_chr1]
            self.pos1 = tmp[idx_pos1]
            self.strand1 = tmp[idx_strand1]
            self.chr2 = tmp[idx_chr2][3:] if tmp[idx_chr2].startswith(("chr", "Chr")) else tmp[idx_chr2]
            self.pos2 = tmp[idx_pos2]
            self.strand2 = tmp[idx_strand2]
            self.gene1 = tmp[idx_gene1]
            self.gene2 = tmp[idx_gene2]
            self.pair_cnt = -1
            self.split_cnt = tmp[idx_split_cnt[0]] + tmp[idx_split_cnt[1]]
            frame_note = {"0": "out-of-frame", "1": "in-frame",
                          "2": "canonical coding start site in tail",
                          "3": "possible 5' UTR fusion in tail"}
            frame = ",".join([frame_note[x] for x in set(tmp[idx_frame].split(",")) if x != ""])
            frame = "NA" if frame == "" else frame
            self.frame = frame

        else:  # Integrate uses this
            self.chr1 = tmp[idx_chr1]
            self.pos1 = tmp[idx_pos1]
            self.strand1 = tmp[idx_strand1] if idx_strand1 != "NA" else "NA"
            self.chr2 = tmp[idx_chr2]
            self.pos2 = tmp[idx_pos2]
            self.strand2 = tmp[idx_strand2] if idx_strand2 != "NA" else "NA"
            self.gene1 = tmp[idx_gene1]
            self.gene2 = tmp[idx_gene2]
            self.pair_cnt = tmp[idx_pair_cnt] if idx_pair_cnt != "NA" else -1
            self.frame = "NA"

        # FIELDS COMMON TO ALL FUSION CALLERS
        self.gene_location1 = tmp[idx_gene_location1] if idx_gene_location1 != "NA" else "NA"
        self.gene_location2 = tmp[idx_gene_location2] if idx_gene_location2 != "NA" else "NA"
        if not hasattr(self, "split_cnt"):  # If not already defined
            self.split_cnt = tmp[idx_split_cnt] if idx_split_cnt != "NA" else -1


class FusionResultFile:
    def __init__(self, result_file):
        self.fusion_results = []
        for line in open(result_file, "r"):
            tmp = line.split()
            # FusionGene	Score	LeftBreakpoint	RightBreakpoint	Gene1Id	Gene2Id	NumSplitReads
            # NumSoftClippedReads	NumPairedReads	ReadNames
            if tmp[0] == "#FusionGene" and tmp[1] == "Score":  # STAR-SEQR header
                self.tool = "DRAGEN"
                self._idx_chr1 = tmp.index("LeftBreakpoint")
                self._idx_chr2 = tmp.index("RightBreakpoint")
                self._idx_pos1 = tmp.index("LeftBreakpoint")
                self._idx_pos2 = tmp.index("RightBreakpoint")
                self._idx_strand1 = tmp.index("LeftBreakpoint")
                self._idx_strand2 = tmp.index("RightBreakpoint")
                # self._idx_split_cnt = tmp.index("NumSplitReads")
                self._idx_split_cnt = "NA"
                # self._idx_pair_cnt = tmp.index("NumPairedReads")
                self._idx_pair_cnt = "NA"
                self._idx_gene1 = tmp.index("#FusionGene")
                self._idx_gene2 = tmp.index("#FusionGene")
                self._idx_gene_location1 = "NA"
                self._idx_gene_location2 = "NA"
                self._idx_frame = "NA"

                # NAME	NREAD_SPANS	NREAD_JXNLEFT	NREAD_JXNRIGHT	FUSION_CLASS	SPLICE_TYPE
                # BRKPT_LEFT	BRKPT_RIGHT	LEFT_SYMBOL	RIGHT_SYMBOL	ANNOT_FORMAT	LEFT_ANNOT
                # RIGHT_ANNOT	DISTANCE	ASSEMBLED_CONTIGS	ASSEMBLY_CROSS_JXN	PRIMERS	ID
                # SPAN_CROSSHOM_SCORE	JXN_CROSSHOM_SCORE	OVERHANG_DIVERSITY	MINFRAG20	MINFRAG35
                # OVERHANG_MEANBQ	SPAN_MEANBQ	JXN_MEANBQ	OVERHANG_BQ15	SPAN_BQ15	JXN_BQ15
                # OVERHANG_MM	SPAN_MM	JXN_MM	OVERHANG_MEANLEN	SPAN_MEANLEN	JXN_MEANLEN	TPM_FUSION
                # TPM_LEFT	TPM_RIGHT	MAX_TRX_FUSION	DISPOSITION
            elif tmp[0] == "NAME" and tmp[1] == "NREAD_SPANS":  # STAR-SEQR header
                self.tool = "STAR-SEQR"
                self._idx_chr1 = tmp.index("BRKPT_LEFT")
                self._idx_chr2 = tmp.index("BRKPT_RIGHT")
                self._idx_pos1 = tmp.index("BRKPT_LEFT")
                self._idx_pos2 = tmp.index("BRKPT_RIGHT")
                self._idx_strand1 = tmp.index("BRKPT_LEFT")
                self._idx_strand2 = tmp.index("BRKPT_RIGHT")
                self._idx_split_cnt = tmp.index("NREAD_JXNLEFT")
                self._idx_pair_cnt = tmp.index("NREAD_SPANS")
                self._idx_gene1 = tmp.index("LEFT_SYMBOL")
                self._idx_gene2 = tmp.index("RIGHT_SYMBOL")
                self._idx_gene_location1 = "NA"
                self._idx_gene_location2 = "NA"
                self._idx_frame = "NA"

                # gene1  gene2   strand1(gene/fusion)    strand2(gene/fusion)    breakpoint1     breakpoint2
                # site1   site2 type     direction1      direction2      split_reads1    split_reads2
                # discordant_mates        coverage1     coverage2        confidence      closest_genomic_breakpoint1
                # closest_genomic_breakpoint2     filters fusion_transcript      reading_frame   peptide_sequence
                # read_identifiers
            elif tmp[0] == "#gene1":  # arriba header
                self.tool = "arriba"
                self._idx_chr1 = tmp.index("breakpoint1")
                self._idx_chr2 = tmp.index("breakpoint2")
                self._idx_pos1 = tmp.index("breakpoint1")
                self._idx_pos2 = tmp.index("breakpoint2")
                self._idx_strand1 = tmp.index("strand1(gene/fusion)")
                self._idx_strand2 = tmp.index("strand2(gene/fusion)")
                self._idx_split_cnt = tmp.index("split_reads1")
                self._idx_pair_cnt = tmp.index("discordant_mates")
                self._idx_gene1 = tmp.index("#gene1")
                self._idx_gene2 = tmp.index("gene2")
                self._idx_gene_location1 = "NA"
                self._idx_gene_location2 = "NA"
                self._idx_frame = tmp.index("reading_frame")

            elif tmp[0] == "#FusionName":  # STAR-Fusion header
                self.tool = "STAR-Fusion"
                self._idx_chr1 = tmp.index("LeftBreakpoint")
                self._idx_chr2 = tmp.index("RightBreakpoint")
                self._idx_pos1 = tmp.index("LeftBreakpoint")
                self._idx_pos2 = tmp.index("RightBreakpoint")
                self._idx_strand1 = tmp.index("LeftBreakpoint")
                self._idx_strand2 = tmp.index("RightBreakpoint")
                self._idx_split_cnt = tmp.index("JunctionReadCount")
                self._idx_pair_cnt = tmp.index("SpanningFragCount")
                self._idx_gene1 = tmp.index("#FusionName")
                self._idx_gene2 = tmp.index("#FusionName")
                self._idx_gene_location1 = "NA"
                self._idx_gene_location2 = "NA"
                try:
                    self._idx_frame = tmp.index("PROT_FUSION_TYPE")
                except ValueError:
                    self._idx_frame = "NA"

            elif tmp[0] == "cluster_id":  # defuse header
                self.tool = "Defuse"

                self._idx_chr1 = tmp.index("gene_chromosome1")
                self._idx_chr2 = tmp.index("gene_chromosome2")
                self._idx_pos1 = tmp.index("genomic_break_pos1")
                self._idx_pos2 = tmp.index("genomic_break_pos2")
                self._idx_strand1 = tmp.index("genomic_strand1")
                self._idx_strand2 = tmp.index("genomic_strand2")
                self._idx_split_cnt = tmp.index("splitr_count")
                self._idx_pair_cnt = tmp.index("span_count")
                # replacing gene names with ENSG numbers
                self._idx_gene1 = tmp.index("gene_name1")
                self._idx_gene2 = tmp.index("gene_name2")
                self._idx_gene_location1 = tmp.index("gene_location1")
                self._idx_gene_location2 = tmp.index("gene_location2")
                self._idx_frame = "NA"

            elif tmp[0] == "FusionID":  # fusionmap header
                self.tool = "FusionMap"

                self._idx_chr1 = tmp.index("Chromosome1")
                self._idx_chr2 = tmp.index("Chromosome2")
                self._idx_pos1 = tmp.index("Position1")
                self._idx_pos2 = tmp.index("Position2")
                self._idx_strand1 = "NA"
                self._idx_strand2 = "NA"
                for i in range(len(tmp)):
                    name = tmp[i]
                    if "UniqueCuttingPositionCount" in name:
                        self._idx_split_cnt = i
                        break
                self._idx_pair_cnt = "NA"
                self._idx_gene1 = tmp.index("KnownGene1")
                self._idx_gene2 = tmp.index("KnownGene2")
                # Frame index
                self._idx_frame = tmp.index("FrameShiftClass")
                self._idx_gene_location1 = "NA"
                self._idx_gene_location2 = "NA"

            elif tmp[0] == "GeneName1":  # ericscript header
                self.tool = "EricScript"

                self._idx_chr1 = tmp.index("chr1")
                self._idx_chr2 = tmp.index("chr2")
                self._idx_pos1 = tmp.index("Breakpoint1")
                self._idx_pos2 = tmp.index("Breakpoint2")
                self._idx_strand1 = "NA"
                self._idx_strand2 = "NA"
                self._idx_split_cnt = tmp.index("crossingreads")
                self._idx_pair_cnt = tmp.index("spanningreads")
                self._idx_gene1 = tmp.index("GeneName1")
                self._idx_gene2 = tmp.index("GeneName2")
                self._idx_gene_location1 = "NA"
                self._idx_gene_location2 = "NA"
                self._idx_frame = "NA"

            elif tmp[0] == "#5P":  # integrate header
                self.tool = "Integrate"

                self._idx_chr1 = tmp.index("Chr1")
                self._idx_chr2 = tmp.index("Chr2")
                self._idx_pos1 = tmp.index("RNA_BK1")
                self._idx_pos2 = tmp.index("RNA_BK2")
                self._idx_strand1 = "NA"
                self._idx_strand2 = "NA"
                self._idx_split_cnt = tmp.index("NUM_SP_RNA")
                self._idx_pair_cnt = tmp.index("NUM_EN_RNA")
                self._idx_gene1 = tmp.index("#5P")
                self._idx_gene2 = tmp.index("3P")
                self._idx_gene_location1 = "NA"
                self._idx_gene_location2 = "NA"
                self._idx_frame = "NA"

            elif tmp[0] == "sample":
                tmp = line.split("\t")
                self.tool = "CICERO"

                self._idx_chr1 = tmp.index("chrA")
                self._idx_chr2 = tmp.index("chrB")
                self._idx_pos1 = tmp.index("posA")
                self._idx_pos2 = tmp.index("posB")
                self._idx_strand1 = tmp.index("ortA")
                self._idx_strand2 = tmp.index("ortB")
                self._idx_split_cnt = (tmp.index("readsA"), tmp.index("readsB"))
                self._idx_pair_cnt = "NA"
                self._idx_gene1 = tmp.index("geneA")
                self._idx_gene2 = tmp.index("geneB")
                self._idx_gene_location1 = tmp.index("featureA")
                self._idx_gene_location2 = tmp.index("featureB")
                self._idx_frame = tmp.index("frame")

            else:
                fusion_line = FusionResult(self.tool, line,
                                           [self._idx_chr1, self._idx_pos1, self._idx_strand1, self._idx_chr2,
                                            self._idx_pos2, self._idx_strand2, self._idx_frame, self._idx_split_cnt,
                                            self._idx_pair_cnt, self._idx_gene1, self._idx_gene2,
                                            self._idx_gene_location1, self._idx_gene_location2])
                if fusion_line.pos1.isdigit() and fusion_line.pos2.isdigit():
                    self.fusion_results.append(fusion_line)


class SampleInfo:
    def __init__(self, line, idxes):
        tmp = line.split()
        # idx_sample, idx_disease, idx_lib, idx_sample_type = idxes
        idx_sample, idx_disease, idx_sample_type = idxes
        self.sample = tmp[idx_sample]
        self.disease = tmp[idx_disease]
        # self.lib = tmp[idx_lib]
        self.sample_type = tmp[idx_sample_type]


class SampleInfoFile:
    def __init__(self, sampleinfo_file):
        self.sampleinfo_dict = {}
        for line in open(sampleinfo_file, "r"):
            tmp = line.split()
            if tmp[0] == "sample":  # header
                self._idx_sample = tmp.index("sample")
                self._idx_disease = tmp.index("disease")
                self._idx_sample_type = tmp.index("sample_type")
            else:
                idxes = [self._idx_sample, self._idx_disease, self._idx_sample_type]
                self.sampleinfo_dict.setdefault(args.sample, SampleInfo(line, idxes))


parser = argparse.ArgumentParser()
parser.add_argument("--sample", help="Sample name")
parser.add_argument('--sample_info_file', action='store', help='Sample infomation file')
parser.add_argument("--tool", help="Tool name")
parser.add_argument("--fusion_result_file", help="Original fusion result file generated by tools")
parser.add_argument("--outdir", help="Output folder")
args = parser.parse_args()

sample = args.sample
sampleinfo = SampleInfoFile(args.sample_info_file)
disease_name = sampleinfo.sampleinfo_dict[sample].disease
# re-assign sample type 
sample_type = sampleinfo.sampleinfo_dict[sample].sample_type
if sample_type in ("NT", "N", "Normal"):
    sample_type = "Normal"
elif sample_type in ("TP", "T", "Tumor"):
    sample_type = "Tumor"
else:
    print >> sys.stderr, "Unknown sample type:", sample_type
    sys.exit(1)

tool = args.tool
fusion_result_file = args.fusion_result_file
fusion_results = FusionResultFile(args.fusion_result_file).fusion_results

print(sample, args.sample_info_file, tool, disease_name, sample_type, args.fusion_result_file, args.outdir)

cff_path = os.path.join(args.outdir, args.sample + "." + args.tool + ".cff")
out_file = open(cff_path, "w")
for fusion in fusion_results:
    print >> out_file, "\t".join(map(str, [fusion.chr1, fusion.pos1, fusion.strand1, fusion.chr2, fusion.pos2,
                                           fusion.strand2, fusion.frame, sample, sample_type, disease_name, args.tool,
                                           fusion.split_cnt, fusion.pair_cnt,
                                           fusion.gene1, fusion.gene_location1, fusion.gene2, fusion.gene_location2]))

out_file.close()
