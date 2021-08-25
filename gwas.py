import argparse
import csv
import fastq2matrix as fm
from collections import defaultdict
from tqdm import tqdm 
from uuid import uuid4


class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'


class vcf_class(fm.vcf.vcf_class):
	def is_annotated(self):
		check = False
		for l in fm.cmd_out(f"bcftools view -h {self.filename}"):
			if "##INFO=<ID=BCSQ" in l:
				check = True
		return check
	def get_genesum(self,outfile=None):
		self.outfile = outfile
		if self.is_annotated()==False:
			quit(bcolors.FAIL + "\nError: VCF file has not been annotated with bcftools csq, please do this before running this script\n" + bcolors.ENDC)
		if self.outfile==None:
			self.outfile = self.prefix+".genesum"
		genesum = defaultdict(lambda:defaultdict(int))
		O = open(self.outfile,"w")
		for l in tqdm(fm.cmd_out("bcftools query -f '[%%SAMPLE\\t%%GT\\t%%TBCSQ\\n]' %(filename)s" % vars(self))):
			row = l.split()
			#por4A	1/1	synonymous|Rv0002|gene1|protein_coding|+|109L|2378G>A	synonymous|Rv0002|gene1|protein_coding|+|109L|2378G>A
			info = row[2].split("|")
			if info[0]=="synonymous": continue
			if info[0][0]=="@": continue
			genesum[info[1]][row[0]]+=1
		for gene in genesum:
			O.write("%s\tNA\tNA\t%s\n" % (gene,"\t".join(str(genesum[gene][s]) for s in self.samples)))
		O.close()
	def get_mean_genotype(self,outfile=None):
		self.outfile = outfile
		if self.outfile==None:
			self.outfile = self.prefix+".geno"
		O = open(self.outfile,"w")
		for l in tqdm(fm.cmd_out("bcftools query -f '%%CHROM\\t%%POS\\t%%REF\\t%%ALT[\\t%%TGT]\\n' %(filename)s" % vars(self))):
			row = l.rstrip().split()
			alts = row[3].split(",")
			for alt in alts:
				ref = "%s/%s" % (row[2],row[2])
				tmp = "%s/%s" % (alt,alt)
				genos = []
				for x in row[4:]:
					if x==ref:
						genos.append("0")
					elif x==tmp:
						genos.append("1")
					else:
						genos.append("NA")
				O.write("%s, %s, %s, %s\n" % (row[0]+"_"+row[1]+"_"+alt,row[2],alt,", ".join(genos)))
		O.close()


def main(args):
	vcf = vcf_class(args.vcf)
	vcf.get_mean_genotype()
	vcf.get_genesum()
	geno_file = vcf.prefix+".geno"
	genesum_file = vcf.prefix+".genesum"

	meta = {}
	for row in csv.DictReader(open(args.pheno_csv)):
		if args.pheno_id_column not in row:
			quit(bcolors.FAIL + f"\nError: Can't find {args.pheno_id_column} in csv columns, set it with --pheno-id-column\n" + bcolors.ENDC)
		meta[row[args.pheno_id_column]] = row

	if args.pheno_name_file:
		phenos = [l.strip() for l in open(args.pheno_name_file)]
	else:
		phenos = args.pheno_names

	
	cmd_file = str(uuid4())
	X = open(cmd_file,"w")
	for pheno in phenos:
		pheno_file = f"{pheno}.pheno"
		if pheno not in row:
			quit(bcolors.FAIL +  f"{pheno} not in CSV file" + bcolors.ENDC)
		P = open(pheno_file,"w")
		P.write("\n".join([meta[s][pheno] if pheno in meta[s] else "NA" for s in vcf.samples]))
		P.close()
		X.write("gemma -p %s -g %s -gk 1 -o %s -maf 0.00005 -miss 0.99 && gemma  -lmm 1 -p %s -g %s  -k output/%s.cXX.txt  -o %s -maf 0.00005 -miss 0.99 && gemma  -lmm 1 -p %s -g %s  -k output/%s.cXX.txt  -o %s.genesum -notsnp\n" % (pheno_file,geno_file,pheno,pheno_file,geno_file,pheno,pheno,pheno_file,genesum_file,pheno,pheno))
	X.close()

	if args.preprocess:
		fm.log("Preprocessing finished\n", True)
	else:
		fm.run_cmd("cat %s | parallel -j %s" % (cmd_file,args.threads))

parser = argparse.ArgumentParser(description='GWAS pipeline',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--vcf', help='BCF file',required = True)
parser.add_argument('--pheno-csv', help='CSV file',required = True)

group = parser.add_mutually_exclusive_group(required = True)
group.add_argument('--pheno-names', nargs="+",help='Columns file')
group.add_argument('--pheno-name-file', help='Columns file')

parser.add_argument('--pheno-id-column',default="id", help='Columns file')
parser.add_argument('--preprocess',default=False,action='store_true', help='Columns file')
parser.add_argument('--threads','-t',default=4,type=int, help='Columns file')
parser.set_defaults(func=main)

args = parser.parse_args()
args.func(args)

