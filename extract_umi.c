/*
 * extract_umi.c
 * A C-code to extract UMI sequences from paired end sequencing data
 * Compilation command
 * gcc -O2 -o extract_umi extract_umi.c -lz
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

#define UMI_LENGTH (9)
const char *usage = "extract_umi in.fastq.R1.gz in.fastq.R2.gz out1.fastq out2.fastq out.umi.fq";
char UMI_R1[UMI_LENGTH] = {0,};
char UMI_R2[UMI_LENGTH] = {0,};
char UMI_QUAL1[UMI_LENGTH] = {0,};
char UMI_QUAL2[UMI_LENGTH] = {0,};
size_t UMI_SZ = sizeof(char) * (UMI_LENGTH - 1);
size_t UMI_SZ2 = sizeof(char) * UMI_LENGTH;

int main(int argc, char **argv)
{
	if(argc < 6){
		fprintf(stderr, "[usage] %s \n\n", usage);
		exit(-1);
	}
	gzFile fp1 = NULL;
	gzFile fp2 = NULL;
	FILE *ofp = NULL;
	FILE *O1 = NULL;
	FILE *O2 = NULL;
	kseq_t *kseq1 = NULL;
	kseq_t *kseq2 = NULL;
	char *seq1 = NULL;
	char *seq2 = NULL;
	char *qual1 = NULL;
	char *qual2 = NULL;
	int l1 = 0;
	int l2 = 0;
	const char *fqname1 = argv[1];
	const char *fqname2 = argv[2];
	const char *O1name = argv[3];
	const char *O2name = argv[4];
	const char *oname = argv[5];
	fp1 = gzopen(fqname1, "r");
	if(NULL == fp1){
		fprintf(stderr, "[extract] failed to open %s\n", fqname1);
		exit(-1);
	}
	fp2 = gzopen(fqname2, "r");
	if(NULL == fp2){
		fprintf(stderr, "[extract] failed to open %s\n", fqname2);
		gzclose(fp1);
		exit(-1);
	}
	ofp = fopen(oname, "w");
	if(NULL == ofp) {
		fprintf(stderr, "[extract] failed to open%s\n", oname);
		gzclose(fp1); gzclose(fp2);
		exit(-1);
	}

	O1 = fopen(O1name, "w");
	if(NULL == O1){
		fprintf(stderr, "[extract] failed to open %s\n", O1name);
		gzclose(fp1); gzclose(fp2); fclose(ofp);
		exit(-1);
	}
	O2 = fopen(O2name, "w");
	if(NULL == O2){
		fprintf(stderr, "[extract] failed to open %s\n", O2name);
		gzclose(fp1); gzclose(fp2); fclose(ofp);fclose(O1);
		exit(-1);

	}
	kseq1 = kseq_init(fp1);
	if(NULL == kseq1){
		fprintf(stderr, "[extract] kseq_init fail - read1\n");
		gzclose(fp1); gzclose(fp2); fclose(ofp); fclose(O1); fclose(O2);
		exit(-1);
	}
	kseq2 = kseq_init(fp2);
	if(NULL == kseq2){
		fprintf(stderr, "[extract] kseq_inti fail - read2\n");
		gzclose(fp1); gzclose(fp2); fclose(ofp); fclose(O1); fclose(O2);
		kseq_destroy(kseq1);
		exit(-1);
	}

	while( ((l1 = kseq_read(kseq1)) >= 0) &&
		((l2 = kseq_read(kseq2)) >= 0)) {
		memset(UMI_R1, 0, UMI_SZ2);
		memset(UMI_R2, 0, UMI_SZ2);
		memset(UMI_QUAL1, 0, UMI_SZ2);
		memset(UMI_QUAL2, 0, UMI_SZ2);
		seq1 = kseq1->seq.s;
		seq2 = kseq2->seq.s;
		qual1 = kseq1->qual.s;
		qual2 = kseq2->qual.s;
		memcpy(UMI_R1, seq1, UMI_SZ);
		memcpy(UMI_R2, seq2, UMI_SZ);
		memcpy(UMI_QUAL1, qual1, UMI_SZ);
		memcpy(UMI_QUAL2, qual2, UMI_SZ);
		fprintf(ofp ,"@%s\n%s-%s\n+\n%s!%s\n", kseq1->name.s, UMI_R1, UMI_R2, UMI_QUAL1, UMI_QUAL2);
		fprintf(O1, "@%s\n%s\n+\n%s\n", kseq1->name.s , seq1 + UMI_LENGTH, qual1 + UMI_LENGTH);
		fprintf(O2, "@%s\n%s\n+\n%s\n", kseq2->name.s, seq2 + UMI_LENGTH, qual2 + UMI_LENGTH);
	}
	gzclose(fp1); gzclose(fp2); fclose(ofp); fclose(O1); fclose(O2);
	kseq_destroy(kseq1); kseq_destroy(kseq2);
	return 0;
}
