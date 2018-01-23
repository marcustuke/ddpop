import itertools
import numpy as np
import sys 
#import pandas as pd
# /gpfs/mrc0/projects/Research_Project-MRC158833/UKBiobank/500K_Genetic_data/direct_genotypes/ukb907_cal_v2.fam /gpfs/mrc0/projects/Research_Project-MRC158833/UKBiobank/500K_Genetic_data/direct_genotypes/ukb_snp_chrX_v2.bim /gpfs/mrc0/projects/Research_Project-MRC158833/UKBiobank/500K_Genetic_data/direct_genotypes/ukb_l2r_chrX_v2.txt /gpfs/mrc0/projects/Research_Project-MRC158833/UKBiobank/500K_Genetic_data/direct_genotypes/ukb_baf_chrX_v2.txt /gpfs/ts0/home/mat206/sandbox/chrX_arms.bed
fam_path=sys.argv[1]
bim_path=sys.argv[2]
l2r_path=sys.argv[3]
baf_path=sys.argv[4]
bed_path=sys.argv[5]
bed_output = open(bed_path + '.outlier_calls.bed','w')
na_output = open(bed_path + '.na_probes.rawcnv','w')
print "Input files: %s %s %s %s %s" % (fam_path, bim_path, l2r_path, baf_path, bed_path)
fam_file = open(fam_path, 'r')
bed_file = open(bed_path, 'r')
bed=bed_file.readlines()
na_lines=dict()
na_lines_rsid=dict()
sid=[]
sex=[]
for fam_ln in fam_file:
    fam = fam_ln.rstrip().split(' ')
    sid.append(fam[1])
    sex.append(fam[4])
l2r_loci = np.zeros([len(bed),len(sid)])
baf_loci = np.zeros([len(bed),len(sid)])
l2r_loci_counts = np.zeros([len(bed),len(sid)])
baf_loci_counts = np.zeros([len(bed),len(sid)])
# Create array of this BED file. Loop through it. Create a list of lists
# Populate counts lists with zeros
probe = 1

print "Processing probes..."
with open(bim_path, 'r') as bim_file:
    with open(l2r_path, 'r') as l2r_file:
        with open(baf_path, 'r') as baf_file:
            for bim_ln, l2r_ln, baf_ln in itertools.izip(bim_file, l2r_file, baf_file):
                bim=bim_ln.rstrip().split('\t')
                for j in range(0, len(bed)):
                    if probe > 1 and bed[j].split('\t')[0] == bim[0] and int(bed[j].split('\t')[1]) <= int(bim[3]) and int(bed[j].rstrip().split('\t')[2]) >= int(bim[3]):
                        l2r=[float(i) if i.isdigit() else i for i in l2r_ln.rstrip().split(' ')]
                        baf=[float(i) if i.isdigit() else i for i in baf_ln.rstrip().split(' ')]
                        for i in range(0, len(sid)):
                            if l2r[i]!='NA':
                                l2r_loci[j,i]+=float(l2r[i])
                                l2r_loci_counts[j,i]+=1
                                # if baf[i]!='NA' and float(baf[i])>0.49 and float(baf[i])<0.51:
                            if baf[i]!='NA':
                                if float(baf[i]) > 0.5:
                                    baf_loci[j,i] += np.exp(2*(1-float(baf[i]))*10)
                                else:
                                    baf_loci[j,i] += np.exp(2*float(baf[i])*10)
                                baf_loci_counts[j,i]+=1
                            # If either L2R or BAF are NA, put in a dict. If key doesn't exist, make it. If it does, append new bim value 
                            if l2r[i]=='NA' or baf[i]=='NA':
                                idx_sid=str(sid[i])
                                val_bim=int(bim[3])
                                val_bim_rsid=bim[1]
                                print "NA: idx_sid = %s, val_bim = %s, val_bim_rsid = %s" % (idx_sid, val_bim, val_bim_rsid)
                                if idx_sid in na_lines:
                                    na_lines[idx_sid].append(val_bim)
                                    na_lines_rsid[idx_sid].append(val_bim_rsid)
                                else:
                                    na_lines[idx_sid] = [val_bim]                          
                                    na_lines_rsid[idx_sid] = [val_bim_rsid]                          
                            if i % 100000 == 0 and i != 0:
                                print 'Sample %s on probe number %s' % (i, probe)
                                probe+=1
                            # if i == 1:
                            #    print 'Probe number %s (not in BED file)' % (probe)
                    else:
                        probe+=1
print na_lines

print "Getting averages..."
line = ''
# Do not add the below to a file. Put them in a 3D array
for i in range(0, len(sid)):
    line = str(sid[i])
    for j in range(0, len(bed)):
        if l2r_loci_counts[j,i]>0:
            l2r_loci[j,i] = l2r_loci[j,i] / l2r_loci_counts[j,i]
        if baf_loci_counts[j,i]>0:
            baf_loci[j,i] = baf_loci[j,i] / baf_loci_counts[j,i]


# R code for outlier detection:
# ss=xf[(xf$V2 > quantile(xf$V2)[4]+(1.5*IQR(xf$V2)) | xf$V2 < quantile(xf$V2)[2]-(1.5*IQR(xf$V2))) & xf$V3 < quantile(xf$V3)[2]-(1.5*IQR(xf$V3)), ]
# ss=xf[xf$V2 > quantile(xf$V2)[4]+l2r_range | xf$V2 < quantile(xf$V2)[2]-l2r_range & xf$V3 < quantile(xf$V3)[2]-(1.5*IQR(xf$V3)), ]

print "Calculating outlier cutoffs..."
l2r_range = 1.5 * np.subtract(*np.percentile(l2r_loci[j,:], [75, 25]))
baf_range = 1.5 * np.subtract(*np.percentile(baf_loci[j,:], [75, 25]))
upper_l2r = np.percentile(l2r_loci[j,:], 75) + l2r_range
lower_l2r = np.percentile(l2r_loci[j,:], 25) - l2r_range
lower_baf = np.percentile(baf_loci[j,:], 25) - baf_range

# Loop through SID inputs, and then BED.
# OUTPUT: BED file - outlier sids
print "Writing files..."
penncnv_sig_prefix="/gpfs/ts0/home/mat206/mat206/500k_CNV_signal_files/all_batches/signals."
for i in range(0, len(sid)):
    line = str(sid[i])
    for j in range(0, len(bed)):
        if (l2r_loci[j,i] >= upper_l2r or l2r_loci[j,i] <= lower_l2r) and baf_loci[j,i] <= lower_baf:
            dup_del, state_cn = '', ''
            if l2r_loci[j,i] >= 0:
                dup_del='dup'
                state_cn='state5,cn=3'
            else: 
                dup_del='del'
                state_cn='state2,cn=1'
            # Print BAF record to a BED file
            # Check SID for presence in na_dict
            if line in na_lines:
                for k in na_lines[line]:
                    # Check probe is within this bed
                    if k >= bed[j].split('\t')[1] and k <= bed[j].split('\t')[2]:
                        # Print PennCNV line k=locus (i.e. chr:k-k) infile.na.rawcnv
                        na_line = 'chr%s:%s-%s numsnp=1 length=1 state2,cn=1 %s%s startsnp=%s endsnp=%s conf=1' % (bed[j].split('\t')[0], k, k, state_cn, penncnv_sig_prefix, sid[i], na_lines_rsid[line], na_lines_rsid[line])
                        na_output.write(na_line)
            # BED file is inbed_chr inbed_pos inbed_end sid l2r baf dup/del infile.outliers.cnv
            bed_line = '%s\t%s\t%s\t%s\t%s\t%s\t%s' % (bed[j].split('\t')[0], bed[j].split('\t')[1], bed[j].rstrip().split('\t')[2], sid[i], l2r_loci[j,i], baf_loci[j,i], dup_del)
            bed_output.write(bed_line)
bed_output.close()
na_output.close()


print "Done."
