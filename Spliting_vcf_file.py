# Executing the python script for adding the depth columns to the VCF file.
# Command: python3 ./data/vcf_depth_wes.py KHAIGPRX1113

####### Depth info spliting in python ######
import pandas as pd
import sys

# read in the VCF file
vcf_file = "/{}/{}_final.vcf"
vcf_data = pd.read_csv(vcf_file.format(sys.argv[1], sys.argv[1]), sep="\t", comment='#', header=None)
print(vcf_data.shape)

# extract the column names from the FORMAT column
format_col = vcf_data[8].str.split(':', expand=True)
format_col.columns = format_col.iloc[0]
format_col = format_col[1:]

# extract the sample data from the SAMPLE column
sample_data = vcf_data[9].str.split(':', expand=True)
sample_data.columns = format_col.columns

# concatenate the original data with the new FORMAT and SAMPLE columns
output_data = pd.concat([vcf_data.iloc[:, :8], sample_data], axis=1)
print(output_data)
dropped = output_data.drop(['RBQ', 'ABQ'], axis=1)
print(dropped)
print(dropped.columns)
print(dropped.shape)

column_names = ['CHROM',  'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO','GT', 'GQ', 'SDP', 'DP', 'RD', 'AD', 'FREQ', 'PVAL', 'RDF', 'RDR', 'ADF', 'ADR']
dropped.columns=column_names
print(dropped.columns)
# write the output data to an Excel file
# output_file = "/run/media/administrator/0c7b01f3-0653-4179-8281-e9396009298b/NGS/KHAIGPRX/main/{}_Depth_vcf.xlsx"
output_file = "/{}/{}_final_depth.tsv"
# output_file = "/media/khserver/Seagate Expansion Drive/Final_VCFs/29_04_2023_VCFs/{}/{}_Depth.vcf"
dropped.to_csv(output_file.format(sys.argv[1], sys.argv[1]), sep="\t", index=False)

# After running the python script 
# Run these shell commands for adding the headers to the output file. 
# cat data/final_vcf/KHAIGPRX1113/KHAIGPRX1113_final.vcf | grep "##" > data/final_vcf/KHAIGPRX1113/KHAIGPRX1113_header.txt
# cat data/final_vcf/KHAIGPRX1113/KHAIGPRX1113_header.txt data/final_vcf/KHAIGPRX1113/KHAIGPRX1113_Depth.vcf > data/final_vcf/KHAIGPRX1113/KHAIGPRX1113_final_Depth.vcf


# apo=pd.read_csv("KHGLBS598_vep.vcf", sep="\t", header=None, comment="#")
# apo.columns = ['CHROM',  'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'SAMPLE']
# # extract the column names from the FORMAT column
# format_col = apo['FORMAT'].str.split(':', expand=True)
# format_col.columns = format_col.iloc[0]
# format_col = format_col[1:]

# # extract the sample data from the SAMPLE column
# sample_data = apo['SAMPLE'].str.split(':', expand=True)
# sample_data.columns = format_col.columns

# # concatenate the original data with the new FORMAT and SAMPLE columns
# output_data = pd.concat([apo.iloc[:, :8], sample_data], axis=1)
# output_data.to_excel("KHGLBS598_vep.xlsx", index=False)
