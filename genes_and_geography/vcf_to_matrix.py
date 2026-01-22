import pysam
import numpy as np
import time
from sklearn import decomposition
import pandas as pd
import sys


pysam.set_verbosity(0)

data_folder = 'data'
vcf_filename = f'{data_folder}/ALL.chr22.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz'
panel_filename = f'{data_folder}/phase1_integrated_calls.20101123.ALL.panel'
csv_output_filename = 'genotype_matrix.csv'
# csv_output_filename = 'big_genotype_matrix.csv'

genotypes = list()
samples = list()
variant_ids = list()
# testing_sample_interval = 100
sample_interval = 100
# sample_interval = 1
approx_total_size = 494328
progress_interval = approx_total_size // 10
print(f'progress_interval: {progress_interval}')
print('sample_interval:', sample_interval)

start_time = time.time()
with pysam.VariantFile(vcf_filename) as vcf_reader:
    for iter, record in enumerate(vcf_reader):
        # Sample across the chromesome (too many bases to read all of them)
        if iter % sample_interval == 0:
            alleles = [record.samples[i].allele_indices for i in record.samples]
            samples = [sample for sample in record.samples]
            variant_ids.append(record.id)
            genotypes.append(alleles)

        if iter % (progress_interval) == 0:
            print(f'Progress: {round(100 * iter / approx_total_size, 2)}%')
        # if iter >= sample_size * sample_interval:
        #     break
end_time = time.time()
print(f'VCF reading time: {end_time - start_time} seconds')
print(f'Number of variants sampled: {len(variant_ids)}')


with open(panel_filename, 'r') as panel_file:
    labels = dict()
    for line in panel_file:
        line = line.strip().split('\t')
        sample_id = line[0]
        population_code = line[1] # Corresponds to countries/region
        labels[sample_id] = population_code

# print(labels)

genotypes = np.array(genotypes)
# print(genotypes.shape)

matrix = np.count_nonzero(genotypes, axis=2)
print(matrix.shape)

# Transpose
matrix = matrix.T
print(f'matrix.shape: {matrix.shape}')

start_time = time.time()
pca = decomposition.PCA(n_components=2)
end_time = time.time()
print(f'PCA initialization time: {end_time - start_time} seconds')

pca.fit(matrix)
print(f'pca.singular_values_: {pca.singular_values_}')
to_plot = pca.transform(matrix)
print(f'to_plot: {to_plot.shape}')

df = pd.DataFrame(matrix, columns=variant_ids, index=samples)
# Map the population code to the dataframe using the sample_ids in labels. Ex: HG00096, HG00097, ...
df['Population code'] = df.index.map(labels)
# df.rename(columns={df.columns[0]: 'Sample'}, inplace=True)

print(df.head()) 
df.to_csv(f'{data_folder}/{csv_output_filename}')