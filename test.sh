# test without voc
python covdist.py dist -p data/test/full.list -r data/NC_045512.2.fasta -o data/test/testout
python covdist.py plot -d data/test/testout/distance.txt -m data/test/metadata.tsv -o data/test/testout

# test with voc
python covdist.py dist -p data/test/full.list -r data/NC_045512.2.fasta -v voc_vcf -o data/test/testout_voc
# pcoa
python /data/yangy34/projects/CoV-Dist/covdist.py plot -d data/test/testout_voc/distance.txt -m data/test/metadata.tsv -v data/test/metadata_voc.tsv -o data/test/testout_voc
# mds
python /data/yangy34/projects/CoV-Dist/covdist.py plot -d data/test/testout_voc/distance.txt -a MDS -m data/test/metadata.tsv -v data/test/metadata_voc.tsv -o data/test/testout_voc
# t-SNE
python /data/yangy34/projects/CoV-Dist/covdist.py plot -d data/test/testout_voc/distance.txt -a TSNE -m data/test/metadata.tsv -v data/test/metadata_voc.tsv -o data/test/testout_voc
