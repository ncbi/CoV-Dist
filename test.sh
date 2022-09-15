python covdist.py dist -p data/test/full.list -r data/NC_045512.2.fasta -o data/test/testout
python covdist.py plot -d data/test/testout/distance.txt -m data/test/metadata.tsv -c date -o data/test/testout
