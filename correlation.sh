cd pqsigrm613/
make keygen_and_sign get_matched_pairs signatures.so
cd -

echo "Generating a keypair and 1000 signatures"
pqsigrm613/keygen_and_sign 1000 pqsigRM-correlation

echo "Finding matched pairs by computing correlations in signature bits"
pqsigrm613/get_matched_pairs pqsigRM-correlation-sigs > pqsigRM-correlation-pairs

echo "Counting the good matched pairs found"
python count_good_matched_pairs.py pqsigRM-correlation-sk pqsigRM-correlation-pairs
