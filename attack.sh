export PYTHONUNBUFFERED=1

cd pqsigrm613/
make keygen_and_sign get_matched_pairs
cd -

echo "Generating a keypair and 100000 signatures"
pqsigrm613/keygen_and_sign 100000 pqsigRM

echo "Finding matched pairs by computing correlations in signature bits"
pqsigrm613/get_matched_pairs pqsigRM-sigs > pqsigRM-pairs

echo "Recovering a permutation that reveals the (U|U+V) structure"
sage find_U_UV.sage pqsigRM-pk pqsigRM-pairs pqsigRM-U_UV_permutation

echo "Verifying the permutation"
python verify_U_UV_permutation.py pqsigRM-sk pqsigRM-U_UV_permutation
