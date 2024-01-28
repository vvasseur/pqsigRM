export PYTHONUNBUFFERED=1

cd pqsigrm613/
make keygen_and_sign signatures.so
cd -

echo "Generating a keypair and 100000 signatures"
pqsigrm613/keygen_and_sign 100000 pqsigRM

echo "Recovering a permutation that reveals two levels of the recursive (U|U+V) structure"
python find_U_UV.py pqsigRM-pk pqsigRM-sigs pqsigRM-U_UV_permutation

echo "Verifying the permutation"
python verify_U_UV_permutation.py pqsigRM-sk pqsigRM-U_UV_permutation
