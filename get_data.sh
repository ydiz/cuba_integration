# tmp=$(sed -n -r  "s/v: \[(.*)\]/\1/p" <<< "$(tail -n 4 output-6916.txt)")
# read -r r z t <<< "$tmp"
#
# echo r: [0, "$r"]
# echo z: [$((-z)), "$z"]
# echo t: [0, "$t"]
#
sed -n -r "s/.*integral = ([0-9.-]+).*/\1/p" output-6916.txt
