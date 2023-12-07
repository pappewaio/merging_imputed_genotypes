infile=${1}
outfile=${2}

head -n1 "$(head -n1 ${infile})" | sed 's/\#//g' | sed 's/\[[0-9]*\]//g' > "${outfile}"
while IFS= read -r line; do
  tail -n+2 "${line}"
done < ${infile} >> "${outfile}"

