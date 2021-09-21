tempfile=$(mktemp)
exec 3>"$tempfile"
exec 4<"$tempfile"
rm "$tempfile"
matlab -nodisplay -nosplash -r "schnackenberg_2dfunc(${1});exit" >&3
mail -s "[fdmsim] schnackenberg ${1} done" yue.liu@maths.ox.ac.uk <&4
