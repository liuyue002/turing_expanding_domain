tempfile=$(mktemp)
exec 3>"$tempfile"
exec 4<"$tempfile"
rm "$tempfile"
matlab -nodisplay -nosplash -r "cdima_2dfunc(${1});exit" >&3
mail -s "[fdmsim] cdima ${1} done" yue.liu@maths.ox.ac.uk <&4
