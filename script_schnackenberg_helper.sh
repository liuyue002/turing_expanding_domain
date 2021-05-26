matlab -nodisplay -nosplash -r "schnackenberg_2dfunc(${1});exit"
mail -s "[fdmsim] schnackenberg ${1} done" yue.liu@maths.ox.ac.uk < /dev/null
