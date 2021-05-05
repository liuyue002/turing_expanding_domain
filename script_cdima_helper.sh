matlab -nodisplay -nosplash -r "cdima_2dfunc(${1});exit"
mail -s "cdima ${1} done" yue.liu@maths.ox.ac.uk < /dev/null
