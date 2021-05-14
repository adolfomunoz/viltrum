RESOLUTION=8
if [ -n ${1} ]; then
	RESOLUTION=${1}
fi

MONTECARLO_SPP="0 "
MONTECARLO_SPP_BINS="0 "
BINS=1011
S=16
B=1
for i in $(seq 1 ${RESOLUTION}); do
	MONTECARLO_SPP+="${S} "
	MONTECARLO_SPP_BINS+="${B} "
	S=$(( ${S}*2 ))
	B=$(( ${B}*2 ))
done
#16 32 64 128 256 512 1024 2048 4096

CV_ITERATIONS=""
CV_ITERATIONS_BINS=""
S=1
B=167
for m in ${MONTECARLO_SPP}; do
	N=$(( (${m}-9)/6 ))
	N=$(( ${N}<0?0:${N} ))
	CV_ITERATIONS+="${N} "
done
for m in ${MONTECARLO_SPP_BINS}; do
	N=$(( (${BINS}*${m}-9)/384 ))
	N=$(( ${N}<0?0:${N} ))
	CV_ITERATIONS_BINS+="${N} "
done

FUNCTIONS="freqband_3 lightmedia_1.5 "
for i in *.jpg; do
	[ -f "$i" ] || break
	FUNCTIONS+="image_${i} "
done

FILENAMES=""
FILENAMES_BINS=""

for f in ${FUNCTIONS}; do
    echo "$f - ${f//_/ } - ${f//./_}"
	../bin/plot-2d -function ${f//_/ } -output "controlvariate-${f//./_}.svg" -iterations 300 -plot-resolution 1000
	inkscape "controlvariate-${f//./_}.svg" --without-gui --export-pdf="controlvariate-${f//./_}.pdf" > /dev/null 2> /dev/null
    ../bin/map-convergence-project-2d -function ${f//_/ }  -bins 1  -montecarlo-error-average 100 -output "map-convergence-${f//./_}.svg" -min-plottable-error 1.e-20 -montecarlo-spp ${MONTECARLO_SPP} -cv-iterations ${CV_ITERATIONS} -min-plottable-error 1.e-10
	inkscape "map-convergence-${f//./_}.svg" --without-gui --export-pdf="map-convergence-${f//./_}.pdf" > /dev/null 2> /dev/null
	FILENAMES+="map-convergence-${f//./_}.txt "
	../bin/map-convergence-project-2d -function ${f//_/ }  -bins ${BINS}  -montecarlo-error-average 100 -output "map-convergence-${f//./_}-${BINS}bins.svg" -min-plottable-error 1.e-20 -montecarlo-spp ${MONTECARLO_SPP_BINS} -cv-iterations ${CV_ITERATIONS_BINS} -min-plottable-error 1.e-10
	inkscape "map-convergence-${f//./_}-${BINS}bins.svg" --without-gui --export-pdf="map-convergence-${f//./_}-${BINS}bins.pdf" > /dev/null 2> /dev/null
	FILENAMES_BINS+="map-convergence-${f//./_}-${BINS}bins.txt "
done

../bin/map-convergence-average -bins 1 -output "map-convergence-average.svg" -montecarlo-spp ${MONTECARLO_SPP} -cv-iterations ${CV_ITERATIONS} -min-plottable-error 1.e-10 -datafiles ${FILENAMES}
inkscape "map-convergence-average.svg" --without-gui --export-pdf="map-convergence-average.pdf" > /dev/null 2> /dev/null
../bin/map-convergence-average -bins ${BINS} -output "map-convergence-average-${BINS}bins.svg" -montecarlo-spp ${MONTECARLO_SPP_BINS} -cv-iterations ${CV_ITERATIONS_BINS} -min-plottable-error 1.e-10 -datafiles ${FILENAMES_BINS}
inkscape "map-convergence-average-${BINS}bins.svg" --without-gui --export-pdf="map-convergence-average-${BINS}bins.pdf" > /dev/null 2> /dev/null