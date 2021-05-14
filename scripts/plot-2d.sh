for i in *.jpg; do
    [ -f "$i" ] || break
    echo "$i"
    ../bin/plot-2d -function image ${i} -output "controlvariate-${i%.jpg}.svg" -iterations 300 -plot-resolution 1000
	inkscape "controlvariate-${i%.jpg}.svg" --without-gui --export-pdf="controlvariate-${i%.jpg}.pdf"
    ../bin/convergence-project-2d -function image ${i} -test-quadrature -cv-ratio 0.5 -number-of-cv-tests 5 -bins 1 -montecarlo-error-average 100 -max-samples 250000 -output "convergence-relative-${i%.jpg}.svg" -plot-resolution 12 -width 200 -spacing 20 -min-plottable-error 1.e-20
	inkscape "convergence-relative-${i%.jpg}.svg" --without-gui --export-pdf="convergence-relative-${i%.jpg}.pdf"
    ../bin/convergence-project-2d -function image ${i} -test-quadrature -adaptive-iterations 2 -number-of-cv-tests 5 -bins 1 -montecarlo-error-average 100 -max-samples 250000 -output "convergence-absolute-${i%.jpg}.svg" -plot-resolution 12 -width 200 -spacing 20 -min-plottable-error 1.e-20
	inkscape "convergence-absolute-${i%.jpg}.svg" --without-gui --export-pdf="convergence-absolute-${i%.jpg}.pdf"
done
