../bin/convergence-project-2d -function step 1 1 -max-samples 200000 -test-order -error-size-weight 0.00001 -number-of-cv-tests 4 -adaptive-iterations 8 -bins 100 -output discontinuity_samples.svg
inkscape discontinuity_samples.svg --without-gui --export-pdf=discontinuity_samples.pdf
../bin/convergence-project-2d -function step 1 1 -max-time 0.2 -test-order -error-size-weight 0.00001 -number-of-cv-tests 4 -adaptive-iterations 8 -bins 100 -plot-resolution 25 -output discontinuity_time.svg
inkscape discontinuity_time.svg --without-gui --export-pdf=discontinuity_time.pdf
../bin/convergence-project-2d -function image flowers.jpg -max-samples 400000 -test-order -error-size-weight 0.00001 -number-of-cv-tests 4 -adaptive-iterations 8 -bins 100 -output flower_samples.svg
inkscape flower_samples.svg --without-gui --export-pdf=flower_samples.pdf
../bin/convergence-project-2d -function image flowers.jpg -max-time 0.4 -test-order -error-size-weight 0.00001 -number-of-cv-tests 4 -adaptive-iterations 8 -bins 100 -plot-resolution 25 -output flower_time.svg
inkscape flower_time.svg --without-gui --export-pdf=flower_time.pdf

