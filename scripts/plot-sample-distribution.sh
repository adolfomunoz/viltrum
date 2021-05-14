../bin/plot-project-2d -function step 1 1 -spp 16 -spp-cv 8 -bins 50 -seed 111 -output samples-discontinuity.svg
../bin/plot-project-2d -function perlin 9 -spp 16 -spp-cv 8 -bins 50 -seed 111 -output samples-perlin.svg
for i in samples-*.svg; do
    [ -f "$i" ] || break
	inkscape "$i" --without-gui --export-pdf="${i%.svg}.pdf"
done
