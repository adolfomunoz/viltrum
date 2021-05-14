../bin/experiment-single-scattering occluder -scale 0.1 -absorption 0.01 -scattering 0.19 -output medium -width 800 -height 400 -max-spp 16 -cv-rate 0.01 -spp-pixel 2

for i in medium-*.hdr; do
    [ -f "$i" ] || break
    luminance-hdr-cli -l "$i" --tmo reinhard02 -o "${i%.hdr}.png"
done
for i in medium-*.svg; do
    [ -f "$i" ] || break
	inkscape "$i" --without-gui --export-pdf="${i%.svg}.pdf"
done
