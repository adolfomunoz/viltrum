../bin/experiment-multiple-scattering arealight occluder -scale 0.2 -absorption 0.1 -scattering 1 -output medium-arealight -width 800 -height 400 -max-spp 32 -cv-rate 0.25 -ground-truth-spp 128


for i in medium-arealight*.hdr; do
    [ -f "$i" ] || break
    luminance-hdr-cli -l "$i" --tmo reinhard02 -o "${i%.hdr}.png"
done
for i in medium-arealight*.svg; do
    [ -f "$i" ] || break
	inkscape "$i" --without-gui --export-pdf="${i%.svg}.pdf"
done