../bin/render4d -scene area-constant -output sphere-04spp -width 800 -height 800 -spp 4 -cv-rate 0.25 
../bin/render4d -scene area-constant -output sphere-16spp -width 800 -height 800 -spp 16 -cv-rate 0.25 
../bin/render4d -scene area-constant -output sphere-64spp -width 800 -height 800 -spp 64 -cv-rate 0.25 

for i in sphere-*.hdr; do
    [ -f "$i" ] || break
    luminance-hdr-cli -l "$i" --tmo reinhard02 -o "${i%.hdr}.png"
done
