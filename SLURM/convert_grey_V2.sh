#!/bin/sh

for i in *.avi;
  do name=`echo "$i" | cut -d'.' -f1`
  echo "$name"
  ffmpeg -i "$i" -vcodec rawvideo "${name}_grey.avi"
  rm "${name}.avi"
done

suffix="_grey"
for i in *.avi;
  do name=`echo "$i" | cut -d'.' -f1`
  echo "$name"
  cname="${name}"
  nname=${cname%$suffix}
  cp "${name}.avi" "${nname}.avi"
  rm "${name}.avi"
done
