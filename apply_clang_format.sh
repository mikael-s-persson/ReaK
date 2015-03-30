#!/bin/bash

TEMP_FNAME=/tmp/acf_names.$RANDOM

git status \
| grep -E "(\.hpp$)|(\.cpp$)|(\.tpp$)|(\.h$)|(\.c$)|(\.inl$)" \
| awk '{
  if(NF==1) {
    print $1 
  } else {
    if($1 != "deleted:")
      print $2;
  }
}' > $TEMP_FNAME

while read line ; do
  clang-format-3.5 -style=file -i $line
done < $TEMP_FNAME

rm $TEMP_FNAME
