#!/bin/bash

TEMP_FNAME=/tmp/acf_names.$RANDOM

# git diff --name-only 08ed2f7052543c37e6333a105392a511d0e38f02 6c1fa58cedfc3f4c6154e29456d4a626a030c915 \
# | grep -E "(\.hpp$)|(\.cpp$)|(\.tpp$)|(\.h$)|(\.c$)|(\.inl$)" \
# | awk '{
#   if(NF==1) {
#     print $1 
#   } else {
#     if($1 != "deleted:")
#       print $2;
#   }
# }' > $TEMP_FNAME

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
