#!/bin/bash
if [ $1 == "alllib" ]
then
assistant-qt4 -collectionFile ReaK_dox_collection.qhc -register ReaK_dox.qch
assistant-qt4 -collectionFile ReaK_dox_collection.qhc
assistant-qt4 -collectionFile ReaK_dox_collection.qhc -unregister ReaK_dox.qch
else
assistant-qt4 -collectionFile ReaK_dox_collection.qhc -register ReaK${1}_dox.qch
assistant-qt4 -collectionFile ReaK_dox_collection.qhc
assistant-qt4 -collectionFile ReaK_dox_collection.qhc -unregister ReaK${1}_dox.qch
fi
