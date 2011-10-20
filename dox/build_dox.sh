#!/bin/bash
build_alllib()
{
doxygen ReaK_Doxyfile
qhelpgenerator html/index.qhp -o ReaK_dox.qch
}
build_core()
{
cd core
doxygen ../ReaKcore_Doxyfile
qhelpgenerator html/index.qhp -o ../ReaKcore_dox.qch
cd ..
}
build_misc()
{
cd misc
doxygen ../ReaKmisc_Doxyfile
qhelpgenerator html/index.qhp -o ../ReaKmisc_dox.qch
cd ..
}
build_math()
{
cd math
doxygen ../ReaKmath_Doxyfile
qhelpgenerator html/index.qhp -o ../ReaKmath_dox.qch
cd ..
}
build_allcore()
{
cd allcore
doxygen ../ReaKallcore_Doxyfile
qhelpgenerator html/index.qhp -o ../ReaKallcore_dox.qch
cd ..
}
build_mbd_kte()
{
cd mbd_kte
doxygen ../ReaKmbd_kte_Doxyfile
qhelpgenerator html/index.qhp -o ../ReaKmbd_kte_dox.qch
cd ..
}
build_ctrl_sys()
{
cd ctrl_sys
doxygen ../ReaKctrl_sys_Doxyfile
qhelpgenerator html/index.qhp -o ../ReaKctrl_sys_dox.qch
cd ..
}
build_pp()
{
cd pp
doxygen ../ReaKpp_Doxyfile
qhelpgenerator html/index.qhp -o ../ReaKpp_dox.qch
cd ..
}
clean_dox()
{
rm -R html
rm -R core/html
rm -R core/pdflatex
rm -R misc/html
rm -R misc/pdflatex
rm -R math/html
rm -R math/pdflatex
rm -R allcore/html
rm -R allcore/pdflatex
rm -R mbd_kte/html
rm -R mbd_kte/pdflatex
rm -R ctrl_sys/html
rm -R ctrl_sys/pdflatex
rm -R pp/html
rm -R pp/pdflatex
}
build_all()
{
clean_dox
build_alllib
build_core
build_misc
build_math
build_allcore
build_mbd_kte
build_ctrl_sys
build_pp
}
case $1 in
  "all") build_all;;
  "core") build_core;;
  "misc") build_misc;;
  "math") build_math;;
  "allcore") build_allcore;;
  "alllib") build_alllib;;
  "mbd_kte") build_mbd_kte;;
  "ctrl_sys") build_ctrl_sys;;
  "pp") build_pp;;
  "clean") clean_dox;;
esac



