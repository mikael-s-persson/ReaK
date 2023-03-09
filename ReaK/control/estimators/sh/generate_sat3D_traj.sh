#!/ bin / bash

BASE_MDL_NAME = "sat3D_airship"

    MDL_DIRECTORY = "models" MDL_FILESTEM = "${MDL_DIRECTORY}/${BASE_MDL_NAME}"

    MAIN_OPTIONS = "--init=${MDL_DIRECTORY}/$1_init.rkx"

    MAIN_OPTIONS =
        "${MAIN_OPTIONS} "
        "--inertia=\"${MDL_FILESTEM}_inertia.rkx\"" MAIN_OPTIONS =
            "${MAIN_OPTIONS} "
            "--Q-matrix=\"${MDL_FILESTEM}_Q.rkx\"" MAIN_OPTIONS =
                "${MAIN_OPTIONS} --R-matrix=\"${MDL_FILESTEM}_R.rkx\""
#MAIN_OPTIONS = "${MAIN_OPTIONS} --R-added=\"${MDL_FILESTEM}_R_added.rkx\""
#MAIN_OPTIONS = \
    "${MAIN_OPTIONS} --IMU-config=\"${MDL_FILESTEM}_IMU_config.rkx\""
    MAIN_OPTIONS =
        "${MAIN_OPTIONS} --xml --protobuf --binary --ssv" MAIN_OPTIONS =
            "${MAIN_OPTIONS} --output=est_results/${BASE_MDL_NAME}/$1"

    SIMTIME_OPTIONS =
        "--time-step=0.01 --start-time=0.0 --end-time=120.0" EXTRA_OPTIONS =
            ""
#EXTRA_OPTIONS = "--gyro"
#EXTRA_OPTIONS = "--IMU"

                echo
            "./estimate_satellite3D --generate-meas-file ${MAIN_OPTIONS} "
            "${SIMTIME_OPTIONS} ${EXTRA_OPTIONS}"./
                estimate_satellite3D-- generate -
            meas - file ${MAIN_OPTIONS} ${SIMTIME_OPTIONS} $ {
  EXTRA_OPTIONS
}
