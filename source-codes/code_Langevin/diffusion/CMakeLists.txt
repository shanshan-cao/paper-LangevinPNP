cmake_minimum_required (VERSION 2.6)
project (diffusion Fortran)

add_subdirectory (src)

# install dependencies
# install (FILES FPcoeff.dat massT.dat mc_glauber.dat parameters_df.dat dNg_over_dt_cD6.dat dNg_over_dt_bD6.dat JetCtl.dat JetData.dat DESTINATION ${PROJECT_NAME})
install (FILES bulk.dat mc_glauber.dat parameters_df.dat dNg_over_dt_cD6.dat dNg_over_dt_bD6.dat JetCtl.dat JetData.dat kappa_tp.dat DESTINATION ${PROJECT_NAME})

