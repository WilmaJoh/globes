import gcm_rocke3d
import numpy as np


gcm_files = np.array([('01_Control_ANN5000', 'ProxCenb04bHR_TL.nc'),
                      ('02_Thermo_ANN4000', 'ProxCenb04eHR_TL.nc'),
                      ('06_Archean-Low_ANN4900-4999', 'ProxCenb04bf_TL_DSA_CLDDIAG.nc'),
                      ('07_Archean-Med_ANN4900-4999', 'ProxCenb04f_TL_CLDDIAG.nc'),
                      ('08_Archean-Med_NoCH4_ANN5002-5101', 'ProxCenb04f3TL_CLDDIAG.nc'),
                      ('09_Archean-High_ANN6002-6101', 'ProxCenb04arcb_TL_DSA_CLDDIAG.nc')])



for (name, suffix) in gcm_files:
    gcm_rocke3d.convertgcm(sim_name = name, file_suffix = suffix, fileout = name + '_OP_psg.dat')