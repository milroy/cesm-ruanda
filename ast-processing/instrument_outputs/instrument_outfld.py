import os
import sys
import glob
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument("--src", dest="srcPath", required=True, type=str)

args = parser.parse_args()
srcPath = args.srcPath

outFiles = ['uwshcu.F90', 'sslt_rebin.F90', 'cloud_diagnostics.F90',
'mo_airglow.F90', 'modal_aero_newnuc.F90', 'modal_aer_opt.F90',
'gw_drag.F90', 'cldwat2m_macro.F90', 'stepon.F90', 'cfc11star.F90',
'mo_setext.F90', 'eddy_diff.F90', 'tropopause.F90',
'cloud_fraction.F90', 'modal_aero_wateruptake.F90', 'mo_photo.F90',
'aerodep_flx.F90', 'aircraft_emit.F90', 'mo_setinv.F90',
'cospsimulator_intr.F90', 'slingo.F90', 'camsrfexch.F90',
'modal_aero_gasaerexch.F90', 'radiation_data.F90', 'noy_ubc.F90',
'cam_history.F90', 'constituents.F90', 'tracer_cnst.F90', 'radsw.F90',
'modal_aero_coag.F90', 'modal_aero_calcsize.F90', 'mo_lightning.F90',
'rad_constituents.F90', 'chemistry.F90', 'tidal_diag.F90',
'constituent_burden.F90', 'mo_setsox.F90', 'cldwat.F90',
'micro_mg_cam.F90', 'cam3_ozone_data.F90', 'mo_waccm_hrates.F90',
'gcr_ionization.F90', 'mo_neu_wetdep.F90', 'mo_extfrc.F90',
'prescribed_ghg.F90', 'tracer_srcs.F90', 'zm_conv_intr.F90',
'mo_airplane.F90', 'aero_model.F90', 'rate_diags.F90',
'sox_cldaero_mod.F90', 'tracers.F90', 'linoz_data.F90', 'mo_sulf.F90',
'conv_water.F90', 'convect_deep.F90', 'mo_imp_sol.F90',
'microp_aero.F90', 'cam_history_support.F90',
'physpkg.F90', 'prescribed_volcaero.F90', 'mo_gas_phase_chemdr.F90',
'prescribed_aero.F90', 'mo_exp_sol.F90', 'mo_cph.F90', 'mo_sad.F90',
'vertical_diffusion.F90', 'mo_chm_diags.F90', 
'subcol.F90', 'aoa_tracers.F90', 'macrop_driver.F90',
'prescribed_strataero.F90', 'radlw.F90', 'prescribed_ozone.F90',
'dp_coupling.F90', 'radiation.F90', 'oldcloud.F90',
'cam_diagnostics.F90', 'cloud_cover_diags.F90', 'convect_shallow.F90',
'lin_strat_chem.F90', 'mo_aurora.F90', 'ndrop.F90',
'aer_rad_props.F90', 'atm_comp_mct.F90']

def readF(path):
	with open(path, 'r') as f:
		for line in f:
			yield line


for file in outFiles:
	walk = readF(srcPath + '/' + file)

	with open('instmt/' + file, 'w') as output:
		prevL = None
		for l in walk:
			output.write(l)

			lstrip = l.strip()
			if not lstrip.startswith('!'):
				if 'call outfld' in lstrip:
					try:
						splits = lstrip.split('call outfld(')[1].split(',')
						name = splits[0].replace('call outfld(', '').replace('\n', '').strip()
					except:
						splits = lstrip.split('call outfld (')[1].split(',')
						name = splits[0].replace('call outfld (', '').replace('\n', '').strip()					
					varname = splits[1].replace('\n', '').strip()
					outFMT = "WRITE(*,'(4A)') 'NAME: ', " + name + ", " + \
							"' FILE: " + file + "', " + " ' LINE: " + lstrip.replace("'", "").replace('\n', '') + "'" + " \n"
					output.write(outFMT)

