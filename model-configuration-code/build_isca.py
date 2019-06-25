import os
import glob
import numpy as np
import shutil
import string
import sys
from isca import SocratesCodeBase, IscaCodeBase, GreyCodeBase, DiagTable, Experiment, Namelist, GFDL_BASE
from isca.util import exp_progress

casename0='SW_absorb_'
radiation_scheme = 'grey'
append=sys.argv[1]
casename=casename0+append
#casename = sys.argv[1]
print('building case:'+casename0)
isca_caseroot = '/glade/u/home/'+os.environ["USER"]+'/Isca_caseroot/'
builddir = isca_caseroot+casename0
print('in:'+builddir)
#executable_name=casename0+'.x'
executable_name='isca.x'
NCORES= int(sys.argv[2])
print('Using '+str(NCORES)+' cores')
#RESOLUTION = 'T42', 40
RESOLUTION = 'T21', 25
todo = sys.argv[3]
debug=False
run_ddt=False
print(len(sys.argv))
if len(sys.argv) >= 5:
	mothercase = sys.argv[4]
else:
	mothercase = None


#Step 1: clean up if required
objdir=builddir.replace('/','_')+'_Isca'
if 'clean' in todo:
	if os.path.exists(os.environ['GFDL_WORK']+'/codebase/'+objdir):
		print('clean up '+os.environ['GFDL_WORK']+'/codebase/'+objdir)
		shutil.rmtree(os.environ['GFDL_WORK']+'/codebase/'+objdir)
	if os.path.exists(os.environ['GFDL_WORK']+'/experiments/'+casename):
		print('clean up '+os.environ['GFDL_WORK']+'/experiments/'+casename)
		shutil.rmtree(os.environ['GFDL_WORK']+'/experiments/'+casename)
	if os.path.exists(builddir):
		print('clean up '+builddir)
		shutil.rmtree(builddir)

if 'reinit' in todo:
	if os.path.exists(os.environ['GFDL_DATA']+casename):
		shutil.rmtree(os.environ['GFDL_DATA']+casename)

#Step 2. Provide the necessary input files and user source code:
if not os.path.exists(builddir):
	os.makedirs(builddir)
if not os.path.exists(builddir+'/Isca'):
	os.symlink(os.environ['GFDL_BASE'],builddir+'/Isca')
#cb = IscaCodeBase.from_directory(builddir+'/Isca',storedir=builddir)  # for saving obj files in Isca_caseroot/
if radiation_scheme == 'grey':
	cb = GreyCodeBase.from_directory(builddir+'/Isca')
	do_grey = True
	do_rrtm = False
	do_socrates = False
if radiation_scheme == 'rrtm':
	cb = IscaCodeBase.from_directory(builddir+'/Isca')
	do_grey = False
	do_rrtm = True
	do_socrates = False
if radiation_scheme == 'socrates':
	cb = SocratesCodeBase.from_directory(builddir+'/Isca')
	do_grey = False
	do_rrtm = False
	do_socrates = True




if 'compile' in todo or 'configure' in todo:
	if mothercase is not None:
		mom_builddir = isca_caseroot+mothercase
		mom_objdir=mom_builddir.replace('/','_')+'_Isca'
		if os.path.exists(os.environ['GFDL_WORK']+'/codebase/'+mom_objdir):
			print('found mother case: '+mothercase+', copying .......')
			#os.system("cp -rf "+os.environ['GFDL_WORK']+'/codebase/'+mom_objdir+'/code '+os.environ['GFDL_WORK']+'/codebase/'+objdir)
			os.system("cp -rf "+os.environ['GFDL_WORK']+'/codebase/'+mom_objdir+'/build '+os.environ['GFDL_WORK']+'/codebase/'+objdir)
	if not os.path.exists(builddir+'/srcmod_all'):
		os.mkdir(builddir+'/srcmod_all')
	if os.path.exists('./isca_srcmod'):
		if os.path.exists(builddir+'/srcmod'):
			shutil.rmtree(builddir+'/srcmod')
		os.system("cp -rf ./isca_srcmod "+builddir+"/srcmod")
	os.system("cp -f "+builddir+"/srcmod/* "+builddir+"/srcmod_all/")
	extra_pathnames_all_tmp = os.listdir(builddir+'/srcmod_all/')
	extra_pathnames_tmp = os.listdir(builddir+'/srcmod/')
	extra_pathnames = [builddir+'/srcmod/' + s for s in extra_pathnames_tmp]
	print('Including the following modified CODEs:')
	print(extra_pathnames)
if 'compile' in todo:
	print('deleting the following OBJ:')
	print(extra_pathnames_all_tmp)
	for s in extra_pathnames_all_tmp[:]:
		ss = os.path.splitext(s)[0]
		for sss in glob.glob(os.environ["GFDL_WORK"]+"/codebase/"+objdir+"/build/isca/"+ss+"*"):
			os.remove(sss)

#Step 3. input files and run length
print('Step 3: setting input files and run length')
runlen=30
totalrun = 600
eachrun = 600
run_day=runlen
run_hour=0
ozonefilename='ozone_uniform7'  # default: 'ozone_1990'
inputfiles = [os.path.join(GFDL_BASE,'input/rrtm_input_files/',ozonefilename+'.nc')]
# copy/link qflux file to builddir
#qflux_file_name = 'merlis_schneider_30_16'
qflux_file_name = 'off'
if qflux_file_name!='off':
	if not os.path.exists(builddir+'/'+qflux_file_name+'.nc'):
		shutil.copy(os.environ['GFDL_BASE']+'/ictp-isca-workshop-2018/experiments/earth_projects/project_5_qfluxes/'+qflux_file_name+'.nc',builddir+'/'+qflux_file_name+'.nc')

#Step 4. Define the diagnostics we want to be output from the model
print('Step 4: setting output variables')
avgflag = True

diag = DiagTable()
diag.add_file('atmos_monthly', 30, 'days', time_units='days')
diag.add_file('atmos_daily', 1, 'days', time_units='days')
if 'debug' in todo:
	debug=True
	run_ddt=True
	avgflag = False
	diag.add_file('atmos_hourly', 1, 'hours', time_units='hours')
	#NCORES=1

#Tell model which diagnostics to write to those files
diag.add_field('dynamics', 'ps', time_avg=avgflag)
diag.add_field('dynamics', 'bk')
diag.add_field('dynamics', 'pk')
diag.add_field('dynamics', 'zsurf')
diag.add_field('dynamics', 'sphum', time_avg=avgflag)
diag.add_field('dynamics', 'ucomp', time_avg=avgflag)
diag.add_field('dynamics', 'vcomp', time_avg=avgflag)
diag.add_field('dynamics', 'vor', time_avg=avgflag)
diag.add_field('dynamics', 'div', time_avg=avgflag)
diag.add_field('dynamics', 'height', time_avg=avgflag)
diag.add_field('dynamics', 'omega', time_avg=avgflag)
diag.add_field('dynamics', 'temp', time_avg=avgflag)

diag.add_field('atmosphere', 'dt_tg_convection', time_avg=avgflag)
diag.add_field('atmosphere', 'dt_tg_diffusion', time_avg=avgflag)
diag.add_field('mixed_layer', 't_surf', time_avg=avgflag)
diag.add_field('mixed_layer', 'flux_r', time_avg=avgflag)
diag.add_field('mixed_layer', 'flux_t', time_avg=avgflag)

if radiation_scheme == 'grey':
	diag.add_field('two_stream', 'flux_sw', time_avg=avgflag)    # 3D for two stream
	diag.add_field('two_stream', 'flux_lw', time_avg=avgflag)
	diag.add_field('two_stream', 'tdt_rad', time_avg=avgflag)
	diag.add_field('two_stream', 'olr', time_avg=avgflag)
	diag.add_field('two_stream', 'swdn_sfc', time_avg=avgflag)
	diag.add_field('two_stream', 'lwdn_sfc', time_avg=avgflag)
	diag.add_field('two_stream', 'lwup_sfc', time_avg=avgflag)
	diag.add_field('two_stream', 'net_lw_surf', time_avg=avgflag)
if radiation_scheme == 'rrtm':
	diag.add_field('rrtm_radiation', 'flux_sw', time_avg=avgflag) # net surface flux for rrtm
	diag.add_field('rrtm_radiation', 'flux_lw', time_avg=avgflag)
	diag.add_field('rrtm_radiation', 'tdt_rad', time_avg=avgflag)
	diag.add_field('rrtm_radiation', 'olr', time_avg=avgflag)
if radiation_scheme == 'socrates':
	diag.add_field('socrates', 'soc_tdt_lw', time_avg=avgflag) 
	diag.add_field('socrates', 'soc_tdt_sw', time_avg=avgflag)
	diag.add_field('socrates', 'soc_tdt_rad', time_avg=avgflag)
	diag.add_field('socrates', 'soc_surf_flux_lw', time_avg=avgflag)
	diag.add_field('socrates', 'soc_surf_flux_sw', time_avg=avgflag)
	diag.add_field('socrates', 'soc_olr', time_avg=avgflag)
	diag.add_field('socrates', 'soc_toa_sw', time_avg=avgflag)


diag.add_field('dynamics', 'ucomp_vcomp', time_avg=avgflag)
diag.add_field('dynamics', 'ucomp_sq', time_avg=avgflag)
diag.add_field('dynamics', 'ucomp_omega', time_avg=avgflag)
diag.add_field('dynamics', 'ucomp_temp', time_avg=avgflag)
diag.add_field('dynamics', 'vcomp_sq', time_avg=avgflag)
diag.add_field('dynamics', 'vcomp_omega', time_avg=avgflag)
diag.add_field('dynamics', 'vcomp_temp', time_avg=avgflag)
diag.add_field('dynamics', 'omega_sq', time_avg=avgflag)
diag.add_field('dynamics', 'omega_temp', time_avg=avgflag)
diag.add_field('dynamics', 'temp_sq', time_avg=avgflag)

diag.add_field('atmosphere', 'rh', time_avg=avgflag)
diag.add_field('mixed_layer','albedo', time_avg=avgflag) 

#Step 5. Define the namelist options, which will get passed to the fortran to configure the model.
#Define values for the 'core' namelist
print('Step 5: setting namelist')
namelist = Namelist({
	'astronomy_nml': {
		'obliq'  : 0.
	},

    'constants_nml': {
        'omega': 29.168e-5  # default:7.2921150e-5
    },

    'main_nml': {
        'days'   : runlen,
        'hours'  : 0,
        'minutes': 0,
        'seconds': 0,
        'dt_atmos':720,
        'current_date' : [0001,1,1,0,0,0],
        'calendar' : 'thirty_day'
    },

    'idealized_moist_phys_nml': {
        'do_damping': False,
        'turb':True,  # turb has to be True for mixed_layer to work
        'mixed_layer_bc':True,
        'do_virtual' :False,
        'do_simple': True,
        'roughness_mom':5e-5, #Ocean roughness lengths 5e-4
        'roughness_heat':5e-5,
        'roughness_moist':5e-5,      
        'two_stream_gray': do_grey, 
        'do_rrtm_radiation': do_rrtm,
        'do_socrates_radiation': do_socrates,
        'convection_scheme':'DRY',
    },

    'vert_turb_driver_nml': {
        'do_mellor_yamada': False,     # default: True
        'do_diffusivity': True,        # default: False
        'do_simple': True,             # default: False
        'constant_gust': 0.0,          # default: 1.0
        'use_tau': False
    },

	'dry_convection_nml': {
		'tau':1800., # restoring time scale in [s]
		'gamma':1. # [0,1], 1 corresponds to dry adiabat, 0 corresponds to constant T
	},

    'mixed_layer_nml': {
        'tconst' : 285.,
        'prescribe_initial_dist':True,
        'evaporation':False,    
        'depth':30.,							# used for heat capacity over ocean, if dynamic_heat_capacity = False
		'land_depth': 30.,						# used for heat capacity over land, if dynamic_heat_capacity = False
		'albedo_choice': 1,						# FEEDBACK A, 6 means changing albedo with bucket_depth
        'albedo_value': 0.0,					# albedo value if albedo_option = 1
		'update_albedo_from_ice': False,
        'do_qflux' : False						# Do not use prescribed qflux formula
    },

    'diffusivity_nml': {
        'do_entrain':False,
        'do_simple': True
    },

    'surface_flux_nml': {
        'use_virtual_temp': False,
        'do_simple': True,
        'old_dtaudv': True,
#		'land_evap_prefactor': 0.2,
    },

    'atmosphere_nml': {
        'idealized_moist_model': True
    },
    
    'damping_driver_nml': {
        'do_rayleigh': True,
        'trayfric': -1.,              # neg. value: time in *days*
        'sponge_pbottom':  150., #Setting the lower pressure boundary for the model sponge layer in Pa.
        'do_conserve_energy': True,         
    },
    
# not really used "sat_vapor_pres_nml"
    'sat_vapor_pres_nml': {
		'tcmin_simple': -273.,
		'show_all_bad_values': True,
		'tcmax_simple': 800.,
        'do_simple':True
    },

    # FMS Framework configuration
    #'diag_manager_nml': {
    #    'mix_snapshot_average_fields': True  # time avg fields are labelled with time in middle of window
    #},

    'fms_nml': {
        'domains_stack_size': 600000                        # default: 0
    },

    'fms_io_nml': {
        'threading_write': 'single',                         # default: multi
        'fileset_write': 'single',                           # default: multi
    },

    'spectral_dynamics_nml': {
		'eddy_sponge_coeff'   :  4.e-4, # 1e-4 ~ 10 day
		'zmu_sponge_coeff'    :  4.e-4,
		'zmv_sponge_coeff'    :  4.e-4,
        'damping_order': 4,             
        'water_correction_limit': 0.,
        'reference_sea_level_press':1.0e5,
        'num_levels':40,
        'valid_range_t':[100.,800.],
        'initial_sphum':[0.],
        'vert_coord_option':'uneven_sigma',
        'surf_res':0.2, #Parameter that sets the vertical distribution of sigma levels
        'scale_heights' : 11.0,
        'exponent':7.0,
        'robert_coeff':0.03,
		'omega_coriolis': 7.2921150e-5  # default 7.2921150e-5
    },

	'spectral_init_cond_nml':{
		'initial_temperature': 300.
	},

	'two_stream_gray_rad_nml': {
	},
	'rrtm_radiation_nml':{
	},
	'socrates_rad_nml':{
	}
})

if radiation_scheme == 'grey':
	namelist['two_stream_gray_rad_nml']['rad_scheme'] = 'FRIERSON' # two-band grey atmosphere
	namelist['two_stream_gray_rad_nml']['solar_constant'] = 1360 # insolation
	namelist['two_stream_gray_rad_nml']['do_seasonal'] = True # calculate zenith angle rather than use the idealized latitude profile
	namelist['two_stream_gray_rad_nml']['solday'] = 90 # no seasonal variation
	namelist['two_stream_gray_rad_nml']['linear_tau'] = 1. # set to 1 to omit the longwave absorptivity dependence on altitude
	namelist['two_stream_gray_rad_nml']['ir_tau_eq'] = 4. # longwave optical depth at equator
	namelist['two_stream_gray_rad_nml']['ir_tau_pole'] = 4. # longwave optical depth at poles
	namelist['two_stream_gray_rad_nml']['atm_abs'] = 4. # shortwave optical depth
	namelist['two_stream_gray_rad_nml']['solar_exponent'] = 1. # set to 1 to omit the shortwave absorptivity dependence on altitude
	namelist['two_stream_gray_rad_nml']['sw_diff'] = 0. # no meridional variation of shortwave optical depth
	namelist['two_stream_gray_rad_nml']['use_daily_mean_insol'] = False # use daily mean insolation?
	namelist['two_stream_gray_rad_nml']['use_exo_diurnal'] = True # allow omega to be used in radiation

if radiation_scheme == 'rrtm':
	namelist['rrtm_radiation_nml']['do_read_ozone'] = True
	namelist['rrtm_radiation_nml']['ozone_file'] = ozonefilename
	namelist['rrtm_radiation_nml']['solr_cnst'] = 1360.
	namelist['rrtm_radiation_nml']['dt_rad'] = 1800.

if radiation_scheme == 'socrates': # need to learn namelist for socrates
	namelist['socrates_rad_nml']['stellar_constant'] = 1370.
	namelist['socrates_rad_nml']['lw_spectral_filename'] = os.path.join(GFDL_BASE,'src/atmos_param/socrates/src/trunk/data/spectra/ga7/sp_lw_ga7_dsa')
	namelist['socrates_rad_nml']['sw_spectral_filename'] = os.path.join(GFDL_BASE,'src/atmos_param/socrates/src/trunk/data/spectra/ga7/sp_sw_ga7_dsa_sun')
	namelist['socrates_rad_nml']['do_read_ozone'] = True
	namelist['socrates_rad_nml']['ozone_file_name'] = ozonefilename
	namelist['socrates_rad_nml']['ozone_field_name'] = ozonefilename
	namelist['socrates_rad_nml']['dt_rad'] = 1800
	namelist['socrates_rad_nml']['store_intermediate_rad'] = True
	namelist['socrates_rad_nml']['chunk_size'] = 16
	namelist['socrates_rad_nml']['use_pressure_interp_for_half_levels'] = False
	namelist['socrates_rad_nml']['tidally_locked'] = False

if debug:
	namelist['main_nml']['days'] = 1
	namelist['main_nml']['hours'] = 0


#Step 6. Compile the fortran code
print('Step 6: compile..............')
if 'compile' in todo:
	cb.compile(extra_pathnames=extra_pathnames, executable_name=executable_name, debug=debug)
	os.system("cp -rf "+os.path.realpath(__file__)+" "+builddir+"/run_isca.py"+append)

#Step 7. Soft links and run script
if 'compile' in todo or 'configure' in todo:
	if os.path.exists(builddir+'/run.sub'+append):
		os.remove(builddir+'/run.sub'+append)
	runsub_template=os.environ['GFDL_BASE']+'/src/extra/python/isca/templates/run.sub'
	with open(runsub_template, 'r') as file :
		filedata = file.read()
	NNODES = np.ceil(NCORES/32.).astype(int)
	print("using"+str(NNODES))
	NCPUS = min(NCORES,32)
	filedata = filedata.replace('_CASENAME_', append)
	filedata = filedata.replace('_NNODES_', str(NNODES))
	filedata = filedata.replace('_NCPUS_', str(NCPUS))
	filedata = filedata.replace('_TODO_', 'run')
	filedata = filedata.replace('_ACCOUNT_', 'UHAR0008')
	filedata = filedata.replace('_TOTALRUN_', str(totalrun))
	filedata = filedata.replace('_EACHRUN_', str(eachrun))
	filedata = filedata.replace('_CURRENTRUN_', str(eachrun))
	filedata = filedata.replace('run_isca.py', 'run_isca.py'+append)
	filedata = filedata.replace('run.sub', 'run.sub'+append)

	with open(builddir+'/run.sub'+append, 'w') as file :
		file.write(filedata)

	with open(runsub_template+'_debug', 'r') as file :
		filedata_debug = file.read()
	filedata_debug = filedata_debug.replace('_CASENAME0_', casename0)
	filedata_debug = filedata_debug.replace('_CASENAME_', append)
	filedata_debug = filedata_debug.replace('_OBJDIR_', os.environ['GFDL_WORK']+'/codebase/'+objdir)
	filedata_debug = filedata_debug.replace('_NNODES_', '1')
	filedata_debug = filedata_debug.replace('_NCPUS_', '1')
	filedata_debug = filedata_debug.replace('_TODO_', 'rundebug')
	filedata_debug = filedata_debug.replace('_ACCOUNT_', 'UHAR0008')
	filedata_debug = filedata_debug.replace('run_isca.py', 'run_isca.py'+append)

	if not os.path.exists(builddir+'/debug'):
		os.makedirs(builddir+'/debug')
	with open(builddir+'/debug/run.sub_debug'+append, 'w') as file :
		file.write(filedata_debug)
	

#	shutil.rmtree(builddir, ignore_errors=True)
	if not os.path.exists(os.environ['GFDL_DATA']+casename):
		os.makedirs(os.environ['GFDL_DATA']+casename)
	if not os.path.exists(builddir+'/data'+append):
		os.symlink(os.environ['GFDL_DATA']+casename,builddir+'/data'+append)
	if not os.path.exists(os.environ['GFDL_WORK']+'/experiment/'+casename+'/run'):
		os.makedirs(os.environ['GFDL_WORK']+'/experiment/'+casename+'/run')
	if not os.path.exists(builddir+'/run'+append):
		os.symlink(os.environ['GFDL_WORK']+'/experiment/'+casename+'/run',builddir+'/run'+append)


#Step 8. Cache everything
if 'compile' in todo or 'archive' in todo:
	print('modified:')
	for fi in os.listdir(builddir+'/srcmod'):
		print(fi)
		os.system("cp -rf "+builddir+"/srcmod/"+fi+" "+isca_caseroot+"/srcmod/"+fi+"_"+casename)
	os.system("cp -rf "+os.path.realpath(__file__)+" "+isca_caseroot+"/srcmod/build_isca.py_"+casename)

#Step 9. Run the fortran code
if 'run' in todo:
	exp = Experiment(casename, codebase=cb)
	exp.clear_rundir()

	exp.diag_table = diag
	exp.namelist = namelist.copy()

	if qflux_file_name!='off':
		inputfiles.append(os.path.join(builddir,qflux_file_name+'.nc'))
		exp.namelist['mixed_layer_nml']['load_qflux'] = True
		exp.namelist['mixed_layer_nml']['time_varying_qflux'] = False
		exp.namelist['mixed_layer_nml']['qflux_file_name'] = qflux_file_name
	else:
		print('no qflux')
		exp.namelist['mixed_layer_nml']['load_qflux'] = False

	exp.inputfiles = inputfiles
	exp.set_resolution(*RESOLUTION)
	for i in range(1,601):
		if i == 1:
			exp.run(1, use_restart=False, num_cores=NCORES,run_ddt=run_ddt)
		else:
			exp.run(i, num_cores=NCORES,run_ddt=run_ddt)
