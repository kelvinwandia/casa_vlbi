# casa_pipescript.py

__rethrow_casa_exceptions = True
context = h_init()
context.set_state('ProjectSummary', 'observatory', 'Karl G. Jansky Very Large Array')
context.set_state('ProjectSummary', 'telescope', 'EVLA')
context.set_state('ProjectStructure', 'recipe_name', 'hifv_cal')
try:
    hifv_importdata(vis=['/raid1/scratch/kelvinw/k2_18b/c_band_x_band/23B-307/observation.60297.645682870374/23B-307.sb44672012.eb44930293.60297.64567559028'],
        # vis=['/raid1/scratch/kelvinw/k2_18b/uncalibrated_vla/23B-307/observation.60230.613748263895/23B-307.sb44594812.eb44691528.60230.613198356485'],
        # ['/raid1/scratch/kelvinw/k2_18b/uncalibrated_vla/23B-307/observation.60239.58857523148/23B-307.sb44594812.eb44725045.60239.588568113424'], 
    session=['default'])
    hifv_hanning(pipelinemode="automatic")
    hifv_flagdata(hm_tbuff='1.5int', intents='*POINTING*,*FOCUS*,*ATMOSPHERE*,*SIDEBAND_RATIO*, \
         *UNKNOWN*, *SYSTEM_CONFIGURATION*, *UNSPECIFIED#UNSPECIFIED*',template=True, \
        #  filetemplate='/raid1/scratch/kelvinw/casa_vlbi/selfcal/vla_flagging_template/s_band_d_config.txt'
        )
    hifv_vlasetjy(pipelinemode="automatic")
    hifv_priorcals(pipelinemode="automatic",apply_tec_correction=True)
    hifv_syspower()
    hifv_testBPdcals(weakbp=False,refantignore='ea01,ea02,ea03,ea04,ea05,ea06,ea07,ea08,ea09,ea10,ea11,ea12,ea13,ea14,ea15,ea16,ea17,ea19,ea20,ea21,ea22,ea23,ea24,ea26,ea27,ea28')
    hifv_checkflag(checkflagmode='bpd-vla')
    hifv_semiFinalBPdcals(weakbp=False,refantignore='ea01,ea02,ea03,ea04,ea05,ea06,ea07,ea08,ea09,ea10,ea11,ea12,ea13,ea14,ea15,ea16,ea17,ea19,ea20,ea21,ea22,ea23,ea24,ea26,ea27,ea28' )
    hifv_checkflag(checkflagmode='allcals-vla')
    hifv_semiFinalBPdcals(weakbp=False,refantignore='ea01,ea02,ea03,ea04,ea05,ea06,ea07,ea08,ea09,ea10,ea11,ea12,ea13,ea14,ea15,ea16,ea17,ea19,ea20,ea21,ea22,ea23,ea24,ea26,ea27,ea28' )
    hifv_solint(refantignore='ea01,ea02,ea03,ea04,ea05,ea06,ea07,ea08,ea09,ea10,ea11,ea12,ea13,ea14,ea15,ea16,ea17,ea19,ea20,ea21,ea22,ea23,ea24,ea26,ea27,ea28')
    hifv_fluxboot(fitorder=2,refantignore='ea01,ea02,ea03,ea04,ea05,ea06,ea07,ea08,ea09,ea10,ea11,ea12,ea13,ea14,ea15,ea16,ea17,ea19,ea20,ea21,ea22,ea23,ea24,ea26,ea27,ea28')
    hifv_finalcals(refantignore='ea01,ea02,ea03,ea04,ea05,ea06,ea07,ea08,ea09,ea10,ea11,ea12,ea13,ea14,ea15,ea16,ea17,ea19,ea20,ea21,ea22,ea23,ea24,ea26,ea27,ea28')
    hifv_applycals(pipelinemode="automatic")
    hifv_checkflag(checkflagmode='target-vla')
    hifv_statwt(pipelinemode="automatic")
    hifv_plotsummary(pipelinemode="automatic")
    hif_makeimlist(intent='PHASE,BANDPASS,TARGET', specmode='cont')
    hif_makeimages(hm_masking='none')
    # hifv_mstransform()
    ## Self calibration
    # hif_checkproductsize(maximsize=640)
    # hif_makeimlist(specmode='cont', datatype='regcal')
    # hif_makeimages(hm_cyclefactor=3.0)
    # hif_selfcal()
    # hif_makeimlist(specmode='cont', datatype='selfcal')
    # hif_makeimages(hm_cyclefactor=3.0)
finally:
    h_save()