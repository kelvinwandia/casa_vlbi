context = h_init()
context.set_state('ProjectSummary', 'observatory', 'Karl G. Jansky Very Large Array')
context.set_state('ProjectSummary', 'telescope', 'EVLA')
try:
    hifv_importdata(vis=['/raid1/scratch/kelvinw/k2_18b/official_pipe_cal/s_band_d_config/23B-307.sb44594812.eb44725045.60239.588568113424/K2-18.ms'])
    hifv_flagdata(intents='*POINTING*,*FOCUS*,*ATMOSPHERE*,*SIDEBAND_RATIO*,*UNKNOWN*, *SYSTEM_CONFIGURATION*, *UNSPECIFIED#UNSPECIFIED*',\
        quack=False, autocorr=False, baseband=False, edgespw=False, clip=False, online=False, shadow=False, scan=True,
        # template=True, filetemplate='/raid1/scratch/kelvinw/casa_vlbi/selfcal/vla_flagging_template/s_band_d_config.txt',
        )
    hif_mstransform(pipelinemode="automatic")
    hif_checkproductsize(maximsize=640)
    hif_makeimlist(specmode='cont')
    hif_makeimages(hm_cyclefactor=3.0)
    hif_selfcal()
    hif_makeimlist(specmode='cont',datatype='selfcal')
    hifv_pbcor(pipelinemode="automatic")
    #hifv_exportdata(imaging_products_only=True)
finally:
    h_save()