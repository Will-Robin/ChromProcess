Search.setIndex({docnames:["API","Classes","Loading","Processing","Utils","Writers","index"],envversion:{"sphinx.domains.c":2,"sphinx.domains.changeset":1,"sphinx.domains.citation":1,"sphinx.domains.cpp":4,"sphinx.domains.index":1,"sphinx.domains.javascript":2,"sphinx.domains.math":2,"sphinx.domains.python":3,"sphinx.domains.rst":2,"sphinx.domains.std":2,sphinx:56},filenames:["API.rst","Classes.rst","Loading.rst","Processing.rst","Utils.rst","Writers.rst","index.rst"],objects:{"ChromProcess.Classes":[[1,0,1,"","AnalysisInformation"],[1,0,1,"","Chromatogram"],[1,0,1,"","DataReport"],[1,0,1,"","ExperimentConditions"],[1,0,1,"","InstrumentCalibration"],[1,0,1,"","MassSpectrum"],[1,0,1,"","Peak"],[1,0,1,"","PeakCollection"],[1,0,1,"","PeakCollectionSeries"]],"ChromProcess.Classes.AnalysisInformation":[[1,1,1,"","write_to_file"]],"ChromProcess.Classes.Chromatogram":[[1,1,1,"","get_mass_spectrum"],[1,1,1,"","ion_chromatogram"],[1,1,1,"","write_peak_collection"],[1,1,1,"","write_to_csv"],[1,1,1,"","write_to_json"]],"ChromProcess.Classes.DataReport":[[1,1,1,"","find_repeat_data_entries"],[1,1,1,"","remove_entries_below_threshold"],[1,1,1,"","remove_repeat_entries"],[1,1,1,"","remove_specific_entries"],[1,1,1,"","write_to_file"]],"ChromProcess.Classes.ExperimentConditions":[[1,1,1,"","write_to_file"]],"ChromProcess.Classes.InstrumentCalibration":[[1,1,1,"","get_info"],[1,1,1,"","modify_boundaries"]],"ChromProcess.Classes.MassSpectrum":[[1,1,1,"","write_to_string"]],"ChromProcess.Classes.Peak":[[1,1,1,"","apply_linear_calibration"],[1,1,1,"","apply_quadratic_calibration"],[1,1,1,"","assign_peak"],[1,1,1,"","calculate_error"],[1,1,1,"","dilution_correction"],[1,1,1,"","get_height"],[1,1,1,"","get_integral"],[1,1,1,"","get_mass_spectrum"],[1,1,1,"","reference_height_to_IS"],[1,1,1,"","reference_integral_to_IS"]],"ChromProcess.Classes.PeakCollection":[[1,1,1,"","add_mass_spectra"],[1,1,1,"","align_peaks_to_IS"],[1,1,1,"","apply_calibrations_to_peaks"],[1,1,1,"","assign_peaks"],[1,1,1,"","calculate_conc_errors"],[1,1,1,"","dilution_correct_peaks"],[1,1,1,"","get_all_assigned_compounds"],[1,1,1,"","get_peak_positions"],[1,1,1,"","reference_heights_to_IS"],[1,1,1,"","reference_integrals_to_IS"],[1,1,1,"","remove_peaks_below_threshold"],[1,1,1,"","write_to_file"]],"ChromProcess.Classes.PeakCollectionSeries":[[1,1,1,"","align_peaks_to_IS"],[1,1,1,"","apply_calibrations"],[1,1,1,"","apply_peak_dilution_factors"],[1,1,1,"","assign_clusters"],[1,1,1,"","assign_peaks"],[1,1,1,"","calculate_conc_errors"],[1,1,1,"","create_peak_series"],[1,1,1,"","get_all_assigned_compounds"],[1,1,1,"","get_peak_clusters"],[1,1,1,"","get_peak_positions"],[1,1,1,"","reference_heights_to_IS"],[1,1,1,"","reference_integrals_to_IS"],[1,1,1,"","remove_peaks_below_threshold"],[1,1,1,"","series_traces_as_dict"],[1,1,1,"","write_data_reports"]],"ChromProcess.Loading.analysis_info":[[2,2,0,"-","analysis_from_csv"],[2,2,0,"-","analysis_from_csv_legacy"]],"ChromProcess.Loading.analysis_info.analysis_from_csv":[[2,3,1,"","analysis_from_csv"]],"ChromProcess.Loading.analysis_info.analysis_from_csv_legacy":[[2,3,1,"","analysis_from_csv_legacy"]],"ChromProcess.Loading.chromatogram.cdf":[[2,2,0,"-","cdf_loading"],[2,2,0,"-","chrom_from_cdf"],[2,2,0,"-","instruments"]],"ChromProcess.Loading.chromatogram.cdf.cdf_loading":[[2,3,1,"","load_from_cdf"]],"ChromProcess.Loading.chromatogram.cdf.chrom_from_cdf":[[2,3,1,"","chrom_from_cdf"]],"ChromProcess.Loading.chromatogram.cdf.instruments":[[2,0,1,"","JEOL"]],"ChromProcess.Loading.chromatogram.cdf.instruments.JEOL":[[2,4,1,"","ACTUAL_SCAN_NUMBER"],[2,4,1,"","AD_COADDITION_FACTOR"],[2,4,1,"","AD_SAMPLING_RATE"],[2,4,1,"","ERROR_LOG"],[2,4,1,"","FLAG_COUNT"],[2,4,1,"","INSTRUMENT_APP_VERSION"],[2,4,1,"","INSTRUMENT_COMMENTS"],[2,4,1,"","INSTRUMENT_FW_VERSION"],[2,4,1,"","INSTRUMENT_ID"],[2,4,1,"","INSTRUMENT_MFR"],[2,4,1,"","INSTRUMENT_MODEL"],[2,4,1,"","INSTRUMENT_NAME"],[2,4,1,"","INSTRUMENT_OS_VERSION"],[2,4,1,"","INSTRUMENT_SERIAL_NO"],[2,4,1,"","INSTRUMENT_SW_VERSI"],[2,4,1,"","INTER_SCAN_TIME"],[2,4,1,"","MASS_RANGE_MAX"],[2,4,1,"","MASS_RANGE_MIN"],[2,4,1,"","MS_INTENSITY_KEY"],[2,4,1,"","MZ_KEY"],[2,4,1,"","POINT_COUNT_KEY"],[2,4,1,"","RESOLUTION"],[2,4,1,"","SCAN_DURATION"],[2,4,1,"","SCAN_INDEX_KEY"],[2,4,1,"","TIC_KEY"],[2,4,1,"","TIME_CONVERSION"],[2,4,1,"","TIME_KEY"],[2,4,1,"","TIME_RANGE_MAX"],[2,4,1,"","TIME_RANGE_MIN"],[2,4,1,"","X_UNIT"],[2,4,1,"","Y_UNIT"]],"ChromProcess.Loading.chromatogram.ion_chromatogram":[[2,2,0,"-","ion_chromatogram_from_peak"],[2,2,0,"-","ion_chromatogram_from_region"]],"ChromProcess.Loading.chromatogram.ion_chromatogram.ion_chromatogram_from_peak":[[2,3,1,"","ion_chromatogram_from_peak"]],"ChromProcess.Loading.chromatogram.ion_chromatogram.ion_chromatogram_from_region":[[2,3,1,"","ion_chromatogram_from_region"]],"ChromProcess.Loading.chromatogram.text":[[2,2,0,"-","chrom_from_csv"],[2,2,0,"-","chrom_from_labsolutions_ascii"],[2,2,0,"-","chrom_from_text"]],"ChromProcess.Loading.chromatogram.text.chrom_from_csv":[[2,3,1,"","chrom_from_csv"]],"ChromProcess.Loading.chromatogram.text.chrom_from_labsolutions_ascii":[[2,3,1,"","chrom_from_labsolutions_ascii"]],"ChromProcess.Loading.chromatogram.text.chrom_from_text":[[2,3,1,"","chrom_from_text"]],"ChromProcess.Loading.custom":[[2,2,0,"-","local_assignments_from_csv"],[2,2,0,"-","mass_spectra_from_csv"],[2,2,0,"-","point_removals_from_csv"]],"ChromProcess.Loading.custom.local_assignments_from_csv":[[2,3,1,"","local_assignments_from_csv"]],"ChromProcess.Loading.custom.mass_spectra_from_csv":[[2,3,1,"","mass_spectra_from_csv"]],"ChromProcess.Loading.custom.point_removals_from_csv":[[2,3,1,"","read_point_removals_file"]],"ChromProcess.Loading.data_report":[[2,2,0,"-","data_report_from_csv"]],"ChromProcess.Loading.data_report.data_report_from_csv":[[2,3,1,"","data_report_from_csv"]],"ChromProcess.Loading.experiment_conditions":[[2,2,0,"-","conditions_from_csv"]],"ChromProcess.Loading.experiment_conditions.conditions_from_csv":[[2,3,1,"","conditions_from_csv"]],"ChromProcess.Loading.instrument_calibration":[[2,2,0,"-","instrument_cal_from_csv"]],"ChromProcess.Loading.instrument_calibration.instrument_cal_from_csv":[[2,3,1,"","instrument_cal_from_csv"]],"ChromProcess.Loading.mass_spectrum":[[2,2,0,"-","mass_spectrum_from_peak"]],"ChromProcess.Loading.mass_spectrum.mass_spectrum_from_peak":[[2,3,1,"","mass_spectrum_from_peak"]],"ChromProcess.Loading.parsers":[[2,2,0,"-","parsers"]],"ChromProcess.Loading.parsers.parsers":[[2,3,1,"","import_file_section"],[2,3,1,"","parse_text_columns"]],"ChromProcess.Loading.peak":[[2,2,0,"-","peak_from_chromatogram"]],"ChromProcess.Loading.peak.peak_from_chromatogram":[[2,3,1,"","peak_from_chromatogram"]],"ChromProcess.Loading.peak_collection":[[2,2,0,"-","peak_collection_from_csv"]],"ChromProcess.Loading.peak_collection.peak_collection_from_csv":[[2,3,1,"","peak_collection_from_csv"]],"ChromProcess.Processing":[[3,2,0,"-","custom"]],"ChromProcess.Processing.chromatogram":[[3,2,0,"-","background_subtraction"],[3,2,0,"-","find_peaks"],[3,2,0,"-","modify_chromatogram"],[3,2,0,"-","stack_chromatograms"]],"ChromProcess.Processing.chromatogram.background_subtraction":[[3,3,1,"","ic_background_subtraction"]],"ChromProcess.Processing.chromatogram.find_peaks":[[3,3,1,"","find_peaks_in_region"]],"ChromProcess.Processing.chromatogram.modify_chromatogram":[[3,3,1,"","add_peaks_to_chromatogram"],[3,3,1,"","integrate_chromatogram_peaks"],[3,3,1,"","internal_standard_integral"]],"ChromProcess.Processing.chromatogram.stack_chromatograms":[[3,3,1,"","get_chrom_time_min_max"],[3,3,1,"","stack_chromatograms"]],"ChromProcess.Processing.peak":[[3,2,0,"-","assign_peak"]],"ChromProcess.Processing.peak.assign_peak":[[3,3,1,"","assign_retention_time"]],"ChromProcess.Utils":[[4,2,0,"-","custom"]],"ChromProcess.Utils.peak_finding":[[4,2,0,"-","pick_peaks"]],"ChromProcess.Utils.peak_finding.pick_peaks":[[4,3,1,"","find_peak_boundaries"],[4,3,1,"","find_peak_boundaries_look_ahead"],[4,3,1,"","find_peaks"],[4,3,1,"","find_peaks_scipy"]],"ChromProcess.Utils.signal_processing":[[4,2,0,"-","deconvolution"],[4,2,0,"-","signal_processing"]],"ChromProcess.Utils.signal_processing.deconvolution":[[4,3,1,"","deconvolute_peak"],[4,3,1,"","deconvolute_region"],[4,3,1,"","fit_gaussian_peaks"]],"ChromProcess.Utils.signal_processing.signal_processing":[[4,3,1,"","savitzky_golay"]],"ChromProcess.Utils.utils":[[4,2,0,"-","clustering"],[4,2,0,"-","error_propagation"],[4,2,0,"-","functions"],[4,2,0,"-","utils"]],"ChromProcess.Utils.utils.clustering":[[4,3,1,"","cluster"],[4,3,1,"","cluster_indices"]],"ChromProcess.Utils.utils.error_propagation":[[4,3,1,"","mult_div_error_prop"],[4,3,1,"","sum_error_prop"]],"ChromProcess.Utils.utils.functions":[[4,3,1,"","inverse_linear"],[4,3,1,"","inverse_quadratic"],[4,3,1,"","inverse_quadratic_standard_error"],[4,3,1,"","linear"],[4,3,1,"","quadratic_function"],[4,3,1,"","residual_squared_error"]],"ChromProcess.Utils.utils.utils":[[4,3,1,"","bin_dictionary"],[4,3,1,"","indices_from_boundary"],[4,3,1,"","is_float"],[4,3,1,"","is_int"],[4,3,1,"","peak_dict_to_spreadsheet"],[4,3,1,"","peak_indices_to_times"]],"ChromProcess.Writers.chromatogram":[[5,2,0,"-","chromatogram_to_csv"],[5,2,0,"-","chromatogram_to_json"],[5,2,0,"-","chromatogram_to_peak_collection"]],"ChromProcess.Writers.chromatogram.chromatogram_to_csv":[[5,3,1,"","chromatogram_to_csv"],[5,3,1,"","write_chromatogram_csv_text"]],"ChromProcess.Writers.chromatogram.chromatogram_to_json":[[5,3,1,"","chromatogram_to_json"],[5,3,1,"","write_chromatogram_json_text"]],"ChromProcess.Writers.chromatogram.chromatogram_to_peak_collection":[[5,3,1,"","chromatogram_to_peak_collection"],[5,3,1,"","write_peak_collection_text"]],"ChromProcess.Writers.data_report":[[5,2,0,"-","data_report_to_csv"]],"ChromProcess.Writers.data_report.data_report_to_csv":[[5,3,1,"","data_report_to_csv"],[5,3,1,"","write_data_report_text"]],"ChromProcess.Writers.general":[[5,2,0,"-","write_header"]],"ChromProcess.Writers.general.write_header":[[5,3,1,"","write_conditions_header"]],"ChromProcess.Writers.mass_spectra":[[5,2,0,"-","mass_spectrum_to_string"]],"ChromProcess.Writers.mass_spectra.mass_spectrum_to_string":[[5,3,1,"","mass_spectrum_to_string_cols"],[5,3,1,"","mass_spectrum_to_string_rows"]],"ChromProcess.Writers.peak":[[5,2,0,"-","peak_to_entry_text"]],"ChromProcess.Writers.peak.peak_to_entry_text":[[5,3,1,"","peak_to_entry_text"]],"ChromProcess.Writers.peak_collection":[[5,2,0,"-","peak_collection_to_csv"]],"ChromProcess.Writers.peak_collection.peak_collection_to_csv":[[5,3,1,"","peak_collection_to_csv"]],"ChromProcess.Writers.peak_collection_series":[[5,2,0,"-","peak_collection_series_to_data_report"]],"ChromProcess.Writers.peak_collection_series.peak_collection_series_to_data_report":[[5,3,1,"","peak_collection_series_to_data_report"]]},objnames:{"0":["py","class","Python class"],"1":["py","method","Python method"],"2":["py","module","Python module"],"3":["py","function","Python function"],"4":["py","attribute","Python attribute"]},objtypes:{"0":"py:class","1":"py:method","2":"py:module","3":"py:function","4":"py:attribute"},terms:{"0":[1,2,3,4],"001":4,"005":4,"025":4,"05":4,"0th":4,"1":[1,2,3,4],"100":4,"10000":4,"1001":4,"12":4,"13":4,"1627":4,"1639":4,"1964":4,"1d":[1,3,4],"1e":4,"1st":4,"2":4,"2d":[1,3],"3":[2,3],"31":4,"36":4,"3rd":4,"4":4,"5":4,"500":[3,4],"60":2,"7":4,"8":4,"9780521880688":4,"class":[0,2,3,5,6],"default":[2,4],"do":[1,2,3],"export":2,"float":[1,2,3,4],"function":[1,2,3,4],"import":4,"int":[2,4],"new":4,"return":[1,2,3,4,5],"true":4,A:[1,2,4],By:4,For:[2,5],IS:4,If:[1,2,3,4],It:4,The:[1,3,4],To:4,_1gaussian:4,a_d_coaddition_factor:2,a_d_sampling_r:2,about:1,abov:4,accord:4,account:1,actual_scan_numb:2,ad:1,ad_coaddition_factor:2,ad_sampling_r:2,adapt:4,add:[1,3],add_mass_spectra:1,add_peaks_to_chromatogram:3,advantag:4,agglom:4,agglomer:1,ahead:4,algorithm:[1,4],align:1,align_peaks_to_i:1,all:[1,2,3,4],allow:4,alwai:4,amplitud:4,an:[1,2,3,4],analysi:[1,2,3],analysis_from_csv:2,analysis_from_csv_legaci:2,analysis_info:2,analysis_inform:[1,5],analysisinform:[1,2,5],analyt:4,andi:2,ani:4,apex:[1,2],api:6,appli:1,apply_calibr:1,apply_calibrations_to_peak:1,apply_linear_calibr:1,apply_peak_dilution_factor:1,apply_quadratic_calibr:1,approach:4,ar:[1,2,4],arbitrai:1,argument:4,arrai:[1,2,3,4],array_lik:4,art:4,assign:[1,2,3],assign_clust:1,assign_peak:[1,3],assign_retention_tim:3,assigned_compound:1,associ:3,atjacob:4,attribut:1,averag:4,axi:4,b:[1,4],background_subtract:3,base:[1,2,3,4,5],baselin:[1,3],baseline_subtract:[1,3],been:1,befor:[1,3],behind:4,below:[1,3,4],better:4,between:[1,2,3,4],beyond:4,bin:[1,4],bin_dictionari:4,bool:[1,3,4],bound:[1,2,4],boundari:[1,2,3,4],c:[1,4],c_set:2,calc:4,calcul:[1,4],calculate_conc_error:1,calculate_error:1,calib:1,calibr:[1,2,4],cambridg:4,can:[1,4],cannot:4,caus:4,cdf:2,cdf_load:2,center:4,centr:4,certain:1,ch1:2,chemistri:4,chrom:2,chrom_from_cdf:2,chrom_from_csv:2,chrom_from_labsolutions_ascii:2,chrom_from_text:2,chrom_stack:3,chromatogram:[1,2,3,4,5],chromatogram_seri:3,chromatogram_to_csv:5,chromatogram_to_json:5,chromatogram_to_peak_collect:5,chromatograph:1,chromatographi:1,chromprocess:[1,2,3,4,5],cluster:[1,4],cluster_indic:4,code:2,collect:[1,5],column:[2,5],com:4,combin:[1,4],comma:2,compound:[1,2],compound_1_nam:2,compound_2_nam:2,compound_nam:1,comput:4,conc_dict:1,concentr:[1,4],concept:4,condit:[1,2,5],conditions_from_csv:2,consid:1,constraint:4,contain:[1,2,3,4],convers:1,convert:[2,4,5],cookbook:4,coordint:2,correct:1,correl:4,correspond:4,could:4,count:[1,2],covari:4,creat:[1,2,3,4,5],create_peak_seri:1,csv:[1,2,5],csv_string:5,current:[3,4],curv:4,custom:2,cut:4,data:[1,2,3,4,5],data_1:4,data_2:4,data_contain:2,data_kei:2,data_report:[2,5],data_report_from_csv:2,data_report_text:5,data_report_to_csv:5,datareport:[1,2,5],dataset:2,decis:3,deconvolut:4,deconvolute_peak:4,deconvolute_region:4,decreas:4,defin:[1,2],delimit:2,deriv:[1,2,4],detail:5,detect:4,detector:2,determin:4,dict:[1,2,3,4],dictionari:[1,2,3,4],diff:4,differ:4,differenti:4,dil:4,dilut:[1,4],dilution_correct:1,dilution_correct_peak:1,dilution_factor:1,directori:[1,5],distanc:4,divid:1,e:4,each:4,earli:4,edit:4,either:4,element:4,end:[1,2,3,4],end_token:2,enough:4,entri:[1,5],equal:4,equat:4,err_dict:1,error:[1,4],error_log:2,error_propag:4,estim:[1,4],etc:2,ewmpti:2,exampl:4,exce:[1,2,3],exp:4,expect:2,experi:[2,5],experiment_condit:2,experimentcondit:[1,2,5],extend:4,extra:1,extract:2,f:4,factor:[1,4],factor_error:1,fall:1,fals:[1,2,3],far:4,featur:4,file:[1,2,5],filenam:[1,2,5],filter:4,find:[1,2,3,4],find_peak:[3,4],find_peak_boundari:4,find_peak_boundaries_look_ahead:4,find_peaks_in_region:3,find_peaks_scipi:4,find_repeat_data_entri:1,first:4,fit:4,fit_gaussian_peak:4,flag_count:2,flanneri:4,fname:2,folder:[1,5],format:[1,2,5],former:4,found:[3,4],fraction:[2,3],frequenc:4,from:[1,2,3,4,5],g:4,gaussian:4,gcm:2,gener:[2,5],get:[1,2,3,4],get_all_assigned_compound:1,get_chrom_time_min_max:3,get_height:1,get_info:1,get_integr:1,get_mass_spectrum:1,get_peak_clust:1,get_peak_posit:1,github:4,given:1,golai:4,gradient:4,grid:4,guess:4,h:4,ha:[3,4],have:[1,4],header:[4,5],header_text:[1,5],heat:3,height:1,held:1,high:4,higher:4,highest:[3,4],histori:4,howev:3,html:4,http:4,ic_background_subtract:3,idea:4,identifi:1,implement:[3,4],import_file_sect:2,includ:[1,2,5],independ:4,index:[4,6],indic:[1,2,4],indices_from_boundari:4,info:3,inform:[1,2,3,5],initi:4,initial_guess:4,inpend:4,insert:2,instrument:2,instrument_app_vers:2,instrument_cal_from_csv:2,instrument_calibr:[1,2],instrument_com:2,instrument_fw_vers:2,instrument_id:2,instrument_mfr:2,instrument_model:2,instrument_nam:2,instrument_os_vers:2,instrument_serial_no:2,instrument_sw_versi:2,instrumentcalibr:[1,2],integ:4,integr:[1,3,4,5],integral_dict:1,integrate_chromatogram_peak:3,inten:1,intens:[1,2,3],intensity_valu:2,inter:3,inter_scan_tim:2,intercept:4,intern:[1,3,4],internal_standard:1,internal_standard_integr:3,interpol:1,interpret:4,inverse_linear:4,inverse_quadrat:4,inverse_quadratic_standard_error:4,io:4,ion:[1,2,3],ion_chromatogram:[1,2],ion_chromatogram_from_peak:2,ion_chromatogram_from_region:2,is_conc:1,is_conc_err:1,is_end:3,is_float:4,is_height:1,is_int:4,is_integr:1,is_set:1,is_start:3,isbn:4,issu:3,item:[2,4],iter:3,its:[1,2,4],j:4,jeol:2,json:5,k:4,kei:[2,4],kind:4,label:4,labsolut:2,least:4,legend:4,len:3,length:4,less:4,librari:2,like:4,line:2,linear:[1,4],linearli:1,linspac:4,list:[1,2,3,4],load:[0,6],load_from_cdf:2,load_m:2,local:3,local_assignments_from_csv:2,longer:4,look:4,look_ahead:4,low:4,lower:[1,2,4],lowerbound:4,lowest:3,lw:4,m:[1,2,3,4],m_z:1,macro:4,magnitud:1,main:4,make:[1,4],manipul:2,map:3,mass:[1,2,3,5],mass_range_max:2,mass_range_min:2,mass_spectra:[2,5],mass_spectra_from_csv:2,mass_spectrum:[1,2,5],mass_spectrum_from_peak:2,mass_spectrum_to_str:5,mass_spectrum_to_string_col:5,mass_spectrum_to_string_row:5,mass_valu:2,massspectrum:[1,2,5],master:4,match:4,matplotlib:4,matrix:4,max_inten:4,max_tim:3,maxim:[1,4],maximum:2,mean:[1,4],measur:4,method:[1,3],might:4,min:2,min_dist:4,min_inten:4,min_tim:3,minim:4,minimum:4,modif:1,modifi:[1,3,4],modified_bound:[1,2],modify_boundari:1,modify_chromatogram:3,modul:[1,6],more:[1,4],move:4,ms:3,ms_intensity_kei:2,ms_list:1,ms_string:[1,5],mult_div_error_prop:4,multipl:4,multipli:1,must:[1,4],mz:[1,2],mz_kei:2,n:4,name:[1,2,3,5],ndarrai:4,neg:4,netcdf4:2,nois:4,noisi:4,none:[1,3,4,5],normal:4,normalis:1,note:[1,4],np:4,num_peak:4,number:4,numer:4,numpi:[1,2,3,4],object:[1,2,3,4,5],obtain:1,odd:4,off:4,older:2,omit:2,ommit:3,one:[1,4],onli:4,oper:[1,4],option:4,order:[1,3,4],ordinal_delimit:2,organis:2,origin:[3,4],other:4,out_log:4,output:[1,4,5],over:[2,4],p:4,page:6,pair:2,paramet:[1,2,3,4,5],parent:[1,2],parent_chromatogram:2,pars:2,parse_text_column:2,parser:2,particularli:4,pass:[1,4],path:[1,2,5],pathlib:[1,2,5],pcov:4,peak:[1,2,3,4,5],peak_collect:[1,2,5],peak_collection_from_csv:2,peak_collection_seri:5,peak_collection_series_to_data_report:5,peak_collection_str:5,peak_collection_to_csv:5,peak_dict:4,peak_dict_to_spreadsheet:4,peak_end:4,peak_featur:4,peak_find:4,peak_from_chromatogram:2,peak_grid_transpos:4,peak_head:4,peak_indic:4,peak_indices_to_tim:4,peak_nam:3,peak_po:1,peak_start:4,peak_to_entry_text:5,peakcollect:[1,2,5],peakcollectionseri:[1,5],peaks_indic:4,peakutil:4,peal:1,perform:3,pick:[3,4],pick_peak:4,picked_peak:4,place:[1,3,5],plot:[3,4],plt:4,po:1,point:[1,2,4],point_count:2,point_count_kei:2,point_delimit:2,point_remov:2,point_removals_from_csv:2,polynomi:4,popt:4,posit:[1,4],possibl:[3,4],potenti:3,pp:4,pre:1,prefer:4,present:[2,3],preserv:4,press:4,prioriti:3,probabl:1,procedur:[1,4],process:[0,6],promin:4,prominc:4,propag:4,pyplot:4,python:[1,2,3,4,5],quadrat:[1,4],quadratic_funct:4,quit:4,r:4,random:4,rang:4,rate:4,rather:4,raw:1,re:3,read:2,read_point_removals_fil:2,readthedoc:4,recip:4,recommend:4,reconsitut:3,reconstitut:3,reduc:4,refactor:[1,4],refer:4,reference_height_to_i:1,reference_heights_to_i:1,reference_integral_to_i:1,reference_integrals_to_i:1,regex:2,region:[2,3,4],region_peak_pick:4,rel:2,reli:[3,4],remain:2,remov:[1,2,3,4],remove_entries_below_threshold:1,remove_list:1,remove_peaks_below_threshold:1,remove_repeat_entri:1,remove_specific_entri:1,repeat:1,repeat_entri:1,report:[1,2,5],requir:[2,4],residu:4,residual_squared_error:4,resolut:2,retent:[1,2,3,4],retention_tim:[1,3],roughli:4,round_digit:2,routin:4,row:[1,5],rt:4,run:4,s:[1,2,4,5],sa2:4,sab:4,sac:4,sampl:1,satisfi:4,save:1,savitzki:4,savitzky_golai:4,savitzkygolai:4,sb2:4,sbc:4,sc2:4,scan_acquisition_tim:2,scan_dur:2,scan_index:2,scan_index_kei:2,scientif:4,scipi:4,scope:4,scrape:2,scraper:[],search:[4,6],second:4,section:2,secur:1,see:[1,4],self:1,separ:2,sequenc:4,seri:[1,3,4,5],series_assigned_compound:1,series_traces_as_dict:1,series_unit:4,series_valu:4,set:[1,2,3,4],shape:4,shimadzu:2,should:[1,4],show:4,sig:4,signal:[1,2,3,4,5],signal_process:4,similar:[1,4],simplifi:4,singl:1,size:4,slope:4,small:4,smooth:4,softwar:2,solut:4,some:[3,4],sort:[1,4],sort_series_ord:4,space:4,specif:2,specifi:[1,4],spectra:[1,2,3],spectrim:1,spectrum:[1,2,5],spectrum_filt:2,spreadsheet:4,sqrt:4,squar:4,stack:3,stack_chromatogram:3,standard:[1,3,4],start:[1,2,3,4],start_token:2,stdev:4,stepwis:4,storag:[1,5],store:1,str:[1,2,3,4,5],string:[2,5],structur:1,substract:1,subtract:3,suit:4,sum:4,sum_error_prop:4,suppli:4,sy2:4,t:4,tabl:5,take:[3,4],techniqu:4,term:4,test:4,teukolski:4,text:[1,2,5],th:4,than:[1,4],thei:1,them:[1,3],thi:[1,2,3,4],thing:4,thre:4,threshold:[1,2,3,4],through:4,tic_kei:2,time:[1,2,3,4,5],time_axi:3,time_convers:2,time_kei:2,time_range_max:2,time_range_min:2,todo:[1,3,4],togeth:1,token:2,top:[1,5],total:[1,2,3],total_intens:2,tree:[3,4],two:[2,4],txt:2,type:[2,4],unchang:2,under:4,unew:4,unit:4,univers:4,until:4,upon:3,upper:[1,2,4],upperbound:4,us:[1,2,4],util:[0,6],valu:[1,2,4],value_dict:4,variabl:4,varianc:4,version:[2,5],vetterl:4,w:4,warrant:4,wast:3,were:4,when:4,where:4,whether:[3,4],which:[1,2,3,4],whose:1,width:4,window:4,window_s:4,wise:[1,5],within:[2,3,4],without:4,wlen:4,work:3,write:[1,5],write_chromatogram_csv_text:5,write_chromatogram_json_text:5,write_conditions_head:5,write_data_report:1,write_data_report_text:5,write_head:5,write_peak_collect:1,write_peak_collection_text:5,write_peak_mass_spectra:2,write_to_csv:1,write_to_fil:1,write_to_json:1,write_to_str:1,writer:[0,1,6],written:[2,4],x:[2,4],x_unit:2,x_valu:2,y:[2,4],y_unit:2,y_valu:2,yhat:4,yield:4,ys:4,ysg:4,z:[1,2,3]},titles:["API","Classes","Loading","Processing","Utils","Writers","Documentation for ChromProcess"],titleterms:{"class":1,api:0,chromprocess:6,content:[0,6],document:6,indic:6,load:2,process:3,tabl:6,util:4,writer:5}})