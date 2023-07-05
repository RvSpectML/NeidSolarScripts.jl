using ArgParse

function parse_commandline_make_manifest_solar()
     s = ArgParseSettings( description = "Make manifest.csv from FITS files in inputdir.")
     @add_arg_table! s begin
         "input"
            help = "Directory (or parent directory) with downloaded FITS files."
            arg_type = String
            required = true
         "output"
            help = "Path for outputs: manifest.csv"
            arg_type = String
            required = true
         "--subdir"
            help = "Subdirectory with downloaded FITS files."
            arg_type = String
            default = ""
#=
         "--pyrohelio"
            help = "Path for pyroheliometer inputs: *.tel"
            arg_type = String
            default = "/gpfs/group/ebf11/default/ebf11/neid_solar/data/pyrohelio/202110_fromChad/"
=#
         "--pyrheliometer"
            help = "Path with pyrheliometer summary statistics: pyrheliometer.csv"
            arg_type = String
         "--root"
            help = "Path to root of data directories"
            arg_type = String
            default = ""
         "--create-continuum"
            help = "Create missing continuum files"
            action = :store_true
         "--calc_order_snr"
            help = "Calculate order SNRs"
            action = :store_true
         "--verbose"
            help = "Verbose outputs"
            action = :store_true
      end

  return parse_args(s)
end


function parse_commandline_make_pyrheliometer_daily()
     s = ArgParseSettings( description = "Make pyrheliometer.csv for specified day.")
     @add_arg_table! s begin
         "manifest_or_date"
            help = "CSV Filename with fields l?filename, object and exptime to get pyrheliometer data for OR Date in YYYYMMDD format to extract pyrheliometer data."
            arg_type = String
            #default = "meta_test.csv" # joinpath(pwd(),"pyrohelio.csv")
            #default = "20210505"
            #default = "20210605"
         "--output"
            help = "Path for outputs: pyrheliometer.csv"
            arg_type = String
            default = joinpath(pwd(),"pyrohelio.csv")
         "--pyrheliometer_dir"
            help = "Path for pyroheliometer inputs: *.tel"
            arg_type = String
            #default = "/gpfs/group/ebf11/default/ebf11/neid_solar/data/pyrohelio/202110_fromChad/"
            default = "/mnt/data_simons/NEID/pyrohelio/"
         "--user"
            help = "NExScI Archive username"
            arg_type = String
         "--password"
            help = "NExScI Archive password"
            arg_type = String
         "--nexsci_login_filename"
            help = "TOML file with NExScI username and password"
            arg_type = String
            default = "nexsci_id.toml"
         "--work_dir"
            help = "Working directory for NExScI query and L0 downloads."
            arg_type = String
            default = "L0"
         "--cookie"
            help = "Filename for NExScI cookie."
            arg_type = String
            default = "neidadmincookie.txt"
         "--query_filename"
            help = "Filename for L0 query results."
            arg_type = String
            default = "meta_l0.csv"
         "--root"
            help = "Path to root of data directories"
            arg_type = String
            default = ""
         "--try_tel_first"
            help = "Try pyrheliometer telemetry file before L0s"
            action = :store_true
         "--verbose"
            help = "Verbose outputs"
            action = :store_true
      end

  return parse_args(s)
end


  function parse_commandline_verify_downloads()
     s = ArgParseSettings( description = "Verify FITS file downloads match contents of meta.csv.")
     @add_arg_table! s begin
         "inputdir"
            help = "Directory with downloaded FITS files."
            arg_type = String
            required = true
         "--object"
            help = "Filter for object name matching (default: empty, so no checking)"
            arg_type = String
            default = ""
#=
         "--root"
             help = "Path to root of data directories"
             arg_type = String
             default = "/storage/group/ebf11/default/pipeline/neid_solar/data"
         "--input"
             help = "Path to inputs FITS files (fits)"
             arg_type = String
             default = "v1.1/L2/"
         "--output"
             help = "Path to daily CCFs (jld2)"
             arg_type = String
             default = "output"
=#
         "--checksums"
            help = "Compute md5 checksums for input FITS files."
            action = :store_true
         "--max_num_spectra"
            help = "Maximum number of spectra to process."
            arg_type = Int
            default = 300
         "--crawl"
            help = "Run on each date's subdirectory."
            action = :store_true
         "--quiet"
            help = "Don't write output to terminal"
            action = :store_true
         "--overrwrite"
            help = "Specify it's ok to overwrite files."
            action = :store_true
      end

     return parse_args(s)
 end

  function parse_commandline_calc_order_ccfs() 
     s = ArgParseSettings( description = "Calculate order CCFs from NEID L2 FITS files.")
     #import_settings!(s, s_files_only, args_only=false)
     add_arg_group!(s, "Files", :argg_files)
     @add_arg_table! s begin
         "manifest"
             help = "Manifest file (CVS) containing FITS files to analyze.\nExpects columns named Filename, bjd, target, and used for continuum continuum_filename."
             arg_type = String
             default = "manifest.csv"
             #required = true
         "output"
             help = "Filename for output CCFs (jld2)"
             arg_type = String
             default = "daily_ccfs.jld2"
         "--param"
             help = "Parameters file"
             arg_type = String
         "--overwrite"
            help = "Specify it's ok to overwrite the output file."
            action = :store_true
         "--verbose"
            help = "Verbose logging."
            action = :store_true
      end
      add_arg_group!(s, "CCF parameters", :argg_ccf_param)
      @add_arg_table! s begin
        "--orders_to_use"
           help = "First and last order _index_ to compute CCFs for"
           nargs = 2
           arg_type = Int64
           #default = [ min_order(NEID2D()), max_order(NEID2D()) ]
           #default = [ first(orders_to_use_default(NEID2D())), last(orders_to_use_default(NEID2D())) ]
        "--mask_scale_factor"
            help = "Specify CCF mask width scale as multiple of NEID default v width (620.953 m/s)"
            arg_type = Float64
            #default = round(lsf_width/default_ccf_mask_v_width(NEID2D()), sigdigits=3)
        "--line_width_50_default"
            help = "Specify default line full width half maximum"
            arg_type = Float64
            #default = 7.9e3 # m/2
        "--variable_mask_scale"
            help = "Vary width of mask to compensate for solar rotation."
            action = :store_true
        "--ccf_mid_vel"
            help = "Middle of velocity range to use for CCFs."
            arg_type = Float64
            #default = 0.0
        "--range_no_mask_change"
            help = "Avoid lines that would result in a mask change due to BC using line width times this factor (TODO: Verify definition)"
            arg_type = Float64
            #default = 5.0
        "--v_step"
            help = "Specify v step for computing CCF "
            arg_type = Float64
            #default = 155.0 # m/2
      end
      add_arg_group!(s, "Line list parameters", :argg_line_list_param)
      @add_arg_table! s begin
         "--line_list_filename"
             help = "Line list filename (input)"
             arg_type = String
             #default = joinpath(pkgdir(NeidSolarScripts),"data","solar_line_list_espresso.csv")
         "--line_list_output_filename"
             help = "Line list filename (output)"
             arg_type = String
         "--recompute_line_weights"
             help = "Force recalculation of SNR-based line weight factors."
             action = :store_true
      end
      add_arg_group!(s, "Continuum normalization parameters", :argg_continuum_param)
      @add_arg_table! s begin
         "--sed_filename"
             help = "Filename for input SED to normalize by"
             arg_type = String
             #default = "/home/eford/Code/RvSpectMLEcoSystem/NeidSolarScripts/data/neidMaster_HR_SmoothLampSED_20210101.fits"
         "--continuum_poly_half_width"
             help = "Half width for window to use for smoothing prior to polynomial fit to continuum."
             arg_type = Int64
             default = 50
         "--quantile_fit_continuum"
             help = "Quantile of rolling window to use prior to polynomial fit to continuum."
             arg_type = Float64
             default = 0.9
         "--order_poly_continuum"
             help = "Order of polynomial to fit after dividing by SED proivded."
             arg_type = Int64
             default = 4
        "--anchors_filename"
             help = "Filename for anchor locations to use in computing continuum to normalize by."
             arg_type = String
             #default = "/home/eford/Code/RvSpectMLEcoSystem/NeidSolarScripts/data/neidMaster_HR_SmoothLampSED_20210101.fits"
        "--anchors_filename_output"
             help = "Filename to write anchor locations to use in computing continuum to normalize by."
             arg_type = String
         "--smoothing_half_width"
             help = "Half width for window to use for smoothing prior to findind local maximum."
             arg_type = Int64
             default = 6
         "--stretch_factor"
             help = "Stretch factor to use for scaling flux in distance calculation for rolling pin continuum normalization"
             arg_type = Float64
             default = 5.0
         "--merging_threshold"
             help = "Threshold (in λ units) for merging nearby continuum anchors."
             arg_type = Float64
             default = 0.25
         "--fwhm_continuum"
            help = "Full width half-max to use for performing continuum normalization."
            arg_type = Float64
            default = Continuum.fwhm_sol/1000
         "--min_rollingpin_r"
            help = "Minimum rolling pin radius for continuum normalization in uints of FWHM."
            arg_type = Float64
            default = 100.0
         "--nu_continuum"
            help = "Exponent for rolling pin radius in continuum normalization"
            arg_type = Float64
            default = 1.0
         "--apply_continuum_normalization"
            help = "Apply continuum normalization."
            action = :store_true
        "--continuum_normalization_individually"
            help = "Calculate continuum normalization for each file rather than averaged spectra."
            action = :store_true
        "--save_continuum_normalization"
            help = "NOT IMPLEMENTED.  Save continuum normalization to file."
            action = :store_true
        "--recompute_continuum_normalization"
            help = "CURRENTLY ONLY OPTION.  Recompute continuum normalization even if continuum file exists."
            action = :store_true
      end
      add_arg_group!(s, "Filter manifest for ", :argg_filter_param)
      @add_arg_table! s begin
         "--target"
            help = "Target field"
            arg_type = String
            default = "Sun"
         "--datestr"
            help = "Filenames that contain date string."
            arg_type = String
        "--max_solar_hour_angle"
           help = "Maximum absolute value of solar hour angle."
           arg_type = Float64
           default = 3.12
        "--max_solar_hour_angle_clean"
           help = "Maximum absolute value of solar hour angle for clean spectrum."
           arg_type = Float64
           default = 2.0
         "--max_airmass"
            help = "Maximum airmass."
            arg_type = Float64
         "--max_airmass_clean"
            help = "Maximum airmass for clean spectrum."
            arg_type = Float64
            default = 2.0
         "--min_expmeter"
            help = "Minimum mean exposure meter flux."
            arg_type = Float64
            default = 6e4
         "--min_expmeter_clean"
            help = "Minimum mean exposure meter flux for clean spectrum."
            arg_type = Float64
            default = 1e5
         "--min_pyrhelio"
            help = "Minimum mean pyrheliometer flux."
            arg_type = Float64
            default = 10^2.95
         "--min_pyrhelio_clean"
            help = "Minimum mean pyrheliometer flux for clean spectrum."
            arg_type = Float64
            default = 10^2.95
         "--max_expmeter_rms_frac"
            help = "Maximum fractional RMS exposure meter flux."
            arg_type = Float64
            default = 0.003
         "--max_expmeter_rms_frac_clean"
            help = "Maximum fractional RMS exposure meter flux for clean spectrum."
            arg_type = Float64
            default = 0.003
         "--max_pyrhelio_rms_frac"
            help = "Maximum fractional RMS pyrheliometer flux."
            arg_type = Float64
            default = 0.0035
         "--max_pyrhelio_rms_frac_clean"
            help = "Maximum fractional RMS pyrheliometer flux for clean spectrum."
            arg_type = Float64
            default = 0.0035
         "--min_expmeter_to_pyrhelio"
            help = "Minimum mean exposure meter flux relative to pyrheliometer flux."
            arg_type = Float64
            default = 0.0
        "--min_expmeter_to_pyrhelio_clean"
            help = "Minimum mean exposure meter flux relative to pyrheliometer flux for clean spectrum."
            arg_type = Float64
            default = 0.0 # 150 from DRP v1.1
#=
         "--min_snr_factor"
            help = "Minimum SNR relative to max_snr."
            arg_type = Float64
         "--min_snr_factor_clean"
            help = "Minimum SNR relative to max_snr for clean spectrum."
            arg_type = Float64
            default = 0.5
         "--max_snr"
            help = "Specify max_snr manually."
            arg_type = Float64
=#
        "--start_time"
            help = "Specify daily start time for CCFs to be processed (HH MM)"
            nargs = 2
            arg_type = Int64
            default = [0, 0]
            #default = [18, 30]
        "--stop_time"
           help = "Specify daily stop time for CCFs to be processed (HH MM)"
           nargs = 2
           arg_type = Int64
           #default = [ min_order(NEID2D()), max_order(NEID2D()) ]
           #default = [22, 30]
           default = [23, 59]
        "--start_time_clean"
            help = "Specify daily start time for CCFs for clean spectrum"
            nargs = 2
            arg_type = Int64
            default = [17, 30]
        "--stop_time_clean"
           help = "Specify daily stop time for CCFs for clean spectrum"
           nargs = 2
           arg_type = Int64
           default = [22, 12]
         "--max_num_spectra"
            help = "Maximum number of spectra to process."
            arg_type = Int
            default = 300  # Enough for one day of NEID
         "--max_num_spectra_clean"
            help = "Maximum number of spectra to include in clean spectrum."
            arg_type = Int
            default = 130  # based on 120 minutes of integration time from 55s exposures
            #default = 65  # based on 60 minutes of integration time from 55s exposures
      end

     return parse_args(s)
 end


  function parse_commandline_daily_report() 
     s = ArgParseSettings( description = "Make dailiy report from daily rvs.")
     #import_settings!(s, s_files_only, args_only=false)
     add_arg_group!(s, "Files", :argg_files)
     @add_arg_table! s begin
         "input"
             help = "Filename for input with daily rvs (csv)"
             arg_type = String
             #default = "daily_rvs.csv"
         "output"
             help = "Filename for output dailiy RVs (toml)"
             arg_type = String
             default = "daily_summary.toml"
#=
         "--template_file"
            help = "Filename with CCF template stored as 'order_ccfs' (jld2) [TODO]"
            arg_type = String
=#
         "--overwrite"
            help = "Specify it's ok to overwrite the output file."
            #default = true
            action = :store_true
      end
      add_arg_group!(s, "Filter rvs for ", :argg_filter_param)
      @add_arg_table! s begin
         "--target"
            help = "Target field"
            arg_type = String
            default = "Sun"
         "--datestr"
            help = "Filenames that contain date string."
            arg_type = String
        "--max_solar_hour_angle"
           help = "Maximum absolute value of solar hour angle."
           arg_type = Float64
           default = 3.12
         "--max_airmass"
            help = "Maximum airmass."
            arg_type = Float64
            default = 2.0
         "--min_expmeter"
            help = "Minimum exposure meter flux relative to model flux."
            arg_type = Float64
            default = 6e4
#=
         "--min_expmeter_factor"
            help = "Minimum exposure meter flux relative to model flux."
            arg_type = Float64
            default = 0.9
         "--max_expmeter_factor"
            help = "Maximum exposure meter flux relative to model flux."
            arg_type = Float64
            default = 1.1
=#
         "--max_expmeter_frac_rms"
            help = "Maximum fractional RMS of exposure meter flux."
            arg_type = Float64
            default = 0.01
#=
         "--min_pyrhelio_factor"
            help = "Minimum pyrheliometer flux relative to model flux."
            arg_type = Float64
            default = 0.9
         "--max_pyrhelio_factor"
            help = "Maximum pyrheliometer flux relative to model flux."
            arg_type = Float64
            default = 1.1
=#
         "--max_pyrhelio_frac_rms"
            help = "Maximum fractional RMS of pyrheliometer meter flux."
            arg_type = Float64
            default = 0.0035
         #=
         "--min_drp_snr"
            help = "Minimum extracted SNR reported by NEID DRP."
            arg_type = Float64
            #default = 0.0
         =#
            #=
         "--start_time"
            help = "Specify daily start time for CCFs to be processed (HH MM)"
            nargs = 2
            arg_type = Int64
            default = [0, 0]
            #default = [18, 30]
        "--stop_time"
           help = "Specify daily stop time for CCFs to be processed (HH MM)"
           nargs = 2
           arg_type = Int64
           #default = [ min_order(NEID2D()), max_order(NEID2D()) ]
           #default = [22, 30]
           default = [23, 59]
           =#
         "--max_num_spectra"
            help = "Maximum number of spectra to process."
            arg_type = Int
            default = 300   # TODO: Increase  after finish testing
      end

     return parse_args(s)
 end

 function parse_commandline_combine_daily_reports()
     s = ArgParseSettings( description = "Combined NEID solar daily reports into one csv file.")
     add_arg_group!(s, "Files", :argg_files)
     @add_arg_table! s begin
         "input_path"
             help = "Path for input files"
             arg_type = String
             default = pwd()
         "output"
             help = "Filename for output (csv)"
             arg_type = String
             default = "summary.csv"
         "--input_filename"
             help = "Filename for daily reports (toml)"
             arg_type = String
             default = "daily_summary.toml"
         "--exclude_filename"
             help = "Filename with dates to exclude (csv)"
             arg_type = String
         "--overwrite"
            help = "Specify it's ok to overwrite the output file."
            #default = true
            action = :store_true
         "--log_info"
            help = "Print info-level messages."
            action = :store_true
         "--log_debug"
            help = "Print debug-level messages."
            action = :store_true
      end

     return parse_args(s)
 end

function parse_commandline_combined_daily_rvs()
     s = ArgParseSettings( description = "Combined NEID solar daily reports into one csv file.")
     add_arg_group!(s, "Files", :argg_files)
     @add_arg_table! s begin
         "input_path"
             help = "Path for input files"
             arg_type = String
             default = pwd()
         "output"
             help = "Filename for output (csv)"
             arg_type = String
             default = "summary.csv"
         "--input_filename"
             help = "Filename for daily tables (csv)"
             arg_type = String
             default = "daily_summary.toml"
         "--exclude_filename"
             help = "Filename with dates to exclude (csv)"
             arg_type = String
         "--overwrite"
            help = "Specify it's ok to overwrite the output file."
            #default = true
            action = :store_true
         "--log_info"
            help = "Print info-level messages."
            action = :store_true
         "--log_debug"
            help = "Print debug-level messages."
            action = :store_true
      end

     return parse_args(s)
 end

