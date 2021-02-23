global max_spectra_to_use = 250
if max_spectra_to_use < 250
   @warn "param.in setting max_spectra_to_use to " * string(max_spectra_to_use)
end

global idx_day_to_use
#global tophap_ccf_mask_scale_factor=1.6

global fits_target_str
@assert fits_target_str == "Solar" || fits_target_str == "Sun"

global linelist_for_ccf_filename = "G2.espresso.mas"

hostname = gethostname()
   if occursin("aci.ics.psu.edu",hostname)
      global ancilary_solar_data_path = "/gpfs/group/ebf11/default/ebf11/neid_solar"
   elseif occursin("nuc8",hostname)  # Eric's home machine :)
      global ancilary_solar_data_path = "/home/eford/Data/SolarSpectra/NEID_solar"
   end
   global ccf_mid_velocity = 0

#global bjd_first_good = 2458745.1296134139
   #global bjd_last_good = 2458745.283
   # Good days identified by Andrea
   good_days = ["20201215", "20201216", "20201218", "20201219", "20201220", "20201226", "20201230", "20210104", "20210109", "20210110", "20210111", "20210112", "20210115", "20210118" ]
   very_good_days = ["20201216", "20201220", "20210109", "20210110", "20210111", "20210112" ]
   function is_good_day(fn::AbstractString)
      m = match(r"neidL1_(\d+)T\d+\.fits", fn)
      if m != Nothing
         substr = m.captures[1]
         return any(map(d->d==substr,good_days))
      else
         return false
      end
   end

global df_files
   global df_files_use = df_files |>
      @filter( _.target == fits_target_str ) |>
      #@filter(bjd_first_good <= _.bjd < bjd_last_good) |>
      @filter( is_good_day(_.Filename) ) |>
      #@take(max_spectra_to_use) |>
      DataFrame
   global df_files_solar_by_day = df_files_use |> @groupby(floor(Int64,_.bjd)) |> @map({obsjd_int=key(_), data=_ } ) |> DataFrame
   if size(df_files_use,1) > max_spectra_to_use
      df_files_use = df_files_solar_by_day[idx_day_to_use,:data] |> @orderby(_.bjd) |> @take(max_spectra_to_use) |> DataFrame
   end
