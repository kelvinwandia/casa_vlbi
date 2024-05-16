--[[
 This is the generic "minimal" AOFlagger strategy, version 2020-06-14
 Author: André Offringa

 It is functionally equal to the default generic strategy, but removes some of the
 tweaking parameters and visualization. It is therefore less easy to tweak, but smaller
 and possibly slightly faster.
]]

aoflagger.require_min_version("3.0")

function execute(input)
  --
  -- Generic settings
  --

  local base_threshold = 0.75-- lower means more sensitive detection
  -- How to flag complex values, options are: phase, amplitude, real, imaginary, complex
  local representation = "amplitude"
  local iteration_count = 5 -- how many iterations to perform?
  local threshold_factor_step = 2.5 -- How much to increase the sensitivity each iteration?
  local transient_threshold_factor = 0.75 -- decreasing this value makes detection of transient RFI more aggressive

  --
  -- End of generic settings
  --

  local inpPolarizations = input:get_polarizations()

  input:clear_mask()

  for ipol, polarization in ipairs(inpPolarizations) do
    local converted_data = input:convert_to_polarization(polarization):convert_to_complex(representation)

    local converted_copy = converted_data:copy()

    for i = 1, iteration_count - 1 do
        local threshold_factor = threshold_factor_step ^ (iteration_count - i)

        local x_threshold_factor = threshold_factor * base_threshold
        -- Runs sumthreshold and detects sharp, line-shaped features in the time-freq domain
        -- Parameters
        -- x_threshold_factor - threshold_factor in time direction
        -- y_threshold_factor - threshold_factor in frequency direction
        -- x_direction (bool) - enable flagging in time direction
        -- y_direction (bool) - enable flagging in frequency direction
        aoflagger.sumthreshold(converted_data, x_threshold_factor, x_threshold_factor * transient_threshold_factor, true, true)

        -- Do timestep & channel flagging
        -- Calculates the root mean square RMS for each timestep and flags channels with outlier
        local chdata = converted_data:copy()
        aoflagger.threshold_timestep_rms(converted_data, 3.5)
        aoflagger.threshold_channel_rms(chdata, 3.0 * threshold_factor, true)
        converted_data:join_mask(chdata)

        -- High pass filtering steps 
        --
        -- high_pass_filter removes the diffuse background in the resized_data
        --
        -- Parameters
            --  data - the data modified in place
            --  xsize (integer): kernel size in time direction
            --  ysize - kernel size in frequency direction
            --  xsigma - Gaussian width in time direction
            --  ysigma - Gaussian size in frequency direction
            --  removes astromomical signals before thresholding
            --  for high resolution data, use a large kernal either in frequency and/or time

        local xsize = 51
        local ysize = 51
        local xsigma = 3.0
        local ysigma = 3.0
        
        --  Downsampling
        --  Decreases the resolution of the data using simple linear binning - which can increase
        --  the speed of data smoothing
        --  Parameters
            -- data - the input data
            -- xfactor(int) - the downsampling in time direction
            -- yfactor(int) - the downsampling in frequency direction
            -- masked(bool) - take flags into account when averaging

        local time_resize_factor = 1 -- allows convolution with a smaller kernel
        local frequency_resize_factor = 1

        converted_data:set_visibilities(converted_copy)
        local resized_data = aoflagger.downsample(converted_data, time_resize_factor, frequency_resize_factor, true)
        aoflagger.low_pass_filter(resized_data, xsize, ysize, xsigma, ysigma)
        aoflagger.upsample(resized_data, converted_data, time_resize_factor, frequency_resize_factor)


        local tmp = converted_copy - converted_data
        tmp:set_mask(converted_data)
        converted_data = tmp

        aoflagger.set_progress((ipol - 1) * iteration_count + i, #inpPolarizations * iteration_count)

    end -- end of iterations

    aoflagger.sumthreshold(converted_data, base_threshold, base_threshold * transient_threshold_factor, true, true)

    if input:is_complex() then
        converted_data = converted_data:convert_to_complex("complex")
    end
    input:set_polarization_data(polarization, converted_data)

    aoflagger.set_progress(ipol, #inpPolarizations)
  end -- end of polarization iterations

    aoflagger.scale_invariant_rank_operator(input, 0.2, 0.2)
    aoflagger.threshold_timestep_rms(input, 4.0)
    input:flag_nans()
end
