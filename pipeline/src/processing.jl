# -*- coding: utf-8 -*-
# @author: max (max@mdotmonar.ch)

include("../../lib/complexity.jl")


using HDF5
using DSP
using JSON
using Statistics
using Plots
using CurveFit

function signal_reconstruct(signal, info)
	return [(s-info[:ADZero])*(info[:ConversionFactor]*(10.0^info[:Exponent])) for s in signal]
end

function group_check(file, list, i=0)
	if i == length(list)
		return true
	end
	path = "/"*join(list[1:i], "/")
	println("Checking path: ", path)
	if list[i+1] in keys(file[path])
		return group_check(file, list, i+1)
	end
	return false
end

function high_pass_filter(signal, nco, order)
	filter = digitalfilter(Highpass(nco), Butterworth(order))
	return filtfilt(filter, signal)
end

function band_pass_filter(signal, nco_low, nco_high, order)
	filter = digitalfilter(Bandpass(nco_low, nco_high), Butterworth(order))
	return filtfilt(filter, signal)
end

# parameters
sampling_rate = 20000
resampling_rate = 250
nyquist_frequency = 0.5 * sampling_rate

function read_raw_normalize_and_average(dataset)
	h5open("./pipeline/processed_data/"*dataset*"_processed.h5", "cw") do processed_file
		# check if average_signal group exists
		if group_check(processed_file, ["average_signal"])
			println("Skipping signal averaging.")
			return
		end

		# open raw file
		file = h5open("./data/"*dataset*"_natural_images.h5", "r")
		stream = read(file, "Data/Recording_0/AnalogStream/Stream_0")
		close(file)

		# open event list file
		csv_file = readlines("./data/event_list_"*dataset*".csv")
		event_lengths = [split(line, ",")[5] for line in csv_file]
		event_lengths = event_lengths[2:end]
		event_lengths = [parse(Int, e_l) for e_l in event_lengths]

		events = [split(line, ",")[2] for line in csv_file]
		events = events[2:end]
		events = [parse(Int, e) for (i, e) in enumerate(events) if event_lengths[i] > 20000*60]

		event_lengths = [e_l for e_l in event_lengths if e_l > 20000*60]
		min_event_length = minimum(event_lengths)


		signal_length = trunc(Int, length(stream["ChannelData"])/252)	# 252 channels
		average_signal = zeros(min_event_length)
		for i in 1:252
			electrode_label = "electrode_"*stream["InfoChannel"][i][:Label]
			println("Processing electrode #"*string(i-1)*": "*electrode_label)

			info = stream["InfoChannel"][i]
			signal_raw = stream["ChannelData"][signal_length*(i-1)+1:signal_length*i]

			signal = signal_reconstruct(signal_raw, info)
			signal = [s*1000 for s in signal] # Convert from V to mV

			# avg repetitions
			avg_repetitions = zeros(min_event_length)
			for event in events
				avg_repetitions = avg_repetitions + signal[event:event+min_event_length-1]
			end
			avg_repetitions = avg_repetitions ./ length(events)

			# normalize
			u = mean(avg_repetitions)
			s = std(avg_repetitions)
			normalized_signal = (avg_repetitions .- u) ./ s

			average_signal = average_signal + normalized_signal
		end
		average_signal = average_signal ./ 252

		processed_file["average_signal/data"] = average_signal
		println("Done.")
	end
end

function filter_signal(dataset, filter)
	# open processed file
	h5open("./pipeline/processed_data/"*dataset*"_processed.h5", "cw") do processed_file
		# check if average_signal group exists
		if group_check(processed_file, ["filtered_signal"])
			println("Skipping signal filtering.")
			return
		end

		println("Filtering signal... ")
		# read average_signal
		average_signal = read(processed_file, "average_signal/data")

		# filter signal
		if filter["type"] == "HPF"
			filtered_signal = high_pass_filter(average_signal, filter["cutoff"]/nyquist_frequency, 5)
		elseif filter["type"] == "BPF"
			filtered_signal = band_pass_filter(average_signal, filter["cutoff_low"]/nyquist_frequency, filter["cutoff_high"]/nyquist_frequency, 5)
		end

		processed_file["filtered_signal/data"] = filtered_signal
		processed_file["filtered_signal/meta/filter_type"] = filter["type"]
		if filter["type"] == "HPF"
			processed_file["filtered_signal/meta/cutoff"] = filter["cutoff"]
		elseif filter["type"] == "BPF"
			processed_file["filtered_signal/meta/cutoff_low"] = filter["cutoff_low"]
			processed_file["filtered_signal/meta/cutoff_high"] = filter["cutoff_high"]
		end
		println("Done.")
	end
end

function resample_signal(dataset, resampling_rate)
	# open processed file
	h5open("./pipeline/processed_data/"*dataset*"_processed.h5", "cw") do processed_file
		# check if average_signal group exists
		if group_check(processed_file, ["signals", "complete"])
			println("Skipping signal to noise ratio selection.")
			return
		end

		println("Resampling signal... ")
		# read filtered_signal
		filtered_signal = read(processed_file, "filtered_signal/data")

		# resample signal
		step = trunc(Int, sampling_rate / resampling_rate)
		resampled_signal = filtered_signal[1:step:end]

		processed_file["signals/complete/data"] = resampled_signal
		processed_file["signals/meta/resampling_rate"] = resampling_rate
		println("Done.")
	end
end

function compute_complexity_curve(dataset, segment, type, m, r, scales)
	# open processed file
	h5open("./pipeline/processed_data/"*dataset*"_processed.h5", "cw") do processed_file
		# check if complexity group exists
		if group_check(processed_file, ["complexity", segment, type, string(r)])
			println("Skipping signal to noise ratio selection.")
			return
		end

		println("Computing complexity curve... ")

		# read signal
		signal = read(processed_file, "signals/"*segment*"/data")

		# compute complexity curve
		if type == "MSE"
			complexity_curve = multiscale_entropy(signal, m, r*std(signal), "sample", scales)
		elseif type == "RCMSE"
			complexity_curve = refined_composite_multiscale_entropy(signal, m, r*std(signal), "sample", scales)
		elseif type == "FMSE"
			complexity_curve = multiscale_entropy(signal, m, r*std(signal), "fuzzy", scales)
		elseif type == "FRCMSE"
			complexity_curve = refined_composite_multiscale_entropy(signal, m, r*std(signal), "fuzzy", scales)
		end

		processed_file["complexity/"*segment*"/"*type*"/"*string(r)*"/data"] = complexity_curve
		processed_file["complexity/"*segment*"/"*type*"/"*string(r)*"/meta/m"] = m
		processed_file["complexity/"*segment*"/"*type*"/"*string(r)*"/meta/r"] = r
		processed_file["complexity/"*segment*"/"*type*"/"*string(r)*"/meta/scales"] = scales
	end
end

function compute_linear_regression(dataset, segment, type, r)
	# open processed file
	h5open("./pipeline/processed_data/"*dataset*"_processed.h5", "cw") do processed_file
		# check if complexity group exists
		if group_check(processed_file, ["complexity", segment, type, string(r), "meta", "linear_fit"])
			println("Skipping linear regression.")
			return
		end

		println("Computing linear regression... ")

		# read complexity curve
		complexity_curve = read(processed_file, "complexity/"*segment*"/"*type*"/"*string(r)*"/data")
		scales = read(processed_file, "complexity/"*segment*"/"*type*"/"*string(r)*"/meta/scales")

		# compute linear regression
		a, b = linear_fit(scales, complexity_curve)

		processed_file["complexity/"*segment*"/"*type*"/"*string(r)*"/meta/linear_fit/a"] = a
		processed_file["complexity/"*segment*"/"*type*"/"*string(r)*"/meta/linear_fit/b"] = b
	end
end