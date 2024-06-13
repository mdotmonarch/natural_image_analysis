# -*- coding: utf-8 -*-
# @author: max (max@mdotmonar.ch)

using HDF5
using Plots
using StatsPlots

include("../../lib/complexity.jl")
include("processing.jl")

function signal_reconstruct(signal, info)
	return [(s-info[:ADZero])*(info[:ConversionFactor]*(10.0^info[:Exponent])) for s in signal]
end

# parameters
filter = Dict(
	"type" => "BPF",
	"cutoff_low" => 0.1, # Hz
	"cutoff_high" => 40 # Hz
)
resampling_rate = 250 # Hz

datasets = [
	"MR-0554",
	"MR-0569_nd4",
	"MR-0596_nd4",
	"MR-0597_nd4",
	"MR-0599_nd4",
	"MR-0622_nd4",
	"MR-0623_nd4",
	"MR-0624_nd4",
	"MR-0625_nd4",
	"MR-0591_nd4",
	"MR-0592_nd4",
	"MR-0620_nd4",
	"MR-0621_nd4",
	"MR-0593_nd4",
	"MR-0552",
	"MR-0629-t1_nd4",
	"MR-0629-t2_nd4",
	"MR-0630_nd4",
	"MR-0548-t1",
	"MR-0548-t2",
	"MR-0568-t1_nd4",
	"MR-0568-t2_nd4",
	"MR-0582_nd4",
	"MR-0587_nd4",
	"MR-0588_nd4",
	"MR-0583_nd4",
	"MR-0575_nd4",
	"MR-0577_nd4",
	"MR-0579_nd4",
	"MR-0585_nd4",
	"MR-0586_nd4",
]

for dataset in datasets
	println("Processing dataset: ", dataset)
	read_raw_normalize_and_average(dataset)
	filter_signal(dataset, filter)
	resample_signal(dataset, resampling_rate)

	compute_complexity_curve(dataset, "complete", "RCMSE", 2, 0.2, [i for i in 1:45])
	compute_linear_regression(dataset, "complete", "RCMSE", 0.2)
end