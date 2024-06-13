# -*- coding: utf-8 -*-
# @author: max (max@mdotmonar.ch)

using HDF5
using Plots
using StatsPlots
using DSP

include("../../lib/complexity.jl")

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
	println("Plotting dataset: ", dataset)
	# if directory does not exist, create it
	if !isdir("./pipeline/plots/"*dataset*"/")
		mkpath("./pipeline/plots/"*dataset*"/")
	end

	# open processed file
	h5open("./pipeline/processed_data/"*dataset*"_processed.h5", "r") do processed_file

		for segment in ["complete"]
			plot(size=(1024, 768))
			signal = read(processed_file, "signals/"*segment*"/data")
			plot!((0:length(signal)-1)./250, signal, label="Response")
			plot!(xlabel="Time", ylabel="Amplitude", title=dataset)
			plot!(ylims=(-3, 3))
			savefig("./pipeline/plots/"*dataset*"/"*dataset*"_signal_"*segment*".png")

			# read entropy curves
			plot(size=(1024, 768))
			for r in ["0.2"]
				rcmse_curve = read(processed_file, "complexity/"*segment*"/RCMSE/"*string(r)*"/data")
				rcmse_scales = read(processed_file, "complexity/"*segment*"/RCMSE/"*string(r)*"/meta/scales")
				plot!(rcmse_scales, rcmse_curve, label="RCMSE_"*string(r))
			end
			plot!(xlabel="Scale", ylabel="Entropy", title=dataset)
			#y label from 0 to 3
			plot!(ylims=(0, 3))
			savefig("./pipeline/plots/"*dataset*"/"*dataset*"_rcmse_"*segment*".png")
		end
	end
end