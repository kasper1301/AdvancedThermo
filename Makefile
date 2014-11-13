

all: experimentalPhaseEnvelope pressureVolume main


experimentalPhaseEnvelope: src/plotPhaseEnvelope.m
	cd src/; matlab -nodisplay -nosplash -r tmp;

pressureVolume: src/pressureVolumePlot.jl
	cd src/; julia pressureVolumePlot.jl

main: main.jl
	cd src; julia main.jl
