set input=square-with-noise.wav
set output=square-with-noise-output.wav
set t0=0
set t1=44100
set gain=24
set sensitivity=12
set smoothing=0
..\..\..\build\Debug\noisereduction_driver.exe -i %input% -o %output% --t0 %t0% --t1 %t1% --noiseGain %gain% --sensitivity %sensitivity% --smoothing %smoothing%
PAUSE
