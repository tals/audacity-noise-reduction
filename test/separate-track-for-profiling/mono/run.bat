set input=square-with-noise.wav
set output=square-with-noise-output.wav
set profile=noise-profile.wav
set gain=24
set sensitivity=12
set smoothing=0
..\..\..\build\Debug\noisereduction_driver.exe -i %input% -o %output% -p %profile% --noiseGain %gain% --sensitivity %sensitivity% --smoothing %smoothing%
PAUSE
