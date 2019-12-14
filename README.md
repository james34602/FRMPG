# FRMPG
## FR minimum phase general filter design

smoothSpectrumImpulseResponse.m is a function can perform N-octave smoothing from impulse response.

minphasefir2.m is modified fir2 function that perform grid interpolation and frequency sampling method to design minimum phase filter.

demo.m show you how to convert a impulse response file from waveform to shrinked interesting points and reconvert it to minimum phase impulse response.
The intermediate result [freq,gV] is extracted interesting point.

demo2.m show you how to craft frequency response from handwritten numeric vectors.

# What's the point?

Headphone frequnecy response equalization should not be generated using simple waveform editing software, or worse, record the "impulse response" from loopback driver, which is quite awful. And possibly create unintentional signal delay, resampling artifacts.
