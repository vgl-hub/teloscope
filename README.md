# teloscope
**Tool Use**

This tool uses a sliding window to detect the occurences of a given string in a sequence. The teloscope can be used to find certain chains of base pairs such as the 'TTAGGG' often found in the telomeres of vertebrates. Given the input of a DNA sequence, the sequence being detected, and specifications for the window, the teloscope will produce a vector that contains the number of times a sequences is in the following window. The best steps to use (from personal experience) would be a step with a quantity of 1 or 3.

**Installation**

Either download one of the releases or "git clone https://github.com/vgl-hub/gfastats.git --recursive" and "make -j" in "gfastats" folder.