# Callmark

`Callmark` is a Matlab GUI tool for marking calls or other sounds in
multi-recording from animals.
It allow user to switch between different channels while editing
the marked call locations.
The marked call locations can then be saved as a separate file
and loaded again for further editing.

## Input
Input to the program can be a `WAV` or a `MAT` file. If a `MAT` file,
it should contain `sig` (sound recording) and `fs` (sampling frequency in Hz).
Below is an example of the structure of an input file that has 32 channels,
4 seconds of recording, sampled at 250 kHz:
```
Name         Size          Bytes       Class     Attributes
  fs           1x1              8      double              
  sig    1000002x32     256000512      double        
```

## GUI Operation

### GUI button functions

*  **Load file**<br>
  Load data file in `WAV` or `MAT` format.
  If this is the first time you load a particular file, you will see
  two file/folder selection windows with prompts
  1) `Select the folder containing the recording files` and
  2) `Select the folder containing *_detect.mat files`.
  These windows won't pop up if you try to load another recording file
  in the same folder.

  After selecting these paths, you will see another window popping up
  for you to select the actual recording file you want to analyze.

* **Detect calls**<br>
  Initial simple energy detection of calls,
  which will also initialize the spectrogram and time series display.
  Note that there is usually a bit of lag after you accept the simple
  threshold detection until the spectrogram shows up, so don't be
  alarmed if nothing seems to be happening for a few seconds.

* **Add mark**<br>
  Add a call mark on where you click.

* **Remove mark**<br>
  Remove all marks enclosed in a dragged rectangular area on the spectrogram.

* **Done & Save**<br>
  Save the marked calls in a new file with default filename:
  `xxx_detect.mat`, where `xxx` is the filename of the loaded data file.
  Therefore, the call marking results and the raw data are kept as
  separate files.

### Other operation notes
* The MATLAB figure tools `Zoom in`, `Zoom out`, and `Pan` are active
in this program so that you can move around easily in the recorded
time series.
* The scroll bar on the right-hand side of the spectrogram is used
to control the coloraxis of spectrogram.
* After saving `xxx_detect.mat`, you can use **Load file** to load
either `xxx_detect.mat` or `xxx.mat` to review what has been saved.
